#!/usr/bin/env Rscript
# build_diagrams.R -- turn a diagram spec (specs/<name>.yaml) into an Excalidraw
# scene (out/<name>.excalidraw). Stage 1 of 2; stage 2 (export.mjs) renders the
# scene to SVG.
#
# The spec carries no pixel coordinates. Layout is computed here on a grid:
#   y band  <- the node's `layer`, in the order the spec's `layers` block gives
#   x column <- longest-path rank over the edge list, so a callee sits to the
#               right of its caller and the pipeline reads left to right
#   row     <- position within (band, column), in spec order
# A node may override its placement with `col:` / `row:`; nothing else about
# position is expressible in the spec.
#
# Output is reproducible: ids are derived from node ids, and the RNG for `seed`
# and `versionNonce` is seeded from the spec name, so rebuilding an unchanged
# spec gives a byte-identical file and an empty git diff.
#
# Usage: Rscript build_diagrams.R [specs/foo.yaml ...]   (default: all specs)

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

# The `kind: workflow` layout engine lives in its own file. It defines no
# elements and draws no random numbers at load time, so sourcing it cannot move
# the architecture kind's output.
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])[1]
  here <- if (is.na(f)) getwd() else dirname(normalizePath(f))
  source(file.path(here, "workflow_layout.R"), local = FALSE)
})

# ---------------------------------------------------------------- geometry ----

NODE_W       <- 300    # node box width
NODE_MIN_H   <- 58     # node box height for a single-line label
LINE_H       <- 21     # text line advance at FONT_SIZE
FONT_SIZE    <- 14
FONT_FAMILY  <- 3      # code face: predictable advance width, suits identifiers
CHAR_W       <- 8.4    # ~0.6 * FONT_SIZE, the monospace advance
COL_GAP      <- 44     # horizontal gap between columns
ROW_GAP      <- 30     # vertical gap between rows inside a band
BAND_PAD_X   <- 26     # band padding around its nodes
BAND_PAD_TOP <- 46     # extra room at the top of a band for its title
BAND_PAD_BOT <- 22
BAND_GAP     <- 34     # vertical gap between bands
ORIGIN       <- c(x = 0, y = 0)

# Status -> stroke / fill. `built: false` additionally dashes the stroke, and
# `emphasis: true` thickens it.
STATUS_STYLE <- list(
  new              = list(stroke = "#1971c2", fill = "#a5d8ff"),
  reused_unchanged = list(stroke = "#495057", fill = "#e9ecef"),
  external         = list(stroke = "#e8590c", fill = "#ffec99")
)
EDGE_STYLE <- list(
  call     = list(stroke = "#6c757d", style = "solid", width = 1),
  data     = list(stroke = "#1971c2", style = "dashed", width = 2),
  dispatch = list(stroke = "#868e96", style = "dotted", width = 1)
)
INFERRED_STROKE <- "#9c36b5"   # edges the record does not state (stated: false)
BAND_FILL       <- "#f8f9fa"
BAND_STROKE     <- "#ced4da"

# An edge label longer than this is moved to the numbered key below the diagram.
# Dense graphs cannot carry twenty mid-air labels without them colliding with
# each other and with the boxes; short ones stay inline because a bound arrow
# label gaps the line it sits on.
LABEL_INLINE_MAX <- 20
BADGE_R          <- 11     # radius of a numbered edge-label badge

# ------------------------------------------------------------------ helpers ---

`%||%` <- function(a, b) if (is.null(a)) b else a

rnd31 <- function() as.integer(floor(stats::runif(1) * 2147483647))

# Excalidraw ids are opaque strings; derive them from spec ids so that a binding
# failure names the node it came from and so rebuilds are stable.
el_id <- function(...) {
  s <- paste(c(...), collapse = "--")
  gsub("[^A-Za-z0-9_-]+", "-", s)
}

# Greedy wrap on whitespace, then hard-break any word longer than the budget.
wrap_label <- function(s, max_chars) {
  words <- strsplit(s, " +")[[1]]
  out <- character(0)
  cur <- ""
  push <- function() if (nzchar(cur)) out <<- c(out, cur)
  for (w in words) {
    while (nchar(w) > max_chars) {
      push(); cur <- ""
      out <- c(out, substr(w, 1, max_chars))
      w <- substr(w, max_chars + 1, nchar(w))
    }
    cand <- if (nzchar(cur)) paste(cur, w) else w
    if (nchar(cand) <= max_chars) {
      cur <- cand
    } else {
      push(); cur <- w
    }
  }
  push()
  if (!length(out)) out <- ""
  out
}

# Fields every element carries. `updated` is taken from the spec's mtime, not
# Sys.time(), to keep the output reproducible.
# `roughness` and `opacity` are appended with the values they replace as
# defaults, so every existing call site -- and the bytes it produces -- is
# unaffected. The workflow kind passes roughness 0 (clean lines) and opacity 55
# (a stage specified but not yet coded).
base_el <- function(id, type, x, y, w, h, stroke, fill,
                    stroke_style = "solid", stroke_width = 1,
                    fill_style = "solid", roundness = 3L, updated,
                    roughness = 1L, opacity = 100L) {
  list(
    id              = id,
    type            = type,
    x               = x,
    y               = y,
    width           = w,
    height          = h,
    angle           = 0,
    strokeColor     = stroke,
    backgroundColor = fill,
    fillStyle       = fill_style,
    strokeWidth     = stroke_width,
    strokeStyle     = stroke_style,
    roughness       = roughness,
    opacity         = opacity,
    groupIds        = list(),
    frameId         = NULL,
    roundness       = if (is.null(roundness)) NULL else list(type = roundness),
    seed            = rnd31(),
    version         = 1L,
    versionNonce    = rnd31(),
    isDeleted       = FALSE,
    boundElements   = list(),
    updated         = updated,
    link            = NULL,
    locked          = FALSE
  )
}

# As for base_el(), the trailing parameters default to the values they replace.
# `auto_resize = FALSE` matters for a generated scene: with TRUE, Excalidraw
# regrows the container when the file is opened, which would silently undo a
# computed layout.
text_el <- function(id, container_id, lines, x, y, w, h, colour, size, updated,
                    align = "center", family = FONT_FAMILY, roughness = 1L,
                    valign = "middle", auto_resize = TRUE, opacity = 100L) {
  txt <- paste(lines, collapse = "\n")
  el <- base_el(id, "text", x, y, w, h, colour, "transparent",
                roundness = NULL, updated = updated,
                roughness = roughness, opacity = opacity)
  el$fontSize      <- size
  el$fontFamily    <- family
  el$text          <- txt
  el$originalText  <- txt
  el$textAlign     <- align
  el$verticalAlign <- valign
  el$containerId   <- container_id
  el$lineHeight    <- 1.25
  el$autoResize    <- auto_resize
  el
}

# -------------------------------------------------------------------- layout --

compute_layout <- function(spec) {
  nodes <- spec$nodes
  ids <- vapply(nodes, `[[`, character(1), "id")
  if (anyDuplicated(ids)) {
    stop("duplicate node ids: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "))
  }

  band_ids <- vapply(spec$layers, `[[`, character(1), "id")
  node_band <- vapply(nodes, `[[`, character(1), "layer")
  bad <- setdiff(node_band, band_ids)
  if (length(bad)) stop("node(s) in undeclared layer: ", paste(bad, collapse = ", "))

  # --- x: longest-path rank over every edge, relaxed iteratively so that a
  # cycle (which this spec does not have, but a later one might) degrades to a
  # capped rank instead of hanging.
  from <- vapply(spec$edges, `[[`, character(1), "from")
  to   <- vapply(spec$edges, `[[`, character(1), "to")
  unknown <- setdiff(c(from, to), ids)
  if (length(unknown)) {
    stop("edge endpoint is not a node: ", paste(unknown, collapse = ", "))
  }
  rank <- setNames(rep(0L, length(ids)), ids)
  for (sweep in seq_along(ids)) {
    changed <- FALSE
    for (e in seq_along(from)) {
      want <- rank[[from[e]]] + 1L
      if (rank[[to[e]]] < want) {
        rank[[to[e]]] <- want
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  # Nodes with no edges at all are reference material (the record names them but
  # nothing in the GSD path calls them). Park them in a column of their own on
  # the right rather than at rank 0 on the left.
  detached <- setdiff(ids, union(from, to))
  if (length(detached)) rank[detached] <- max(rank) + 1L

  col <- rank[ids]
  for (i in seq_along(nodes)) {
    if (!is.null(nodes[[i]]$col)) col[i] <- as.integer(nodes[[i]]$col)
  }

  # --- label wrap and box height
  max_chars <- floor((NODE_W - 24) / CHAR_W)
  lines <- lapply(nodes, function(n) wrap_label(n$label %||% n$id, max_chars))
  h <- vapply(lines, function(l) max(NODE_MIN_H, length(l) * LINE_H + 22), numeric(1))

  # --- row within (band, column), in spec order unless overridden
  row <- integer(length(ids))
  for (b in band_ids) {
    for (cc in sort(unique(col[node_band == b]))) {
      sel <- which(node_band == b & col == cc)
      ord <- seq_along(sel)
      given <- vapply(nodes[sel], function(n) as.integer(n$row %||% NA_integer_), integer(1))
      if (any(!is.na(given))) ord <- order(ifelse(is.na(given), ord * 100L, given))
      row[sel[ord]] <- seq_along(sel) - 1L
    }
  }

  # --- pixels. Column x is global so bands align vertically. Row pitch is per
  # (band, row) so a tall three-line label does not overlap the row below.
  col_x <- setNames(numeric(0), character(0))
  for (cc in sort(unique(col))) {
    col_x[as.character(cc)] <- ORIGIN[["x"]] + cc * (NODE_W + COL_GAP)
  }

  x <- col_x[as.character(col)]
  y <- numeric(length(ids))
  bands <- list()
  cursor <- ORIGIN[["y"]]
  for (b in band_ids) {
    sel <- which(node_band == b)
    if (!length(sel)) next
    nrow_b <- max(row[sel]) + 1L
    row_h <- vapply(seq_len(nrow_b) - 1L, function(r) {
      hh <- h[sel][row[sel] == r]
      if (length(hh)) max(hh) else NODE_MIN_H
    }, numeric(1))
    row_y <- cursor + BAND_PAD_TOP + c(0, cumsum(row_h + ROW_GAP))[seq_len(nrow_b)]
    y[sel] <- row_y[row[sel] + 1L]
    band_h <- BAND_PAD_TOP + sum(row_h) + ROW_GAP * (nrow_b - 1L) + BAND_PAD_BOT
    bands[[b]] <- list(
      id = b,
      label = spec$layers[[which(band_ids == b)]]$label %||% b,
      x = min(x[sel]) - BAND_PAD_X,
      y = cursor,
      w = max(x[sel] + NODE_W) - min(x[sel]) + 2 * BAND_PAD_X,
      h = band_h
    )
    cursor <- cursor + band_h + BAND_GAP
  }
  # Square the bands off to a common left edge and width so they read as bands.
  bx <- min(vapply(bands, `[[`, numeric(1), "x"))
  bw <- max(vapply(bands, function(b) b$x + b$w, numeric(1))) - bx
  bands <- lapply(bands, function(b) { b$x <- bx; b$w <- bw; b })

  list(
    nodes = data.frame(
      id = ids, band = node_band, col = as.integer(col), row = row,
      x = as.numeric(x), y = y, w = NODE_W, h = h,
      stringsAsFactors = FALSE
    ),
    lines = lines,
    bands = bands,
    canvas = c(x = bx, y = ORIGIN[["y"]], w = bw, h = cursor - BAND_GAP)
  )
}

# ------------------------------------------------------------------- arrows ---

# Straight arrow between the two box centres, clipped by Excalidraw's own
# binding at render time (hence focus 0 / gap 4 rather than hand-clipped points).
arrow_geometry <- function(a, b) {
  ax <- a$x + a$w / 2; ay <- a$y + a$h / 2
  bx <- b$x + b$w / 2; by <- b$y + b$h / 2
  list(x = ax, y = ay, dx = bx - ax, dy = by - ay)
}

# --------------------------------------------------------------- scene build --

build_scene <- function(spec, spec_path) {
  set.seed(sum(utf8ToInt(basename(spec_path))))
  updated <- as.numeric(file.mtime(spec_path)) * 1000
  updated <- round(updated)

  lay <- compute_layout(spec)
  nd <- lay$nodes
  rownames(nd) <- nd$id
  elements <- list()

  # Z-order is array order, and it is the main thing standing between a dense
  # call graph and an unreadable one:
  #   bands -> arrows -> nodes and their text -> badges -> legend and key
  # Arrows under the boxes means a long arrow crossing the middle of the picture
  # never obscures a node's label; the binding gap keeps its arrowhead visible.

  # 1. layer bands first, so they paint behind everything
  for (b in lay$bands) {
    bid <- el_id("band", b$id)
    band <- base_el(bid, "rectangle", b$x, b$y, b$w, b$h,
                    BAND_STROKE, BAND_FILL, stroke_style = "dashed",
                    fill_style = "solid", roundness = 2L, updated = updated)
    tid <- el_id("band", b$id, "label")
    band$boundElements <- list(list(id = tid, type = "text"))
    elements <- c(elements, list(band))
    # Band title sits in the padding strip above the first row of nodes; it is a
    # free text element (left-aligned) rather than the band's centred label.
    ttl <- text_el(tid, bid, b$label, b$x + 18, b$y + 12, b$w - 36, 24,
                   "#495057", 16, updated, align = "left")
    ttl$verticalAlign <- "top"
    elements <- c(elements, list(ttl))
  }

  # 2. edges. Long labels become a numbered badge plus a key entry; short ones
  # stay bound to the arrow.
  badges <- list()
  key <- list()
  inline_labels <- list()
  for (e in spec$edges) {
    kind <- e$kind %||% "call"
    es <- EDGE_STYLE[[kind]]
    if (is.null(es)) stop("edge ", e$from, " -> ", e$to, " has unknown kind '", kind, "'")
    inferred <- identical(e$stated, FALSE)
    geo <- arrow_geometry(nd[e$from, ], nd[e$to, ])
    aid <- el_id("e", e$from, e$to)
    arr <- base_el(aid, "arrow", geo$x, geo$y, abs(geo$dx), abs(geo$dy),
                   if (inferred) INFERRED_STROKE else es$stroke, "transparent",
                   stroke_style = if (inferred) "dotted" else es$style,
                   stroke_width = es$width, roundness = 2L, updated = updated)
    arr$points       <- list(c(0, 0), c(geo$dx, geo$dy))
    arr$lastCommittedPoint <- NULL
    arr$startBinding <- list(elementId = el_id("n", e$from), focus = 0, gap = 4)
    arr$endBinding   <- list(elementId = el_id("n", e$to),   focus = 0, gap = 4)
    arr$startArrowhead <- NULL
    arr$endArrowhead   <- "arrow"
    arr$elbowed        <- FALSE

    lbl <- e$label
    inline <- !is.null(lbl) && nchar(lbl) <= LABEL_INLINE_MAX
    if (inline) {
      # Free text on a white plate, not text bound to the arrow: Excalidraw
      # renders a bound label with its container, so a bound label inherits the
      # arrow's z-order and ends up half-hidden under a box whatever position
      # its own element takes in the array.
      ln <- wrap_label(lbl, 22)
      lw <- max(nchar(ln)) * 7.2 + 10
      lh <- length(ln) * 15 + 6
      elements <- c(elements, list(arr))
      inline_labels[[length(inline_labels) + 1L]] <- list(
        id = el_id("l", e$from, e$to), lines = ln,
        x = geo$x + geo$dx / 2 - lw / 2, y = geo$y + geo$dy / 2 - lh / 2,
        w = lw, h = lh
      )
    } else {
      elements <- c(elements, list(arr))
      if (!is.null(lbl)) {
        n <- length(key) + 1L
        badges[[length(badges) + 1L]] <- list(n = n, geo = geo, edge = e)
        key[[n]] <- list(n = n, from = e$from, to = e$to, label = lbl)
      }
    }
  }

  # 3. nodes
  for (i in seq_along(spec$nodes)) {
    n <- spec$nodes[[i]]
    g <- nd[n$id, ]
    st <- STATUS_STYLE[[n$status]]
    if (is.null(st)) stop("node '", n$id, "' has unknown status '", n$status, "'")
    dashed <- identical(n$built, FALSE)
    nid <- el_id("n", n$id)
    tid <- el_id("t", n$id)
    box <- base_el(nid, "rectangle", g$x, g$y, g$w, g$h,
                   st$stroke, st$fill,
                   stroke_style = if (dashed) "dashed" else "solid",
                   stroke_width = if (isTRUE(n$emphasis)) 4 else 2,
                   roundness = 3L, updated = updated)
    box$boundElements <- list(list(id = tid, type = "text"))
    elements <- c(elements, list(box))
    ln <- lay$lines[[i]]
    th <- length(ln) * LINE_H
    elements <- c(elements, list(text_el(
      tid, nid, ln,
      g$x + 12, g$y + (g$h - th) / 2, g$w - 24, th,
      "#1e1e1e", FONT_SIZE, updated
    )))
  }

  # 4. inline arrow labels and badges for the keyed ones, both on top of the
  # nodes so neither can be painted over.
  for (il in inline_labels) {
    plate <- base_el(el_id(il$id, "plate"), "rectangle", il$x, il$y, il$w, il$h,
                     "transparent", "#ffffff", stroke_width = 1,
                     roundness = NULL, updated = updated)
    txt <- text_el(il$id, NULL, il$lines, il$x + 5, il$y + 3,
                   il$w - 10, il$h - 6, "#495057", 12, updated)
    txt$containerId <- NULL
    elements <- c(elements, list(plate), list(txt))
  }

  # badges: nudge along the arrow if one lands on another already placed.
  # A badge may not sit on another badge, on a node box, or on an inline label:
  # slide it along its own arrow until it finds clear air, and fall back to the
  # midpoint if the whole line is congested.
  occupied <- c(
    lapply(seq_len(nrow(nd)), function(i) c(nd$x[i], nd$y[i], nd$w[i], nd$h[i])),
    lapply(inline_labels, function(l) c(l$x, l$y, l$w, l$h))
  )
  hits <- function(cx, cy) {
    any(vapply(occupied, function(r) {
      cx < r[1] + r[3] && r[1] < cx + 2 * BADGE_R &&
        cy < r[2] + r[4] && r[2] < cy + 2 * BADGE_R
    }, logical(1)))
  }
  fracs <- c(0.5, 0.62, 0.38, 0.72, 0.28, 0.82, 0.18, 0.9, 0.1)
  for (b in badges) {
    g <- b$geo
    frac <- 0.5
    for (f in fracs) {
      cx <- g$x + g$dx * f - BADGE_R
      cy <- g$y + g$dy * f - BADGE_R
      if (!hits(cx, cy)) {
        frac <- f
        break
      }
    }
    cx <- g$x + g$dx * frac - BADGE_R
    cy <- g$y + g$dy * frac - BADGE_R
    occupied[[length(occupied) + 1L]] <- c(cx, cy, 2 * BADGE_R + 4, 2 * BADGE_R + 4)
    bid <- el_id("badge", b$n)
    tid <- el_id("badge", b$n, "text")
    circ <- base_el(bid, "ellipse", cx, cy, 2 * BADGE_R, 2 * BADGE_R,
                    "#495057", "#ffffff", stroke_width = 1,
                    roundness = NULL, updated = updated)
    circ$boundElements <- list(list(id = tid, type = "text"))
    elements <- c(elements, list(circ))
    elements <- c(elements, list(text_el(
      tid, bid, as.character(b$n),
      cx + 3, cy + BADGE_R - 8, 2 * BADGE_R - 6, 16, "#495057", 12, updated
    )))
  }

  # 5. legend and the label key, below the last band
  elements <- c(elements, legend_elements(lay, updated))
  elements <- c(elements, key_elements(key, lay, updated))

  # Bind arrows back onto the shapes they touch: Excalidraw needs the shape's
  # boundElements to list every arrow bound to it, not just its label.
  by_id <- setNames(seq_along(elements), vapply(elements, `[[`, character(1), "id"))
  for (e in spec$edges) {
    aid <- el_id("e", e$from, e$to)
    for (endpoint in c(e$from, e$to)) {
      k <- by_id[[el_id("n", endpoint)]]
      elements[[k]]$boundElements <- c(elements[[k]]$boundElements,
                                       list(list(id = aid, type = "arrow")))
    }
  }

  list(
    scene = list(
      type     = "excalidraw",
      version  = 2L,
      source   = "multigrain/dev/diagrams",
      elements = elements,
      appState = list(gridSize = 20L, viewBackgroundColor = "#ffffff"),
      files    = setNames(list(), character(0))
    ),
    layout = lay
  )
}

# The numbered key: one line per edge label too long to sit on its arrow.
key_elements <- function(key, lay, updated) {
  if (!length(key)) return(list())
  x0 <- lay$canvas[["x"]] + 1020
  y0 <- lay$canvas[["y"]] + lay$canvas[["h"]] + 26
  out <- list(text_el(el_id("keyhead"), NULL, "Edge labels", x0, y0 - 24, 400, 20,
                      "#495057", 14, updated, align = "left"))
  out[[1]]$containerId <- NULL
  per_col <- ceiling(length(key) / 2)
  for (k in key) {
    i <- k$n - 1L
    cx <- x0 + (i %/% per_col) * 700
    cy <- y0 + (i %% per_col) * 32
    lines <- wrap_label(sprintf("%d  %s -> %s: %s", k$n, k$from, k$to, k$label), 88)
    el <- text_el(el_id("key", k$n), NULL, lines, cx, cy, 690,
                  length(lines) * 15, "#495057", 11, updated, align = "left")
    el$containerId <- NULL
    out <- c(out, list(el))
  }
  out
}

legend_elements <- function(lay, updated) {
  items <- list(
    list(kind = "box", text = "new on gsd-build",            status = "new",              built = TRUE),
    list(kind = "box", text = "reused unchanged",            status = "reused_unchanged", built = TRUE),
    list(kind = "box", text = "external dependency",         status = "external",         built = TRUE),
    list(kind = "box", text = "dashed border: not yet coded (P3-P5)", status = "new",      built = FALSE),
    list(kind = "line", text = "call",                       style = "call"),
    list(kind = "line", text = "object handed on",           style = "data"),
    list(kind = "line", text = "S3 dispatch",                style = "dispatch"),
    list(kind = "line", text = "call site not stated by the record", style = "inferred")
  )
  x0 <- lay$canvas[["x"]] + 4
  y0 <- lay$canvas[["y"]] + lay$canvas[["h"]] + 26
  out <- list()
  for (i in seq_along(items)) {
    it <- items[[i]]
    cx <- x0 + ((i - 1) %/% 4) * 460
    cy <- y0 + ((i - 1) %% 4) * 34
    kid <- el_id("legend", i)
    if (identical(it$kind, "box")) {
      st <- STATUS_STYLE[[it$status]]
      out <- c(out, list(base_el(
        kid, "rectangle", cx, cy, 34, 20, st$stroke, st$fill,
        stroke_style = if (isFALSE(it$built)) "dashed" else "solid",
        stroke_width = 2, roundness = 3L, updated = updated
      )))
    } else {
      inferred <- identical(it$style, "inferred")
      es <- if (inferred) list(stroke = INFERRED_STROKE, style = "dotted", width = 1)
            else EDGE_STYLE[[it$style]]
      ln <- base_el(kid, "arrow", cx, cy + 10, 34, 0, es$stroke, "transparent",
                    stroke_style = es$style, stroke_width = es$width,
                    roundness = 2L, updated = updated)
      ln$points <- list(c(0, 0), c(34, 0))
      ln$startBinding <- NULL
      ln$endBinding <- NULL
      ln$startArrowhead <- NULL
      ln$endArrowhead <- "arrow"
      ln$elbowed <- FALSE
      out <- c(out, list(ln))
    }
    lab <- text_el(el_id("legend", i, "text"), NULL, it$text,
                   cx + 46, cy, 400, 20, "#495057", 13, updated, align = "left")
    lab$containerId <- NULL
    out <- c(out, list(lab))
  }
  out
}

# --------------------------------------------------------------- validation ---

# Re-read the written file and check the invariants that make a scene openable:
# bindings resolve, container/text references are reciprocal, node boxes do not
# overlap. Band rectangles and the legend are excluded from the overlap test --
# bands contain nodes by design.
validate_scene <- function(path, lay) {
  s <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  els <- s$elements
  ids <- vapply(els, `[[`, character(1), "id")
  if (anyDuplicated(ids)) {
    stop("FAIL duplicate element ids: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "))
  }
  idx <- setNames(seq_along(els), ids)
  problems <- character(0)
  notes <- character(0)

  for (el in els) {
    # every binding target exists
    for (f in c("startBinding", "endBinding")) {
      tgt <- el[[f]]$elementId
      if (!is.null(tgt) && !(tgt %in% ids)) {
        problems <- c(problems, sprintf("%s: %s -> missing element '%s'", el$id, f, tgt))
      }
    }
    # a bound text element points back at its container
    for (be in el$boundElements %||% list()) {
      if (!(be$id %in% ids)) {
        problems <- c(problems, sprintf("%s: boundElements -> missing element '%s'",
                                        el$id, be$id))
        next
      }
      child <- els[[idx[[be$id]]]]
      if (identical(be$type, "text") && !identical(child$containerId, el$id)) {
        problems <- c(problems, sprintf(
          "%s: bound text '%s' has containerId '%s'", el$id, be$id,
          child$containerId %||% "<null>"))
      }
    }
    # a contained text element is listed by its container
    if (!is.null(el$containerId)) {
      if (!(el$containerId %in% ids)) {
        problems <- c(problems, sprintf("%s: containerId -> missing element '%s'",
                                        el$id, el$containerId))
      } else {
        parent <- els[[idx[[el$containerId]]]]
        listed <- vapply(parent$boundElements %||% list(),
                         function(b) identical(b$id, el$id), logical(1))
        if (!any(listed)) {
          problems <- c(problems, sprintf(
            "%s: containerId '%s' does not list it in boundElements",
            el$id, el$containerId))
        }
      }
    }
    # arrow points are relative, first point at the origin
    if (identical(el$type, "arrow")) {
      p1 <- unlist(el$points[[1]])
      if (!identical(as.numeric(p1), c(0, 0))) {
        problems <- c(problems, sprintf("%s: first arrow point is not [0,0]", el$id))
      }
    }
  }

  # node boxes must not overlap
  nd <- lay$nodes
  for (i in seq_len(nrow(nd) - 1L)) {
    for (j in seq(i + 1L, nrow(nd))) {
      a <- nd[i, ]; b <- nd[j, ]
      if (a$x < b$x + b$w - 1e-9 && b$x < a$x + a$w - 1e-9 &&
          a$y < b$y + b$h - 1e-9 && b$y < a$y + a$h - 1e-9) {
        problems <- c(problems, sprintf("node boxes overlap: %s and %s", a$id, b$id))
      }
    }
  }

  # Workflow kind only: a free caption should not sit on a box. Text widths for
  # a proportional face are estimated, so a marginal overlap is reported as a
  # note rather than failing the build; a caption more than half buried is an
  # error, because that is a layout fault rather than a metrics rounding.
  if (!is.null(lay$captions) && length(lay$captions)) {
    for (cp in lay$captions) {
      area <- max(cp$w * cp$h, 1e-9)
      for (i in seq_len(nrow(nd))) {
        ox <- min(cp$x + cp$w, nd$x[i] + nd$w[i]) - max(cp$x, nd$x[i])
        oy <- min(cp$y + cp$h, nd$y[i] + nd$h[i]) - max(cp$y, nd$y[i])
        if (ox <= 0 || oy <= 0) next
        frac <- (ox * oy) / area
        msg <- sprintf("caption '%s' covers %.0f%% of box %s",
                       cp$id, 100 * frac, nd$id[i])
        if (frac > 0.5) problems <- c(problems, msg) else notes <- c(notes, msg)
      }
    }
  }

  if (length(notes)) {
    cat(sprintf("    note: %s\n", notes), sep = "")
  }
  if (length(problems)) {
    stop("scene validation failed for ", path, ":\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }
  length(els)
}

# -------------------------------------------------------------------- driver --

build_one <- function(spec_path, out_dir) {
  spec <- yaml::read_yaml(spec_path)
  kind <- spec$kind %||% "architecture"
  need <- switch(kind,
    architecture = c("layers", "nodes", "edges"),
    workflow     = "stages",
    stop(spec_path, ": unknown kind '", kind, "'")
  )
  for (f in need) {
    if (is.null(spec[[f]])) stop(spec_path, ": spec has no `", f, "` block")
  }
  built <- switch(kind,
    architecture = build_scene(spec, spec_path),
    workflow     = build_workflow_scene(spec, spec_path)
  )
  name <- sub("\\.ya?ml$", "", basename(spec_path))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(out_dir, paste0(name, ".excalidraw"))
  jsonlite::write_json(built$scene, out, auto_unbox = TRUE, null = "null",
                       digits = NA, pretty = 2)
  n <- validate_scene(out, built$layout)
  counts <- switch(kind,
    architecture = c(length(spec$nodes), length(spec$edges)),
    workflow     = c(nrow(built$layout$nodes), length(built$layout$edges))
  )
  cat(sprintf("  %-28s %3d nodes, %3d edges, %3d elements -> %s\n",
              basename(spec_path), counts[[1]], counts[[2]], n, out))
  invisible(out)
}

main <- function(args) {
  here <- tryCatch({
    a <- commandArgs(trailingOnly = FALSE)
    dirname(normalizePath(sub("^--file=", "", a[grepl("^--file=", a)])[1]))
  }, error = function(e) getwd())
  specs <- if (length(args)) args else
    list.files(file.path(here, "specs"), pattern = "\\.ya?ml$", full.names = TRUE)
  if (!length(specs)) stop("no specs found under ", file.path(here, "specs"))
  cat("stage 1: spec -> excalidraw\n")
  for (s in specs) build_one(normalizePath(s), file.path(here, "out"))
  invisible(NULL)
}

# Run only when invoked as a script (Rscript build_diagrams.R ...), so that the
# file can be sourced to get at its functions -- the test below does that.
if (any(grepl("^--file=", commandArgs(trailingOnly = FALSE)))) {
  main(commandArgs(trailingOnly = TRUE))
}
