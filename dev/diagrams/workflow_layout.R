# workflow_layout.R -- the `kind: workflow` layout engine, sourced by
# build_diagrams.R. Kept in its own file so that the architecture engine stays
# readable; it defines no elements at load time and draws no random numbers at
# file scope, so sourcing it cannot move the architecture kind's output.
#
# This reproduces the idiom of dev/diagrams/optimisation_pipeline.excalidraw: a
# vertical spine that forks into two lanes and re-merges, purple decision
# diamonds, clean lines (roughness 0) in Excalifont, and -- unlike the
# architecture kind -- no text bound to any arrow. Every caption is a free grey
# text placed beside its connector.
#
# The constants below are measured from that file rather than invented. Two are
# easy to get wrong: the intra-lane invariant is the 15px *gap*, not a fixed
# pitch (box heights vary from 37 to 50 with their text), and a lane-to-spine
# merge is a 4-point stub-across-drop, while 3-point routes are reserved for a
# decision diamond's side exit.

WF_ROLE <- list(
  api       = list(stroke = "#1e88e5", fill = "#e3f2fd", width = 2,   style = "solid",  shape = "rectangle"),
  subsystem = list(stroke = "#fb8c00", fill = "#fff3e0", width = 2.5, style = "dashed", shape = "rectangle"),
  helper    = list(stroke = "#43a047", fill = "#e8f5e9", width = 2,   style = "solid",  shape = "rectangle"),
  decision  = list(stroke = "#8e24aa", fill = "#f3e5f5", width = 2,   style = "solid",  shape = "diamond"),
  fallback  = list(stroke = "#e65100", fill = "#fff3e0", width = 2,   style = "solid",  shape = "rectangle"),
  terminal  = list(stroke = "#e53935", fill = "#fce4ec", width = 2,   style = "solid",  shape = "rectangle")
)
WF_ORIGIN    <- c(x = 110, y = 131)
WF_SPINE_X   <- 500
WF_LANE_X    <- c(spine = 500, left = 280, right = 760, gutter = 1010)
WF_CHANNEL_X <- c(left = 190, right = 1040)
WF_W <- list(spine = 360, header = 340, step = 280, result = 320,
             gutter = 254, decision = 280)
WF_GAP_LANE  <- 15   # header -> step -> ... -> result; the lane invariant
WF_GAP_SPINE <- 50   # between consecutive spine blocks
WF_GAP_STAGE <- 55   # across a lane boundary
WF_STUB      <- 20   # short vertical stub before a sideways run
WF_MIN_DROP  <- 24
WF_PAD_X     <- 10
WF_PAD_Y     <- 8
WF_FS <- list(spine = 16, spine_multi = 13, header = 16, step = 11,
              result = 12, decision = 14, caption = 12)
WF_MIN_H <- list(spine = 42, header = 37, step = 37, result = 42, decision = 80)
WF_CHAR_F <- 0.52    # Excalifont is proportional; measured 0.45 to 0.55
WF_FAMILY <- 5L      # Excalifont
WF_INK    <- "#1e1e1e"
WF_ANNOT  <- "#495057"
WF_DIM    <- 55L     # opacity for a stage specified but not yet coded
WF_CAP_DX <- 8
WF_CAP_DY <- 4

# Lane centres sit on the 20px lattice; a box's left edge is centre - width/2
# and must NOT be re-snapped. Snapping it moved a 340-wide header off its own
# lane centre, which made every header-to-step edge read as a column change.
wf_snap   <- function(x) x
wf_text_w <- function(lines, size) max(nchar(lines)) * WF_CHAR_F * size
wf_text_h <- function(lines, size) length(lines) * size * 1.25

wf_role <- function(role) {
  r <- WF_ROLE[[role]]
  if (is.null(r)) {
    stop("unknown role '", role, "'; expected one of ",
         paste(names(WF_ROLE), collapse = ", "), call. = FALSE)
  }
  r
}

wf_box_h <- function(lines, size, kind) {
  th <- wf_text_h(lines, size)
  if (identical(kind, "decision")) {
    max(WF_MIN_H$decision, 2.6 * th)
  } else {
    max(WF_MIN_H[[kind]], th + 2 * WF_PAD_Y)
  }
}

# One box: title, then body lines, then the status line. The status is the last
# line of the same text element, so it inherits the box colour and size. That is
# deliberate: it leaves the five role colours doing exactly one job.
wf_norm_item <- function(it, place, kind, default_w) {
  if (is.null(it$id)) stop("a ", kind, " box has no id", call. = FALSE)
  lines <- c(it$title %||% it$id,
             unlist(it$body %||% character(0)),
             it$status %||% character(0))
  font <- it$font %||%
    (if (identical(kind, "spine") && length(lines) > 1L) WF_FS$spine_multi
     else WF_FS[[kind]])
  role <- it$role %||% "helper"
  list(id = it$id, place = place, kind = kind, role = role,
       shape = wf_role(role)$shape, lines = lines, font = font,
       width = it$width %||% default_w, height = it$height,
       built = !identical(it$built, FALSE))
}

# ------------------------------------------------------------------- layout ----
# One y-cursor walks the flat stage list. A lane is a compound stage: a wider
# header, narrow steps at a constant gap, then a wider result -- the "hat and
# foot" that brackets a lane without any frame or background rectangle.
# `align:` lets a stage share another stage's top y in a different column, and
# cursor <- max(...) afterwards folds that branch back into the trunk.
wf_layout <- function(spec) {
  implicit <- !identical(spec$implicit_edges, FALSE)
  items <- list()
  imp <- list()
  ytop <- list()
  cursor <- WF_ORIGIN[["y"]]

  for (st in spec$stages) {
    if (is.null(st$id)) stop("every stage needs an id", call. = FALSE)
    place <- st$place %||% "spine"
    if (!(place %in% names(WF_LANE_X))) {
      stop("stage '", st$id, "': unknown place '", place, "'", call. = FALSE)
    }
    cx <- WF_LANE_X[[place]]
    top <- cursor
    if (!is.null(st$align)) {
      if (is.null(ytop[[st$align]])) {
        stop("stage '", st$id, "': align target '", st$align,
             "' is not a stage placed before it", call. = FALSE)
      }
      top <- ytop[[st$align]]
    }

    if (!is.null(st$steps)) {
      chain <- character(0)
      yy <- top
      parts <- c(
        list(list(it = st$header, kind = "header")),
        lapply(st$steps, function(s) list(it = s, kind = "step")),
        list(list(it = st$result, kind = "result"))
      )
      for (p in parts) {
        if (is.null(p$it)) next
        bx <- wf_norm_item(p$it, place, p$kind, WF_W[[p$kind]])
        bx$h <- bx$height %||% wf_box_h(bx$lines, bx$font, p$kind)
        bx$x <- wf_snap(cx - bx$width / 2)
        bx$y <- yy
        items <- c(items, list(bx))
        chain <- c(chain, bx$id)
        ytop[[bx$id]] <- yy
        yy <- yy + bx$h + WF_GAP_LANE
      }
      bottom <- yy - WF_GAP_LANE
      if (implicit && length(chain) > 1L) {
        for (k in seq_len(length(chain) - 1L)) {
          imp <- c(imp, list(list(from = chain[[k]], to = chain[[k + 1L]])))
        }
      }
      gap_after <- st$gap %||% WF_GAP_STAGE
    } else {
      kind <- if (identical(st$role, "decision")) "decision" else "spine"
      dw <- if (identical(kind, "decision")) {
        WF_W$decision
      } else if (identical(place, "spine")) {
        WF_W$spine
      } else if (identical(place, "gutter")) {
        WF_W$gutter
      } else {
        WF_W$step
      }
      bx <- wf_norm_item(st, place, kind, dw)
      bx$h <- bx$height %||% wf_box_h(bx$lines, bx$font, kind)
      bx$x <- wf_snap(cx - bx$width / 2)
      bx$y <- top
      items <- c(items, list(bx))
      bottom <- top + bx$h
      gap_after <- st$gap %||%
        (if (identical(place, "spine")) WF_GAP_SPINE else WF_GAP_STAGE)
    }
    ytop[[st$id]] <- top
    cursor <- max(cursor, bottom + gap_after)
  }

  nd <- data.frame(
    id    = vapply(items, `[[`, character(1), "id"),
    place = vapply(items, `[[`, character(1), "place"),
    kind  = vapply(items, `[[`, character(1), "kind"),
    role  = vapply(items, `[[`, character(1), "role"),
    shape = vapply(items, `[[`, character(1), "shape"),
    x     = vapply(items, `[[`, numeric(1), "x"),
    y     = vapply(items, `[[`, numeric(1), "y"),
    w     = vapply(items, `[[`, numeric(1), "width"),
    h     = vapply(items, `[[`, numeric(1), "h"),
    font  = vapply(items, `[[`, numeric(1), "font"),
    built = vapply(items, `[[`, logical(1), "built"),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(nd$id)) {
    stop("duplicate box ids: ",
         paste(unique(nd$id[duplicated(nd$id)]), collapse = ", "), call. = FALSE)
  }

  # Only intra-lane links are implicit. Everything that changes column is
  # declared, so that an `align`ed branch can never acquire an inferred edge.
  # An explicit edge with the same endpoints wins over the implicit one.
  explicit <- spec$edges %||% list()
  ekey <- function(e) paste(e$from, e$to, sep = "\r")
  seen <- vapply(explicit, ekey, character(1))
  imp <- Filter(function(e) !(ekey(e) %in% seen), imp)
  edges <- c(explicit, imp)
  for (e in edges) {
    for (endp in c(e$from, e$to)) {
      if (!(endp %in% nd$id)) {
        stop("edge endpoint is not a box: '", endp, "'", call. = FALSE)
      }
    }
  }

  list(nodes = nd, items = items, edges = edges, kind = "workflow",
       canvas = c(x = min(nd$x), y = WF_ORIGIN[["y"]],
                  w = max(nd$x + nd$w) - min(nd$x),
                  h = max(nd$y + nd$h) - WF_ORIGIN[["y"]]))
}

# --------------------------------------------------------- connector router ----
# Points are relative to the arrow's own x/y with the first point at [0,0].
# Four cases, matching the distribution observed in the source diagram.
wf_route <- function(e, a, b) {
  ax <- a$x + a$w / 2
  bx <- b$x + b$w / 2
  a_bot <- a$y + a$h
  b_top <- b$y
  drop <- b_top - a_bot
  same <- abs(ax - bx) < 1

  forced <- !is.null(e$route) || !is.null(e$via_x)

  if (!forced && same && drop > 0) {
    # plain vertical, same column, adjacent
    list(x = ax, y = a_bot, points = list(c(0, 0), c(0, drop)),
         anchor = list(kind = "v", x = ax, y = a_bot + drop / 2, side = 1))
  } else if (!forced && identical(a$shape, "diamond") && !same && drop >= WF_MIN_DROP) {
    # a decision's side exit: leave the side vertex, run across, then drop.
    # The second point stops 10px shy of the third so the bend reads as a curve
    # under roundness 2, as it does in the source.
    dir <- sign(bx - ax)
    sx <- if (dir > 0) a$x + a$w else a$x
    sy <- a$y + a$h / 2
    run <- bx - sx
    list(x = sx, y = sy,
         points = list(c(0, 0), c(run - 10 * dir, 0), c(run, b_top - sy)),
         anchor = list(kind = "h", x = sx + run / 2, y = sy, side = dir))
  } else if (!forced && !same && drop >= WF_STUB + WF_MIN_DROP) {
    # stub down, across, then drop: the lane/spine merge and split
    d <- bx - ax
    list(x = ax, y = a_bot,
         points = list(c(0, 0), c(0, WF_STUB), c(d, WF_STUB), c(d, drop)),
         anchor = list(kind = "h", x = ax + d / 2, y = a_bot + WF_STUB,
                       side = sign(d)))
  } else if (!forced && !same && drop > 0) {
    # a small column change over a short drop: jog straight across rather than
    # sweeping out to a channel, which would be absurd for a 15px gap
    d <- bx - ax
    list(x = ax, y = a_bot,
         points = list(c(0, 0), c(d, drop / 2), c(d, drop)),
         anchor = list(kind = "v", x = bx, y = a_bot + drop / 2, side = 1))
  } else {
    # long haul: sweep out to a side channel and back in
    via <- e$via_x %||% switch(
      e$route %||% "auto",
      left = WF_CHANNEL_X[["left"]],
      right = WF_CHANNEL_X[["right"]],
      gutter = WF_LANE_X[["gutter"]],
      if (max(ax, bx) <= WF_SPINE_X) WF_CHANNEL_X[["left"]]
      else WF_CHANNEL_X[["right"]]
    )
    d <- bx - ax
    v <- via - ax
    list(x = ax, y = a_bot,
         points = list(c(0, 0), c(v, WF_STUB), c(v, drop - WF_STUB), c(d, drop)),
         anchor = list(kind = "v", x = via, y = a_bot + drop / 2,
                       side = if (via > WF_SPINE_X) 1 else -1))
  }
}

wf_caption_geom <- function(anchor, lines, size) {
  w <- wf_text_w(lines, size) + 6
  h <- wf_text_h(lines, size)
  if (identical(anchor$kind, "v")) {
    x <- if (anchor$side < 0) anchor$x - WF_CAP_DX - w else anchor$x + WF_CAP_DX
    y <- anchor$y - h / 2
  } else {
    x <- anchor$x - w / 2
    y <- anchor$y + WF_CAP_DY
  }
  list(x = x, y = y, w = w, h = h)
}

wf_rect_hit <- function(a, b, eps = 1e-9) {
  a$x < b$x + b$w - eps && b$x < a$x + a$w - eps &&
    a$y < b$y + b$h - eps && b$y < a$y + a$h - eps
}

# Nudge a caption along its connector's free axis until it clears every box and
# every caption already placed. Same idiom as the architecture kind's badges.
wf_place_caption <- function(geom, anchor, boxes, placed) {
  for (k in c(0, 0.6, -0.6, 1.2, -1.2, 1.9, -1.9)) {
    g <- geom
    g$y <- geom$y + k * (geom$h + 5)
    if (!identical(anchor$kind, "v") && k != 0) {
      # a run caption may also shift a little along the run
      g$x <- geom$x + sign(k) * 12
    }
    clash <- any(vapply(boxes, function(b) wf_rect_hit(g, b), logical(1))) ||
      any(vapply(placed, function(p) wf_rect_hit(g, p), logical(1)))
    if (!clash) return(g)
  }
  geom
}

# ---------------------------------------------------------------- the scene ----
# Z-order, as for the architecture kind: arrows, then boxes and their bound text,
# then captions and annotations, then the legend. Arrows under the boxes means a
# long haul crossing the picture never obscures a label.
build_workflow_scene <- function(spec, spec_path) {
  set.seed(sum(utf8ToInt(basename(spec_path))))
  updated <- round(as.numeric(file.mtime(spec_path)) * 1000)

  lay <- wf_layout(spec)
  nd <- lay$nodes
  rownames(nd) <- nd$id
  by_id <- setNames(lay$items, vapply(lay$items, `[[`, character(1), "id"))
  elements <- list()
  captions <- list()

  box_rects <- lapply(seq_len(nrow(nd)), function(i) {
    list(x = nd$x[i], y = nd$y[i], w = nd$w[i], h = nd$h[i])
  })
  annots <- spec$annotations %||% list()
  obstacles <- c(box_rects, lapply(annots, function(an) {
    at <- unlist(an$at)
    lines <- unlist(an$text)
    list(x = at[[1]], y = at[[2]],
         w = wf_text_w(lines, WF_FS$caption) + 6,
         h = wf_text_h(lines, WF_FS$caption))
  }))

  # 1. arrows, with their captions collected for later placement
  for (e in lay$edges) {
    a <- nd[e$from, ]
    b <- nd[e$to, ]
    a$shape <- nd[e$from, "shape"]
    rt <- wf_route(e, a, b)
    aid <- el_id("a", e$from, e$to)
    dim_edge <- !nd[e$from, "built"] && !nd[e$to, "built"]
    dxs <- vapply(rt$points, function(p) p[[1]], numeric(1))
    dys <- vapply(rt$points, function(p) p[[2]], numeric(1))
    arr <- base_el(aid, "arrow", rt$x, rt$y,
                   diff(range(dxs)), diff(range(dys)),
                   WF_ANNOT, "transparent", stroke_width = 2,
                   roundness = 2L, updated = updated, roughness = 0L,
                   opacity = if (dim_edge) WF_DIM else 100L)
    arr$points <- rt$points
    arr$lastCommittedPoint <- NULL
    arr$startBinding <- list(elementId = el_id("n", e$from), focus = 0, gap = 1)
    arr$endBinding <- list(elementId = el_id("n", e$to), focus = 0, gap = 1)
    arr$startArrowhead <- NULL
    arr$endArrowhead <- "arrow"
    arr$elbowed <- FALSE
    elements <- c(elements, list(arr))

    txt <- c(if (!is.null(e$branch)) e$branch else character(0),
             unlist(e$caption %||% character(0)))
    if (length(txt)) {
      captions[[length(captions) + 1L]] <- list(
        id = el_id("c", e$from, e$to), lines = txt, anchor = rt$anchor
      )
    }
  }

  # 2. boxes and their bound text
  for (i in seq_len(nrow(nd))) {
    it <- by_id[[nd$id[i]]]
    st <- wf_role(it$role)
    op <- if (it$built) 100L else WF_DIM
    nid <- el_id("n", it$id)
    tid <- el_id("t", it$id)
    shp <- base_el(nid, st$shape, nd$x[i], nd$y[i], nd$w[i], nd$h[i],
                   st$stroke, st$fill, stroke_style = st$style,
                   stroke_width = st$width,
                   roundness = if (identical(st$shape, "diamond")) 2L else 3L,
                   updated = updated, roughness = 0L, opacity = op)
    shp$boundElements <- list(list(id = tid, type = "text"))
    elements <- c(elements, list(shp))
    th <- wf_text_h(it$lines, it$font)
    elements <- c(elements, list(text_el(
      tid, nid, it$lines,
      nd$x[i] + WF_PAD_X, nd$y[i] + (nd$h[i] - th) / 2,
      nd$w[i] - 2 * WF_PAD_X, th,
      WF_INK, it$font, updated,
      family = WF_FAMILY, roughness = 0L, auto_resize = FALSE, opacity = op
    )))
  }

  # 3. captions, de-collided against the boxes and each other
  placed <- list()
  cap_boxes <- list()
  for (cp in captions) {
    g <- wf_place_caption(
      wf_caption_geom(cp$anchor, cp$lines, WF_FS$caption),
      cp$anchor, obstacles, placed
    )
    placed[[length(placed) + 1L]] <- g
    cap_boxes[[length(cap_boxes) + 1L]] <- c(g, list(id = cp$id))
    el <- text_el(cp$id, NULL, cp$lines, g$x, g$y, g$w, g$h,
                  WF_ANNOT, WF_FS$caption, updated, align = "left",
                  family = WF_FAMILY, roughness = 0L, valign = "top",
                  auto_resize = FALSE)
    el$containerId <- NULL
    elements <- c(elements, list(el))
  }

  # 4. free annotations, placed exactly where the spec asks.
  # `at` is a positional [x, y] pair, not a mapping: YAML 1.1 reads a bare `y:`
  # key as the boolean true, so `at: {x: 700, y: 140}` silently loses the y.
  for (an in annots) {
    lines <- unlist(an$text)
    at <- unlist(an$at)
    if (length(at) != 2L || !is.numeric(at)) {
      stop("annotation '", an$id %||% "?",
           "': `at` must be a pair of numbers, [x, y]", call. = FALSE)
    }
    el <- text_el(el_id("an", an$id %||% paste(at, collapse = "-")), NULL,
                  lines, at[[1]], at[[2]],
                  wf_text_w(lines, WF_FS$caption) + 6,
                  wf_text_h(lines, WF_FS$caption),
                  WF_ANNOT, WF_FS$caption, updated, align = "left",
                  family = WF_FAMILY, roughness = 0L, valign = "top",
                  auto_resize = FALSE)
    el$containerId <- NULL
    elements <- c(elements, list(el))
  }

  # 5. legend
  if (!identical(spec$legend, FALSE)) {
    elements <- c(elements, wf_legend_elements(lay, updated))
  }

  # Bind arrows back onto the shapes they touch.
  idx <- setNames(seq_along(elements), vapply(elements, `[[`, character(1), "id"))
  for (e in lay$edges) {
    aid <- el_id("a", e$from, e$to)
    for (endp in c(e$from, e$to)) {
      k <- idx[[el_id("n", endp)]]
      elements[[k]]$boundElements <- c(elements[[k]]$boundElements,
                                       list(list(id = aid, type = "arrow")))
    }
  }

  lay$captions <- cap_boxes
  list(
    scene = list(
      type = "excalidraw", version = 2L,
      source = "multigrain/dev/diagrams",
      elements = elements,
      appState = list(gridSize = 20L, viewBackgroundColor = "#ffffff"),
      files = setNames(list(), character(0))
    ),
    layout = lay
  )
}

wf_legend_elements <- function(lay, updated) {
  rows <- list(
    list(role = "api",       text = "public API, the functions a user calls"),
    list(role = "subsystem", text = "internal subsystem, or its result object"),
    list(role = "helper",    text = "leaf helper called inside a subsystem"),
    list(role = "decision",  text = "branch"),
    list(role = "terminal",  text = "returned object"),
    list(role = "api",       text = "faded: specified but not yet coded (P3 to P5)",
         built = FALSE)
  )
  x0 <- lay$canvas[["x"]]
  y0 <- lay$canvas[["y"]] + lay$canvas[["h"]] + 46
  out <- list(text_el(el_id("wf", "legendhead"), NULL, "Colour marks the role.",
                      x0, y0 - 26, 420, 18, WF_ANNOT, 13, updated,
                      align = "left", family = WF_FAMILY, roughness = 0L,
                      valign = "top", auto_resize = FALSE))
  out[[1]]$containerId <- NULL
  for (i in seq_along(rows)) {
    r <- rows[[i]]
    st <- wf_role(r$role)
    op <- if (identical(r$built, FALSE)) WF_DIM else 100L
    cx <- x0 + ((i - 1L) %/% 3L) * 470
    cy <- y0 + ((i - 1L) %% 3L) * 32
    out <- c(out, list(base_el(
      el_id("wf", "legend", i), st$shape, cx, cy, 38, 22,
      st$stroke, st$fill, stroke_style = st$style, stroke_width = st$width,
      roundness = if (identical(st$shape, "diamond")) 2L else 3L,
      updated = updated, roughness = 0L, opacity = op
    )))
    lab <- text_el(el_id("wf", "legend", i, "t"), NULL, r$text,
                   cx + 50, cy + 3, 400, 18, WF_ANNOT, WF_FS$caption, updated,
                   align = "left", family = WF_FAMILY, roughness = 0L,
                   valign = "top", auto_resize = FALSE)
    lab$containerId <- NULL
    out <- c(out, list(lab))
  }
  out
}
