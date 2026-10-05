// export.mjs -- render an Excalidraw scene to SVG. Stage 2 of 2; stage 1
// (build_diagrams.R) produces the scene from a spec.
//
//   node export.mjs out/foo.excalidraw out/foo.svg [--png out/foo.png]
//
// @excalidraw/excalidraw 0.18 ships only an ESM build whose module graph pulls
// in ~25 bare specifiers (react, react/jsx-runtime, clsx, roughjs, pako,
// lodash.throttle, ...), most of them CJS-only. An import map cannot resolve
// that, so the entry is bundled to an IIFE with esbuild and served, with the
// package's own font directory, from a throwaway localhost server -- module
// loading and font loading both need an http origin, and file:// gives neither.
//
// exportToSvg needs a DOM (it measures text and drives roughjs), hence
// Playwright rather than a headless shim. Nothing is fetched from the network:
// every off-origin request is aborted and then reported as a failure.

import { createServer } from "node:http";
import { readFile, writeFile, mkdir, mkdtemp, rm } from "node:fs/promises";
import { existsSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import * as esbuild from "esbuild";
import { chromium } from "playwright";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const PKG = path.join(HERE, "node_modules", "@excalidraw", "excalidraw");
const DIST = path.join(PKG, "dist", "prod");

const MIME = {
  ".js": "text/javascript; charset=utf-8",
  ".mjs": "text/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8",
  ".html": "text/html; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".woff2": "font/woff2",
  ".woff": "font/woff",
  ".ttf": "font/ttf",
  ".otf": "font/otf",
  ".png": "image/png",
  ".svg": "image/svg+xml"
};

function parseArgs(argv) {
  const positional = [];
  let png = null;
  for (let i = 0; i < argv.length; i++) {
    if (argv[i] === "--png") png = argv[++i];
    else positional.push(argv[i]);
  }
  const [input, output] = positional;
  if (!input || !output) {
    throw new Error(
      "usage: node export.mjs <scene.excalidraw> <out.svg> [--png <out.png>]"
    );
  }
  return { input, output, png };
}

// Confirm the dist layout rather than assuming it: the entry point moved
// between 0.17 and 0.18, so a missing file should say so plainly.
function resolveEntry() {
  const pkgJson = path.join(PKG, "package.json");
  if (!existsSync(pkgJson)) {
    throw new Error(`@excalidraw/excalidraw is not installed under ${PKG}`);
  }
  const entry = path.join(DIST, "index.js");
  if (!existsSync(entry)) {
    throw new Error(
      `expected the excalidraw browser bundle at ${entry}; ` +
        `dist layout has changed, re-check the package's "exports" field`
    );
  }
  if (!existsSync(path.join(DIST, "fonts"))) {
    throw new Error(`expected a local font directory at ${path.join(DIST, "fonts")}`);
  }
  return entry;
}

// Bundle `import { exportToSvg } from "@excalidraw/excalidraw"` to a single
// IIFE that hangs the API off window. react/react-dom resolve from the local
// node_modules; nothing is fetched.
async function bundleExcalidraw() {
  const stub = `
    import * as Excalidraw from "@excalidraw/excalidraw";
    window.__excalidraw = Excalidraw;
  `;
  const result = await esbuild.build({
    stdin: { contents: stub, resolveDir: HERE, loader: "js" },
    bundle: true,
    format: "iife",
    platform: "browser",
    target: "chrome120",
    write: false,
    minify: true,
    logLevel: "silent",
    define: { "process.env.NODE_ENV": '"production"', global: "window" },
    loader: { ".css": "empty", ".woff2": "dataurl", ".woff": "dataurl", ".ttf": "dataurl" }
  });
  return result.outputFiles[0].text;
}

function pageHtml() {
  return `<!doctype html>
<html><head><meta charset="utf-8"><title>export</title>
<style>html,body{margin:0;padding:0;background:#fff}</style></head>
<body><div id="root"></div>
<script>window.EXCALIDRAW_ASSET_PATH = "/excalidraw/dist/prod/";</script>
<script src="/bundle.js"></script>
</body></html>`;
}

// Serves three things and nothing else: the page, the bundle, and the
// excalidraw package directory (for fonts).
async function startServer(bundle) {
  const served = [];
  const missed = [];
  const server = createServer(async (req, res) => {
    const url = new URL(req.url, "http://localhost");
    served.push(url.pathname);
    const send = (code, body, type) => {
      if (code === 404) missed.push(url.pathname);
      res.writeHead(code, { "Content-Type": type || "text/plain" });
      res.end(body);
    };
    try {
      if (url.pathname === "/" || url.pathname === "/index.html") {
        return send(200, pageHtml(), MIME[".html"]);
      }
      if (url.pathname === "/bundle.js") {
        return send(200, bundle, MIME[".js"]);
      }
      if (url.pathname.startsWith("/excalidraw/")) {
        const rel = decodeURIComponent(url.pathname.slice("/excalidraw/".length));
        const file = path.resolve(PKG, rel);
        // never serve outside the package directory
        if (!file.startsWith(path.resolve(PKG))) return send(403, "forbidden");
        if (!existsSync(file)) return send(404, "not found");
        return send(200, await readFile(file), MIME[path.extname(file)] || "application/octet-stream");
      }
      return send(404, "not found");
    } catch (err) {
      return send(500, String(err));
    }
  });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  const { port } = server.address();
  return { server, port, served, missed };
}

async function main() {
  const { input, output, png } = parseArgs(process.argv.slice(2));
  resolveEntry();

  const scene = JSON.parse(await readFile(input, "utf8"));
  if (!Array.isArray(scene.elements) || scene.elements.length === 0) {
    throw new Error(`${input} has no elements`);
  }

  const bundle = await bundleExcalidraw();
  const { server, port, served, missed } = await startServer(bundle);
  const origin = `http://127.0.0.1:${port}`;

  const browser = await chromium.launch({ headless: true });
  const consoleErrors = [];
  const offOrigin = [];
  try {
    const page = await browser.newPage({ viewport: { width: 1600, height: 1200 } });
    page.on("console", (m) => {
      if (m.type() === "error" || m.type() === "warning") {
        consoleErrors.push(`[${m.type()}] ${m.text()}`);
      }
    });
    page.on("pageerror", (e) => consoleErrors.push(`[pageerror] ${e.message}`));
    // Anything not served by us is a bug (a CDN font, say); abort and report.
    await page.route("**", (route) => {
      const u = route.request().url();
      if (u.startsWith(origin)) return route.continue();
      offOrigin.push(u);
      return route.abort();
    });

    await page.goto(`${origin}/`, { waitUntil: "load" });
    await page.waitForFunction(
      () => typeof window.__excalidraw?.exportToSvg === "function",
      null,
      { timeout: 30000 }
    );

    // Excalidraw's restore() drops elements it considers malformed, and does so
    // silently -- the same pass exportToSvg runs internally. Count before and
    // after, so a bad element is a build failure rather than a gap in the SVG.
    const kept = await page.evaluate((sceneJson) => {
      const s = JSON.parse(sceneJson);
      const r = window.__excalidraw.restore(
        { elements: s.elements, appState: s.appState || {}, files: s.files || {} },
        null,
        null
      );
      // Also report any element whose geometry restore() changed: a container
      // regrown to fit its text moves every box below it and sends bound
      // arrows off on hooks, which is invisible in the JSON but obvious in the
      // SVG. Generated scenes must survive restore unchanged.
      const before = new Map(s.elements.map((e) => [e.id, e]));
      const moved = [];
      for (const e of r.elements) {
        const b = before.get(e.id);
        if (!b) continue;
        for (const k of ["x", "y", "width", "height"]) {
          if (Math.abs((b[k] ?? 0) - (e[k] ?? 0)) > 0.5) {
            moved.push(`${e.id}.${k}: ${Number(b[k]).toFixed(1)} -> ${Number(e[k]).toFixed(1)}`);
            break;
          }
        }
      }
      return { kept: r.elements.length, moved };
    }, JSON.stringify(scene));
    if (kept.kept !== scene.elements.length) {
      throw new Error(
        `excalidraw dropped ${scene.elements.length - kept.kept} of ` +
          `${scene.elements.length} elements as malformed during restore`
      );
    }
    if (kept.moved.length) {
      console.log(
        `    note: restore moved ${kept.moved.length} element(s), e.g.\n` +
          kept.moved.slice(0, 8).map((m) => `      ${m}`).join("\n")
      );
    }

    const svg = await page.evaluate(async (sceneJson) => {
      const s = JSON.parse(sceneJson);
      const node = await window.__excalidraw.exportToSvg({
        elements: s.elements,
        appState: {
          ...(s.appState || {}),
          exportBackground: true,
          exportWithDarkMode: false,
          exportEmbedScene: false
        },
        files: s.files || {}
      });
      if (!node) return null;
      return new XMLSerializer().serializeToString(node);
    }, JSON.stringify(scene));

    if (!svg) {
      throw new Error(
        "exportToSvg returned nothing" +
          (consoleErrors.length ? `\npage console:\n  ${consoleErrors.join("\n  ")}` : "")
      );
    }

    await mkdir(path.dirname(path.resolve(output)), { recursive: true });
    await writeFile(output, svg, "utf8");

    if (png) {
      // Rasterise the SVG we just wrote, as a cheap visual check.
      await page.setContent(
        `<html><body style="margin:0;background:#fff">${svg}</body></html>`
      );
      const el = await page.$("svg");
      await mkdir(path.dirname(path.resolve(png)), { recursive: true });
      await el.screenshot({ path: png, scale: "css" });
    }

    if (offOrigin.length) {
      throw new Error(
        `export tried to fetch ${offOrigin.length} off-origin resource(s), ` +
          `which were blocked:\n  ${[...new Set(offOrigin)].join("\n  ")}`
      );
    }

    const kb = (Buffer.byteLength(svg, "utf8") / 1024).toFixed(0);
    console.log(
      `  ${path.basename(input)} -> ${output} (${kb} kB${png ? `, ${png}` : ""})`
    );
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((err) => {
  console.error(`export failed: ${err.message}`);
  process.exit(1);
});
