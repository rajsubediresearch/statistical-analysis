# BioStat Explorer — Shinylive Deployment Guide

## What this is
A browser-based biostatistics tool built with R Shiny + Shinylive.
Runs entirely client-side — no server, no cost, deploys free on GitHub Pages.

## Models included
- Simple & Multiple Linear Regression
- Simple & Multiple Logistic Regression
- Ordinal Logistic Regression (MASS::polr)
- Multinomial Logistic Regression (nnet::multinom)
- Independent t-test
- One-way ANOVA

## Features
- Two built-in sample datasets (Epidemiology & Clinical Trial)
- Upload your own CSV/TXT
- Auto-generated plain-English interpretations
- Coefficient tables with CIs
- Visualization tab (histogram, scatter, boxplot, density, bar)
- Causal inference disclaimer

---

## Deployment Steps

### Step 1 — Install R packages locally

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "MASS", "nnet", "shinylive"))
```

### Step 2 — Export the app to Shinylive format

```r
library(shinylive)

# Run from the parent directory of biostats-tool/
shinylive::export(
  appdir = "biostats-tool",
  destdir = "biostats-tool-site"
)
```

This creates a `biostats-tool-site/` folder with static HTML/JS/WASM files.

### Step 3 — Deploy to GitHub Pages

**Option A — Separate repo (recommended)**
1. Create a new repo: `biostat-explorer` on GitHub
2. Push the contents of `biostats-tool-site/` to the `main` branch
3. Go to repo Settings → Pages → Deploy from branch → main / root
4. Your tool will be live at:
   `https://rajsubediresearch.github.io/biostat-explorer/`
   And once your custom domain is set: `https://rajsubediresearch.com/biostat-explorer/`

**Option B — Subfolder of existing site**
1. Copy `biostats-tool-site/` contents into a `biostat-explorer/` folder
   in your main `rajsubediresearch.github.io` repo
2. Push to main
3. Live at: `https://rajsubediresearch.com/biostat-explorer/`

### Step 4 — Add to your main site
In your `rajsubediresearch.github.io` index, add a card linking to:
`/biostat-explorer/`

---

## Local testing (before deploying)

```r
# Test the Shiny app locally first
shiny::runApp("biostats-tool")

# Then test the Shinylive export locally
shinylive::export("biostats-tool", "biostats-tool-site")
httpuv::runStaticServer("biostats-tool-site")
```

---

## Notes on package support in WebR/Shinylive
- `shiny`, `bslib`, `DT`, `ggplot2`, `MASS`, `nnet` — all supported ✅
- `readxl` — NOT available in WebR; Excel uploads require local Shiny
  (CSV fallback is included in the app)
- First load takes ~10-20 seconds as WebR downloads R packages to browser

---

## Adding a DOI (Zenodo)
Once deployed, archive a release on Zenodo for a citable DOI badge,
just like your DAG Builder and Age Standardization tools.
