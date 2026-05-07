# Run this script to test the app locally before deploying

# 1. Install required packages (run once)
pkgs <- c("shiny", "bslib", "DT", "ggplot2", "MASS", "nnet", "readxl")
missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing) > 0) install.packages(missing)

# 2. Run the app locally
shiny::runApp("app.R", launch.browser = TRUE)

# ── When ready to deploy to Shinylive ────────────────────────────────────────
# install.packages("shinylive")   # install once
#
# From PARENT directory of this folder:
# shinylive::export(
#   appdir  = "biostats-tool",
#   destdir = "biostats-tool-site"
# )
#
# Then push biostats-tool-site/ to GitHub and enable Pages.
