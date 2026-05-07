# BioStat Explorer

**An interactive, browser-based biostatistics tool for education and research.**  
Built with R Shiny · Deployable via Shinylive on GitHub Pages · No server required · Free

> Developed by [Raj Subedi](https://rajsubediresearch.com) · PhD Student, Epidemiology · Georgia State University

---

## Overview

BioStat Explorer lets users upload their own dataset (CSV, TXT, or Excel), select variables, choose a statistical model, and instantly receive results — including a coefficient table, confidence intervals, and **auto-generated plain-English interpretations** of the findings.

Designed for epidemiology students, public health researchers, and anyone learning applied biostatistics without needing to write code.

---

## Statistical Methods

| Method | Function |
|--------|----------|
| Simple Linear Regression | `lm()` |
| Multiple Linear Regression | `lm()` |
| Simple Logistic Regression | `glm(..., binomial)` |
| Multiple Logistic Regression | `glm(..., binomial)` |
| Ordinal Logistic Regression | `MASS::polr()` |
| Multinomial Logistic Regression | `nnet::multinom()` |
| Independent t-test | `t.test()` |
| One-way ANOVA | `aov()` |

---

## Features

-  **Flexible data input** — upload CSV, TXT, or Excel files; or use built-in sample datasets
-  **Auto-interpretation** — results are narrated in plain English with variable names embedded (e.g., *"Smoking is associated with 68% higher odds of disease, adjusting for age and BMI"*)
-  **Coefficient tables** — estimates, standard errors, p-values, and 95% CIs
-  **Analysis plots** — contextually appropriate visualizations per model type (scatter + regression line, violin plots, stacked bar charts, etc.)
-  **Explore tab** — standalone histogram, density, scatter, boxplot, and bar chart builder
-  **Causal inference disclaimer** — reminds users that statistical associations require careful interpretation
-  **No installation needed** — deployable as a static site via Shinylive + GitHub Pages

---

## Built-in Sample Datasets

Two simulated datasets are included for immediate exploration:

- **Epidemiology Dataset (n=300)** — age, sex, BMI, smoking, vaccination, systolic BP, disease outcome, severity, transport mode
- **Clinical Trial Dataset (n=250)** — treatment group, age, weight, cholesterol, BP change, recovery, outcome category

---

## Run Locally

```r
# Install dependencies (run once)
install.packages(c("shiny", "bslib", "DT", "ggplot2", "MASS", "nnet", "readxl"))

# Launch the app
shiny::runApp("app.R")
```

---

## Deploy to GitHub Pages (Shinylive)

```r
install.packages("shinylive")

# Export to static site (run from parent directory)
shinylive::export(appdir = ".", destdir = "docs")
```

Then in your GitHub repo: **Settings → Pages → Deploy from branch → main / docs**

The app will be live at `https://yourusername.github.io/statistical-analysis/`

---

## Dependencies

```r
library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(MASS)
library(nnet)
library(readxl)
```

---

## Citation

If you use this tool in teaching or research, please cite:

> Subedi, R. (2026). *BioStat Explorer: An interactive biostatistics teaching tool*. GitHub. https://github.com/rajsubediresearch/statistical-analysis

---

## Author

**Raj Subedi**  
PhD Student · Epidemiology · Georgia State University  
2CI Fellow · Graduate Research Assistant, Chowell Lab  
🌐 [rajsubediresearch.com](https://rajsubediresearch.com) · [Google Scholar](https://scholar.google.com) · [GitHub](https://github.com/rajsubediresearch)
