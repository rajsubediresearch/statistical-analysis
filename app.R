library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(MASS)
library(nnet)
library(readxl)
library(gridExtra)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !identical(a, "")) a else b

# ── Descriptive statistics helper ─────────────────────────────────────────────
make_desc_stats <- function(df, vars = NULL) {
  if (!is.null(vars) && length(vars) > 0) df <- df[, vars, drop = FALSE]

  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  cat_df <- df[, !sapply(df, is.numeric), drop = FALSE]

  out <- list()

  if (ncol(num_df) > 0) {
    num_tbl <- do.call(rbind, lapply(names(num_df), function(v) {
      x <- na.omit(num_df[[v]])
      data.frame(
        Variable = v,
        N        = length(x),
        Missing  = sum(is.na(num_df[[v]])),
        Mean     = round(mean(x), 3),
        SD       = round(sd(x), 3),
        Median   = round(median(x), 3),
        IQR      = round(IQR(x), 3),
        Min      = round(min(x), 3),
        Max      = round(max(x), 3),
        check.names = FALSE
      )
    }))
    out$numeric <- num_tbl
  }

  if (ncol(cat_df) > 0) {
    cat_tbl <- do.call(rbind, lapply(names(cat_df), function(v) {
      x  <- cat_df[[v]]
      tb <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
      do.call(rbind, lapply(names(tb), function(lv) {
        data.frame(
          Variable   = v,
          Level      = ifelse(is.na(lv) | lv == "NA", "<NA>", lv),
          N          = as.integer(tb[[lv]]),
          `% (valid)` = round(100 * as.integer(tb[[lv]]) / sum(!is.na(x)), 1),
          check.names = FALSE
        )
      }))
    }))
    out$categorical <- cat_tbl
  }

  out
}

# ── Sample datasets ───────────────────────────────────────────────────────────
make_epi_data <- function() {
  set.seed(42)
  n <- 300
  age        <- round(rnorm(n, 45, 15))
  sex        <- factor(sample(c("Male","Female"), n, replace=TRUE))
  bmi        <- round(rnorm(n, 26, 5), 1)
  smoking    <- factor(sample(c("Never","Former","Current"), n, replace=TRUE,
                              prob=c(.5,.3,.2)))
  vaccinated <- rbinom(n, 1, .65)
  sbp        <- round(110 + 0.4*age + 3*(sex=="Male") + 0.8*bmi +
                        5*(smoking=="Current") + rnorm(n,0,10))
  disease    <- rbinom(n, 1, plogis(-3 + 0.03*age + 0.5*(smoking=="Current") -
                                      0.8*vaccinated))
  severity   <- cut(sbp, breaks=c(-Inf,120,140,Inf),
                    labels=c("Normal","Pre-hypertensive","Hypertensive"),
                    ordered=TRUE)
  transport  <- factor(sample(c("Walk","Bike","Car","Transit"), n,
                               replace=TRUE, prob=c(.2,.15,.4,.25)))
  data.frame(age, sex, bmi, smoking, vaccinated, sbp, disease,
             severity, transport, stringsAsFactors=FALSE)
}

make_clinical_data <- function() {
  set.seed(99)
  n <- 250
  treatment   <- factor(sample(c("Control","DrugA","DrugB"), n, replace=TRUE))
  age         <- round(rnorm(n, 55, 12))
  weight      <- round(rnorm(n, 75, 15), 1)
  cholesterol <- round(180 + 0.3*age + 2*(treatment=="Control") + rnorm(n,0,20))
  bp_change   <- round(-5*(treatment=="DrugA") - 8*(treatment=="DrugB") +
                         0.1*age + rnorm(n,0,5), 1)
  recovered   <- rbinom(n, 1, plogis(-1 + 1.2*(treatment=="DrugA") +
                                       1.8*(treatment=="DrugB") - 0.02*age))
  outcome_cat <- cut(bp_change, breaks=c(-Inf,-5,0,Inf),
                     labels=c("Improved","Stable","Worsened"), ordered=TRUE)
  data.frame(treatment, age, weight, cholesterol, bp_change,
             recovered, outcome_cat, stringsAsFactors=FALSE)
}

SAMPLE_DATASETS <- list(
  "Epidemiology Dataset (n=300)"   = make_epi_data(),
  "Clinical Trial Dataset (n=250)" = make_clinical_data()
)

# ── Interpretation helpers ────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < .001) "< 0.001" else paste0("= ", round(p, 3))
}

interp_lm <- function(exposure, outcome, covariates, coef, ci_lo, ci_hi, p) {
  adj <- if (length(covariates) > 0)
    paste0(", adjusting for ", paste(covariates, collapse=", ")) else ""
  dir <- if (coef >= 0) "increase" else "decrease"
  if (p < .05) {
    paste0("<b>", exposure, "</b> is associated with a <b>", round(abs(coef),3),
           " unit ", dir, "</b> in <b>", outcome, "</b> (\u03b2 = ", round(coef,3),
           ", 95% CI: ", round(ci_lo,3), " to ", round(ci_hi,3),
           ", p ", fmt_p(p), ")", adj, ".")
  } else {
    paste0("<b>", exposure, "</b> was <b>not significantly associated</b> with <b>",
           outcome, "</b> (\u03b2 = ", round(coef,3),
           ", 95% CI: ", round(ci_lo,3), " to ", round(ci_hi,3),
           ", p ", fmt_p(p), ")", adj, ".")
  }
}

interp_logistic <- function(exposure, outcome, covariates, or, ci_lo, ci_hi, p) {
  adj <- if (length(covariates) > 0)
    paste0(", adjusting for ", paste(covariates, collapse=", ")) else ""
  dir <- if (or >= 1) "higher" else "lower"
  pct <- round(abs(1 - or) * 100)
  if (p < .05) {
    paste0("<b>", exposure, "</b> is associated with <b>", pct, "% ", dir,
           " odds</b> of <b>", outcome, "</b> (OR = ", round(or,3),
           ", 95% CI: ", round(ci_lo,3), " to ", round(ci_hi,3),
           ", p ", fmt_p(p), ")", adj, ".")
  } else {
    paste0("<b>", exposure, "</b> was <b>not significantly associated</b> with the odds of <b>",
           outcome, "</b> (OR = ", round(or,3),
           ", 95% CI: ", round(ci_lo,3), " to ", round(ci_hi,3),
           ", p ", fmt_p(p), ")", adj, ".")
  }
}

interp_ttest <- function(g1, g2, outcome, diff, ci_lo, ci_hi, p) {
  sig <- if (p < .05) "significantly " else "not significantly "
  paste0("Mean <b>", outcome, "</b> was <b>", sig, "different</b> between <b>",
         g1, "</b> and <b>", g2, "</b> (mean difference = ", round(diff,3),
         ", 95% CI: ", round(ci_lo,3), " to ", round(ci_hi,3),
         ", p ", fmt_p(p), ").")
}

interp_anova <- function(outcome, groupvar, p) {
  sig <- if (p < .05)
    paste0("at least one group differs significantly in mean <b>", outcome, "</b>")
  else
    paste0("no significant difference in mean <b>", outcome,
           "</b> was detected across groups")
  paste0("The one-way ANOVA indicates that <b>", sig,
         "</b> across levels of <b>", groupvar, "</b> (p ", fmt_p(p), ").")
}

interp_ordinal <- function(exposure, outcome, covariates, or, ci_lo, ci_hi, p) {
  adj <- if (length(covariates) > 0)
    paste0(", adjusting for ", paste(covariates, collapse=", ")) else ""
  dir <- if (or >= 1) "higher" else "lower"
  paste0("<b>", exposure, "</b> is associated with <b>", round(or,3), " times ",
         dir, " odds</b> of being in a more severe category of <b>", outcome,
         "</b> (OR = ", round(or,3), ", 95% CI: ", round(ci_lo,3), " to ",
         round(ci_hi,3), ", p ", fmt_p(p), ")", adj, ".")
}

# ── Analysis plot helper ──────────────────────────────────────────────────────
make_analysis_plot <- function(df, mt, outcome, exposure) {
  tryCatch({
    base <- ggplot(df) + theme_minimal(base_size=13) +
      theme(plot.title=element_text(face="bold", size=13))
    if (mt %in% c("lm_simple","lm_multi")) {
      if (is.numeric(df[[exposure]])) {
        base + aes_string(x=exposure, y=outcome) +
          geom_point(alpha=.5, color="#2c7bb6") +
          geom_smooth(method="lm", se=TRUE, color="#e74c3c") +
          labs(title=paste("Linear fit:", outcome, "~", exposure))
      } else {
        base + aes_string(x=exposure, y=outcome, fill=exposure) +
          geom_boxplot(alpha=.7) + theme(legend.position="none") +
          labs(title=paste(outcome, "by", exposure))
      }
    } else if (mt %in% c("logit_simple","logit_multi")) {
      if (is.numeric(df[[exposure]])) {
        base + aes_string(x=exposure, y=outcome) +
          geom_jitter(height=.05, alpha=.4, color="#2c7bb6") +
          geom_smooth(method="glm", method.args=list(family=binomial),
                      se=TRUE, color="#e74c3c") +
          labs(title=paste("Logistic fit:", outcome, "~", exposure))
      } else {
        tmp <- as.data.frame(table(df[[exposure]], df[[outcome]]))
        names(tmp) <- c("Group","Outcome","Count")
        ggplot(tmp, aes(x=Group, y=Count, fill=Outcome)) +
          geom_bar(stat="identity", position="fill", alpha=.8) +
          scale_y_continuous(labels=scales::percent) +
          labs(title=paste("Outcome proportion by", exposure), y="Proportion") +
          theme_minimal(base_size=13)
      }
    } else if (mt == "ttest") {
      base + aes_string(x=exposure, y=outcome, fill=exposure) +
        geom_violin(alpha=.5) + geom_boxplot(width=.2, alpha=.8) +
        theme(legend.position="none") +
        labs(title=paste(outcome, "by", exposure))
    } else if (mt == "anova") {
      base + aes_string(x=exposure, y=outcome, fill=exposure) +
        geom_boxplot(alpha=.7) +
        stat_summary(fun=mean, geom="point", shape=18, size=3, color="#e74c3c") +
        theme(legend.position="none") +
        labs(title=paste(outcome, "by", exposure), caption="Red diamond = group mean")
    } else if (mt == "ordinal") {
      tmp <- as.data.frame(table(df[[exposure]], df[[outcome]]))
      names(tmp) <- c("Exposure","Severity","Count")
      ggplot(tmp, aes(x=Exposure, y=Count, fill=Severity)) +
        geom_bar(stat="identity", position="fill", alpha=.8) +
        scale_y_continuous(labels=scales::percent) +
        labs(title=paste("Severity by", exposure), y="Proportion") +
        theme_minimal(base_size=13)
    } else {
      ggplot() + annotate("text", x=.5, y=.5,
        label="No auto-plot for this model type.", size=5, color="gray50") +
        theme_void()
    }
  }, error=function(e) {
    ggplot() + annotate("text", x=.5, y=.5,
      label=paste("Plot error:", e$message), size=4, color="red") + theme_void()
  })
}

# ── Post-hoc Tukey helper ─────────────────────────────────────────────────────
make_tukey_table <- function(aov_mod) {
  tk  <- TukeyHSD(aov_mod)
  mat <- tk[[1]]
  data.frame(
    Comparison   = rownames(mat),
    Difference   = round(mat[, "diff"],  4),
    `CI Low`     = round(mat[, "lwr"],   4),
    `CI High`    = round(mat[, "upr"],   4),
    `p adjusted` = sapply(mat[, "p adj"], function(p)
      if (p < .001) "<0.001" else as.character(round(p, 4))),
    check.names  = FALSE
  )
}

# ── Regression diagnostics helpers ───────────────────────────────────────────
make_lm_diag_plots <- function(mod) {
  df_diag <- data.frame(
    fitted    = fitted(mod),
    resid     = residuals(mod),
    std_resid = rstandard(mod),
    leverage  = hatvalues(mod),
    cooks     = cooks.distance(mod),
    obs       = seq_along(fitted(mod))
  )
  p1 <- ggplot(df_diag, aes(fitted, resid)) +
    geom_point(alpha = .5, color = "#2c7bb6") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#e74c3c") +
    geom_smooth(method = "loess", se = FALSE, color = "#f39c12", linewidth = .8) +
    labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
    theme_minimal(base_size = 12)

  p2 <- ggplot(df_diag, aes(sample = std_resid)) +
    stat_qq(alpha = .5, color = "#2c7bb6") +
    stat_qq_line(color = "#e74c3c") +
    labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Std. residuals") +
    theme_minimal(base_size = 12)

  p3 <- ggplot(df_diag, aes(fitted, sqrt(abs(std_resid)))) +
    geom_point(alpha = .5, color = "#2c7bb6") +
    geom_smooth(method = "loess", se = FALSE, color = "#f39c12", linewidth = .8) +
    labs(title = "Scale-Location", x = "Fitted values",
         y = expression(sqrt("|Std. residuals|"))) +
    theme_minimal(base_size = 12)

  p4 <- ggplot(df_diag, aes(leverage, std_resid, size = cooks)) +
    geom_point(alpha = .6, color = "#2c7bb6") +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#e74c3c") +
    geom_vline(xintercept = 2 * mean(df_diag$leverage),
               linetype = "dashed", color = "#f39c12") +
    labs(title = "Leverage vs Std. Residuals", x = "Leverage",
         y = "Std. residuals", size = "Cook's D") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

make_logit_diag_plots <- function(mod) {
  df_diag <- data.frame(
    fitted    = fitted(mod),
    resid_dev = residuals(mod, type = "deviance"),
    resid_prs = residuals(mod, type = "pearson"),
    leverage  = hatvalues(mod),
    cooks     = cooks.distance(mod),
    obs       = seq_along(fitted(mod)),
    y         = mod$y
  )

  # Calibration: robust binning — at least 10 obs per bin
  n        <- nrow(df_diag)
  n_bins   <- max(3, min(10, floor(n / 15)))
  df_diag$bin <- tryCatch({
    cuts <- unique(quantile(df_diag$fitted,
                            probs = seq(0, 1, length.out = n_bins + 1),
                            na.rm = TRUE))
    if (length(cuts) < 3) cuts <- c(0, 0.5, 1)
    cut(df_diag$fitted, breaks = cuts, include.lowest = TRUE)
  }, error = function(e) factor(rep(1, n)))

  cal <- tryCatch({
    tmp <- aggregate(cbind(fitted, y) ~ bin, data = df_diag, FUN = mean, na.rm = TRUE)
    tmp[complete.cases(tmp), ]
  }, error = function(e) data.frame(fitted = numeric(0), y = numeric(0)))

  p1 <- if (nrow(cal) >= 2) {
    ggplot(cal, aes(fitted, y)) +
      geom_point(size = 3, color = "#2c7bb6") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#e74c3c") +
      geom_smooth(method = "loess", se = FALSE, color = "#f39c12",
                  linewidth = .7, na.rm = TRUE) +
      labs(title = "Calibration plot", x = "Mean predicted probability",
           y = "Observed proportion") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_minimal(base_size = 12)
  } else {
    ggplot() +
      annotate("text", x = .5, y = .5, label = "Not enough variation\nfor calibration plot",
               size = 4, color = "gray50") + theme_void() +
      labs(title = "Calibration plot")
  }

  p2 <- ggplot(df_diag, aes(sample = resid_dev)) +
    stat_qq(alpha = .5, color = "#2c7bb6") +
    stat_qq_line(color = "#e74c3c") +
    labs(title = "Normal Q-Q (deviance resid.)",
         x = "Theoretical quantiles", y = "Deviance residuals") +
    theme_minimal(base_size = 12)

  p3 <- ggplot(df_diag, aes(obs, cooks)) +
    geom_segment(aes(xend = obs, yend = 0), color = "#2c7bb6", alpha = .6) +
    geom_point(color = "#2c7bb6", size = 1.5) +
    geom_hline(yintercept = 4 / n,
               linetype = "dashed", color = "#e74c3c") +
    labs(title = "Cook's Distance", x = "Observation", y = "Cook's D",
         caption = "Dashed line: 4/n threshold") +
    theme_minimal(base_size = 12)

  p4 <- ggplot(df_diag, aes(leverage, resid_prs)) +
    geom_point(alpha = .5, color = "#2c7bb6") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#e74c3c") +
    geom_vline(xintercept = 2 * mean(df_diag$leverage, na.rm = TRUE),
               linetype = "dashed", color = "#f39c12") +
    labs(title = "Leverage vs Pearson Residuals",
         x = "Leverage", y = "Pearson residuals") +
    theme_minimal(base_size = 12)

  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  theme = bs_theme(bootswatch="flatly", primary="#2c7bb6"),

  tags$style(HTML("
    body { margin:0; padding:0; background:#f5f6fa; }
    .topnav {
      background:#2c7bb6;
      display:flex; align-items:center;
      padding:0 1.2rem; height:50px;
      box-shadow:0 2px 4px rgba(0,0,0,0.15);
    }
    .topnav .brand {
      color:white; font-weight:700; font-size:1.15rem;
      margin-right:1.5rem; white-space:nowrap;
    }
    .topnav .nav-links { display:flex; gap:4px; flex:1; }
    .topnav .nav-links a {
      color:rgba(255,255,255,0.85);
      padding:6px 14px; border-radius:5px;
      text-decoration:none; font-size:0.92rem;
      transition:background 0.15s;
      cursor:pointer;
    }
    .topnav .nav-links a:hover,
    .topnav .nav-links a.active {
      background:rgba(255,255,255,0.22); color:white;
    }
    .topnav .brand-link {
      color:rgba(255,255,255,0.75); font-size:0.82rem;
      text-decoration:none; white-space:nowrap;
    }
    .topnav .brand-link:hover { color:white; }
    .page-wrap {
      display:flex; min-height:calc(100vh - 50px);
    }
    .left-sidebar {
      width:290px; min-width:290px;
      background:white;
      border-right:1px solid #dee2e6;
      padding:1.2rem 1rem;
      overflow-y:auto;
    }
    .right-content {
      flex:1; padding:1.2rem 1.4rem;
      overflow-x:hidden;
    }
    .card { border:1px solid #dee2e6; border-radius:8px;
            background:white; margin-bottom:1rem; }
    .card-header {
      background:#f8f9fa; border-bottom:1px solid #dee2e6;
      padding:0.6rem 1rem; font-weight:600; border-radius:8px 8px 0 0;
    }
    .card-body { padding:1rem; }
  ")),

  # ── Top nav ────────────────────────────────────────────────────────────────
  tags$div(class="topnav",
    tags$span("BioStat Explorer", class="brand"),
    tags$div(class="nav-links",
      actionLink("go_data",     "📂 Data"),
      actionLink("go_analysis", "📊 Analysis"),
      actionLink("go_plots",    "📈 Plots")
    ),
    tags$a("Raj Subedi | rajsubediresearch.com",
           href="https://rajsubediresearch.com",
           target="_blank", class="brand-link")
  ),

  # ── Page body ──────────────────────────────────────────────────────────────
  tags$div(class="page-wrap",
    tags$div(class="left-sidebar", uiOutput("sidebar_ui")),
    tags$div(class="right-content", uiOutput("main_ui"))
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  current_tab <- reactiveVal("data")
  observeEvent(input$go_data,     { current_tab("data") })
  observeEvent(input$go_analysis, { current_tab("analysis") })
  observeEvent(input$go_plots,    { current_tab("plots") })

  # ── Dataset ──────────────────────────────────────────────────────────────
  # Sheet names for Excel files
  excel_sheets <- reactive({
    req(input$user_file)
    ext <- tolower(tools::file_ext(input$user_file$name))
    if (ext %in% c("xlsx","xls")) {
      tryCatch(readxl::excel_sheets(input$user_file$datapath),
               error=function(e) "Sheet1")
    } else NULL
  })

  output$sheet_selector_ui <- renderUI({
    sheets <- excel_sheets()
    req(sheets)
    selectInput("excel_sheet", "Select sheet:", choices=sheets)
  })

  # Dynamic upload options — separator for CSV/TXT, sheet for Excel
  output$upload_options_ui <- renderUI({
    req(input$user_file)
    ext <- tolower(tools::file_ext(input$user_file$name))
    if (ext %in% c("xlsx","xls")) {
      uiOutput("sheet_selector_ui")
    } else {
      radioButtons("sep", "Separator:",
                   choices=c(Comma=",", Tab="\t", Semicolon=";"),
                   selected=",", inline=TRUE)
    }
  })

  active_data <- reactive({
    src <- input$data_source %||% "sample"
    if (src == "upload") {
      req(input$user_file)
      ext <- tolower(tools::file_ext(input$user_file$name))
      if (ext %in% c("xlsx","xls")) {
        sheet <- input$excel_sheet %||% 1
        tryCatch(
          as.data.frame(readxl::read_excel(
            input$user_file$datapath,
            sheet = sheet,
            na    = c("","NA","N/A")
          )),
          error=function(e) {
            validate(need(FALSE, paste("Could not read Excel file:", e$message)))
          }
        )
      } else {
        read.csv(input$user_file$datapath,
                 sep              = input$sep %||% ",",
                 stringsAsFactors = FALSE,
                 na.strings       = c("","NA","N/A"))
      }
    } else {
      nm <- input$sample_choice %||% names(SAMPLE_DATASETS)[1]
      SAMPLE_DATASETS[[nm]]
    }
  })

  num_cols <- reactive({ names(Filter(is.numeric, active_data())) })
  cat_cols <- reactive({ names(Filter(Negate(is.numeric), active_data())) })
  all_cols <- reactive({ names(active_data()) })
  ord_cols <- reactive({
    df <- active_data()
    names(df)[sapply(df, function(x) is.factor(x) && is.ordered(x))]
  })

  # ── Sidebar ───────────────────────────────────────────────────────────────
  output$sidebar_ui <- renderUI({
    switch(current_tab(),

      "data" = tagList(
        h5("Load Data", class="fw-bold mb-3"),
        radioButtons("data_source", NULL,
                     choices=c("Use sample dataset"="sample",
                               "Upload my own"="upload")),
        conditionalPanel("input.data_source == 'sample'",
          selectInput("sample_choice", "Choose dataset:",
                      choices=names(SAMPLE_DATASETS))
        ),
        conditionalPanel("input.data_source == 'upload'",
          fileInput("user_file", "Upload CSV, TXT, or Excel",
                    accept=c(".csv",".txt",".xlsx",".xls")),
          uiOutput("upload_options_ui")
        ),
        hr(),
        uiOutput("data_summary_ui")
      ),

      "analysis" = tagList(
        h5("Model Setup", class="fw-bold mb-3"),
        tags$div(
          tags$label("Statistical Method:", class="control-label"),
          tags$select(
            id = "model_type",
            class = "form-control",
            tags$optgroup(label = "Descriptive",
              tags$option(value = "desc", "Descriptive Statistics")
            ),
            tags$optgroup(label = "Regression",
              tags$option(value = "lm_simple",   "Simple Linear Regression"),
              tags$option(value = "lm_multi",    "Multiple Linear Regression"),
              tags$option(value = "logit_simple","Simple Logistic Regression"),
              tags$option(value = "logit_multi", "Multiple Logistic Regression"),
              tags$option(value = "ordinal",     "Ordinal Logistic Regression"),
              tags$option(value = "multinom",    "Multinomial Logistic Regression")
            ),
            tags$optgroup(label = "Group Comparison",
              tags$option(value = "ttest", "Independent t-test"),
              tags$option(value = "anova", "One-way ANOVA")
            )
          )
        ),
        hr(),
        uiOutput("var_selectors"),
        hr(),
        conditionalPanel(
          "input.model_type !== 'desc'",
          checkboxInput("show_plot",
                        "📈 Show analysis plot below results",
                        value=FALSE)
        ),
        br(),
        uiOutput("run_btn_ui")
      ),

      "plots" = tagList(
        h5("Plot Options", class="fw-bold mb-3"),
        uiOutput("plot_controls"),
        br(),
        actionButton("make_plot", "Generate Plot",
                     class="btn-primary w-100")
      )
    )
  })

  # ── Main panel ────────────────────────────────────────────────────────────
  output$main_ui <- renderUI({
    switch(current_tab(),

      "data" = tagList(
        tags$div(class="card",
          tags$div(class="card-header", "Data Preview"),
          tags$div(class="card-body", DTOutput("data_preview"))
        )
      ),

      "analysis" = tagList(
        tags$div(class="card",
          tags$div(class="card-header", "Results"),
          tags$div(class="card-body", uiOutput("results_ui"))
        ),
        uiOutput("posthoc_card"),
        uiOutput("diagnostics_card"),
        uiOutput("analysis_plot_card")
      ),

      "plots" = tagList(
        tags$div(class="card",
          tags$div(class="card-header", "Visualization"),
          tags$div(class="card-body",
                   plotOutput("main_plot", height="450px"))
        )
      )
    )
  })

  # ── Run button ────────────────────────────────────────────────────────────
  output$run_btn_ui <- renderUI({
    lbl <- if (!is.null(input$model_type) && input$model_type == "desc")
      "▶ Compute Statistics" else "▶ Run Analysis"
    actionButton("run_model", lbl, class = "btn-primary w-100")
  })

  output$desc_num_table <- renderDT({
    req(model_result(), input$model_type == "desc")
    res <- model_result()$result
    req(res$numeric)
    datatable(res$numeric, rownames = FALSE,
              options = list(pageLength = 20, dom = "t", scrollX = TRUE))
  })

  output$desc_cat_table <- renderDT({
    req(model_result(), input$model_type == "desc")
    res <- model_result()$result
    req(res$categorical)
    datatable(res$categorical, rownames = FALSE,
              options = list(pageLength = 30, dom = "tp", scrollX = TRUE))
  })

  # ── Data tab outputs ──────────────────────────────────────────────────────
  output$data_summary_ui <- renderUI({
    df <- active_data()
    tagList(
      tags$small(class="text-muted",
        paste0("Rows: ", nrow(df), " | Columns: ", ncol(df))), br(),
      tags$small(class="text-muted",
        paste0("Numeric: ", sum(sapply(df, is.numeric)),
               " | Categorical: ", sum(!sapply(df, is.numeric))))
    )
  })

  output$data_preview <- renderDT(
    datatable(active_data(),
              options=list(pageLength=10, scrollX=TRUE),
              rownames=FALSE)
  )

  # ── Variable selectors ───────────────────────────────────────────────────
  output$var_selectors <- renderUI({
    mt <- input$model_type %||% "lm_simple"
    nc <- num_cols(); ac <- all_cols(); cc <- cat_cols()

    def_out <- if (length(nc) >= 1) nc[1] else ac[1]
    def_exp <- if (length(nc) >= 2) nc[2] else {
      others <- ac[ac != def_out]; if (length(others) > 0) others[1] else ac[1]
    }
    def_grp <- if (length(cc) >= 1) cc[1] else ac[1]
    def_ord <- if (length(ord_cols()) > 0) ord_cols()[1] else ac[1]

    switch(mt,
      desc         = tagList(
        tags$small(class = "text-muted",
          "Select variables to summarise, or leave blank for all columns."),
        br(), br(),
        selectInput("desc_vars", "Variables (optional):",
                    choices = ac, multiple = TRUE, selectize = TRUE)
      ),
      lm_simple    = tagList(
        selectInput("outcome",  "Outcome (continuous):", nc, selected=def_out),
        selectInput("exposure", "Exposure:", ac, selected=def_exp)),
      lm_multi     = tagList(
        selectInput("outcome",    "Outcome (continuous):", nc, selected=def_out),
        selectInput("exposure",   "Main exposure:", ac, selected=def_exp),
        selectInput("covariates", "Covariates (ctrl+click):", ac, multiple=TRUE)),
      logit_simple = tagList(
        selectInput("outcome",  "Outcome (binary 0/1):", ac, selected=def_out),
        selectInput("exposure", "Exposure:", ac, selected=def_exp)),
      logit_multi  = tagList(
        selectInput("outcome",    "Outcome (binary 0/1):", ac, selected=def_out),
        selectInput("exposure",   "Main exposure:", ac, selected=def_exp),
        selectInput("covariates", "Covariates (ctrl+click):", ac, multiple=TRUE)),
      ordinal      = tagList(
        selectInput("outcome",    "Outcome (ordered factor):",
                    c(ord_cols(), ac), selected=def_ord),
        selectInput("exposure",   "Main exposure:", ac, selected=def_exp),
        selectInput("covariates", "Covariates (ctrl+click):", ac, multiple=TRUE)),
      multinom     = tagList(
        selectInput("outcome",    "Outcome (categorical):", ac, selected=def_grp),
        selectInput("exposure",   "Main exposure:", ac, selected=def_exp),
        selectInput("covariates", "Covariates (ctrl+click):", ac, multiple=TRUE)),
      ttest        = tagList(
        selectInput("outcome",  "Continuous variable:", nc, selected=def_out),
        selectInput("exposure", "Grouping variable (2 groups):", ac, selected=def_grp)),
      anova        = tagList(
        selectInput("outcome",  "Continuous variable:", nc, selected=def_out),
        selectInput("exposure", "Grouping variable:", ac, selected=def_grp))
    )
  })

  # ── Model fitting ─────────────────────────────────────────────────────────
  model_result <- eventReactive(input$run_model, {
    df  <- active_data()
    mt  <- input$model_type

    # Descriptive stats path — no outcome/exposure needed
    if (mt == "desc") {
      return(list(.__type__ = "desc",
                  result    = make_desc_stats(df, input$desc_vars)))
    }

    out <- input$outcome
    exp <- input$exposure
    cov <- input$covariates

    validate(
      need(!is.null(out) && !is.null(exp), "Select outcome and exposure."),
      need(out != exp,  "Outcome and exposure must be different."),
      need(!out %in% cov, "Outcome cannot also be a covariate.")
    )

    # Wrap names in backticks to handle spaces/special characters in column names
    bt  <- function(x) paste0('`', x, '`')
    rhs <- if (length(cov) > 0 &&
               mt %in% c("lm_multi","logit_multi","ordinal","multinom"))
      paste(c(bt(exp), sapply(cov, bt)), collapse=" + ") else bt(exp)
    fml <- as.formula(paste(bt(out), "~", rhs))

    tryCatch(switch(mt,
      lm_simple    = , lm_multi    = lm(fml, data=df),
      logit_simple = , logit_multi = glm(fml, data=df, family=binomial),
      ordinal      = polr(fml, data=df, Hess=TRUE),
      multinom     = multinom(fml, data=df, trace=FALSE),
      ttest = {
        grps <- unique(na.omit(df[[exp]]))
        validate(need(length(grps)==2,
          "t-test requires exactly 2 groups."))
        t.test(df[[out]] ~ df[[exp]])
      },
      anova = aov(fml, data=df)
    ), error=function(e) validate(need(FALSE, paste("Error:", e$message))))
  })

  # ── Results UI ────────────────────────────────────────────────────────────
  output$results_ui <- renderUI({
    req(model_result())
    mt  <- input$model_type
    mod <- model_result()

    # ── Descriptive statistics output ───────────────────────────────────────
    if (mt == "desc") {
      res <- mod$result
      cards <- list()
      if (!is.null(res$numeric)) {
        cards <- c(cards, list(tagList(
          h6("Continuous Variables", class = "fw-bold"),
          DTOutput("desc_num_table"), br()
        )))
      }
      if (!is.null(res$categorical)) {
        cards <- c(cards, list(tagList(
          h6("Categorical Variables", class = "fw-bold"),
          DTOutput("desc_cat_table")
        )))
      }
      return(tagList(cards))
    }

    out <- input$outcome
    exp <- input$exposure
    cov <- if (!is.null(input$covariates) &&
               mt %in% c("lm_multi","logit_multi","ordinal","multinom"))
             input$covariates else character(0)

    result <- tryCatch({
      if (mt %in% c("lm_simple","lm_multi")) {
        s  <- summary(mod); cf <- coef(s); ci <- confint(mod)
        erows <- grep(paste0("^", exp), rownames(cf), value=TRUE)
        interps <- lapply(erows, function(rn)
          interp_lm(rn, out, cov, cf[rn,1], ci[rn,1], ci[rn,2], cf[rn,4]))
        list(interps=interps,
             extra=paste0("R\u00b2 = ", round(s$r.squared,4),
                          " | Adj R\u00b2 = ", round(s$adj.r.squared,4)))

      } else if (mt %in% c("logit_simple","logit_multi")) {
        s  <- summary(mod); cf <- coef(s)
        ci <- suppressMessages(confint(mod))
        or <- exp(cf[,1]); ci_or <- exp(ci)
        erows <- grep(paste0("^", exp), rownames(cf), value=TRUE)
        interps <- lapply(erows, function(rn)
          interp_logistic(rn, out, cov,
                          or[rn], ci_or[rn,1], ci_or[rn,2], cf[rn,4]))
        list(interps=interps, extra=NULL)

      } else if (mt == "ordinal") {
        cf <- coef(summary(mod))
        ci <- suppressMessages(confint(mod))
        or <- exp(cf[,1]); ci_or <- exp(ci)
        erows <- grep(paste0("^", exp), rownames(cf), value=TRUE)
        interps <- lapply(erows, function(rn) {
          p <- 2*pnorm(-abs(cf[rn,3]))
          interp_ordinal(rn, out, cov,
                         or[rn], ci_or[rn,1], ci_or[rn,2], p)
        })
        list(interps=interps, extra=NULL)

      } else if (mt == "multinom") {
        cf  <- summary(mod)$coefficients
        se  <- summary(mod)$standard.errors
        interps <- unlist(lapply(rownames(cf), function(lv) {
          z <- cf[lv,]/se[lv,]; p <- 2*pnorm(-abs(z))
          erows <- grep(paste0("^", exp), names(z), value=TRUE)
          lapply(erows, function(rn) {
            or    <- exp(cf[lv,rn])
            ci_lo <- exp(cf[lv,rn] - 1.96*se[lv,rn])
            ci_hi <- exp(cf[lv,rn] + 1.96*se[lv,rn])
            paste0("<em>Level: ", lv, " vs reference —</em> ",
                   interp_logistic(rn, out, cov, or, ci_lo, ci_hi, p[rn]))
          })
        }), recursive=FALSE)
        list(interps=interps, extra=NULL)

      } else if (mt == "ttest") {
        tt   <- mod
        diff <- tt$estimate[1] - tt$estimate[2]
        list(interps=list(interp_ttest(
               names(tt$estimate)[1], names(tt$estimate)[2],
               out, diff, tt$conf.int[1], tt$conf.int[2], tt$p.value)),
             extra=paste0("t = ", round(tt$statistic,4),
                          " | df = ", round(tt$parameter,1),
                          " | p ", fmt_p(tt$p.value)))

      } else if (mt == "anova") {
        s <- summary(mod)[[1]]; p <- s[["Pr(>F)"]][1]
        list(interps=list(interp_anova(out, exp, p)), extra=NULL)
      }
    }, error=function(e) {
      list(interps=list(paste("Could not compute:", e$message)), extra=NULL)
    })

    tagList(
      if (!is.null(result$extra))
        div(class="alert alert-info mb-2",
            HTML(paste0("<strong>Model fit:</strong> ", result$extra))),
      h6("Coefficient Table", class="fw-bold"),
      DTOutput("coef_table"),
      hr(),
      h6("📝 Interpretation", class="fw-bold"),
      div(class="alert alert-success",
        lapply(result$interps, function(i) tagList(HTML(i), br(), br()))
      ),
      div(class="alert alert-warning",
        tags$small("⚠️ This tool reports statistical associations only. ",
                   "Causal interpretation requires careful consideration of ",
                   "study design, confounding, and biological plausibility.")
      )
    )
  })

  output$coef_table <- renderDT({
    req(model_result())
    mt <- input$model_type
    if (mt == "desc") return(datatable(data.frame()))
    mod <- model_result()
    tbl <- tryCatch({
      if (mt %in% c("lm_simple","lm_multi")) {
        cf <- coef(summary(mod)); ci <- confint(mod)
        data.frame(Term=rownames(cf), Estimate=round(cf[,1],4),
                   SE=round(cf[,2],4), `t value`=round(cf[,3],4),
                   `p value`=sapply(cf[,4], function(p)
                     if(p<.001)"<0.001" else as.character(round(p,4))),
                   `CI Low`=round(ci[,1],4), `CI High`=round(ci[,2],4),
                   check.names=FALSE)
      } else if (mt %in% c("logit_simple","logit_multi")) {
        cf <- coef(summary(mod)); ci <- suppressMessages(confint(mod))
        or <- exp(cf[,1]); ci_or <- exp(ci)
        data.frame(Term=rownames(cf), OR=round(or,4),
                   `CI Low`=round(ci_or[,1],4), `CI High`=round(ci_or[,2],4),
                   `p value`=sapply(cf[,4], function(p)
                     if(p<.001)"<0.001" else as.character(round(p,4))),
                   check.names=FALSE)
      } else if (mt == "ordinal") {
        cf <- coef(summary(mod)); ci <- suppressMessages(confint(mod))
        or <- exp(cf[,1]); ci_or <- exp(ci)
        data.frame(Term=rownames(cf), OR=round(or,4),
                   `CI Low`=round(ci_or[,1],4), `CI High`=round(ci_or[,2],4),
                   check.names=FALSE)
      } else if (mt == "multinom") {
        as.data.frame(round(exp(summary(mod)$coefficients),4))
      } else if (mt == "ttest") {
        tt <- mod
        data.frame(Group=names(tt$estimate), Mean=round(tt$estimate,4),
                   check.names=FALSE)
      } else if (mt == "anova") {
        as.data.frame(round(summary(mod)[[1]],4))
      }
    }, error=function(e) data.frame(Error=e$message))
    datatable(tbl, rownames=TRUE,
              options=list(pageLength=20, dom="t", scrollX=TRUE))
  })

  output$analysis_plot_card <- renderUI({
    req(input$run_model > 0, isTRUE(input$show_plot))
    tags$div(class="card",
      tags$div(class="card-header", "Analysis Plot"),
      tags$div(class="card-body",
               plotOutput("analysis_plot", height="380px"))
    )
  })

  output$analysis_plot <- renderPlot({
    req(input$run_model, isTRUE(input$show_plot), model_result())
    make_analysis_plot(active_data(), input$model_type,
                       input$outcome, input$exposure)
  })

  # ── Post-hoc Tukey (ANOVA only) ───────────────────────────────────────────
  output$posthoc_card <- renderUI({
    req(model_result())
    if (input$model_type != "anova") return(NULL)
    tags$div(class = "card",
      tags$div(class = "card-header", "🔍 Post-hoc Tests (Tukey HSD)"),
      tags$div(class = "card-body",
        p(class = "text-muted",
          tags$small("Pairwise comparisons with Bonferroni-corrected p-values via Tukey's Honestly Significant Difference.")),
        DTOutput("posthoc_table")
      )
    )
  })

  output$posthoc_table <- renderDT({
    req(model_result(), input$model_type == "anova")
    tbl <- tryCatch(
      make_tukey_table(model_result()),
      error = function(e) data.frame(Error = e$message)
    )
    datatable(tbl, rownames = FALSE,
              options = list(pageLength = 20, dom = "tp", scrollX = TRUE)) |>
      formatStyle("p adjusted",
        backgroundColor = styleEqual("<0.001", "#d4edda"),
        color            = styleEqual("<0.001", "#155724"))
  })

  # ── Regression diagnostics ────────────────────────────────────────────────
  output$diagnostics_card <- renderUI({
    req(model_result())
    mt <- input$model_type
    if (!mt %in% c("lm_simple", "lm_multi", "logit_simple", "logit_multi")) return(NULL)
    label <- if (mt %in% c("lm_simple", "lm_multi"))
      "🔬 Linear Regression Diagnostics"
    else
      "🔬 Logistic Regression Diagnostics"

    tags$div(class = "card",
      tags$div(class = "card-header", label),
      tags$div(class = "card-body",
        uiOutput("diag_stats_ui"),
        hr(),
        plotOutput("diag_plots", height = "520px")
      )
    )
  })

  output$diag_stats_ui <- renderUI({
    req(model_result())
    mt  <- input$model_type
    mod <- model_result()

    if (mt %in% c("lm_simple", "lm_multi")) {
      # ── Linear regression key diagnostics ──────────────────────────────────
      s        <- summary(mod)
      res      <- residuals(mod)
      std_res  <- rstandard(mod)
      lev      <- hatvalues(mod)
      cooks    <- cooks.distance(mod)
      n        <- length(res)
      p_terms  <- length(coef(mod))
      lev_thr  <- 2 * p_terms / n
      cook_thr <- 4 / n

      n_high_lev   <- sum(lev > lev_thr, na.rm = TRUE)
      n_high_cook  <- sum(cooks > cook_thr, na.rm = TRUE)
      n_outlier    <- sum(abs(std_res) > 2, na.rm = TRUE)

      # Shapiro-Wilk on residuals (max 5000 obs)
      sw <- tryCatch({
        sw_samp <- if (n > 5000) sample(res, 5000) else res
        st <- shapiro.test(sw_samp)
        paste0("W = ", round(st$statistic, 4), ", p ", fmt_p(st$p.value),
               if (st$p.value > .05) " — residuals appear normal"
               else " — possible non-normality")
      }, error = function(e) "Could not compute (n may be too large)")

      # Breusch-Pagan-style: cor of |resid| with fitted
      bp_note <- tryCatch({
        ct <- cor.test(abs(res), fitted(mod))
        paste0("r(|resid|, fitted) = ", round(ct$estimate, 3),
               ", p ", fmt_p(ct$p.value),
               if (ct$p.value > .05) " — no strong evidence of heteroscedasticity"
               else " — possible heteroscedasticity")
      }, error = function(e) "Could not compute")

      div(class = "alert alert-info mb-0",
        tags$strong("Key diagnostic values"), br(), br(),
        tags$table(class = "table table-sm table-borderless mb-1",
          style = "font-size:0.88rem;",
          tags$tbody(
            tags$tr(
              tags$td(tags$b("R²")),
              tags$td(round(s$r.squared, 4)),
              tags$td(tags$b("Adj. R²")),
              tags$td(round(s$adj.r.squared, 4)),
              tags$td(tags$b("Residual SE")),
              tags$td(round(s$sigma, 4))
            ),
            tags$tr(
              tags$td(tags$b("Resid. range")),
              tags$td(paste0(round(min(res),3), " to ", round(max(res),3))),
              tags$td(tags$b("|Std.resid| > 2")),
              tags$td(paste0(n_outlier, " obs (", round(100*n_outlier/n,1), "%)")),
              tags$td(tags$b("F-statistic")),
              tags$td(paste0(round(s$fstatistic[1],3), " (p ",
                             fmt_p(pf(s$fstatistic[1], s$fstatistic[2],
                                      s$fstatistic[3], lower.tail=FALSE)), ")"))
            ),
            tags$tr(
              tags$td(tags$b("High leverage")),
              tags$td(paste0(n_high_lev, " obs (threshold: ", round(lev_thr,3), ")")),
              tags$td(tags$b("High Cook's D")),
              tags$td(paste0(n_high_cook, " obs (threshold: ", round(cook_thr,3), ")")),
              tags$td(tags$b("Max Cook's D")),
              tags$td(round(max(cooks, na.rm=TRUE), 4))
            ),
            tags$tr(
              tags$td(tags$b("Normality (S-W)")),
              tags$td(colspan = "3", sw),
              tags$td(tags$b("Homoscedasticity")),
              tags$td(colspan = "1", bp_note)
            )
          )
        )
      )

    } else if (mt %in% c("logit_simple", "logit_multi")) {
      # ── Logistic regression key diagnostics ────────────────────────────────
      null  <- tryCatch(update(mod, . ~ 1), error = function(e) NULL)
      dev   <- mod$deviance
      df_r  <- mod$df.residual
      mcf   <- if (!is.null(null)) round(1 - dev / null$deviance, 4) else NA
      aic   <- round(AIC(mod), 2)
      bic   <- round(BIC(mod), 2)
      n     <- nobs(mod)
      cooks <- cooks.distance(mod)
      lev   <- hatvalues(mod)
      p_t   <- length(coef(mod))
      lev_thr  <- 2 * p_t / n
      cook_thr <- 4 / n
      n_high_lev  <- sum(lev > lev_thr, na.rm = TRUE)
      n_high_cook <- sum(cooks > cook_thr, na.rm = TRUE)

      # Hosmer-Lemeshow
      hl_txt <- tryCatch({
        y_hat <- fitted(mod)
        y_obs <- mod$y
        cuts  <- unique(quantile(y_hat, probs=seq(0,1,length.out=11), na.rm=TRUE))
        if (length(cuts) < 3) stop("degenerate")
        grp   <- cut(y_hat, breaks=cuts, include.lowest=TRUE)
        obs1  <- tapply(y_obs, grp, sum)
        exp1  <- tapply(y_hat, grp, sum)
        ng    <- tapply(y_obs, grp, length)
        x2    <- sum((obs1-exp1)^2 / (exp1*(1-exp1/ng)), na.rm=TRUE)
        df_hl <- length(obs1) - 2
        p_hl  <- pchisq(x2, df=df_hl, lower.tail=FALSE)
        list(stat=paste0("χ²(", df_hl, ") = ", round(x2,3)),
             p   =fmt_p(p_hl),
             note=if(p_hl>.05) "good fit" else "possible misfit")
      }, error=function(e) list(stat="—", p="—", note="could not compute"))

      # Discrimination: concordance (c-stat / AUC approximation)
      c_stat <- tryCatch({
        y_hat <- fitted(mod); y_obs <- mod$y
        pairs <- outer(y_hat[y_obs==1], y_hat[y_obs==0], "-")
        conc  <- sum(pairs > 0, na.rm=TRUE)
        ties  <- sum(pairs == 0, na.rm=TRUE)
        total <- sum(y_obs==1) * sum(y_obs==0)
        round((conc + 0.5*ties) / total, 4)
      }, error=function(e) NA)

      div(class = "alert alert-info mb-0",
        tags$strong("Key diagnostic values"), br(), br(),
        tags$table(class = "table table-sm table-borderless mb-1",
          style = "font-size:0.88rem;",
          tags$tbody(
            tags$tr(
              tags$td(tags$b("AIC")),        tags$td(aic),
              tags$td(tags$b("BIC")),        tags$td(bic),
              tags$td(tags$b("McFadden R²")),tags$td(mcf)
            ),
            tags$tr(
              tags$td(tags$b("Residual deviance")),
              tags$td(paste0(round(dev,2), " on ", df_r, " df")),
              tags$td(tags$b("C-statistic (AUC)")),
              tags$td(if(!is.na(c_stat)) paste0(c_stat,
                if(c_stat>=.8) " — good discrimination"
                else if(c_stat>=.7) " — acceptable"
                else " — poor discrimination") else "—"),
              tags$td(tags$b("n")), tags$td(n)
            ),
            tags$tr(
              tags$td(tags$b("Hosmer-Lemeshow")),
              tags$td(hl_txt$stat),
              tags$td(tags$b("p")),
              tags$td(paste0(hl_txt$p, " — ", hl_txt$note)),
              tags$td(tags$b("High Cook's D")),
              tags$td(paste0(n_high_cook, " obs (threshold: ", round(cook_thr,3), ")"))
            ),
            tags$tr(
              tags$td(tags$b("High leverage")),
              tags$td(paste0(n_high_lev, " obs (threshold: ", round(lev_thr,3), ")")),
              tags$td(tags$b("Max Cook's D")),
              tags$td(round(max(cooks,na.rm=TRUE),4)),
              tags$td(tags$b("Max leverage")),
              tags$td(round(max(lev,na.rm=TRUE),4))
            )
          )
        )
      )
    }
  })

  output$diag_plots <- renderPlot({
    req(model_result())
    mt  <- input$model_type
    mod <- model_result()
    if (mt %in% c("lm_simple", "lm_multi")) {
      ps <- make_lm_diag_plots(mod)
    } else if (mt %in% c("logit_simple", "logit_multi")) {
      ps <- make_logit_diag_plots(mod)
    } else return(NULL)
    grid.arrange(ps$p1, ps$p2, ps$p3, ps$p4, ncol = 2)
  })

  # ── Explore plots ─────────────────────────────────────────────────────────
  output$plot_controls <- renderUI({
    nc <- num_cols(); ac <- all_cols()
    tagList(
      selectInput("plot_type", "Plot type:",
                  choices=c("Histogram","Box plot","Scatter plot",
                            "Bar chart","Density plot")),
      selectInput("plot_x", "X variable:", ac),
      conditionalPanel(
        "input.plot_type == 'Scatter plot' || input.plot_type == 'Box plot'",
        selectInput("plot_y", "Y variable:", nc)
      ),
      conditionalPanel(
        "input.plot_type != 'Histogram' && input.plot_type != 'Density plot'",
        selectInput("plot_fill", "Color by (optional):", c("None", ac))
      )
    )
  })

  output$main_plot <- renderPlot({
    input$make_plot
    isolate({
      req(input$plot_x)
      df  <- active_data(); pt <- input$plot_type
      xv  <- input$plot_x;  yv <- input$plot_y
      clr <- if (!is.null(input$plot_fill) &&
                 input$plot_fill != "None") input$plot_fill else NULL
      p <- ggplot(df, aes_string(x=xv)) +
        theme_minimal(base_size=14) +
        theme(plot.title=element_text(face="bold"))
      if (pt=="Histogram")
        p + geom_histogram(fill="#2c7bb6",color="white",bins=30,alpha=.8) +
          labs(title=paste("Histogram of",xv), y="Count")
      else if (pt=="Density plot")
        p + geom_density(fill="#2c7bb6",alpha=.5) +
          labs(title=paste("Density of",xv), y="Density")
      else if (pt=="Box plot") {
        req(yv)
        ggplot(df,aes_string(x=xv,y=yv,fill=if(!is.null(clr))clr else NULL)) +
          geom_boxplot(alpha=.7) + theme_minimal(base_size=14) +
          labs(title=paste("Box plot:",yv,"by",xv))
      } else if (pt=="Scatter plot") {
        req(yv)
        ggplot(df,aes_string(x=xv,y=yv,color=if(!is.null(clr))clr else NULL)) +
          geom_point(alpha=.6,size=2) +
          geom_smooth(method="lm",se=TRUE,color="#e74c3c") +
          theme_minimal(base_size=14) +
          labs(title=paste("Scatter:",xv,"vs",yv))
      } else {
        p + geom_bar(fill="#2c7bb6",alpha=.8) +
          labs(title=paste("Bar chart of",xv), y="Count")
      }
    })
  })
}

shinyApp(ui, server)
