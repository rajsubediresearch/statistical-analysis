library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(MASS)
library(nnet)
library(readxl)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !identical(a, "")) a else b

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
        selectInput("model_type", "Statistical Method:",
                    choices=c(
                      "Simple Linear Regression"        = "lm_simple",
                      "Multiple Linear Regression"      = "lm_multi",
                      "Simple Logistic Regression"      = "logit_simple",
                      "Multiple Logistic Regression"    = "logit_multi",
                      "Ordinal Logistic Regression"     = "ordinal",
                      "Multinomial Logistic Regression" = "multinom",
                      "Independent t-test"              = "ttest",
                      "One-way ANOVA"                   = "anova"
                    )),
        hr(),
        uiOutput("var_selectors"),
        hr(),
        checkboxInput("show_plot",
                      "📈 Show analysis plot below results",
                      value=FALSE),
        br(),
        actionButton("run_model", "▶ Run Analysis",
                     class="btn-primary w-100")
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
    out <- input$outcome
    exp <- input$exposure
    cov <- if (!is.null(input$covariates)) input$covariates else character(0)

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
    mt <- input$model_type; mod <- model_result()
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
