# =====================================================================
# Spatial Statistics Dashboard - COVID-19 & Community Health Centres, South Africa
# Methods: (1) Global Moran's I  (2) Local Moran's I (LISA)  (3) CSR tests on CHC point pattern
# =====================================================================
library(dplyr); library(ggplot2)
library(sf); library(leaflet); library(spdep); library(spatstat)
# Load these LAST: spatstat defines its own card(), which would mask bslib::card()
library(shiny); library(bslib)

# ---------------------------------------------------------------- DATA
needed <- file.path("data", c("ts.rds", "pop.rds", "prov.rds", "chc.rds"))
if (!all(file.exists(needed)))
  stop("Missing: ", paste(needed[!file.exists(needed)], collapse = ", "),
       ". Working directory is ", getwd(),
       ". Run 00_prepare_data.R from the app folder first.")
ts   <- readRDS("data/ts.rds")
pop  <- readRDS("data/pop.rds")
prov <- readRDS("data/prov.rds")          # sf, EPSG:4326, column 'province'
chc  <- readRDS("data/chc.rds")           # sf points, EPSG:4326

# Project to metres so neighbours, centroids and distances are valid
crs_sa <- "+proj=aea +lat_1=-24 +lat_2=-33 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m +no_defs"
prov_p <- st_transform(prov, crs_sa)
chc_p  <- st_transform(chc, crs_sa)
cent   <- st_coordinates(st_centroid(st_geometry(prov_p)))

chc_n <- st_drop_geometry(chc) %>% count(province = Province, name = "n_chc")

# Panel of indicators for every province x date
panel <- ts %>%
  left_join(pop, by = "province") %>%
  left_join(chc_n, by = "province") %>%
  mutate(cases_raw   = cases,
         cases_100k  = cases / pop * 1e5,
         deaths_100k = deaths / pop * 1e5,
         cfr         = ifelse(cases > 0, deaths / cases * 100, NA_real_),
         rec_ratio   = ifelse(cases > 0, recoveries / cases * 100, NA_real_),
         chc_100k    = n_chc / pop * 1e5)

indicators <- c("Cumulative cases per 100k" = "cases_100k",
                "Cumulative deaths per 100k" = "deaths_100k",
                "Case fatality ratio (%)" = "cfr",
                "Recovery ratio (%)" = "rec_ratio",
                "CHCs per 100k people (static)" = "chc_100k",
                "Cumulative cases, raw count (not population-adjusted)" = "cases_raw")

# ---- point pattern objects (units = km) -----------------------------
g_km <- st_geometry(prov_p) * 0.001
st_crs(g_km) <- NA
prov_km <- st_sf(province = prov_p$province, geometry = g_km)

make_owin <- function(sfc_km) {
  u <- st_union(sfc_km) %>% st_simplify(dTolerance = 3, preserveTopology = TRUE) %>% st_cast("POLYGON")
  a <- as.numeric(st_area(u)); u <- u[a >= 0.01 * max(a)]       # drop tiny islands
  orient <- function(m, ccw) {                                    # spatstat: exterior ccw, holes cw
    m <- m[-nrow(m), , drop = FALSE]
    s <- 0.5 * sum(m[, 1] * c(m[-1, 2], m[1, 2]) - c(m[-1, 1], m[1, 1]) * m[, 2])
    if ((s > 0) != ccw) m <- m[nrow(m):1, , drop = FALSE]
    list(x = m[, 1], y = m[, 2])
  }
  rings <- list()
  for (i in seq_along(u)) for (r in seq_along(u[[i]]))
    rings[[length(rings) + 1]] <- orient(u[[i]][[r]], ccw = (r == 1))
  tryCatch(owin(poly = rings), error = function(e) NULL)
}
win <- make_owin(g_km)
win_ok <- !is.null(win)
if (!win_ok) win <- owin(range(cent[, 1] / 1000) + c(-500, 500), range(cent[, 2] / 1000) + c(-500, 500))
xy <- st_coordinates(chc_p) / 1000
X  <- suppressWarnings(ppp(xy[, 1], xy[, 2], window = win))
n_dropped <- nrow(xy) - npoints(X)

# ------------------------------------------------------------ HELPERS
nb_edges <- function(nb, cent) {
  do.call(rbind, lapply(seq_along(nb), function(i) {
    j <- nb[[i]]; j <- j[j > 0]
    if (length(j)) data.frame(x = as.numeric(cent[i, 1]), y = as.numeric(cent[i, 2]), xend = as.numeric(cent[j, 1]), yend = as.numeric(cent[j, 2]))
  }))
}
moran_text <- function(I, EI, p, alpha) {
  if (p < alpha) {
    if (I > EI) paste0("p < ", alpha, ": reject H0. I is above E[I], so similar values tend to be neighbours (positive spatial autocorrelation / clustering).")
    else paste0("p < ", alpha, ": reject H0. I is below E[I], so neighbours tend to have dissimilar values (negative spatial autocorrelation).")
  } else paste0("p >= ", alpha, ": fail to reject H0. No statistical evidence against a spatially random pattern. With only 9 provinces the test has low power, so this is not proof of randomness.")
}
cluster_cols <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6", "High-Low" = "#fdae61",
                  "Low-High" = "#abd9e9", "Not significant" = "#d9d9d9")

latest <- panel[panel$Date == max(panel$Date), ]

abbr <- c("Eastern Cape" = "EC", "Free State" = "FS", "Gauteng" = "GP", "KwaZulu-Natal" = "KZN",
          "Limpopo" = "LP", "Mpumalanga" = "MP", "Northern Cape" = "NC", "North West" = "NW",
          "Western Cape" = "WC")
prov_cols <- setNames(c("#1f77b4", "#e07b39", "#2ca02c", "#d62728", "#9467bd",
                        "#8c564b", "#e377c2", "#17a2b8", "#b5a300"), sort(unique(panel$province)))

# Consistent ggplot look for every chart
theme_dash <- function(base = 13) {
  theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold", colour = "#1f4e79"),
          plot.subtitle = element_text(colour = "#5a6572"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e6ebf1"),
          axis.title = element_text(colour = "#3a4654"),
          legend.position = "bottom", legend.title = element_text(face = "bold"))
}

# Base maps shared by both leaflet maps
add_basemaps <- function(m) {
  # Key-free base maps (CARTO tiles now need an API key, so they are not used)
  m %>% addProviderTiles("Esri.WorldGrayCanvas", group = "Light grey") %>%
    addProviderTiles("OpenStreetMap", group = "Streets") %>%
    addProviderTiles("Esri.WorldTopoMap", group = "Topographic") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite")
}
base_groups <- c("Light grey", "Streets", "Topographic", "Satellite")

# Coloured HTML tables (heat-map cells, badges)
html_tbl <- function(df, style_fun = NULL) {
  hdr <- tags$tr(lapply(names(df), tags$th))
  rows <- lapply(seq_len(nrow(df)), function(i) tags$tr(lapply(seq_along(df), function(j) {
    st <- if (is.null(style_fun)) NULL else style_fun(i, j, df[[j]][i], names(df)[j])
    tags$td(as.character(df[[j]][i]), style = st)
  })))
  div(class = "table-wrap", tags$table(class = "dash-table", tags$thead(hdr), tags$tbody(rows)))
}

css <- "
body{background:#eef2f7}
.navbar{background:linear-gradient(90deg,#0f2f4f,#1f4e79) !important;box-shadow:0 2px 8px rgba(0,0,0,.25)}
.navbar-brand,.navbar .navbar-brand{color:#fff !important;font-weight:700;letter-spacing:.2px}
.navbar .nav-link{color:#cfe0f1 !important;font-weight:500}
.navbar .nav-link.active,.navbar .nav-link:hover{color:#fff !important;border-bottom:3px solid #e07b39}
.card{border:none !important;border-radius:14px !important;box-shadow:0 2px 10px rgba(15,47,79,.10)}
.card-header{background:#fff !important;font-weight:700;color:#1f4e79;border-bottom:3px solid #e07b39 !important;border-radius:14px 14px 0 0 !important}
.bslib-sidebar-layout > .sidebar{background:#fff !important;box-shadow:2px 0 8px rgba(15,47,79,.08)}
.hero{background-image:radial-gradient(circle at 12% 18%,rgba(255,255,255,.16),transparent 42%),
 radial-gradient(circle at 90% 85%,rgba(224,123,57,.45),transparent 38%),
 linear-gradient(120deg,#0f2f4f,#1f4e79 55%,#2e75b6);
 color:#fff;padding:2.6rem 2.2rem;border-radius:16px;margin-bottom:1.3rem;box-shadow:0 6px 18px rgba(15,47,79,.25)}
.hero h1{font-weight:800;margin-bottom:.5rem} .hero p{font-size:1.1rem;max-width:780px;opacity:.95}
.intro{background:#fff;border-left:6px solid #e07b39;padding:.8rem 1.1rem;border-radius:10px;margin-bottom:1rem;box-shadow:0 1px 6px rgba(15,47,79,.08)}
.method-card .card-header{font-size:1.05rem}
.step{background:#fff;border-radius:12px;padding:1rem 1.2rem;box-shadow:0 1px 6px rgba(15,47,79,.08);display:flex;gap:.9rem;align-items:center}
.step-n{background:#e07b39;color:#fff;font-weight:800;border-radius:50%;min-width:2.2rem;height:2.2rem;display:flex;align-items:center;justify-content:center}
.bslib-value-box.vb-blue{background:linear-gradient(135deg,#1f4e79,#2e75b6) !important;color:#fff !important}
.bslib-value-box.vb-red{background:linear-gradient(135deg,#a4262c,#d9534f) !important;color:#fff !important}
.bslib-value-box.vb-green{background:linear-gradient(135deg,#1b7f3b,#4cb86a) !important;color:#fff !important}
.bslib-value-box.vb-orange{background:linear-gradient(135deg,#d9642a,#f0a04b) !important;color:#fff !important}
.bslib-value-box.vb-teal{background:linear-gradient(135deg,#0e7c86,#35b6bf) !important;color:#fff !important}
.bslib-value-box.vb-grey{background:linear-gradient(135deg,#475569,#7b8da3) !important;color:#fff !important}
.bslib-value-box.vb-amber{background:linear-gradient(135deg,#b7791f,#ecc94b) !important;color:#222 !important}
.table-wrap{overflow-x:auto}
.dash-table{width:100%;border-collapse:separate;border-spacing:0;font-size:.9rem}
.dash-table th{background:#1f4e79;color:#fff;padding:.5rem .6rem;text-align:center;position:sticky;top:0}
.dash-table th:first-child{border-radius:8px 0 0 0;text-align:left} .dash-table th:last-child{border-radius:0 8px 0 0}
.dash-table td{padding:.4rem .6rem;border-bottom:1px solid #e6ebf1}
.dash-table tbody tr:nth-child(even){background:#f6f9fc} .dash-table tbody tr:hover{background:#e8f0f8}
.app-footer{text-align:center;color:#6b7785;font-size:.85rem;padding:1rem}
"
intro <- function(...) div(class = "intro", ...)
vbox <- function(title, value, cls, ...) value_box(title, value, ..., theme = "primary", class = cls)

# ----------------------------------------------------------------- UI
ui <- page_navbar(
  title = "Spatial Statistics: COVID-19 & Health Centres in South Africa",
  id = "nav", selected = "home", header = tags$style(HTML(css)),
  footer = div(class = "app-footer", "Spatial Statistics group project  |  Data: provincial COVID-19 counts, Community Health Centres, GADM boundaries, Stats SA population (verify)"),
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1f4e79", base_font = "Segoe UI, Helvetica, Arial, sans-serif"),
  sidebar = sidebar(id = "side", width = 320, open = FALSE,
    conditionalPanel("input.nav != 'pp'",
      h6("Lattice data (provinces)"),
      selectInput("ind", "Indicator", indicators),
      sliderInput("date", "Date (cumulative to)", min = min(ts$Date), max = max(ts$Date),
                  value = max(ts$Date), timeFormat = "%d %b %Y", animate = TRUE),
      hr(), h6("Spatial weights W"),
      selectInput("nbtype", "Neighbourhood", c("Queen contiguity" = "queen", "Rook contiguity" = "rook",
                                               "k nearest neighbours" = "knn")),
      conditionalPanel("input.nbtype == 'knn'", sliderInput("k", "k", 1, 5, 2, step = 1)),
      selectInput("wstyle", "Weights style", c("Row-standardised (W)" = "W", "Binary (B)" = "B")),
      selectInput("alt", "Alternative hypothesis", c("Two-sided" = "two.sided",
                  "Positive autocorrelation (greater)" = "greater", "Negative (less)" = "less")),
      sliderInput("alpha", "Significance level", 0.01, 0.10, 0.05, step = 0.01)
    ),
    conditionalPanel("input.nav == 'pp'",
      h6("Point pattern (CHC locations)"),
      radioButtons("ppview", "View", c("Intensity (kernel density)" = "kde",
                                       "Quadrat test" = "quad",
                                       "Distance functions (Monte Carlo)" = "fun")),
      conditionalPanel("input.ppview == 'kde'", sliderInput("bw", "Bandwidth sigma (km)", 20, 250, 80, step = 10),
                       checkboxInput("showpts", "Show points", TRUE)),
      conditionalPanel("input.ppview == 'quad'", sliderInput("nx", "Quadrats: columns", 2, 8, 4),
                       sliderInput("ny", "Quadrats: rows", 2, 8, 3)),
      conditionalPanel("input.ppview == 'fun'",
        selectInput("fn", "Function", c("K  (Ripley)" = "K", "L  (variance-stabilised K)" = "L",
                                        "G  (nearest-neighbour)" = "G", "F  (empty space)" = "F")),
        sliderInput("nsim", "CSR simulations", 19, 199, 39, step = 20),
        actionButton("run", "Run simulation", class = "btn-primary"))
    )
  ),
  nav_panel("Home", value = "home",
    div(class = "hero", h1("COVID-19 and Health Centres in South Africa"),
        p("An interactive spatial statistics dashboard. Explore where COVID-19 burden clusters across provinces, and whether Community Health Centres are spread randomly across the country."),
        actionButton("go_ov", "Explore the map", class = "btn-light btn-lg")),
    layout_columns(
      vbox("Confirmed cases", format(sum(latest$cases), big.mark = " "), "vb-blue", p(format(max(latest$Date), "%d %b %Y"))),
      vbox("Deaths", format(sum(latest$deaths), big.mark = " "), "vb-red"),
      vbox("Recoveries", format(sum(latest$recoveries), big.mark = " "), "vb-green"),
      vbox("Health centres", nrow(chc), "vb-orange", p("with coordinates"))),
    h4("How to use it", style = "margin-top:1.2rem;color:#1f4e79;font-weight:700"),
    layout_columns(col_widths = c(4, 4, 4),
      div(class = "step", span(class = "step-n", "1"), span("Pick an indicator and date in the sidebar")),
      div(class = "step", span(class = "step-n", "2"), span("Choose how provinces are neighbours (spatial weights W)")),
      div(class = "step", span(class = "step-n", "3"), span("Open a method tab and read the statistic"))),
    h4("Choose a spatial method", style = "margin-top:1.2rem;color:#1f4e79;font-weight:700"),
    layout_columns(col_widths = c(4, 4, 4),
      card(class = "method-card", card_header("1. Global Moran's I"),
           p("Are provinces with similar COVID-19 values also neighbours? One number for the whole country."),
           p(em("Data: provinces (lattice)")), actionButton("go_gm", "Open", class = "btn-primary")),
      card(class = "method-card", card_header("2. Local Moran's I (LISA)"),
           p("Which provinces form hot spots, cold spots or spatial outliers?"),
           p(em("Data: provinces (lattice)")), actionButton("go_lm", "Open", class = "btn-primary")),
      card(class = "method-card", card_header("3. Point pattern (CSR)"),
           p("Are health centre locations random, clustered or regular? Density maps, quadrat test, K/L/G/F functions."),
           p(em("Data: CHC locations (points)")), actionButton("go_pp", "Open", class = "btn-primary")))),
  nav_panel("Map explorer", value = "ov", uiOutput("vb"),
    layout_columns(col_widths = c(8, 4),
      card(card_header("Map"), checkboxInput("showchc", "Show Community Health Centres", FALSE),
           leafletOutput("map", height = 430)),
      card(card_header("Through time"), plotOutput("tsplot", height = 430)))),
  nav_panel("1. Global Moran's I", value = "gm",
    intro(strong("Question: "), "is the pattern of the chosen indicator across provinces spatially random? ", strong("H0: "), "I = E[I] (no spatial autocorrelation)."),
    layout_columns(col_widths = c(12), uiOutput("gm_vb")),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Moran scatterplot"), plotOutput("gm_scatter", height = 360)),
      card(card_header("Permutation test (Monte Carlo)"), sliderInput("mcsim", "Permutations", 99, 4999, 999, step = 100),
           plotOutput("gm_mc", height = 300))),
    navset_card_tab(
      nav_panel("Neighbour structure", plotOutput("nbplot", height = 340)),
      nav_panel("Weights matrix W", uiOutput("wtab")),
      nav_panel("Moran's I through time", plotOutput("gm_time", height = 300)))),
  nav_panel("2. Local Moran's I (LISA)", value = "lm",
    intro(strong("Question: "), "which individual provinces are significant clusters or outliers? ", strong("H0 (per province): "), "no local spatial association."),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Cluster map"), leafletOutput("lmap", height = 430)),
      card(card_header("Local statistics"), selectInput("padj", "p-value adjustment",
             c("none", "bonferroni", "fdr" = "BH")), uiOutput("ltab"))),
    card(card_header("How to read this"), uiOutput("ltxt"))),
  nav_panel("3. Point pattern (CSR)", value = "pp",
    intro(strong("Question: "), "are health centre locations consistent with Complete Spatial Randomness? ", strong("H0: "), "homogeneous Poisson process."), uiOutput("pp_head"), uiOutput("pp_ui")),
  nav_panel("About", value = "ab", card(uiOutput("about")))
)

# ------------------------------------------------------------- SERVER
server <- function(input, output, session) {

  observeEvent(input$go_ov, nav_select("nav", "ov"))
  observeEvent(input$go_gm, nav_select("nav", "gm"))
  observeEvent(input$go_lm, nav_select("nav", "lm"))
  observeEvent(input$go_pp, nav_select("nav", "pp"))
  observeEvent(input$nav, sidebar_toggle("side", open = !(input$nav %in% c("home", "ab"))))

  dat <- reactive({
    d <- panel[panel$Date == input$date, ]
    out <- prov_p %>% left_join(d, by = "province")
    out$value <- out[[input$ind]]
    out
  })
  v <- reactive({ x <- dat()$value; validate(need(!anyNA(x) && sd(x) > 0, "Indicator is missing or constant on this date.")); x })

  nb <- reactive({
    switch(input$nbtype,
      queen = poly2nb(prov_p, queen = TRUE, snap = 1000),
      rook  = poly2nb(prov_p, queen = FALSE, snap = 1000),
      knn   = knn2nb(knearneigh(cent, k = input$k)))
  })
  lw <- reactive(nb2listw(nb(), style = input$wstyle, zero.policy = TRUE))

  # ---------- Overview
  output$vb <- renderUI({
    d <- panel[panel$Date == input$date, ]
    layout_columns(
      vbox("Confirmed cases", format(sum(d$cases), big.mark = " "), "vb-blue"),
      vbox("Deaths", format(sum(d$deaths), big.mark = " "), "vb-red"),
      vbox("Recoveries", format(sum(d$recoveries), big.mark = " "), "vb-green"),
      vbox("Health centres mapped", nrow(chc), "vb-orange"))
  })
  output$map <- renderLeaflet({
    d <- st_transform(dat(), 4326); pal <- colorNumeric("YlOrRd", d$value)
    m <- leaflet(d) %>% add_basemaps() %>%
      addPolygons(fillColor = ~pal(value), fillOpacity = 0.65, color = "#ffffff", weight = 1.5, group = "Provinces",
                  highlightOptions = highlightOptions(weight = 3, color = "#1f4e79", bringToFront = TRUE),
                  label = ~sprintf("%s: %s", province, round(value, 2))) %>%
      addLegend(pal = pal, values = ~value, title = names(indicators)[indicators == input$ind], position = "bottomright")
    overlays <- "Provinces"
    if (isTRUE(input$showchc)) {
      m <- m %>% addCircleMarkers(data = chc, radius = 4, stroke = TRUE, color = "#ffffff", weight = 1,
                                  fillColor = "#1f4e79", fillOpacity = 0.85, label = ~Name, group = "Health centres")
      overlays <- c(overlays, "Health centres")
    }
    m %>% addLayersControl(baseGroups = base_groups, overlayGroups = overlays,
                           options = layersControlOptions(collapsed = TRUE))
  })
  output$tsplot <- renderPlot({
    lab <- names(indicators)[indicators == input$ind]
    now <- panel[panel$Date == input$date, ]
    ggplot(panel, aes(Date, .data[[input$ind]], colour = province)) +
      geom_line(linewidth = 1, alpha = 0.9) +
      geom_vline(xintercept = input$date, linetype = 2, colour = "#e07b39", linewidth = 0.8) +
      geom_point(data = now, size = 3.5) +
      scale_colour_manual(values = prov_cols, name = NULL) +
      scale_x_date(date_labels = "%b") +
      guides(colour = guide_legend(nrow = 3)) +
      labs(x = NULL, y = lab, title = "Province trends", subtitle = "Dashed line = selected date") + theme_dash()
  })

  # ---------- 1. Global Moran
  gm <- reactive(moran.test(v(), lw(), alternative = input$alt, zero.policy = TRUE))
  output$gm_vb <- renderUI({
    m <- gm(); I <- unname(m$estimate[1]); EI <- unname(m$estimate[2])
    tagList(
      layout_columns(
        vbox("Moran's I", sprintf("%.3f", I), "vb-blue"),
        vbox("E[I] = -1/(n-1)", sprintf("%.3f", EI), "vb-grey"),
        vbox("z-score", sprintf("%.2f", unname(m$statistic)), "vb-teal"),
        vbox("p-value", sprintf("%.4f", m$p.value), if (m$p.value < input$alpha) "vb-green" else "vb-amber")),
      p(strong("H0:"), " the spatial pattern is random (I = E[I]).  ", strong("Conclusion: "),
        moran_text(I, EI, m$p.value, input$alpha)))
  })
  output$gm_scatter <- renderPlot({
    z <- as.numeric(scale(v())); lz <- lag.listw(lw(), z, zero.policy = TRUE)
    quad <- ifelse(z > 0 & lz > 0, "High-High", ifelse(z < 0 & lz < 0, "Low-Low",
            ifelse(z > 0 & lz < 0, "High-Low", "Low-High")))
    df <- data.frame(z, lz, quad, province = dat()$province)
    L <- max(abs(c(z, lz))) * 1.2
    ggplot(df, aes(z, lz)) +
      annotate("rect", xmin = 0, xmax = L, ymin = 0, ymax = L, fill = cluster_cols["High-High"], alpha = 0.10) +
      annotate("rect", xmin = -L, xmax = 0, ymin = -L, ymax = 0, fill = cluster_cols["Low-Low"], alpha = 0.10) +
      annotate("rect", xmin = 0, xmax = L, ymin = -L, ymax = 0, fill = cluster_cols["High-Low"], alpha = 0.15) +
      annotate("rect", xmin = -L, xmax = 0, ymin = 0, ymax = L, fill = cluster_cols["Low-High"], alpha = 0.20) +
      annotate("text", x = L * 0.95, y = L * 0.95, label = "High-High", hjust = 1, fontface = "bold", colour = "#a4262c") +
      annotate("text", x = -L * 0.95, y = -L * 0.95, label = "Low-Low", hjust = 0, fontface = "bold", colour = "#1f5f99") +
      annotate("text", x = L * 0.95, y = -L * 0.95, label = "High-Low", hjust = 1, fontface = "bold", colour = "#b9661f") +
      annotate("text", x = -L * 0.95, y = L * 0.95, label = "Low-High", hjust = 0, fontface = "bold", colour = "#4a8db5") +
      geom_hline(yintercept = 0, colour = "grey55") + geom_vline(xintercept = 0, colour = "grey55") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#0f2f4f", linewidth = 1.1) +
      geom_point(aes(fill = quad), shape = 21, size = 4.5, colour = "white", stroke = 1.2) +
      geom_text(aes(label = abbr[province]), vjust = -1.2, size = 3.6, fontface = "bold", colour = "#1f4e79") +
      scale_fill_manual(values = cluster_cols, guide = "none") +
      coord_cartesian(xlim = c(-L, L), ylim = c(-L, L)) +
      labs(x = "Standardised value (z)", y = "Spatial lag of z (neighbour average)",
           title = "Moran scatterplot", subtitle = "Slope of the dark line = Moran's I when W is row-standardised") + theme_dash()
  })
  output$gm_mc <- renderPlot({
    set.seed(1)
    mc <- moran.mc(v(), lw(), nsim = input$mcsim, alternative = input$alt, zero.policy = TRUE)
    sims <- mc$res[seq_len(input$mcsim)]; obs <- unname(mc$statistic)
    ggplot(data.frame(sims), aes(sims)) +
      geom_histogram(bins = 30, fill = "#5b8fc4", colour = "white", alpha = 0.95) +
      geom_vline(xintercept = obs, colour = "#e07b39", linewidth = 1.5) +
      annotate("label", x = obs, y = Inf, vjust = 1.4, label = sprintf("Observed I = %.3f", obs),
               fill = "#e07b39", colour = "white", fontface = "bold") +
      labs(x = "Moran's I under random re-allocation of values", y = "Count",
           title = "Permutation test", subtitle = sprintf("Monte Carlo p = %.4f", mc$p.value)) + theme_dash()
  })
  output$nbplot <- renderPlot({
    e <- nb_edges(nb(), cent)
    nodes <- data.frame(x = cent[, 1], y = cent[, 2], lab = abbr[prov_p$province])
    ggplot() +
      geom_sf(data = dat(), aes(fill = value), colour = "white", linewidth = 0.8) +
      geom_segment(data = e, aes(x = x, y = y, xend = xend, yend = yend), colour = "#0f2f4f", linewidth = 0.9, alpha = 0.85) +
      geom_point(data = nodes, aes(x, y), shape = 21, size = 5, fill = "white", colour = "#0f2f4f", stroke = 1.3) +
      geom_text(data = nodes, aes(x, y, label = lab), size = 2.6, fontface = "bold", colour = "#0f2f4f") +
      scale_fill_distiller(palette = "YlOrRd", direction = 1, name = names(indicators)[indicators == input$ind]) +
      labs(title = "Neighbour structure", subtitle = "Lines join each province to its neighbours") +
      theme_void(base_size = 13) +
      theme(plot.title = element_text(face = "bold", colour = "#1f4e79"), legend.position = "bottom")
  })
  output$wtab <- renderUI({
    Wm <- listw2mat(lw()); nm <- prov_p$province
    df <- data.frame(Province = nm, stringsAsFactors = FALSE)
    for (k in seq_along(nm)) df[[abbr[nm[k]]]] <- ifelse(Wm[, k] == 0, ".", formatC(Wm[, k], format = "f", digits = 2))
    df[["Row sum"]] <- formatC(rowSums(Wm), format = "f", digits = 2)
    html_tbl(df, function(i, j, val, col) {
      if (j == 1 || col == "Row sum") return("font-weight:700;color:#1f4e79;")
      w <- suppressWarnings(as.numeric(val))
      if (is.na(w) || w == 0) return("color:#b8c2cc;text-align:center;")
      paste0("background:rgba(31,78,121,", round(0.25 + 0.75 * min(w, 1), 2), ");color:", if (w > 0.5) "#fff" else "#0f2f4f",
             ";text-align:center;font-weight:700;")
    })
  })
  output$gm_time <- renderPlot({
    l <- lw(); dts <- sort(unique(panel$Date))
    res <- lapply(dts, function(dt) {
      d <- panel[panel$Date == dt, ]; x <- d[[input$ind]][match(prov_p$province, d$province)]
      if (anyNA(x) || sd(x) == 0) return(NULL)
      m <- moran.test(x, l, alternative = input$alt, zero.policy = TRUE)
      data.frame(Date = dt, I = unname(m$estimate[1]), EI = unname(m$estimate[2]), p = m$p.value)
    })
    r <- bind_rows(res); validate(need(nrow(r) > 0, "No valid dates."))
    r$sig <- ifelse(r$p < input$alpha, "Significant", "Not significant")
    ggplot(r, aes(Date, I)) +
      geom_ribbon(aes(ymin = EI, ymax = I), fill = "#5b8fc4", alpha = 0.18) +
      geom_hline(aes(yintercept = EI), linetype = 2, colour = "grey45") +
      geom_vline(xintercept = input$date, colour = "#e07b39", linetype = 3, linewidth = 0.9) +
      geom_line(colour = "#1f4e79", linewidth = 1.1) + geom_point(aes(colour = sig), size = 2.2) +
      scale_colour_manual(values = c("Significant" = "#d7191c", "Not significant" = "#8a96a3"), name = NULL) +
      scale_x_date(date_labels = "%b") +
      labs(x = NULL, y = "Moran's I", title = "Spatial autocorrelation through time",
           subtitle = "Dashed = E[I]; orange dotted = selected date") + theme_dash()
  })

  # ---------- 2. Local Moran
  lisa <- reactive({
    x <- v(); l <- lw(); z <- as.numeric(scale(x)); lz <- lag.listw(l, z, zero.policy = TRUE)
    lm <- localmoran(x, l, alternative = "two.sided", zero.policy = TRUE)
    p <- p.adjust(lm[, 5], method = input$padj)
    quad <- ifelse(z > 0 & lz > 0, "High-High", ifelse(z < 0 & lz < 0, "Low-Low",
            ifelse(z > 0 & lz < 0, "High-Low", "Low-High")))
    cl <- ifelse(p < input$alpha, quad, "Not significant")
    data.frame(Province = prov_p$province, Value = x, Ii = lm[, 1], z = lm[, 4], p = p, Quadrant = quad, Cluster = cl)
  })
  output$lmap <- renderLeaflet({
    r <- lisa(); d <- st_transform(prov_p, 4326); d$Cluster <- factor(r$Cluster, levels = names(cluster_cols))
    pal <- colorFactor(cluster_cols, levels = names(cluster_cols))
    leaflet(d) %>% add_basemaps() %>%
      addPolygons(fillColor = ~pal(Cluster), fillOpacity = 0.7, color = "#ffffff", weight = 1.8, group = "LISA clusters",
                  highlightOptions = highlightOptions(weight = 3, color = "#1f4e79", bringToFront = TRUE),
                  label = ~sprintf("%s: %s", province, Cluster)) %>%
      addLegend(colors = unname(cluster_cols), labels = names(cluster_cols), title = "LISA cluster", position = "bottomright") %>%
      addLayersControl(baseGroups = base_groups, overlayGroups = "LISA clusters",
                       options = layersControlOptions(collapsed = TRUE))
  })
  output$ltab <- renderUI({
    r <- lisa(); a <- input$alpha
    df <- data.frame(Province = r$Province, Value = formatC(r$Value, format = "f", digits = 2),
                     `Local I` = formatC(r$Ii, format = "f", digits = 3), z = formatC(r$z, format = "f", digits = 2),
                     p = formatC(r$p, format = "f", digits = 4), Cluster = r$Cluster,
                     stringsAsFactors = FALSE, check.names = FALSE)
    html_tbl(df, function(i, j, val, col) {
      if (col == "Province") return("font-weight:700;color:#1f4e79;")
      if (col == "Cluster") {
        txt <- if (val %in% c("High-High", "Low-Low")) "#fff" else "#222"
        return(paste0("background:", cluster_cols[[val]], ";color:", txt, ";font-weight:700;text-align:center;"))
      }
      if (col == "Local I") return(paste0("font-weight:700;color:", if (as.numeric(val) > 0) "#d7191c" else "#2c7bb6", ";"))
      if (col == "p" && as.numeric(val) < a) return("font-weight:800;color:#1b7f3b;")
      NULL
    })
  })
  output$ltxt <- renderUI(tagList(
    p("Local Moran's I decomposes the global statistic into one value per province. A province is a ",
      strong("High-High"), " cluster if it and its neighbours are above average, ", strong("Low-Low"),
      " if both are below, and a spatial outlier (", strong("High-Low / Low-High"),
      ") if it differs from its neighbours. Only provinces with p below the chosen level are coloured."),
    p(strong("Caution: "), "with n = 9 each province has only 1-5 neighbours, so local p-values are very unstable, and 9 simultaneous tests inflate false positives. Compare p-value adjustments.")))

  # ---------- 3. Point pattern
  output$pp_head <- renderUI(tagList(
    p(sprintf("%d health centres with valid coordinates are analysed (%d fell outside the simplified window).", npoints(X), n_dropped)),
    if (!win_ok) p(strong("Warning: "), "province window could not be built; a bounding box is used instead.")))
  output$pp_ui <- renderUI({
    switch(input$ppview,
      kde = card(card_header("Kernel intensity estimate"), plotOutput("kde", height = 520),
                 p("Intensity lambda(s) = expected CHCs per unit area. Bandwidth sigma controls smoothing: small = bumpy, large = smooth.")),
      quad = layout_columns(col_widths = c(7, 5), card(card_header("Quadrat counts"), plotOutput("quadplot", height = 460)),
                            card(card_header("Chi-square test of CSR"), uiOutput("quadtxt"))),
      fun = layout_columns(col_widths = c(8, 4), card(card_header("Observed vs CSR envelope"), plotOutput("funplot", height = 460)),
                           card(card_header("How to read"), uiOutput("funtxt"))))
  })
  output$kde <- renderPlot({
    d <- as.data.frame(density(X, sigma = input$bw))
    p <- ggplot() + geom_raster(data = d, aes(x, y, fill = value * 1e4)) +
      geom_sf(data = prov_km, fill = NA, colour = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "plasma", name = "CHCs per 10 000 km^2") +
      labs(title = "Where health centres are concentrated") +
      theme_void(base_size = 13) + theme(plot.title = element_text(face = "bold", colour = "#1f4e79"), legend.position = "bottom", legend.key.width = unit(2, "cm"))
    if (isTRUE(input$showpts)) p <- p + geom_point(data = as.data.frame(X), aes(x, y), size = 0.6, colour = "white", alpha = 0.7)
    p
  })
  qt <- reactive(suppressWarnings(quadrat.test(X, nx = input$nx, ny = input$ny, method = "Chisq")))  # warning is shown in the UI text instead
  output$quadplot <- renderPlot(plot(qt(), main = ""))
  output$quadtxt <- renderUI({
    t <- qt(); cnt <- as.vector(quadratcount(X, nx = input$nx, ny = input$ny))
    tagList(
      p(strong("H0:"), " CHCs follow a homogeneous Poisson process (CSR)."),
      p(sprintf("X-squared = %.1f, df = %d, p = %.4g", unname(t$statistic), unname(t$parameter), t$p.value)),
      p(sprintf("Approx. variance-to-mean ratio of counts = %.2f (CSR is about 1; >1 clustered).", var(cnt) / mean(cnt))),
      p(strong("Note: "), "If the grid is fine, some quadrats have expected counts below 5 and the chi-square approximation is unreliable. Use fewer quadrats (smaller sliders)."),
      p(strong("Note: "), "With a fine grid, some quadrats have expected counts below 5 and the chi-square approximation is unreliable; use fewer quadrats."),
      p(em("Quadrats cut by the border have smaller area; the chi-square test accounts for this, the VMR does not. Expected counts below 5 weaken the test; try fewer quadrats."), p(em("Results depend on quadrat size (scale).")))
    )
  })
  fsim <- eventReactive(input$run, {
    withProgress(message = "Simulating CSR patterns...", {
      set.seed(42)
      switch(input$fn,
        K = envelope(X, Kest, nsim = input$nsim, correction = "border", verbose = FALSE),
        L = envelope(X, Lest, nsim = input$nsim, correction = "border", verbose = FALSE),
        G = envelope(X, Gest, nsim = input$nsim, correction = "km", verbose = FALSE),
        F = envelope(X, Fest, nsim = input$nsim, correction = "km", verbose = FALSE))
    })
  })
  output$funplot <- renderPlot({
    e <- fsim(); fn <- isolate(input$fn)
    if (fn == "L") plot(e, . - r ~ r, main = "L(r) - r", legend = FALSE, shadecol = "#cfe0f1") else plot(e, main = paste0(fn, "(r)"), legend = FALSE, shadecol = "#cfe0f1")
  })
  output$funtxt <- renderUI(tagList(
    p("Black = observed. Grey band = range of ", isolate(input$nsim), " simulated CSR patterns."),
    tags$ul(
      tags$li(strong("K / L: "), "above the band = more clustering than CSR at that distance; below = more regular."),
      tags$li(strong("G: "), "above the band = shorter nearest-neighbour distances (clustering)."),
      tags$li(strong("F: "), "reversed: below the band = large empty gaps (clustering).")),
    p(strong("Caveat: "), "CHCs are sited where people live, so the process is almost certainly not homogeneous. Rejection of CSR here mostly reflects varying intensity (first-order), not attraction between centres (second-order). Distances are in km; edge effects handled by border/Kaplan-Meier corrections."),
    p(em("Click 'Run simulation' in the sidebar after changing settings."))))

  output$about <- renderUI(tagList(
    h4("About this dashboard"),
    p("Data: provincial cumulative COVID-19 cases, deaths and recoveries (5 Mar - 4 Oct 2020) and 287 Community Health Centre locations (one without coordinates)."),
    p("Population denominators are Stats SA 2020 mid-year estimates entered manually (data/population_2020.csv) and should be verified."),
    p("Boundaries: GADM level 1. Methods: Global and Local Moran's I (spdep), CSR tests, kernel intensity and K/L/G/F functions (spatstat)."),
    p(strong("Limitations: "), "9 provinces is a very small lattice (low power, MAUP, ecological fallacy); the CHC list may be incomplete; reporting and testing differ by province.")))
}

shinyApp(ui, server)
