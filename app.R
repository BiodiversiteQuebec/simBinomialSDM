# Load necessary libraries
library(shiny)
library(geoR)
library(terra)
library(sf)
library(INLA)
library(inlabru)
library(mgcv)
library(shinyjs)
library(shinycssloaders)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 8180)

# Define UI
ui <- fluidPage(
  tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML("
      h2 {
        font-size: 26px;
      }
      h3 {
        font-size: 20px;
      }
      .shiny-input-container {
        font-size: 12px;
      }
      .form-control.shiny-bound-input, .selectize-input {
              height: 2.5vh;
              width: 50%;
      }"
      ))
  ),
  titlePanel("Simulate a binomial data surface from predictors and effort and find back with a purely spatial model"),
  sidebarLayout(
    sidebarPanel(
      h3("Simulation Grid"),
      div(numericInput("epsg", "EPSG:", value = 6623), style = "display: none;"),
      numericInput("simres", "Resolution for simulations\n(smaller values take longer)", value = 250),
      h3("Aggregation Grid"),
      numericInput("aggres", "Resolution for data aggregation\n(combine observation in cells)", value = 300),
      h3("INLA"),
      numericInput("pedge", "Proportion of mesh edge as\na fraction of the study region", value = 0.02),
      numericInput("prior_range", "Prior range \n Pr(range < prior) = prob", value = 1000),
      numericInput("prior_sd", "Prior sd \n Pr(sd > prior) = prob", value = 1),
      numericInput("prior_prob", "Prior prob", value = 0.1),
      h3("GAM"),
      numericInput("k", "Allowed spline complexity\n(higher values take longer to run)", value = 20),
      h3("__________________"),
      actionButton("run", "RUN MODELS"),
      width = 2
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabset",
        #tags$head(tags$style(".shiny-plot-output{height:80vh !important;}")),
        tabPanel("Predictors", 
                 fluidRow(
                   column(4, offset = 0, align = "left", 
                          h3("Gradient Predictor"),
                          numericInput("b1", "Y Gradient (b1):", value = -0.00000006)
                   ),
                   column(4, offset = 0, align = "left", 
                          h3("Continuous Predictor"),
                          numericInput("b3", "Coefficient (b3):", value = 1),
                          numericInput("vpred", "Variance of field (vpred)", value = 1),
                          numericInput("rpred", "Range of field (rpred)", value = 5000)
                   ),
                   column(4, offset = 0, align = "left", 
                          h3("Binary Predictor"),
                          numericInput("b2", "Coefficient (b2):", value = 2),
                          numericInput("vbin", "Variance of field (vbin):", value = 1),
                          numericInput("rbin", "Range of field (rbin):", value = 1000),
                          numericInput("th", "Quantile for binarization (th):", value = 0.5)
                   )
                 ),
                 fluidRow(
                   column(4, offset = 0, align = "center", 
                          actionButton("update_gradient", "UPDATE")
                   ),
                   column(4, offset = 0, align = "center", 
                          actionButton("update_gaussian", "UPDATE")
                   ),
                   column(4, offset = 0, align = "center", 
                          actionButton("update_binary", "UPDATE")
                   )
                 ),
                 fluidRow(
                   column(4, offset = 0, 
                          plotOutput("plot_gradient") 
                   ),
                   column(4, offset = 0, 
                          plotOutput("plot_gaussian")
                   ),
                   column(4, offset = 0,
                          plotOutput("plot_binary")
                   )
                 )#,
                 #tableOutput('show_inputs')
        ),
        tabPanel("Effort and Observations", 
                 fluidRow(
                   column(4, offset = 0, align = "left", 
                          h3("Effort field and gradient"),
                          numericInput("beff", "Y Gradient (beff):", value = -0.0006),
                          numericInput("meff", "Mean effort (meff):", value = 2),
                          numericInput("veff", "Variance of field (veff):", value = 2),
                          numericInput("reff", "Range of field (reff):", value = 200)
                          
                   ),
                   column(4, offset = 0, align = "left", 
                          h3("Sprinkled low effort"),
                          numericInput("trunc", "Remove low effort below (trunc):", value = 10),
                          numericInput("ns", "Number of low effort locations (ns):", value = 20),
                          numericInput("nmax", "Max low effort sampled in 1:nmax:", value = 5),
                          
                   ),
                   column(4, offset = 0, align = "left", 
                          h3("Model Parameters"),
                          numericInput("b0", "Intercept (b0):", value = -1)#,
                          #actionButton("update_observations", "Update")
                   )
                 ),
                 fluidRow(
                   column(12, offset = 0, align = "center", 
                          actionButton("update_effort", "UPDATE")
                   )
                 ),
                 fluidRow(
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_effort")
                   ),
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_aggeffort")
                   ),
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_obs")
                   ),
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_raw")
                   ) 
                 )
        ),
        tabPanel("Truth and Predictions", 
                 fluidRow(
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_truth"),
                          plotOutput("plot_aggeffort2")
                          
                   ),
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_preds"),
                          plotOutput("plot_sd")
                          
                   ),
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_gampreds"),
                          plotOutput("plot_low")
                   ),
                   column(3, offset = 0, align = "center", 
                          plotOutput("plot_raw2"),
                          plotOutput("plot_high")
                   )
                 )
        )
      ),
      width = 10
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  logit <- function(x){
    log(x / (1 - x))
  }
  
  invlogit <- function(x){
    exp(x) / (exp(x) + 1)
  }
  
  clear_preds <- function(){
    output$plot_truth <- NULL
    output$plot_preds <- NULL
    output$plot_gampreds <- NULL
    output$plot_sd <- NULL
    output$plot_low <- NULL
    output$plot_high <- NULL
  }
  
  clear_obs <- function(){
    output$plot_obs <- output$plot_obs2 <- NULL
    output$plot_raw <- output$plot_raw2 <- NULL
    output$plot_effort <- output$plot_effort2 <- NULL
    output$plot_aggeffort <- output$plot_aggeffort2 <- NULL
  }
  
  clear_predictors <- function(){
    output$plot_gaussian <- NULL
    output$plot_gradient <- NULL
    output$plot_binary <- NULL
  }
  
  
  rv <- reactiveValues(
          gaussian = NULL,
          binary = NULL,
          gradient = NULL,
          dat = NULL,
          xy = NULL,
          rr = NULL,
          truth = NULL
        )
  
  
  
  params <- reactive({
    x <- reactiveValuesToList(input)
    
    r <- rast(resolution = input$simres, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = paste0("epsg:",input$epsg))
    grid <- xyFromCell(r, 1:ncell(r))
    
    region <- ext(r) |> vect() |> st_as_sf()
    st_crs(region) <- input$epsg
    
    rv$rr <- disagg(r, 5)
    
    agg <- rast(resolution = input$aggres, ext = ext(region), crs = crs(region))

    x <- c(x, r = r, list(grid = grid, region = region, agg = agg, gaussian = NULL, gradient = NULL, binary = NULL))
    x$mar <- c(2,1,4,4)
    x$textpos <- matrix(c(c(ext(r)[1] + diff(ext(r)[1:2])/2), ext(r)[4]), ncol = 2) * c(1.02, 1.02) 
    x$textadj <- c(0.5, 0)
    
    edge <- min(abs(c(diff(st_bbox(region)[c(3,1)]),diff(st_bbox(region)[c(4,2)]))))
    #pedge <- 0.02
    edge <- edge * x$pedge
    
    mesh <- fm_mesh_2d_inla(
      boundary = region, max.edge = c(edge, 3 * edge),
      cutoff = edge, offset = c(edge, 3 * edge),
      crs = fm_crs(region)
    )
    smesh <- mesh |> fm_as_sfc() |> st_as_sf()
    
    x$mesh <- mesh
    x$smesh <- smesh

    x

  })
  
  #output$show_inputs <- renderText({
  #  input$grid
  #})
  
  
  observeEvent(ignoreInit = TRUE, list(
    input$simres,
    input$epsg
  ),
  {
    clear_preds()
    clear_obs()
    clear_predictors()
  })
  
  observeEvent(ignoreInit = TRUE, list(
    input$aggres
  ),
  {
    clear_preds()
    clear_obs()
  })
  
  
  
  observeEvent(input$update_gaussian, {
    
    clear_obs()
    clear_preds()
    
    epsg <- params()$epsg
    simres <- params()$simres
    b3 <- params()$b3
    vpred <- params()$vpred
    rpred <- params()$rpred
    grid <- params()$grid
    r <- params()$r
    mar <- params()$mar
    textpos <- params()$textpos
    textadj <- params()$textadj

    gf1 <- grf(grid = grid, cov.pars = c(vpred, rpred)) 
    gaussian <- data.frame(gf1$coords, val = gf1$data) |> 
      rast(type = "xyz") |>
      disagg(fact = 5, method = "bilinear")
    crs(gaussian) <- "epsg:6623"
    
    rv$gaussian <- gaussian
    
    output$plot_gaussian <- renderPlot({
      plot(invlogit(b3 * gaussian), mar = mar)
      text(textpos, label = "Gaussian Field Predictor (GF)\n(from prob = 0.5)", xpd = TRUE, adj = textadj)
    })
    
    
  })
  
  
  observeEvent(input$update_gradient, {
    
    clear_obs()
    clear_preds()
    
    epsg <- params()$epsg
    simres <- params()$simres
    grid <- params()$grid
    r <- params()$r
    b1 <- params()$b1
    mar <- params()$mar
    textpos <- params()$textpos
    textadj <- params()$textadj
    
    gradient <- invlogit(setValues(rv$rr, crds(rv$rr)[, 2])^2 * b1)
    
    rv$gradient <- gradient
    
    output$plot_gradient <- renderPlot({
      plot(gradient, mar = mar)
      text(textpos, label = "Gradient Predictor\n(from prob = 0.5)", xpd = TRUE, adj = textadj)
    })
    
  })
  
  
  observeEvent(input$update_binary, {
    
    clear_obs()
    clear_preds()
    
    epsg <- params()$epsg
    simres <- params()$simres
    grid <- params()$grid
    r <- params()$r
    b2 <- params()$b2
    vbin <- params()$vbin
    rbin <- params()$rbin
    th <- params()$th
    mar <- params()$mar
    textpos <- params()$textpos
    textadj <- params()$textadj
    
    gf1 <- grf(grid = grid, cov.pars = c(vbin, rbin))
    binary <- data.frame(gf1$coords,val = gf1$data) |> 
      rast(type = "xyz") |>
      disagg(fact = 5, method = "bilinear")
    binary <- ifel(binary < global(binary, fun = function(i){quantile(i, th)})[1, 1], 1, 0)  
    crs(binary) <- "epsg:6623"
    
    rv$binary <- binary
    
    output$plot_binary <- renderPlot({
      plot(invlogit((b2 * binary) - (b2 / 2)), mar = mar)
      text(textpos, label = "Binary Predictor\n(from prob = 0.5)", xpd = TRUE, adj = textadj)
    })
    
  })
  
  observeEvent(input$update_effort, {
    
    showPageSpinner()
    
    clear_preds()
    
    epsg <- params()$epsg
    simres <- params()$simres
    grid <- params()$grid
    r <- params()$r
    agg <- params()$agg
    aggres <- params()$aggres
    mar <- params()$mar
    textpos <- params()$textpos
    textadj <- params()$textadj
    b0 <- params()$b0
    b1 <- params()$b1
    b2 <- params()$b2
    b3 <- params()$b3
    beff <- params()$beff
    meff <- params()$meff
    veff <- params()$veff
    reff <- params()$reff
    trunc <- params()$trunc
    ns <- params()$ns
    nmax <- params()$nmax
    smesh <- params()$smesh
    
    #rr <- disagg(r, 5)
    rr <- rv$rr# <- rr
    
    gf2 <- grf(grid = grid, cov.pars = c(veff, reff))
    effort <- data.frame(gf2$coords, val = gf2$data + (beff * gf2$coords[, 2]) + meff) |> 
      rast(type = "xyz") |>
      disagg(fact = 5, method = "bilinear") |>
      exp() |>
      round()
    effort[effort < trunc] <- 0
    sprinkle <- setValues(effort, sample(c(rep(0, ncell(effort) - ns), sample(1:nmax, size = ns, replace = TRUE)))) 
    sprinklewhere <- data.frame(xyFromCell(sprinkle, 1:ncell(sprinkle)), n = values(sprinkle)[, 1]) |> st_as_sf(coords = c("x", "y"), crs = epsg)
    sprinklewhere <- sprinklewhere[sprinklewhere$n > 0, ]
    effort <- effort + sprinkle
  #  
    
    ### Models ########################################
    lp <- b0 + 
      b1 * xyFromCell(rr, 1:ncell(rr))[, 2]^2 + 
      b3 * values(rv$gaussian)[, 1] +
      b2 * values(rv$binary)[, 1]
    p <- invlogit(lp)
    truth <- setValues(rr, p)
    rv$truth <- truth
    
    nobs <- rbinom(ncell(rr), size = values(effort)[, 1], prob = values(truth)[, 1])
    props <- nobs / values(effort)[, 1]
    raw <- setValues(truth, props)
    obs <- setValues(rr, nobs)
    obs <- data.frame(xyFromCell(obs, 1:ncell(obs)), n = nobs)
    obs$effort <- values(effort)[, 1]
    obs <- st_as_sf(obs, coords = c("x", "y"))
    dat <- obs
    dat <- cbind(dat, as.data.frame(st_coordinates(dat)))
    st_crs(dat) <- epsg
    obs <- obs[obs$n > 0, ]
    
    r1 <- rasterize(dat["n"], agg, field = "n", fun = sum)
    r2 <- rasterize(dat["effort"], agg, field = "effort", fun = sum)
    xy <- data.frame(xyFromCell(r1, 1:ncell(r1)), n = values(r1)[, 1]) |> 
      cbind(effort = values(r2)[, 1]) |>
      st_as_sf(coords = c("x", "y"), crs = epsg)
    xy <- cbind(xy, as.data.frame(st_coordinates(xy)))
    
    rv$dat <- dat
    rv$xy <- xy
    
    output$plot_obs <- output$plot_obs2 <- renderPlot({
      plot(raw, col = "white", legend = TRUE, mar = mar)
      par(mar = mar)
      text(textpos, label = "Observations", xpd = TRUE, adj = textadj)
      plot(smesh, border = adjustcolor("black", 0.05), add = TRUE)
      plot(st_geometry(obs), pch = 16, cex = 0.75, col = adjustcolor("black", 0.25), add = TRUE)
    })
    
    output$plot_raw <- output$plot_raw2 <- renderPlot({
      plot(r1 / r2, mar = mar)
      text(textpos, label = "Aggregated Raw Proportions Modeled", xpd = TRUE, adj = textadj)
    })
    
    
    output$plot_effort <- renderPlot({
      eff <- effort
      eff[eff == 0] <- NA
      plot(eff, mar = mar)
      text(textpos, label = paste("Raw Effort (", global(effort, "sum")[1, 1], "checklists)"), xpd = TRUE, adj = textadj)
    })
    
    output$plot_aggeffort <- output$plot_aggeffort2 <- renderPlot({
      r3 <- r2
      r3[r3 < 1] <- NA
      plot(r3, mar = mar)
      text(textpos, label = "Aggregated Effort\n(with low effort as numbers)", xpd = TRUE, adj = textadj)
      par(mar = mar)
      text(st_coordinates(sprinklewhere), label = sprinklewhere$n, cex = 0.6, col = adjustcolor("black", 0.5))
    })
    
    hidePageSpinner()
    
  })
  
  
  observeEvent(input$run, {
    
    showPageSpinner()
    
    updateTabsetPanel(
      session, "tabset",
      selected = "Truth and Predictions"
    )
    
    epsg <- params()$epsg
    simres <- params()$simres
    aggres <- params()$aggres
    b0 <- params()$b0
    b1 <- params()$b1
    b3 <- params()$b3
    vpred <- params()$vpred
    rpred <- params()$rpred
    pedge <- params()$pedge
    region <- params()$region
    mar <- params()$mar
    textpos <- params()$textpos
    textadj <- params()$textadj
    k <- params()$k
    prior_range <- params()$prior_range
    prior_sd <- params()$prior_sd
    prior_prob <- params()$prior_prob
    mesh <- params()$mesh
    
    matern <- inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(prior_sd, prior_prob),
                                  prior.range = c(prior_range, prior_prob)
    )
    
    comps <- ~ Intercept(1) + field(geometry, model = matern) 
    
    bru_options(bru_verbose = TRUE)
    
    fit <- bru(
      comps,
      like(
        family = "binomial", 
        data = rv$xy,
        formula = n ~ Intercept + field,
        E = NULL,
        weights = NULL,
        Ntrials = rv$xy$effort,
        options = list(control.inla = list(int.strategy = "eb"), control.predictor = list(link = 1))
      )
    )
    
    predictions <- predict(
      fit, rv$dat,
      ~ Intercept + field
    )
    pred <- predictions
    pred$mean <- invlogit(pred$mean)
    preds <- rasterize(pred["mean"], rv$rr, field = "mean", fun = mean)
    preds <- mask(preds, vect(region))
    
    unc <- rasterize(predictions["sd"], rv$rr, field = "sd", fun = mean)
    unc <- mask(unc, vect(region))
    
    low <- invlogit(rasterize(predictions["q0.025"], rv$rr, field = "q0.025", fun = mean))
    low <- mask(low, vect(region))
    
    high <- invlogit(rasterize(predictions["q0.975"], rv$rr, field = "q0.975", fun = mean))
    high <- mask(high, vect(region))
    
    gm <- gam(cbind(n, effort) ~ te(X, Y, bs = "ds", k = k, m = c(1, 0.5)), data = rv$xy, family = binomial(link = "logit"), method = "REML")
    p <- predict(gm, newdata = rv$dat, type = "response")
    gampreds <- setValues(preds, as.vector(p))
    
    ### Plot Outputs ###
    
    zlim <- range(c(values(rv$truth)[, 1], values(preds)[, 1], values(gampreds)[, 1]), na.rm = TRUE)
    
    output$plot_truth <- renderPlot({
      plot(rv$truth, mar = mar, range = zlim)
      text(textpos, label = "Truth\n(prob = GF + binary + gradient)", xpd = TRUE, adj = textadj)
    })
    
    output$plot_preds <- renderPlot({
      plot(preds, mar = mar, range = zlim)
      text(textpos, label = "Predictions\n(from a purely spatial model)", xpd = TRUE, adj = textadj)
    })
    
    output$plot_gampreds<- renderPlot({
      plot(gampreds, mar = mar, range = zlim)
      text(textpos, label = "Predictions\n(from a GAM)", xpd = TRUE, adj = textadj)
    })
    
    output$plot_sd <- renderPlot({
      plot(unc, mar = mar)
      text(textpos, label = "Predictions sd\n(from the purely spatial model)", xpd = TRUE, adj = textadj)
    })
    
    output$plot_low <- renderPlot({
      plot(low, mar = mar)
      text(textpos, label = "low CI\n(from the purely spatial model)", xpd = TRUE, adj = textadj)
    })
    
    output$plot_high <- renderPlot({
      plot(high, mar = mar)
      text(textpos, label = "high CI\n(from the purely spatial model)", xpd = TRUE, adj = textadj)
    })
    
    hidePageSpinner()
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
