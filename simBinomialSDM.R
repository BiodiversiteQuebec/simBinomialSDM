
library(geoR)
library(terra)
library(sf)
library(INLA)
library(inlabru)
library(mgcv)

logit <- function(x){
  log(p / (1 - p))
}
  
invlogit <- function(x){
  exp(x) / (exp(x) + 1)
}  

### Simulation grid
epsg = 6623
nc = 2000 # number of grid cells

### Data aggregation grid
res = 300 # resolution

### Model
b0 = -1 # intercept
b1 = -0.00000006 # y gradient

### Continuous habitat model
b3 = 1 # coef
vpred = 0.005 # variance of field
rpred = 5000 # range of field

### Binary habitat model
b2 = 0.005 # coef
vbin = 1 # variance of field
rbin = 500 # range of field
th = 0.75 # quantile for binarization

### Effort model
beff = -0.0006 # y gradient
meff = 2 # mean effort
veff = 2 # variance of field
reff = 200 # range of field
trunc = 10 # remove low effort below trunc
ns = 60 # nb of low effort locations
nmax = 10 # max low effort sampled in 1:nmax



#r <- rast(resolution = 200, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000)
#grid <- xyFromCell(r, 1:ncell(r))

#gf1 <- grf(grid = grid, cov.pars = c(vpred, rpred))


### Gaussian field predictors #################
gf1 <- grf(nc, grid = "reg", cov.pars = c(vpred, rpred), xlims = c(0, 10000), ylims = c(0, 10000)) 
#image(gf1)
predictors <- data.frame(gf1$coords,val = gf1$data) |> 
     rast(type = "xyz") |>
     disagg(fact = 5, method = "bilinear")
crs(predictors) <- "epsg:6623"


### Binary predictors #########################
gf1 <- grf(nc, grid = "reg", cov.pars = c(vbin, rbin), xlims = c(0, 10000), ylims = c(0, 10000))
binary <- data.frame(gf1$coords,val = gf1$data) |> 
  rast(type = "xyz") |>
  disagg(fact = 5, method = "bilinear")
binary <- ifel(binary < global(binary, fun = function(i){quantile(i, th)})[1, 1], 0, 1)  
crs(binary) <- "epsg:6623"


### Gradient predictor ########################
gradient <- invlogit(setValues(predictors, crds(predictors)[, 2])^2 * b1)


### Effort field ##############################
gf2 <- grf(nc, grid = "reg", cov.pars = c(veff, reff), xlims = c(0, 10000), ylims = c(0, 10000))
#image(gf2)
effort <- data.frame(gf2$coords,val = gf2$data + (beff * gf2$coords[, 2])+ meff) |> 
  rast(type = "xyz") |>
  disagg(fact = 5, method = "bilinear") |>
  exp() |>
  round()
# sprinkle low effort
effort[effort < trunc] <- 0
sprinkle <- setValues(effort, sample(c(rep(0, ncell(effort) - ns), sample(1:nmax, size = ns, replace = TRUE))))  
sprinklewhere <- data.frame(xyFromCell(sprinkle, 1:ncell(sprinkle)), n = values(sprinkle)[, 1]) |> st_as_sf(coords = c("x", "y"), crs = epsg)
sprinklewhere <- sprinklewhere[sprinklewhere$n > 0, ]
sprinklewhere
effort <- effort + sprinkle


#plot(effort)
#rev(sort(table(values(effort)[, 1])))


lp <- b0 + 
      b1 * xyFromCell(predictors, 1:ncell(predictors))[, 2]^2 + 
      b3 * values(predictors)[, 1] +
      b2 * values(binary)[, 1]
p <- invlogit(lp)
truth <- setValues(predictors, p)
#plot(truth)

nobs <- rbinom(ncell(predictors), size = values(effort)[, 1] , prob = values(truth)[, 1])
props <- nobs / values(effort)[, 1]

raw <- setValues(truth, props)
#plot(raw)

obs <- setValues(predictors, nobs)
obs <- data.frame(xyFromCell(obs, 1:ncell(obs)), n = nobs)
obs$effort <- values(effort)[, 1]
obs <- st_as_sf(obs, coords = c("x", "y"))
dat <- obs
dat <- cbind(dat, as.data.frame(st_coordinates(dat)))
st_crs(dat) <- epsg
obs <- obs[obs$n > 0, ]


#######################################
### model #############################
#######################################
region <- ext(predictors) |> vect() |> st_as_sf()
st_crs(region) <- epsg

edge <- min(abs(c(diff(st_bbox(region)[c(3,1)]),diff(st_bbox(region)[c(4,2)]))))
pedge <- 0.02
edge <- edge * pedge

mesh <- fm_mesh_2d_inla(
  boundary = region, max.edge = c(edge, 3 * edge), # km inside and outside
  cutoff = edge, offset = c(edge, 3 * edge),
  crs = fm_crs(region)
) # cutoff is min edge
smesh <- mesh |> fm_as_sfc() |> st_as_sf()


matern <- inla.spde2.pcmatern(mesh,
                      prior.sigma = c(1, 0.1),
                      prior.range = c(1000, 0.1)
)

comps <- ~ Intercept(1) + field(geometry, model = matern) 

# aggregate data or not

r <- rast(resolution = res, ext = ext(region))
#r <- predictors
r1 <- rasterize(dat["n"], r, field = "n", fun = sum)
r2 <- rasterize(dat["effort"], r, field = "effort", fun = sum)
xy <- data.frame(xyFromCell(r1, 1:ncell(r1)), n = values(r1)[, 1]) |> 
        cbind(effort = values(r2)[, 1]) |>
        st_as_sf(coords = c("x", "y"), crs = epsg)
xy <- cbind(xy, as.data.frame(st_coordinates(xy)))
#xy <- dat


bru_options(bru_verbose = TRUE)

fit <- bru(
  comps,
  like(
    family = "binomial", 
    data = xy,
    formula = n ~ Intercept +
      #latitude + latitude2 +
      #temp + temp2 +
      field,
    E = NULL,
    weights = NULL,#s2$pres/s2$counts, #ifelse(s2$counts>250,250,s2$counts)/250,
    Ntrials = xy$effort,
    options = list(control.inla = list(int.strategy = "eb"), control.predictor = list(link = 1))
  )
)

predictions <- predict(
  fit, dat,
  ~ Intercept + 
    #latitude + latitude2 + 
    #temp + temp2 +
    field
)
pred <- predictions
pred$mean <- invlogit(pred$mean)
preds <- rasterize(pred["mean"], predictors, field = "mean", fun = mean)
preds <- mask(preds, vect(region))

unc <- rasterize(predictions["sd"], predictors, field = "sd", fun = mean)
unc <- mask(unc, vect(region))


gm <- gam(cbind(n, effort) ~ te(X, Y, bs = "ds", k = 15, m = c(1, 0.5)), data = xy, family = binomial(link = "logit"), method = "REML")
p <- predict(gm, newdata = dat, type = "response")
gampreds <- setValues(preds, as.vector(p))


#####################################
### Plots ###########################
#####################################

png("simGaussianFields.png", width = 12, height = 9, units = "in", res = 500)

par(mfrow = c(3, 4))
mar <- c(1,1,2,4)
textpos <- matrix(c(c(ext(predictors)[1] + diff(ext(predictors)[1:2])/2), ext(predictors)[4]), ncol = 2) * c(1.02, 1.02) 
textadj <- c(0.5, 0)

plot(invlogit(b3 * predictors), mar = mar); text(textpos, label = "Gaussian field predictor (GF)\n(from prob = 0.5)", xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

plot(invlogit((b2 * binary) - (b2/2)), mar = mar);text(textpos, label = "Binary predictor\n(from prob = 0.5)", xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

plot(gradient, mar = mar);text(textpos, label = "Gradient predictor\n(from prob = 0.5)", xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

plot(raw, mar = mar, col = "white", legend = FALSE);text(textpos, label = "Observations and INLA mesh", xpd = TRUE, adj = textadj)
par(mar = mar)
plot(smesh, border = adjustcolor("black", 0.05), add = TRUE)
plot(st_geometry(obs), pch = 16, cex = 0.75, col = adjustcolor("black", 0.25), add = TRUE)

eff <- effort
eff[eff == 0] <- NA
plot(eff, mar = mar);text(textpos, label = paste("Effort (",global(effort,"sum")[1,1],"checklists)"), xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

plot(raw, mar = mar); text(textpos, label = "Raw proportions", xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

r3 <- r2
r3[r3 < 1] <- NA
plot(r3, mar = mar); text(textpos, label = "Aggregated effort modeled\n(with low effort as numbers)", xpd = TRUE, adj = textadj)
text(st_coordinates(sprinklewhere), label = sprinklewhere$n, cex = 0.6, col = adjustcolor("black", 0.5))

plot(r1/r2, mar = mar); text(textpos, label = "Aggregated raw proportions modeled", xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

zlim <- range(c(values(truth)[, 1], values(preds)[, 1], values(gampreds)[, 1]), na.rm = FALSE)
plot(truth, mar = mar, range = zlim); text(textpos, label = "Truth\n(prob = GF + binary + gradient)", xpd = TRUE, adj = textadj)
#plot(st_geometry(obs), cex = 0.5, lwd = 0.2, add = TRUE)

plot(preds, mar = mar, range = zlim); text(textpos, label = "Predictions\n(from a purely spatial model)", xpd = TRUE, adj = textadj)

plot(gampreds, mar = mar, range = zlim); text(textpos, label = "Predictions\n(from a GAM)", xpd = TRUE, adj = textadj)

plot(unc, mar = mar); text(textpos, label = "Predictions sd\n(from the purely spatial model)", xpd = TRUE, adj = textadj)

dev.off()
system("xdg-open simGaussianFields.png")


#x <- seq(0, 100000, by = 1)
#y <- invlogit(-0.000000002 * x^2 + 1)
#plot(x, y, type = "l")


library(geodata)
canada <- gadm("GADM", country = "CAN", level = 1) |>
  st_as_sf()  # transform a SpatVector to an sf data.frame

