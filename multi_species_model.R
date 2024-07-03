# Multi-species occupancy modeling of ground-dwelling mammals in central Laos: a case study for monitoring in tropical forests 

# Journal: Wildlife Biology

# Authors: Ioannis Alexiou, Camille NZ Coudrat, Jürgen Niedballa, Andreas Wilting and Andrew Tilker

# R script author: Ioannis Alexiou

# 23/06/2024




 
# Necessary files
# 
# The following files are needed to run this R script:
# 1. CTtable.csv (information on camera traps)
# 2. recorded_community.csv (information on detected taxa)
# 3. covariates.csv (table of site covariates)
# 4. prediction.tif (covariates raster)
# 5. aoi.shp (the area of interest vector)
# 
# Those files are available from the authors upon request.

# Working directory: please include all input files (csv, tif etc) at the same directory as this R script



# R packages 

library(camtrapR)
library(DT)
library(knitr)
library(ggplot2)
library(sf)
library(coda)
library(raster)
library(sp)
library(purrr)
library(terra)

# Camera trap station information file

camtraps <- read.csv("CTtable.csv", row.names=1)

# Camera trap station operation matrix

camop <- cameraOperation(CTtable      = camtraps,
                         stationCol   = "Station_correct",
                         cameraCol = "Camera",
                         setupCol     = "Setup_date",
                         retrievalCol = "Retrieval_date",
                         hasProblems  = TRUE,
                         byCamera = FALSE,
                         camerasIndependent = FALSE,
                         allCamsOn = FALSE)

# List of recorded species

recordTable <- read.csv("recorded_community.csv", row.names=1)

# Detection histories

DetHist_list <- lapply(unique(recordTable$Species), FUN = function(x) {
  detectionHistory(
    recordTable         = recordTable,
    camOp                = camop,
    stationCol           = "combined",
    speciesCol           = "Species",
    recordDateTimeCol    = "DateTimeOriginal",
    species              = x,     
    occasionLength       = 20,
    day1                 = "station",
    datesAsOccasionNames = FALSE,
    includeEffort        = TRUE,
    scaleEffort          = TRUE,
    timeZone             = "Asia/Ho_Chi_Minh")})

# assign species names to list of detection histories

names(DetHist_list) <- unique(recordTable$Species)

# List of detected histories

ylist <- lapply(DetHist_list, FUN = function(x) x$detection_history)

# Covariates 

covtable <- read.csv("covariates.csv")

# Data list

data_list <- list(ylist    = ylist,
                  siteCovs = covtable,
                  obsCovs  = list(effort = DetHist_list[[1]]$effort))

# Prepare model

modelFile_jags <- tempfile(fileext = ".txt")

model1_jags <- communityModel(data_list,
                              occuCovs = list(ranef = c("Elevation", 
                                                        "Remoteness", 
                                                        "TRI",
                                                        "Canopy_height",
                                                        "Elevation_squared")),
                              detCovsObservation = list(ranef = c("effort", "zone")),
                              intercepts = list(det = "ranef", occu = "ranef"),
                              model = "RN",
                              modelFile = modelFile_jags,
                              keyword_quadratic = "_squared")

summary(model1_jags)

# Run model

out_ahm_jags <- fit(model1_jags, 
                    n.iter = 60000, 
                    n.burnin = 30000,
                    thin = 10,
                    chains = 3,
                    cores = 3,      
                    quiet = F
)

# Analyze model results

# Create mcmc.list object

out_ahm_jags <- mcmc.list(out_ahm_jags)

# Summarize results

out_ahm_jags_summary <- summary(out_ahm_jags)

DT::datatable(round(out_ahm_jags_summary$statistics, 2))

# Estimate Gelman-Rubin convergence statistic

gd <- gelman.diag(out_ahm_jags,  multivariate = FALSE)

# Results data frame

stats <- as.data.frame(cbind(out_ahm_jags_summary$statistics, out_ahm_jags_summary$quantiles, gd$psrf))

# Plot beta coefficients

plot_coef(model1_jags,
          out_ahm_jags,
          "state",
          combine=T,
          level = c(outer = 0.95, inner = 0.75)) + theme_minimal()

# Plot response curves

plot_effects(model1_jags, 
             out_ahm_jags,
             submodel = "state")

plot_effects(model1_jags, 
             out_ahm_jags,
             submodel = "det")

# Spatial predictions

# Load stack prediction

stack_prediction <- rast("prediction.tif")

# Load AOI (area of interest) vector

aoi <- vect("aoi.shp")

# Convert AOI to a raster

r_aoi <- rast(aoi, nrow=100, ncol=100)

# PAO (Percentage of Area Occupied) estimates

predictions_pao <- predict(object    = model1_jags, 
                           mcmc.list = out_ahm_jags,
                           x         = stack_prediction,
                           type      = "pao", 
                           aoi       = r_aoi,
                           draws = 10)

# Plot PAO

ggplot2::ggplot(predictions_pao$pao_df, aes(x = Species , y = PAO)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# species occupancy estimates

predictions_psi <- predict(object    = model1_jags, 
                           mcmc.list = out_ahm_jags,
                           x         = stack_prediction,
                           type      = "psi")

# Plot species occupancy

plot(predictions_psi$mean, zlim = c(0,1), 
     col = hcl.colors(100), 
     maxnl = 9,   # plotting only the first 9 species for space reasons
     asp = 1)  

# Species richness estimates

predictions_rich <- predict(object   = model1_jags, 
                            mcmc.list = out_ahm_jags,
                            x         = stack_prediction,
                            type      = "richness")

# Plot species richness

plot(predictions_rich, col = hcl.colors(100), asp = 1)

# sessionInfo()

# R version 4.1.0 (2021-05-18)
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] purrr_0.3.4    coda_0.19-4    ggplot2_3.4.2  knitr_1.40     DT_0.25       
# [6] sf_1.0-12      raster_3.6-20  sp_1.6-0       terra_1.7-23   camtrapR_2.3.0
# 
# loaded via a namespace (and not attached):
#   [1] xfun_0.33           RcppNumerical_0.5-0 tidyselect_1.1.2    splines_4.1.0      
# [5] lattice_0.20-44     colorspace_2.1-0    vctrs_0.6.1         generics_0.1.3     
# [9] htmltools_0.5.2     mgcv_1.8-35         utf8_1.2.3          rlang_1.1.1        
# [13] e1071_1.7-13        pillar_1.9.0        withr_3.0.0         glue_1.6.2         
# [17] DBI_1.2.2           secr_4.5.10         lifecycle_1.0.4     stringr_1.5.1      
# [21] munsell_0.5.1       gtable_0.3.4        htmlwidgets_1.5.4   codetools_0.2-18   
# [25] fastmap_1.1.0       parallel_4.1.0      class_7.3-19        fansi_1.0.4        
# [29] Rcpp_1.0.10         KernSmooth_2.23-20  scales_1.2.1        classInt_0.4-9     
# [33] RcppParallel_5.1.7  abind_1.4-5         digest_0.6.27       stringi_1.7.12     
# [37] dplyr_1.0.10        grid_4.1.0          cli_3.2.0           tools_4.1.0        
# [41] magrittr_2.0.1      proxy_0.4-27        tibble_3.2.1        pkgconfig_2.0.3    
# [45] MASS_7.3-54         Matrix_1.5-1        data.table_1.14.8   timechange_0.2.0   
# [49] lubridate_1.9.2     assertthat_0.2.1    rstudioapi_0.14     R6_2.5.1           
# [53] units_0.8-1         nlme_3.1-152        compiler_4.1.0     
