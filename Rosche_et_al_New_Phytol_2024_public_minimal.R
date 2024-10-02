#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#> Herbarium specimens reveal a cryptic invasion of tetraploid Centaurea stoebe in Europe
#>
#> Study on spatio-temporal range dynamics of polyploid complexes across sympatric
#> distribution of cytotypes of Centhaurea stoebe improve our understanding on how niche dynamics,
#> habitat type preferences and colonization success can be related to a plant´s ploidy level.
#>
#> This script contains the code to produce all results presented in the manuscript.
#>
#> AUTHORS:
#> Olivier Brönnimann (initial code on PCA, MCA, niche/range dynamics),
#> Christoph Rosche (GAMs, BRTs),
#> Julian Selke (niche/range dynamics, functions)
#>
#> Corresponding author: christoph.rosche@botanik.uni-halle.de


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

set.seed(2024)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Setup                                                                                         ####
# installing/loading packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

packages_needed <- c(
  "ade4",               # PCA
  "conflicted",         # manage conflicted functions
  "cowplot",            # plot arrangement
  "data.table",         # data manipulation
  "dismo",              # mess model
  "docstring",          # function documentation without building a package
  "dplyr",              # data manipulation
  "dynRB",              # estimate niche overlap
  "ecospat",            # dispersal route calculation
  "emojifont",          # up and down arrows
  "ggblend",            # plotting; color overlap
  "ggnewscale",         # multiple color scales in ggplots
  "ggplot2",            # plotting
  "ggpubr",             # plot arrangement
  "ggrepel",            # plots
  "ggtext",             # text style in plots
  "grid",               # plotting
  "gridExtra",          # plot arrangement
  "lme4",               # statistics
  "MASS",               # statistics
  "metR",               # spatial data plotting
  "mgcv",               # statistics
  "multcomp",           # statistics
  "nlme",               # statistics
  "patchwork",          # plot arrangement
  "plyr",               # data manipulation
  "prismatic",          # color manipulation
  "raster",             # spatial data manipulation
  "readxl",             # reading excel data
  "rnaturalearth",      # world plot
  "rnaturalearthdata",  # world plot
  "sf",                 # spatial data handling / plotting
  "skimr",              # data summary
  "sloop",              # object information
  "sp",                 # spatial data handling / plotting
  "stringr",            # string manipulation
  "terra",              # spatial data manipulation
  "tibble",             # data manipulation
  "tidyr",              # data manipulation
  "tidyterra",          # spatial data manipulation / plotting
  "tidyverse",          # data manipulation
  "vegan",              # RDA
  "viridis",            # colors
  "wesanderson",        # colors
  NULL                  # dummy
)

# check if packages need to be installed
new_packages <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]

# value greater than zero coerces to TRUE; install missing packages
if(length(new_packages)) install.packages(new_packages)

# load all packages (and confirm success)
all(unlist(lapply(packages_needed, require, character.only = TRUE)))

# cleanup
rm(new_packages, packages_needed)

# use functions from respective packages
conflicted::conflicts_prefer(
  docstring::`?`,
  dplyr::mutate,
  dplyr::filter,
  dplyr::select
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# File import ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cen <- readxl::read_xlsx("Dataset_S1_Rosche_2024_New_Phytologist_Centaurea.xlsx",
                         sheet = 1,
                         na = c("NA", ""))

areal <- raster::shapefile("Areal.shp")

areal_poly <- sf::st_read("Areal.shp")

# if does not exist, file will be created later in the script
if (file.exists("clim_europe.tif")) clim <- terra::rast("clim_europe.tif")

varimp_all_rec <- read.csv("varimp_all_recent_.csv", header = TRUE, sep = ";", dec = ".")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Global variables ####
# definition of parameters used throughout the script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

TIME_SPAN  = 1800:2022
QUANT_HIGH = 0.9
QUANT_LOW  = 0.1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Functions ####
# convenience functions for this script; inline documentation with 'docstring' package
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

show_progress <- function() {
  # NOT AVAILABLE VIA DOCSTRING!
  #' Show Loop Progress
  #'
  #' @description Print progress bar, percentage, current iteration, and, optionally, estimation
  #' of remaining time until completion of loop to console.
  #'
  #' @param i integer. The loop counter. Expected to start at 1.
  #' @param n integer. The maximum loop count, i.e. length of the data.
  #' @param names A character vector of length n. The names of the dimension
  #' being looped over, e.g. row.names. Defaults to iteration number.
  #' @param time logical. Indication of whether to print an estimate of time to
  #' complete remaining iterations. See Details section for implementation.
  #' @param extra_poly integer. Degree of the polynomial used for estimation of
  #' remaining time. Default is 0, i.e. mean duration of loops times
  #' number of remaining loops is used for time estimation. If set to -1, the
  #' default is up to loop 50 and quadratic estimation for remaining loops.
  #' @details The time estimation is implemented by calculating a time difference
  #' of the current and the last k iterations where k is i-1. Values are stored
  #' in a cache for estimation of remaining time. The cache is cleared when the
  #' loop is finished, i.e. when \eqn{i == n}, when the the loop is not finished
  #' (i.e. \eqn{i == n} was not reached) but \eqn{i == 1}, or when the
  #' function is called in a different environment. However, the cache is not
  #' cleared if the last call did not complete the loop and the current call
  #' does not start at \eqn{i == 1}. This may be the case if the loop is
  #' continued, but a falsely set \eqn{i = x | x > 1} as start value will not
  #' be recognized as an error.

  # convenience function
  pretty_time <- function(x) {
    if (x > 90 * 60 * 24) return(paste(round(x / (60 * 60 * 24)), "days"))
    if (x > 90 * 60) return(paste(round(x / (60 * 60)), "hours"))
    if (x > 90) return(paste(round(x / 60), "minutes"))
    return(paste(round(x), "seconds"))
  }
  # data to enclose by function
  time_cache   <- c()
  time_stamp   <- 0
  loop_count   <- 0
  env_sentinel <- rlang::caller_env()
  # sentinel is set when i == n, i.e. loop completed;
  # if loop terminates in error, cache is cleared on next call
  completed <- FALSE
  # closure to return
  clsr <- function(i, n, names = NA_character_, time = FALSE, extra_poly = 0) {
    if (identical(names,  NA_character_)) names <- as.character(1:n)
    if (!is.numeric(i) || !is.numeric(n) || i <= 0 || n <= 0) {
      stop("Both arguments 'i' and 'n' must be positive integers.")
    }
    if (i > n) {
      stop("Index 'i' must not be greater than length of data 'n'.")
    }
    if (!i == as.integer(i)) {
      warning("Removing decimal places for indexing.")
      i <- as.integer(i)
    }
    if (!is.character(names)) {
      stop("Argument 'names' must be a character vector.")
    }
    # when restarting after error/abort or called in another environment
    if ((completed == FALSE & i == 1) |
        !identical(env_sentinel, rlang::caller_env())) {
      # cleanup
      time_cache   <<- c()
      time_stamp   <<- 0
      loop_count   <<- 0
      env_sentinel <<- rlang::caller_env()
    }
    loop_count <<- loop_count + 1
    if (loop_count == 1 & i != 1) stop("Index 'i' must start at 1.")
    cat("\014\n")
    extra <- nchar("||100%")
    width <- options()$width
    step  <- round(i / n * (width - extra))
    text  <- sprintf(
      "|%s%s|% 3s%%",
      strrep("=", step),
      strrep(" ", width - step - extra),
      round(i / n * 100)
    )
    if (time == TRUE) {
      if (i == 1) {
        # message for first timediff calculation
        cat("Calculating remaining time...\n")
        time_stamp <<- Sys.time()
      } else {
        now <- Sys.time()
        loop_dur <- difftime(now, time_stamp, units = "secs")
        # smoothing over cached values
        time_cache <<- c(time_cache, loop_dur)
        # switch from linear prediction to polynomial when sufficient data are collected
        if (length(unique(time_cache)) <= extra_poly |
            (extra_poly == -1 & i < 50) |
            extra_poly == 0) {
          time_remaining <- mean(time_cache) * ((n - i) + 1)
        } else {
          # adjust degree after switching to quadratic estimation
          if (extra_poly == -1) extra_poly <- 2
          time_data <- data.frame(iterations = 1:length(time_cache),
                                  durations = time_cache)
          pred_data <- data.frame(iterations = (i + 1):n)
          pred_data$predictions <-
            stats::predict(lm(durations ~ poly(iterations, degree = extra_poly, raw = TRUE),
                              data = time_data),
                           newdata = pred_data)
          time_remaining <- sum(pred_data$predictions)
        }
        cat(
          "Approximately",
          pretty_time(time_remaining),
          "remaining.\n"
        )
        # save for next iteration
        time_stamp <<- Sys.time()
      }
    }
    cat(text)
    cat("\nProcessing ", names[i], "\n")
    if (i == n) {
      # cleanup
      completed    <<- TRUE
      time_cache   <<- c()
      loop_count   <<- 0
      env_sentinel <<- rlang::caller_env()
      cat("\r\n", crayon::green$bold(">>> COMPLETED"), "\n\n")
    } else {
      completed <<- FALSE
    }
  }
  return(clsr)
}

# assign return of higher order function, i.e. function showing progress, to desired name
show_progress <- show_progress()



compute_niche_time <- function(data_info = cen.info,
                               data_scores = cen.scores,
                               qhigh = QUANT_HIGH,
                               qlow = QUANT_LOW,
                               timespan = TIME_SPAN) {
  #' Compute niche sizes
  #'
  #' @description Computation of niche breadth through time based on PCA on BIOCLIM variables
  #' extracted respective coordinates. Quantile values for individual years are computed for a
  #' presentation more robust and smooth.
  #'
  #' @param data_info data.frame. Data on habitat, range, etc. Column names must not be changed!
  #' @param data_scores data.frame. Data on principal components.
  #' @param qlow numeric. 0 <= qlow <= qhigh <= 1
  #' @param qhigh numeric. 0 <= qlow <= qhigh <= 1
  #' @param timespan numeric. Vector of years to include.
  #' @return data.frame
  #' @export
  #' @import tidyr, dplyr, stringr
  arr <- array(dim = c(length(timespan), 3, 4, 2, 2, 2),
               dimnames = list(timespan,
                               c("all", "cwn", "native"),
                               c("xnq", "xrq", "xnq.test", "xrq.test"),
                               c(1, 2),
                               c(qlow, qhigh),
                               c("diploid", "tetraploid")))

  loop_data <- expand_grid(c("all", "cwn", "native"),
                           c("diploid", "tetraploid"),
                           1:2,
                           c(qlow, qhigh),
                           TIME_SPAN) %>%
    tidyr::unite(data = ., col = "names", colnames(.), remove = FALSE)
  loop_count <- 0
  for (sub_range in c("all", "cwn", "native")) {
    for (ploidy in c("diploid", "tetraploid")) {
      for (pc in 1:2) {
        for (q in c(qlow, qhigh)) {
          for (y in TIME_SPAN) {
            loop_count <- loop_count + 1
            show_progress(i = loop_count,
                          n = nrow(loop_data),
                          names = loop_data$names,
                          time = TRUE)

            if (sub_range == "all") {
              s <- rep(TRUE, nrow(data_info))
            } else if (sub_range == "cwn") {
              s <- data_info$Study_range == "yes" & data_info$range == "expanded"
            } else if (sub_range == "native") {
              s <- data_info$range == "native"
            }
            n       <- data_info$ploidy_level == ploidy
            year    <- data_info$year < y
            natural <- data_info$habitat == "natural"
            ruderal <- data_info$habitat == "ruderal"
            xn      <- data_scores[which(s & n & year & natural), pc]
            xr      <- data_scores[which(s & n & year & ruderal), pc]
            if (length(xn) == 0) {
              xnq <- NA
            } else {
              xnq <- quantile(xn, q)
            }
            if (length(xr) == 0) {
              xrq <- NA
            } else {
              xrq <- quantile(xr, q)
            }
            if (length(xr) == 0) xrq <- NA
            if (!is.na(xnq) & !is.na(xrq)) {
              repxnq <- repxrq <- c()
              for (i in 1:1000) {
                repxnq <- c(repxnq, quantile(xn[sample(1:length(xn), length(xn), replace = TRUE)],
                                             q))
                repxrq <- c(repxrq, quantile(xr[sample(1:length(xr), length(xr), replace = TRUE)],
                                             q))
              }
              xnq.test <- (xnq > quantile(repxrq, qlow)) * (xnq < quantile(repxrq, qhigh))
              xrq.test <- (xrq > quantile(repxnq, qlow)) * (xrq < quantile(repxnq, qhigh))
            } else {
              xnq.test <- 1
              xrq.test <- 1
            }
            range_index    <- which(c("all", "cwn", "native") == sub_range)
            TIME_SPAN_index <- which(timespan == y)
            quantile_index <- which(c(qlow, qhigh) == q)
            ploidy_index   <- which(c("diploid", "tetraploid") == ploidy)
            arr[TIME_SPAN_index, range_index, 1, pc, quantile_index, ploidy_index] <- xnq
            arr[TIME_SPAN_index, range_index, 2, pc, quantile_index, ploidy_index] <- xrq
            arr[TIME_SPAN_index, range_index, 3, pc, quantile_index, ploidy_index] <- xnq.test
            arr[TIME_SPAN_index, range_index, 4, pc, quantile_index, ploidy_index] <- xrq.test
          }
        }
      }
    }
  }
  tmp_arr <- arr %>%
    as.data.frame.table() %>%
    tidyr::pivot_wider(names_from = Var5, values_from = Freq) %>%
    `names<-`(c("year", "range", "type", "PC", "ploidy", "q0.1", "q0.9")) %>%
    dplyr::mutate(year = as.numeric(as.character(year)),
                  type = stringr::str_replace(type, "xnq", "natural"),
                  type = stringr::str_replace(type, "xrq", "ruderal"),
                  type = as.factor(type),
                  PC = as.factor(paste("PC", PC)),
                  m_type = gsub("\\.test", "", .$type))
  arr <- merge(tmp_arr %>%
                 dplyr::filter(stringr::str_detect(string = .$type,
                                                   pattern = "test",
                                                   negate = TRUE)),
               tmp_arr %>%
                 dplyr::filter(stringr::str_detect(string = .$type,
                                                   pattern = "test",
                                                   negate = FALSE)),
               by = c("year", "range", "PC", "ploidy", "m_type")) %>%
    `names<-`(c("year", "range", "PC", "ploidy", "m_type", "habitat",
                "q_low", "q_high", "test", "t_low", "t_high"))
  return(arr)
}

# convenience function
info <- function(object) {
  #' Print object information
  #'
  #' @description Prints object information:
  #' 1. human readable size in memory
  #' 2. object type in relation to R's OOP system
  #' 3. storage mode
  #' 4. class within the OOP system
  #' 5. environment where the object's name is bound
  #' 6. names of the object's attributes
  #'
  #' @param object any object
  #' @return NULL
  #' @export
  #' @import sloop
  message(
    paste0(
      "size:\t\t", format(object.size(object), unit = "auto"),
      "\nobj.type:\t", sloop::otype(object),
      "\nstg.mode:\t", storage.mode(object),
      "\nclass:\t\t", paste(class(object), collapse = ", "),
      "\nnamespace:\t", paste(find(deparse(substitute(object))), collapse = ", "),
      "\nattributes:\t", paste(names(attributes(object)), collapse = ", ")
    )
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data import ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cen <- cen %>%
  mutate(across(c(herbarium_voucher,
                  country,
                  habitat,
                  range,
                  EUNIS,
                  ploidy_level), ~ as.factor(.x))) %>%
  mutate(year = as.numeric(str_extract(collection_date, "^\\d{4}")))

skimr::skim(cen)

info(areal)

info(areal_poly)

areal_df <- areal %>% .fortify()
head(areal_df)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# EXECUTE ONLY ONCE
# aggregate climate variables, crop to land area, and to store them in a single file


if (!file.exists("clim_europe.tif")) {
  e <- raster::extent(-12, 58, 35, 71)
  # IMPORTANT: Apparently, this source link below is not stable. However, the files should still be
  #            available on https://envicloud.wsl.ch/
  # source: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
  # -> climatologies -> 1981-2010 -> bio
  # download date: 29 Sept. 2023
  files <- paste0("../envicloud/chelsa/chelsa_V2/GLOBAL/",
                  "climatologies/1981-2010/bio/CHELSA_bio", 1:19, "_1981-2010_V.2.1.tif")
  clim <- terra::rast()
  # source: https://www.sciencebase.gov/catalog/item/508fece8e4b0a1b43c29ca22
  # -> Attached Files -> download shapefiles / click 'download all'
  # download date: 29 Sept. 2023
  land <- terra::vect("../official/wwf_terr_ecos.shp") %>% terra::aggregate()

  for(i in 1:length(files)) {
    temp <- terra::rast(files[i]) %>% terra::crop(ext(e)) %>% terra::mask(land)
    clim <- c(clim, temp)
  }

  terra::writeRaster(clim, filename = "clim_europe.tif")
  clim <- terra::rast("clim_europe.tif")

}

info(clim)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Proportion of tetraploids to diploids through time, space and habitat type ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . GAMs: (1)	Accounting for spatial autocorrelation ####

gamA <- gam(species ~ s(year),
            method = "REML", family = "binomial", data = cen)

gamB <- gam(species ~ s(year) + s(country, bs = "re"),
            method = "REML", family = "binomial", data = cen)

gamC <- gam(species ~ s(year) + s(latitude, longitude, bs = "tp"),
            method = "REML", family = "binomial", data = cen)

gamD <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
            method = "REML", family = "binomial", data = cen)

summary(gamA)
gam.check(gamA)

summary(gamB)
gam.check(gamB)

summary(gamC)
gam.check(gamC)

summary(gamD)
gam.check(gamD)

AIC(gamA, gamB, gamC, gamD)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . GAMs: (2)	Different smooths between both ranges ####

gamE <- gam(species ~ s(year, by = range) + range + s(latitude, longitude, bs = "sos"),
            method = "REML", family = "binomial", data = cen)

gamF <- gam(species ~ s(year) + range + s(latitude, longitude, bs = "sos"),
            method = "REML", family = "binomial", data = cen)

summary(gamE)
gam.check(gamE)

summary(gamF)
gam.check(gamF)

AIC(gamD, gamE, gamF)

plot.gam(gamD)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . GAMs: (3)	Different smooths between habitat types in the native range ####

cen_native_HT <- cen %>% filter(range == "native" & !is.na(habitat))

gamG <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
            method = "REML", family = "binomial", data = cen_native_HT)

gamH <- gam(species ~ s(year) + habitat + s(latitude, longitude, bs = "sos"),
            method = "REML", family = "binomial", data = cen_native_HT)

gamI <- gam(species ~ s(year, by = habitat) + habitat + s(latitude, longitude, bs = "sos"),
            method = "REML", family = "binomial", data = cen_native_HT)

summary(gamG)
gam.check(gamG)

summary(gamH)
gam.check(gamH)

summary(gamI)
gam.check(gamI)

AIC(gamG, gamH, gamI)

##### significance of the predicts in the habitats of the native range

cen_native_HT_natural <- cen_native_HT %>% filter(habitat == "natural")

gam_native_natural <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
                          family = "binomial",
                          method = "REML",
                          data = cen_native_HT_natural)

summary(gam_native_natural)

gam_native_ruderal <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
                          family = "binomial",
                          method = "REML",
                          data = cen_native_HT_natural)

summary(gam_native_ruderal)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . GAMs: (4)	Different smooths between habitat types in the expanded range ####

cen_expanded_HT <- cen %>% filter(range == "expanded" & !is.na(habitat))

gamJ <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
            method = "REML",
            family = "binomial",
            data = cen_expanded_HT)

gamK <- gam(species ~ s(year) + habitat + s(latitude, longitude, bs = "sos"),
            method = "REML",
            family = "binomial",
            data = cen_expanded_HT)

gamL <- gam(species ~ s(year, by = habitat) + habitat + s(latitude, longitude, bs = "sos"),
            method = "REML",
            family = "binomial",
            data = cen_expanded_HT)


summary(gamJ)
gam.check(gamJ)

summary(gamK)
gam.check(gamK)

summary(gamL)
gam.check(gamL)

AIC(gamJ, gamK, gamL)

##### significance of the predicts in the habitats of the expanded range


gam_expanded_natural <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
                            family = "binomial",
                            method = "REML",
                            data = cen_expanded_HT %>% filter(habitat == "natural"))

summary(gam_expanded_natural)


gam_expanded_ruderal <- gam(species ~ s(year) + s(latitude, longitude, bs = "sos"),
                            family = "binomial",
                            method = "REML",
                            data = cen_expanded_HT %>% filter(habitat == "ruderal"))

summary(gam_expanded_ruderal)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Geographical range size and the climatic niche breadth of tetraploids through ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . Dispersal routes ####

set.seed(2024)

pt <- st_cast(areal_poly, "POINT") %>%
  unlist() %>%
  .[. > 0] %>%
  matrix(ncol = 2, byrow = TRUE) %>%
  .[sample(1:nrow(.), 500), ]

plot(pt)

NWE4x <- cen %>%
  filter(ploidy_level == "tetraploid" & range == "expanded") %>%
  select(longitude, latitude, year)

pt <- data.frame(cbind(pt, 1850))
names(pt) <- names(NWE4x)
d <- rbind(pt, NWE4x)

set.seed(2024)
rte <- ecospat.mdr(data = d,
                   xcol = 1,
                   ycol = 2,
                   datecol = 3,
                   mode = "min",
                   rep = 50,
                   mean.date.error = 2,
                   fixed.sources.rows = which(d$year == 1850))

mca <- rte[[1]]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . MESS ####

# min/max long/lat
e <- ext(c(quantile(NWE4x$longitude, c(0, 1)) + c(-1, 1),
           quantile(NWE4x$latitude,  c(0, 1)) + c(-1, 1)))

e[1] <- -5

clim10 <- crop(clim, e)

Native4x <- cen %>%
  dplyr::filter(ploidy_level == "tetraploid" & range == "native") %>%
  dplyr::select(longitude, latitude, year)

clim10.native4x <- terra::extract(clim10, Native4x[, 1:2])
clim10.native4x <- clim10.native4x[, -1]

clim10 <- raster::stack(clim10)
MESS <- dismo::mess(clim10, clim10.native4x, full = FALSE)

mess_df <- as.data.frame(as(MESS, "SpatialPixelsDataFrame"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . Niche dynamics analysis ####

# sample climate for pca
set.seed(2024)

# PCA scores
cen.env <- terra::extract(clim, cbind(cen$longitude, cen$latitude))
cen.na <- is.na(cen.env[, 1])
# TRUE == 1 & FALSE == 0
sum(cen.na) # 46 cen removed
cen.env <- cen.env[!cen.na, ]
cen.info <- cen[!cen.na, ]
names(cen.env) <- names(clim)

pca <- rda(cen.env, scale = TRUE)
biplot(pca, display = 'species', scaling = 'species')

cen.scores <- predict(pca, cen.env, type = "wa", scaling = 1)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . Niche overlap ####

cen_no <- cen %>% filter(Study_range == "yes")
cen.env <- terra::extract(clim, cbind(cen_no$longitude, cen_no$latitude))
r <- dynRB_VPa(cbind(cen %>% filter(Study_range == "yes") %>% pull(ploidy_level), cen.env) %>% drop_na(),
               steps = 201, pca.corr = TRUE)
r$result

cen_no_4 <- cen %>% filter(Study_range == "yes" & ploidy_level == "tetraploid")
cen.env <- terra::extract(clim, cbind(cen_no_4$longitude, cen_no_4$latitude))
r4 <- dynRB_VPa(cbind(cen %>% filter(Study_range == "yes" & ploidy_level == "tetraploid") %>% pull(range), cen.env) %>% drop_na(),
                steps = 201, pca.corr = TRUE)
r4$result



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . . . . . Scale data ####

# data presented in final figures are 10% quantiles, i.e., 0.1 and 0.9
side_df <- expand.grid(1:2,
                       c("natural", "ruderal"),
                       c("all", "cwn", "native"),
                       c(0.1, 0.05, 0.025, 0),
                       c("diploid", "tetraploid")) %>%
  `names<-`(c("PC", "habitat", "range", "qlow", "ploidy")) %>%
  mutate(qhigh = case_when(qlow == 0.1 ~ 0.9,
                           qlow == 0.05 ~ 0.95,
                           qlow == 0.025 ~ 0.975,
                           qlow == 0 ~ 1,
                           .default = NA),
         qlow_val = NA, qhigh_val = NA)

qlows = c(0.1, 0.05, 0.025, 0)
qhighs = c(0.9, 0.95, 0.975, 1)

for(q in 1:length(qlows)) {
  for(pc in c(1, 2)) {
    for(ploidy in c("diploid", "tetraploid")) {
      for(hab in c("natural", "ruderal")) {
        for(ran in c("all", "cwn", "native")) {
          if (ran == "all") {
            s <- cen.info$year %in% TIME_SPAN
          } else if (ran == "cwn") {
            s <- cen.info$Study_range == "yes" &
              cen.info$range == "expanded" &
              cen.info$year %in% TIME_SPAN
          } else if (ran == "native") {
            s <- cen.info$range == "native" & cen.info$year %in% TIME_SPAN
          }

          dat <- cen.scores[!is.na(cen.info$ploidy_level == ploidy & cen.info$habitat == hab & s) &
                              cen.info$ploidy_level == ploidy & cen.info$habitat == hab & s,
                            pc]
          side_df[side_df$PC == pc &
                    side_df$habitat == hab &
                    side_df$qlow == qlows[q] &
                    side_df$range == ran &
                    side_df$ploidy == ploidy,
                  c("qlow_val", "qhigh_val")] <- quantile(dat, c(qlows[q], qhighs[q]))
        }
      }
    }
  }
}


df_10 <- compute_niche_time(data_info   = cen.info,
                            data_scores = cen.scores,
                            qhigh       = QUANT_HIGH,
                            qlow        = QUANT_LOW,
                            timespan    = TIME_SPAN)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# . . . . Range dynamics analysis ####

clim20 <- terra::ifel(clim[[1]] > -20, TRUE, FALSE)

aggr <- 10
yrange <- seq(min(cen$year), max(cen$year))

all_dat <- list(NA, NA ,NA)
names(all_dat) <- c("all", "cwn", "native")
loop_count <- 1

for (subset in c("all", "cwn", "native")) {
  cat("\n")
  if (subset == "all")  s <- 1:nrow(cen)
  else if (subset == "cwn") s <- which(cen$Study_range == "yes" & cen$range == "expanded")
  else if (subset == "native") s <- which(cen$range == "native")
  else stop("unkown subset")

  grid <- aggregate(clim20, aggr, fun = "max", na.rm = TRUE)
  grid <- grid - 1
  grid <- rep(grid, 4)
  names(grid) <- c("4xn","4xr","2xn","2xr")

  ygrid <- grid
  cells <- 0
  ncells <- cbind(yrange, NA, NA, NA, NA)
  colnames(ncells) <- c("yrange", "4xn", "4xr", "2xn", "2xr")

  for (ploidy in c("diploid", "tetraploid")) {
    for (habitat in c("natural", "ruderal")) {
      cat("\n")
      for (y in yrange){
        keep <- base::intersect(which(cen$year < y &
                                        cen$ploidy_level == ploidy &
                                        cen$habitat == habitat),
                                s)

        o <- terra::vect(cbind(cen[keep, "longitude"], cen[keep, "latitude"]),
                         geom = c("longitude", "latitude"))
        cells0 <- cells
        cells <- terra::extract(grid, o, cells = T)$cell

        #store
        p <- ifelse(ploidy == "diploid", "2x", "4x")
        h <- ifelse(habitat == "natural", "n", "r")
        ph <- paste0(p, h)
        grid[[which(names(grid) == ph)]][cells] <- 1
        ncells[which(ncells[, 1] == y), which(colnames(ncells) == ph)] <-
          as.numeric(terra::global(grid[[which(names(grid) == ph)]], "sum", na.rm = TRUE))

        ygrid[[which(names(ygrid) == ph)]][setdiff(cells, cells0)] <- y
        cat("\r")
        cat("\r", loop_count, "/ 12  --  Range: ",
            paste(subset, " -  Ploidy: ", ploidy, " -  Habitat: ", habitat, " -- "),
            paste(round(1:length(yrange)/length(yrange), 3)[y - min(yrange) + 1] * 100, "% "))
      }
      loop_count <- loop_count + 1
    }
  }
  names(ygrid) <- c("tetraploid - natural",
                    "tetraploid - ruderal",
                    "diploid - natural",
                    "diploid - ruderal")
  all_dat[[subset]] <- ncells
}

all_dat_p <- purrr::map_df(all_dat, ~ as.data.frame(.x), .id = "id")

all_dat_p %<>%
  dplyr::select(-yrange) %>%
  group_by(id) %>%
  nest() %>%
  mutate(max = max(unlist(data))) %>%
  dplyr::select(-data) %>%
  ungroup() %>%
  merge(all_dat_p, by = base::intersect(names(.), names(all_dat_p)), all.y = TRUE)

all_dat_p %<>%
  pivot_longer(c(`4xn`, `4xr`, `2xn`, `2xr`)) %>%
  mutate(ploidy = ifelse(str_detect(name, "2"), "diploid", "tetraploid"),
         habitat = ifelse(str_detect(name, "n"), "natural", "ruderal"))

str(all_dat_p)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Q4 - Determinants of the current spread of tetraploids in the expanded range ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# these data are not publicly available
if (file.exists("data_10km.rds")) {

  # extracting rds data for of several spatio-temporally explicit predictors
  t10 <- readRDS("data_10km.rds")
  dat10 <- t10$data

  summary(dat10)
  dat_spec <- cen[cen$year > 1988,]
  summary(dat_spec)
  dat <- merge(dat_spec, dat10, by = "sample_id", all = TRUE)
  summary(dat)

  # variable transformation
  dat <- dat %>%
    mutate(connectivity = 1 / travel,
           across(c(tetraplois, longitude, latitude), ~ as.numeric(.x)),
           habitat = as.factor(habitat),
           impervious = log(impervious + 0.001),
           urban = logit(urban),
           across(c(population, travel, connectivity), ~ log(.x)),
           across(c(road, rail), ~ log(.x / 1000)))

  summary(dat)

  # run brt current_spread
  gbm_current_spread <- gbm.step(dat = dat,
                                 # columns of explanatory variables
                                 gbm.x = c(4, 5, 9, 14:15, 17, 19:23),
                                 # column of response variable
                                 gbm.y = 3,
                                 distribution = "bernoulli",
                                 tree.complexity = 1,
                                 n.minobsinnode = 10,
                                 learning.rate = 0.01,
                                 bag.fraction = 0.5,
                                 tolerance.method = "fixed",
                                 tolerance = 0.01)

  summary(gbm_current_spread)

  gbm.plot(gbm_current_spread)
  gbm_current_spread.simp <- gbm.simplify(gbm_current_spread)
  gbm_current_spread.simp4 <- gbm.step(dat = dat,
                                       gbm.x = gbm_current_spread.simp$pred.list[[2]],
                                       tree.complexity = 1,
                                       gbm.y = 3,
                                       distribution = "binomial",
                                       learning.rate = 0.01,
                                       n.minobsinnode = 10)

  summary(gbm_current_spread.simp4)
  gbm.plot(gbm_current_spread.simp4)

  # summary
  varimp_all_current_spread <- summary(gbm_current_spread)
  varimp_all_current_spread <- varimp_all_current_spread[, 1:2]

  varimp_all_rec <- varimp_all_current_spread

}


# input
varimp_all_rec <- varimp_all_rec %>%
  mutate(direction = c(-1, -1, 1, 1, 1, 1, 0, 0, 1, -1, 0, 0)) %>%
  mutate(direction = case_when(direction == -1 ~ "\u21D8", # down arrow
                               direction == 1 ~ "\u21D7", # up arrow
                               .default = "n.s."))
skimr::skim(varimp_all_rec)


varimp_all_rec <- varimp_all_rec %>%
  mutate(var = factor(var)) %>%
  # names for plotting
  mutate(var = plyr::revalue(var, c("latitude" = "Latitude",
                                    "longitude" = "Longitude",
                                    "impervious" = "Impervious cover",
                                    "rail" = "Railway density",
                                    "distance" = "Spatial distance",
                                    "urban" = "Urban/rural ratio",
                                    "PC2" = "Temperature",
                                    "PC1" = "Precipitation ",
                                    "MESS" = "Climatic distance",
                                    "road" = "Road density",
                                    "population" = "Population density",
                                    "connectivity" = "Connectivity index"))) %>%
  mutate(categ = factor(c(1, 1, 2, 2, 2, 3, 1, 3, 4, 4, 3, 4),
                        levels = 1:4,
                        labels = c("Space", "Urban",
                                   "Dispersal", "Climate")))

summary(varimp_all_rec)

varimp_all_rec <- varimp_all_rec %>%
  mutate(ypos_lab = ifelse(rel.inf > 10, rel.inf - 1, rel.inf - 0.8),
         lab = paste0(round(rel.inf, 1), " %"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# session info ####
# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")





