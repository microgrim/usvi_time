# 06_permanova.R

# Resource allocation time ------------------------------------------------


if(file.exists(paste0(getwd(), "/", "00_resource_allocation.R"))){
  cat("Preparing resource allocations.")
  source(paste0(getwd(), "/", "00_resource_allocation.R"), local = FALSE,
         echo = TRUE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
  try(f_projectpath())
} else {
  cat("Preparing resource allocations.")
  
  #load packages
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(BiocManager)
  library(BiocParallel)
  library(data.table)
  library(cli)
  library(furrr)
  library(progressr)
  
  #determine multithreading capability
  if(grepl("arch64", Sys.getenv("R_PLATFORM"))){
    print("Detected Mac, using parallelly...")
    nthreads <- parallelly::availableCores(omit = 1) - 1
    future::plan(multisession, workers = nthreads)
    options(future.globals.maxSize = 10e9)
  } else {
    if(grepl("x86_64", Sys.getenv("R_PLATFORM"))){
      print("Detected Windows")
      nthreads <- parallelly::availableCores(omit = 1) - 1
      future::plan(sequential)
      options(future.globals.maxSize = 10e9)
    } else {
      print("Using data.table")
      nthreads <- data.table::getDTthreads()
      future::plan(sequential)
      options(future.globals.maxSize = 10e9)
    }
  }
  
  if(nthreads > 4){
    bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                  workers = nthreads - 2,
                                                  stop.on.error = TRUE,
                                                  RNGseed = 48105,
                                                  progressbar = TRUE)
  } else {
    bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                  workers = nthreads,
                                                  stop.on.error = TRUE,
                                                  RNGseed = 48105,
                                                  progressbar = TRUE)
  }
  register(bpparam_multi)
  
  
  
  #set project paths
  if(grepl("arch64", Sys.getenv("R_PLATFORM")) | grepl("usvi_temporal", getwd())){
    projectpath <- getwd()
  } else {
    if(grepl("vortex", getwd())){
      # projectpath <- "/proj/omics/apprill"
      projectpath <- "/vortexfs1/home/sharon.grim/projects/apprill/usvi_time"
    } else {
      projectpath <- "/user/sharon.grim/projects/apprill/usvi_time"
    }
  }
  
}

# Load additional packages ------------------------------------------------

library(tidyverse)
library(vegan)
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)
library(scales)
library(PERMANOVA)



# Custom functions --------------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_custom_functions.R"))){
  cat("Loading custom functions.")
  source(paste0(getwd(), "/", "00_custom_functions.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 

no_progress <- function() {} #your placeholder function for progress reporting

idx_samples <- function(x) grep("samples", names(x), value = TRUE)
coalesce2 <- function(x, y, sep = ".") ifelse(x == y, coalesce(x, y, sep = sep), paste0(x, "_vs_", y))

set.alpha <- 0.05
set.seed(48105)


# find existing processed files -------------------------------------------
to_import <- c("usvi_prok_asvs.df", "usvi_prok_asvs.taxa", "usvi_prok_decontam_idx", "usvi_hobo_light_temp_filtered.df")


for(file in to_import){
  if(!exists(file, envir = .GlobalEnv)){
    namevar <- file
    if(file.exists(paste0(projectpath, "/", namevar, ".tsv", ".gz"))){
      cli::cli_alert_info("Importing this dataset: {namevar}")
      temp_df <- readr::read_delim(paste0(projectpath, "/", namevar, ".tsv", ".gz"),
                                   col_names = TRUE, show_col_types = FALSE, delim = "\t", num_threads = nthreads)
      assign(paste0(namevar), temp_df, envir = .GlobalEnv)
      rm(temp_df)
      to_import <- to_import[grep(namevar, to_import, value = FALSE, invert = TRUE)]
    } else {
      cli::cli_alert_warning("Please prepare this dataset: {namevar}")
    }
  }
}

if(file.exists(paste0(projectpath, "/", "usvi_metadata_tidy.txt"))){
  metadata <- readr::read_delim(paste0(projectpath, "/", "usvi_metadata_tidy.txt"),
                                delim = "\t",
                                quote = "",
                                col_names = TRUE,
                                show_col_types = FALSE,
                                num_threads = nthreads)
} else {
  cli::cli_alert_warning("Please pre-process the metadata in the prior script step.")
}

if(file.exists(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"))){
  metabolomics_sample_metadata <- readr::read_delim(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
} else {
  cli::cli_alert_warning("Please process the metabolomics sample data previously.")
}

if(ncol(usvi_prok_decontam_idx) == 1){
  usvi_prok_decontam_idx <- usvi_prok_decontam_idx %>%
    tidyr::separate_wider_delim(1, names = c("asv_id", "keep"), delim = " ", cols_remove = TRUE, too_few = "align_start")
}

usvi_prok_asvs.taxa <- usvi_prok_asvs.taxa %>%
  dplyr::left_join(., usvi_prok_decontam_idx, by = join_by(asv_id)) %>%
  dplyr::filter(keep == TRUE) %>%
  dplyr::select(-keep) %>%
  droplevels

#replace NA in taxonomy with last known level

usvi_prok_filled.taxa.df <- usvi_prok_asvs.taxa %>%
  dplyr::mutate(Phylum = coalesce(Phylum, Domain)) %>%
  dplyr::mutate(Class = coalesce(Class, Phylum)) %>%
  dplyr::mutate(Order = coalesce(Order, Class)) %>%
  dplyr::mutate(Family = coalesce(Family, Order)) %>%
  dplyr::mutate(Genus = coalesce(Genus, Family)) %>%
  dplyr::mutate(Species = coalesce(Species, Genus)) %>%
  dplyr::mutate(across(everything(), ~factor(.x))) %>%
  dplyr::relocate(asv_id) %>%
  dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                          grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                          .default = Phylum)) %>%
  
  droplevels


usvi_metabolomics.df <- readr::read_delim(paste0(projectpath, "/", "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_wideFormat_noMetadata.csv"),
                                          col_names = TRUE, show_col_types = FALSE, delim = ",", num_threads = nthreads)
colnames(usvi_metabolomics.df)[1] <- "metab_deriv_label"


#there are samples "CINAR_BC_81A" and "CINAR_BC_81B" in the metabolomics dataset 
#and in the metadata, there are two DNA samples associated with "Deriv_81": Metab_219 (LB_seagrass dawn) and Metab_319 (tektite dawn)

usvi_metabolomics_long.df <- readr::read_delim(paste0(projectpath, "/", "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_longFormat_wMetadata.csv"), 
                                               col_select = c(2:last_col()),
                                               col_names = TRUE, show_col_types = FALSE, delim = ",", num_threads = nthreads) %>%
  dplyr::mutate(sample_id = paste0("Metab_", DNA_no))





#llight:

#it looks like at MIS, we used HOBOs deployed at multiple depths to 23m below surface to measure light
#Hobos record in Lux or Lumens, which is not easily translated to PAR
#to get the full light field, we also used multiple light profilers: Hyperspectral, LiCor, C-OPs
#the CONVERSION of the reported Hobo light in Lux, to uE was as follows:

#at 23m depth, only blue-green PAR was available as determined by the Blackbird hyperspectral profiler. 
#in July 2016, Blackbird measured ~75uE of available light almost entirely in Blue-Green wavelengths
#HOBOs measured light at 23m depth on those same days as ~996 lux
#so the light-profile/depth specific conversion was:
#Lux/54/0.245142572 = available PAR at 23m depth
#in this case, 996/54/.245142572 = 75 uE

#in MIS the attenuation coefficients averaged -0.132840405 (range: -0.161106824, -0.101399784)
#calculate depths and assign sampling_time windows before subsetting for dawn and peak_photo
# metabolomics_sample_metadata %>%
#   dplyr::mutate(day = lubridate::ymd(sampling_date)) %>%
#   dplyr::group_by(site, day, sampling_time) %>%
#   dplyr::summarise(depth = mean(depth))

if(!exists("usvi_hobo_light_temp_filtered.df", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_hobo_light_temp", ".tsv.gz"))){
    usvi_hobo_light_temp <- readr::read_delim(paste0(projectpath, "/", "usvi_hobo_light_temp", ".tsv.gz"))
    
    interval <- lubridate::as.interval(duration(hours = 4), ymd_hms("2021-01-21 11:00:00", tz = "America/Virgin")) 
    interval_vector_peak <- list(interval,
                                 lubridate::int_shift(interval, days(1)),
                                 lubridate::int_shift(interval, days(2)),
                                 lubridate::int_shift(interval, days(3)),
                                 lubridate::int_shift(interval, days(4)),
                                 lubridate::int_shift(interval, days(5)))
    
    interval <- lubridate::as.interval(duration(hours = 4), ymd_hms("2021-01-21 04:00:00", tz = "America/Virgin")) 
    interval_vector_dawn <- list(interval,
                                 lubridate::int_shift(interval, days(1)),
                                 lubridate::int_shift(interval, days(2)),
                                 lubridate::int_shift(interval, days(3)),
                                 lubridate::int_shift(interval, days(4)),
                                 lubridate::int_shift(interval, days(5)))
    
    usvi_hobo_light_temp_filtered.df <- usvi_hobo_light_temp %>%
      dplyr::mutate(date_ast = lubridate::force_tz(date_ast, tzone = "America/Virgin")) %>%
      dplyr::filter(grepl(paste0(c("2021-01-22", "2021-01-23", "2021-01-24", "2021-01-25", "2021-01-26"), collapse = "|"), day)) %>%
      dplyr::select(site, date_ast, day, time, temp, lux, lumens) %>%
      dplyr::rowwise(.) %>%
      dplyr::mutate(sampling_time = dplyr::case_when(date_ast %within% interval_vector_peak ~ "peak_photo",
                                                     date_ast %within% interval_vector_dawn ~ "dawn",
                                                     .default = NA)) %>%
      dplyr::mutate(sampling_time = factor(sampling_time)) %>%
      dplyr::left_join(., metabolomics_sample_metadata %>%
                         dplyr::rowwise(.) %>%
                         dplyr::mutate(day = lubridate::ymd(sampling_date)) %>%
                         dplyr::group_by(site, day, sampling_time) %>%
                         dplyr::summarise(depth = mean(depth)), by = join_by(site, day, sampling_time)) %>%
      dplyr::group_by(site, day) %>%
      tidyr::fill(depth, .direction = "downup") %>%
      dplyr::mutate(PAR = (lux/54)) %>%
      dplyr::mutate(PAR = PAR/exp(-0.13*(depth - 23))) %>%
      tidyr::pivot_longer(., cols = c("temp", "lumens", "lux", "PAR"),
                          # tidyr::pivot_longer(., cols = c("temp", "lumens", "lux"),
                          names_to = "parameter",
                          values_to = "value") %>%
      dplyr::group_by(site) %>%
      dplyr::slice_head(prop = 0.9) %>%
      droplevels
    
    readr::write_delim(usvi_hobo_light_temp_filtered.df, paste0(projectpath, "/", "usvi_hobo_light_temp_filtered.df", ".tsv.gz"),
                       delim = "\t", col_names = TRUE, num_threads = nthreads)
  }
}

usvi_light_temp_metadata.df <- bind_rows(
  (usvi_hobo_light_temp_filtered.df %>%
     dplyr::group_by(site, day, sampling_time, parameter) %>%
     dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>% dplyr::slice_max(n = 3, by = c(site, day, sampling_time,parameter), order_by = max, with_ties = FALSE) %>% #do you want to randomly pick the top 3 in each?
     droplevels)) %>%
  tidyr::drop_na(.) %>%
  dplyr::left_join(., metabolomics_sample_metadata %>%
                     dplyr::rowwise(.) %>%
                     dplyr::mutate(day = lubridate::ymd(sampling_date)) %>%
                     dplyr::group_by(site, day, sampling_time, sampling_day) %>%
                     dplyr::summarise(depth = mean(depth)), by = join_by(site, day, sampling_time)) %>%
  dplyr::arrange(site, sampling_day, sampling_time, parameter) %>%
  dplyr::mutate(replicate = rep(c(rep(1:3, 4),
                                  rep(1:3, 4),
                                  rep(1:3, 4)), 10)) %>%
  tidyr::pivot_wider(., id_cols = c("site", "sampling_day", "sampling_time", "replicate", "depth"),
                     names_from = "parameter",
                     values_from = "max") %>%
  # dplyr::rowwise(.) %>%
  # dplyr::mutate(PAR2 = (lux/54)) %>%
  # dplyr::mutate(PAR2 = PAR2/exp(-0.13*(depth - 23))) %>%
  # dplyr::ungroup(.) %>%
  dplyr::select(-replicate) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  dplyr::mutate(replicate = rep(c(sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3),
                                  sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3),
                                  sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3)), 2)) %>%
  dplyr::arrange(site, sampling_day, sampling_time, replicate) %>%
  droplevels


#if you did not calculate PAR for each specific depth and sampling_time window at each site and time, subset for the dawn and peak photo windows:
{
  # usvi_hobo_light_temp_expanded.df <- usvi_hobo_light_temp %>%
  #   dplyr::filter(grepl(paste0(c("2021-01-22", "2021-01-23", "2021-01-24", "2021-01-25", "2021-01-26"), collapse = "|"), day)) %>%
  #   dplyr::select(site, date_ast, day, time, temp, lux, lumens) %>%
  #   dplyr::left_join(., metabolomics_sample_metadata %>%
  #                      dplyr::mutate(day = lubridate::ymd(sampling_date)) %>%
  #                      dplyr::group_by(site, day) %>%
  #                      dplyr::summarise(depth = mean(depth))) %>%
  #   dplyr::group_by(site, day) %>%
  #   tidyr::fill(depth, .direction = "downup") %>%
  #   dplyr::mutate(PAR = (lux/54)) %>%
  #   dplyr::mutate(PAR = PAR/exp(-0.13*(depth - 23))) %>%
  #   tidyr::pivot_longer(., cols = c("temp", "lumens", "lux", "PAR"),
  #                       # tidyr::pivot_longer(., cols = c("temp", "lumens", "lux"),
  #                       names_to = "parameter",
  #                       values_to = "value") %>%
  #   dplyr::group_by(site) %>%
  #   dplyr::slice_head(prop = 0.9) %>%
  #   droplevels
  # 
  # # temp_df <- bind_rows(
  # usvi_light_temp_metadata.df <- bind_rows(
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-22 05:00:00"), ymd_hms("2021-01-22 08:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>% dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>% #do you want to randomly pick the top 3 in each?
  #      dplyr::mutate(sampling_day = "Day1") %>%
  #      dplyr::mutate(sampling_time = "dawn") %>%
  #      droplevels), 
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-22 11:00:00"), ymd_hms("2021-01-22 14:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day1") %>%
  #      dplyr::mutate(sampling_time = "peak_photo") %>%
  #      droplevels),
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-23 05:00:00"), ymd_hms("2021-01-23 08:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day2") %>%
  #      dplyr::mutate(sampling_time = "dawn") %>%
  #      droplevels), 
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-23 11:00:00"), ymd_hms("2021-01-23 14:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day2") %>%
  #      dplyr::mutate(sampling_time = "peak_photo") %>%
  #      droplevels),
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-24 05:00:00"), ymd_hms("2021-01-24 08:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day3") %>%
  #      dplyr::mutate(sampling_time = "dawn") %>%
  #      droplevels), 
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-24 11:00:00"), ymd_hms("2021-01-24 14:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day3") %>%
  #      dplyr::mutate(sampling_time = "peak_photo") %>%
  #      droplevels),
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-25 05:00:00"), ymd_hms("2021-01-25 08:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day4") %>%
  #      dplyr::mutate(sampling_time = "dawn") %>%
  #      droplevels), 
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-25 11:00:00"), ymd_hms("2021-01-25 14:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day4") %>%
  #      dplyr::mutate(sampling_time = "peak_photo") %>%
  #      droplevels),
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-26 05:00:00"), ymd_hms("2021-01-26 08:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day5") %>%
  #      dplyr::mutate(sampling_time = "dawn") %>%
  #      droplevels), 
  #   (usvi_hobo_light_temp_expanded.df %>%
  #      dplyr::filter(date_ast %within% interval(ymd_hms("2021-01-26 11:00:00"), ymd_hms("2021-01-26 14:00:00"))) %>%
  #      dplyr::group_by(site, day, parameter) %>%
  #      # dplyr::summarise(max = max(value, na.rm = TRUE)) %>%
  #      dplyr::reframe(max = pmax.int(value, na.rm = TRUE)) %>%
  #      dplyr::slice_max(n = 3, by = c(site, day, parameter), order_by = max, with_ties = FALSE) %>%
  #      dplyr::mutate(sampling_day = "Day5") %>%
  #      dplyr::mutate(sampling_time = "peak_photo") %>%
  #      droplevels)
  # ) %>%
  #   dplyr::left_join(., metabolomics_sample_metadata %>%
  #                      dplyr::group_by(site, sampling_day, sampling_time) %>%
  #                      dplyr::summarise(depth = mean(depth)),
  #                    by = join_by(site, sampling_day, sampling_time)) %>%
  #   dplyr::arrange(site, sampling_day, sampling_time, parameter) %>%
  #   dplyr::mutate(replicate = rep(c(rep(1:3, 4),
  #                                   rep(1:3, 4),
  #                                   rep(1:3, 4)), 10)) %>%
  #   tidyr::pivot_wider(., id_cols = c("site", "sampling_day", "sampling_time", "replicate", "depth"),
  #                      names_from = "parameter",
  #                      values_from = "max") %>%
  #   dplyr::rowwise(.) %>%
  #   dplyr::mutate(PAR = (lux/54)) %>%
  #   dplyr::mutate(PAR = PAR/exp(-0.13*(depth - 23))) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::select(-replicate) %>%
  #   dplyr::arrange(site, sampling_time, sampling_day) %>%
  #   dplyr::mutate(replicate = rep(c(sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3),
  #                                   sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3),
  #                                   sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3), sample(3, 3)), 2)) %>%
  #   dplyr::arrange(site, sampling_day, sampling_time, replicate) %>%
  #   droplevels
}


temp_g <- print(
  ggplot(data = usvi_hobo_light_temp_filtered.df,
         aes(x = date_ast, y = value, fill = site, group = interaction(site, parameter)))
  + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_grid(parameter~site, labeller = labeller(site = site_lookup),
               scales = "free_y")
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

# Global labellers/renamers -----------------------------------------------

model_dispersion_lookup <- data.frame(v1 = c("manual_dispersion", 
                                             "manual_dawn", "manual_photo", "auto_dawn", "auto_photo",
                                             "auto_dispersion"),
                                      label = c("All LB seagrass, manual dispersion",
                                                "Dawn LB seagrass, manual dispersion",
                                                "Peak photo LB seagrass, manual dispersion",
                                                "Dawn LB seagrass, auto dispersion",
                                                "Peak photo LB seagrass, auto dispersion",
                                                "All LB seagrass, auto dispersion")) %>%
  tibble::deframe(.)


group_labels_lookup <- c("Tektite_d" = "Tektite dawn",
                         "Tektite_p" = "Tektite peak photo",
                         "Yawzi_d" = "Yawzi dawn",
                         "Yawzi_p" = "Yawzi peak photo",
                         "LB_seagrass_d" = "LB seagrass dawn",
                         "LB_seagrass_p" = "LB seagrass peak photo",
                         "Tektite" = "Tektite",
                         "Yawzi" = "Yawzi",
                         "LB_seagrass" = "LB seagrass")

group_labels_colors <- viridisLite::turbo(length(group_labels_lookup)) %>%
  setNames(., names(group_labels_lookup))

site_lookup <- data.frame(site = c("LB_seagrass", "Tektite", "Yawzi", "control_extraction", "control_pcr", "control_seq"),
                          label = c("Lameshur Bay seagrass", "Tektite Reef", "Yawzi Reef",
                                    "Control (DNA Extraction)", "Control (PCR)", "Control (Sequencing)")) %>%
  tibble::deframe(.)
site_colors <- pals::kelly(22)[6:(5+length(site_lookup))] %>%
  # site_colors <- viridisLite::cividis(n = length(site_lookup), direction = 1) %>%
  setNames(., names(site_lookup))
sampling_time_lookup <- data.frame(sampling_time = c("dawn", "peak_photo"),
                                   label = c("Dawn", "Peak photosynthesis")) %>%
  tibble::deframe(.)
sampling_time_colors <- pals::ocean.haline(n = length(sampling_time_lookup)) %>%
  setNames(., names(sampling_time_lookup))
sampling_day_lookup <- data.frame(sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5"),
                                  label = c("20210122", "20210123", "20210124", "20210125", "20210126")) %>%
  tibble::deframe(.)
sampling_day_colors <- pals::ocean.thermal(n = length(sampling_day_lookup)) %>%
  setNames(., names(sampling_day_lookup))

# usvi_genera_relabel <- usvi_prok_filled.taxa.df %>%
#   dplyr::mutate(across(everything(), ~stringr::str_replace_all(.x, " clade", ""))) %>%
#   dplyr::rowwise(.) %>%
#   dplyr::mutate(first = dplyr::case_when((Order == Genus) ~ Class,
#                                          .default = Order)) %>%
#   dplyr::mutate(first = dplyr::case_when((Class == Order) ~ NA,
#                                          (Class != Phylum) ~ Phylum,
#                                          grepl(paste0(c("Synechococcales"), collapse = "|"), first) ~ NA,
#                                          .default = first)) %>%
#   dplyr::mutate(second = dplyr::case_when((Order == Genus) ~ Order,
#                                           .default = Genus)) %>%
#   dplyr::mutate(taxa_label = dplyr::case_when(!is.na(first) ~ paste0(first, "; ", second),
#                                               .default = second)) %>%
#   dplyr::select(asv_id, taxa_label) %>%
#   tibble::deframe(.)


usvi_genera_relabel <- usvi_prok_asvs.taxa %>%
  dplyr::select(-sequence) %>%
  dplyr::mutate(final_taxonomy = across(c(Genus:Domain)) %>% purrr::reduce(coalesce)) %>%
  dplyr::mutate(taxonomy = dplyr::case_when(
    (is.na(Phylum)) ~ paste(Domain), #Domain;specific
    (is.na(Class)) ~ paste(Domain, final_taxonomy, sep = ";"), #Domain;specific
    (stringr::str_length(Class) < 5 | grepl("(^[0-9])", Class) | is.na(Order)) ~ paste(Phylum, final_taxonomy, sep = ";"),
    (grepl("SAR11 clade", Order)) ~ paste0(Class, ";", "SAR11 ", final_taxonomy), #Class;specific
    (!is.na(Order) & !(Order == Class)) ~ paste(Class, final_taxonomy, sep = ";"), #Class;specific
    (!is.na(Order)) ~ paste(Order, final_taxonomy, sep = ";"))) %>% #Order;specific
  dplyr::relocate(contains("_id"), taxonomy) %>%
  dplyr::select(-c(final_taxonomy)) %>%
  droplevels %>%
  # dplyr::arrange(taxonomy) %>%
  # dplyr::mutate(across(contains("_id"), ~factor(.x, levels = unique(.x)))) %>%
  # droplevels
  dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>% dplyr::mutate(taxonomy = paste0(asv_id, ": ", taxonomy)) %>% dplyr::select(asv_id, taxonomy) %>% tibble::deframe(.)

global_labeller <- labeller(
  # model2 = model_dispersion_lookup,
  model = model_dispersion_lookup,
  Combo = group_labels_lookup,
  site = site_lookup,
  sampling_day = sampling_day_lookup,
  sampling_time = sampling_time_lookup,
  asv_id = usvi_genera_relabel,
  .multi_line = TRUE,
  .default = label_wrap_gen(25, multi_line = TRUE)
  # .default = label_value
)

# Prepare data for permanova ----------------------------------------------

drop <- c("CINAR_BC_73")
keep <- c("sample_id", "metab_deriv_label", "sample_type", "site", "sampling_time", "sampling_day", "sampling_date", "depth", 
          "fcm_Prochlorococcus", "fcm_Synechococcus", "fcm_Picoeukaryotes", "fcm_Unpigmented_cells",
          "replicate", "grouping", "PAR", "lumens", "lux", "temp")

usvi_sus_metabolites_idx <- usvi_metabolomics_long.df %>%
  dplyr::arrange(LOD) %>%
  dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
  dplyr::arrange(metabolites) %>%
  dplyr::filter(is.na(LOD)) %>%
  droplevels
usvi_metabolomics.tbl <- usvi_metabolomics.df %>%
  dplyr::select(!c(usvi_sus_metabolites_idx[["metabolites"]])) %>%
  dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  t() %>%
  as.data.frame(.)


usvi_selected_metadata <- metabolomics_sample_metadata %>%
  dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::select(intersect(colnames(metabolomics_sample_metadata), keep)) %>%
  # dplyr::select(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, starts_with("fcm")) %>%
  # dplyr::select(-fcm_label) %>%
  dplyr::mutate(across(c(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site), ~factor(.x))) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(site_type = dplyr::case_when(grepl("LB", site) ~ "seagrass",
                                             grepl("Yawzi|Tektite", site) ~ "reef",
                                             .default = NA)) %>%
  dplyr::mutate(rownames = sample_id) %>% tibble::column_to_rownames(var = "rownames") %>%
  dplyr::mutate(grouping = interaction(site, sampling_time)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_id, replicate)) %>%
  dplyr::left_join(., usvi_light_temp_metadata.df, by = join_by(site, sampling_day, sampling_time, replicate)) %>%
  tidyr::drop_na(.) %>%
  droplevels



#for permanova input matrices, twy two options:
#the ASV table at genus-level
#then the ASV table

usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
  dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus) %>%
  dplyr::distinct(Domain, Phylum, Class, Order, Family, Genus, .keep_all = TRUE) %>%
  droplevels

usvi_sw_genus.tbl <- usvi_prok_asvs.df %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  otu_to_taxonomy(., usvi_prok_filled.taxa.df, level = "Genus") %>%
  dplyr::left_join(., usvi_sw_genus.taxa.df,
                   by = join_by(Domain, Phylum, Class, Order, Family, Genus)) %>%
  dplyr::relocate(asv_id) %>%
  dplyr::select(-c(Domain, Phylum, Class, Order, Family, Genus)) %>%
  droplevels %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.)


usvi_sw_asv.tbl <- usvi_prok_asvs.df %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  droplevels %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.) %>%
  dplyr::select(colnames(usvi_sw_genus.tbl))

meta.microb <- metabolomics_sample_metadata %>%
  dplyr::select(intersect(colnames(metabolomics_sample_metadata), keep)) %>%
  dplyr::left_join(., usvi_selected_metadata %>%
                     dplyr::select(sample_id, site, site_type, grouping, replicate, PAR, lumens, lux, temp), multiple = "all", relationship = "many-to-many") %>%
# meta.microb <- usvi_selected_metadata %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = colnames(usvi_sw_genus.tbl))) %>%
  dplyr::arrange(sample_id) %>%
  tidyr::drop_na(.) %>%
  dplyr::mutate(grouping2 = paste0(grouping, ".", sampling_day)) %>%
  dplyr::mutate(across(c(sample_id, sample_type, site, sampling_time, sampling_day, metab_deriv_label, site_type, grouping, grouping2), ~factor(.x))) %>%
  dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                grouping = fct_relevel(grouping, "LB_seagrass.dawn"),
                grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day1")) %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  droplevels

meta.metab <- metabolomics_sample_metadata %>%
  # dplyr::select(intersect(colnames(metabolomics_sample_metadata), keep)) %>%
  dplyr::left_join(., usvi_selected_metadata %>%
                     dplyr::select(sample_id, site, site_type, grouping, replicate, PAR, lumens, lux, temp), multiple = "all", relationship = "many-to-many") %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  tidyr::fill(PAR, lumens, lux, temp, .direction = "down") %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(metab_deriv_label %in% colnames(usvi_metabolomics.tbl)) %>%
  dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
  dplyr::mutate(grouping2 = paste0(grouping, ".", sampling_day)) %>%
  dplyr::select(sample_id, colnames(meta.microb)) %>%
  dplyr::mutate(metab_deriv_label = factor(metab_deriv_label, levels = colnames(usvi_metabolomics.tbl))) %>%
  droplevels %>%
  dplyr::arrange(metab_deriv_label) %>%
  tidyr::drop_na(.) %>%
  dplyr::mutate(across(c(sample_id, sample_type, site, sampling_time, sampling_day, metab_deriv_label, site_type, grouping, grouping2), ~factor(.x))) %>%
  dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                grouping = fct_relevel(grouping, "LB_seagrass.dawn"),
                grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day2")) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  droplevels


dist_usvi_asv.d <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_genus.d <- usvi_sw_genus.tbl %>%
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)

dist_usvi_asv.mat <- as.matrix(dist_usvi_asv.d)
dist_usvi_genus.mat <- as.matrix(dist_usvi_genus.d)

dist_usvi_metab.d <- usvi_metabolomics.tbl %>%
  dplyr::select(rownames(meta.metab)) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  t() %>%
  vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_metab.mat <- as.matrix(dist_usvi_metab.d)

#stat checking:
{
  
# #which transformation to use?
# 
# #check homogeneity via anova(betadisper):
# #for asvs
# #site:sampling_time: Pr(>F) = 1.401e-05
# #site: Pr(>F) =1.704e-07
# #sampling_time: Pr(>F) = 0.2268
# #sampling_day: 0.5866
# 
# vegan::betadisper(dist_usvi_asv.d, type = "median",
#                   # interaction(meta.microb$site, meta.microb$sampling_time)) %>%
#   meta.microb$sampling_day) %>%
#   anova(.)
# 
# #for metabolomics
# #site:sampling_time: Pr(>F) = 0.008554
# #site: Pr(>F) =0.006649
# #sampling_time: Pr(>F) = 0.5337
# #sampling_day:  Pr(>F) = 0.9971
# 
# vegan::betadisper(dist_usvi_metab.d, type = "median",
#                   # interaction(meta.metab$site, meta.metab$sampling_time)) %>%
#                   meta.metab$sampling_day) %>%
#   anova(.)


# #visualizing dispersions in the ASV matrix based on site*sampling_time:
# betadisp_dist_usvi_asv.df <- vegan::betadisper(dist_usvi_asv.d, type = "median",
#                                                interaction(usvi_selected_metadata$site, usvi_selected_metadata$sampling_time)) %>%
#                                       # usvi_selected_metadata$sampling_time) %>%
#   # purrr::pluck("distances") %>% tibble::enframe(value = "dissimilarity", name = NULL) %>%
#   purrr::pluck("distances") %>%  tibble::enframe(value = "dissimilarity", name = "sample_id") %>% droplevels %>%
#   dplyr::left_join(usvi_selected_metadata %>%
#               dplyr::select(sample_id, site, sampling_time)) %>%
#   # tibble::column_to_rownames(var = "sample_id") %>%
#   droplevels
# 
# print(
#   ggplot(data = betadisp_dist_usvi_asv.df)
#   + geom_violin(aes(x = site, y = dissimilarity, fill = sampling_time), color = "white", draw_quantiles = c(0.25, 0.5, 0.75))
#   + theme_dark()
#   + scale_shape_manual(name = "Sampling site and time",
#                        values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#   + scale_fill_manual(name = "Sampling time",
#                       values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                override.aes = list(color = "black", shape = 21, size = 3)),
#            color = "none")
# )
  
}

# Comparing physicochemical and cotninuous variables using manova ---------
# 
# usvi_physicochem.tbl <- meta.microb %>%
#   dplyr::select(depth, PAR, temp)
# usvi_fcm.tbl <- meta.microb %>%
#   dplyr::select(starts_with("fcm_"))
# 
# 
# # #using manova doesn't seem as suitable:
# # temp_manova_fit <- manova(as.matrix(usvi_physicochem.tbl) ~ meta.microb$site*meta.microb$sampling_time )
# # summary.aov(temp_manova_fit)
# 
# 
# # #MRPP evaluates whether there is a significant difference between two or more GROUPs of sampling units
# # #the continuous data can be the input, but the variables to test are categorical
# # 
# # mrpp.lighttemp <- mrpp(as.matrix(usvi_physicochem.tbl), interaction(meta.microb$site, meta.microb$sampling_time), distance = "bray")
# # meandist(vegdist(usvi_physicochem.tbl, "bray"), interaction(meta.microb$site, meta.microb$sampling_time))
# # 
# # mrpp.fcm <- mrpp(as.matrix(usvi_fcm.tbl), interaction(meta.microb$site, meta.microb$sampling_time), distance = "bray")
# # meandist(vegdist(usvi_fcm.tbl, "bray"), interaction(meta.microb$site, meta.microb$sampling_time))
# # 
# # mrpp.asv <- mrpp(dist_usvi_asv.d, interaction(meta.microb$site, meta.microb$sampling_time), distance = "bray")
# # meandist(dist_usvi_asv.d, interaction(meta.microb$site, meta.microb$sampling_time))
# 
# 
# ##trying nested anova a la: https://rcompanion.org/rcompanion/d_07.html
# # library(nlme)
# # library(multcomp)
# # 
# # model_lme <- lme(fcm_Prochlorococcus ~ (site*sampling_time)/sampling_day,
# #                  random = ~1|sampling_day, data = meta.microb, method = "REML")
# # # K <- diag(length(coef(model_lme)))[-1,]
# # # rownames(K) <- names(coef(model_lme))[-1]
# # K <- diag(length(coef(model_lme)))
# # rownames(K) <- names(coef(model_lme))
# # 
# # anova.lme(model_lme,
# #           type="sequential",
# #           adjustSigma = FALSE)
# # model_lme_posthoc <- glht(model_lme, linfct = K)
# # model_lme_summary <- summary(model_lme_posthoc,
# #               test=adjusted("fdr"))
# # 
# # cld(model_lme_posthoc,
# #     level=0.05,
# #     decreasing=TRUE)
# 
# aov(fcm_Prochlorococcus ~ site*sampling_time + sampling_day/sampling_time, data = meta.microb) %>% summary
# #one at a time: calculate Tukey HSD from the flow counts and other continuous data, compared to site*sampling_time
# TukeyHSD(aov(fcm_Prochlorococcus ~ site*sampling_time + (site:sampling_day)/sampling_time, data = meta.microb))
# TukeyHSD(aov(fcm_Prochlorococcus ~ (site/sampling_time) + (sampling_day/sampling_time), data = meta.microb))
# TukeyHSD(aov(fcm_Prochlorococcus ~ site + site*(sampling_day/sampling_time), data = meta.microb))
# 
# #for all the continuous variables:
# continuous_idx <- c("depth", "fcm_Prochlorococcus", "fcm_Synechococcus", "fcm_Picoeukaryotes", "fcm_Unpigmented_cells",
#                     "PAR", "lumens", "lux", "temp")
# 
# temp_list <- NULL
# for(i in continuous_idx){
# # for(i in continuous_idx[1:2]){
#   var1 <- rlang::parse_expr(i)
#   design_formula <- paste0(var1, "~", "site*sampling_time*sampling_day")
#   design_formula <- rlang::parse_expr(design_formula)
#   temp_tuk <- list(TukeyHSD(aov(eval(design_formula), data = meta.microb)))
#   temp_list <- c(temp_list, temp_tuk) %>%
#     setNames(., c(names(temp_list), var1))
# }
# usvi_cts_param.list <- temp_list
# 


# Permanova via adonis2 ---------------------------------------------------
#please see here for why we should not include continuous variables in adonis2's implementation of permanova:
#https://learninghub.primer-e.com/books/should-i-use-primer-or-r/chapter/3-permanova-vs-adonis2-in-r/export/html

#we have a Crossed Model, where Site has 3 levels, Sampling_time has 2 levels, and Sampling_day has 5 levels

# #testing out adonis2's ability to handle strata:
# dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
# Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
# Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
# Y <- data.frame(Agropyron, Schizachyrium)
# adonis2(Y ~ NO3, data = dat, permutations = 199)
# ## Correct with strata
# with(dat, adonis2(Y ~ NO3, data = dat, permutations = 199, strata = field))



adonis_asv.res <- vegan::adonis2(data = meta.microb, method = "bray", permutations = 999,
                                 # formula = dist_usvi_asv.d ~ site*sampling_time, #pr site:sampling_time = 0.226
                                 # formula = dist_usvi_asv.d ~ (site*sampling_time) + site:sampling_time:sampling_day, #pr = 0.072 for site:sampling_time
                                 # formula = dist_usvi_asv.d ~ (site*sampling_time)/sampling_day, #site:sampling_time pr = 0.085
                                 # formula = dist_usvi_asv.d ~ site + sampling_time + site:sampling_time + site:sampling_time:sampling_day, #site:sampling_time pr = 0.081
                                 # formula = dist_usvi_asv.d ~ site/sampling_day/sampling_time, #site:sampling_time pr = 0.002
                                 formula = dist_usvi_asv.d ~ sampling_day*sampling_time,
                                   parallel = nthreads, by = "terms")
# summary(adonis_asv.res)
adonis_asv.res <- with(meta.microb, adonis2(data = meta.microb, method = "bray", permutations = 999,
                                            # dist_usvi_asv.d ~ site*sampling_time*sampling_day, strata = sampling_day, #Pr for site:sampling_time = 0.102; site:sampling_time:sampling_day: 0.134
                                            # dist_usvi_asv.d ~ site*sampling_time*sampling_day, strata = site, #Pr for site:sampling_time = 0.038; site:sampling_time:sampling_day: 0.051
                                            # dist_usvi_asv.d ~ sampling_time*sampling_day*site, strata = site, #Pr for sampling_time:site = 0.064; sampling_time:sampling_day:site: 0.053
                                            dist_usvi_asv.d ~ sampling_time/sampling_day, strata = site, #Pr for sampling_time:sampling_day = 0.152; sampling_time = 0.001, sampling_day = 0.001
                                            parallel = nthreads, by = "terms"))

adonis_asv.res


adonis_genus.res <- vegan::adonis2(data = meta.microb, method = "bray", permutations = 999,
                                 # formula = dist_usvi_asv.d ~ site*sampling_time,
                                 # formula = dist_usvi_asv.d ~ site*sampling_time*sampling_day + PAR*temp,
                                 # formula = dist_usvi_asv.d ~  fcm_Synechococcus + fcm_Prochlorococcus + fcm_Picoeukaryotes + fcm_Unpigmented_cells,
                                 formula = dist_usvi_genus.d ~ site*sampling_time*sampling_day,
                                 # formula = dist_usvi_asv.d ~ fcm_Synechococcus*fcm_Prochlorococcus*fcm_Picoeukaryotes*fcm_Unpigmented_cells,
                                 parallel = nthreads, by = "terms")
summary(adonis_genus.res)


# Kosher to use continuous variables in permanova? ------------------------

# 
# sim_fcm.tbl <- usvi_selected_metadata %>%
#   dplyr::select(starts_with("fcm_")) %>%
#   droplevels
# 
# dummy_metadata <- usvi_selected_metadata %>%
#   # dplyr::left_join(., metadata %>%
#   #                    dplyr::select(sample_id, replicate)) %>%
#   # dplyr::left_join(., usvi_light_temp_metadata.df, by = join_by(site, sampling_day, sampling_time, replicate)) %>%
#   # dplyr::select(sample_id, starts_with("fcm_"), PAR, temp) %>%
#   tibble::column_to_rownames(var = "sample_id") %>%
#   droplevels
# 
# # 
# # temp_idx <- usvi_selected_metadata %>%
# #   dplyr::select(fcm_Prochlorococcus) %>% tibble::deframe(.)
# # 
# # #which transformation to use?
# # dist_usvi_asv.df <- vegan::betadisper(dist_usvi_asv.d, type = "median",
# #                                       usvi_selected_metadata$grouping) %>%
# #   # purrr::pluck("distances") %>% tibble::enframe(value = "dissimilarity", name = NULL) %>%
# #   purrr::pluck("distances") %>%  tibble::enframe(value = "dissimilarity", name = "sample_id") %>% droplevels %>%
# #   dplyr::left_join(usvi_selected_metadata %>%
# #               dplyr::select(sample_id, contains("fcm_"))) %>%
# #   tibble::column_to_rownames(var = "sample_id")
# #   
# # 
# # print(
# #   ggplot(data = dist_usvi_asv.df)
# #       + geom_histogram(aes(x = fcm_Prochlorococcus)))
# # 
# # temp_df <- dist_usvi_asv.df %>%
# #   dplyr::select(contains("fcm_")) %>%
# #   apply(., 1, function(x) car::bcPower(x, lambda = 0.5)) %>%
# #   # apply(., 2, function(x) MASS::boxcox(x, plotit = FALSE))  %>%
# #   t() %>% as.data.frame(.) %>% 
# #   # dplyr::reframe(across(contains("fcm_"), ~quantile(.x, probs = seq(0, 1, 0.25)))) %>%
# #   droplevels
# # 
# # MASS::boxcox(dissimilarity ~ fcm_Synechococcus + fcm_Prochlorococcus + fcm_Picoeukaryotes + fcm_Unpigmented_cells, 
# #              plotit = FALSE,
# #              data = dist_usvi_asv.df)
# # 
# # plot(transform_boxcox(), xlim = c(0, 10))
# # plot(transform_boxcox(1), xlim = c(0, 10))
# 
# g1 <- (ggplot(data = dummy_metadata)
#       # + geom_point(aes(x = PAR, y = temp, shape = site, fill = sampling_time))
#       # + geom_point(aes(x = PAR, y = fcm_Prochlorococcus, shape = site, fill = sampling_time))
#       # + geom_point(aes(x = PAR, y = fcm_Synechococcus, shape = site, fill = sampling_time))
#       + geom_point(aes(x = temp, y = fcm_Prochlorococcus, shape = site, fill = sampling_time), size = 3)
#       # + scale_x_continuous(transform = "log10", name = "PAR")
#       # + scale_y_continuous(transform = "log10", name = "cell counts")
#       + scale_shape_manual(name = "Sampling site and time",
#                            values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#       + scale_fill_manual(name = "Sampling time",
#                           values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#       + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                    override.aes = list(color = "black", shape = 21, size = 3)),
#                color = "none")
#       )
# g2 <-  (ggplot(data = dummy_metadata)
#         # + geom_point(aes(x = PAR, y = temp, shape = site, fill = sampling_time))
#         # + geom_point(aes(x = PAR, y = fcm_Prochlorococcus, shape = site, fill = sampling_time))
#         # + geom_point(aes(x = PAR, y = fcm_Synechococcus, shape = site, fill = sampling_time))
#         + geom_point(aes(x = PAR, y = fcm_Prochlorococcus, shape = site, fill = sampling_time), size = 3)
#         # + scale_x_continuous(transform = "log10", name = "PAR")
#         # + scale_y_continuous(transform = "log10", name = "cell counts")
#         + scale_shape_manual(name = "Sampling site and time",
#                              values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#         + scale_fill_manual(name = "Sampling time",
#                             values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#         + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                      override.aes = list(color = "black", shape = 21, size = 3)),
#                  color = "none")
# )
# g3 <- (ggplot(data = dummy_metadata)
#        # + geom_point(aes(x = PAR, y = temp, shape = site, fill = sampling_time))
#        # + geom_point(aes(x = PAR, y = fcm_Prochlorococcus, shape = site, fill = sampling_time))
#        # + geom_point(aes(x = PAR, y = fcm_Synechococcus, shape = site, fill = sampling_time))
#        + geom_point(aes(x = temp, y = PAR, shape = site, fill = sampling_time), size = 3)
#        # + scale_x_continuous(transform = "log10", name = "PAR")
#        # + scale_y_continuous(transform = "log10", name = "cell counts")
#        + scale_shape_manual(name = "Sampling site and time",
#                             values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#        + scale_fill_manual(name = "Sampling time",
#                            values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#        + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                     override.aes = list(color = "black", shape = 21, size = 3)),
#                 color = "none")
# )
# 
# gpatch <- g1 + g2 + g3 + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Continuous data",
#                                                                                                  tag_levels = "A")
# ggsave(paste0(projectpath, "/", "modeled_PAR_fcmPro-", Sys.Date(), ".png"),
#        gpatch,
#        width = 16, height = 8, units = "in")
# 
# # lm(fcm_Prochlorococcus ~ lux*temp, dummy_metadata) %>%
# #   summary(.)
# # lm(fcm_Prochlorococcus ~ lumens*temp, dummy_metadata) %>%
# #   summary(.)
# lm(fcm_Prochlorococcus ~ PAR*temp, dummy_metadata) %>%
#   summary(.)
# 
# lm(fcm_Synechococcus ~ PAR*temp, dummy_metadata) %>%
#   summary(.)
# # lm(fcm_Synechococcus ~ lumens*temp, dummy_metadata) %>%
# #   summary(.)
# 
# lm(fcm_Prochlorococcus ~ site*sampling_time, dummy_metadata) %>% residuals() %>% hist(., col="darkgray")
# 
# # adonis_asv.sim_res <- vegan::adonis2(data = dummy_metadata, method = "bray", permutations = 999,
# #                                      # formula = sim_fcm.tbl[, 1] ~  PAR, #Pro counts are correlated and significant with PAR
# #                                  # formula = sim_fcm.tbl[, 1] ~  PAR*temp, #Pro counts are correlated and significant with PAR, temp, and PAR:temp
# #                                  # formula = sim_fcm.tbl ~ site, #none of the flow counts are signifiantly correlated with PAR; they are weakly sig with lumens, lux, and temp
# #                                  # formula = sim_fcm.tbl ~ site*sampling_time + PAR*temp, #none of the flow counts are signifiantly correlated with PAR, temp, or PAR:temp
# #                                  formula = sim_fcm.tbl ~ site*sampling_time, #none of the flow counts are signifiantly correlated with PAR, temp, or PAR:temp
# #                                  parallel = nthreads, by = "terms")
# # 
# # adonis_asv.sim_res
# 
# lm(fcm_Prochlorococcus ~ PAR, dummy_metadata) %>% summary(.)
# lm(fcm_Prochlorococcus ~ temp, dummy_metadata) %>% summary(.)
# lm(fcm_Prochlorococcus ~ site, dummy_metadata) %>% summary(.)
# lm_1 <- lm(fcm_Prochlorococcus ~ PAR, dummy_metadata)
# lm_2 <- lm(fcm_Prochlorococcus ~ temp, dummy_metadata)
# lm_3 <- lm(fcm_Prochlorococcus ~ PAR*temp, dummy_metadata)
# lm_4 <- lm(fcm_Prochlorococcus ~ site, dummy_metadata)
# # lm(fcm_Prochlorococcus ~ PAR*temp, dummy_metadata) %>% summary(.)
# 
# anova(lm_1)
# anova(lm_2)
# anova(lm_3)
# anova(lm_4)
# 
# #residuals should be normally distributed
# hist(residuals(lm_4),
#      col="darkgray")
# 
# adonis_asv.sim_res2 <- vegan::adonis2(data = dummy_metadata, method = "bray", permutations = 999,
#                                      # formula = sim_fcm.tbl[, 1] ~  PAR, #Pro counts are correlated and significant with PAR
#                                      # formula = sim_fcm.tbl[, 1] ~  PAR*temp, #Pro counts are correlated and significant with PAR, temp, and PAR:temp
#                                      # formula = sim_fcm.tbl ~ site, #none of the flow counts are signifiantly correlated with PAR; they are weakly sig with lumens, lux, and temp
#                                      # formula = sim_fcm.tbl ~ site*sampling_time + PAR*temp, #none of the flow counts are signifiantly correlated with PAR, temp, or PAR:temp
#                                      formula = sim_fcm.tbl ~ site*sampling_time*sampling_day, #none of the flow counts are signifiantly correlated with PAR, temp, or PAR:temp
#                                      parallel = nthreads, by = "terms")
# adonis_asv.sim_res2



# Try PERMANOVA -----------------------------------------------------------



#try PERMANOVA package
library(PERMANOVA)

pn_usvi_asv.mat <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  PERMANOVA::IniTransform(., transform = 4) #column centering: columns are ASVs, remove the column means

pn_usvi_asv.bc <- pn_usvi_asv.mat %>%
  PERMANOVA::DistContinuous(., coef = "Bray_Curtis")

#for PERMANOVA,
#can't use in option 'group' one of the Effects' levels
pn_usvi_asv.res <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc, 
                                        group = meta.microb$grouping, #6 contrasts: site.sampling_time
                                        seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv.res

pn_usvi_asv_day.res <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc, 
                                        group = meta.microb$grouping2, #30 contrasts: site.sampling_time.sampling_day
                                        # group = meta.microb$grouping, #6 contrasts: site.sampling_time
                                        seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv_day.res


#going deeper into PERMANOVA package:

#try buildingthe contrast matrix from scratch:
contrast_mat_order <- meta.microb %>%
  dplyr::arrange(site, sampling_time) %>%
  dplyr::select(grouping) %>%
  dplyr::mutate(label = grouping) %>%
  dplyr::mutate(label = as.numeric(label)) %>%
  droplevels

treatments <- factor(contrast_mat_order[["label"]], labels = unique(contrast_mat_order[["grouping"]]))
contrasts(treatments) <- meta.microb[rownames(contrast_mat_order),] %>%
  dplyr::select(site, sampling_time, grouping) %>%
  dplyr::distinct(.) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::pivot_wider(., id_cols = NULL,
                     names_from = c("grouping"),
                     values_from = "value",
                     values_fill = 0) %>%
  dplyr::mutate(site = as.numeric(site) - 1,
                # grouping =as.numeric(grouping),
                sampling_time = as.numeric(sampling_time) -1) %>%
  droplevels %>%
  dplyr::distinct(.) %>%
  as.matrix(.)


#build 0 intercept model and use limma to make contrast matrix
mod_mat <- model.matrix(~0 + treatments)
colnames(mod_mat) <- levels(treatments)
rownames(mod_mat) <- rownames(contrast_mat_order)

library(limma)
cmtx <- limma::makeContrasts( "LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
                              "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
                              "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo",
                              levels= mod_mat )

pn_C <- t(cmtx) #this is so it is amenable to PERMANOVA functions

interaction_effects <- data.frame(label = c("LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
                                            "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
                                            "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo"),
                                  to_test = c("sampling_time", "sampling_time", "sampling_time",
                                              "site", "site", "site",
                                              "site", "site", "site")) %>%
  dplyr::arrange(rownames(pn_C)) %>%
  dplyr::mutate(relabel = factor(label),
                retest = factor(to_test)) %>%
  dplyr::mutate(retest = as.numeric(retest)) %>%
  dplyr::mutate(relabel = as.numeric(relabel))
# pn_Effects <- factor(interaction_effects[["relabel"]], labels = unique(interaction_effects[["label"]]))
pn_Effects <- factor(interaction_effects[["retest"]], labels = unique(interaction_effects[["to_test"]]))

data_mat <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(contrast_mat_order)) %>%
  t() %>%
  PERMANOVA::IniTransform(., transform = 4) %>% #column centering: columns are ASVs, remove the column means
  PERMANOVA::DistContinuous(., coef = "Bray_Curtis")

pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(data_mat,
                                         group = contrast_mat_order$grouping, 
                                         C = pn_C, 
                                         Effects = pn_Effects,
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
summary(pn_usvi_asv.res2)


#try MANOVA?
mn_usvi_asv.res <- PERMANOVA::MANOVA(data_mat[["Data"]], Group = contrast_mat_order$grouping, 
                  C = pn_C, 
                  Effects = pn_Effects,
                  InitialTransform = 4, Contrasts = TRUE)
summary(mn_usvi_asv.res)




# pn_usvi_asv.contrast <- tibble::tribble(
#   ~LB_seagrass, ~Tektite, ~Yawzi,
#   
# )
# pn_usvi_asv.contrast <- contr.sum(levels(meta.microb$grouping), contrasts = FALSE, sparse = FALSE)
# colnames(pn_usvi_asv.contrast) <- levels(meta.microb$grouping)
cH <- contr.sum(levels(meta.microb$grouping), contrasts = TRUE, sparse = FALSE)
# apply(cH, 2, sum)
# crossprod(cH)


pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_time", "sampling_day")], MaxOrderIter = 3)
#what does the contrast matrix look like?
# #pulling out the contrasts matrix doesn't work as well:
# pn_C <- pn_usvi_asv.contrast[["Contrasts"]]
# colnames(pn_C) <- levels(pn_usvi_asv.contrast[["Groups"]])
# rownames(pn_C) <- pn_usvi_asv.contrast[["Effects"]]
# apply(pn_C, 2, sum)
# # pn_Effects <- factor(c(1:6))
# # levels(pn_Effects) <- levels(pn_usvi_asv.contrast[["Effects"]])
# pn_Effects <- levels(pn_usvi_asv.contrast[["Effects"]])

pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
                                         group = meta.microb$grouping2, #site.sampling_time.sampling_day
                                         # group = meta.microb$site, #can't use in group one of the Effects levels
                                         C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
summary(pn_usvi_asv.res2)
# Effects 
# Explained     Residual df Num df Denom       F-exp p-value p-value adj.
# site                       -9.682771e+37 1.147236e+35      2       58  -24476.244   0.886        0.993
# sampling_time               7.077761e+37 1.147236e+35      1       58   35782.528   0.106        0.318
# sampling_day                4.293153e+38 1.147236e+35      4       58   54261.463   0.038        0.228
# site*sampling_time         -5.729072e+38 1.147236e+35      2       58 -144820.284   0.993        0.993
# site*sampling_day           8.990167e+37 1.147236e+35      8       58    5681.367   0.193        0.386
# sampling_time*sampling_day -8.095247e+37 1.147236e+35      4       58  -10231.639   0.867        0.993
# Total                       4.106122e+38 1.147236e+35     21       58    9885.249   0.067        0.067

pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_time")], MaxOrderIter = 3)
pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
                                         group = meta.microb$grouping, #site.sampling_time
                                         C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
# Effects 
# Explained     Residual df Num df Denom      F-exp p-value p-value adj.
# site               -2.840383e+38 1.384244e+39      2       82  -8.412945   0.563        0.658
# sampling_time      -3.995956e+38 1.384244e+39      1       82 -23.671288   0.658        0.658
# site*sampling_time -2.898833e+38 1.384244e+39      2       82  -8.586069   0.621        0.658
# Total              -9.735171e+38 1.384244e+39      5       82 -11.533863   0.845        0.845

#when using ConstructContrasts to build the effects levels, they seem only to be able to design pairs:
pn_usvi_asv.contrast[["Effects"]]
#sampling_time sampling_day site sampling_time*sampling_day sampling_time*site sampling_day*site


pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("grouping", "sampling_day")], MaxOrderIter = 3)
pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
                                         group = meta.microb$grouping2, 
                                         C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
# Effects 
# Explained     Residual df Num df Denom      F-exp p-value p-value adj.
# grouping              -5.784676e+38 1.147236e+35      5       58 -58490.336   0.968        0.968
# sampling_day           4.082913e+38 1.147236e+35      4       58  51604.217   0.042        0.105
# grouping*sampling_day  5.807885e+38 1.147236e+35     20       58  14681.252   0.070        0.105
# Total                  4.106122e+38 1.147236e+35     29       58   7158.284   0.067        0.067



pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("sampling_time", "sampling_day", "site")], MaxOrderIter = 3)
# #pulling out the contrasts matrix doesn't work as well:
# pn_C <- pn_usvi_asv.contrast[["Contrasts"]]
# colnames(pn_C) <- levels(pn_usvi_asv.contrast[["Groups"]])
# rownames(pn_C) <- pn_usvi_asv.contrast[["Effects"]]
# # pn_Effects <- factor(c(1:6))
# # levels(pn_Effects) <- levels(pn_usvi_asv.contrast[["Effects"]])
# pn_Effects <- levels(pn_usvi_asv.contrast[["Effects"]])
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          seed = 48105, nperm = 999, PostHoc = "fdr",
#                                          # C = pn_C, Effects = pn_Effects,
#                                          group = meta.microb$grouping)

contrasts(meta.microb$grouping, contrasts = FALSE)

pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
                                         group = meta.microb$grouping2,
                                         C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
# Effects 
# Explained     Residual df Num df Denom       F-exp p-value p-value adj.
# sampling_time              -3.875917e+38 1.147236e+35      1       58 -195951.925   0.984        0.984
# sampling_day                3.719104e+38 1.147236e+35      4       58   47006.019   0.053        0.152
# site                       -3.006977e+38 1.147236e+35      2       58  -76010.779   0.965        0.984
# sampling_time*sampling_day  4.008087e+38 1.147236e+35      4       58   50658.489   0.031        0.152
# sampling_time*site         -2.898182e+38 1.147236e+35      2       58  -73260.636   0.967        0.984
# sampling_day*site           3.384847e+38 1.147236e+35      8       58   21390.660   0.076        0.152
# Total                       4.106122e+38 1.147236e+35     21       58    9885.249   0.067        0.067



pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_time")], MaxOrderIter = 3)
# pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_day")], MaxOrderIter = 3)
# crossprod(pn_usvi_asv.contrast[["Contrasts"]])
# apply(pn_usvi_asv.contrast[["Contrasts"]], 2, sum)

pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
                                         group = meta.microb$grouping, 
                                         C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
                                         seed = 48105, nperm = 999, PostHoc = "fdr")

pn_usvi_asv.res <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
                                         group = meta.microb$sampling_time, 
                                        # C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
                                         seed = 48105, nperm = 999, PostHoc = "fdr")


summary(pn_usvi_asv.res2)


# mn_usvi_asv.res <- PERMANOVA::MANOVA(Y = pn_usvi_asv.bc[["D"]], InitialTransform = 4,
#                      Group = meta.microb$grouping, 
#                      C = pn_usvi_asv.contrast[["Contrasts"]], 
#                      # M = 
#                      Contrasts = TRUE)
# summary(mn_usvi_asv.res)


#construct a matrix for contrasts:
# treatments <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6),labels=c("LB_seagrass.dawn", "LB_seagrass.peak_photo", "Tektite.dawn", "Tektite.peak_photo", "Yawzi.dawn", "Yawzi.peak_photo"))
# contrasts(treatments) <- cbind(site=c(0,0,1,1,2,2), sampling_time=c(0,1,0,1,0,1), `LB_seagrass.peak_photo` = c(0,1,0,0,0,0), `Tektite.peak_photo`=c(0,0,0,1,0,0), `Yawzi.peak_photo`=c(0,0,0,0,0,1))
# mod_mat <- model.matrix(~treatments)
# colnames(mod_mat) <- c("Intercept","site","sampling_time","LB_seagrass.peak_photo","Tektite.peak_photo","Yawzi.peak_photo")
# # rownames(mod_mat) <- as.character(treatments) #this is the order of the rows, but assigning rownames is not helpful

#include the days, not that it really matters:
# treatments <- factor(c(rep(1, 5),rep(2, 5),rep(3, 5),rep(4, 5),rep(5, 5),rep(6,5)), labels=c("LB_seagrass.dawn", "LB_seagrass.peak_photo", "Tektite.dawn", "Tektite.peak_photo", "Yawzi.dawn", "Yawzi.peak_photo"))
# contrasts(treatments) <- cbind(site=c(0,0,1,1,2,2), sampling_time=c(0,1,0,1,0,1), `LB_seagrass.peak_photo` = c(0,1,0,0,0,0), `Tektite.peak_photo`=c(0,0,0,1,0,0), `Yawzi.peak_photo`=c(0,0,0,0,0,1))

# treatments <- factor(as.numeric(meta.microb[["grouping"]]), labels = unique(meta.microb[["grouping"]]))
contrast_mat_order <- meta.microb %>%
  dplyr::arrange(site, sampling_time) %>%
  dplyr::select(grouping) %>%
  dplyr::mutate(label = grouping) %>%
  dplyr::mutate(label = as.numeric(label)) %>%
  droplevels
treatments <- factor(contrast_mat_order[["label"]], labels = unique(contrast_mat_order[["grouping"]]))

contrasts(treatments) <- meta.microb[rownames(contrast_mat_order),] %>%
# contrasts(treatments) <- meta.microb %>%
#   tibble::rownames_to_column(var = "sample_id") %>%
  # dplyr::arrange(site, sampling_time, sample_id) %>%
  dplyr::select(site, sampling_time, grouping) %>%
  dplyr::distinct(.) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::pivot_wider(., id_cols = NULL,
                     names_from = c("grouping"),
                     values_from = "value",
                     values_fill = 0) %>%
  dplyr::mutate(site = as.numeric(site) - 1,
                # grouping =as.numeric(grouping),
                sampling_time = as.numeric(sampling_time) -1) %>%
  droplevels %>%
  dplyr::distinct(.) %>%
  as.matrix(.)


#0 intercept
mod_mat <- model.matrix(~0 + treatments)
colnames(mod_mat) <- levels(treatments)
rownames(mod_mat) <- rownames(contrast_mat_order)

library(limma)
cmtx <- limma::makeContrasts( "LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
                              "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
                              "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo",
                              levels= mod_mat )
# data_mat <- dist_usvi_asv.mat[rownames(contrast_mat_order),]
data_mat <- usvi_sw_asv.tbl %>%
  dplyr::select(rownames(contrast_mat_order)) %>%
  apply(., 2, function(x) log2(x + 1))

fit <- limma::lmFit(data_mat, mod_mat) %>%
# fit <- limma::lmFit(data_mat, ) %>%
  limma::contrasts.fit(., contrasts=cmtx, coefficients=NULL)
fit2 <- eBayes(fit)
topTable(fit2,coef=2)
volcanoplot(fit2,coef=2,highlight=2)

limma.res <- decideTests(fit2, method="nestedF", adjust.method="fdr", p.value = 0.05)
summary(limma.res)



pn_Effects <- factor(colnames(cmtx)) %>%as.numeric
levels(pn_Effects) <- colnames(cmtx)
# pn_Effects <- factor(c(1, 2, 3))
# levels(pn_Effects) <- c("site", "sampling_time", "interaction")

data_mat <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(contrast_mat_order)) %>%
  t() %>%
  PERMANOVA::IniTransform(., transform = 4) %>% #column centering: columns are ASVs, remove the column means
  PERMANOVA::DistContinuous(., coef = "Bray_Curtis")

cmtx <- ConstructContrasts(meta.microb[rownames(contrast_mat_order), c("site", "sampling_time")], MaxOrderIter=3)
pn_usvi_asv.res3 <- PERMANOVA::PERMANOVA(data_mat,
                                         group = contrast_mat_order$grouping, 
                                         C = cmtx[["Contrasts"]], Effects=cmtx[["Effects"]],
                                         # C = t(cmtx), 
                                         # Effects = pn_Effects,
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
summary(pn_usvi_asv.res3)




# #construct a 5row x 6col matrix for each contrast to compare:
# mod_mat <- model.matrix(~0 + levels(meta.microb$grouping))
# colnames(mod_mat) <- c("LB_seagrass.dawn", "LB_seagrass.peak_photo", "Tektite.dawn", "Tektite.peak_photo", "Yawzi.dawn", "Yawzi.peak_photo")
# # mod_mat <- model.matrix(~levels(meta.microb$grouping))
# # colnames(mod_mat) <- c("Intercept", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo", "Yawzi.peak_photo")
# 
# rownames(mod_mat) <- c("LB_seagrass.dawn", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo", "Yawzi.peak_photo")
# 
# c("LB_seagrass.dawn", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo", "Yawzi.peak_photo")
# cH <- tibble::tribble(
#   ~"LB_seagrass.dawn", ~"Tektite.dawn", ~"Yawzi.dawn", ~"LB_seagrass.peak_photo", ~"Tektite.peak_photo", ~"Yawzi.peak_photo",
#   1, 0, 0, 0, 0, 0,
#   0, 1, 0, 0, 0, 0,
#   0, 0, 1, 0, 0, 0,
#   0, 0, 0, 1, 0, 0,
#   0, 0, 0, 0, 1, 0
# ) %>%
#   as.matrix(.)
# rownames(cH) <- c("LB_seagrass.dawn", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo")
# 
# # 
# pn_usvi_asv.contrast <- contr.treatment(levels(meta.microb$grouping), base = 1, contrasts = TRUE, sparse = FALSE)

# mod_mat <- model.matrix(~0 + site + sampling_time + sampling_day + site:sampling_time:sampling_day,
#                         data = meta.microb)
# if(length(which(colSums(mod_mat) == 0)) > 0){ #none of these should be 0--if so, you need to fix your design/metadata coding
#   idx <- names(colSums(mod_mat))[which(colSums(mod_mat) == 0)]
#   idx <- cli::cli_vec(idx)
#   cli::cli_alert_warning("One or more groups in the comparison does not have enough replicates:", wrap = TRUE)
#   cli::cli_bullets_raw(idx)
# }

# condition_filter <- list(`site` = c("LB_seagrass", "Yawzi", "Tektite"),
#                          `sampling_time` = c("dawn", "peak_photo"),
#                          `sampling_day` = c("Day1", "Day2", "Day3", "Day4", "Day5"))
# 
# condition_idx_names <- purrr::map(condition_filter,
#                                   ~unlist(., recursive = FALSE)) %>%
#   purrr::reduce(interaction, sep = "_") %>%
#   levels %>%
#   as.character
# 
# #brute force: list the model vectors
# for(i in c("LB_seagrass", "Yawzi", "Tektite")){
#   namevar <- i
#   for(j in c("dawn", "peak_photo")){
#     lightvar <- j
#     if(is.null(dim(mod_mat[meta.microb$site == namevar & meta.microb$sampling_time == lightvar,]))){
#       temp_vec <- mod_mat[meta.microb$site == namevar & meta.microb$sampling_time == lightvar, ]
#     } else {
#       temp_vec <- colMeans(mod_mat[meta.microb$site == namevar & meta.microb$sampling_time == lightvar, ])
#     }
#     assign(paste0("contrast_", namevar, "_", lightvar), temp_vec, envir = .GlobalEnv)
#     rm(temp_vec)
#     rm(lightvar)
#   }
#   rm(namevar)
# }
# 
# contrast_id <- condition_idx_names %>%
#   as.data.frame() %>%
#   dplyr::mutate(baseline = "baseline_vec") %>%
#   dplyr::mutate(across(everything(), as.character))

# pn_usvi.effects <- factor(c(1,2,3))
# levels(pn_usvi.effects) <- c("site", "sampling_time", "sampling_day")




###tngent:
#using Limma to detect SDA ASVs between pairs

contrast_mat_order <- meta.microb %>%
  dplyr::arrange(site, sampling_time) %>%
  dplyr::select(grouping) %>%
  dplyr::mutate(label = grouping) %>%
  dplyr::mutate(label = as.numeric(label)) %>%
  droplevels
treatments <- factor(contrast_mat_order[["label"]], labels = unique(contrast_mat_order[["grouping"]]))

contrasts(treatments) <- meta.microb[rownames(contrast_mat_order),] %>%
  dplyr::select(site, sampling_time, grouping) %>%
  dplyr::distinct(.) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::pivot_wider(., id_cols = NULL,
                     names_from = c("grouping"),
                     values_from = "value",
                     values_fill = 0) %>%
  dplyr::mutate(site = as.numeric(site) - 1,
                sampling_time = as.numeric(sampling_time) -1) %>%
  droplevels %>%
  dplyr::distinct(.) %>%
  as.matrix(.)


#0 intercept
mod_mat <- model.matrix(~0 + treatments)
colnames(mod_mat) <- levels(treatments)
rownames(mod_mat) <- rownames(contrast_mat_order)

library(limma)
cmtx <- limma::makeContrasts( "LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
                              "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
                              "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo",
                              levels= mod_mat )

data_mat <- usvi_sw_asv.tbl %>%
  dplyr::select(rownames(contrast_mat_order)) %>%
  apply(., 2, function(x) log2(x + 1))

fit <- limma::lmFit(data_mat, mod_mat) %>%
  limma::contrasts.fit(., contrasts=cmtx, coefficients=NULL)
fit2 <- eBayes(fit)
topTable(fit2,coef=2)
volcanoplot(fit2,coef=2,highlight=2)

limma.res <- decideTests(fit2, method="nestedF", adjust.method="fdr", p.value = 0.05)
summary(limma.res)



