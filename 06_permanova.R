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
      dplyr::mutate(site = factor(site, levels = names(site_lookup))) %>%
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
  dplyr::mutate(site = factor(site, levels = names(site_lookup))) %>%
  droplevels %>%
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
# 
# 
g1_par <- print(
  ggplot(data = usvi_hobo_light_temp_filtered.df %>%
           dplyr::filter(grepl("PAR", parameter)),
         aes(x = date_ast, y = value, fill = site, color = site, group = interaction(site, parameter)))
  + geom_point(shape = 19, size = 1)
  + geom_line(show.legend = FALSE)
  + theme_bw()
  + facet_grid(site~., labeller = labeller(site = site_lookup),
               scales = "fixed")
  + scale_y_continuous(name = expression(paste("PAR (\U00B5 mol photons", ~m^-2 ~s^-1, ")")))
  + scale_discrete_manual(aesthetics = c("color"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
                          drop = TRUE)
  + scale_discrete_manual(aesthetics = c("fill"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
                          drop = TRUE)
  + theme(axis.text.x = element_text(angle = 90),
          panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          legend.key = element_blank(),
          legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
          legend.text = element_text(size = 12, colour = "grey30"))
)
g1_temp <- print(
  ggplot(data = usvi_hobo_light_temp_filtered.df %>%
           dplyr::filter(grepl("temp", parameter)),
         aes(x = date_ast, y = value, fill = site, color = site, group = interaction(site, parameter)))
  + geom_point(shape = 19, size = 1)
  + geom_line(show.legend = FALSE)
  + theme_bw()
  + facet_grid(site~., labeller = labeller(site = site_lookup),
               scales = "fixed")
  + scale_y_continuous(name = expression(paste("Temperature (˚C)")))
  + scale_discrete_manual(aesthetics = c("color"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
                          drop = TRUE)
  + scale_discrete_manual(aesthetics = c("fill"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
                          drop = TRUE)
  + theme(axis.text.x = element_text(angle = 90),
          panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          legend.key = element_blank(),
          legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
          legend.text = element_text(size = 12, colour = "grey30"))
)

gpatch <- g1_par + g1_temp + patchwork::plot_annotation(title = "Physicochemical parameters", subtitle ="Measured in USVI sites via HOBO loggers", tag_levels = "A")
gpatch
ggsave(paste0(projectpath, "/", "usvi_light_temp-", Sys.Date(), ".png"),
       gpatch,
       width = 10, height = 8, units = "in")


# Global labellers/renamers -----------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_global_labellers.R"))){
  cat("Loading project-wide labellers and lookups.")
  source(paste0(getwd(), "/", "00_global_labellers.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 

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
keep <- c("sample_id", "metab_deriv_label", "sample_type", "site", "sampling_time", "sampling_day", "sampling_date", "depth", "site_type",
          "fcm_Prochlorococcus", "fcm_Synechococcus", "fcm_Picoeukaryotes", "fcm_Unpigmented_cells",
          "replicate", "grouping", "PAR", "lumens", "lux", "temp")

#as a reminder, here are the sus metabolites without LODs reported:
#2'deoxyguanosine, HMP, adenosine, inosine, pyridoxine
if(!exists("usvi_sus_metabolites_idx", envir = .GlobalEnv)){
  usvi_sus_metabolites_idx <- usvi_metabolomics_long.df %>%
    dplyr::arrange(LOD) %>%
    dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
    dplyr::arrange(metabolites) %>%
    dplyr::filter(is.na(LOD)) %>%
    droplevels
}

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
  dplyr::select(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, starts_with("fcm")) %>%
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
  # tidyr::drop_na(.) %>% #removed sample Metab_306 due to lack of FCM measurements
  droplevels

#permanova at the ASV level:

#old: using genus-level
{
  # usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
  #   dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus) %>%
  #   dplyr::distinct(Domain, Phylum, Class, Order, Family, Genus, .keep_all = TRUE) %>%
  #   droplevels
  # 
  # usvi_sw_genus.tbl <- usvi_prok_asvs.df %>%
  #   dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  #   droplevels %>%
  #   tidyr::pivot_wider(., id_cols = "asv_id",
  #                      names_from = "sample_ID",
  #                      values_from = "counts",
  #                      values_fill = 0) %>%
  #   otu_to_taxonomy(., usvi_prok_filled.taxa.df, level = "Genus") %>%
  #   dplyr::left_join(., usvi_sw_genus.taxa.df,
  #                    by = join_by(Domain, Phylum, Class, Order, Family, Genus)) %>%
  #   dplyr::relocate(asv_id) %>%
  #   dplyr::select(-c(Domain, Phylum, Class, Order, Family, Genus)) %>%
  #   droplevels %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   as.data.frame(.)
}



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
  # dplyr::select(colnames(usvi_sw_genus.tbl)) %>%
  droplevels

# meta.microb <- metabolomics_sample_metadata %>%
#   dplyr::select(intersect(colnames(metabolomics_sample_metadata), keep)) %>%
#   dplyr::left_join(., usvi_selected_metadata %>%
#                      dplyr::select(sample_id, site, site_type, grouping, replicate, PAR, lumens, lux, temp), multiple = "all", relationship = "many-to-many") %>%
meta.microb <- usvi_selected_metadata %>%
    dplyr::select(intersect(colnames(usvi_selected_metadata), keep)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = colnames(usvi_sw_asv.tbl))) %>%
  dplyr::arrange(sample_id) %>%
  # tidyr::drop_na(.) %>%
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
  dplyr::mutate(metab_deriv_label = factor(metab_deriv_label, levels = colnames(usvi_metabolomics.tbl))) %>%
  dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
  dplyr::mutate(grouping2 = paste0(grouping, ".", sampling_day)) %>%
  dplyr::select(sample_id, colnames(meta.microb)) %>%
  dplyr::mutate(metab_deriv_label = factor(metab_deriv_label, levels = colnames(usvi_metabolomics.tbl))) %>%
  droplevels %>%
  dplyr::arrange(metab_deriv_label) %>%
  # tidyr::drop_na(.) %>%
  dplyr::mutate(across(c(sample_id, sample_type, site, sampling_time, sampling_day, metab_deriv_label, site_type, grouping, grouping2), ~factor(.x))) %>%
  dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                grouping = fct_relevel(grouping, "LB_seagrass.dawn"),
                grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day2")) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  droplevels
dist_usvi_metab.d <- usvi_metabolomics.tbl %>%
  dplyr::select(rownames(meta.metab)) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  t() %>%
  vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_metab.mat <- as.matrix(dist_usvi_metab.d)


dist_usvi_asv.d <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE)

dist_usvi_asv.mat <- as.matrix(dist_usvi_asv.d)


dist_usvi_asv_log2.d <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE)

dist_usvi_asv_log2.mat <- as.matrix(dist_usvi_asv_log2.d)


# dist_usvi_genus.d <- usvi_sw_genus.tbl %>%
#   dplyr::select(rownames(meta.microb)) %>%
#   t() %>%
#   vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE)

# dist_usvi_genus.mat <- as.matrix(dist_usvi_genus.d)


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
{
  

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
}

# Permanova via adonis2 ---------------------------------------------------
#please see here for why we should not include continuous variables in adonis2's implementation of permanova:
#https://learninghub.primer-e.com/books/should-i-use-primer-or-r/chapter/3-permanova-vs-adonis2-in-r/export/html

#we have a Crossed Model, where Site has 3 levels, Sampling_time has 2 levels, and Sampling_day has 5 levels
#with this data set we really ought to use type III SS
#but adonis only does type I SS.

#SS(total) for microbiomes: 1.49
(1/nrow(dist_usvi_asv.mat))*(sum((dist_usvi_asv.mat[lower.tri(dist_usvi_asv.mat, diag = FALSE)])^2))

#per-site Sum of Squares
dist_usvi_asv.ss <- meta.microb %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x]) %>%
  map(., ~sum((.x[lower.tri(.x, diag = FALSE)])^2)/nrow(.x))
#LB: 0.4070378
#Tektite: 0.1798429
#Yawzi: 0.1220277


#this is the final one we should use:
#in this output table, we want the results for:
#site
#site:sampling_time
#site:sampling_day
#site:sampling_time:sampling_day

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
                                 # strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time*site*sampling_day,
                                 parallel = nthreads, by = "terms"))
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1  0.11015 0.07360 26.1100 0.0001 ***
#   site                             2  0.78347 0.52345 92.8533 0.0001 ***
#   sampling_day                     4  0.20673 0.13812 12.2504 0.0001 ***
#   sampling_time:site               2  0.01104 0.00738  1.3089 0.3095    
# sampling_time:sampling_day       4  0.04129 0.02759  2.4467 0.0524 .  
# site:sampling_day                8  0.08502 0.05681  2.5192 0.0175 *  
#   sampling_time:site:sampling_day  8  0.01013 0.00677  0.3002 0.9451    
# Residual                        59  0.24891 0.16630                   
# Total                           88  1.49675 1.00000                 

#do samplign day before ssite:
with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
                                 # strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time*sampling_day*site,
                                 parallel = nthreads, by = "terms"))
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1  0.11015 0.07360 26.1100 0.0001 ***
#   sampling_day                     4  0.20072 0.13410 11.8943 0.0001 ***
#   site                             2  0.78948 0.52746 93.5656 0.0001 ***
#   sampling_time:sampling_day       4  0.04200 0.02806  2.4889 0.0469 *  
#   sampling_time:site               2  0.01033 0.00690  1.2244 0.3417    
# sampling_day:site                8  0.08502 0.05681  2.5192 0.0189 *  
#   sampling_time:sampling_day:site  8  0.01013 0.00677  0.3002 0.9451    
# Residual                        59  0.24891 0.16630                   
# Total                           88  1.49675 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

temp_res <- with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
                                             # strata = site,
                                             formula = dist_usvi_asv.d ~ sampling_time*sampling_day*site,
                                             parallel = nthreads, by = "terms"))
permanova_asv_estimate_variation.df <- temp_res %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::mutate(term = as.character(term)) %>%
  dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                      .default = NA))
permanova_asv_estimate_variation.df <- permanova_asv_estimate_variation.df %>%
  dplyr::mutate(permanova_asv_estimate_variation.df %>%
                     dplyr::filter(!grepl("Total", term)) %>%
                     dplyr::summarise(TotalEMS = sum(MS))) %>%
  dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
                                          .default = NA)) %>%
  dplyr::mutate(sq_root = MS^(1/2))
# term Df SumOfSqs      R2       F Pr(>F)      MS TotalEMS V_term sq_root
# [1,]    4  1  0.11015 0.07360 26.1100 0.0001 0.11015  0.58685 18.770 0.33189
# [2,]    2  4  0.20072 0.13410 11.8943 0.0001 0.05018  0.58685  8.551 0.22401
# [3,]    8  2  0.78948 0.52746 93.5656 0.0001 0.39474  0.58685 67.264 0.62828
# [4,]    5  4  0.04200 0.02806  2.4889 0.0509 0.01050  0.58685  1.789 0.10247
# [5,]    7  2  0.01033 0.00690  1.2244 0.3353 0.00517  0.58685  0.880 0.07187
# [6,]    3  8  0.08502 0.05681  2.5192 0.0172 0.01063  0.58685  1.811 0.10309
# [7,]    6  8  0.01013 0.00677  0.3002 0.9454 0.00127  0.58685  0.216 0.03559
# [8,]    1 59  0.24891 0.16630                0.00422  0.58685  0.719 0.06495
# [9,]    9 88  1.49675 1.00000                         0.58685             

#how does it look with log2(x+1) relative abundance as basis of diss matrix?
with(meta.microb, vegan::adonis2(data = meta.microb, 
                                 # method = "horn", 
                                 permutations = 9999,
                                 # strata = site,
                                 formula = dist_usvi_asv_log2.d ~ sampling_time*sampling_day*site,
                                 parallel = nthreads, by = "terms"))
# Df SumOfSqs       R2        F Pr(>F)    
# sampling_time                    1  0.06863  0.04266  36.7012 0.0001 ***
#   sampling_day                     4  0.18506  0.11504  24.7411 0.0001 ***
#   site                             2  1.22081  0.75892 326.4262 0.0001 ***
#   sampling_time:sampling_day       4  0.02061  0.01281   2.7552 0.0457 *  
#   sampling_time:site               2 -0.01211 -0.00753  -3.2375 1.0000    
# sampling_day:site                8  0.01191  0.00741   0.7964 0.5860    
# sampling_time:sampling_day:site  8  0.00337  0.00209   0.2250 0.9675    
# Residual                        59  0.11033  0.06859                    
# Total                           88  1.60861  1.00000                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 





#Site has a significant difference via P-value and F-ratio
#since it has the highest pseudo F-ratio, frst investigate site
#run pair-wise comparisons of Site: LB vs Tektite, LB vs Yawzi, Tektite vs Yawzi
#results: 


#per-site average dissimilarity for ASVs:

dist_usvi_asv.df <- meta.microb %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "site")

  
dist_usvi_asv.df %>%
  dplyr::group_by(site) %>%
  # dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
  #                  min_dist = min(dissimilarity, na.rm = TRUE),
  #                  max_dist = max(dissimilarity, na.rm = TRUE)) 
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
                   min_dist = min(dissimilarity, na.rm = TRUE),
                   max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_asv.df[["dissimilarity"]]))))


# site        avg_dist min_dist max_dist avg_dist_rescaled min_dist_rescaled max_dist_rescaled
# <chr>          <dbl>    <dbl>    <dbl>             <dbl>             <dbl>             <dbl>
#   1 LB_seagrass   0.138   0.00266    0.482             0.286           0.00553             1    
# 2 Tektite       0.0837  0.00273    0.423             0.174           0.00568             0.878
# 3 Yawzi         0.0650  0.00179    0.376             0.135           0.00372             0.780

#is it different with log-transformed?

meta.microb %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv_log2.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "site") %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# site        avg_dist
# <chr>          <dbl>
#   1 LB_seagrass   0.119 
# 2 Tektite       0.0649
# 3 Yawzi         0.0369

#what about if we used bray-curtis instead of MH?
#still LB is higher than Tektite and higher than Yawzi

dist_usvi_asv_bc.mat <- usvi_sw_asv.tbl %>% 
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE) %>%
  as.matrix(.)

meta.microb %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv_bc.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "site") %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# # A tibble: 3 × 2
# site        avg_dist
# <chr>          <dbl>
#   1 LB_seagrass    0.320
# 2 Tektite        0.224
# 3 Yawzi          0.180


#approach 1: Kruskal-Wallis testing
{
  #Tektite vs Yawzi:
  dist_usvi_asv.df %>%
    dplyr::filter(!grepl("LB", site)) %>%
    kruskal.test(dissimilarity ~ site, .)
  # Kruskal-Wallis rank sum test
  # 
  # data:  dissimilarity by site
  # Kruskal-Wallis chi-squared = 22.837, df = 1, p-value = 1.763e-06
  
  #Yawzi vs LB:
  dist_usvi_asv.df %>%
    dplyr::filter(!grepl("Tektite", site)) %>%
    kruskal.test(dissimilarity ~ site, .)
  # Kruskal-Wallis rank sum test
  # 
  # data:  dissimilarity by site
  # Kruskal-Wallis chi-squared = 153.81, df = 1, p-value < 2.2e-16
  
  #Tektite vs LB:
  dist_usvi_asv.df %>%
    dplyr::filter(!grepl("Yawzi", site)) %>%
    kruskal.test(dissimilarity ~ site, .)
  # Kruskal-Wallis rank sum test
  # 
  # data:  dissimilarity by site
  # Kruskal-Wallis chi-squared = 77.509, df = 1, p-value < 2.2e-16
  
  #all 3:
  dist_usvi_asv.df %>%
    kruskal.test(dissimilarity ~ site, .)
  # data:  dissimilarity by site
  # Kruskal-Wallis chi-squared = 167.83, df = 2, p-value < 2.2e-16
  
  
  #which is identical to:
  
  # temp_list <- meta.microb %>%
  #   split(., f = .$site) %>%
  #   map(., ~.x %>%
  #         tibble::rownames_to_column(var = "sample_id") %>%
  #         dplyr::select(sample_id) %>%
  #         tibble::deframe(.)) %>%
  #   map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  #   map(., ~.x[lower.tri(.x, diag = FALSE)])
  # 
  # kruskal.test(temp_list)
  #Kruskal-Wallis chi-squared = 167.83, df = 2, p-value < 2.2e-16
  # 
  }


#interlude: KW testing metabolites
{
temp_list <- meta.metab %>% #metadata file where rownames are "CINAR_BC_*"
  split(., f = .$grouping) %>% #grouping looks like: "Yawzi.dawn"
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~usvi_metabolomics.tbl %>% #this is a wide dataframe where rownames are the metabolite names, and columsn are sample names "CINAR_BC_*"
        apply(., 2, function(x) log2(x + 1)) %>%
        as.data.frame(.) %>%
        dplyr::select(all_of(.x)) %>%
        tibble::rownames_to_column(var = "metabolite") %>%
        tidyr::pivot_longer(., cols = !c("metabolite"),
                            names_to = "metab_deriv_label",
                            values_to = "concentration")) %>%
  bind_rows(., .id = "grouping") %>%
  dplyr::left_join(., meta.metab %>%
                     tibble::rownames_to_column(var = "metab_deriv_label") %>%
                     dplyr::select(metab_deriv_label, site))
temp_kw_res <- temp_list %>%
  split(., f = .$metabolite) %>%
  map(., ~.x %>%
        split(., f = .$site) %>%
        map(., ~.x %>% 
              kruskal.test(concentration ~ grouping, .)))


temp_kw_res <- temp_list %>%
  split(., f = .$metabolite) %>%
  map(., ~.x %>%
        split(., f = .$site) %>%
        map(., ~.x %>% 
              kruskal.test(concentration ~ grouping, .))) %>%
  map_depth(., 2, ~.x %>%
              purrr::pluck("p.value")) %>%
  bind_rows(., .id = "metabolite") %>%
  tidyr::pivot_longer(., cols = !"metabolite",
                      names_to = "site.dawn_vs_afternoon",
                      values_to = "p.value")


q_cutoff <- p.adjust(temp_kw_res[["p.value"]], method = "fdr") %>%
  na.omit(.) %>%
  ashr::qval.from.lfdr(.) %>%
  quantile(., probs = q_value, na.rm = TRUE, names = FALSE, type = 7)
temp_kw_res <- temp_kw_res %>%
  dplyr::mutate(filtered_p.value = dplyr::case_when(p.value > 0.05 ~ NA,
                                               .default = p.value)) %>%
  dplyr::mutate(adj_p.value = dplyr::case_when(p.value > q_cutoff ~ NA,
                                           .default = p.value))
length(na.omit(temp_kw_res[["filtered_p.value"]]))
length(na.omit(temp_kw_res[["adj_p.value"]]))

  
# temp_list <- usvi_metabolomics.tbl %>%
#   dplyr::select(rownames(meta.metab)) %>%
#   tibble::rownames_to_column(var = "metabolite") %>%
#   split(., f = .$metabolite) %>%
#   map(., ~.x %>%
#         tibble::column_to_rownames(var = "metabolite"))

kruskal.test(temp_list)
}

#approach 2: pairwise t-tests
{
  # 
  # t.test(temp_list[[1]], temp_list[[2]])
  # # pairwise.t.test(c(temp_list[["LB_seagrass"]], temp_list[["Tektite"]]),
  # #                 c(rep("LB_seagrass", length(temp_list[["LB_seagrass"]])),
  # #                   rep("Tektite", length(temp_list[["Tektite"]]))),
  # #                 pool.sd = FALSE, paired = FALSE,
  # #                 p.adjust.method = "fdr")
  # # # LB_seagrass
  # # # Tektite <2e-16     
  # # # 
  # # # P value adjustment method: fdr 
  # 
  
  temp_res <- pairwise.t.test((dist_usvi_asv.df %>%
                                 split(., f = .$site) %>%
                                 map(., ~.x %>%
                                       dplyr::select(dissimilarity) %>%
                                       tibble::deframe(.)) %>%
                                 purrr::reduce(c)),
                              (dist_usvi_asv.df %>%
                                 split(., f = .$site) %>%
                                 map(., ~.x %>%
                                       dplyr::select(site) %>%
                                       tibble::deframe(.)) %>%
                                 purrr::reduce(c)),
                              pool.sd = FALSE, paired = FALSE, alternative = "two.sided",
                              p.adjust.method = "BH")
  
  
  # LB_seagrass Tektite
  # Tektite < 2e-16     -      
  #   Yawzi   < 2e-16     7.4e-05
  # 
  # P value adjustment method: fdr 
  
  
  t.test(dist_usvi_asv.df %>%
           dplyr::filter(grepl("Tektite", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         dist_usvi_asv.df %>%
           dplyr::filter(grepl("Yawzi", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         conf.level = 0.99)
  
  # t = 3.9811, df = 854.51, p-value = 7.442e-05
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.00948756 0.02793962
  # sample estimates:
  #   mean of x  mean of y 
  # 0.08370106 0.06498747 
  
  t.test(dist_usvi_asv.df %>%
           dplyr::filter(grepl("LB", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         dist_usvi_asv.df %>%
           dplyr::filter(grepl("Tektite", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.))
  # t = 8.8941, df = 739.58, p-value < 2.2e-16
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.04228695 0.06624246
  # sample estimates:
  #   mean of x  mean of y 
  # 0.13796577 0.08370106 
  
  t.test(dist_usvi_asv.df %>%
           dplyr::filter(grepl("LB", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         dist_usvi_asv.df %>%
           dplyr::filter(grepl("Yawzi", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.))
  # t = 12.434, df = 685.06, p-value < 2.2e-16
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.06145430 0.08450229
  # sample estimates:
  #   mean of x  mean of y 
  # 0.13796577 0.06498747 
  
}
# resample_depth <- dist_usvi_metab.df %>%
#   split(., f = .$site) %>%
#   map(., ~round(nrow(.x)/10)) %>%
#   purrr::reduce(mean)
# 
# length_v <- c("LB_seagrass", "Yawzi", "Tektite")
# ttest_res <- combn(length_v, 2) %>%
#   t() %>%
#   as.data.frame() %>%
#   setNames(., c("pair1", "pair2")) %>%
#   dplyr::mutate(p.value = NA)
# for(i in seq_len(nrow(ttest_res))){
#   var1 <- ttest_res[i, 1]
#   var2 <- ttest_res[i, 2]
#   temp_i_a <- dist_usvi_metab.df %>%
#     dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
#     dplyr::select(dissimilarity) %>%
#     tibble::deframe(.) %>% sample(., resample_depth, replace = FALSE)
#   temp_i_b <- dist_usvi_metab.df %>%
#     dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
#     dplyr::select(dissimilarity) %>%
#     tibble::deframe(.) %>% sample(., resample_depth, replace = FALSE)
#   ttest_res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
# }

#approach 3: bootstrap resampling the dissimilarity values in each site, and pairwise t-test them
#generalized function to grab dissimilarities
F_ttest_gn <- function(dataset, variable, f){
  #d is a dataframe containing all the measurements by site
  #variable is your factor you want to examine pairwise the levels
  variable <- rlang::parse_expr(variable)
  temp_vars <- unique(dataset[[variable]])
  
  resample_depth <- dataset %>%
    split(., f = .[[variable]]) %>%
    map(., ~round(nrow(.x)/10)) %>%
    purrr::reduce(mean)
  
  ttest_res <- combn(temp_vars, 2) %>%
    t() %>%
    as.data.frame() %>%
    setNames(., c("pair1", "pair2")) %>%
    dplyr::mutate(p.value = NA)
  for(i in seq_len(nrow(ttest_res))){
    var1 <- ttest_res[i, 1]
    var2 <- ttest_res[i, 2]
    temp_i_a <- dataset %>%
      dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
      dplyr::select(dissimilarity) %>%
      tibble::deframe(.) %>% sample(., resample_depth, replace = FALSE)
    temp_i_b <- dataset %>%
      dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
      dplyr::select(dissimilarity) %>%
      tibble::deframe(.) %>% sample(., resample_depth, replace = FALSE)
    ttest_res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
  }

  #list the resulting p-values from t-testing LB vs Yawzi, LB vs Tektite, and Yawzi vs Tektite
  ttest_res %>% dplyr::select(p.value) %>% tibble::deframe(.)
}

#testing it:
# F_ttest_gn(dataset=dist_usvi_asv.df, variable = "site")

dist_usvi_asv_sites.boot <- boot::boot(dist_usvi_asv.df, F_ttest_gn, variable = "site", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

#for reference: 
q_value <- 0.10
dist_usvi_asv_sites.boot.summary.df <- unique(dist_usvi_asv.df[["site"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_sites.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_asv_sites.boot$t %>%
                                 as.data.frame() %>%
                                 dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                 t())))

# dist_usvi_asv_sites.boot$t0[1] #0.03060117
# dist_usvi_asv_sites.boot$t0[2] #3.239403e-05
# dist_usvi_asv_sites.boot$t0[3] #0.4459671
# 
# 
# quantile(dist_usvi_asv_sites.boot$t[,1], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.0001492198
# quantile(dist_usvi_asv_sites.boot$t[,2], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #3.01793e-06
# quantile(dist_usvi_asv_sites.boot$t[,3], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.01629039


#older method, not generalized to any variable:
{
  # F_ttest <- function(d, f){
  #   #d is a dataframe containing all the measurements by site
  #   
  #   # length_v <- unique(d$site)
  #   length_v <- c("LB_seagrass", "Yawzi", "Tektite")
  #   
  #   temp_i_a <- d %>%
  #     # dplyr::filter(grepl(length_v[1], site)) %>%
  #     dplyr::filter(grepl("LB", site)) %>%
  #     dplyr::select(dissimilarity) %>%
  #     # tibble::deframe(.)
  #     tibble::deframe(.) %>% sample(., 25, replace = FALSE)
  #   temp_i_b <- d %>%
  #     # dplyr::filter(grepl(length_v[2], site)) %>%
  #     dplyr::filter(grepl("Yawzi", site)) %>%
  #     dplyr::select(dissimilarity) %>%
  #     # tibble::deframe(.)
  #     tibble::deframe(.) %>% sample(., 30, replace = FALSE)
  #   temp_i_c <- d %>%
  #     # dplyr::filter(grepl(length_v[3], site)) %>%
  #     dplyr::filter(grepl("Tektite", site)) %>%
  #     dplyr::select(dissimilarity) %>%
  #     # tibble::deframe(.)
  #     tibble::deframe(.) %>% sample(., 30, replace = FALSE)
  #   
  #   #list the resulting p-values from t-testing LB vs Yawzi, LB vs Tektite, and Yawzi vs Tektite
  #   c(t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value, 
  #     t.test(temp_i_a, temp_i_c, conf.level = 0.95)$p.value, 
  #     t.test(temp_i_b, temp_i_c, conf.level = 0.95)$p.value)
  # }
  # temp_boot_asv <- boot::boot(dist_usvi_asv.df, F_ttest, sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
  # 
  # q_value <- 0.10
  # quantile(temp_boot_asv$t[,1], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.0001014252
  # quantile(temp_boot_asv$t[,2], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.001321693
  # quantile(temp_boot_asv$t[,3], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.02430411
  # # temp_boot_asv$t0[1] <= quantile(temp_boot_asv$t[,1], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #LB and Yawzi are not sig diff
  # # temp_boot_asv$t0[2] <= quantile(temp_boot_asv$t[,2], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #LB and Tektite are not sig diff
  # # temp_boot_asv$t0[3] <= quantile(temp_boot_asv$t[,3], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #so Tektite and Yawzi are not sig diff in the metablomes
  # 
  # temp_boot_asv$t0[1]
  # temp_boot_asv$t0[2]
  # temp_boot_asv$t0[3]
}





# Test the metabolomes in PERMANOVA ---------------------------------------


#reepat for the metabolomes:

#SS(total) for metabolomes:  2.739826
(1/nrow(dist_usvi_metab.mat))*(sum((dist_usvi_metab.mat[lower.tri(dist_usvi_metab.mat, diag = FALSE)])^2))

#per-site Sum of Squares
meta.metab %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x]) %>%
  map(., ~sum((.x[lower.tri(.x, diag = FALSE)])^2)/nrow(.x))
#LB: 0.2827906
#Tektite: 0.9711734
#Yawzi: 0.7658972

with(meta.metab, vegan::adonis2(data = meta.metab, method = "bray", permutations = 9999,
                                # strata = site,
                                formula = dist_usvi_metab.d ~ sampling_time*site*sampling_day,
                                parallel = nthreads, by = "terms"))
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# vegan::adonis2(formula = dist_usvi_metab.d ~ sampling_time * site * sampling_day, data = meta.metab, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1  0.10903 0.03979  4.4923 0.0095 ** 
#   site                             2  0.72391 0.26422 14.9139 0.0001 ***
#   sampling_day                     3  0.09925 0.03622  1.3631 0.2067    
# sampling_time:site               2  0.08861 0.03234  1.8256 0.1068    
# sampling_time:sampling_day       3  0.09200 0.03358  1.2636 0.2636    
# site:sampling_day                6  0.31689 0.11566  2.1762 0.0134 *  
#   sampling_time:site:sampling_day  6  0.12092 0.04413  0.8304 0.6394    
# Residual                        49  1.18922 0.43405                   
# Total                           72  2.73983 1.00000    


#does it matter if we write the formula as sampling_time*sampling_day*site instead?

with(meta.metab, vegan::adonis2(data = meta.metab, method = "bray", permutations = 9999,
                                # strata = site,
                                formula = dist_usvi_metab.d ~ sampling_time*sampling_day*site,
                                parallel = nthreads, by = "terms"))
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1  0.10903 0.03979  4.4923 0.0091 ** 
#   sampling_day                     3  0.10704 0.03907  1.4702 0.1650    
# site                             2  0.71612 0.26137 14.7533 0.0001 ***
#   sampling_time:sampling_day       3  0.09027 0.03295  1.2399 0.2706    
# sampling_time:site               2  0.09034 0.03297  1.8611 0.0936 .  
# sampling_day:site                6  0.31689 0.11566  2.1762 0.0104 *  
#   sampling_time:sampling_day:site  6  0.12092 0.04413  0.8304 0.6446    
# Residual                        49  1.18922 0.43405                   
# Total                           72  2.73983 1.00000    


temp_res <- with(meta.metab, vegan::adonis2(data = meta.metab, method = "bray", permutations = 9999,
                                            # strata = site,
                                            formula = dist_usvi_metab.d ~ sampling_time*sampling_day*site,
                                            parallel = nthreads, by = "terms"))
permanova_metab_estimate_variation.df <- temp_res %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::mutate(term = as.character(term)) %>%
  dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                      .default = NA))
  
permanova_metab_estimate_variation.df <- permanova_metab_estimate_variation.df %>%
  dplyr::mutate(permanova_metab_estimate_variation.df %>%
                  dplyr::filter(!grepl("Total", term)) %>%
                  dplyr::summarise(TotalEMS = sum(MS))) %>%
  dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
                                          .default = NA)) %>%
  dplyr::mutate(sq_root = MS^(1/2))


permanova_metab_estimate_variation.df %>% 
  dplyr::summarise(sum(sq_root, na.rm = TRUE))

# term Df SumOfSqs      R2       F Pr(>F)      MS TotalEMS V_term sq_root
# [1,]    4  1  0.10903 0.03979  4.4923 0.0087 0.10903  0.67526 16.146 0.33019
# [2,]    2  3  0.10704 0.03907  1.4702 0.1655 0.03568  0.67526  5.284 0.18889
# [3,]    8  2  0.71612 0.26137 14.7533 0.0001 0.35806  0.67526 53.025 0.59838
# [4,]    5  3  0.09027 0.03295  1.2399 0.2700 0.03009  0.67526  4.456 0.17347
# [5,]    7  2  0.09034 0.03297  1.8611 0.0980 0.04517  0.67526  6.689 0.21253
# [6,]    3  6  0.31689 0.11566  2.1762 0.0127 0.05282  0.67526  7.821 0.22982
# [7,]    6  6  0.12092 0.04413  0.8304 0.6479 0.02015  0.67526  2.985 0.14196
# [8,]    1 49  1.18922 0.43405                0.02427  0.67526  3.594 0.15579
# [9,]    9 72  2.73983 1.00000                         0.67526      


#what is the interaction of terms within each factor?
dist_usvi_metab.df <- meta.metab %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "site")
dist_usvi_metab.df %>%
  dplyr::group_by(site) %>%
  # dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
  #                  min_dist = min(dissimilarity, na.rm = TRUE),
  #                  max_dist = max(dissimilarity, na.rm = TRUE))
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
                   min_dist = min(dissimilarity, na.rm = TRUE),
                   max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_metab.df[["dissimilarity"]]))))
# site        avg_dist min_dist max_dist avg_dist_rescaled min_dist_rescaled max_dist_rescaled
# <chr>          <dbl>    <dbl>    <dbl>             <dbl>             <dbl>             <dbl>
#   1 LB_seagrass    0.151   0.0568    0.334             0.236            0.0887             0.522
# 2 Tektite        0.257   0.0270    0.640             0.401            0.0421             1    
# 3 Yawzi          0.219   0.0357    0.602             0.342            0.0558             0.940

#Tektite vs Yawzi:
dist_usvi_metab.df %>%
  dplyr::filter(!grepl("LB", site)) %>%
  kruskal.test(dissimilarity ~ site, .)
# Kruskal-Wallis chi-squared = 15.546, df = 1, p-value = 8.052e-05

#Yawzi vs LB:
dist_usvi_metab.df %>%
  dplyr::filter(!grepl("Tektite", site)) %>%
  kruskal.test(dissimilarity ~ site, .)
# Kruskal-Wallis chi-squared = 28.413, df = 1, p-value = 9.801e-08

#Tektite vs LB:
dist_usvi_metab.df %>%
  dplyr::filter(!grepl("Yawzi", site)) %>%
  kruskal.test(dissimilarity ~ site, .)
# Kruskal-Wallis chi-squared = 119.69, df = 1, p-value < 2.2e-16

#all 3:
dist_usvi_metab.df %>%
  kruskal.test(dissimilarity ~ site, .)
# data:  dissimilarity by site
# Kruskal-Wallis chi-squared = 106.24, df = 2, p-value < 2.2e-16


pairwise.t.test((dist_usvi_metab.df %>%
                   split(., f = .$site) %>%
                   map(., ~.x %>%
                         dplyr::select(dissimilarity) %>%
                         tibble::deframe(.)) %>%
                   purrr::reduce(c)),
                (dist_usvi_metab.df %>%
                   split(., f = .$site) %>%
                   map(., ~.x %>%
                         dplyr::select(site) %>%
                         tibble::deframe(.)) %>%
                   purrr::reduce(c)),
                pool.sd = FALSE, paired = FALSE, alternative = "two.sided",
                p.adjust.method = "BH")

# LB_seagrass Tektite
# Tektite < 2e-16     -      
#   Yawzi   5.7e-16     0.00022
# 
# P value adjustment method: BH 



#old method:
{
  #make a random resampling
  # 
  # temp_responses <- (dist_usvi_metab.df %>%
  #     split(., f = .$site) %>%
  #     map(., ~.x %>%
  #           dplyr::select(dissimilarity) %>%
  #           tibble::deframe(.)) %>%
  #     purrr::reduce(c)) %>%
  #   sort(.)
  # #LB: sample(temp_responses, 253)
  # #Tektite and Yawzi: sample(temp_responses, 300)
  # sample(temp_responses, 253)
  # hist(temp_responses)
  # quantile(temp_responses, probs = seq(0, 1, 0.1), names = FALSE)
  # # [1] 0.02697357 0.08573037 0.10930988 0.13063755 0.15480734 0.18241717 0.22346676 0.25924514 0.30878149 0.37856812
  # # [11] 0.63994859
  # LaplacesDemon::p.interval(temp_responses, prob = 0.95)
  # # Lower     Upper
  # # [1,] 0.05489936 0.4568425
  # # attr(,"Probability.Interval")
  # # [1] 0.9495897
  # temp_boot <- boot::boot(dist_usvi_metab.df, F_ttest, sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
  
  
  # #bootstrap p-values approximated by:
  # sum(abs(temp_boot$t[,1]-1) > abs(temp_boot$t0[1]-1))/(1+temp_boot$R) #LB vs Yawzi
  # sum(abs(temp_boot$t[,2]-1) > abs(temp_boot$t0[2]-1))/(1+temp_boot$R) #LB vs Tektite
  # sum(abs(temp_boot$t[,3]-1) > abs(temp_boot$t0[3]-1))/(1+temp_boot$R) #Yawzi vs Tektite
  
  #look at the distribution of p-values obtained by bootstrapping, then slice at the 10% quantile
  q_value <- 0.10
  quantile(temp_boot$t[,1], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.0004332786
  quantile(temp_boot$t[,2], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #2.12744e-06
  quantile(temp_boot$t[,3], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #0.01809638
  temp_boot$t0[1] <= quantile(temp_boot$t[,1], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #LB and Yawzi are not sig diff
  temp_boot$t0[2] <= quantile(temp_boot$t[,2], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #LB and Tektite are not sig diff
  temp_boot$t0[3] <= quantile(temp_boot$t[,3], probs = q_value, na.rm = TRUE, names = FALSE,type = 7) #so Tektite and Yawzi are not sig diff in the metablomes
  
  
  
  #old way of simulating p-values from t.testing random samples of the responses
  {
    # modeled_ttest_p <- matrix(nrow = 9999, ncol = 1)
    # for(i in seq_len(9999)){
    #   temp_i_a <- sample(temp_responses, 253)
    #   temp_i_b <- sample(temp_responses, 300)
    #   temp_res <- t.test(temp_i_a, temp_i_b,
    #                      conf.level = 0.95)
    #   modeled_ttest_p[i, 1] <- temp_res$p.value
    # }
    # hist(modeled_ttest_p)
    # # padj_cutoff <- modeled_ttest_p %>%
    # #   as.vector(.) %>%
    # #   p.adjust(., method = "BH") %>% #multiple testing corrections, "BH" is an alias for "fdr" accordiong to p.adjust()
    # #   unlist
    # 
    # # padj_cutoff <- ashr::qval.from.lfdr((dist_usvi_metab.df %>%
    # # split(., f = .$site) %>%
    # # map(., ~.x %>%
    # # dplyr::select(dissimilarity) %>%
    # # tibble::deframe(.)) %>%
    # # purrr::reduce(c))) %>%
    # # unlist %>%
    # #   quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)
    # q_value <- 0.01
    # padj_cutoff <- modeled_ttest_p %>%
    #   as.vector(.) %>%
    #   # p.adjust(., method = "BH") %>% #multiple testing corrections, "BH" is an alias for "fdr" accordiong to p.adjust()
    #   # na.omit(.) %>%
    #   ashr::qval.from.lfdr(.) %>%
    #   unlist %>%
    #   quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)
    # 
  }
}


dist_usvi_metab_sites.boot <- boot::boot(dist_usvi_metab.df, F_ttest_gn, variable = "site", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

#for reference: 
q_value <- 0.10
dist_usvi_metab_sites.boot.summary.df <- unique(dist_usvi_metab.df[["site"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_sites.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_metab_sites.boot$t %>%
                                         setNames(., c("V1", "V2", "V3")) %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))





temp_res <- t.test(dist_usvi_metab.df %>%
         dplyr::filter(grepl("Tektite", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab.df %>%
         dplyr::filter(grepl("Yawzi", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       conf.level = 0.95)
# t = 3.7198, df = 597.53, p-value = 0.0002182
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01782635 0.05770426
# sample estimates:
#   mean of x mean of y 
# 0.2568152 0.2190499 

t.test(dist_usvi_metab.df %>%
         dplyr::filter(grepl("LB", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab.df %>%
         dplyr::filter(grepl("Tektite", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.))
# t = -13.479, df = 427.79, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.12144943 -0.09053627
# sample estimates:
#   mean of x mean of y 
# 0.1508223 0.2568152 

t.test(dist_usvi_metab.df %>%
         dplyr::filter(grepl("LB", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab.df %>%
         dplyr::filter(grepl("Yawzi", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.))
# t = -8.4801, df = 421.71, p-value = 3.832e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.08404212 -0.05241297
# sample estimates:
#   mean of x mean of y 
# 0.1508223 0.2190499 


#Summarizing intra-site distances via pariwise t-tests:


length_v <- c("LB_seagrass", "Yawzi", "Tektite")
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(p.value = NA, t = NA)
ttest_usvi_metab_sites.res <- ttest_res
ttest_usvi_asv_sites.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_metab_sites.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic
  ttest_usvi_metab_sites.res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
  
  temp_i_c <- dist_usvi_asv.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_asv_sites.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  ttest_usvi_asv_sites.res[i, "p.value"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$p.value
}


# Look at within-site within-sampling time --------------------------------


#heads up: the bootstrapped pairwise t-tests for 'dawn' vs 'peak_photo' microbial samples yielded:
# pair1      pair2   p.value       t  q.value
# 1  dawn peak_photo 0.5313026 1.52252  0.07870059

#for metabolics samples:
# pair1      pair2    p.value        t  q.value
# 1  dawn peak_photo 0.002738135 -3.12219 0.02773236

meta.microb %>%
  split(., f = .$sampling_time) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "sampling_time") %>% 
  dplyr::group_by(sampling_time) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# sampling_time avg_dist
# <chr>            <dbl>
#   1 dawn             0.145
# 2 peak_photo       0.138

meta.metab %>%
  split(., f = .$sampling_time) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "sampling_time") %>% 
  dplyr::group_by(sampling_time) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# sampling_time avg_dist
# <chr>            <dbl>
#   1 dawn             0.240
# 2 peak_photo       0.259


#bootstrap and save the results for 'dawn' vs 'aftenroon":
dist_usvi_metab_time.df <- meta.metab %>%
  split(., f = .$sampling_time) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "sampling_time")

dist_usvi_metab_time.boot <- boot::boot(dist_usvi_metab_time.df, F_ttest_gn, variable = "sampling_time", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_metab_time.boot.summary.df <- unique(dist_usvi_metab_time.df[["sampling_time"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_time.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_metab_time.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))
dist_usvi_asv_time.df <- meta.microb %>%
  split(., f = .$sampling_time) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "sampling_time")

dist_usvi_asv_time.boot <- boot::boot(dist_usvi_asv_time.df, F_ttest_gn, variable = "sampling_time", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_asv_time.boot.summary.df <- unique(dist_usvi_asv_time.df[["sampling_time"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_time.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_asv_time.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))

length_v <- unique(dist_usvi_metab_time.df[["sampling_time"]]) %>% as.character(.)
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA)
ttest_usvi_metab_time.res <- ttest_res
ttest_usvi_asv_time.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab_time.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab_time.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_metab_time.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic
  
  
  temp_i_c <- dist_usvi_asv_time.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv_time.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_asv_time.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  
}

ttest_usvi_metab_time.res <- ttest_usvi_metab_time.res %>%
  dplyr::left_join(., dist_usvi_metab_time.boot.summary.df)
ttest_usvi_asv_time.res <- ttest_usvi_asv_time.res %>%
  dplyr::left_join(., dist_usvi_asv_time.boot.summary.df)

ttest_usvi_metab_time.res

ttest_usvi_asv_time.res

# Site by time ------------------------------------------------------------


#so look at site x time:
dist_usvi_asv_grouping.df <- meta.microb %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "grouping")

dist_usvi_asv_grouping.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
                   min_dist = min(dissimilarity, na.rm = TRUE),
                   max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_asv_grouping.df[["dissimilarity"]]))))
# # A tibble: 6 × 2
# grouping               avg_dist
# <chr>                     <dbl>
#   1 LB_seagrass.dawn         0.120 
# 2 LB_seagrass.peak_photo   0.133 
# 3 Tektite.dawn             0.0837
# 4 Tektite.peak_photo       0.0859
# 5 Yawzi.dawn               0.0737
# 6 Yawzi.peak_photo         0.0487



dist_usvi_asv_grouping.boot <- boot::boot(dist_usvi_asv_grouping.df, F_ttest_gn, variable = "grouping", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_asv_grouping.boot.summary.df <- unique(dist_usvi_asv_grouping.df[["grouping"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_grouping.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_asv_grouping.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))


dist_usvi_metab_grouping.df <- meta.metab %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "grouping")


dist_usvi_metab_grouping.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
                   min_dist = min(dissimilarity, na.rm = TRUE),
                   max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_metab_grouping.df[["dissimilarity"]]))))
# grouping               avg_dist
# <chr>                     <dbl>
#   1 LB_seagrass.dawn          0.135
# 2 LB_seagrass.peak_photo    0.137
# 3 Tektite.dawn              0.269
# 4 Tektite.peak_photo        0.230
# 5 Yawzi.dawn                0.139
# 6 Yawzi.peak_photo          0.279


dist_usvi_metab_grouping.boot <- boot::boot(dist_usvi_metab_grouping.df, F_ttest_gn, variable = "grouping", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_metab_grouping.boot.summary.df <- unique(dist_usvi_metab_grouping.df[["grouping"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_grouping.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_metab_grouping.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))

length_v <- unique(dist_usvi_metab_grouping.df[["grouping"]]) %>% as.character(.)
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA)
ttest_usvi_metab_grouping.res <- ttest_res
ttest_usvi_asv_grouping.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab_grouping.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab_grouping.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_metab_grouping.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic

  
  temp_i_c <- dist_usvi_asv_time.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv_time.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_asv_grouping.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  
}

ttest_usvi_metab_grouping.res <- ttest_usvi_metab_grouping.res %>%
  dplyr::left_join(., dist_usvi_metab_grouping.boot.summary.df)
ttest_usvi_asv_grouping.res <- ttest_usvi_asv_grouping.res %>%
  dplyr::left_join(., dist_usvi_asv_grouping.boot.summary.df)

ttest_usvi_metab_grouping.res
# pair1                  pair2          t      p.value     q.value
# 1        LB_seagrass.dawn           Tektite.dawn -6.9524786 0.0506089059 0.007521892
# 2        LB_seagrass.dawn             Yawzi.dawn -0.3562355 0.5036596644 0.112459468
# 3        LB_seagrass.dawn LB_seagrass.peak_photo -0.2535510 0.2574257136 0.122252021
# 4        LB_seagrass.dawn     Tektite.peak_photo -6.5144061 0.0008105861 0.008802235
# 5        LB_seagrass.dawn       Yawzi.peak_photo -8.0718672 0.0780585121 0.003003581
# 6            Tektite.dawn             Yawzi.dawn  6.2649761 0.4928625541 0.008544316
# 7            Tektite.dawn LB_seagrass.peak_photo  6.8032034 0.0177295113 0.008336287
# 8            Tektite.dawn     Tektite.peak_photo  1.7979570 0.2630125324 0.091231116
# 9            Tektite.dawn       Yawzi.peak_photo -0.4077902 0.3254864876 0.122948678
# 10             Yawzi.dawn LB_seagrass.peak_photo  0.1504834 0.2874147171 0.100162518
# 11             Yawzi.dawn     Tektite.peak_photo -5.5162225 0.4307662395 0.010928800
# 12             Yawzi.dawn       Yawzi.peak_photo -7.2114314 0.1001830464 0.004208829
# 13 LB_seagrass.peak_photo     Tektite.peak_photo -6.3111660 0.0160331044 0.009203213
# 14 LB_seagrass.peak_photo       Yawzi.peak_photo -7.9059676 0.0013519500 0.003524839
# 15     Tektite.peak_photo       Yawzi.peak_photo -2.3857111 0.5707757935 0.060426384
ttest_usvi_asv_grouping.res
# pair1                  pair2          t     p.value     q.value
# 1        LB_seagrass.dawn           Tektite.dawn  3.2816177 0.033883002 0.037833663
# 2        LB_seagrass.dawn             Yawzi.dawn  4.1371680 0.309493280 0.019404630
# 3        LB_seagrass.dawn LB_seagrass.peak_photo -1.0499835 0.148447780 0.101599542
# 4        LB_seagrass.dawn     Tektite.peak_photo  2.9333845 0.505118351 0.036089589
# 5        LB_seagrass.dawn       Yawzi.peak_photo  7.2450485 0.207276782 0.003949202
# 6            Tektite.dawn             Yawzi.dawn  1.0605164 0.305865644 0.101223199
# 7            Tektite.dawn LB_seagrass.peak_photo -4.4325617 0.091043498 0.014618039
# 8            Tektite.dawn     Tektite.peak_photo -0.2215124 0.313041380 0.115449047
# 9            Tektite.dawn       Yawzi.peak_photo  4.4715409 0.090764553 0.022191146
# 10             Yawzi.dawn LB_seagrass.peak_photo -5.2640002 0.002083231 0.007848393
# 11             Yawzi.dawn     Tektite.peak_photo -1.2110278 0.115193619 0.101249962
# 12             Yawzi.dawn       Yawzi.peak_photo  3.1137935 0.049286459 0.053565491
# 13 LB_seagrass.peak_photo     Tektite.peak_photo  4.0355111 0.554597544 0.014515163
# 14 LB_seagrass.peak_photo       Yawzi.peak_photo  8.4670159 0.003359854 0.001419172
# 15     Tektite.peak_photo       Yawzi.peak_photo  4.3318484 0.383426594 0.029717709

#in the microbiomes, significant differences between:
#Yawzi dawn vs afternoon
#Lameshur dawn vs Tektite dawn
#Lameshur afternoon vs Yawzi afternoon
#Yawzi dawn vs Lameshur afternoon
#weakly, Tekeite dawn vs Lameshur afternoon
#weakly, Tektite dawn vs Yawzi afternoon


#in metabolomes:
#Lameshur afternoon vs Tekeite afternoon
#Lameshur afternoonv s Yawzi afternoon
#Lameshur dawn vs Tekeite afternoon
#Tekeite dawn vs Lameshur afternoon
#somewhat, Lameshur dawn vs Tektite dawn
#somewhat, Lameshur dawn vs Yawzi afternoon


# Pairwise t-tests for day and site ---------------------------------------

#note that day was significant for microbiomes, but not metabolomes:

dist_usvi_day.df <- meta.microb %>%
  split(., f = .$sampling_day) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "sampling_day")
dist_usvi_day.df %>% 
  dplyr::group_by(sampling_day) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# sampling_day avg_dist
# <chr>           <dbl>
#   1 Day1            0.107
# 2 Day2            0.176
# 3 Day3            0.140
# 4 Day4            0.128
# 5 Day5            0.135

dist_usvi_day.df <- dist_usvi_day.df %>%
  dplyr::mutate(type = "asv") %>%
  bind_rows(., (meta.metab %>%
  split(., f = .$sampling_day) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "sampling_day") %>%
    dplyr::mutate(type = "metab")))


dist_usvi_day.df %>%
  dplyr::group_by(type, sampling_day) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))  %>%
  dplyr::left_join(., dist_usvi_day.df %>%
                     dplyr::group_by(type) %>%
                     dplyr::summarise(max_dist = max(dissimilarity, na.rm = TRUE))) %>%
  dplyr::mutate(avg_dist_rescale = avg_dist/max_dist)
  # dplyr::mutate(avg_dist_rescaled = avg_dist/max_dist)
# sampling_day avg_dist
# <chr>           <dbl>
#   1 Day2            0.259
# 2 Day3            0.255
# 3 Day4            0.250
# 4 Day5            0.253


#repeat for day and site:

dist_usvi_asv_site_day.df <- meta.microb %>%
  dplyr::mutate(site_day = paste0(site, ".", sampling_day)) %>%
  split(., f = .$site_day) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "site_day")

dist_usvi_asv_site_day.df %>%
  dplyr::group_by(site_day) %>%
  # dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
  #                  min_dist = min(dissimilarity, na.rm = TRUE),
  #                  max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_asv_site_day.df[["dissimilarity"]]))))
# site_day         avg_dist avg_dist_rescaled
# <chr>               <dbl>             <dbl>
#   1 LB_seagrass.Day1   0.147             0.373 
# 2 LB_seagrass.Day2   0.0884            0.224 
# 3 LB_seagrass.Day3   0.110             0.279 
# 4 LB_seagrass.Day4   0.0924            0.235 
# 5 LB_seagrass.Day5   0.119             0.303 
# 6 Tektite.Day1       0.0703            0.178 
# 7 Tektite.Day2       0.102             0.260 
# 8 Tektite.Day3       0.0459            0.117 
# 9 Tektite.Day4       0.0601            0.153 
# 10 Tektite.Day5       0.0486            0.123 
# 11 Yawzi.Day1         0.0270            0.0685
# 12 Yawzi.Day2         0.0874            0.222 
# 13 Yawzi.Day3         0.0527            0.134 
# 14 Yawzi.Day4         0.0657            0.167 
# 15 Yawzi.Day5         0.0508            0.129 



dist_usvi_asv_site_day.boot <- boot::boot(dist_usvi_asv_site_day.df, F_ttest_gn, variable = "site_day", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_asv_site_day.boot.summary.df <- unique(dist_usvi_asv_site_day.df[["site_day"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_site_day.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_asv_site_day.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))


dist_usvi_metab_site_day.df <- meta.metab %>%
  dplyr::mutate(site_day = paste0(site, ".", sampling_day)) %>%
  split(., f = .$site_day) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "metab_deriv_label") %>%
        dplyr::select(metab_deriv_label) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_metab.mat[.x, .x])  %>%
  map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
  map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
  bind_rows(., .id = "site_day")


dist_usvi_metab_site_day.df %>%
  dplyr::group_by(site_day) %>%
  # dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
  #                  min_dist = min(dissimilarity, na.rm = TRUE),
  #                  max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_metab_site_day.df[["dissimilarity"]]))))
# site_day         avg_dist avg_dist_rescaled
# <chr>               <dbl>             <dbl>
#   1 LB_seagrass.Day2    0.112             0.197
# 2 LB_seagrass.Day3    0.103             0.180
# 3 LB_seagrass.Day4    0.112             0.195
# 4 LB_seagrass.Day5    0.212             0.370
# 5 Tektite.Day2        0.321             0.561
# 6 Tektite.Day3        0.283             0.495
# 7 Tektite.Day4        0.205             0.359
# 8 Tektite.Day5        0.234             0.409
# 9 Yawzi.Day2          0.209             0.365
# 10 Yawzi.Day3          0.160             0.280
# 11 Yawzi.Day4          0.179             0.313
# 12 Yawzi.Day5          0.213             0.372


dist_usvi_metab_site_day.boot <- boot::boot(dist_usvi_metab_site_day.df, F_ttest_gn, variable = "site_day", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_metab_site_day.boot.summary.df <- unique(dist_usvi_metab_site_day.df[["site_day"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_site_day.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q.value" = (dist_usvi_metab_site_day.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))

length_v <- unique(dist_usvi_asv_site_day.df[["site_day"]]) %>% as.character(.)
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA)
ttest_usvi_metab_site_day.res <- ttest_res
ttest_usvi_asv_site_day.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab_site_day.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab_site_day.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  if(length(temp_i_a) > 0 & length(temp_i_b) > 0){
    ttest_usvi_metab_site_day.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic  
  } else {
    ttest_usvi_metab_site_day.res[i, "t"] <- NA
  }

  temp_i_c <- dist_usvi_asv_site_day.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv_site_day.df %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  if(length(temp_i_c) > 0 & length(temp_i_d) > 0){
    ttest_usvi_asv_site_day.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  } else {
    ttest_usvi_asv_site_day.res[i, "t"] <- NA
  }
}

ttest_usvi_metab_site_day.res <- ttest_usvi_metab_site_day.res %>%
  tidyr::drop_na(t) %>%
  dplyr::left_join(., dist_usvi_metab_site_day.boot.summary.df)
ttest_usvi_asv_site_day.res <- ttest_usvi_asv_site_day.res %>%
  tidyr::drop_na(t) %>%
  dplyr::left_join(., dist_usvi_asv_site_day.boot.summary.df)


ttest_usvi_asv_site_day.res %>%
  dplyr::filter(p.value < 0.1) %>%
  # dplyr::summarise(num_SDA = length(t))
  droplevels
#105 results, of which 6 have bootstrapped p-values < 0.1
#none are between the same days across sites, or between days at the same site

ttest_usvi_metab_site_day.res %>%
  dplyr::filter(p.value < 0.1) %>%
  # dplyr::summarise(num_SDA = length(t))
  droplevels
#66 results, of which 7 have bootstrapped p-values < 0.1
#2 comparisons are between days at LB
#2 comparisons are between Yawzi and LB on days 3 and 4
#1 comparison is between Yawzi and Tektite on day 4


# adonis2 per site --------------------------------------------------------


#let's subset for individual sites, to constrain the site differences

site.meta.metab <- metabolomics_sample_metadata %>%
  # dplyr::select(intersect(colnames(metabolomics_sample_metadata), keep)) %>%
  dplyr::left_join(., usvi_selected_metadata %>%
                     dplyr::select(sample_id, site, site_type, grouping, replicate, PAR, lumens, lux, temp), multiple = "all", relationship = "many-to-many") %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  tidyr::fill(PAR, lumens, lux, temp, .direction = "down") %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(metab_deriv_label %in% colnames(usvi_metabolomics.tbl)) %>%
  dplyr::mutate(metab_deriv_label = factor(metab_deriv_label, levels = colnames(usvi_metabolomics.tbl))) %>%
  dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
  dplyr::mutate(grouping2 = paste0(grouping, ".", sampling_day)) %>%
  dplyr::select(sample_id, colnames(meta.microb)) %>%
  dplyr::mutate(metab_deriv_label = factor(metab_deriv_label, levels = colnames(usvi_metabolomics.tbl))) %>%
  droplevels %>%
  dplyr::arrange(metab_deriv_label) %>%
  # tidyr::drop_na(.) %>%
  dplyr::mutate(across(c(sample_id, sample_type, site, sampling_time, sampling_day, metab_deriv_label, site_type, grouping, grouping2), ~factor(.x))) %>%
  dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                grouping = fct_relevel(grouping, "LB_seagrass.dawn"),
                grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day2")) %>%
  dplyr::filter(grepl("Tektite", site)) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  droplevels
dist_usvi_site_metab.d <- usvi_metabolomics.tbl %>%
  dplyr::select(rownames(site.meta.metab)) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  t() %>%
  vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)

with(site.meta.metab, vegan::adonis2(data = site.meta.metab, method = "bray", permutations = 9999,
                                # strata = site,
                                formula = dist_usvi_site_metab.d ~ sampling_time*sampling_day,
                                parallel = nthreads, by = "terms"))
#for LB seagrass:
# vegan::adonis2(formula = dist_usvi_site_metab.d ~ sampling_time * sampling_day, data = site.meta.metab, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2      F Pr(>F)    
# sampling_time               1 0.056188 0.19869 7.2475 0.0002 ***
#   sampling_day                3 0.069641 0.24626 2.9942 0.0037 ** 
#   sampling_time:sampling_day  3 0.040670 0.14381 1.7486 0.0814 .  
# Residual                   15 0.116292 0.41123                  
# Total                      22 0.282791 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#for yawzi:
# vegan::adonis2(formula = dist_usvi_site_metab.d ~ sampling_time * sampling_day, data = site.meta.metab, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2      F Pr(>F)   
# sampling_time               1  0.07629 0.09961 3.6698 0.0445 * 
#   sampling_day                3  0.25303 0.33037 4.0571 0.0062 **
#   sampling_time:sampling_day  3  0.08317 0.10859 1.3335 0.2589   
# Residual                   17  0.35341 0.46144                 
# Total                      24  0.76590 1.00000                 
# ---

#for Tektite:
# vegan::adonis2(formula = dist_usvi_site_metab.d ~ sampling_time * sampling_day, data = site.meta.metab, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2      F Pr(>F)
# sampling_time               1  0.06751 0.06951 1.5951 0.1969
# sampling_day                3  0.09062 0.09331 0.7137 0.6296
# sampling_time:sampling_day  3  0.09353 0.09630 0.7366 0.6175
# Residual                   17  0.71951 0.74087              
# Total                      24  0.97117 1.00000      


#calculating estimates of variation in the models:
#start with Tektite, since it is so variable.

temp_res <- with(site.meta.metab, vegan::adonis2(data = site.meta.metab, method = "bray", permutations = 9999,
                                     # strata = site,
                                     formula = dist_usvi_site_metab.d ~ sampling_time*sampling_day,
                                     parallel = nthreads, by = "terms"))

permanova_tektite_metab_estimate_variation.df <- temp_res %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::mutate(term = as.character(term)) %>%
  dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                      .default = NA))

permanova_tektite_metab_estimate_variation.df <- permanova_tektite_metab_estimate_variation.df %>%
  dplyr::mutate(permanova_tektite_metab_estimate_variation.df %>%
                  dplyr::filter(!grepl("Total", term)) %>%
                  dplyr::summarise(TotalEMS = sum(MS))) %>%
  dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
                                          .default = NA)) %>%
  dplyr::mutate(sq_root = MS^(1/2))


permanova_tektite_metab_estimate_variation.df %>% 
  dplyr::summarise(sum(sq_root, na.rm = TRUE))
#term   V_term
#sampling_time    39.42973
#sampling_day   17.64227
#sampling_time:sampling_day   18.20849
#Residual   24.71950
#Total



# How to get type III -----------------------------------------------------



#calculating F-ratios for types II and II SS

# options(contrasts = c("contr.sum","contr.poly"))
# model <- lm(time ~ topic * sys, data=search)
# drop1(model, .~., test="F")

# #in lm, can use the distance matrix as input, but for car::Anova we need one response variable for each sample
## also, lm on the distance matrix generates a "mlm" object which can't be subjected to drop1()
# model <- lm(dist_usvi_metab.mat ~ sampling_time*sampling_day*site - 1, data = meta.metab,
#             contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))


temp_dist.df <- vegan::betadisper(dist_usvi_metab.d, type = "median",
                  sqrt.dist = TRUE,
                    meta.metab$grouping2) %>%
  purrr::pluck("distances") %>%  
  tibble::enframe(value = "dissimilarity", name = "metab_deriv_label")  %>%
  dplyr::left_join(., meta.metab %>%
                     dplyr::select(site, sampling_time, sampling_day) %>%
                     tibble::rownames_to_column(var = "metab_deriv_label"),
                   by = join_by(metab_deriv_label)) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  droplevels


model <- lm(dissimilarity ~ sampling_time*sampling_day*site, data = temp_dist.df,
                # contrasts=list(sampling_time="contr.helmert", sampling_day="contr.helmert", site="contr.helmert"))
                contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))
drop1(model, .~., test="F")
# Single term deletions
# 
# Model:
#   dissimilarity ~ sampling_time * sampling_day * site
# Df Sum of Sq     RSS     AIC F value   Pr(>F)   
# <none>                                       0.60252 -302.19                    
# sampling_time                    1  0.004907 0.60743 -303.60  0.3991 0.530510   
# sampling_day                     3  0.022186 0.62471 -305.55  0.6014 0.617173   
# site                             2  0.130498 0.73302 -291.88  5.3064 0.008203 **
#   sampling_time:sampling_day       3  0.027702 0.63022 -304.91  0.7510 0.527070   
# sampling_time:site               2  0.030140 0.63266 -302.62  1.2256 0.302428   
# sampling_day:site                6  0.036791 0.63931 -309.86  0.4987 0.806255   
# sampling_time:sampling_day:site  6  0.029381 0.63190 -310.71  0.3982 0.876589   


car::Anova(lm(dissimilarity ~ sampling_time*sampling_day*site, data = temp_dist.df,
              contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum")), 
           type=3, contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))

# Anova Table (Type III tests)
# 
# Response: dissimilarity
# Sum Sq Df  F value    Pr(>F)    
# (Intercept)                     3.7653  1 306.2155 < 2.2e-16 ***
#   sampling_time                   0.0049  1   0.3991  0.530510    
# sampling_day                    0.0222  3   0.6014  0.617173    
# site                            0.1305  2   5.3064  0.008203 ** 
#   sampling_time:sampling_day      0.0277  3   0.7510  0.527070    
# sampling_time:site              0.0301  2   1.2256  0.302428    
# sampling_day:site               0.0368  6   0.4987  0.806255    
# sampling_time:sampling_day:site 0.0294  6   0.3982  0.876589    
# Residuals                       0.6025 49                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# car::Anova(lm(dissimilarity ~ sampling_time*sampling_day*site - 1, data = temp_dist.df,
#               contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum")), 
#            type=3, contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))
# # Anova Table (Type III tests)
# # 
# # Response: dissimilarity
# # Sum Sq Df  F value    Pr(>F)    
# # sampling_time                   3.7662  2 153.1435 < 2.2e-16 ***
# #   sampling_day                    0.0222  3   0.6014  0.617173    
# # site                            0.1305  2   5.3064  0.008203 ** 
# #   sampling_time:sampling_day      0.0277  3   0.7510  0.527070    
# # sampling_time:site              0.0301  2   1.2256  0.302428    
# # sampling_day:site               0.0368  6   0.4987  0.806255    
# # sampling_time:sampling_day:site 0.0294  6   0.3982  0.876589    
# # Residuals                       0.6025 49                       
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 
# temp_mat <- usvi_metabolomics.tbl %>%
#   dplyr::select(rownames(meta.metab)) %>%
#   apply(., 2, function(x) log2(x + 1)) %>%
#   t()
# # pco(., dist = "bray", sqrt.dist = TRUE)
# # model <- vegan::dbrda(temp_mat ~ sampling_time*sampling_day*site, data = meta.metab, sqrt.dist = TRUE, dist = "bray")
# # rda(temp_mat ~ sampling_time*sampling_day*site, data = meta.metab)
# # # model <- vegan::dbrda(dist_usvi_metab.d ~ sampling_time*sampling_day*site, data = meta.metab, sqrt.dist = TRUE)
# # summary(model)                  
# # 
# 
# temp_res2 <- vegan::dbrda(temp_mat ~ sampling_time*sampling_day*site, data = meta.metab, sqrt.dist = TRUE, dist = "bray")
# drop1(temp_res2, scope=formula(temp_res2), test="perm")
# 
temp_res3 <- vegan::cca(temp_mat ~ sampling_time*sampling_day*site, data = meta.metab, dist = "bray", by = "terms")
# drop1(temp_res3, scope=formula(temp_res3), test="perm")
# 

vcov(temp_res3, complete=TRUE)
Anova(temp_res3, type=3, test.statistic = "F", vcov. = vcov(temp_res3, complete = TRUE), singular.ok = TRUE)
# Analysis of Deviance Table (Type III tests)
# 
# Response: temp_mat
# Df       F    Pr(>F)    
# sampling_time                    1  3.6461   0.06207 .  
# sampling_day                     3  0.3548   0.78583    
# site                             2 69.8436 4.509e-15 ***
#   sampling_time:sampling_day       3  1.1660   0.33226    
# sampling_time:site               2  1.3846   0.26005    
# sampling_day:site                6  1.0927   0.38007    
# sampling_time:sampling_day:site  6  1.9672   0.08858 .  
# Residuals                       49                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


temp_mat <- usvi_metabolomics.tbl %>%
  dplyr::select(rownames(meta.metab)) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  t()
# temp_mod <- lm(temp_mat ~ sampling_time*sampling_day*site - 1, data = meta.metab,
#                             contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))
# # #can't use the btw-sample distance matrix as input for car::Manova(lm(.)) or car::Anova(lm(.))
# # temp_mod <- lm(dist_usvi_metab.mat ~ sampling_time*sampling_day*site, data = meta.metab,
# #                contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))
# 
# temp_mod <- glm(dissimilarity ~ sampling_time*sampling_day*site, data = temp_dist.df)
glm(dissimilarity ~ sampling_time*sampling_day*site, data = temp_dist.df) %>%
  Anova(., type=3, contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"),
        test.statistic = "F", singular.ok = TRUE)
# Analysis of Deviance Table (Type III tests)
# 
# Response: dissimilarity
# Error estimate based on Pearson residuals 
# 
# Sum Sq Df F values   Pr(>F)   
# sampling_time                   0.00491  1   0.3991 0.530510   
# sampling_day                    0.02219  3   0.6014 0.617173   
# site                            0.13050  2   5.3064 0.008203 **
#   sampling_time:sampling_day      0.02770  3   0.7510 0.527070   
# sampling_time:site              0.03014  2   1.2256 0.302428   
# sampling_day:site               0.03679  6   0.4987 0.806255   
# sampling_time:sampling_day:site 0.02938  6   0.3982 0.876589   
# Residuals                       0.60252 49                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# temp_mat2 <- usvi_metabolomics.tbl %>%
#   dplyr::select(rownames(meta.metab)) %>%
#   apply(., 2, function(x) log2(x + 1)) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "metabolite") %>%
#   tidyr::pivot_longer(., cols = !c("metabolite"),
#                       names_to = "metab_deriv_label",
#                       values_to = "conc") %>%
#   split(., f = .$metabolite) %>%
#   map(., ~.x %>%
#         dplyr::select(metab_deriv_label, conc) %>%
#         droplevels %>%
#         tibble::deframe(.)) %>%
#   purrr::reduce(c)
# 
# temp_group <- usvi_metabolomics.tbl %>%
#   dplyr::select(rownames(meta.metab)) %>%
#   apply(., 2, function(x) log2(x + 1)) %>%
#   t() %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "metab_deriv_label") %>%
#   dplyr::left_join(., meta.metab %>%
#                      tibble::rownames_to_column(var = "metab_deriv_label")) %>%
#   dplyr::select(metab_deriv_label, site, sampling_time, sampling_day) %>%
#   tidyr::pivot_longer(., cols = !c("metab_deriv_label"),
#                       names_to = "group",
#                       values_to = "value") %>% 
#   split(., f = .$group) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         dplyr::distinct(metab_deriv_label, value) %>%
#         dplyr::mutate(value = as.character(value)) %>%
#         tibble::deframe(.) %>% simplify(.))
# # tibble::deframe(.)) %>% purrr::reduce(c)
# 
# # temp_mat2 <- bind_cols(temp_mat, meta.metab[, c("site", "sampling_time", "sampling_day")])

manova(temp_mat ~ sampling_time*sampling_day*site , data = bind_cols(temp_mat, meta.metab),
       contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum")) %>%
  Manova(., type=3, contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"),
         test.statistic = "Wilks", singular.ok = TRUE)
# Type III MANOVA Tests: Wilks test statistic
# Df  test stat approx F num Df  den Df  Pr(>F)   
# sampling_time                    2 5.7000e-09  270.507     98  2.0000 0.00369 **
#   sampling_day                     3 4.4900e-08    7.532    147  3.9293 0.03133 * 
#   site                             2 2.2131e-06   13.698     98  2.0000 0.07035 . 
# sampling_time:sampling_day       3 2.5835e-06    1.929    147  3.9293 0.28088   
# sampling_time:site               2 1.4059e-05    5.422     98  2.0000 0.16812   
# sampling_day:site                6 0.0000e+00    2.685    294 14.9605 0.01516 * 
#   sampling_time:sampling_day:site  6 2.0000e-10    2.201    294 14.9605 0.03943 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

temp_ss_res <- manova(temp_mat ~ sampling_time*sampling_day*site -1 , data = bind_cols(temp_mat, meta.metab),
       contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum")) %>%
  Anova(., type = 3, test.statistic = "Wilks", singular.ok = TRUE)
# temp_ss_res %>%
#   summary(., univariate=TRUE, multivariate=FALSE,
#         p.adjust.method=TRUE)
print(temp_ss_res)


lm(temp_mat ~ sampling_time*sampling_day*site -1 , data = bind_cols(temp_mat, meta.metab),
       contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum")) %>%
  Manova(., type = 3, test.statistic = "Wilks", singular.ok = TRUE)

#doesn'tw ork: 
# nlme::lme(temp_mat ~ sampling_time + sampling_day + site, data = bind_cols(temp_mat, meta.metab))




# temp_res2 <- manova(dist_usvi_metab.mat ~ sampling_time*sampling_day*site, data = meta.metab, contrasts=list(sampling_time="contr.poly", sampling_day="contr.poly", site="contr.poly"))
# car::Manova(temp_res2, type = 3)
# 
# car::Manova(manova(dist_usvi_metab.mat ~ sampling_time*sampling_day*site, data = meta.metab, contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum")),
#            type=3, contrasts=list(sampling_time="contr.sum", sampling_day="contr.sum", site="contr.sum"))
# 
# drop1(temp_res2, scope=formula(temp_res2), test="perm")
# 
# 
# with(meta.metab, vegan::adonis2(data = meta.metab, method = "bray", permutations = 9999,
#                                      # strata = site,
#                                 formula = temp_mat ~ sampling_time*sampling_day*site, 
#                                      parallel = nthreads, by = "terms"))
# # Permutation test for adonis under reduced model
# # Terms added sequentially (first to last)
# # Permutation: free
# # Number of permutations: 9999
# # 
# # vegan::adonis2(formula = temp_mat ~ sampling_time * sampling_day * site, data = meta.metab, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
# # Df SumOfSqs      R2       F Pr(>F)    
# # sampling_time                    1  0.10903 0.03979  4.4923 0.0081 ** 
# #   sampling_day                     3  0.10704 0.03907  1.4702 0.1660    
# # site                             2  0.71612 0.26137 14.7533 0.0001 ***
# #   sampling_time:sampling_day       3  0.09027 0.03295  1.2399 0.2642    
# # sampling_time:site               2  0.09034 0.03297  1.8611 0.1005    
# # sampling_day:site                6  0.31689 0.11566  2.1762 0.0098 ** 
# #   sampling_time:sampling_day:site  6  0.12092 0.04413  0.8304 0.6489    
# # Residual                        49  1.18922 0.43405                   
# # Total                           72  2.73983 1.00000                   
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 

(vegan::adonis2(data = temp_dist.df, method = "bray", permutations = 9999,
                # strata = site,
                formula = temp_dist.df$dissimilarity ~ sampling_time*sampling_day*site, 
                parallel = nthreads, by = "terms"))
# vegan::adonis2(formula = temp_dist.df$dissimilarity ~ sampling_time * sampling_day * site, data = temp_dist.df, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2      F Pr(>F)   
# sampling_time                    1  0.06039 0.02045 1.5547 0.2113   
# sampling_day                     3  0.13526 0.04580 1.1608 0.3189   
# site                             2  0.46074 0.15600 5.9309 0.0036 **
#   sampling_time:sampling_day       3  0.08235 0.02788 0.7067 0.5611   
# sampling_time:site               2  0.10370 0.03511 1.3348 0.2623   
# sampling_day:site                6  0.10166 0.03442 0.4362 0.8874   
# sampling_time:sampling_day:site  6  0.10598 0.03588 0.4547 0.8700   
# Residual                        49  1.90328 0.64445                 
# Total                           72  2.95337 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# Old stuff ---------------------------------------------------------------




#calculate within-sampling time by sampling day:


# 
# pairwise.t.test(c(temp_list[["LB_seagrass"]], temp_list[["Tektite"]], temp_list[["Yawzi"]]),
#                 c(rep("LB_seagrass", length(temp_list[["LB_seagrass"]])),
#                   rep("Tektite", length(temp_list[["Tektite"]])),
#                   rep("Yawzi", length(temp_list[["Yawzi"]]))),
#                 pool.sd = FALSE, paired = FALSE,
#                 p.adjust.method = "fdr")
# # LB_seagrass Tektite
# # Tektite < 2e-16     -
# #   Yawzi   < 2e-16     7.4e-05
# #
# # P value adjustment method: fdr




#per-site Sum of Squares
dist_usvi_asv.ss <- meta.microb %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x]) %>%
  map(., ~sum((.x[lower.tri(.x, diag = FALSE)])^2)/nrow(.x))
#LB: 0.4070378
#Tektite: 0.1798429
#Yawzi: 0.1220277


#per-day sum of squares
meta.microb %>%
  split(., f = .$sampling_day) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x]) %>%
  map(., ~sum((.x[lower.tri(.x, diag = FALSE)])^2)/nrow(.x))

#within-site SS for sampling time
meta.microb %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x]) %>%
  map(., ~sum((.x[lower.tri(.x, diag = FALSE)])^2)/nrow(.x))

# meta.microb2 <- meta.microb %>%
#   # dplyr::mutate(site = factor(site, levels = c("none", "LB_seagrass", "Tektite", "Yawzi")))
#   dplyr::mutate(site = fct_relevel(site, "LB_seagrass"))
# levels(meta.microb2$site)
# 
with(meta.microb, vegan::adonis2(data = meta.microb2, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site_type:site,
                                 # formula = dist_usvi_asv.d ~ site/sampling_time,
                                 parallel = nthreads, by = "onedf"))



# with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
#                                  strata = site,
#                                  formula = dist_usvi_asv.d ~ site + sampling_time + sampling_day + site/(sampling_time*sampling_day), #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
#                                  parallel = nthreads, by = "terms"))
# # Df SumOfSqs      R2       F Pr(>F)    
# # site                             2  0.78784 0.52637 93.3716  0.001 ***
# #   sampling_time                    1  0.10578 0.07067 25.0733  0.001 ***
# #   sampling_day                     4  0.20673 0.13812 12.2504  0.001 ***
# #   site:sampling_time               2  0.01104 0.00738  1.3089  0.296    
# # site:sampling_day                8  0.08363 0.05587  2.4778  0.043 *  
# #   site:sampling_time:sampling_day 12  0.05281 0.03529  1.0432  0.421    
# # Residual                        59  0.24891 0.16630                   
# # Total                           88  1.49675 1.00000  
# 
# #sampling time within site: 0.3232394
# (0.01104/2)/((1.49675-0.01104)/(89-2))



#adonis2 output of an F statistics:
# (35564/3)/((161648)/(64)) = 4.69
#what it should have been using PRIMER: 2.8086
(35564/3)/((16883)/(4))
#this is: [SS(A)/(a-1)]/[SS(B:A)/(N - (B:A))]


with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site/sampling_time/sampling_day,
                                 parallel = nthreads, by = "terms"))
# vegan::adonis2(formula = dist_usvi_asv.d ~ site/sampling_time/sampling_day, data = meta.microb, permutations = 999, method = "horn", by = "terms", parallel = nthreads, strata = site)
# Df SumOfSqs      R2       F Pr(>F)    
# site                             2  0.78784 0.52637 93.3716  0.001 ***
#   site:sampling_time               3  0.11606 0.07754  9.1698  0.001 ***
#   site:sampling_time:sampling_day 24  0.34394 0.22979  3.3969  0.003 ** 
#   Residual                        59  0.24891 0.16630                   
# Total                           88  1.49675 1.00000  

#F statistic hand-calculations:

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "margin"))
# vegan::adonis2(formula = dist_usvi_asv.d ~ site, data = meta.microb, permutations = 999, method = "horn", by = "margin", parallel = nthreads)
# Df SumOfSqs      R2      F Pr(>F)    
# site      2  0.78784 0.52637 47.788  0.001 ***
#   Residual 86  0.70891 0.47363                  
# Total    88  1.49675 1.00000    

#this is the F-ratio for just site: (0.78784/2)/(0.70891/86)

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site + sampling_time + sampling_day, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "margin"))
# Df SumOfSqs      R2      F Pr(>F)    
# site           2  0.78948 0.52746 80.661  0.001 ***
#   sampling_time  1  0.10558 0.07054 21.575  0.001 ***
#   sampling_day   4  0.20673 0.13812 10.561  0.001 ***
#   Residual      81  0.39640 0.26484                  
# Total         88  1.49675 1.00000      


with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site/sampling_time + site/sampling_day, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "margin"))
# Df SumOfSqs      R2      F Pr(>F)    
# site:sampling_time  3  0.11627 0.07768 9.1198  0.001 ***
#   site:sampling_day  12  0.29112 0.19450 5.7088  0.001 ***
#   Residual           71  0.30173 0.20159                  
# Total              88  1.49675 1.00000    

#for just sampling_time within site: 2.414431
(0.11627/3)/((1.49675-0.11627)/(89-3))


#for just sampling_day within site: 1.549414
(0.29112/12)/((1.49675-0.29112)/(89-12))

#calculate residual of SS before interaciton terms: 0.39640
#SS: 1.49675
#SS(site): 0.78948
#SS(sampling_time): 0.10558
#SS(sampling_day): 0.20673
#0.39640


#note that in type III SS, the sum of the SS of the interaction terms won't add up to this remainder.

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time:sampling_day, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "margin"))
# vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time:sampling_day, data = meta.microb, permutations = 999, method = "horn", by = "terms", parallel = nthreads, strata = site)
# Df SumOfSqs      R2      F Pr(>F)    
# sampling_time:sampling_day  9  0.34792 0.23245 2.6584  0.001 ***
#   Residual                   79  1.14883 0.76755                  
# Total                      88  1.49675 1.00000        


with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site/(sampling_time + sampling_day), #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "margin"))
# Df SumOfSqs      R2      F Pr(>F)    
# site:sampling_time  3  0.11627 0.07768 9.1198  0.001 ***
#   site:sampling_day  12  0.29112 0.19450 5.7088  0.001 ***
#   Residual           71  0.30173 0.20159                  
# Total              88  1.49675 1.00000            


#SS(sampling_time:sampling_day): 0.34792
#SS(site:sampling_time): 0.11627
#SS(site:sampling_day): 0.29112

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ site + site/(sampling_time*sampling_day), #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "terms"))

# Df SumOfSqs      R2       F Pr(>F)    
# site                             2  0.78784 0.52637 93.3716  0.001 ***
#   site:sampling_time               3  0.11606 0.07754  9.1698  0.001 ***
#   site:sampling_day               12  0.29112 0.19450  5.7505  0.001 ***
#   site:sampling_time:sampling_day 12  0.05281 0.03529  1.0432  0.429    
# Residual                        59  0.24891 0.16630                   
# Total                           88  1.49675 1.00000     



with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time + sampling_day/sampling_time, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "terms"))
# Df SumOfSqs      R2      F Pr(>F)    
# sampling_time               1  0.11015 0.07360 7.5748  0.001 ***
#   sampling_day                4  0.20072 0.13410 3.4507  0.001 ***
#   sampling_time:sampling_day  4  0.03705 0.02475 0.6369  0.282    
# Residual                   79  1.14883 0.76755                  
# Total                      88  1.49675 1.00000      


with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "terms"))

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 # strata = site,
                                 formula = dist_usvi_asv.d ~ site + sampling_time + sampling_day + site/(sampling_time:sampling_day), #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 parallel = nthreads, by = "terms"))


adonis_asv.res <- vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 # formula = dist_usvi_asv.d ~ site*sampling_time, #pr site:sampling_time = 0.226
                                 # formula = dist_usvi_asv.d ~ (site*sampling_time) + site:sampling_time:sampling_day, #pr = 0.072 for site:sampling_time
                                 # formula = dist_usvi_asv.d ~ (site*sampling_time)/sampling_day, #site:sampling_time pr = 0.085
                                 # formula = dist_usvi_asv.d ~ site + sampling_time + site:sampling_time + site:sampling_time:sampling_day, #site:sampling_time pr = 0.081
                                 # formula = dist_usvi_asv.d ~ site/sampling_day/sampling_time, #site:sampling_time pr = 0.002
                                 # formula = dist_usvi_asv.d ~ sampling_time*sampling_day, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                 formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + sampling_time:sampling_day + sampling_time:site + sampling_day:site + sampling_time:sampling_day:site, #pr sampling_day: 0.013; sampling_time: 0.005; sampling_day:sampling_time = 0.635
                                   parallel = nthreads, by = "terms")

# adonis_asv.res <- with(meta.microb, adonis2(data = meta.microb, method = "horn", permutations = 999,
#                                             strata = site, #use strata for nested (e.g., block) designs
#                                             formula = dist_usvi_asv.d ~ sampling_time:sampling_day + sampling_time + sampling_day, #Pr for sampling_time = 0.001, sampling_day = 0.001; sampling_time:sampling_day = 0.282
#                                             parallel = nthreads, by = "terms"))

# adonis_asv.res <- with(meta.microb, adonis2(data = meta.microb, method = "horn", permutations = 999,
#                                             # strata = sampling_time,formula = dist_usvi_asv.d ~ site*sampling_day,  #Pr: site=0.001, sampling_time=0.001, site:sampling_time=0.530
#                                             strata = sampling_day,formula = dist_usvi_asv.d ~ site*sampling_time,  #Pr: site=0.001, sampling_time=0.001, site:sampling_time=0.530
#                                             parallel = nthreads, by = "terms"))

adonis_asv.res
summary(adonis_asv.res)

# adonis_genus.res <- vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
#                                  # formula = dist_usvi_asv.d ~ site*sampling_time,
#                                  # formula = dist_usvi_asv.d ~ site*sampling_time*sampling_day + PAR*temp,
#                                  # formula = dist_usvi_asv.d ~  fcm_Synechococcus + fcm_Prochlorococcus + fcm_Picoeukaryotes + fcm_Unpigmented_cells,
#                                  formula = dist_usvi_genus.d ~ site*sampling_time*sampling_day,
#                                  # formula = dist_usvi_asv.d ~ fcm_Synechococcus*fcm_Prochlorococcus*fcm_Picoeukaryotes*fcm_Unpigmented_cells,
#                                  parallel = nthreads, by = "terms")
# summary(adonis_genus.res)


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

# pn_usvi_asv.mh <- pn_usvi_asv.mat %>%
pn_usvi_asv.mh <- usvi_sw_asv.tbl %>%
  dplyr::select(rownames(meta.microb)) %>%
  t() %>%
  vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE) %>%
  as.matrix(.)

# PERMANOVA::PerMANOVA.Simple(pn_usvi_asv.mh, meta.microb$site, nperm = 999, seed = 48105)

pn_usvi_asv.mh.list <- list(t(usvi_sw_asv.tbl), pn_usvi_asv.mh, "Morisita-Horn") %>%
  setNames(., c("Data", "D", "Coefficient"))
pn_usvi_asv.res3 <- PERMANOVA::PERMANOVA(pn_usvi_asv.mh.list, 
                                         group = meta.microb$grouping, #6 contrasts: site.sampling_time
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv.res1 <- PERMANOVA::PERMANOVA(pn_usvi_asv.mh.list, 
                                         group = meta.microb$site, 
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.mh.list, 
                                         group = meta.microb$sampling_time,
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv.res4 <- PERMANOVA::PERMANOVA(pn_usvi_asv.mh.list, 
                                         group = meta.microb$grouping2, #30 contrasts: site.sampling_time.sampling_day
                                         seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv.res1$pvalue #0.001
pn_usvi_asv.res2$pvalue #0.01
pn_usvi_asv.res3$pvalue #0.001
pn_usvi_asv.res4$pvalue #0.001



#try bray-curtis
pn_usvi_asv.bc <- pn_usvi_asv.mat %>%
  PERMANOVA::DistContinuous(., coef = "Bray_Curtis")

#for PERMANOVA,
#can't use in option 'group' one of the Effects' levels
pn_usvi_asv.res <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc, 
                                        group = meta.microb$grouping, #6 contrasts: site.sampling_time
                                        seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv.res$pvalue #0.017

pn_usvi_asv_day.res <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc, 
                                        group = meta.microb$grouping2, #30 contrasts: site.sampling_time.sampling_day
                                        # group = meta.microb$grouping, #6 contrasts: site.sampling_time
                                        seed = 48105, nperm = 999, PostHoc = "fdr")
pn_usvi_asv_day.res$pvalue #0.02


# #going deeper into PERMANOVA package:
# 
# #try buildingthe contrast matrix from scratch:
# contrast_mat_order <- meta.microb %>%
#   dplyr::arrange(site, sampling_time) %>%
#   dplyr::select(grouping) %>%
#   dplyr::mutate(label = grouping) %>%
#   dplyr::mutate(label = as.numeric(label)) %>%
#   droplevels
# 
# treatments <- factor(contrast_mat_order[["label"]], labels = unique(contrast_mat_order[["grouping"]]))
# contrasts(treatments) <- meta.microb[rownames(contrast_mat_order),] %>%
#   dplyr::select(site, sampling_time, grouping) %>%
#   dplyr::distinct(.) %>%
#   dplyr::mutate(value = 1) %>%
#   tidyr::pivot_wider(., id_cols = NULL,
#                      names_from = c("grouping"),
#                      values_from = "value",
#                      values_fill = 0) %>%
#   dplyr::mutate(site = as.numeric(site) - 1,
#                 # grouping =as.numeric(grouping),
#                 sampling_time = as.numeric(sampling_time) -1) %>%
#   droplevels %>%
#   dplyr::distinct(.) %>%
#   as.matrix(.)
# 
# 
# #build 0 intercept model and use limma to make contrast matrix
# mod_mat <- model.matrix(~0 + treatments)
# colnames(mod_mat) <- levels(treatments)
# rownames(mod_mat) <- rownames(contrast_mat_order)
# 
# library(limma)
# cmtx <- limma::makeContrasts( "LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
#                               "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
#                               "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo",
#                               levels= mod_mat )
# 
# pn_C <- t(cmtx) #this is so it is amenable to PERMANOVA functions
# 
# interaction_effects <- data.frame(label = c("LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
#                                             "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
#                                             "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo"),
#                                   to_test = c("sampling_time", "sampling_time", "sampling_time",
#                                               "site", "site", "site",
#                                               "site", "site", "site")) %>%
#   dplyr::arrange(rownames(pn_C)) %>%
#   dplyr::mutate(relabel = factor(label),
#                 retest = factor(to_test)) %>%
#   dplyr::mutate(retest = as.numeric(retest)) %>%
#   dplyr::mutate(relabel = as.numeric(relabel))
# # pn_Effects <- factor(interaction_effects[["relabel"]], labels = unique(interaction_effects[["label"]]))
# pn_Effects <- factor(interaction_effects[["retest"]], labels = unique(interaction_effects[["to_test"]]))
# 
# data_mat <- usvi_sw_asv.tbl %>% 
#   dplyr::select(rownames(contrast_mat_order)) %>%
#   t() %>%
#   PERMANOVA::IniTransform(., transform = 4) %>% #column centering: columns are ASVs, remove the column means
#   PERMANOVA::DistContinuous(., coef = "Bray_Curtis")
# 
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(data_mat,
#                                          group = contrast_mat_order$grouping, 
#                                          C = pn_C, 
#                                          Effects = pn_Effects,
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# summary(pn_usvi_asv.res2)
# 
# 
# #try MANOVA?
# mn_usvi_asv.res <- PERMANOVA::MANOVA(data_mat[["Data"]], Group = contrast_mat_order$grouping, 
#                   C = pn_C, 
#                   Effects = pn_Effects,
#                   InitialTransform = 4, Contrasts = TRUE)
# summary(mn_usvi_asv.res)
# 
# 
# 
# 
# # pn_usvi_asv.contrast <- tibble::tribble(
# #   ~LB_seagrass, ~Tektite, ~Yawzi,
# #   
# # )
# # pn_usvi_asv.contrast <- contr.sum(levels(meta.microb$grouping), contrasts = FALSE, sparse = FALSE)
# # colnames(pn_usvi_asv.contrast) <- levels(meta.microb$grouping)
# cH <- contr.sum(levels(meta.microb$grouping), contrasts = TRUE, sparse = FALSE)
# # apply(cH, 2, sum)
# # crossprod(cH)
# 
# 
# pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_time", "sampling_day")], MaxOrderIter = 3)
# #what does the contrast matrix look like?
# # #pulling out the contrasts matrix doesn't work as well:
# # pn_C <- pn_usvi_asv.contrast[["Contrasts"]]
# # colnames(pn_C) <- levels(pn_usvi_asv.contrast[["Groups"]])
# # rownames(pn_C) <- pn_usvi_asv.contrast[["Effects"]]
# # apply(pn_C, 2, sum)
# # # pn_Effects <- factor(c(1:6))
# # # levels(pn_Effects) <- levels(pn_usvi_asv.contrast[["Effects"]])
# # pn_Effects <- levels(pn_usvi_asv.contrast[["Effects"]])
# 
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          group = meta.microb$grouping2, #site.sampling_time.sampling_day
#                                          # group = meta.microb$site, #can't use in group one of the Effects levels
#                                          C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# summary(pn_usvi_asv.res2)
# # Effects 
# # Explained     Residual df Num df Denom       F-exp p-value p-value adj.
# # site                       -9.682771e+37 1.147236e+35      2       58  -24476.244   0.886        0.993
# # sampling_time               7.077761e+37 1.147236e+35      1       58   35782.528   0.106        0.318
# # sampling_day                4.293153e+38 1.147236e+35      4       58   54261.463   0.038        0.228
# # site*sampling_time         -5.729072e+38 1.147236e+35      2       58 -144820.284   0.993        0.993
# # site*sampling_day           8.990167e+37 1.147236e+35      8       58    5681.367   0.193        0.386
# # sampling_time*sampling_day -8.095247e+37 1.147236e+35      4       58  -10231.639   0.867        0.993
# # Total                       4.106122e+38 1.147236e+35     21       58    9885.249   0.067        0.067
# 
# pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_time")], MaxOrderIter = 3)
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          group = meta.microb$grouping, #site.sampling_time
#                                          C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# # Effects 
# # Explained     Residual df Num df Denom      F-exp p-value p-value adj.
# # site               -2.840383e+38 1.384244e+39      2       82  -8.412945   0.563        0.658
# # sampling_time      -3.995956e+38 1.384244e+39      1       82 -23.671288   0.658        0.658
# # site*sampling_time -2.898833e+38 1.384244e+39      2       82  -8.586069   0.621        0.658
# # Total              -9.735171e+38 1.384244e+39      5       82 -11.533863   0.845        0.845
# 
# #when using ConstructContrasts to build the effects levels, they seem only to be able to design pairs:
# pn_usvi_asv.contrast[["Effects"]]
# #sampling_time sampling_day site sampling_time*sampling_day sampling_time*site sampling_day*site
# 
# 
# pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("grouping", "sampling_day")], MaxOrderIter = 3)
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          group = meta.microb$grouping2, 
#                                          C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# # Effects 
# # Explained     Residual df Num df Denom      F-exp p-value p-value adj.
# # grouping              -5.784676e+38 1.147236e+35      5       58 -58490.336   0.968        0.968
# # sampling_day           4.082913e+38 1.147236e+35      4       58  51604.217   0.042        0.105
# # grouping*sampling_day  5.807885e+38 1.147236e+35     20       58  14681.252   0.070        0.105
# # Total                  4.106122e+38 1.147236e+35     29       58   7158.284   0.067        0.067
# 
# 
# 
# pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("sampling_time", "sampling_day", "site")], MaxOrderIter = 3)
# # #pulling out the contrasts matrix doesn't work as well:
# # pn_C <- pn_usvi_asv.contrast[["Contrasts"]]
# # colnames(pn_C) <- levels(pn_usvi_asv.contrast[["Groups"]])
# # rownames(pn_C) <- pn_usvi_asv.contrast[["Effects"]]
# # # pn_Effects <- factor(c(1:6))
# # # levels(pn_Effects) <- levels(pn_usvi_asv.contrast[["Effects"]])
# # pn_Effects <- levels(pn_usvi_asv.contrast[["Effects"]])
# # pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
# #                                          seed = 48105, nperm = 999, PostHoc = "fdr",
# #                                          # C = pn_C, Effects = pn_Effects,
# #                                          group = meta.microb$grouping)
# 
# contrasts(meta.microb$grouping, contrasts = FALSE)
# 
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          group = meta.microb$grouping2,
#                                          C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# # Effects 
# # Explained     Residual df Num df Denom       F-exp p-value p-value adj.
# # sampling_time              -3.875917e+38 1.147236e+35      1       58 -195951.925   0.984        0.984
# # sampling_day                3.719104e+38 1.147236e+35      4       58   47006.019   0.053        0.152
# # site                       -3.006977e+38 1.147236e+35      2       58  -76010.779   0.965        0.984
# # sampling_time*sampling_day  4.008087e+38 1.147236e+35      4       58   50658.489   0.031        0.152
# # sampling_time*site         -2.898182e+38 1.147236e+35      2       58  -73260.636   0.967        0.984
# # sampling_day*site           3.384847e+38 1.147236e+35      8       58   21390.660   0.076        0.152
# # Total                       4.106122e+38 1.147236e+35     21       58    9885.249   0.067        0.067
# 
# 
# 
# pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_time")], MaxOrderIter = 3)
# # pn_usvi_asv.contrast <- PERMANOVA::ConstructContrasts(meta.microb[, c("site", "sampling_day")], MaxOrderIter = 3)
# # crossprod(pn_usvi_asv.contrast[["Contrasts"]])
# # apply(pn_usvi_asv.contrast[["Contrasts"]], 2, sum)
# 
# pn_usvi_asv.res2 <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          group = meta.microb$grouping, 
#                                          C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# 
# pn_usvi_asv.res <- PERMANOVA::PERMANOVA(pn_usvi_asv.bc,
#                                          group = meta.microb$sampling_time, 
#                                         # C = pn_usvi_asv.contrast[["Contrasts"]], Effects = pn_usvi_asv.contrast[["Effects"]],
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# 
# 
# summary(pn_usvi_asv.res2)
# 
# 
# # mn_usvi_asv.res <- PERMANOVA::MANOVA(Y = pn_usvi_asv.bc[["D"]], InitialTransform = 4,
# #                      Group = meta.microb$grouping, 
# #                      C = pn_usvi_asv.contrast[["Contrasts"]], 
# #                      # M = 
# #                      Contrasts = TRUE)
# # summary(mn_usvi_asv.res)
# 
# 
# #construct a matrix for contrasts:
# # treatments <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6),labels=c("LB_seagrass.dawn", "LB_seagrass.peak_photo", "Tektite.dawn", "Tektite.peak_photo", "Yawzi.dawn", "Yawzi.peak_photo"))
# # contrasts(treatments) <- cbind(site=c(0,0,1,1,2,2), sampling_time=c(0,1,0,1,0,1), `LB_seagrass.peak_photo` = c(0,1,0,0,0,0), `Tektite.peak_photo`=c(0,0,0,1,0,0), `Yawzi.peak_photo`=c(0,0,0,0,0,1))
# # mod_mat <- model.matrix(~treatments)
# # colnames(mod_mat) <- c("Intercept","site","sampling_time","LB_seagrass.peak_photo","Tektite.peak_photo","Yawzi.peak_photo")
# # # rownames(mod_mat) <- as.character(treatments) #this is the order of the rows, but assigning rownames is not helpful
# 
# #include the days, not that it really matters:
# # treatments <- factor(c(rep(1, 5),rep(2, 5),rep(3, 5),rep(4, 5),rep(5, 5),rep(6,5)), labels=c("LB_seagrass.dawn", "LB_seagrass.peak_photo", "Tektite.dawn", "Tektite.peak_photo", "Yawzi.dawn", "Yawzi.peak_photo"))
# # contrasts(treatments) <- cbind(site=c(0,0,1,1,2,2), sampling_time=c(0,1,0,1,0,1), `LB_seagrass.peak_photo` = c(0,1,0,0,0,0), `Tektite.peak_photo`=c(0,0,0,1,0,0), `Yawzi.peak_photo`=c(0,0,0,0,0,1))
# 
# # treatments <- factor(as.numeric(meta.microb[["grouping"]]), labels = unique(meta.microb[["grouping"]]))
# contrast_mat_order <- meta.microb %>%
#   dplyr::arrange(site, sampling_time) %>%
#   dplyr::select(grouping) %>%
#   dplyr::mutate(label = grouping) %>%
#   dplyr::mutate(label = as.numeric(label)) %>%
#   droplevels
# treatments <- factor(contrast_mat_order[["label"]], labels = unique(contrast_mat_order[["grouping"]]))
# 
# contrasts(treatments) <- meta.microb[rownames(contrast_mat_order),] %>%
# # contrasts(treatments) <- meta.microb %>%
# #   tibble::rownames_to_column(var = "sample_id") %>%
#   # dplyr::arrange(site, sampling_time, sample_id) %>%
#   dplyr::select(site, sampling_time, grouping) %>%
#   dplyr::distinct(.) %>%
#   dplyr::mutate(value = 1) %>%
#   tidyr::pivot_wider(., id_cols = NULL,
#                      names_from = c("grouping"),
#                      values_from = "value",
#                      values_fill = 0) %>%
#   dplyr::mutate(site = as.numeric(site) - 1,
#                 # grouping =as.numeric(grouping),
#                 sampling_time = as.numeric(sampling_time) -1) %>%
#   droplevels %>%
#   dplyr::distinct(.) %>%
#   as.matrix(.)
# 
# 
# #0 intercept
# mod_mat <- model.matrix(~0 + treatments)
# colnames(mod_mat) <- levels(treatments)
# rownames(mod_mat) <- rownames(contrast_mat_order)
# 
# library(limma)
# cmtx <- limma::makeContrasts( "LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
#                               "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
#                               "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo",
#                               levels= mod_mat )
# # data_mat <- dist_usvi_asv.mat[rownames(contrast_mat_order),]
# data_mat <- usvi_sw_asv.tbl %>%
#   dplyr::select(rownames(contrast_mat_order)) %>%
#   apply(., 2, function(x) log2(x + 1))
# 
# fit <- limma::lmFit(data_mat, mod_mat) %>%
# # fit <- limma::lmFit(data_mat, ) %>%
#   limma::contrasts.fit(., contrasts=cmtx, coefficients=NULL)
# fit2 <- eBayes(fit)
# topTable(fit2,coef=2)
# volcanoplot(fit2,coef=2,highlight=2)
# 
# limma.res <- decideTests(fit2, method="nestedF", adjust.method="fdr", p.value = 0.05)
# summary(limma.res)
# 
# 
# 
# pn_Effects <- factor(colnames(cmtx)) %>%as.numeric
# levels(pn_Effects) <- colnames(cmtx)
# # pn_Effects <- factor(c(1, 2, 3))
# # levels(pn_Effects) <- c("site", "sampling_time", "interaction")
# 
# data_mat <- usvi_sw_asv.tbl %>% 
#   dplyr::select(rownames(contrast_mat_order)) %>%
#   t() %>%
#   PERMANOVA::IniTransform(., transform = 4) %>% #column centering: columns are ASVs, remove the column means
#   PERMANOVA::DistContinuous(., coef = "Bray_Curtis")
# 
# cmtx <- ConstructContrasts(meta.microb[rownames(contrast_mat_order), c("site", "sampling_time")], MaxOrderIter=3)
# pn_usvi_asv.res3 <- PERMANOVA::PERMANOVA(data_mat,
#                                          group = contrast_mat_order$grouping, 
#                                          C = cmtx[["Contrasts"]], Effects=cmtx[["Effects"]],
#                                          # C = t(cmtx), 
#                                          # Effects = pn_Effects,
#                                          seed = 48105, nperm = 999, PostHoc = "fdr")
# summary(pn_usvi_asv.res3)
# 
# 
# 
# 
# # #construct a 5row x 6col matrix for each contrast to compare:
# # mod_mat <- model.matrix(~0 + levels(meta.microb$grouping))
# # colnames(mod_mat) <- c("LB_seagrass.dawn", "LB_seagrass.peak_photo", "Tektite.dawn", "Tektite.peak_photo", "Yawzi.dawn", "Yawzi.peak_photo")
# # # mod_mat <- model.matrix(~levels(meta.microb$grouping))
# # # colnames(mod_mat) <- c("Intercept", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo", "Yawzi.peak_photo")
# # 
# # rownames(mod_mat) <- c("LB_seagrass.dawn", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo", "Yawzi.peak_photo")
# # 
# # c("LB_seagrass.dawn", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo", "Yawzi.peak_photo")
# # cH <- tibble::tribble(
# #   ~"LB_seagrass.dawn", ~"Tektite.dawn", ~"Yawzi.dawn", ~"LB_seagrass.peak_photo", ~"Tektite.peak_photo", ~"Yawzi.peak_photo",
# #   1, 0, 0, 0, 0, 0,
# #   0, 1, 0, 0, 0, 0,
# #   0, 0, 1, 0, 0, 0,
# #   0, 0, 0, 1, 0, 0,
# #   0, 0, 0, 0, 1, 0
# # ) %>%
# #   as.matrix(.)
# # rownames(cH) <- c("LB_seagrass.dawn", "Tektite.dawn", "Yawzi.dawn", "LB_seagrass.peak_photo", "Tektite.peak_photo")
# # 
# # # 
# # pn_usvi_asv.contrast <- contr.treatment(levels(meta.microb$grouping), base = 1, contrasts = TRUE, sparse = FALSE)
# 
# # mod_mat <- model.matrix(~0 + site + sampling_time + sampling_day + site:sampling_time:sampling_day,
# #                         data = meta.microb)
# # if(length(which(colSums(mod_mat) == 0)) > 0){ #none of these should be 0--if so, you need to fix your design/metadata coding
# #   idx <- names(colSums(mod_mat))[which(colSums(mod_mat) == 0)]
# #   idx <- cli::cli_vec(idx)
# #   cli::cli_alert_warning("One or more groups in the comparison does not have enough replicates:", wrap = TRUE)
# #   cli::cli_bullets_raw(idx)
# # }
# 
# # condition_filter <- list(`site` = c("LB_seagrass", "Yawzi", "Tektite"),
# #                          `sampling_time` = c("dawn", "peak_photo"),
# #                          `sampling_day` = c("Day1", "Day2", "Day3", "Day4", "Day5"))
# # 
# # condition_idx_names <- purrr::map(condition_filter,
# #                                   ~unlist(., recursive = FALSE)) %>%
# #   purrr::reduce(interaction, sep = "_") %>%
# #   levels %>%
# #   as.character
# # 
# # #brute force: list the model vectors
# # for(i in c("LB_seagrass", "Yawzi", "Tektite")){
# #   namevar <- i
# #   for(j in c("dawn", "peak_photo")){
# #     lightvar <- j
# #     if(is.null(dim(mod_mat[meta.microb$site == namevar & meta.microb$sampling_time == lightvar,]))){
# #       temp_vec <- mod_mat[meta.microb$site == namevar & meta.microb$sampling_time == lightvar, ]
# #     } else {
# #       temp_vec <- colMeans(mod_mat[meta.microb$site == namevar & meta.microb$sampling_time == lightvar, ])
# #     }
# #     assign(paste0("contrast_", namevar, "_", lightvar), temp_vec, envir = .GlobalEnv)
# #     rm(temp_vec)
# #     rm(lightvar)
# #   }
# #   rm(namevar)
# # }
# # 
# # contrast_id <- condition_idx_names %>%
# #   as.data.frame() %>%
# #   dplyr::mutate(baseline = "baseline_vec") %>%
# #   dplyr::mutate(across(everything(), as.character))
# 
# # pn_usvi.effects <- factor(c(1,2,3))
# # levels(pn_usvi.effects) <- c("site", "sampling_time", "sampling_day")
# 
# 
# 
# 
# ###tngent:
# #using Limma to detect SDA ASVs between pairs
# 
# contrast_mat_order <- meta.microb %>%
#   dplyr::arrange(site, sampling_time) %>%
#   dplyr::select(grouping) %>%
#   dplyr::mutate(label = grouping) %>%
#   dplyr::mutate(label = as.numeric(label)) %>%
#   droplevels
# treatments <- factor(contrast_mat_order[["label"]], labels = unique(contrast_mat_order[["grouping"]]))
# 
# contrasts(treatments) <- meta.microb[rownames(contrast_mat_order),] %>%
#   dplyr::select(site, sampling_time, grouping) %>%
#   dplyr::distinct(.) %>%
#   dplyr::mutate(value = 1) %>%
#   tidyr::pivot_wider(., id_cols = NULL,
#                      names_from = c("grouping"),
#                      values_from = "value",
#                      values_fill = 0) %>%
#   dplyr::mutate(site = as.numeric(site) - 1,
#                 sampling_time = as.numeric(sampling_time) -1) %>%
#   droplevels %>%
#   dplyr::distinct(.) %>%
#   as.matrix(.)
# 
# 
# #0 intercept
# mod_mat <- model.matrix(~0 + treatments)
# colnames(mod_mat) <- levels(treatments)
# rownames(mod_mat) <- rownames(contrast_mat_order)
# 
# library(limma)
# cmtx <- limma::makeContrasts( "LB_seagrass.peak_photo-LB_seagrass.dawn", "Tektite.peak_photo-Tektite.dawn", "Yawzi.peak_photo-Yawzi.dawn",
#                               "Tektite.dawn-LB_seagrass.dawn", "Yawzi.dawn-LB_seagrass.dawn", "Yawzi.dawn-Tektite.dawn",
#                               "Tektite.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-LB_seagrass.peak_photo", "Yawzi.peak_photo-Tektite.peak_photo",
#                               levels= mod_mat )
# 
# data_mat <- usvi_sw_asv.tbl %>%
#   dplyr::select(rownames(contrast_mat_order)) %>%
#   apply(., 2, function(x) log2(x + 1))
# 
# fit <- limma::lmFit(data_mat, mod_mat) %>%
#   limma::contrasts.fit(., contrasts=cmtx, coefficients=NULL)
# fit2 <- eBayes(fit)
# topTable(fit2,coef=2)
# volcanoplot(fit2,coef=2,highlight=2)
# 
# limma.res <- decideTests(fit2, method="nestedF", adjust.method="fdr", p.value = 0.05)
# summary(limma.res)
# 


