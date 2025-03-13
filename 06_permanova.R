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
    dplyr::mutate(contrast = paste0(pair1, ":", pair2)) %>%
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
  # ttest_res %>% dplyr::select(p.value) %>% tibble::deframe(.)
  ttest_res %>% dplyr::select(contrast, p.value) %>% tibble::deframe(.)
}

F_ttest_rn <- function(dataset, mle){
  #mle is a required component for the ran.gen argument in boot::boot
  #but know that it is functionally equivalent to "resample_depth" calculated in F_ttest_gn
  #dataset is agrouped dataset, by "variable", so that you correctly resample
  resample_depth <- mle
  
  out <- dataset %>%
    dplyr::slice_sample(., n = resample_depth) %>%
    droplevels
  out
}



F_ttest_para <- function(dataset, variable){ 
  #this version is for parametric simulation using resampling and correctly reports the full sampleset p-value
  #d is a dataframe containing all the measurements by site
  #variable is your factor you want to examine pairwise the levels
  variable <- rlang::parse_expr(variable)
  temp_vars <- unique(dataset[[variable]])
  
  ttest_res <- combn(temp_vars, 2) %>%
    t() %>%
    as.data.frame() %>%
    setNames(., c("pair1", "pair2")) %>%
    dplyr::mutate(p.value = NA)
  for(i in seq_len(nrow(ttest_res))){
    var1 <- ttest_res[i, 1]
    var2 <- ttest_res[i, 2]
    temp_i_a <- dataset %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
      dplyr::select(dissimilarity) %>%
      tibble::deframe(.)
    temp_i_b <- dataset %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
      dplyr::select(dissimilarity) %>%
      tibble::deframe(.)
    ttest_res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
  }
  
  ttest_res %>% dplyr::select(p.value) %>% tibble::deframe(.)
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

# Read in metabolites -----------------------------------------------------

if(file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
  temp_list <- readr::read_rds(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
  list2env(temp_list, envir = .GlobalEnv)
  rm(temp_list)
} else {
  usvi_metabolomics_long.df <- readr::read_delim(paste0(projectpath, "/", 
                                                        # "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_longFormat_wMetadata.csv"), 
                                                        "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_longFormat_outliersRmvd_wMetadata_20250204.csv"),
                                                 col_select = c(2:last_col()),
                                                 col_names = TRUE, show_col_types = FALSE, delim = ",", num_threads = nthreads) %>%
    dplyr::select(-"Sample_ID") %>%
    dplyr::mutate(sample_id = paste0("Metab_", DNA_no)) %>%
    dplyr::select(metabolites, adaptedDervLabel, Site, Date, Time, starts_with("sampl", ignore.case = TRUE), concentration, LOD, LOQ, contains("flag")) %>%
    dplyr::mutate(LODflag = dplyr::case_when((LODflag == FALSE) ~ 0, (LODflag == TRUE) ~ 1, .default = NA)) %>% #"LODflag" == FALSE means that the metabolite was measured at a concentration above the LOD (so it's good)
    droplevels
  
  usvi_metabolomics.df <- usvi_metabolomics_long.df %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration) %>%
    tidyr::pivot_wider(., id_cols = "adaptedDervLabel",
                       names_from = "metabolites",
                       values_fill = NA,
                       values_from = "concentration") %>%
    dplyr::rename(metab_deriv_label = "adaptedDervLabel")
  
  usvi_metabolomics_old_long.df <- readr::read_delim(paste0(projectpath, "/", 
                                                            "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_longFormat_wMetadata.csv"),
                                                     # "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_longFormat_outliersRmvd_wMetadata_20250204.csv"),
                                                     col_select = c(2:last_col()),
                                                     col_names = TRUE, show_col_types = FALSE, delim = ",", num_threads = nthreads) %>%
    dplyr::select(-"Sample_ID") %>%
    dplyr::mutate(sample_id = paste0("Metab_", DNA_no)) %>%
    dplyr::select(metabolites, adaptedDervLabel, Site, Date, Time, starts_with("sampl", ignore.case = TRUE), concentration, LOD, LOQ, contains("flag")) %>%
    dplyr::mutate(LODflag = dplyr::case_when((LODflag == TRUE) ~ 0, (LODflag == FALSE) ~ 1, .default = NA)) %>% #"LODflag" == TRUE means that the metabolite was measured at a concentration above the LOD (so it's good)
    droplevels
  
  usvi_metab_cinar_bc_73.df <- usvi_metabolomics_old_long.df %>%
    dplyr::ungroup(.) %>%
    dplyr::rowwise(.) %>%
    dplyr::filter(grepl("CINAR_BC_73", adaptedDervLabel)) %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::filter(metabolites %in% unique(usvi_metabolomics_long.df[["metabolites"]])) %>%
    droplevels %>%
    dplyr::right_join(., usvi_metabolomics_long.df %>%
                        dplyr::distinct(metabolites, LOD, LOQ) %>% 
                        # dplyr::rename(`LOD.new` = "LOD") %>%
                        droplevels,
                      by = join_by(metabolites), suffix = c(".old", ".new")) %>%
    dplyr::ungroup(.) %>%
    tidyr::fill(adaptedDervLabel, .direction = "down") %>%
    dplyr::mutate(LODflag = dplyr::case_when(is.na(LODflag) ~ 1, .default = LODflag),
                  concentration = dplyr::case_when((LODflag == 1) ~ LOD/10, .default = concentration)) %>%
    droplevels
  
  if(!file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
    temp_list <- list(usvi_metabolomics.df, usvi_metabolomics_long.df, usvi_metabolomics_old_long.df, usvi_metab_cinar_bc_73.df) %>%
      setNames(., c("usvi_metabolomics.df", "usvi_metabolomics_long.df", "usvi_metabolomics_old_long.df", "usvi_metab_cinar_bc_73.df"))
    readr::write_rds(temp_list, paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
    rm(temp_list)
  }
}




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
# g1_par <- print(
#   ggplot(data = usvi_hobo_light_temp_filtered.df %>%
#            dplyr::filter(grepl("PAR", parameter)),
#          aes(x = date_ast, y = value, fill = site, color = site, group = interaction(site, parameter)))
#   + geom_point(shape = 19, size = 1)
#   + geom_line(show.legend = FALSE)
#   + theme_bw()
#   + facet_grid(site~., labeller = labeller(site = site_lookup),
#                scales = "fixed")
#   + scale_y_continuous(name = expression(paste("PAR (\U00B5 mol photons", ~m^-2 ~s^-1, ")")))
#   + scale_discrete_manual(aesthetics = c("color"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
#                           drop = TRUE)
#   + scale_discrete_manual(aesthetics = c("fill"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
#                           drop = TRUE)
#   + theme(axis.text.x = element_text(angle = 90),
#           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#           panel.grid = element_blank(),
#           axis.title.x = element_blank(),
#           legend.position = "none",
#           legend.key = element_blank(),
#           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           legend.text = element_text(size = 12, colour = "grey30"))
# )
# g1_temp <- print(
#   ggplot(data = usvi_hobo_light_temp_filtered.df %>%
#            dplyr::filter(grepl("temp", parameter)),
#          aes(x = date_ast, y = value, fill = site, color = site, group = interaction(site, parameter)))
#   + geom_point(shape = 19, size = 1)
#   + geom_line(show.legend = FALSE)
#   + theme_bw()
#   + facet_grid(site~., labeller = labeller(site = site_lookup),
#                scales = "fixed")
#   + scale_y_continuous(name = expression(paste("Temperature (˚C)")))
#   + scale_discrete_manual(aesthetics = c("color"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
#                           drop = TRUE)
#   + scale_discrete_manual(aesthetics = c("fill"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),
#                           drop = TRUE)
#   + theme(axis.text.x = element_text(angle = 90),
#           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#           panel.grid = element_blank(),
#           axis.title.x = element_blank(),
#           legend.position = "none",
#           legend.key = element_blank(),
#           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           legend.text = element_text(size = 12, colour = "grey30"))
# )
# 
# gpatch <- g1_par + g1_temp + patchwork::plot_annotation(title = "Physicochemical parameters", subtitle ="Measured in USVI sites via HOBO loggers", tag_levels = "A")
# gpatch
# ggsave(paste0(projectpath, "/", "usvi_light_temp-", Sys.Date(), ".png"),
#        gpatch,
#        width = 10, height = 8, units = "in")


#look at specifically 5-7am and 12-2pm

interval <- lubridate::as.interval(duration(hours = 1), ymd_hms("2021-01-21 13:00:00", tz = "America/Virgin")) 
interval_vector_peak <- list(interval,
                             lubridate::int_shift(interval, days(1)),
                             lubridate::int_shift(interval, days(2)),
                             lubridate::int_shift(interval, days(3)),
                             lubridate::int_shift(interval, days(4)),
                             lubridate::int_shift(interval, days(5)))

interval <- lubridate::as.interval(duration(hours = 1), ymd_hms("2021-01-21 06:00:00", tz = "America/Virgin")) 
interval_vector_dawn <- list(interval,
                             lubridate::int_shift(interval, days(1)),
                             lubridate::int_shift(interval, days(2)),
                             lubridate::int_shift(interval, days(3)),
                             lubridate::int_shift(interval, days(4)),
                             lubridate::int_shift(interval, days(5)))


usvi_hobo_light_temp_sampling.df <- usvi_hobo_light_temp_filtered.df %>%
  dplyr::filter(grepl("temp|PAR", parameter)) %>%
  tidyr::drop_na(sampling_time) %>%
  dplyr::mutate(sampled_sw = dplyr::case_when(date_ast %within% interval_vector_peak ~ 1,
                                                 date_ast %within% interval_vector_dawn ~ 1,
                                                 .default = NA)) %>%
  dplyr::mutate(sampled_sw = sampled_sw*min(value,na.rm = TRUE), .by = parameter) %>%
  droplevels

g1_par <- print(
  ggplot(data = usvi_hobo_light_temp_filtered.df %>%
           dplyr::filter(grepl("PAR", parameter)),
         aes(x = date_ast, y = value, group = interaction(site, parameter), color = "grey"))
  + geom_smooth(show.legend = FALSE, method = "loess", span = 0.1)
  + geom_point(data = usvi_hobo_light_temp_sampling.df %>%
                 dplyr::filter(grepl("PAR", parameter)),
               aes(x = date_ast, y = value, fill = site, color = site, group = interaction(site, parameter)), shape = 19, size = 1)
  + geom_crossbar(data = usvi_hobo_light_temp_sampling.df %>%
                     dplyr::filter(grepl("PAR", parameter)) %>%
                     tidyr::drop_na(sampled_sw),
                   aes(xmin= date_ast, xmax = date_ast, y = sampled_sw, group = interaction(day, site, parameter, sampling_time)), color = "black", linewidth = 1, alpha = 1) 
  + theme_bw()
  + facet_grid(site~., labeller = labeller(site = site_lookup),
               scales = "free_y")
  + scale_y_continuous(name = expression(paste("PAR (\U00B5 mol photons", ~m^-2 ~s^-1, ") at sampling depth")))
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
  ggplot(data = usvi_hobo_light_temp_sampling.df %>%
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

usvi_metabolomics.tbl <- usvi_metabolomics_long.df %>%
  dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
  # dplyr::bind_rows(., usvi_metab_cinar_bc_73.df %>%
  #                    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
  #                    droplevels) %>%
  dplyr::rename(metabolite = "metabolites", metab_deriv_label = "adaptedDervLabel") %>%
  dplyr::filter(LODflag == 0) %>%
  dplyr::select(-LODflag) %>%
  dplyr::filter(!(metab_deriv_label %in% drop)) %>%
  tidyr::pivot_wider(., id_cols = "metab_deriv_label",
                     names_from = "metabolite",
                     # values_fill = 0,
                     values_from = "concentration") %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!(intersect(usvi_sus_metabolites_idx[["metabolites"]], colnames(.)))) %>%
  t() %>%
  as.data.frame(.)


usvi_selected_metadata <- metabolomics_sample_metadata %>%
  dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::filter(!grepl("Day1", sampling_day)) %>%
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

usvi_sw_asv.tbl <- usvi_prok_asvs.df %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  dplyr::filter(!grepl("Metab_280", sample_ID)) %>% #drop the sample corresponding to CINAR_BC_73: Metab_280
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
                grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day2")) %>%
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



# Permanova via adonis2 ---------------------------------------------------
#please see here for why we should not include continuous variables in adonis2's implementation of permanova:
#https://learninghub.primer-e.com/books/should-i-use-primer-or-r/chapter/3-permanova-vs-adonis2-in-r/export/html

#we have a Crossed Model, where Site has 3 levels, Sampling_time has 2 levels, and Sampling_day has 5 levels
#with this data set we really ought to use type III SS
#but adonis only does type I SS.

#what it used to be with Day1 samples:
{
  meta.microb_full <- metabolomics_sample_metadata %>%
    dplyr::select(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
    # dplyr::filter(!grepl("Day1", sampling_day)) %>%
    dplyr::mutate(across(c(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site), ~factor(.x))) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    dplyr::mutate(site_type = dplyr::case_when(grepl("LB", site) ~ "seagrass",
                                               grepl("Yawzi|Tektite", site) ~ "reef",
                                               .default = NA)) %>%
    # dplyr::mutate(sample_id = factor(sample_id, levels = colnames(usvi_sw_asv_full.tbl))) %>%
    dplyr::arrange(sample_id) %>%
    tidyr::drop_na(.) %>%
    dplyr::mutate(grouping = interaction(site, sampling_time)) %>%
    dplyr::mutate(grouping2 = paste0(grouping, ".", sampling_day)) %>%
    dplyr::mutate(across(c(sample_id, sample_type, site, sampling_time, sampling_day, metab_deriv_label, site_type, grouping, grouping2), ~factor(.x))) %>%
    dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                  grouping = fct_relevel(grouping, "LB_seagrass.dawn"),
                  grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day2")) %>%
    dplyr::mutate(rownames = sample_id) %>%
    tibble::column_to_rownames(var = "rownames") %>%
    droplevels
  
  usvi_sw_asv_full.tbl <- usvi_prok_asvs.df %>%
    dplyr::filter(sample_ID %in% meta.microb_full[["sample_id"]]) %>%
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
  
  dist_usvi_asv_full.mh.d <- usvi_sw_asv_full.tbl %>%
    dplyr::select(rownames(meta.microb_full)) %>%
    t() %>%
    vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE)
  
  temp_permanova_asv_est_var.df <- with(meta.microb_full, vegan::adonis2(data = meta.microb_full, method = "horn", permutations = 9999,
                                                                         # strata = site,
                                                                         formula = dist_usvi_asv_full.mh.d ~ sampling_time*sampling_day*site,
                                                                         parallel = nthreads, by = "terms"))
  temp_permanova_asv_est_var.df <- temp_permanova_asv_est_var.df %>%
    tibble::rownames_to_column(var = "term") %>%
    dplyr::mutate(term = as.character(term)) %>%
    dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                        .default = NA))
  temp_permanova_asv_est_var.df <- temp_permanova_asv_est_var.df %>%
    dplyr::mutate(temp_permanova_asv_est_var.df %>%
                    dplyr::filter(!grepl("Total", term)) %>%
                    dplyr::summarise(TotalEMS = sum(MS))) %>%
    dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
                                            .default = NA)) %>%
    dplyr::select(-TotalEMS) %>%
    dplyr::mutate(sq_root = MS^(1/2)) %>%
    dplyr::mutate(perc_variation = sq_root/(sum(sq_root, na.rm = TRUE))*100)
}

# #SS(total) for microbiomes: 1.239223
# (1/nrow(dist_usvi_asv.mat))*(sum((dist_usvi_asv.mat[lower.tri(dist_usvi_asv.mat, diag = FALSE)])^2))

#per-site Sum of Squares
dist_usvi_asv.ss <- meta.microb %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::select(sample_id) %>%
        tibble::deframe(.)) %>%
  map(., ~dist_usvi_asv.mat[.x, .x]) %>%
  map(., ~sum((.x[lower.tri(.x, diag = FALSE)])^2)/nrow(.x))
#LB: 0.2625108
#Tektite: 0.1379405
#Yawzi: 0.09376283


#this is the final one we should use:
#in this output table, we want the results for:
#site
#site:sampling_time
#site:sampling_day
#site:sampling_time:sampling_day



#do samplign day before ssite:
with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
                                 # strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time*sampling_day*site,
                                 parallel = nthreads, by = "terms"))
# Df SumOfSqs       R2       F Pr(>F)    
# sampling_time                    1  0.07469  0.06027 15.6865 0.0002 ***
#   sampling_day                     3  0.10838  0.08746  7.5874 0.0002 ***
#   site                             2  0.74539  0.60150 78.2766 0.0001 ***
#   sampling_time:sampling_day       3  0.04151  0.03350  2.9060 0.0464 *  
#   sampling_time:site               2 -0.02104 -0.01698 -2.2093 1.0000    
# sampling_day:site                6  0.07101  0.05730  2.4857 0.0358 *  
#   sampling_time:sampling_day:site  6 -0.00449 -0.00363 -0.1573 0.9979    
# Residual                        47  0.22378  0.18058                   
# Total                           70  1.23922  1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#note that now we have negative SS for two interactions: sampling_time:site, and sampling_time:sampling_day:site
#see here for what to do about negative SS: https://learninghub.primer-e.com/books/permanova-for-primer-guide-to-software-and-statistical-methods/page/129-pooling-or-excluding-terms

#is it negative if we included the Metab_280 sample?
#tldr: yes
#using BC distance metric does not generate neagtive SS for sampling_time:sampling_day and sampling_time:sampling_day:site
#but with Bray, sampling_time:sampling_day is not flagged as significant, whereas sampling_day:site is more strongly significant.

{
  # dist_usvi_asv_full.mh.d <- usvi_sw_asv_full.tbl %>%
  #   dplyr::select(rownames(meta.microb_full)) %>%
  #   t() %>%
  #   vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE)
# 
#   
  # 
  # #try a different distance metric:
  ###try BC
  # dist_usvi_asv_full.bc.d <- usvi_sw_asv_full.tbl %>% 
  #   dplyr::select(rownames(meta.microb_full)) %>%
  #   t() %>%
  #   vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)
  # 
  # 
  # with(meta.microb_full, vegan::adonis2(data = meta.microb_full, method = "bray", permutations = 9999,
  #                                       # strata = site,
  #                                       formula = dist_usvi_asv_full.bc.d ~ sampling_time*sampling_day*site,
  #                                       parallel = nthreads, by = "terms"))
  # # Df SumOfSqs       R2       F Pr(>F)    
  # # sampling_time                    1  0.06931  0.05341 14.7839 0.0004 ***
  # #   sampling_day                     3  0.11603  0.08941  8.2490 0.0001 ***
  # #   site                             2  0.79357  0.61150 84.6288 0.0001 ***
  # #   sampling_time:sampling_day       3  0.04353  0.03354  3.0945 0.0351 *  
  # #   sampling_time:site               2 -0.01896 -0.01461 -2.0219 1.0000    
  # # sampling_day:site                6  0.06956  0.05360  2.4728 0.0391 *  
  # #   sampling_time:sampling_day:site  6 -0.00035 -0.00027 -0.0126 0.9915    
  # # Residual                        48  0.22505  0.17342                   
  # # Total                           71  1.29774  1.00000                   
  # # ---
  # #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # # vegan::adonis2(formula = dist_usvi_asv_full.bc.d ~ sampling_time * sampling_day * site, data = meta.microb_full, permutations = 9999, method = "bray", by = "terms", parallel = nthreads)
  # # Df SumOfSqs      R2       F Pr(>F)    
  # # sampling_time                    1   0.1475 0.02954  7.0551 0.0024 ** 
  # #   sampling_day                     3   0.3526 0.07061  5.6214 0.0001 ***
  # #   site                             2   2.8066 0.56198 67.1123 0.0001 ***
  # #   sampling_time:sampling_day       3   0.1224 0.02450  1.9509 0.0819 .  
  # # sampling_time:site               2   0.0559 0.01119  1.3361 0.2333    
  # # sampling_day:site                6   0.3268 0.06543  2.6046 0.0055 ** 
  # #   sampling_time:sampling_day:site  6   0.1787 0.03577  1.4240 0.1593    
  # # Residual                        48   1.0037 0.20097                   
  # # Total                           71   4.9942 1.00000                   
  # # ---
  # #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # 
  # 
  # 
  
  # dist_usvi_asv.bc.d <- usvi_sw_asv.tbl %>% 
  #   dplyr::select(rownames(meta.microb)) %>%
  #   t() %>%
  #   vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE)
  # 
  # temp_res <- with(meta.microb, vegan::adonis2(data = meta.microb_full, method = "bray", permutations = 9999,
  #                                              # strata = site,
  #                                              # formula = dist_usvi_asv.bc.d ~ sampling_time*sampling_day*site,
  #                                              formula = dist_usvi_asv_full.bc.d ~ sampling_time*sampling_day*site,
  #                                              parallel = nthreads, by = "terms"))
  # permanova_asv_estimate_variation.bc.df <- temp_res %>%
  #   tibble::rownames_to_column(var = "term") %>%
  #   dplyr::mutate(term = as.character(term)) %>%
  #   dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
  #                                       .default = NA))
  # permanova_asv_estimate_variation.bc.df <- permanova_asv_estimate_variation.bc.df %>%
  #   dplyr::mutate(permanova_asv_estimate_variation.bc.df %>%
  #                   dplyr::filter(!grepl("Total", term)) %>%
  #                   dplyr::summarise(TotalEMS = sum(MS))) %>%
  #   dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
  #                                           .default = NA)) %>%
  #   dplyr::mutate(sq_root = MS^(1/2))
  # 
}



#if we manually specify the terms (not subtracting them):
{
  #if we set sampling_time:sampling_day:site to 0
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time * sampling_day * site - sampling_time:sampling_day:site, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs       R2       F Pr(>F)    
  # sampling_time               1  0.07469  0.06027 18.0516 0.0001 ***
  #   sampling_day                3  0.10838  0.08746  8.7313 0.0001 ***
  #   site                        2  0.74539  0.60150 90.0782 0.0001 ***
  #   sampling_time:sampling_day  3  0.04151  0.03350  3.3442 0.0266 *  
  #   sampling_time:site          2 -0.02104 -0.01698 -2.5424 1.0000    
  # sampling_day:site           6  0.07101  0.05730  2.8604 0.0203 *  
  #   Residual                   53  0.21929  0.17695                   
  # Total                      70  1.23922  1.00000                   
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  #if we set both sampling_time:sampling_day:site = 0, and sampling_time:site = 0:
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time * sampling_day * site - sampling_time:site - sampling_time:sampling_day:site, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs      R2        F Pr(>F)    
  # sampling_time               1  0.07469 0.06027  20.6686 0.0001 ***
  #   sampling_day                3  0.10838 0.08746   9.9971 0.0001 ***
  #   site                        2  0.74539 0.60150 103.1374 0.0001 ***
  #   sampling_time:sampling_day  3  0.04151 0.03350   3.8290 0.0130 *  
  #   sampling_day:site           6  0.07051 0.05690   3.2521 0.0089 ** 
  #   Residual                   55  0.19875 0.16038                    
  # Total                      70  1.23922 1.00000                    
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + sampling_time:sampling_day + sampling_day:site, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs      R2        F Pr(>F)    
  # sampling_time               1  0.07469 0.06027  20.6686 0.0001 ***
  #   sampling_day                3  0.10838 0.08746   9.9971 0.0001 ***
  #   site                        2  0.74539 0.60150 103.1374 0.0001 ***
  #   sampling_time:sampling_day  3  0.04151 0.03350   3.8290 0.0171 *  
  #   sampling_day:site           6  0.07051 0.05690   3.2521 0.0093 ** 
  #   Residual                   55  0.19875 0.16038                    
  # Total                      70  1.23922 1.00000                    
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + sampling_time:sampling_day + sampling_day:site/sampling_time, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs       R2       F Pr(>F)    
  # sampling_time                    1  0.07469  0.06027 15.6865 0.0001 ***
  #   sampling_day                     3  0.10838  0.08746  7.5874 0.0001 ***
  #   site                             2  0.74539  0.60150 78.2766 0.0001 ***
  #   sampling_time:sampling_day       3  0.04151  0.03350  2.9060 0.0458 *  
  #   sampling_day:site                6  0.07051  0.05690  2.4682 0.0404 *  
  #   sampling_time:sampling_day:site  8 -0.02503 -0.02020 -0.6572 1.0000    
  # Residual                        47  0.22378  0.18058                   
  # Total                           70  1.23922  1.00000                   
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + sampling_time/(sampling_day:site), data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs      R2       F Pr(>F)    
  # sampling_time                    1  0.07469 0.06027 15.6865 0.0003 ***
  #   sampling_day                     3  0.10838 0.08746  7.5874 0.0004 ***
  #   site                             2  0.74539 0.60150 78.2766 0.0001 ***
  #   sampling_time:sampling_day:site 17  0.08699 0.07020  1.0747 0.4061    
  # Residual                        47  0.22378 0.18058                   
  # Total                           70  1.23922 1.00000                   
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_time:sampling_day)/site, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs      R2       F Pr(>F)    
  # sampling_time                    1  0.07469 0.06027 15.6865 0.0005 ***
  #   sampling_day                     3  0.10838 0.08746  7.5874 0.0002 ***
  #   site                             2  0.74539 0.60150 78.2766 0.0001 ***
  #   sampling_time:sampling_day       3  0.04151 0.03350  2.9060 0.0451 *  
  #   sampling_time:sampling_day:site 14  0.04548 0.03670  0.6823 0.7729    
  # Residual                        47  0.22378 0.18058                   
  # Total                           70  1.23922 1.00000                   
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  # vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
  # Df SumOfSqs      R2       F Pr(>F)    
  # sampling_time                    1  0.07469 0.06027 15.6865 0.0003 ***
  #   sampling_day                     3  0.10838 0.08746  7.5874 0.0003 ***
  #   site                             2  0.74539 0.60150 78.2766 0.0001 ***
  #   sampling_day:site                6  0.06906 0.05572  2.4173 0.0419 *  
  #   sampling_time:sampling_day:site 11  0.01793 0.01447  0.3424 0.9538    
  # Residual                        47  0.22378 0.18058                   
  # Total                           70  1.23922 1.00000                   
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
}


#instead of changing the distance metric, pool the terms smartly (setting V(sampling_time:sampling_day:site) = 0 and V(sampling_time:site) = 0), by adding them to the next lowest term, in this case it was sampling_time:sampling_day
#

with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
                                 # strata = site,
                                 formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
                                 parallel = nthreads, by = "terms"))

# vegan::adonis2(formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time, data = meta.microb, permutations = 9999, method = "horn", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1  0.07469 0.06027 15.6865 0.0003 ***
#   sampling_day                     3  0.10838 0.08746  7.5874 0.0003 ***
#   site                             2  0.74539 0.60150 78.2766 0.0001 ***
#   sampling_day:site                6  0.06906 0.05572  2.4173 0.0464 *  
#   sampling_time:sampling_day:site 11  0.01793 0.01447  0.3424 0.9479    
# Residual                        47  0.22378 0.18058                   
# Total                           70  1.23922 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Calculate estimates of variation ----------------------------------------

if(any(grepl("sites-", list.files(projectpath, pattern = "usvi_permanova_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_permanova_sites.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
 
permanova_asv_estimate_variation.df <- with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
                                             # strata = site,
                                             formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
                                             parallel = nthreads, by = "terms"))
permanova_asv_estimate_variation.df <- permanova_asv_estimate_variation.df %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::mutate(term = as.character(term)) %>%
  dplyr::mutate(term = dplyr::case_when(term == "sampling_time:sampling_day:site" ~ "pooled",
                                        .default = term)) %>%
  dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                      .default = NA))
permanova_asv_estimate_variation.df <- permanova_asv_estimate_variation.df %>%
  dplyr::mutate(permanova_asv_estimate_variation.df %>%
                     dplyr::filter(!grepl("Total", term)) %>%
                     dplyr::summarise(TotalEMS = sum(MS))) %>%
  dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
                                          .default = NA)) %>%
  dplyr::select(-c(TotalEMS, R2)) %>%
  dplyr::mutate(sq_root = MS^(1/2)) %>%
  dplyr::mutate(perc_variation = sq_root/(sum(sq_root, na.rm = TRUE))*100)
# term Df SumOfSqs       F Pr(>F)      MS V_term sq_root perc_variation
# [1,]    4  1  0.07469 15.6865 0.0002 0.07469 14.896 0.27329         21.177
# [2,]    2  3  0.10838  7.5874 0.0005 0.03613  7.205 0.19007         14.728
# [3,]    6  2  0.74539 78.2766 0.0001 0.37270 74.330 0.61049         47.306
# [4,]    3  6  0.06906  2.4173 0.0419 0.01151  2.295 0.10728          8.313
# [5,]    5 11  0.01793  0.3424 0.9532 0.00163  0.325 0.04038          3.129
# [6,]    1 47  0.22378                0.00476  0.950 0.06900          5.347
# [7,]    7 70  1.23922                                                     

#the column "sq root" typically represents the % variation, however if the sum of the terms' sqrts is not 100, then it's not as accurate
#normalize to 100% by summing the sqrts of the terms' estimates

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
  bind_rows(., .id = "site") %>%
  dplyr::group_by(site)

  
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
#   1 LB_seagrass   0.127   0.00266    0.424             0.299           0.00628             1    
# 2 Tektite       0.0807  0.00273    0.423             0.190           0.00645             0.997
# 3 Yawzi         0.0634  0.00258    0.376             0.149           0.00607             0.886

# 
# dist_usvi_asv.df %>%
#   dplyr::group_by(site) %>%
#   TukeyHSD(.)

# #is it different with log-transformed?
# 
# meta.microb %>%
#   split(., f = .$site) %>%
#   map(., ~.x %>%
#         tibble::rownames_to_column(var = "sample_id") %>%
#         dplyr::select(sample_id) %>%
#         tibble::deframe(.)) %>%
#   map(., ~dist_usvi_asv_log2.mat[.x, .x])  %>%
#   map(., ~.x[lower.tri(.x, diag = FALSE)]) %>%
#   map(., ~.x %>% tibble::enframe(., name = NULL, value = "dissimilarity")) %>%
#   bind_rows(., .id = "site") %>%
#   dplyr::group_by(site) %>%
#   dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# # site        avg_dist
# # <chr>          <dbl>
# #   1 LB_seagrass   0.119 
# # 2 Tektite       0.0600
# # 3 Yawzi         0.0390

#what about if we used bray-curtis instead of MH?
#still LB is higher than Tektite and higher than Yawzi


#approach 1: Kruskal-Wallis testing
{
  # #Tektite vs Yawzi:
  # dist_usvi_asv.df %>%
  #   dplyr::filter(!grepl("LB", site)) %>%
  #   kruskal.test(dissimilarity ~ site, .)
  # # Kruskal-Wallis rank sum test
  # # 
  # # data:  dissimilarity by site
  # # Kruskal-Wallis chi-squared = 10.624, df = 1, p-value = 0.001116
  # 
  # #Yawzi vs LB:
  # dist_usvi_asv.df %>%
  #   dplyr::filter(!grepl("Tektite", site)) %>%
  #   kruskal.test(dissimilarity ~ site, .)
  # # Kruskal-Wallis rank sum test
  # # 
  # # data:  dissimilarity by site
  # # Kruskal-Wallis chi-squared = 94.392, df = 1, p-value < 2.2e-16
  # 
  # #Tektite vs LB:
  # dist_usvi_asv.df %>%
  #   dplyr::filter(!grepl("Yawzi", site)) %>%
  #   kruskal.test(dissimilarity ~ site, .)
  # # Kruskal-Wallis rank sum test
  # # 
  # # data:  dissimilarity by site
  # # Kruskal-Wallis chi-squared = 48.111, df = 1, p-value = 4.028e-12
  # 
  # #all 3:
  # dist_usvi_asv.df %>%
  #   kruskal.test(dissimilarity ~ site, .)
  # # Kruskal-Wallis rank sum test
  # # 
  # # data:  dissimilarity by site
  # # Kruskal-Wallis chi-squared = 100.96, df = 2, p-value < 2.2e-16
  # 
  # 
  # #which is identical to:
  # 
  # # temp_list <- meta.microb %>%
  # #   split(., f = .$site) %>%
  # #   map(., ~.x %>%
  # #         tibble::rownames_to_column(var = "sample_id") %>%
  # #         dplyr::select(sample_id) %>%
  # #         tibble::deframe(.)) %>%
  # #   map(., ~dist_usvi_asv.mat[.x, .x])  %>%
  # #   map(., ~.x[lower.tri(.x, diag = FALSE)])
  # # 
  # # kruskal.test(temp_list)
  # #Kruskal-Wallis chi-squared = 167.83, df = 2, p-value < 2.2e-16
  # # 
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
                     dplyr::select(metab_deriv_label, site)) %>%
    tidyr::drop_na(.)
temp_kw_res <- temp_list %>%
  split(., f = .$metabolite) %>%
  map(., ~.x %>%
        split(., f = .$site) %>%
        map(., ~.x %>%
                    dplyr::select(grouping, concentration) %>%
                    tidyr::drop_na(.) %>%
              droplevels) %>%
        map_if(., ~(length(unique(.x[["grouping"]])) > 1), ~.x %>% kruskal.test(concentration ~ grouping, .), .else = NULL)) %>%
  map_depth(., 2, ~.x %>%
              purrr::pluck("p.value")) %>%
  bind_rows(., .id = "metabolite") %>%
  tidyr::pivot_longer(., cols = !"metabolite",
                      names_to = "site.dawn_vs_afternoon",
                      values_to = "p.value") %>%
  tidyr::drop_na(.)
  

# temp_kw_res <- temp_list %>%
#   split(., f = .$metabolite) %>%
#   map(., ~.x %>%
#         split(., f = .$site) %>%
#         map(., ~.x %>% 
#               kruskal.test(concentration ~ grouping, .))) %>%
#   map_depth(., 2, ~.x %>%
#               purrr::pluck("p.value")) %>%
#   bind_rows(., .id = "metabolite") %>%
#   tidyr::pivot_longer(., cols = !"metabolite",
#                       names_to = "site.dawn_vs_afternoon",
#                       values_to = "p.value")


q_cutoff <- p.adjust(temp_kw_res[["p.value"]], method = "fdr") %>%
  na.omit(.) %>%
  ashr::qval.from.lfdr(.) %>%
  quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
temp_kw_res <- temp_kw_res %>%
  dplyr::mutate(filtered_p.value = dplyr::case_when(p.value > 0.05 ~ NA,
                                               .default = p.value)) %>%
  dplyr::mutate(padj_10 = dplyr::case_when(p.value <= q_cutoff[2] ~ p.value,
                                           .default = NA)) %>%
  dplyr::mutate(padj_05 = dplyr::case_when(p.value <= q_cutoff[1] ~ p.value,
                                           .default = NA))
length(na.omit(temp_kw_res[["filtered_p.value"]])) #26
length(na.omit(temp_kw_res[["padj_10"]])) #20
length(na.omit(temp_kw_res[["padj_05"]])) #15

  
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
  
  
  # LB_seagrass     Tektite
  # Tektite 2.906243e-10          NA
  # Yawzi   7.636926e-19 0.003527418
  
  
  t.test(dist_usvi_asv.df %>%
           dplyr::filter(grepl("Tektite", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         dist_usvi_asv.df %>%
           dplyr::filter(grepl("Yawzi", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         conf.level = 0.99)
  
  # t = 2.9305, df = 539.46, p-value = 0.003527
  # alternative hypothesis: true difference in means is not equal to 0
  # 99 percent confidence interval:
  #   0.002043565 0.032619614
  # sample estimates:
  #   mean of x  mean of y 
  # 0.08070424 0.06337265 
  
  t.test(dist_usvi_asv.df %>%
           dplyr::filter(grepl("LB", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         dist_usvi_asv.df %>%
           dplyr::filter(grepl("Tektite", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.))
  # t = 6.5023, df = 494.4, p-value = 1.937e-10
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.03229326 0.06025932
  # sample estimates:
  #   mean of x  mean of y 
  # 0.12698053 0.08070424 
  
  t.test(dist_usvi_asv.df %>%
           dplyr::filter(grepl("LB", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.),
         dist_usvi_asv.df %>%
           dplyr::filter(grepl("Yawzi", site)) %>%
           dplyr::select(dissimilarity) %>%
           tibble::deframe(.))
  # t = 9.4029, df = 458.51, p-value < 2.2e-16
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.05031426 0.07690149
  # sample estimates:
  #   mean of x  mean of y 
  # 0.12698053 0.06337265 
  
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


#testing it:

set.seed(48105)
# dist_usvi_asv_sites.boot <- boot::boot(dist_usvi_asv.df, F_ttest_gn, variable = "site", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)

resample_depth <- dist_usvi_asv.df %>%
    split(., f = .$site) %>%
    map(., ~round(nrow(.x)/10)) %>%
    purrr::reduce(mean)
dist_usvi_asv_sites.boot <- boot::boot(dist_usvi_asv.df, F_ttest_para, variable = "site", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)



# (dist_usvi_asv_sites.boot$t %>%
#     as.data.frame() %>%
#     dplyr::reframe(across(everything(), \(x) quantile(x, probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE, type = 7))))
#     
#     # dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>% t())

#for reference: 
# q_value <- 0.10
dist_usvi_asv_sites.boot.summary.df <- unique(dist_usvi_asv.df[["site"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_sites.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value appears to be random based on repeating this bootstrpaping procedure.
  # bind_cols(., tibble::enframe(dist_usvi_asv_sites.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q05" = (dist_usvi_asv_sites.boot$t %>%
                                 as.data.frame() %>%
                                 dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                 t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_asv_sites.boot$t %>%
                                         as.data.frame() %>%
                                         dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                         t())))
# hist(dist_usvi_asv_sites.boot$t[,1])

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


t.test(dist_usvi_asv.df %>%
         dplyr::filter(grepl("Yawzi", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_asv.df %>%
         dplyr::filter(grepl("Tektite", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       conf.level = 0.95)


#reepat for the metabolomes:

#SS(total) for metabolomes:  2.554602
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
#LB: 0.2570766
#Tektite: 0.9269946
#Yawzi: 0.6523068

with(meta.metab, vegan::adonis2(data = meta.metab, method = "bray", permutations = 9999,
                                # strata = site,
                                formula = dist_usvi_metab.d ~ sampling_time*sampling_day*site,
                                parallel = nthreads, by = "terms"))
# Df SumOfSqs      R2       F Pr(>F)
# sampling_time                    1  0.09094 0.03560  3.9751 0.0174 *
#   sampling_day                     3  0.08985 0.03517  1.3091 0.2456
# site                             2  0.71497 0.27987 15.6259 0.0001 ***
#   sampling_time:sampling_day       3  0.08584 0.03360  1.2507 0.2708
# sampling_time:site               2  0.05069 0.01984  1.1079 0.3426
# sampling_day:site                6  0.29794 0.11663  2.1705 0.0119 *
#   sampling_time:sampling_day:site  6  0.10337 0.04046  0.7530 0.7189
# Residual                        49  1.12100 0.43882
# Total                           72  2.55460 1.00000
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


permanova_metab_estimate_variation.df <- with(meta.metab, vegan::adonis2(data = meta.metab, method = "bray", permutations = 9999,
                                            # strata = site,
                                            formula = dist_usvi_metab.d ~ sampling_time*sampling_day*site,
                                            parallel = nthreads, by = "terms"))
permanova_metab_estimate_variation.df <- permanova_metab_estimate_variation.df %>%
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
  dplyr::select(-c(TotalEMS, R2)) %>%
  dplyr::mutate(sq_root = MS^(1/2)) %>%
  dplyr::mutate(perc_variation = sq_root/(sum(sq_root, na.rm = TRUE))*100)


# term Df SumOfSqs       F Pr(>F)      MS V_term sq_root perc_variation
# [1,]    4  1  0.09094  3.9751 0.0155 0.09094 14.618 0.30156        15.8199
# [2,]    2  3  0.08985  1.3091 0.2464 0.02995  4.814 0.17306         9.0786
# [3,]    8  2  0.71497 15.6259 0.0001 0.35748 57.464 0.59790        31.3656
# [4,]    5  3  0.08584  1.2507 0.2704 0.02861  4.599 0.16915         8.8737
# [5,]    7  2  0.05069  1.1079 0.3532 0.02535  4.074 0.15920         8.3518
# [6,]    3  6  0.29794  2.1705 0.0113 0.04966  7.982 0.22284        11.6900
# [7,]    6  6  0.10337  0.7530 0.7280 0.01723  2.769 0.13126         6.8856
# [8,]    1 49  1.12100                0.02288  3.678 0.15125         7.9347
# [9,]    9 72  2.55460                   

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
  bind_rows(., .id = "site") %>%
  dplyr::group_by(site)
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
#   1 LB_seagrass    0.143   0.0535    0.325             0.236            0.0883             0.537
# 2 Tektite        0.251   0.0294    0.605             0.414            0.0485             1    
# 3 Yawzi          0.203   0.0350    0.529             0.335            0.0577             0.873

# dist_usvi_metab.df %>%
#   dplyr::group_by(site) %>%
#   dplyr::summarise(min_dist = min(dissimilarity, na.rm = TRUE),
#                    max_dist = max(dissimilarity, na.rm = TRUE)) %>%
#   tidyr::pivot_longer(., cols = contains("dist"),
#                       names_to = "metric",
#                       values_to = "dist") %>%
#   dplyr::mutate(dist = dplyr::case_when(grepl("min", metric) ~ -(dist),
#                                             .default = dist)) %>%
#   dplyr::group_by(site) %>%
#   dplyr::summarise(delta_beta = sum(dist, na.rm = TRUE)) %>%
#   droplevels

# #Tektite vs Yawzi:
# dist_usvi_metab.df %>%
#   dplyr::filter(!grepl("LB", site)) %>%
#   kruskal.test(dissimilarity ~ site, .)
# 
# # data:  dissimilarity by site
# # Kruskal-Wallis chi-squared = 24.227, df = 1, p-value = 8.563e-07
# 
# #Yawzi vs LB:
# dist_usvi_metab.df %>%
#   dplyr::filter(!grepl("Tektite", site)) %>%
#   kruskal.test(dissimilarity ~ site, .)
# # Kruskal-Wallis chi-squared = 24.928, df = 1, p-value = 5.952e-07
# 
# #Tektite vs LB:
# dist_usvi_metab.df %>%
#   dplyr::filter(!grepl("Yawzi", site)) %>%
#   kruskal.test(dissimilarity ~ site, .)
# # Kruskal-Wallis chi-squared = 123.09, df = 1, p-value < 2.2e-16
# 
# #all 3:
# dist_usvi_metab.df %>%
#   kruskal.test(dissimilarity ~ site, .)
# # data:  dissimilarity by site
# # Kruskal-Wallis chi-squared = 112.63, df = 2, p-value < 2.2e-16
# 
# 
# pairwise.t.test((dist_usvi_metab.df %>%
#                    split(., f = .$site) %>%
#                    map(., ~.x %>%
#                          dplyr::select(dissimilarity) %>%
#                          tibble::deframe(.)) %>%
#                    purrr::reduce(c)),
#                 (dist_usvi_metab.df %>%
#                    split(., f = .$site) %>%
#                    map(., ~.x %>%
#                          dplyr::select(site) %>%
#                          tibble::deframe(.)) %>%
#                    purrr::reduce(c)),
#                 pool.sd = FALSE, paired = FALSE, alternative = "two.sided",
#                 p.adjust.method = "BH")
# 
# # LB_seagrass Tektite
# # Tektite < 2e-16     -      
# #   Yawzi   2.0e-14     8.6e-07
# # 
# # P value adjustment method: BH 


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

set.seed(48105)
#try with strata?
# temp_df <- dist_usvi_metab.df %>%
#   # dplyr::mutate(series = as.character(site)) %>%
#   dplyr::mutate(series = dplyr::case_when(site == "LB_seagrass" ~ 1,
#                                           site == "Yawzi" ~ 2,
#                                           site == "Tektite" ~ 3)) %>%
# dplyr::group_by(site)
# 
# temp.boot <- boot::boot(temp_df, F_ttest_gn, variable = "site", 
#                         strata = temp_df[["series"]], 
#                         sim = "permutation",
#                         R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
# # hist(temp.boot$t[,1], breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1))
# # hist(temp.boot$t[,2], breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1))
# # hist(temp.boot$t[,3], breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1))
# 
# temp.boot.summary.df <- unique(dist_usvi_metab.df[["site"]]) %>%
#   combn(., 2) %>%
#   t() %>%
#   as.data.frame() %>%
#   setNames(., c("pair1", "pair2")) %>%
#   bind_cols(., tibble::enframe(temp.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value appears to be random based on repeating this bootstrpaping procedure.
#   # bind_cols(., tibble::enframe(dist_usvi_metab_sites.boot$t0, name = NULL, value = "p.value")) %>%
#   # bind_cols(., data.frame("q.value" = (dist_usvi_metab_sites.boot$t %>%
#   #                                        setNames(., c("V1", "V2", "V3")) %>%
#   #                                        as.data.frame() %>%
#   #                                        dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
#   #                                        t())))
#   bind_cols(., data.frame("q05" = (temp.boot$t %>%
#                                      as.data.frame() %>%
#                                      dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
#                                      t()))) %>%
#   bind_cols(., data.frame("q10" = (temp.boot$t %>%
#                                      as.data.frame() %>%
#                                      dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
#                                      t())))

# 
# temp_df %>%
#   # dplyr::group_by({variable2}) %>%
#   dplyr::slice_sample(., n = 25) %>%
#   droplevels

# #know that mle = resample_depth
# temp.boot <- boot::boot(temp_df, F_ttest_para, variable = "site", 
#                        R = 999, sim = "parametric", 
#                        # strata = temp_df[["series"]], 
#                        ran.gen = F_ttest_rn, mle = 25)


# dist_usvi_metab_sites.boot <- boot::boot(dist_usvi_metab.df, F_ttest_gn, variable = "site", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
resample_depth <- dist_usvi_metab.df %>%
  split(., f = .$site) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_metab_sites.boot <- boot::boot(dist_usvi_metab.df, F_ttest_para, variable = "site", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)



#for reference: 
q_value <- 0.10
dist_usvi_metab_sites.boot.summary.df <- unique(dist_usvi_metab.df[["site"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_sites.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value appears to be random based on repeating this bootstrpaping procedure.
  # bind_cols(., tibble::enframe(dist_usvi_metab_sites.boot$t0, name = NULL, value = "p.value")) %>%
  # bind_cols(., data.frame("q.value" = (dist_usvi_metab_sites.boot$t %>%
  #                                        setNames(., c("V1", "V2", "V3")) %>%
  #                                        as.data.frame() %>%
  #                                        dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
  #                                        t())))
bind_cols(., data.frame("q05" = (dist_usvi_metab_sites.boot$t %>%
                                   as.data.frame() %>%
                                   dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                   t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_metab_sites.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t())))

# hist(dist_usvi_metab_sites.boot$t[,1], breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1))
# hist(dist_usvi_metab_sites.boot$t[,2], breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1))
# hist(dist_usvi_metab_sites.boot$t[,3], breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1))


t.test(dist_usvi_metab.df %>%
         dplyr::filter(grepl("Tektite", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab.df %>%
         dplyr::filter(grepl("Yawzi", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       conf.level = 0.95)
# t = 4.9744, df = 596.93, p-value = 8.569e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02899711 0.06683073
# sample estimates:
#   mean of x mean of y 
# 0.250585  0.202671 

t.test(dist_usvi_metab.df %>%
         dplyr::filter(grepl("LB", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab.df %>%
         dplyr::filter(grepl("Tektite", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.))
# t = -13.896, df = 430.39, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.12281951 -0.09238019
# sample estimates:
#   mean of x mean of y 
# 0.1429851 0.2505850 

t.test(dist_usvi_metab.df %>%
         dplyr::filter(grepl("LB", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab.df %>%
         dplyr::filter(grepl("Yawzi", site)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.))
# t = -7.9733, df = 440.05, p-value = 1.345e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.07439820 -0.04497366
# sample estimates:
#   mean of x mean of y 
# 0.1429851 0.2026710 


#Summarizing intra-site distances via pariwise t-tests:

length_v <- unique(dist_usvi_metab.df[["site"]]) %>% as.character(.)
# length_v <- c("LB_seagrass", "Tektite", "Yawzi")
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA)  %>%
  dplyr::mutate(p.value = NA)
ttest_usvi_metab_sites.res <- ttest_res
ttest_usvi_asv_sites.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_metab_sites.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic
  ttest_usvi_metab_sites.res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
  
  temp_i_c <- dist_usvi_asv.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_asv_sites.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  ttest_usvi_asv_sites.res[i, "p.value"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$p.value
}

ttest_usvi_metab_sites.res <- ttest_usvi_metab_sites.res %>%
  dplyr::left_join(., dist_usvi_metab_sites.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))
ttest_usvi_asv_sites.res <- ttest_usvi_asv_sites.res %>%
  dplyr::left_join(., dist_usvi_asv_sites.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))

readr::write_delim(ttest_usvi_metab_sites.res, paste0(projectpath, "/", "permanova_ttest_usvi_metab_sites.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)
readr::write_delim(ttest_usvi_asv_sites.res, paste0(projectpath, "/", "permanova_ttest_usvi_asv_sites.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)

save(dist_usvi_metab_sites.boot, dist_usvi_metab_sites.boot.summary.df, 
     permanova_metab_estimate_variation.df, meta.metab, dist_usvi_metab.df,
     dist_usvi_asv_sites.boot, dist_usvi_asv_sites.boot.summary.df, 
     permanova_asv_estimate_variation.df, meta.microb, dist_usvi_asv.df,
     ttest_usvi_metab_sites.res, ttest_usvi_asv_sites.res,
     file = paste0(projectpath, "/", "usvi_permanova_sites-", Sys.Date(), ".RData"))

}


# Look at within-site within-sampling time --------------------------------


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
# 2 peak_photo       0.145

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
#   1 dawn             0.231
# 2 peak_photo       0.250


if(any(grepl("time-", list.files(projectpath, pattern = "usvi_permanova_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_permanova_time.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
  
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
  bind_rows(., .id = "sampling_time") %>%
  dplyr::group_by(sampling_time)

t.test(dist_usvi_metab_time.df %>%
         dplyr::filter(grepl("dawn", sampling_time)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       dist_usvi_metab_time.df %>%
         dplyr::filter(grepl("peak", sampling_time)) %>%
         dplyr::select(dissimilarity) %>%
         tibble::deframe(.),
       conf.level = 0.95)

# dist_usvi_metab_time.boot <- boot::boot(dist_usvi_metab_time.df, F_ttest_gn, variable = "sampling_time", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
resample_depth <- dist_usvi_metab_time.df %>%
  split(., f = .$sampling_time) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_metab_time.boot <- boot::boot(dist_usvi_metab_time.df, F_ttest_para, variable = "sampling_time", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)


hist(dist_usvi_metab_time.boot$t[,1])
# q_value <- 0.10
dist_usvi_metab_time.boot.summary.df <- unique(dist_usvi_metab_time.df[["sampling_time"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  # bind_cols(., tibble::enframe(dist_usvi_metab_time.boot$t0, name = NULL, value = "p.value")) %>%
  # bind_cols(., data.frame("q.value" = (dist_usvi_metab_time.boot$t %>%
  #                                        as.data.frame() %>%
  #                                        dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
  #                                        t())))
  bind_cols(., tibble::enframe(dist_usvi_metab_time.boot$t0, name = NULL, value = "boot.p.value")) %>%
  bind_cols(., data.frame("q05" = (dist_usvi_metab_time.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_metab_time.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
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
  bind_rows(., .id = "sampling_time") %>%
  dplyr::group_by(sampling_time)


# dist_usvi_asv_time.boot <- boot::boot(dist_usvi_asv_time.df, F_ttest_gn, variable = "sampling_time", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
resample_depth <- dist_usvi_asv_time.df %>%
  split(., f = .$sampling_time) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_asv_time.boot <- boot::boot(dist_usvi_asv_time.df, F_ttest_para, variable = "sampling_time", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)


hist(dist_usvi_asv_time.boot$t[,1])

# q_value <- 0.10
dist_usvi_asv_time.boot.summary.df <- unique(dist_usvi_asv_time.df[["sampling_time"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_time.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value appears to be random based on repeating this bootstrpaping procedure.
  # bind_cols(., tibble::enframe(dist_usvi_asv_time.boot$t0, name = NULL, value = "p.value")) %>%
  # bind_cols(., data.frame("q.value" = (dist_usvi_asv_time.boot$t %>%
  #                                        as.data.frame() %>%
  #                                        dplyr::reframe(across(everything(), quantile, probs = q_value, na.rm = TRUE, names = FALSE, type = 7)) %>%
  #                                        t())))
  bind_cols(., data.frame("q05" = (dist_usvi_asv_time.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_asv_time.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t())))

length_v <- unique(dist_usvi_metab_time.df[["sampling_time"]]) %>% as.character(.)
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA) %>%
  dplyr::mutate(`p.value` = NA)
ttest_usvi_metab_time.res <- ttest_res
ttest_usvi_asv_time.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab_time.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab_time.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_metab_time.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic
  ttest_usvi_metab_time.res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
  
  temp_i_c <- dist_usvi_asv_time.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv_time.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_asv_time.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  ttest_usvi_asv_time.res[i, "p.value"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$p.value
}

ttest_usvi_metab_time.res <- ttest_usvi_metab_time.res %>%
  dplyr::left_join(., dist_usvi_metab_time.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))
ttest_usvi_asv_time.res <- ttest_usvi_asv_time.res %>%
  dplyr::left_join(., dist_usvi_asv_time.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))

ttest_usvi_metab_time.res

ttest_usvi_asv_time.res

readr::write_delim(ttest_usvi_metab_time.res, paste0(projectpath, "/", "permanova_ttest_usvi_metab_time.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)
readr::write_delim(ttest_usvi_asv_time.res, paste0(projectpath, "/", "permanova_ttest_usvi_asv_time.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)



save(dist_usvi_metab_time.boot, dist_usvi_metab_time.boot.summary.df, dist_usvi_metab_time.df,
     dist_usvi_asv_time.boot, dist_usvi_asv_time.boot.summary.df, dist_usvi_asv_time.df,
     ttest_usvi_metab_time.res, ttest_usvi_asv_time.res,
     file = paste0(projectpath, "/", "usvi_permanova_time-", Sys.Date(), ".RData"))

}
# Site by time ------------------------------------------------------------

if(any(grepl("site_time-", list.files(projectpath, pattern = "usvi_permanova_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_permanova_site_time.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
  
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
  bind_rows(., .id = "grouping") %>%
  dplyr::group_by(grouping)

dist_usvi_asv_grouping.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
                   min_dist = min(dissimilarity, na.rm = TRUE),
                   max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_asv_grouping.df[["dissimilarity"]]))))
# grouping               avg_dist min_dist max_dist avg_dist_rescaled min_dist_rescaled max_dist_rescaled
# <chr>                     <dbl>    <dbl>    <dbl>             <dbl>             <dbl>             <dbl>
#   1 LB_seagrass.dawn         0.0916  0.00266    0.270             0.229           0.00665             0.673
# 2 LB_seagrass.peak_photo   0.156   0.00873    0.401             0.390           0.0218              1    
# 3 Tektite.dawn             0.0750  0.00273    0.291             0.187           0.00682             0.725
# 4 Tektite.peak_photo       0.0899  0.00377    0.329             0.224           0.00941             0.821
# 5 Yawzi.dawn               0.0749  0.00259    0.278             0.187           0.00646             0.693
# 6 Yawzi.peak_photo         0.0444  0.00500    0.198             0.111           0.0125              0.493

dist_usvi_asv_grouping.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                time = stringr::str_split_i(grouping, "\\.", 2)) %>%
  dplyr::mutate(across(c(site, time), ~factor(.x))) %>%
  dplyr::mutate(avg_dist = dplyr::case_when(time == "dawn" ~ -(avg_dist),
                                            .default = avg_dist)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(delta_beta = sum(avg_dist, na.rm = TRUE)) %>%
  droplevels
# # A tibble: 3 × 2
# site        delta_beta
# <fct>            <dbl>
#   1 LB_seagrass     0.0645
# 2 Tektite         0.0149
# 3 Yawzi          -0.0305

set.seed(48105)
# dist_usvi_asv_grouping.boot <- boot::boot(dist_usvi_asv_grouping.df, F_ttest_gn, variable = "grouping", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
resample_depth <- dist_usvi_asv_grouping.df %>%
  split(., f = .$grouping) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_asv_grouping.boot <- boot::boot(dist_usvi_asv_grouping.df, F_ttest_para, variable = "grouping", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)


# q_value <- 0.10
dist_usvi_asv_grouping.boot.summary.df <- unique(dist_usvi_asv_grouping.df[["grouping"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  # bind_cols(., tibble::enframe(dist_usvi_asv_grouping.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_grouping.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value appears to be random based on repeating this bootstrpaping procedure.
  bind_cols(., data.frame("q05" = (dist_usvi_asv_grouping.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_asv_grouping.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
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
  bind_rows(., .id = "grouping") %>%
  dplyr::group_by(grouping)


dist_usvi_metab_grouping.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
                   min_dist = min(dissimilarity, na.rm = TRUE),
                   max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_metab_grouping.df[["dissimilarity"]]))))
# grouping               avg_dist min_dist max_dist avg_dist_rescaled min_dist_rescaled max_dist_rescaled
# <chr>                     <dbl>    <dbl>    <dbl>             <dbl>             <dbl>             <dbl>
#   1 LB_seagrass.dawn          0.128   0.0535    0.259             0.220            0.0919             0.446
# 2 LB_seagrass.peak_photo    0.133   0.0556    0.284             0.229            0.0956             0.488
# 3 Tektite.dawn              0.260   0.0414    0.582             0.447            0.0711             1    
# 4 Tektite.peak_photo        0.232   0.0294    0.498             0.399            0.0505             0.857
# 5 Yawzi.dawn                0.134   0.0350    0.332             0.230            0.0601             0.571
# 6 Yawzi.peak_photo          0.259   0.0686    0.529             0.445            0.118              0.909

dist_usvi_metab_grouping.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                time = stringr::str_split_i(grouping, "\\.", 2)) %>%
  dplyr::mutate(across(c(site, time), ~factor(.x))) %>%
  dplyr::mutate(avg_dist = dplyr::case_when(time == "dawn" ~ -(avg_dist),
                                            .default = avg_dist)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(delta_beta = sum(avg_dist, na.rm = TRUE)) %>%
  droplevels
# # A tibble: 3 × 2
# site        delta_beta
# <fct>            <dbl>
#   1 LB_seagrass    0.00551
# 2 Tektite       -0.0277 
# 3 Yawzi          0.125  

set.seed(48105)
resample_depth <- dist_usvi_metab_grouping.df %>%
  split(., f = .$grouping) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_metab_grouping.boot <- boot::boot(dist_usvi_metab_grouping.df, F_ttest_para, variable = "grouping", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_metab_grouping.boot.summary.df <- unique(dist_usvi_metab_grouping.df[["grouping"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_grouping.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value appears to be random based on repeating this bootstrpaping procedure.
  # bind_cols(., tibble::enframe(dist_usvi_metab_grouping.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q05" = (dist_usvi_metab_grouping.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_metab_grouping.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t())))

length_v <- unique(dist_usvi_metab_grouping.df[["grouping"]]) %>% as.character(.)
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA) %>%
  dplyr::mutate(p.value = NA)
ttest_usvi_metab_grouping.res <- ttest_res
ttest_usvi_asv_grouping.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab_grouping.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab_grouping.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_metab_grouping.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic
  ttest_usvi_metab_grouping.res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value
  
  temp_i_c <- dist_usvi_asv_grouping.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv_grouping.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  ttest_usvi_asv_grouping.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
  ttest_usvi_asv_grouping.res[i, "p.value"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$p.value
  
}

ttest_usvi_metab_grouping.res <- ttest_usvi_metab_grouping.res %>%
  # dplyr::left_join(., dist_usvi_metab_grouping.boot.summary.df)
  dplyr::left_join(., dist_usvi_metab_grouping.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))
ttest_usvi_asv_grouping.res <- ttest_usvi_asv_grouping.res %>%
  # dplyr::left_join(., dist_usvi_asv_grouping.boot.summary.df)
  dplyr::left_join(., dist_usvi_asv_grouping.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))

ttest_usvi_metab_grouping.res
# pair1                  pair2           t      p.value         q05         q10
# 1        LB_seagrass.dawn           Tektite.dawn -7.18817569 5.879049e-02 0.002838329 0.006115816
# 2        LB_seagrass.dawn             Yawzi.dawn -0.47572203 1.629540e-02 0.059171237 0.114628017
# 3        LB_seagrass.dawn LB_seagrass.peak_photo -0.53043075 6.372982e-01 0.070538418 0.121354236
# 4        LB_seagrass.dawn     Tektite.peak_photo -7.07125354 3.284103e-07 0.002657340 0.006489542
# 5        LB_seagrass.dawn       Yawzi.peak_photo -7.61852425 1.505827e-02 0.001502684 0.003824924
# 6            Tektite.dawn             Yawzi.dawn  6.39916002 1.423647e-01 0.003158393 0.006926526
# 7            Tektite.dawn LB_seagrass.peak_photo  6.81924084 7.541640e-01 0.003044016 0.007031463
# 8            Tektite.dawn     Tektite.peak_photo  1.30189416 9.267313e-01 0.053462899 0.101644100
# 9            Tektite.dawn       Yawzi.peak_photo  0.04674720 8.979607e-01 0.066708965 0.126135985
# 10             Yawzi.dawn LB_seagrass.peak_photo  0.02861654 5.370420e-01 0.042630277 0.097260366
# 11             Yawzi.dawn     Tektite.peak_photo -6.00281601 5.546100e-01 0.003977986 0.008832722
# 12             Yawzi.dawn       Yawzi.peak_photo -6.71707758 2.745982e-01 0.002222292 0.005185253
# 13 LB_seagrass.peak_photo     Tektite.peak_photo -6.59460114 6.061675e-02 0.003388207 0.008130165
# 14 LB_seagrass.peak_photo       Yawzi.peak_photo -7.21474369 4.199997e-01 0.002016677 0.004734910
# 15     Tektite.peak_photo       Yawzi.peak_photo -1.31343777 9.405063e-01 0.050585101 0.099311294
ttest_usvi_asv_grouping.res
# pair1                  pair2            t    p.value          q05         q10
# 1        LB_seagrass.dawn           Tektite.dawn  1.410729901 0.30252946 0.0391881551 0.079333501
# 2        LB_seagrass.dawn             Yawzi.dawn  1.391956585 0.94932343 0.0375353927 0.083367831
# 3        LB_seagrass.dawn LB_seagrass.peak_photo -4.270729785 0.00933399 0.0074522812 0.017894323
# 4        LB_seagrass.dawn     Tektite.peak_photo  0.137144000 0.57207990 0.0656889120 0.123106245
# 5        LB_seagrass.dawn       Yawzi.peak_photo  4.945466223 0.05298388 0.0081880989 0.017390905
# 6            Tektite.dawn             Yawzi.dawn  0.005695843 0.22500250 0.0819694534 0.132361170
# 7            Tektite.dawn LB_seagrass.peak_photo -5.354195172 0.16193637 0.0029496823 0.007127690
# 8            Tektite.dawn     Tektite.peak_photo -1.166492963 0.84264153 0.0657565985 0.109702051
# 9            Tektite.dawn       Yawzi.peak_photo  3.175916699 0.44402908 0.0388328193 0.062521063
# 10             Yawzi.dawn LB_seagrass.peak_photo -5.301768123 0.12500506 0.0027059068 0.007034590
# 11             Yawzi.dawn     Tektite.peak_photo -1.154447771 0.16792411 0.0631370396 0.106698918
# 12             Yawzi.dawn       Yawzi.peak_photo  3.087296629 0.85116696 0.0404019377 0.067588469
# 13 LB_seagrass.peak_photo     Tektite.peak_photo  4.172437573 0.05602462 0.0060249461 0.013647421
# 14 LB_seagrass.peak_photo       Yawzi.peak_photo  8.288106126 0.02527052 0.0007432683 0.001861703
# 15     Tektite.peak_photo       Yawzi.peak_photo  4.238559560 0.36369284 0.0207406480 0.036856611

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
readr::write_delim(ttest_usvi_metab_grouping.res, paste0(projectpath, "/", "permanova_ttest_usvi_metab_site_time.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)
readr::write_delim(ttest_usvi_asv_grouping.res, paste0(projectpath, "/", "permanova_ttest_usvi_asv_site_time.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)

save(dist_usvi_metab_grouping.boot, dist_usvi_metab_grouping.boot.summary.df, dist_usvi_metab_grouping.df,
     dist_usvi_asv_grouping.boot, dist_usvi_asv_grouping.boot.summary.df, dist_usvi_asv_grouping.df,
     ttest_usvi_metab_grouping.res, ttest_usvi_asv_grouping.res,
     file = paste0(projectpath, "/", "usvi_permanova_site_time-", Sys.Date(), ".RData"))

}

# Pairwise t-tests for day and site ---------------------------------------


if(any(grepl("site_day-", list.files(projectpath, pattern = "usvi_permanova_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_permanova_site_day.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
  
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
  bind_rows(., .id = "sampling_day") %>%
  dplyr::group_by(sampling_day)
dist_usvi_day.df %>% 
  dplyr::group_by(sampling_day) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))
# sampling_day avg_dist
# <chr>           <dbl>
#   1 Day2            0.176
# 2 Day3            0.140
# 3 Day4            0.128
# 4 Day5            0.135

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
    dplyr::mutate(type = "metab"))) %>%
  dplyr::group_by(sampling_day)


dist_usvi_day.df %>%
  dplyr::group_by(type, sampling_day) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE))  %>%
  dplyr::left_join(., dist_usvi_day.df %>%
                     dplyr::group_by(type) %>%
                     dplyr::summarise(max_dist = max(dissimilarity, na.rm = TRUE))) %>%
  dplyr::mutate(avg_dist_rescale = avg_dist/max_dist)
# type  sampling_day avg_dist max_dist avg_dist_rescale
# <chr> <chr>           <dbl>    <dbl>            <dbl>
#   1 asv   Day2            0.176    0.604            0.291
# 2 asv   Day3            0.140    0.604            0.232
# 3 asv   Day4            0.128    0.604            0.212
# 4 asv   Day5            0.135    0.604            0.223
# 5 metab Day2            0.247    0.605            0.408
# 6 metab Day3            0.246    0.605            0.407
# 7 metab Day4            0.242    0.605            0.401
# 8 metab Day5            0.247    0.605            0.408


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
  bind_rows(., .id = "site_day") %>%
  dplyr::group_by(site_day)

dist_usvi_asv_site_day.df %>%
  dplyr::group_by(site_day) %>%
  # dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
  #                  min_dist = min(dissimilarity, na.rm = TRUE),
  #                  max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_asv_site_day.df[["dissimilarity"]]))))
# site_day         avg_dist avg_dist_rescaled
# <chr>               <dbl>             <dbl>
#   1 LB_seagrass.Day2   0.0884             0.224
# 2 LB_seagrass.Day3   0.110              0.279
# 3 LB_seagrass.Day4   0.0924             0.235
# 4 LB_seagrass.Day5   0.119              0.303
# 5 Tektite.Day2       0.102              0.260
# 6 Tektite.Day3       0.0459             0.117
# 7 Tektite.Day4       0.0601             0.153
# 8 Tektite.Day5       0.0486             0.123
# 9 Yawzi.Day2         0.0874             0.222
# 10 Yawzi.Day3         0.0527             0.134
# 11 Yawzi.Day4         0.0657             0.167
# 12 Yawzi.Day5         0.0508             0.129



# dist_usvi_asv_site_day.boot <- boot::boot(dist_usvi_asv_site_day.df, F_ttest_gn, variable = "site_day", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
resample_depth <- dist_usvi_asv_site_day.df %>%
  split(., f = .$site_day) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_asv_site_day.boot <- boot::boot(dist_usvi_asv_site_day.df, F_ttest_para, variable = "site_day", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)

# q_value <- 0.10
dist_usvi_asv_site_day.boot.summary.df <- unique(dist_usvi_asv_site_day.df[["site_day"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_asv_site_day.boot$t0, name = NULL, value = "boot.p.value")) %>% #this value is essentially random...
  # bind_cols(., tibble::enframe(dist_usvi_asv_site_day.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q05" = (dist_usvi_asv_site_day.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_asv_site_day.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
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
  bind_rows(., .id = "site_day") %>%
  dplyr::group_by(site_day)


dist_usvi_metab_site_day.df %>%
  dplyr::group_by(site_day) %>%
  # dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE),
  #                  min_dist = min(dissimilarity, na.rm = TRUE),
  #                  max_dist = max(dissimilarity, na.rm = TRUE)) %>%
  dplyr::summarise(avg_dist = mean(dissimilarity, na.rm = TRUE)) %>%
  dplyr::mutate(across(contains("dist"), list(rescaled = ~.x/max(dist_usvi_metab_site_day.df[["dissimilarity"]]))))
# site_day         avg_dist avg_dist_rescaled
# <chr>               <dbl>             <dbl>
# 1 LB_seagrass.Day2   0.109              0.197
# 2 LB_seagrass.Day3   0.0974             0.177
# 3 LB_seagrass.Day4   0.103              0.188
# 4 LB_seagrass.Day5   0.203              0.368
# 5 Tektite.Day2       0.317              0.575
# 6 Tektite.Day3       0.275              0.499
# 7 Tektite.Day4       0.199              0.361
# 8 Tektite.Day5       0.228              0.413
# 9 Yawzi.Day2         0.176              0.320
# 10 Yawzi.Day3         0.148              0.268
# 11 Yawzi.Day4         0.173              0.314
# 12 Yawzi.Day5         0.206              0.374


# dist_usvi_metab_site_day.boot <- boot::boot(dist_usvi_metab_site_day.df, F_ttest_gn, variable = "site_day", sim = "permutation", R = 9999, stype = "f", parallel = "multicore", ncpus = nthreads)
resample_depth <- dist_usvi_metab_site_day.df %>%
  split(., f = .$site_day) %>%
  map(., ~round(nrow(.x)/10)) %>%
  purrr::reduce(mean)
dist_usvi_metab_site_day.boot <- boot::boot(dist_usvi_metab_site_day.df, F_ttest_para, variable = "site_day", sim = "parametric", R = 9999, ran.gen = F_ttest_rn, mle = resample_depth, parallel = "multicore", ncpus = nthreads)


# q_value <- 0.10
dist_usvi_metab_site_day.boot.summary.df <- unique(dist_usvi_metab_site_day.df[["site_day"]]) %>%
  combn(., 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  bind_cols(., tibble::enframe(dist_usvi_metab_site_day.boot$t0, name = NULL, value = "boot.p.value")) %>%
  # bind_cols(., tibble::enframe(dist_usvi_metab_site_day.boot$t0, name = NULL, value = "p.value")) %>%
  bind_cols(., data.frame("q05" = (dist_usvi_metab_site_day.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t()))) %>%
  bind_cols(., data.frame("q10" = (dist_usvi_metab_site_day.boot$t %>%
                                     as.data.frame() %>%
                                     dplyr::reframe(across(everything(), quantile, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7)) %>%
                                     t())))

length_v <- unique(dist_usvi_asv_site_day.df[["site_day"]]) %>% as.character(.)
ttest_res <- combn(length_v, 2) %>%
  t() %>%
  as.data.frame() %>%
  setNames(., c("pair1", "pair2")) %>%
  dplyr::mutate(t = NA) %>%
  dplyr::mutate(`p.value` = NA)
ttest_usvi_metab_site_day.res <- ttest_res
ttest_usvi_asv_site_day.res <- ttest_res
for(i in seq_len(nrow(ttest_res))){
  var1 <- ttest_res[i, 1]
  var2 <- ttest_res[i, 2]
  temp_i_a <- dist_usvi_metab_site_day.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_b <- dist_usvi_metab_site_day.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  if(length(temp_i_a) > 0 & length(temp_i_b) > 0){
    ttest_usvi_metab_site_day.res[i, "t"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$statistic  
    ttest_usvi_metab_site_day.res[i, "p.value"] <- t.test(temp_i_a, temp_i_b, conf.level = 0.95)$p.value  
  } else {
    ttest_usvi_metab_site_day.res[i, "t"] <- NA
  }

  temp_i_c <- dist_usvi_asv_site_day.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var1, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  temp_i_d <- dist_usvi_asv_site_day.df %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(if_any(everything(), ~grepl(var2, .x))) %>%
    dplyr::select(dissimilarity) %>%
    tibble::deframe(.)
  if(length(temp_i_c) > 0 & length(temp_i_d) > 0){
    ttest_usvi_asv_site_day.res[i, "t"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$statistic
    ttest_usvi_asv_site_day.res[i, "p.value"] <- t.test(temp_i_c, temp_i_d, conf.level = 0.95)$p.value
  } else {
    ttest_usvi_asv_site_day.res[i, "t"] <- NA
  }
}

ttest_usvi_metab_site_day.res <- ttest_usvi_metab_site_day.res %>%
  # tidyr::drop_na(t) %>%
  dplyr::left_join(., dist_usvi_metab_site_day.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))
ttest_usvi_asv_site_day.res <- ttest_usvi_asv_site_day.res %>%
  # tidyr::drop_na(t) %>%
  dplyr::left_join(., dist_usvi_asv_site_day.boot.summary.df) %>%
  dplyr::mutate(sig = dplyr::case_when(p.value < q05 ~ p.value, .default = NA))


readr::write_delim(ttest_usvi_metab_site_day.res, paste0(projectpath, "/", "permanova_ttest_usvi_metab_site_day.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)
readr::write_delim(ttest_usvi_asv_site_day.res, paste0(projectpath, "/", "permanova_ttest_usvi_asv_site_day.res-", Sys.Date(), ".tsv"),
                   delim = "\t", col_names = TRUE)

ttest_usvi_asv_site_day.res %>%
  dplyr::filter(p.value < 0.05) %>%
  # dplyr::summarise(num_SDA = length(t))
  droplevels

ttest_usvi_metab_site_day.res %>%
  dplyr::filter(p.value < 0.05) %>%
  # dplyr::summarise(num_SDA = length(t))
  droplevels
#2 comparisons are between days at LB
#2 comparisons are between Yawzi and LB on days 3 and 4
#1 comparison is between Yawzi and Tektite on day 4
save(dist_usvi_metab_site_day.boot, dist_usvi_metab_site_day.boot.summary.df, dist_usvi_metab_site_day.df,
     dist_usvi_asv_site_day.boot, dist_usvi_asv_site_day.boot.summary.df, dist_usvi_asv_site_day.df,
     ttest_usvi_metab_site_day.res, ttest_usvi_asv_site_day.res,
     file = paste0(projectpath, "/", "usvi_permanova_site_day-", Sys.Date(), ".RData"))

}
###STOP HERE: 20250208


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





# Segue: BIOS cruise 2023 July --------------------------------------------

# ctd_data.df <- readr::read_delim("~/projects/kujawinski/biosscope/ctd_50060_data.txt", 
#                                  # comment = "%",
#                                  quote = "\"", delim = ";", trim_ws = TRUE,
#                                  skip = 1,
#                                  col_names = TRUE, show_col_types = FALSE) %>%
#   droplevels
# new_names <- colnames(ctd_data.df)
# ctd_data.df <- ctd_data.df %>%
#   dplyr::rename(original = `% Cruise/Cast ID`) %>%
#   dplyr::mutate(original = stringr::str_squish(original)) %>%
#   dplyr::select(original) %>%
#   tidyr::separate_wider_delim(original, names = new_names, delim = " ", too_few = "align_start", too_many = "merge") %>%
#   droplevels 
# 
# ctd_metadata.df <- readr::read_delim("~/projects/kujawinski/biosscope/ctd_50060_meta.txt",
#                                  # comment = "%",
#                                  quote = "\"", delim = ";", trim_ws = TRUE,
#                                  skip = 1,
#                                  col_names = TRUE, show_col_types = FALSE) %>%
#   droplevels
# new_names <- colnames(ctd_metadata.df)[1:12]
# new_names[1] <- "CastID"
# ctd_metadata.df <- ctd_metadata.df %>%
#   dplyr::rename(original = `% Cruise/Cast ID`) %>%
#   dplyr::mutate(original = stringr::str_squish(original)) %>%
#   dplyr::select(original) %>%
#   tidyr::separate_wider_delim(original, names = new_names, delim = " ", too_few = "align_start", too_many = "merge") %>%
#   dplyr::mutate(across(everything(), ~as.numeric(.x))) %>%
#   dplyr::mutate(time_in = lubridate::date_decimal(`Decimal Year-in`, tz = "Atlantic/Bermuda")) %>%
#   dplyr::mutate(time_out = lubridate::date_decimal(`Decimal Year-out`, tz = "Atlantic/Bermuda")) %>%
#   droplevels
# 
# temp_df <- ctd_data.df %>%
#   dplyr::select(`% Cruise/Cast ID`, `Depth(m)`, `Flu (rfu)`, `PAR (uE/m^2)`) %>%
#   setNames(., c("CastID", "Depth", "Flu", "PAR")) %>%
#   tidyr::drop_na(.) %>%
#   dplyr::mutate(across(everything(), ~as.numeric(.x))) %>%
#   dplyr::left_join(., ctd_metadata.df %>%
#                      dplyr::select(CastID, time_in) %>%
#                      droplevels, by = join_by(CastID)) %>%
#   dplyr::mutate(Flu = Flu*10000) %>%
#   # dplyr::mutate(PAR = (PAR + 1)) %>%
#   dplyr::filter(if_all(everything(), ~.x > -1)) %>%
#   dplyr::mutate(CastID = factor(CastID),
#                 time_in = factor(time_in)) %>%
#   tidyr::pivot_longer(., cols = !c("CastID", "Depth", "time_in"),
#                       names_to = "parameter",
#                       values_to = 'value') %>%
#   droplevels
# 
# 
# temp_g <- print(ggplot(data= temp_df %>%
#                dplyr::filter(Depth <= 50), aes(x = value, y = Depth, group = interaction(CastID, parameter)))
#       + geom_path(aes(color = parameter), linewidth = 1)
#       + geom_point(aes(fill = parameter), shape = 21) 
#       + scale_discrete_manual(aesthetics = c("fill", "color"),
#                               values = viridis::viridis(option = "viridis", begin = 0.5, end = 1, n = length(unique(temp_df[["parameter"]]))))
#       + theme_bw()
#       + facet_grid(.~time_in, scales = "free")
#       + scale_x_continuous(name = "PAR", 
#                            labels = scaleFUN1,
#                            # transform = scales::transform_pseudo_log(sigma = 10, base = 10),
#                            transform = "log10",
#                            sec.axis = sec_axis(~./10000, name = "Flu",  labels = scales::label_number(accuracy = 0.001)))
#       + scale_y_continuous(n.breaks = 10, name = "Depth (m)", transform = "reverse", expand = expansion(mult = c(0, 0.1)))
#       + theme(axis.text.x = element_text(angle = -90))
#       )
# 
# ggsave(paste0("~/projects/kujawinski/biosscope/ctd_50060_par_flu", ".png"),
#        temp_g,
#        width = 10, height = 6, units = "in")

# temp_df2 <- ctd_data.df %>%
#   dplyr::select(`% Cruise/Cast ID`, `Depth(m)`, `Flu (rfu)`, `PAR (uE/m^2)`) %>%
#   setNames(., c("CastID", "Depth", "Flu", "PAR")) %>%
#   tidyr::drop_na(.) %>%
#   dplyr::mutate(across(everything(), ~as.numeric(.x))) %>%
#   dplyr::left_join(., ctd_metadata.df %>%
#                      dplyr::select(CastID, time_in) %>%
#                      droplevels, by = join_by(CastID)) %>%
#   dplyr::filter(if_all(everything(), ~.x > -1)) %>%
#   dplyr::mutate(CastID = factor(CastID),
#                 time_in = factor(time_in)) %>%
#   tidyr::pivot_longer(., cols = !c("CastID", "Depth", "time_in"),
#                       names_to = "parameter",
#                       values_to = 'value') %>%
#   droplevels
# print(ggplot(data= temp_df2 %>%
#                dplyr::filter(Depth <= 50), aes(x = value, y = Depth, group = interaction(CastID, parameter)))
#       + geom_path(aes(color = time_in), linewidth = 1)
#       + geom_point(aes(fill = time_in), shape = 21) 
#       + theme_bw()
#       # + facet_grid(.~CastID, scales = "free")
#       + facet_grid(.~parameter, scales = "free", space = "fixed")
#       + scale_x_continuous(name = "measurement")
#       + scale_discrete_manual(aesthetics = c("fill", "color"),
#                               values = viridis::viridis(option = "turbo", n = length(unique(temp_df[["time_in"]]))))
#       + scale_y_continuous(n.breaks = 10, name = "Depth (m)", transform = "reverse", expand = expansion(mult = c(0, 0.1)))
#       + theme(axis.text.x = element_text(angle = -90))
# )
