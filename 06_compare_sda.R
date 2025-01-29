# 06_compare_sda.R

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
library(DESeq2)
library(edgeR)
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)
library(ggVennDiagram)


# Custom functions --------------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_custom_functions.R"))){
  cat("Loading custom functions.")
  source(paste0(getwd(), "/", "00_custom_functions.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 


set.alpha <- 0.05
set.seed(48105)


# find existing processed files -------------------------------------------
to_import <- c("usvi_prok_asvs.df", "usvi_prok_asvs.taxa", "usvi_prok_decontam_idx")

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


# usvi_metabolomics.df <- readr::read_delim(paste0(projectpath, "/", "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_wideFormat_noMetadata.csv"),
#                                           col_names = TRUE, show_col_types = FALSE, delim = ",", num_threads = nthreads)
# colnames(usvi_metabolomics.df)[1] <- "metab_deriv_label"
# 
# 
# #there are samples "CINAR_BC_81A" and "CINAR_BC_81B" in the metabolomics dataset 
# #and in the metadata, there are two DNA samples associated with "Deriv_81": Metab_219 (LB_seagrass dawn) and Metab_319 (tektite dawn)
# 
# usvi_metabolomics_long.df <- readr::read_delim(paste0(projectpath, "/", "USVI2021_CINARtemporal_BzCl_Exometabolite_QCd_longFormat_wMetadata.csv"), 
#                                                col_select = c(2:last_col()),
#                                                col_names = TRUE, show_col_types = FALSE, delim = ",", num_threads = nthreads) %>%
#   dplyr::mutate(sample_id = paste0("Metab_", DNA_no))
# 
# # long metabolomics dataset from Brianna confirms that CINAR_BC_81A is the BC sample associated with Tektite Metab_319
# # and CINAR_BC_81B is the BC sample associated with LB_seagrass Metab_219
# 
# 
drop <- c("CINAR_BC_73")

usvi_selected_metadata <- metabolomics_sample_metadata %>%
  dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::select(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
  dplyr::mutate(across(c(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site), ~factor(.x))) %>%
  # dplyr::select(sample_id, metab_deriv_label) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(site_type = dplyr::case_when(grepl("LB", site) ~ "seagrass",
                                             grepl("Yawzi|Tektite", site) ~ "reef",
                                             .default = NA)) %>%
  # dplyr::mutate(rownames = metab_deriv_label) %>% tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::mutate(rownames = sample_id) %>% tibble::column_to_rownames(var = "rownames") %>%
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
  contrast = contrast_labels_lookup,
  enriched = enriched_labels_lookup,
  .multi_line = TRUE,
  .default = label_wrap_gen(25, multi_line = TRUE)
  # .default = label_value
)




# Read in results from radEmu ---------------------------------------------
#since RadEmu was only at the genera-level, don't read it in.
# 
# 
# #read in the genera that were significant by site through rademu:
# if(any(grepl("genera_reef_highlow", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
#   # temp_file <- list.files(projectpath, pattern = "usvi_rademu_genera_reef_highlow.*.RData")[1]
#   temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_rademu_genera_reef_highlow.*.RData"))
#   load(paste0(projectpath, "/", temp_file))
#   rm(temp_file)
# }
# 
# # rademu_reef_res.df <- score_lowreef_res_list %>%
# #   bind_rows(.) %>%
# #   dplyr::mutate(group = "low_reef") %>%
# #   bind_rows(., score_highreef_res_list %>%
# #               bind_rows(.) %>%
# #               dplyr::mutate(group = "high_reef")) %>%
# #   dplyr::rename(asv_id = "category") %>%
# #   droplevels
# 
# if(any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
#   # temp_file <- list.files(projectpath, pattern = "usvi_rademu_genera_site_time.*.RData")[1]
#   temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_rademu_genera_site_time.*.RData"))
#   load(paste0(projectpath, "/", temp_file))
#   rm(temp_file)
# }
# 
# rademu_res.df <- bind_rows((score_1_res_list %>%
#                               bind_rows(.) %>% #these were genera that were found to be higher in samples collected at dawn
#                               dplyr::mutate(group = "high_dawn") %>%
#                               dplyr::filter(pval < 0.05) %>%
#                               tidyr::drop_na(score_stat)),
#                            (score_2_res_list %>%
#                               bind_rows(.) %>% # these were genera that were found to be higher in samples collected at peak photo
#                               dplyr::mutate(group = "low_dawn") %>%
#                               dplyr::filter(pval < 0.05) %>%
#                               tidyr::drop_na(score_stat)),
#                            (score_4_res_list %>%
#                               bind_rows(.) %>% #these were genera that were found to be higher in samples collected in seagrass
#                               dplyr::mutate(group = "high_seagrass") %>%
#                               dplyr::filter(pval < 0.05) %>%
#                               tidyr::drop_na(score_stat)),
#                            (score_5_res_list %>%
#                               bind_rows(.) %>% #these were genera that were found to be higher in samples collected in reef sites
#                               dplyr::mutate(group = "low_seagrass") %>%
#                               dplyr::filter(pval < 0.05) %>%
#                               tidyr::drop_na(score_stat)),
#                            (score_lowreef_res_list %>%
#                               bind_rows(.) %>%
#                               dplyr::mutate(group = "low_reef") %>%
#                               dplyr::filter(pval < 0.05) %>%
#                               tidyr::drop_na(score_stat)),
#                            (score_highreef_res_list %>%
#                               bind_rows(.) %>%
#                               dplyr::mutate(group = "high_reef") %>%
#                               dplyr::filter(pval < 0.05) %>%
#                               tidyr::drop_na(score_stat))) %>%
#   dplyr::rename(asv_id = "category") %>%
#   droplevels %>%
#   dplyr::rename(contrast = "covariate") %>%
#   dplyr::mutate(contrast = stringr::str_remove_all(contrast, "site") %>%
#                   stringr::str_replace_all(., ":sampling_time", "_")) %>%
#   dplyr::mutate(contrast = dplyr::case_when(grepl("reef", group) ~ paste0(contrast, " - LB_seagrass_all"),
#                                             grepl("dawn", group) ~ paste0(contrast, " - all_peak_photo"),
#                                             grepl("seagrass", group) ~ paste0(contrast, " - LB_seagrass_all"),
#                                             .default = contrast)) %>%
#   dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
#   dplyr::mutate(commonid = stringr::str_split_i(group, "_", 2)) %>%
#   dplyr::arrange(commonid, asv_id) %>%
#   dplyr::group_by(commonid) %>%
#   dplyr::distinct(asv_id, group, contrast, .keep_all = TRUE) %>%
#   droplevels
# 
# #for the temporal nuances:
# #subtract all genera that were just significantly higher/lower in seagrass regardless of time
# # temp_list <- rademu_res.df %>%
# #   split(., f = .$commonid) %>%
# #   map(., ~.x %>%
# #         droplevels)
# rademu_seagrass_res.df <- rademu_res.df %>%
#   dplyr::filter(commonid == "seagrass") %>%
#   dplyr::distinct(asv_id, group, .keep_all = TRUE) %>%
#   dplyr::arrange(asv_id) %>%
#   droplevels
# 
# rademu_reef_res.df <- rademu_res.df %>%
#   dplyr::filter(commonid == "reef") %>%
#   dplyr::distinct(asv_id, group, .keep_all = TRUE) %>%
#   dplyr::arrange(asv_id) %>%
#   droplevels
# 
# rademu_time_res.df <- rademu_res.df %>%
#   dplyr::filter(commonid == "dawn") %>%
#   dplyr::distinct(asv_id, group, .keep_all = TRUE) %>%
#   dplyr::arrange(asv_id) %>%
#   droplevels
# 
# #these are the genera that were shared between high/low reef and high/low seagrass (23):
# rademu_site_genera_idx <- intersect(rademu_reef_res.df[["asv_id"]], rademu_seagrass_res.df[["asv_id"]])
# # union(rademu_reef_res.df[["asv_id"]], rademu_seagrass_res.df[["asv_id"]])
# 
# # #these are the genera that were not shared between high/low reef and high/low seagrass (13):
# # setdiff(rademu_reef_res.df[["asv_id"]], rademu_seagrass_res.df[["asv_id"]])
# 
# 
# 
# #now find the difference between the genera flagged as spatiotemporally different, and the site different:
# # setdiff(rademu_reef_res.df[["asv_id"]], rademu_time_res.df[["asv_id"]])
# setdiff(rademu_time_res.df[["asv_id"]], rademu_site_genera_idx)
# # [1] "ASV_00008" "ASV_00012" "ASV_00215" "ASV_00227" "ASV_00243" "ASV_00282" "ASV_01134"
# 




# Read in results from DESeq2 ---------------------------------------------

#proceed with the results at the ASV-level

if(any(grepl("asvs_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_asvs_site_time.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
}
#usvi_deseq_asvs_res.list, usvi_deseq_asvs_abund_filtered.df, usvi_deseq_asvs_res.df

# if(!any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
#   temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_genera_site_time.*.RData"))
#   load(paste0(projectpath, "/", temp_file))
#   rm(temp_file)
# }
# #usvi_deseq_genera_res.list, usvi_deseq_genera_abund_filtered.df, usvi_deseq_genera_abund.df, usvi_deseq_genera_res.df

if(any(grepl("asvs_intrasite_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_asvs_intrasite_time-.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
  #loads in: usvi_deseq_asvs_intrasite_res.list, usvi_deseq_asvs_intrasite_abund.df, usvi_deseq_asvs_intrasite_res.df
}


# Read in results from corncob --------------------------------------------

#read in the genera that were significant via corncob:
# if(any(grepl("cc_dt_sda", list.files(projectpath, pattern = ".*ps_usvi_filtered-.*.RData")))){
#   temp_file <- data.table::last(list.files(projectpath, pattern = "cc_dt_sda_ps_usvi_filtered*.RData"))
#   load(paste0(projectpath, "/", temp_file))
#   rm(temp_file)
# }

if(any(grepl("cc_dt_usvi_summary.df-", list.files(projectpath, pattern = "^cc_dt_usvi_summary.*.tsv")))){
# if(file.exists(paste0(projectpath, "/", "cc_dt_usvi_summary.df", ".tsv"))){
    temp_file <- data.table::last(list.files(projectpath, pattern = "^cc_dt_usvi_summary.df-.*.tsv"))
    cc_dt_usvi_summary.df <- readr::read_delim(paste0(projectpath, "/", temp_file),
  # cc_dt_usvi_summary.df <- readr::read_delim(paste0(projectpath, "/", "cc_dt_usvi_summary.df", ".tsv"),
                                             delim = "\t", col_names = TRUE, show_col_types = FALSE)
}


# Compare SDA results from the different methods --------------------------

#Genus level--old

{
# 
# 
# #since rademu couldn't proceed at the ASV level, use just corncob results at the agglomerated genera-level 
# 
# corncob_res.df <- cc_dt_usvi_summary.df %>%
#   dplyr::filter(grepl("agg_genus", taxon_resolution)) %>%
#   dplyr::mutate(test_type = "corncob") %>%
#   dplyr::rename(padj = "p_adj") %>%
#   dplyr::rename(model = "test") %>%
#   dplyr::select(contrast, hold, variable, asv_id, p_fdr, padj, taxonomy, test_type, model) %>%
#   dplyr::distinct(asv_id, hold, variable, .keep_all = TRUE) %>%
#   # tidyr::separate_longer_delim(variable, delim = " - ") %>%
#   dplyr::mutate(pair1 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\1", variable), "_", hold),
#                                          .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\1", variable)))) %>%
#   dplyr::mutate(pair2 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\3", variable), "_", hold),
#                                          .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)))) %>%
#   # dplyr::mutate(contrast = paste0(pair1, " - ", pair2)) %>%
#   dplyr::select(contrast, hold, variable, asv_id, padj, p_fdr, test_type, model) %>%
#   # dplyr::mutate(sampling_time = dplyr::case_when(grepl("dawn|peak", contrast) ~ stringr::str_extract(contrast, "dawn|peak_photo"),
#   #                                                .default = "all")) %>%
#   # dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
#   # dplyr::mutate(Combo = dplyr::case_when(grepl("Yawzi_|Tektite_", contrast) ~ stringr::str_extract(contrast, "Yawzi_.|Tektite_."),
#   #                                        .default = site)) %>%
#   # dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
#   # dplyr::mutate(baseline = dplyr::case_when(grepl("all_peak_photo", baseline) ~ "all_peak_photo",
#   #                                           grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
#   #                                           .default = "LB_seagrass_all")) %>%
#   droplevels
# 
# 
# usvi_sda_genera_compare.df <- rademu_res.df %>%
#   dplyr::mutate(test_type = "rademu") %>%
#   dplyr::rename(estimate = "score_stat") %>%
#   dplyr::mutate(model = "auto_fit") %>%
#   dplyr::mutate(sampling_time = dplyr::case_when(grepl("dawn|peak", contrast) ~ stringr::str_extract(contrast, "dawn|peak_photo"),
#                                                  .default = "all")) %>%
#   dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
#   dplyr::mutate(Combo = dplyr::case_when(grepl("Yawzi_|Tektite_", contrast) ~ stringr::str_extract(contrast, "Yawzi_.|Tektite_."),
#                                          .default = site)) %>%
#   dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
#   dplyr::mutate(baseline = dplyr::case_when(grepl("all_peak_photo", baseline) ~ "all_peak_photo",
#                                             grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
#                                             .default = "LB_seagrass_all")) %>%
#   bind_rows(., (usvi_deseq_genera_abund.df %>%
#                   # bind_rows(., (usvi_deseq_genera_abund_filtered.df %>%
#                   dplyr::select(-padj_qcutoff) %>%
#                   dplyr::mutate(test_type = "deseq"))) %>%
#   dplyr::mutate(pair1 = stringr::str_split_i(contrast, " - ", 1),
#                 pair2 = stringr::str_split_i(contrast, " - ", 2)) %>%
#   bind_rows(., corncob_res.df)
# 
# 
# 
# #which genera are flagged as significantly differentially abundant in all tests?
# shared_sda_genera_idx <- usvi_sda_genera_compare.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
#   dplyr::group_by(asv_id) %>%
#   dplyr::summarise(num_results = length(test_type)) %>%
#   # dplyr::filter(num_results > 1) %>%
#   dplyr::filter(num_results > 2) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::distinct(asv_id) %>%
#   unlist %>%
#   as.character
# # droplevels
# 
# 
# #which are not?
# selfish_sda_genera_idx <- usvi_sda_genera_compare.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
#   dplyr::group_by(asv_id) %>%
#   dplyr::summarise(num_results = length(test_type)) %>%
#   dplyr::filter(num_results <= 1) %>%
#   # dplyr::filter(num_results <= 2) %>%
#   dplyr::distinct(asv_id) %>%
#   dplyr::left_join(., usvi_sda_genera_compare.df %>%
#                      dplyr::distinct(asv_id, test_type, .keep_all = FALSE),
#                    by = join_by(asv_id)) %>%
#   dplyr::select(test_type, asv_id) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::distinct(asv_id, .keep_all = TRUE) %>%
#   tibble::deframe(.)
# # unlist %>%
# # as.character
# 
# usvi_sda_genera_compare.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(asv_id, test_type, pval, padj, p_fdr, .keep_all = FALSE) %>%
#   dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
#   dplyr::select(asv_id, test_type, p_value) %>%
#   dplyr::distinct(asv_id, test_type, .keep_all = TRUE) %>%
#   tidyr::pivot_wider(., id_cols = "asv_id",
#                      names_from = "test_type",
#                      values_from = "p_value") %>%
#   dplyr::mutate(across(c("rademu", "deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
#                                                                           .default = NA))) %>%
#   tidyr::unite(test_type, c("rademu", "deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
#   droplevels
# 
# semiself_sda_genera_idx <- usvi_sda_genera_compare.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
#   dplyr::group_by(asv_id) %>%
#   dplyr::summarise(num_results = length(test_type)) %>%
#   # dplyr::filter(num_results <= 1) %>%
#   dplyr::filter(num_results <= 2) %>%
#   dplyr::distinct(asv_id) %>%
#   dplyr::anti_join(., tibble::enframe(selfish_sda_genera_idx, name = "test_type", value = "asv_id")) %>%
#   dplyr::left_join(., (usvi_sda_genera_compare.df %>%
#                          dplyr::ungroup(.) %>%
#                          dplyr::distinct(asv_id, test_type, pval, padj, p_fdr, .keep_all = FALSE) %>%
#                          dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
#                          dplyr::select(asv_id, test_type, p_value) %>%
#                          dplyr::distinct(asv_id, test_type, .keep_all = TRUE) %>%
#                          tidyr::pivot_wider(., id_cols = "asv_id",
#                                             names_from = "test_type",
#                                             values_from = "p_value") %>%
#                          dplyr::mutate(across(c("rademu", "deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
#                                                                                                  .default = NA))) %>%
#                          tidyr::unite(test_type, c("rademu", "deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
#                          droplevels),
#                    # dplyr::left_join(., usvi_sda_genera_compare.df %>%
#                    #                    dplyr::select(asv_id, test_type) %>%
#                    #                    tidyr::pivot_wider(., id_cols = "asv_id",
#                    #                                       names_from = "metric",
#                    #                                       values_from = "test_type") %>%
#                    #                    dplyr::distinct(asv_id, test_type, .keep_all = FALSE),
#                    by = join_by(asv_id)) %>%
#   dplyr::select(test_type, asv_id) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::distinct(asv_id, .keep_all = TRUE) %>%
#   tibble::deframe(.)
# 
# #save the results as a table listing which genera are significant in which comparison/test
# 
# usvi_sda_genera_compare_summary.df <- usvi_sda_genera_compare.df %>%
#   dplyr::distinct(asv_id, contrast, test_type) %>%
#   dplyr::arrange(asv_id) %>%
#   droplevels %>%
#   tidyr::pivot_wider(., id_cols = "asv_id",
#                      names_from = c(contrast),
#                      values_from = c(test_type),
#                      values_fill = NA) %>%
#   droplevels %>%
#   dplyr::left_join(., selfish_sda_genera_idx %>%
#                      tibble::enframe(name = "test_type.x", value = "asv_id")) %>%
#   dplyr::left_join(., semiself_sda_genera_idx %>%
#                      tibble::enframe(name = "test_type.z", value = "asv_id")) %>%
#   dplyr::left_join(., data.frame(asv_id = shared_sda_genera_idx,
#                                  test_type.y = "all")) %>%
#   tidyr::unite("test_type", c("test_type.x", "test_type.z", "test_type.y"), remove = TRUE, na.rm = TRUE) %>%
#   dplyr::left_join(., usvi_sw_genus.taxa.df, by = join_by(asv_id)) %>%
#   droplevels
# 
# readr::write_delim(usvi_sda_genera_compare_summary.df, paste0(projectpath, "/", "usvi_sda_genera_compare_summary.tsv"),
#                    delim = "\t", col_names = TRUE, na = "")
# if(!any(grepl("rds$", list.files(projectpath, pattern = "usvi_sda_genera_compare.*.rds")))){
#   readr::write_rds(usvi_sda_genera_compare.df, 
#                    paste0(projectpath, "/", "usvi_sda_genera_compare.rds"))
#   #also save the results with the pvals and estimates
#   readr::write_delim(usvi_sda_genera_compare.df, 
#                      paste0(projectpath, "/", "usvi_sda_genera_compare.tsv"),
#                      delim = "\t", col_names = TRUE, na = "")
# }
# 
# 
# 
# # Plot the genera that were SDA across all 3 methods ----------------------
# 
# shared_sda_genera_idx_cluster <- usvi_sda_genera_compare.df %>%
#   dplyr::filter(asv_id %in% shared_sda_genera_idx) %>%
#   dplyr::distinct(contrast, asv_id, .keep_all = FALSE) %>%
#   dplyr::group_by(asv_id) %>%
#   dplyr::summarise(num_results = length(contrast)) %>%
#   dplyr::arrange(num_results) %>%
#   dplyr::mutate(grouping = c(rep(1:3, each = 5, length.out = length(shared_sda_genera_idx)))) %>%
#   dplyr::arrange(asv_id) %>%
#   split(., f = .$grouping) %>%
#   map(., ~.x %>%
#         dplyr::select(asv_id, num_results) %>%
#         droplevels)
# 
# rm(list = apropos("g_sda_.*", mode = "list"))
# 
# for(i in seq_along(shared_sda_genera_idx_cluster)){
#   namevar <- shared_sda_genera_idx_cluster[[i]] %>%
#     purrr::pluck("asv_id")
#   temp_df2 <- usvi_sda_genera_compare.df %>%
#     dplyr::filter(asv_id %in% namevar) %>%
#     droplevels
#   g <- print(
#     ggplot(data = temp_df2)
#     + theme_bw() 
#     + geom_boxplot(aes(y = relabund, x = sampling_time, group = sampling_time), 
#                    alpha = 0.6, show.legend = FALSE, color = "black",
#                    position = position_dodge(0.8))
#     + geom_point(aes(y = relabund, x = sampling_time, group = sampling_time, fill = sampling_time, shape = site),
#                  alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
#                  position = position_jitter(width = 0.2, seed = 48105))
#     + scale_shape_manual(name = "Sampling site and time",
#                          values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#     + scale_fill_manual(name = "Sampling time",
#                         values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#     + scale_x_discrete(labels = sampling_time_lookup, name = "Sampling time")
#     + scale_y_continuous(name = "Relative abundance (%)")
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                  override.aes = list(color = "black", shape = 21, size = 3)),
#              color = "none")
#     # + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
#     + theme(axis.text.x = element_blank(),
#             strip.text.y = element_text(size = rel(0.7)),
#             axis.title.y = element_blank())
#     + facet_grid(asv_id ~ site,
#                  drop = TRUE,
#                  # labeller = labeller(site = site_lookup, asv_id = stringr::str_wrap(usvi_genera_relabel)),
#                  labeller = global_labeller,
#                  scales = "free")
#     + ggtitle(paste0("Group ", i))
#   )
#   assign(paste0("g_sda_", i), g, envir = .GlobalEnv, inherits = TRUE)
#   rm(g)
#   rm(namevar)
# }
# 
# gpatch <- apropos("g_sda_.*", mode = "list") %>%
#   lapply(., get) %>%
#   # purrr::reduce(., `+ theme(legend.position = "none")`) %>%
#   purrr::reduce(., `+`)
# g2_sda_genera <- gpatch + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")
# g2_sda_genera
# 
# if(!any(grepl("genera_relabund", list.files(projectpath, pattern = "usvi_sda_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_sda_genera_relabund-", Sys.Date(), ".png"),
#          g2_sda_genera,
#          width = 15, height = 10, units = "in")
# }

}

# Plot the relative abundances of ASVs found sig --------------------------


corncob_asv_res.df <- cc_dt_usvi_summary.df %>%
  dplyr::filter(grepl("asv", taxon_resolution)) %>%
  dplyr::filter(p_fdr <= 0.05) %>%
  dplyr::mutate(test_type = "corncob") %>%
  dplyr::rename(padj = "p_adj") %>%
  dplyr::rename(model = "test") %>%
  dplyr::select(contrast, hold, variable, asv_id, p_fdr, padj, q_cutoff, test_type, model) %>%
  dplyr::distinct(asv_id, hold, variable, .keep_all = TRUE) %>%
  dplyr::mutate(contrast = dplyr::case_when(grepl("Tektite - Yawzi", contrast) ~ stringr::str_replace(contrast, "Tektite - Yawzi", "Yawzi - Tektite"),
                                            .default = contrast)) %>%
  dplyr::mutate(variable = dplyr::case_when(grepl("Tektite - Yawzi", variable) ~ "Yawzi - Tektite",
                                            .default = variable)) %>%
  dplyr::mutate(pair1 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\1", variable), "_", hold),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\1", variable)))) %>%
  dplyr::mutate(pair2 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\3", variable), "_", hold),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)))) %>%
  # dplyr::select(contrast, hold, variable, asv_id, padj, p_fdr, test_type, model) %>%
  droplevels


#stick with "manual_" dispersion model, and filter for only those where baseline is time specific
usvi_deseq_asvs_filtered_res.df <- usvi_deseq_asvs_res.df %>%
  dplyr::filter(grepl("manual_", model)) %>%
  dplyr::filter(grepl("photo|dawn|reef", model)) %>%
  dplyr::mutate(sampling_time = stringr::str_extract(contrast, "dawn|peak_photo")) %>%
  dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi_|Tektite_") %>%
                  stringr::str_remove(., "_")) %>%
  dplyr::mutate(hold = stringr::str_remove_all(model, "manual_|auto_")) %>%
  dplyr::mutate(hold = dplyr::case_when(grepl("lameshur", hold) ~ "LB_seagrass",
                                        grepl("photo", hold) ~ "peak_photo",
                                        grepl("yawzi|tektite|lb", hold) ~ stringr::str_to_sentence(hold),
                                        grepl("reef", model) ~ sampling_time,
                                        .default = hold)) %>%
  dplyr::mutate(variable = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ stringr::str_remove_all(contrast, paste0("_", hold)),
                                            (hold %in% names(site_lookup)) ~ stringr::str_remove_all(contrast, paste0(site, "_")),
                                            .default = NA)) %>%
  dplyr::mutate(pair1 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\1", contrast)),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\1", contrast)))) %>%
  dplyr::mutate(pair2 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)))) %>%
  dplyr::mutate(contrast = paste0(hold, " (", variable, ")")) %>%
  droplevels
usvi_deseq_asvs_filtered_res.df <- usvi_deseq_asvs_filtered_res.df %>%
  bind_rows(., (usvi_deseq_asvs_intrasite_res.df %>%
                  dplyr::filter(grepl("manual_", model)) %>%
                  dplyr::mutate(site = stringr::str_extract(contrast, "LB_seagrass_|Yawzi_|Tektite_") %>%
                                  stringr::str_remove(., "_$")) %>%
                  dplyr::mutate(hold = site) %>%
                  dplyr::mutate(pair1 = dplyr::case_when((hold %in% names(site_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\1", contrast)) %>% stringr::str_remove_all(., paste0(site, "_")),
                                                         .default = NA)) %>%
                  dplyr::mutate(pair2 = dplyr::case_when((hold %in% names(site_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>% stringr::str_remove_all(., paste0(site, "_")),
                                                         .default = NA)) %>%
                  dplyr::mutate(variable = stringr::str_remove_all(contrast, paste0(hold, "_"))) %>%
                  dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
                  dplyr::mutate(contrast = paste0(hold, " (", variable, ")")) %>%
                  droplevels)) %>%
  dplyr::select(contrast, hold, variable, pair1, pair2,  asv_id, baseMean, contains("log2FoldChange"), lfcSE, stat, pvalue, padj, padj_qcutoff, model) %>%
  droplevels

usvi_sda_asvs_compare.df <- bind_rows(
  (usvi_deseq_asvs_filtered_res.df %>% #using the table of all ASVs that were SDA by DESeq2, we have 312 combined from both methods
     tidyr::drop_na(padj_qcutoff) %>%
     dplyr::select(-padj_qcutoff) %>%
     dplyr::mutate(test_type = "deseq"))) %>%
  dplyr::mutate(pair1 = gsub("([[:print:]]+)( - )(.*)$", "\\1", variable),
                pair2 = gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)) %>%
  bind_rows(., corncob_asv_res.df %>%
              dplyr::mutate(pair1 = gsub("([[:print:]]+)( - )(.*)$", "\\1", variable),
                            pair2 = gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)) %>%
              droplevels) %>%
  droplevels



#which ASVs are flagged as significantly differentially abundant in both test types?
#using all DESeq2 results + corncob: 181
#initially picked up 192 because 11 ASVs were SDA in deseq in one contrast, and SDA in corncob in a different contrast, but not the same contrast pair.
# usvi_sda_asvs_compare.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(test_type, contrast, asv_id, .keep_all = FALSE) %>%
#   dplyr::summarise(num_results = length(unique(test_type)), .by = c("asv_id", "contrast")) %>%
#   dplyr::filter(num_results > 1) %>%
#   droplevels

shared_sda_asvs_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, contrast, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id, contrast) %>%
  dplyr::summarise(num_results = length(test_type)) %>%
  dplyr::filter(num_results > 1) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id) %>%
  unlist %>%
  as.character


#howeve,r a lot of these are low abundance:
# temp_df <- usvi_sda_asvs_compare.df %>%
#   dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>%
#   droplevels %>%
#   dplyr::arrange(desc(baseMean)) %>%
#   dplyr::distinct(asv_id, .keep_all = TRUE) %>%
#   dplyr::select(asv_id, baseMean)
# hist(temp_df[["baseMean"]], breaks = 10)

#181 ASVs have a baseMean at least 10
#103 ASVs have a baseMean of 100 or more
#36 ASVs have a baseMean of 500+
#17 ASVs have a baseMean of 1000 or more
# usvi_sda_asvs_compare.df %>% dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>% droplevels %>% dplyr::ungroup(.) %>% dplyr::arrange(desc(baseMean)) %>% dplyr::distinct(asv_id, baseMean) %>% dplyr::ungroup(.) %>%
#   dplyr::filter((baseMean >= 500)) %>%
#   dplyr::distinct(asv_id)

shared_sda_asvs_abund_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>%
  droplevels %>%
  dplyr::distinct(asv_id, baseMean) %>%
  dplyr::ungroup(.) %>%
  dplyr::filter((baseMean >= 500)) %>%
  dplyr::distinct(asv_id) %>%
  unlist %>%
  as.character
  


#which are not?
#using all DESeq2 results + corncob: 215
##  specifically, 189 ASVs were identified as SDA only through DESeq2
##  and 26 ASVs identified as SDA only through corncob

selfish_sda_asvs_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(!(asv_id %in% shared_sda_asvs_idx)) %>%
  dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id, test_type) %>%
    dplyr::summarise(num_results = length(test_type)) %>%
    dplyr::filter(num_results <= 1) %>%
    dplyr::select(test_type, asv_id) %>%
    dplyr::arrange(asv_id) %>%
    dplyr::distinct(asv_id, .keep_all = TRUE) %>%
    tibble::deframe(.)



# usvi_sda_asvs_compare.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(asv_id, test_type, padj, p_fdr, .keep_all = FALSE) %>%
#   dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
#   dplyr::select(asv_id, test_type, p_value) %>%
#   dplyr::distinct(asv_id, test_type, .keep_all = TRUE) %>%
#   tidyr::pivot_wider(., id_cols = "asv_id",
#                      names_from = "test_type",
#                      values_from = "p_value") %>%
#   dplyr::mutate(across(c("deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
#                                                                 .default = NA))) %>%
#   tidyr::unite(test_type, c("deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
#   droplevels

#save the results as a table listing which genera are significant in which comparison/test

usvi_sda_asvs_compare_summary.df <- usvi_sda_asvs_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(contrast, hold, variable, asv_id, test_type,  padj, p_fdr, .keep_all = FALSE) %>%
  dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
  tidyr::drop_na(p_value) %>%
  tidyr::pivot_wider(., id_cols = c(contrast, hold, variable, asv_id),
                     names_from = "test_type",
                     values_from = "p_value") %>%
  # tidyr::drop_na(p_value) %>%
  dplyr::mutate(across(c("deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
                                                                .default = NA))) %>%
  tidyr::unite(test_type, c("deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = c(contrast),
                     values_from = c(test_type),
                     values_fill = NA) %>%
  droplevels %>%
    dplyr::left_join(., selfish_sda_asvs_idx %>%
                       tibble::enframe(name = "test_type.x", value = "asv_id")) %>%
    dplyr::left_join(., data.frame(asv_id = shared_sda_asvs_idx,
                                   test_type.y = "all")) %>%
    tidyr::unite("test_type", c("test_type.x", "test_type.y"), remove = TRUE, na.rm = TRUE) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(-sequence), by = join_by(asv_id)) %>%
  droplevels
# usvi_sda_asvs_compare_summary.df %>%
#   dplyr::select(asv_id, `peak_photo (Yawzi - Tektite)`) %>%
#   tidyr::drop_na(.) %>%
#   dplyr::filter(grepl("^deseq_corncob", `peak_photo (Yawzi - Tektite)`))

# usvi_sda_asvs_compare_summary.df <- usvi_sda_asvs_compare.df %>%
#   dplyr::distinct(asv_id, contrast, test_type) %>%
#   dplyr::arrange(asv_id) %>%
#   droplevels %>%
#   tidyr::pivot_wider(., id_cols = "asv_id",
#                      names_from = c(contrast),
#                      values_from = c(test_type),
#                      values_fill = NA) %>%
#   droplevels %>%
#   dplyr::left_join(., selfish_sda_asvs_idx %>%
#                      tibble::enframe(name = "test_type.x", value = "asv_id")) %>%
#   dplyr::left_join(., data.frame(asv_id = shared_sda_asvs_idx,
#                                  test_type.y = "all")) %>%
#   tidyr::unite("test_type", c("test_type.x", "test_type.y"), remove = TRUE, na.rm = TRUE) %>%
#   dplyr::left_join(., usvi_prok_filled.taxa.df %>%
#                      dplyr::select(-sequence), by = join_by(asv_id)) %>%
#   droplevels

if(!any(grepl("sda_asvs_compare_summary", list.files(projectpath, pattern = "usvi_.*.tsv")))){
  readr::write_delim(usvi_sda_asvs_compare_summary.df, paste0(projectpath, "/", "usvi_sda_asvs_compare_summary-", Sys.Date(), ".tsv"),
                     delim = "\t", col_names = TRUE, na = "")
}


# Review SDA ASVs across methods ------------------------------------------

#there are 396 ASVs identified through either DESeq2 and/or corncob, as SDA between groups of samples
#which of them are most interesting?


usvi_sda_asvs_relabund.df <- usvi_prok_asvs.df %>%
  dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = !c("asv_id"),
                      names_to = "sample_id",
                      values_to = "relabund") %>%
  dplyr::filter(asv_id %in% unique(usvi_sda_asvs_compare.df[["asv_id"]])) %>%
  dplyr::left_join(., usvi_selected_metadata %>%
                     dplyr::select(sample_id, sampling_time, site),
                   by = join_by(sample_id), relationship = "many-to-many", multiple = "all") %>%
  droplevels

#for now (20250110)
#plot only those SDA ASVs from the contrasts:
#Yawzi dawn vs afternoon, Tektite dawn vs afternoon, LB dawn vs afternoon


{
  
  # usvi_sda_genera_shared_samples.df <- usvi_sda_genera_compare.df %>%
  #   dplyr::filter(asv_id %in% shared_sda_genera_idx) %>%
  #   droplevels %>%
  #   dplyr::distinct(asv_id, contrast, test_type, group, baseline, sampling_time, site, pair1, pair2, .keep_all = FALSE) %>%
  #   dplyr::mutate(sampling_time = stringr::str_replace_all(sampling_time, "all", "dawn,peak_photo")) %>%
  #   tidyr::separate_longer_delim(sampling_time, delim = ",") %>%
  #   dplyr::rename(pair1_site ="site",
  #                 pair1_sampling_time = "sampling_time") %>%
  #   # dplyr::mutate(pair2_site = "LB_seagrass",
  #   #               pair2_sampling_time = stringr::str_remove_all(pair2, "LB_seagrass_")) %>%
  #   # dplyr::mutate(pair2_sampling_time = stringr::str_replace_all(pair2_sampling_time, "all", "dawn,peak_photo")) %>%
  #   # tidyr::separate_longer_delim(pair2_sampling_time, delim = ",") %>%
  #   dplyr::left_join(., usvi_selected_metadata %>%
  #                      dplyr::select(sample_id, sampling_time, site),
  #                    by = join_by("pair1_site" == "site", "pair1_sampling_time" == "sampling_time"), relationship = "many-to-many", multiple = "all") %>%
  #   # dplyr::bind_rows(., usvi_selected_metadata %>%
  #   #                    dplyr::filter(grepl("LB", site)) %>%
  #   #                    dplyr::select(sample_id, sampling_time) %>%
  #   #                    dplyr::rename(pair2_sampling_time = "sampling_time")) %>%
  #   droplevels
  
  # temp_df2 <- usvi_sw_genus.tbl %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   t() %>%
  #   as.data.frame(.) %>%
  #   tibble::rownames_to_column(var = "sample_id") %>%
  #   tidyr::pivot_longer(., cols = !c("sample_id"),
  #                       names_to = "asv_id",
  #                       values_to = "relabund") %>%
  #   dplyr::filter(asv_id %in% shared_sda_genera_idx) %>%
  #   dplyr::left_join(., usvi_selected_metadata %>%
  #                      dplyr::select(sample_id, sampling_time, site),
  #                    by = join_by(sample_id), relationship = "many-to-many", multiple = "all") %>%
  #   droplevels
  # usvi_sda_genera_shared_samples_baseline.df <- temp_df2 %>%
  #   dplyr::filter(grepl("LB", site)) %>%
  #   dplyr::filter(asv_id %in% shared_sda_genera_idx) %>%
  #   dplyr::rename(pair2_sampling_time = "sampling_time",
  #                 pair2_site = "site") %>%
  #   dplyr::mutate(pair2 = "dawn,peak_photo,all") %>%
  #   tidyr::separate_longer_delim(pair2, delim = ",") %>%
  #   dplyr::mutate(pair2 = paste0("LB_seagrass_", pair2)) %>%
  #   droplevels
  # 
  # temp_list <- usvi_sda_genera_shared_samples.df %>%
  #   dplyr::filter(asv_id == shared_sda_genera_idx[[1]]) %>%
  #   dplyr::filter(pair1 == "Yawzi_dawn") %>%
  #   droplevels %>%
  #   dplyr::left_join(., temp_df2, by = join_by(sample_id, "pair1_site" == "site", "pair1_sampling_time" == "sampling_time", asv_id),
  #                    relationship = "many-to-many", multiple = "all") %>%
  #   # # dplyr::bind_rows(., temp_df2 %>%
  #   #                    dplyr::filter(asv_id == shared_sda_genera_idx[[1]]) %>%
  #   #                    dplyr::filter(grepl("LB", site)) %>%
  #   #                    droplevels %>%
  #   #                    dplyr::mutate(test_type = "rademu,deseq") %>%
  #   #                    tidyr::separate_longer_delim(test_type, delim = ",") %>%
  #   #                    droplevels) %>%
  #   # dplyr::select(-group) %>%
  #   dplyr::distinct(., .keep_all = TRUE) %>%
  #   # dplyr::full_join(., usvi_sda_genera_shared_samples_baseline.df %>%
  #   #                    dplyr::filter(asv_id == shared_sda_genera_idx[[1]]) %>%
  #   #                    droplevels) %>%
  #   droplevels
  # print(ggplot(data = temp_list,
  #              aes(x = sample_id, y = relabund, fill = asv_id))
  #       + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
  #       + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
  #                  position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
  #       + theme_bw() 
  #       + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  #       + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
  #               strip.text.y = element_text(angle = 0))
  #       + facet_wrap(test_type~contrast))
  # 
  # print(ggplot(data = temp_list,
  #              aes(x = pair1, y = relabund, fill = asv_id))
  #       + geom_point(shape = 21)
  #       + theme_bw() 
  #       + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  #       + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
  #               strip.text.y = element_text(angle = 0))
  #       + facet_wrap(.~test_type))
  
  # temp_df2 <- usvi_sw_genus.tbl %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   t() %>%
  #   as.data.frame(.) %>%
  #   tibble::rownames_to_column(var = "sample_id") %>%
  #   tidyr::pivot_longer(., cols = !c("sample_id"),
  #                       names_to = "asv_id",
  #                       values_to = "relabund") %>%
  #   dplyr::filter(asv_id %in% selfish_sda_genera_idx) %>%
  #   dplyr::left_join(., usvi_sda_genera_compare.df %>%
  #                      dplyr::distinct(asv_id, site, sampling_time, .keep_all = FALSE) %>%
  #                      droplevels,
  #                    by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
  #   droplevels
  
  
  # temp_df1 <- usvi_sda_genera_compare.df %>%
  #   dplyr::distinct(asv_id, contrast, test_type, estimate, log2FoldChange_MMSE) %>%
  #   tidyr::unite("value", c(estimate, log2FoldChange_MMSE), remove = TRUE, na.rm = TRUE) %>%
  #   tidyr::nest(., .by = test_type, .key = "genera") %>%
  #   # dplyr::distinct(asv_id, contrast, test_type, pval, padj) %>%
  #   # droplevels %>%
  #   # tidyr::pivot_wider(., id_cols = "asv_id",
  #   #                    names_from = c(contrast),
  #   #                    values_from = c(test_type),
  #   #                    # names_from = c(test_type, contrast),
  #   #                    # values_from = c(pval, padj),
  #   #                    values_fill = NA) %>%
  #   droplevels
  
  
  
  
  
  # temp_df2 <- usvi_sda_genus_relabund.df %>%
  #   dplyr::filter(asv_id == shared_sda_genera_idx[1]) %>%
  #   droplevels
  
  # for(i in seq_len(5)){
  # for(i in seq_along(shared_sda_genera_idx)){
  #   namevar <- shared_sda_genera_idx[i]
  #   temp_df2 <- usvi_sda_genus_relabund.df %>%
  #     dplyr::filter(asv_id == namevar) %>%
  #     droplevels
  #   g <- print(
  #     ggplot(data = temp_df2)
  #     + theme_bw() 
  #     + geom_boxplot(aes(y = relabund, x = sampling_time, group = sampling_time), 
  #                    alpha = 0.6, show.legend = FALSE, color = "black",
  #                    position = position_dodge(0.8))
  #     + geom_point(aes(y = relabund, x = sampling_time, group = sampling_time, fill = sampling_time, shape = site),
  #                  alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
  #                  position = position_jitter(width = 0.2, seed = 48105))
  #     + scale_shape_manual(name = "Sampling site and time",
  #                          values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
  #     + scale_fill_manual(name = "Sampling time",
  #                         values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #     + scale_x_discrete(labels = sampling_time_lookup, name = "Sampling time")
  #     + scale_y_continuous(name = "Relative abundance (%)")
  #     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
  #                                  override.aes = list(color = "black", shape = 21, size = 3)),
  #              color = "none")
  #     # + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
  #     + theme(axis.text.x = element_blank(),
  #             axis.title.y = element_blank(),
  #             strip.text.y = element_text(angle = 0))
  #     + facet_grid(. ~ site,
  #                  drop = TRUE,
  #                  scales = "fixed", labeller = global_labeller)
  #     + ggtitle(paste0("Genera ", namevar))
  #   )
  #   assign(paste0("g_sda_", i), g, envir = .GlobalEnv, inherits = TRUE)
  #   rm(g)
  #   rm(namevar)
  # }
  
}
#181 ASVs were flagged in both methods, as SDA. how do they look?
# #the number of contrasts that an ASV was found to be significant, can vary between 1 and 9
# length(unique(usvi_sda_asvs_compare.df[["contrast"]]))
## [1] 9

temp_rankabund_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::arrange(desc(baseMean), .by_group = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, baseMean) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  # tibble::deframe(.)
  tibble::deframe(.) %>% quantile(., probs = seq(0.75, 0.25, -0.25), names = FALSE) %>% setNames(., c("top", "middle", "bottom"))

temp_list <- usvi_sda_asvs_compare.df %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::arrange(desc(baseMean), .by_group = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, baseMean) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::mutate(quantile_group = dplyr::case_when((baseMean >= temp_rankabund_idx[1]) ~ "top",
                                                  (baseMean < temp_rankabund_idx[1]) & (baseMean >= temp_rankabund_idx[2]) ~ "middle",
                                                  (baseMean < temp_rankabund_idx[2]) & (baseMean >= temp_rankabund_idx[3]) ~ "bottom",
                                                  (baseMean < temp_rankabund_idx[3]) ~ "rare",
                                                  .default = NA)) %>%
  dplyr::select(asv_id, quantile_group) %>%
  dplyr::mutate(quantile_group = factor(quantile_group, levels = c("top", "middle", "bottom", "rare"))) %>%
  dplyr::arrange(quantile_group) %>%
  dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  dplyr::right_join(., usvi_sda_asvs_compare.df %>%
                      dplyr::filter(!grepl("dispersion", model)) %>%
                        dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>%
                        dplyr::distinct(contrast, asv_id, .keep_all = FALSE) %>%
                        dplyr::group_by(asv_id) %>%
                        dplyr::summarise(num_results = length(contrast)) %>%
                        dplyr::arrange(num_results) %>%
                      droplevels, by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
  split(., f = .$quantile_group) %>%
  map(., ~.x %>%
        dplyr::select(asv_id, num_results) %>%
        droplevels)

temp_rankabund_idx <- temp_list %>%
  map(., ~quantile(.x$num_results, probs = seq(0.8, 0.2, -0.25), names = FALSE) %>% trunc(.) %>% unique(.))

shared_sda_asvs_all_relabund_list <- map(names(temp_list),
                                   ~temp_list[[.x]] %>%
                                     dplyr::mutate(grouping = dplyr::case_when((num_results >= temp_rankabund_idx[[.x]][1]) ~ "A",
                                                                                     (num_results < temp_rankabund_idx[[.x]][1]) & (num_results >= temp_rankabund_idx[[.x]][2])  ~ "B",
                                                                               (num_results < temp_rankabund_idx[[.x]][2]) ~ "C",
                                                                                     .default = NA)) %>%
                                     dplyr::select(asv_id, grouping) %>%
                                     dplyr::arrange(grouping) %>%
                                     dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
                                     droplevels) %>%
  setNames(., names(temp_list)) %>%
  map(., ~.x %>%
        dplyr::inner_join(., usvi_sda_asvs_compare.df, by = join_by(asv_id), multiple = "all", relationship = "many-to-many")) %>%
  map(., ~.x %>%
        dplyr::mutate(across(c(asv_id, grouping, contrast, hold, variable, model, test_type), ~factor(.x)))) %>%
  map(., ~.x %>% 
        dplyr::distinct(asv_id, grouping, contrast, hold, variable, pair1, pair2) %>%
        droplevels %>%
        split(., f = .$hold) %>%
        map(., ~.x %>%
              dplyr::mutate(sampling_time = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ hold,
                                                             .default = NA),
                            site= dplyr::case_when((hold %in% names(site_lookup)) ~ hold,
                                                   .default = NA)) %>%
              droplevels) %>%
        map(., ~.x %>%
              tidyr::pivot_longer(., cols = c(pair1, pair2),
                                  names_to = NULL,
                                  values_to = "value") %>%
              dplyr::mutate(sampling_time = dplyr::case_when(is.na(sampling_time) ~ value,
                                                             .default = sampling_time),
                            site = dplyr::case_when(is.na(site) ~ value,
                                                    .default = site)) %>%
              dplyr::select(-value) %>%
              dplyr::distinct(.)) %>%
        bind_rows(., .id = NULL) %>%
        dplyr::left_join(., usvi_sda_asvs_relabund.df, 
                         by = join_by(asv_id, site, sampling_time), relationship = "many-to-many", multiple = "all") %>%
        tidyr::drop_na(.) %>%
        dplyr::mutate(pair1 = gsub("([[:print:]]+)( - )(.*)$", "\\1", variable),
                      pair2 = gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)) %>%
        tidyr::pivot_longer(., cols = c(pair1, pair2),
                            names_to = NULL,
                            values_to = "pair") %>%
        dplyr::rowwise(.) %>%
        dplyr::filter((pair %in% sampling_time) | (pair %in% site)) %>%
        # dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
        dplyr::mutate(hold2 = dplyr::case_when((hold %in% names(site_lookup)) ~ "spatial", 
                                               .default = "temporal")) %>%
        dplyr::distinct(.) %>%
        droplevels) %>%
  map(., ~.x %>%
        dplyr::mutate(across(c(asv_id, grouping, contrast, hold, variable, sampling_time, site, pair, hold2, sample_id), ~factor(.x))) %>%
        dplyr::mutate(pair = factor(pair, levels = c(names(site_lookup), names(sampling_time_lookup))),
                      sampling_time = factor(sampling_time, levels = names(sampling_time_lookup)),
                      site = factor(site, levels = names(site_lookup))) %>%
        droplevels)

temp_list <- shared_sda_asvs_all_relabund_list %>%
  map(., ~.x %>%
        dplyr::left_join(., (usvi_sda_asvs_compare_summary.df %>%
                               dplyr::select(!c("test_type", c(Domain:Species))) %>%
                               tidyr::pivot_longer(., cols = !c(asv_id),
                                                   names_to = "contrast",
                                                   values_to = "significance") %>%
                               dplyr::mutate(significance = dplyr::case_when(grepl("_", significance) ~ "B",
                                                                             grepl("^deseq", significance) ~ "D",
                                                                             grepl("^corncob", significance) ~ "C",
                                                                             .default = NA)) %>%
                               # dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
                               tidyr::drop_na(.) %>%
                               droplevels),
                         by = join_by(asv_id, contrast)) %>%
        dplyr::mutate(label_y = ceiling(relabund)*1.1) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(desc(label_y)) %>%
        dplyr::select(asv_id, contrast, significance, label_y) %>% 
        dplyr::distinct(asv_id, contrast, significance, .keep_all = TRUE) %>% 
        tidyr::drop_na(.) %>%
        droplevels) %>%
  map(., ~.x %>%
        dplyr::mutate(across(c(asv_id, contrast, significance), ~factor(.x))) %>%
        droplevels)

shared_sda_asvs_all_relabund_list <- map(names(shared_sda_asvs_all_relabund_list),
                  ~temp_list[[.x]] %>%
                    dplyr::right_join(., shared_sda_asvs_all_relabund_list[[.x]] %>%
                                        # dplyr::filter(asv_id %in% namevar) %>%
                                        droplevels,
                                      by = join_by(asv_id, contrast), relationship = "one-to-many", multiple = "first") %>%
                    droplevels) %>%
  setNames(., names(shared_sda_asvs_all_relabund_list))

#filter for only those 181 ASVs flagged as SDA in both methods:

#old (lol 2024-12-05) method:
{
  # temp_df <- usvi_sda_asvs_compare_summary.df %>%
  #                                dplyr::select(!c("test_type", c(Domain:Species))) %>%
  #                                tidyr::pivot_longer(., cols = !c(asv_id),
  #                                                    names_to = "contrast",
  #                                                    values_to = "significance") %>%
  #                                dplyr::mutate(significance = dplyr::case_when(grepl("_", significance) ~ "B",
  #                                                                              grepl("^deseq", significance) ~ "D",
  #                                                                              grepl("^corncob", significance) ~ "C",
  #                                                                              .default = NA)) %>%
  #                                tidyr::drop_na(.) %>%
  #                                droplevels %>%
  #   dplyr::filter(grepl("B", significance)) %>%
  #         droplevels
  #
  # temp_df %>%
  #   dplyr::group_by(asv_id) %>%
  #   dplyr::summarise(num_results = length(contrast)) %>%
  #   dplyr::arrange(desc(num_results))
  
  # usvi_sda_asvs_both_list <- usvi_sda_asvs_compare_summary.df %>%
  #   dplyr::select(!c("test_type", c(Domain:Species))) %>%
  #   tidyr::pivot_longer(., cols = !c(asv_id),
  #                       names_to = "contrast",
  #                       values_to = "significance") %>%
  #   dplyr::mutate(significance = dplyr::case_when(grepl("_", significance) ~ "B",
  #                                                 grepl("^deseq", significance) ~ "D",
  #                                                 grepl("^corncob", significance) ~ "C",
  #                                                 .default = NA)) %>%
  #   tidyr::drop_na(.) %>%
  #   droplevels %>%
  #   dplyr::filter(grepl("B", significance)) %>%
  #   split(., f = .$contrast) %>%
  #   map(., ~.x %>%
  #         droplevels %>%
  #         dplyr::inner_join(., usvi_sda_asvs_compare.df, by = join_by(asv_id, contrast), multiple = "all", relationship = "many-to-many")) %>%
  #   map(., ~.x %>%
  #         dplyr::mutate(across(c(asv_id, contrast, hold, variable, model, test_type), ~factor(.x)))) %>%
  #   map(., ~.x %>% 
  #         dplyr::distinct(asv_id, contrast, hold, variable, pair1, pair2) %>%
  #         droplevels %>%
  #         split(., f = .$hold) %>%
  #         map(., ~.x %>%
  #               dplyr::mutate(sampling_time = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ hold,
  #                                                              .default = NA),
  #                             site= dplyr::case_when((hold %in% names(site_lookup)) ~ hold,
  #                                                    .default = NA)) %>%
  #               droplevels) %>%
  #         map(., ~.x %>%
  #               tidyr::pivot_longer(., cols = c(pair1, pair2),
  #                                   names_to = NULL,
  #                                   values_to = "value") %>%
  #               dplyr::mutate(sampling_time = dplyr::case_when(is.na(sampling_time) ~ value,
  #                                                              .default = sampling_time),
  #                             site = dplyr::case_when(is.na(site) ~ value,
  #                                                     .default = site)) %>%
  #               dplyr::select(-value) %>%
  #               dplyr::distinct(.)) %>%
  #         bind_rows(., .id = NULL) %>%
  #         dplyr::left_join(., usvi_sda_asvs_relabund.df, 
  #                          by = join_by(asv_id, site, sampling_time), relationship = "many-to-many", multiple = "all") %>%
  #         tidyr::drop_na(.) %>%
  #         dplyr::mutate(pair1 = gsub("([[:print:]]+)( - )(.*)$", "\\1", variable),
  #                       pair2 = gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)) %>%
  #         tidyr::pivot_longer(., cols = c(pair1, pair2),
  #                             names_to = NULL,
  #                             values_to = "pair") %>%
  #         dplyr::rowwise(.) %>%
  #         dplyr::filter((pair %in% sampling_time) | (pair %in% site)) %>%
  #         # dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
  #         dplyr::mutate(hold2 = dplyr::case_when((hold %in% names(site_lookup)) ~ "spatial", 
  #                                                .default = "temporal")) %>%
  #         dplyr::distinct(.) %>%
  #         droplevels) %>%
  #   map(., ~.x %>%
  #         dplyr::mutate(across(c(asv_id, contrast, hold, variable, sampling_time, site, pair, hold2, sample_id), ~factor(.x))) %>%
  #         dplyr::mutate(pair = factor(pair, levels = c(names(site_lookup), names(sampling_time_lookup))),
  #                       sampling_time = factor(sampling_time, levels = names(sampling_time_lookup)),
  #                       site = factor(site, levels = names(site_lookup))) %>%
  #         dplyr::mutate(group_label = dplyr::case_when(grepl("dawn", sampling_time) ~ "_d",
  #                                                 grepl("photo", sampling_time) ~ "_p",
  #                                                .default = NA) %>%
  #                         paste0(site, .)) %>%
  #         dplyr::mutate(group_label = factor(group_label, levels = names(group_labels_lookup))) %>%
  #         droplevels)
  }

temp_df <- shared_sda_asvs_all_relabund_list %>%
  map(., ~.x %>%
        dplyr::mutate(across(!"relabund", ~factor(.x))) %>%
        dplyr::arrange(asv_id, contrast, significance) %>%
        dplyr::distinct(asv_id, contrast, significance)) %>%
  bind_rows(., .id = NULL) %>%
  dplyr::filter(grepl("B", significance)) %>%
  droplevels


shared_sda_asvs_idx_filtered_list <- shared_sda_asvs_all_relabund_list %>%
  map(., ~.x %>%
        dplyr::select(asv_id, contrast, site, sampling_time, sample_id, relabund) %>%
        dplyr::mutate(across(!"relabund", ~factor(.x))) %>%
        ungroup %>%      
        dplyr::left_join(., temp_df, by = join_by(asv_id, contrast), relationship = "many-to-many", multiple = "all") %>%
        dplyr::filter(grepl("B", significance)) %>%
        droplevels) %>%
# temp_list2 <- shared_sda_asvs_idx_filtered_list %>%
  bind_rows(., .id = NULL) %>%
  dplyr::distinct(asv_id, contrast, significance, .keep_all = FALSE) %>%
  droplevels

shared_sda_asvs_idx_filtered_list <- shared_sda_asvs_idx_filtered_list %>%
  dplyr::left_join(., bind_rows(shared_sda_asvs_all_relabund_list, .id = "quantile_group") %>%
                     dplyr::select(-significance)) %>%
  dplyr::mutate(group_label = dplyr::case_when(grepl("dawn", sampling_time) ~ "_d",
                                               grepl("photo", sampling_time) ~ "_p",
                                               .default = NA) %>%
                  paste0(site, .)) %>%
  dplyr::mutate(group_label = factor(group_label, levels = names(group_labels_lookup))) %>%
  dplyr::mutate(quantile_group = factor(quantile_group, levels = c("top", "middle", "bottom", "rare"))) %>%
  split(., f = .$contrast) %>%
  map(., ~.x %>%
        droplevels)

shared_sda_asvs_idx_filtered_quantile_list <- shared_sda_asvs_idx_filtered_list %>%
  bind_rows(., .id = "contrast") %>%
  split(., f = .$quantile_group) %>%
  map(., ~.x %>%
        droplevels)


#make a summary table of the range of relative abundances of SDA ASVs
usvi_sda_asvs_shared_summary.df <- usvi_sda_asvs_compare_summary.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>%
  droplevels %>%
  dplyr::select(-c(test_type:Species)) %>%
  tidyr::pivot_longer(cols = !c("asv_id"),
                      names_to = "contrast",
                      values_to = "sig") %>%
  dplyr::mutate(contrast = factor(contrast, levels = unique(.[["contrast"]]))) %>%
  dplyr::mutate(obs = dplyr::case_when((sig == "deseq_corncob") ~ 1,
                                       .default = 0)) %>%
  dplyr::group_by(contrast, sig) %>%
  dplyr::summarise(num_SDA = sum(obs), .groups = "keep") %>%
  dplyr::group_by(contrast) %>%
  tidyr::complete(fill = list(sig= "deseq_corncob", num_sda = 0)) %>%
  dplyr::filter(grepl("deseq_corncob", sig)) %>%
  dplyr::group_by(contrast, sig) %>%
  dplyr::summarise(num_SDA = sum(num_SDA), .groups = "keep") %>%
  dplyr::mutate(hold = dplyr::case_when(grepl("^dawn|^peak", contrast) ~ "temporal",
                                        .default = "spatial")) %>%
  dplyr::mutate(variable = dplyr::case_when(hold == "temporal" ~ "spatial",
                                            .default = "temporal")) %>%
  dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
  dplyr::relocate(hold, variable) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(-sig) %>%
  droplevels

usvi_sda_asvs_shared_relabund.df <- shared_sda_asvs_idx_filtered_list %>%
  # usvi_sda_asvs_shared_relabund.df <- shared_sda_asvs_all_relabund_list %>%
  map(., ~.x %>%
        dplyr::select(asv_id, contrast, site, sampling_time, sample_id, relabund) %>%
        dplyr::mutate(across(!"relabund", ~factor(.x))) %>%
        ungroup %>%      
        dplyr::left_join(., temp_df, by = join_by(asv_id, contrast), relationship = "many-to-many", multiple = "all") %>%
        droplevels %>%
        dplyr::group_by(asv_id, contrast, significance, site, sampling_time) %>%
        dplyr::summarise(mean_relabund = mean(relabund, na.rm = TRUE),
                         sd_relabund = sd(relabund,  na.rm = TRUE),
                         .groups = "keep") %>%
        dplyr::filter(grepl("B", significance)) %>%
        droplevels) %>% bind_rows(., .id = NULL) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, contrast, site, sampling_time, contains("relabund")) %>%
  droplevels


usvi_sda_asvs_shared_summary.df <- usvi_sda_asvs_shared_summary.df %>%
  dplyr::left_join(., (usvi_sda_asvs_shared_relabund.df %>%
                         dplyr::group_by(contrast) %>%
                         # dplyr::summarise(num_SDA_ASVs = length(unique(asv_id))) %>%
                         dplyr::summarise(across(contains("mean_relabund"), list(min = min, median = median, mean = mean, max = max))) %>%
                         dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
                         droplevels), by = join_by(contrast))
  
if(!file.exists(paste0(projectpath, "/", "usvi_shared_sda_asvs_filtered_relabund-.*.RData"))){
  save(list = c("shared_sda_asvs_idx_filtered_list", "shared_sda_asvs_all_relabund_list"),
       file = paste0(projectpath, "/", "usvi_shared_sda_asvs_filtered_relabund-", Sys.Date(), ".RData"))
}

if(!any(grepl("sda_asvs_shared_summary", list.files(projectpath, pattern = "usvi_.*.tsv")))){
  readr::write_delim(usvi_sda_asvs_shared_summary.df, paste0(projectpath, "/", "usvi_sda_asvs_shared_summary-", Sys.Date(), ".tsv"),
                     delim = "\t", col_names = TRUE, na = "")
}


# Plot the relabund of ASVs -----------------------------------------------

#testing it out:
{
  
# ## temp_df1 <- shared_sda_asvs_idx_list[[1]] %>%
# temp_df1 <- shared_sda_asvs_idx_filtered_list[[1]] %>%
#   # dplyr::filter(grepl("A", grouping)) %>%
#   droplevels
# 
# temp_g1 <- (
#   ggplot(data = temp_df1 %>%
#            droplevels, aes(y = relabund, x = asv_id, group = interaction(asv_id,group_label)))
#   + theme_bw()
#   + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
#                 alpha = 1, show.legend = FALSE, color = "grey",
#                 position = position_dodge(0.8))
#   # + geom_point(aes(fill = sampling_time, shape = site),
#   + geom_point(aes(fill = group_label, shape = site),
#                alpha = 0.8, color = "black", size = 1, show.legend = TRUE,
#                position = position_jitter(width = 0.2, seed = 48105))
#   # + geom_text(aes(label = significance, x = 1.5, vjust = "outward", hjust = "center", group = contrast,
#   #                 y = label_y),
#   #             colour = "black", fontface = "bold")
#   + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#   # + scale_fill_manual(name = "Sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_discrete_manual(aesthetics = c("fill"), 
#                           values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
#   # + scale_x_discrete(labels = c(site_lookup, sampling_time_lookup), name = NULL)
#   + scale_x_discrete(labels = usvi_genera_relabel, name = NULL)
#   + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.5)))
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical",
#                                override.aes = list(color = "black", shape = 21, size = 3)),
#            shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical",
#                                 override.aes = list(color = "black", size = 3)),
#            color = "none")
#   + theme(strip.text.y = element_text(size = rel(0.7), angle = 0),
#           strip.text.x = element_text(size = rel(0.7), angle = 0),
#           axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 3),
#           axis.title.y = element_blank())
#   + facet_wrap(. ~ quantile_group,
#                drop = TRUE,
#                labeller = global_labeller,
#                scales = "free")
#   # + ggtitle(paste0("Highest recurrence across contrasts of observed SDA ASVs in ", namevar, " relative abundance ranks"))
# )
# temp_g1
# 
# temp_df2 <- shared_sda_asvs_idx_filtered_quantile_list[[1]] %>%
#   dplyr::filter(grepl("A", grouping)) %>%
#   droplevels
# 
# temp_g2 <- (
#   ggplot(data = temp_df2 %>%
#            droplevels, aes(y = relabund, x = pair, group = interaction(pair,group_label)))
#   + theme_bw()
#   + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
#                 alpha = 1, show.legend = FALSE, color = "grey",
#                 position = position_dodge(0.8))
#   # + geom_point(aes(fill = sampling_time, shape = site),
#   + geom_point(aes(fill = group_label, shape = site),
#                alpha = 0.8, color = "black", size = 1, show.legend = TRUE,
#                position = position_jitter(width = 0.2, seed = 48105))
#   # + geom_text(aes(label = significance, x = 1.5, vjust = "outward", hjust = "center", group = contrast,
#   #                 y = label_y),
#   #             colour = "black", fontface = "bold")
#   + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#   # + scale_fill_manual(name = "Sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_discrete_manual(aesthetics = c("fill"),
#                           values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
#   # + scale_x_discrete(labels = c(site_lookup, sampling_time_lookup), name = NULL)
#   + scale_x_discrete(labels = usvi_genera_relabel, name = NULL)
#   + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.5)))
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical",
#                                override.aes = list(color = "black", shape = 21, size = 3)),
#            shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical",
#                                 override.aes = list(color = "black", size = 3)),
#            color = "none")
#   + theme(strip.text.y = element_text(size = rel(0.7), angle = 0),
#           strip.text.x = element_text(size = rel(0.7), angle = 0),
#           axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 3),
#           axis.title.y = element_blank())
#   + facet_grid(asv_id ~ contrast,
#                drop = TRUE,
#                labeller = global_labeller,
#                scales = "free")
#   # + ggtitle(paste0("Highest recurrence across contrasts of observed SDA ASVs in ", namevar, " relative abundance ranks"))
# )
# temp_g2
  
}

#plotting ASVs, grouped by the frequency they were found SDA in the 9 contrasts, then grouped 1-4 by their mean abundance in the dataset
rm(list = apropos("g2_sda_.*", mode = "list"))
# for(i in c(1)){
for(i in seq_along(shared_sda_asvs_idx_filtered_quantile_list)){
  
  namevar <- names(shared_sda_asvs_idx_filtered_quantile_list[i])
  temp_df1 <- shared_sda_asvs_idx_filtered_quantile_list[[i]] %>%
    # dplyr::filter(grepl("top", quantile_group)) %>%
    dplyr::filter(grepl("A", grouping)) %>%
    droplevels
  temp_df2 <- shared_sda_asvs_idx_filtered_quantile_list[[i]] %>%
    # dplyr::filter(grepl("middle", quantile_group)) %>%
    dplyr::filter(grepl("B", grouping)) %>%
    droplevels
  temp_df3 <- shared_sda_asvs_idx_filtered_quantile_list[[i]] %>%
    # dplyr::filter(grepl("bottom", quantile_group)) %>%
    dplyr::filter(grepl("C", grouping)) %>%
    droplevels
  # temp_df4 <- shared_sda_asvs_idx_filtered_list[[i]] %>%
  #   dplyr::filter(grepl("rare", quantile_group)) %>%
  #   # dplyr::filter(grepl("C", grouping)) %>%
  #   droplevels

  temp_g1 <- (
    ggplot(data = temp_df1 %>% droplevels, 
           aes(y = relabund, x = pair, group = interaction(pair,group_label)))
           # aes(y = relabund, x = pair, group = pair))
    + theme_bw() 
    + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
                  alpha = 1, show.legend = FALSE, color = "grey",
                  position = position_dodge(0.8))
    + geom_point(aes(fill = group_label, shape = site),
    # + geom_point(aes(fill = sampling_time, shape = site), 
                 alpha = 0.8, color = "black", size = 1, show.legend = TRUE, 
                 position = position_jitter(width = 0.2, seed = 48105))
    # + geom_text(aes(label = significance, x = 1.5, vjust = "outward", hjust = "center", group = contrast,
    #                 y = label_y), colour = "black", fontface = "bold")
    + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
    # + scale_fill_manual(name = "Sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
    # + scale_x_discrete(labels = c(site_lookup, sampling_time_lookup), name = NULL)
    + scale_discrete_manual(aesthetics = c("fill"), 
                            values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
    + scale_x_discrete(labels = usvi_genera_relabel, name = NULL)
    + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.5)))
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
                                  override.aes = list(color = "black", size = 3)),
             color = "none")
    + theme(strip.text.y = element_text(size = rel(0.7), angle = 0), 
            strip.text.x = element_text(size = rel(0.7), angle = 0), 
            axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 3),
            axis.title.y = element_blank())
    + facet_grid(asv_id ~ contrast,
                 drop = TRUE,
                 labeller = global_labeller,
                 scales = "free")
    + ggtitle(paste0("Highest recurrence across contrasts of observed SDA ASVs in ", namevar, " relative abundance ranks"))
  )
  temp_g2 <- (
    ggplot(data = temp_df2 %>% droplevels, 
           aes(y = relabund, x = pair, group = interaction(pair,group_label)))
           # aes(y = relabund, x = pair, group = pair))
    + theme_bw() 
    + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
                  alpha = 1, show.legend = FALSE, color = "grey",
                  position = position_dodge(0.8))
    + geom_point(aes(fill = group_label, shape = site),
                 # + geom_point(aes(fill = sampling_time, shape = site), 
                 alpha = 0.8, color = "black", size = 1, show.legend = TRUE, 
                 position = position_jitter(width = 0.2, seed = 48105))
    # + geom_text(aes(label = significance, x = 1.5, vjust = "outward", hjust = "center", group = contrast,
    #                 y = label_y), colour = "black", fontface = "bold")
    + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
    # + scale_fill_manual(name = "Sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
    # + scale_x_discrete(labels = c(site_lookup, sampling_time_lookup), name = NULL)
    + scale_discrete_manual(aesthetics = c("fill"), 
                            values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
    + scale_x_discrete(labels = usvi_genera_relabel, name = NULL)
    + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.5)))
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
                                  override.aes = list(color = "black", size = 3)),
             color = "none")
    + theme(strip.text.y = element_text(size = rel(0.7), angle = 0), 
            strip.text.x = element_text(size = rel(0.7), angle = 0), 
            axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 3),
            axis.title.y = element_blank())
    + facet_grid(asv_id ~ contrast,
                 drop = TRUE,
                 labeller = global_labeller,
                 scales = "free")
    + ggtitle(paste0("Mid recurrence across contrasts of observed SDA ASVs in ", namevar, " relative abundance ranks"))
  )
  temp_g3 <- (
    ggplot(data = temp_df3 %>% droplevels, 
           aes(y = relabund, x = pair, group = interaction(pair,group_label)))
           # aes(y = relabund, x = pair, group = pair))
    + theme_bw() 
    + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
                  alpha = 1, show.legend = FALSE, color = "grey",
                  position = position_dodge(0.8))
    + geom_point(aes(fill = group_label, shape = site),
                 # + geom_point(aes(fill = sampling_time, shape = site), 
                 alpha = 0.8, color = "black", size = 1, show.legend = TRUE, 
                 position = position_jitter(width = 0.2, seed = 48105))
    # + geom_text(aes(label = significance, x = 1.5, vjust = "outward", hjust = "center", group = contrast,
    #                 y = label_y), colour = "black", fontface = "bold")
    + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
    # + scale_fill_manual(name = "Sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
    # + scale_x_discrete(labels = c(site_lookup, sampling_time_lookup), name = NULL)
    + scale_discrete_manual(aesthetics = c("fill"), 
                            values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
    + scale_x_discrete(labels = usvi_genera_relabel, name = NULL)
    + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.5)))
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
                                  override.aes = list(color = "black", size = 3)),
             color = "none")
    + theme(strip.text.y = element_text(size = rel(0.7), angle = 0), 
            strip.text.x = element_text(size = rel(0.7), angle = 0), 
            axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 3),
            axis.title.y = element_blank())
    + facet_grid(asv_id ~ contrast, 
                 drop = TRUE,
                 labeller = global_labeller,
                 scales = "free")
    + ggtitle(paste0("Low recurrence across contrasts of observed SDA ASVs in ", namevar, " relative abundance ranks"))
  )
  assign(paste0("g2_sda_", i, "_", namevar, "_Aubiquitous"), temp_g1, envir = .GlobalEnv, inherits = TRUE)
  assign(paste0("g2_sda_", i, "_", namevar, "_Bfrequent"), temp_g2, envir = .GlobalEnv, inherits = TRUE)
  assign(paste0("g2_sda_", i, "_", namevar, "_Cnotfrequent"), temp_g3, envir = .GlobalEnv, inherits = TRUE)
  rm(temp_g1)
  rm(temp_g2)
  rm(temp_g3)
  rm(namevar)
}


gpatch <- apropos("g2_sda_1.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`)
g2_sda_asvs_g1 <- (gpatch & theme(legend.position = "none")) + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")

gpatch <- apropos("g2_sda_2.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`)
g2_sda_asvs_g2 <- (gpatch & theme(legend.position = "none")) + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")

gpatch <- apropos("g2_sda_3.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`)
g2_sda_asvs_g3 <- (gpatch & theme(legend.position = "none")) + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")

gpatch <- apropos("g2_sda_4.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`)
g2_sda_asvs_g4 <- (gpatch & theme(legend.position = "none")) + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")


if(!any(grepl("asvs_relabund_(A|B|C|D)", list.files(projectpath, pattern = "usvi_sda_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund_Atop-", Sys.Date(), ".png"),
         g2_sda_asvs_g1,
         width = 24, height = 15, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund_Bmid-", Sys.Date(), ".png"),
         g2_sda_asvs_g2,
         width = 24, height = 15, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund_Clow-", Sys.Date(), ".png"),
         g2_sda_asvs_g3,
         width = 24, height = 15, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund_Drare-", Sys.Date(), ".png"),
         g2_sda_asvs_g4,
         width = 24, height = 15, units = "in")
}

#old method:
{

# rm(list = apropos("g_sda_.*", mode = "list"))
# # for(i in c(1)){
# for(i in seq_along(shared_sda_asvs_idx_cluster)){
#   namevar <- temp_list[[i]] %>%
#     dplyr::filter(grepl("A", grouping)) %>%
#     purrr::pluck("asv_id")
#   temp_df2 <- usvi_sda_asvs_relabund.df %>%
#     dplyr::filter(asv_id %in% namevar) %>%
#     droplevels
#   g <- print(
#     ggplot(data = temp_df2)
#     + theme_bw() 
#     # + geom_boxplot(aes(y = relabund, x = sampling_time, group = sampling_time), 
#     #                alpha = 0.6, show.legend = FALSE, color = "black", outliers = TRUE, outlier.shape = NA,
#     #                position = position_dodge(0.8))
#     + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
#                   aes(y = relabund, x = sampling_time, group = sampling_time), 
#                   alpha = 0.6, show.legend = FALSE, color = "grey", 
#                   position = position_dodge(0.8))
#     + geom_point(aes(y = relabund, x = sampling_time, group = sampling_time, fill = sampling_time, shape = site),
#                  alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
#                  position = position_jitter(width = 0.2, seed = 48105))
#     + scale_shape_manual(name = "Sampling site and time",
#                          values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#     + scale_fill_manual(name = "Sampling time",
#                         values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#     + scale_x_discrete(labels = sampling_time_lookup, name = "Sampling time")
#     + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.1)))
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                  override.aes = list(color = "black", shape = 21, size = 3)),
#              color = "none")
#     # + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
#     + theme(axis.text.x = element_blank(),
#             strip.text.y = element_text(size = rel(0.7)),
#             axis.title.y = element_blank())
#     + facet_grid(asv_id ~ site,
#                  drop = TRUE,
#                  # labeller = labeller(site = site_lookup, asv_id = stringr::str_wrap(usvi_genera_relabel)),
#                  labeller = global_labeller,
#                  scales = "free")
#     + ggtitle(paste0("Group ", i))
#   )
#   assign(paste0("g_sda_", i), g, envir = .GlobalEnv, inherits = TRUE)
#   rm(g)
#   rm(namevar)
# }

# gpatch <- apropos("g_sda_.*", mode = "list") %>%
#   lapply(., get) %>%
#   # purrr::reduce(., `+ theme(legend.position = "none")`) %>%
#   purrr::reduce(., `+`)
# g2_sda_asvs <- gpatch + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")
# # g2_sda_asvs
# 
# if(!any(grepl("asvs_relabund", list.files(projectpath, pattern = "usvi_sda_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund-", Sys.Date(), ".png"),
#          g2_sda_asvs,
#          width = 15, height = 24, units = "in")
# }
# 
# 
# for(i in seq_along(shared_sda_asvs_idx_cluster)){
#   namevar <- paste0("g_sda_", i)
#   temp_g2_sda_asvs <- get(namevar, inherits = TRUE)
#   temp_g2_sda_asvs <- temp_g2_sda_asvs + theme(legend.position = "none")
#   ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund_group", i, "-", Sys.Date(), ".png"),
#          temp_g2_sda_asvs,
#          width = 8, height = 16, units = "in")
# }
  
}

# Plot mean relative abundance of SDA ASVs --------------------------------

#Instead of plotting 25 points per group, calculate the mean relative abundance of that ASV in that group
#plot point+sterr bars for each SDA ASV 

usvi_sda_asvs_both_list <- shared_sda_asvs_idx_filtered_list
usvi_sda_asvs_both_mean_list <- usvi_sda_asvs_both_list %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::group_by(asv_id, contrast, hold, variable, sampling_time, site, pair, hold2, group_label) %>%
        dplyr::summarise(mean = mean(relabund, na.rm = TRUE),
                         sd = sd(relabund, na.rm = TRUE), .groups = "keep") %>%
        dplyr::left_join(., shared_sda_asvs_idx_cluster %>%
                           bind_rows(., .id = "rank") %>%
                           dplyr::select(asv_id, rank) %>%
                           dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
                           dplyr::mutate(rank = factor(rank, levels = c("top", "middle", "bottom", "rare"))), by = join_by(asv_id)) %>%
        droplevels)

#arrange the asvs by whether they are more abudant in first or second part of variable
abundrank_labeller <- (c(`0` = "Abundant and SDA in ",
                                    `1` = "Overall rare and SDA in "))
#test it out:
{

  temp_df1 <- usvi_sda_asvs_both_mean_list[[5]]
  temp_df1 <- temp_df1 %>%
    dplyr::left_join(., (temp_df1 %>%
                           dplyr::ungroup(.) %>%
                           dplyr::select(asv_id, pair, mean) %>%
                           dplyr::group_by(asv_id) %>%
                           dplyr::slice_max(mean, na_rm = TRUE) %>%
                           dplyr::mutate(enriched = pair) %>%
                           dplyr::ungroup(.) %>%
                           dplyr::mutate(abundrank = dplyr::case_when(mean <= mean(.[["mean"]]) ~ 1, .default = 0)) %>%
                           dplyr::select(asv_id, enriched, abundrank)), by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
    dplyr::arrange(abundrank, enriched, asv_id) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(grouping = recode(enriched, !!!c(group_labels_lookup, sampling_time_lookup, site_lookup))) %>%
    dplyr::mutate(grouping = recode(abundrank, !!!abundrank_labeller) %>%
                    paste0(., grouping)) %>%
    dplyr::arrange(enriched, abundrank, asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]])))

# 
#   temp_g1 <- (
#     ggplot(data = temp_df1 %>%
#              dplyr::filter(abundrank == 0) %>%
#              # dplyr::filter(enriched == unique(.[["enriched"]])[1]) %>%
#              droplevels,
#            aes(x = rev(asv_id), y = mean, group = interaction(asv_id, group_label)))
#     + theme_bw()
#     + geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd), color = group_label), width = 0.3,
#                     position = position_dodge(width = 0.4, preserve = "total"),
#                     show.legend = FALSE)
#     + geom_point(aes(fill = group_label, shape = site), na.rm = TRUE,
#                  position = position_dodge(width = 0.4, preserve = "total"),
#                  alpha = 1, color = "black", size = 2, show.legend = TRUE)
#     + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#     + scale_discrete_manual(aesthetics = c("fill", "color"),
#                             values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
#     # + scale_x_discrete(labels = usvi_genera_relabel,
#     + scale_x_discrete(labels = NULL,
#                        expand = expansion(mult = c(0.1)),
#                        name = NULL)
#     + scale_y_continuous(name = "Relative abundance (%)",
#                          expand = expansion(mult = c(0.01,0.1)))
#     + geom_hline(aes(yintercept = 0), color = "black")
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Variables",  direction = "vertical",
#                                  override.aes = list(color = "black", shape = 21, size = 3)),
#              shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical",
#                                   override.aes = list(color = "black", size = 3)),
#              color = "none")
#     + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 5),
#             # strip.text.y = element_text(size = rel(0.8), angle = 0),  #facet_wrap
#             strip.text.y = element_blank(),
#             strip.background = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
#             panel.grid.major.y = element_line(linetype = "solid", color = "grey80", linewidth = 0.5),
#             # axis.text.y = element_text(angle = -90),
#             # strip.text.y.right = element_blank(),
#             # strip.text.x = element_text(size = rel(0.8), angle = 0), strip.text.y= element_blank(), #facet_wrap
#             # strip.text.y = element_text(size = rel(0.8), angle = -90), strip.text.x = element_blank(), #facet_grid
#             legend.position = "none")
#     # + facet_wrap(.~abundrank, ncol = 1, scales = "free", axes = "all", strip.position = "right",
#     # + facet_grid(row = vars(enriched), cols = vars(abundrank), scales = "free", space = "free", axes = "all",
#     + facet_grid(row = vars(grouping), scales = "free", space = "free", axes = "all",
#     # + facet_grid(row = vars(abundrank), scales = "free", space = "free", axes = "all",
#                  drop = TRUE, shrink = FALSE,
#                  labeller = global_labeller)
#     + coord_flip()
#     # + ggtitle(paste0("Observed SDA ASVs in ", namevar))
#   )
#   temp_g1
#   temp_g2 <- (
#     ggplot(data = temp_df1 %>%
#              dplyr::filter(abundrank == 1) %>%
#              # dplyr::filter(enriched == unique(.[["enriched"]])[2]) %>%
#              droplevels,
#            aes(x = rev(asv_id), y = mean, group = interaction(asv_id, group_label)))
#     + theme_bw()
#     + geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd), color = group_label), width = 0.3,
#                     position = position_dodge(width = 0.4, preserve = "total"),
#                     show.legend = FALSE)
#     + geom_point(aes(fill = group_label, shape = site), na.rm = TRUE,
#                  position = position_dodge(width = 0.4, preserve = "total"),
#                  alpha = 1, color = "black", size = 2, show.legend = TRUE)
#     + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#     + scale_discrete_manual(aesthetics = c("fill", "color"),
#                             values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
#     # + scale_x_discrete(labels = usvi_genera_relabel,
#     + scale_x_discrete(labels = NULL,
#                        expand = expansion(mult = c(0.1)),
#                        name = NULL)
#     + scale_y_continuous(name = "Relative abundance (%)",
#                          expand = expansion(mult = c(0.01,0.1)))
#     + geom_hline(aes(yintercept = 0), color = "black")
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Variables",  direction = "vertical",
#                                  override.aes = list(color = "black", shape = 21, size = 3)),
#              shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical",
#                                   override.aes = list(color = "black", size = 3)),
#              color = "none")
#     + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 5),
#             # strip.text.y = element_text(size = rel(0.8), angle = 0),  #facet_wrap
#             strip.text.y = element_blank(),
#             strip.background = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
#             panel.grid.major.y = element_line(linetype = "solid", color = "grey80", linewidth = 0.5),
#             # axis.text.y = element_text(angle = -90),
#             # strip.text.y.right = element_blank(),
#             # strip.text.x = element_text(size = rel(0.8), angle = 0), strip.text.y= element_blank(), #facet_wrap
#             # strip.text.y = element_text(size = rel(0.8), angle = -90), strip.text.x = element_blank(), #facet_grid
#             legend.position = "none")
#     # + facet_wrap(.~abundrank, ncol = 2, scales = "free", axes = "all", strip.position = "right",
#     # + facet_grid(row = vars(enriched), cols = vars(abundrank), scales = "free", space = "free", axes = "all",
#     + facet_grid(row = vars(grouping), scales = "free", space = "free", axes = "all",
#     # + facet_grid(col = vars(abundrank), scales = "free", space = "free", axes = "all",
#                  drop = TRUE, shrink = FALSE,
#                  labeller = global_labeller)
#     + coord_flip()
#     # + ggtitle(paste0("Observed SDA ASVs in ", namevar))
#   )
#   temp_g2
# temp_g1 | temp_g2
  }

rm(list = apropos("g3_sda_.*", mode = "list"))
# for(i in c(1)){
for(i in seq_along(usvi_sda_asvs_both_mean_list)){
  
  namevar <- names(usvi_sda_asvs_both_mean_list[i]) %>%
    gsub("LB_seagrass", "Lameshur Bay", .) %>%
    gsub("peak_photo", "afternoon", .)
  temp_df1 <- usvi_sda_asvs_both_mean_list[[i]] %>%
    droplevels
  temp_df1 <- temp_df1 %>%
    dplyr::left_join(., (temp_df1 %>%
                           dplyr::ungroup(.) %>%
                           dplyr::select(asv_id, pair, mean) %>%
                           dplyr::group_by(asv_id) %>%
                           dplyr::slice_max(mean, na_rm = TRUE) %>%
                           dplyr::mutate(enriched = pair) %>%
                           dplyr::ungroup(.) %>%
                           dplyr::mutate(abundrank = dplyr::case_when(mean <= mean(.[["mean"]]) ~ 1, .default = 0)) %>%
                           dplyr::select(asv_id, enriched, abundrank)), by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
    dplyr::arrange(abundrank, enriched, asv_id) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(grouping = recode(enriched, !!!c(group_labels_lookup, sampling_time_lookup, site_lookup))) %>%
    dplyr::mutate(grouping = recode(abundrank, !!!abundrank_labeller) %>%
                    paste0(., grouping)) %>%
    dplyr::arrange(enriched, abundrank, asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]])))
  
  namevar2 <- paste0(unique(temp_df1$hold2), "_", unique(temp_df1$hold))
  if(nrow(temp_df1) > 1){
    temp_g1 <- (
      ggplot(data = temp_df1 %>%
               droplevels,
             aes(x = rev(asv_id), y = mean, group = interaction(asv_id, group_label)))
      + theme_bw() 
      + geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd), color = group_label), width = 0.3,
                      position = position_dodge(width = 0.4, preserve = "total"),
                      show.legend = FALSE)
      + geom_point(aes(fill = group_label, shape = site), na.rm = TRUE,
                   position = position_dodge(width = 0.4, preserve = "total"),
                   alpha = 1, color = "black", size = 2, show.legend = TRUE)
      + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
      + scale_discrete_manual(aesthetics = c("fill", "color"), 
                              values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
      + scale_x_discrete(labels = usvi_genera_relabel, expand = expansion(add = 1), name = NULL)
      + scale_y_continuous(name = "Relative abundance (%)", labels = scaleFUN2,
                           expand = expansion(mult = c(0.01,0.1)))
      + geom_hline(aes(yintercept = 0), color = "black")
      + guides(fill = guide_legend(order = 1, ncol = 1, title = "Variables",  direction = "vertical", 
                                   override.aes = list(color = "black", shape = 21, size = 3)),
               shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
                                    override.aes = list(color = "black", size = 3)),
               color = "none")
      + theme(strip.text.y = element_text(size = rel(0.8), angle = 0),
              strip.background = element_blank(),
              legend.position = "none",
              panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
              panel.grid.major.y = element_line(linetype = "solid", color = "grey80", linewidth = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = rel(0.8)),
              axis.text.y = element_text(angle = 0, size = rel(0.8)))
      # + facet_wrap(vars(grouping), ncol = 1,
      # + facet_wrap(.~enriched, ncol = 2,
      + facet_grid(rows = vars(grouping), cols = NULL, space = "free", axes = "all",
                   drop = TRUE, shrink = TRUE,
                   labeller = global_labeller,
                   scales = "free")
      + coord_flip()
      + ggtitle(paste0(namevar))
    )
    assign(paste0("g3_sda_", i, "_", namevar2), temp_g1, envir = .GlobalEnv, inherits = TRUE)
  rm(temp_g1)
  }
  rm(namevar)
  rm(temp_df1)
}

rm(list = apropos("g4_sda_.*", mode = "list"))
# for(i in c(1)){
for(i in seq_along(usvi_sda_asvs_both_mean_list)){
  
  namevar <- names(usvi_sda_asvs_both_mean_list[i]) %>%
    gsub("LB_seagrass", "Lameshur Bay", .) %>%
    gsub("peak_photo", "afternoon", .)
  temp_df1 <- usvi_sda_asvs_both_mean_list[[i]] %>%
    droplevels
  temp_df1 <- temp_df1 %>%
    dplyr::left_join(., (temp_df1 %>%
                           dplyr::ungroup(.) %>%
                           dplyr::select(asv_id, pair, mean) %>%
                           dplyr::group_by(asv_id) %>%
                           dplyr::slice_max(mean, na_rm = TRUE) %>%
                           dplyr::mutate(enriched = pair) %>%
                           dplyr::ungroup(.) %>%
                           dplyr::mutate(abundrank = dplyr::case_when(mean <= mean(.[["mean"]]) ~ 1, .default = 0)) %>%
                           dplyr::select(asv_id, enriched, abundrank)), by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
    dplyr::arrange(abundrank, enriched, asv_id) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(grouping = recode(enriched, !!!c(group_labels_lookup, sampling_time_lookup, site_lookup))) %>%
    dplyr::mutate(grouping = recode(abundrank, !!!abundrank_labeller) %>%
                    paste0(., grouping)) %>%
    dplyr::arrange(enriched, abundrank, asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]])))
  
  namevar2 <- paste0(unique(temp_df1$hold2), "_", unique(temp_df1$hold))
  if(nrow(temp_df1) > 1){
    temp_g1 <- (
      ggplot(data = temp_df1 %>%
               dplyr::filter(abundrank == 0) %>%
               droplevels,
             aes(x = rev(asv_id), y = mean, group = interaction(asv_id, group_label)))
      + theme_bw() 
      + geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd), color = group_label), width = 0.3,
                      position = position_dodge(width = 0.4, preserve = "total"),
                      show.legend = FALSE)
      + geom_point(aes(fill = group_label, shape = site), na.rm = TRUE,
                   position = position_dodge(width = 0.4, preserve = "total"),
                   alpha = 1, color = "black", size = 2, show.legend = TRUE)
      + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
      + scale_discrete_manual(aesthetics = c("fill", "color"), 
                              values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
      + scale_x_discrete(labels = usvi_genera_relabel, expand = expansion(add = 1), name = NULL)
      + scale_y_continuous(name = "Relative abundance (%)", labels = scaleFUN2,
                           expand = expansion(mult = c(0.01,0.1)))
      + geom_hline(aes(yintercept = 0), color = "black")
      + guides(fill = guide_legend(order = 1, ncol = 1, title = "Variables",  direction = "vertical", 
                                   override.aes = list(color = "black", shape = 21, size = 3)),
               shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
                                    override.aes = list(color = "black", size = 3)),
               color = "none")
      + theme(strip.text.y = element_text(size = rel(0.8), angle = 0),
              strip.background = element_blank(),
              legend.position = "none",
              panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
              panel.grid.major.y = element_line(linetype = "solid", color = "grey80", linewidth = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = rel(0.8)),
              axis.text.y = element_text(angle = 0, size = rel(0.8)))
      # + facet_wrap(vars(grouping), ncol = 1,
      # + facet_wrap(.~enriched, ncol = 2,
      + facet_grid(rows = vars(grouping), cols = NULL, space = "free", axes = "all",
                   drop = TRUE, shrink = TRUE,
                   labeller = global_labeller,
                   scales = "free")
      + coord_flip()
      + ggtitle(paste0(namevar))
    )
    temp_g2 <- (
      ggplot(data = temp_df1 %>%
               dplyr::filter(abundrank == 1) %>%
               droplevels,
             aes(x = rev(asv_id), y = mean, group = interaction(asv_id, group_label)))
      + theme_bw() 
      + geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd), color = group_label), width = 0.3,
                      position = position_dodge(width = 0.4, preserve = "total"),
                      show.legend = FALSE)
      + geom_point(aes(fill = group_label, shape = site), na.rm = TRUE,
                   position = position_dodge(width = 0.4, preserve = "total"),
                   alpha = 1, color = "black", size = 2, show.legend = TRUE)
      + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
      + scale_discrete_manual(aesthetics = c("fill", "color"), 
                              values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
      + scale_x_discrete(labels = usvi_genera_relabel, expand = expansion(add = 1), name = NULL)
      + scale_y_continuous(name = "Relative abundance (%)", labels = scaleFUN2,
                           expand = expansion(mult = c(0.01,0.1)))
      + geom_hline(aes(yintercept = 0), color = "black")
      + guides(fill = guide_legend(order = 1, ncol = 1, title = "Variables",  direction = "vertical", 
                                   override.aes = list(color = "black", shape = 21, size = 3)),
               shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
                                    override.aes = list(color = "black", size = 3)),
               color = "none")
      + theme(strip.text.y = element_text(size = rel(0.8), angle = 0),
              strip.background = element_blank(),
              legend.position = "none",
              panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
              panel.grid.major.y = element_line(linetype = "solid", color = "grey80", linewidth = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = rel(0.8)),
              axis.text.y = element_text(angle = 0, size = rel(0.8)))
      # + facet_wrap(vars(grouping), ncol = 1,
      # + facet_wrap(.~enriched, ncol = 2,
      + facet_grid(rows = vars(grouping), cols = NULL, space = "free", axes = "all",
                   drop = TRUE, shrink = TRUE,
                   labeller = global_labeller,
                   scales = "free")
      + coord_flip()
      # + ggtitle(paste0(namevar))
    )
    
    temp_g <- temp_g1 | temp_g2
    assign(paste0("g4_sda_", i, "_", namevar2), temp_g, envir = .GlobalEnv, inherits = TRUE)
    if(!any(grepl(paste0(namevar2, Sys.Date(), collapse = "&"), list.files(projectpath, pattern = "usvi_sda_asvs_g4_sda_.*.png")))){
      ggsave(paste0(projectpath, "/", "usvi_sda_asvs_g4_sda_", i, "_", namevar2 , "-", Sys.Date(), ".png"),
             temp_g,
             width = 16, height = 10, units = "in")
    }
    rm(temp_g1)
    rm(temp_g)
    rm(temp_g2)
  }
  rm(namevar)
  rm(temp_df1)
}

#23: Tektite, 22: Lameshur, 21: Yawzi
group_label_site_shape <- tibble::enframe(group_labels_lookup, name = "group_label", value = "label") %>%
  dplyr::mutate(site = stringr::str_extract_all(group_label, paste0(names(site_lookup), collapse = "|")) %>%
                  as.character(.)) %>%
  dplyr::mutate(site = dplyr::case_when(grepl("LB", site) ~ 22,
                                        grepl("Yawzi", site) ~ 21,
                                        grepl("Tektite", site) ~ 23)) %>%
  dplyr::select(group_label, site) %>%
  tibble::deframe(.)
g4_legend <- print(ggplot(data = data.frame(
  (bind_rows(usvi_sda_asvs_both_mean_list, .id = NULL) %>%
     dplyr::ungroup(.) %>%
     dplyr::distinct(group_label, site) %>%
     droplevels)),
  aes(x = 1, y = 1, fill = group_label, shape = group_label)) 
  + geom_blank() 
  + geom_point() 
  + theme_void()
  + guides(fill = guide_legend(order = 1, ncol = 2, title = "Sample group",  direction = "vertical", 
                               override.aes = list(color = "black", size = 4, stroke = 1)),
           shape = guide_legend(order = 1, ncol = 1, title = "Sample group", direction = "vertical"),
           color = "none")
  + theme(legend.background = element_rect(fill = NA, colour = "grey30"))
  + scale_discrete_manual(aesthetics = c("fill"), 
                          values = group_labels_colors, labels = c(group_labels_lookup), breaks = names(group_labels_lookup))
  + scale_discrete_manual(aesthetics = c("shape"), values = group_label_site_shape, 
                          labels = (group_labels_lookup), 
                          breaks = names(group_labels_lookup),
                          drop = TRUE)) %>%
  g_legend()

gpatch_layout <- "
  AABBCC
  AABBCC
  AABBCC
"
gpatch1 <- apropos("g3_sda_.*_temporal_dawn.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`) 
gpatch1 <- (gpatch1 + patchwork::plot_layout(tag_level = "keep", design = gpatch_layout))  + patchwork::plot_annotation(title = "Spatially SDA ASVs at dawn", tag_levels = list(c("A", "B", "C")))
# gpatch1

gpatch2 <- apropos("g3_sda_.*_temporal_peak.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`) 
gpatch2 <- (gpatch2 + patchwork::plot_layout(tag_level = "keep", design = gpatch_layout))  + patchwork::plot_annotation(title = "Spatially SDA ASVs in afternoon", tag_levels = list(c("A", "B", "C")))
# gpatch2

gpatch3_layout <- "
  AAA
  AAA
  AAA
  #B#
"

gpatch3 <- apropos("g3_sda_.*_spatial.*", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `+`) 
gpatch3 <- (gpatch3 / g4_legend + patchwork::plot_layout(tag_level = "keep", design = gpatch3_layout))  + patchwork::plot_annotation(title = "Temporally SDA ASVs in each site", tag_levels = list(c("A", "B", "")))
# gpatch3


if(!any(grepl(paste0("temporal|spatial", Sys.Date(), collapse = "&"), list.files(projectpath, pattern = "usvi_sda_asvs_g3_sda_.*.png")))){
ggsave(paste0(projectpath, "/", "usvi_sda_asvs_g3_sda_", "spatial", "-", Sys.Date(), ".png"),
       gpatch3,
       width = 16, height = 8, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_g3_sda_", "temporal_afternoon", "-", Sys.Date(), ".png"),
         gpatch2,
         width = 16, height = 16, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_g3_sda_", "temporal_dawn", "-", Sys.Date(), ".png"),
         gpatch1,
         width = 16, height = 16, units = "in")
}



# Look at the 36 most abundant ASVs that were SDA -------------------------
temp_rankabund_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_abund_idx) %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::arrange(desc(baseMean), .by_group = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, baseMean) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  # tibble::deframe(.)
  tibble::deframe(.) %>% quantile(., probs = seq(0.75, 0.25, -0.25), names = FALSE) %>% setNames(., c("top", "middle", "bottom"))

shared_sda_asvs_abund_idx_cluster <- usvi_sda_asvs_compare.df %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::arrange(desc(baseMean), .by_group = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, baseMean) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::mutate(quantile_group = dplyr::case_when((baseMean >= temp_rankabund_idx[2]) ~ "top",
                                                  (baseMean < temp_rankabund_idx[2]) ~ "bottom",
                                                  .default = NA)) %>%
  dplyr::select(asv_id, quantile_group) %>%
  dplyr::mutate(quantile_group = factor(quantile_group, levels = c("top", "middle", "bottom", "rare"))) %>%
  dplyr::arrange(quantile_group) %>%
  dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  droplevels %>%
  dplyr::right_join(., usvi_sda_asvs_compare.df %>%
                      dplyr::filter(!grepl("dispersion", model)) %>%
                      dplyr::filter(asv_id %in% shared_sda_asvs_abund_idx) %>%
                      dplyr::distinct(contrast, asv_id, .keep_all = FALSE) %>%
                      dplyr::group_by(asv_id) %>%
                      dplyr::summarise(num_results = length(contrast)) %>%
                      dplyr::arrange(num_results) %>%
                      droplevels, by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
  split(., f = .$quantile_group) %>%
  map(., ~.x %>%
        dplyr::select(asv_id, num_results) %>%
        droplevels)

shared_sda_asvs_abund_cluster_list <- shared_sda_asvs_abund_idx_cluster %>%
  map(., ~.x %>%
      dplyr::inner_join(., usvi_sda_asvs_compare.df, by = join_by(asv_id), multiple = "all", relationship = "many-to-many")) %>%
  map(., ~.x %>%
        dplyr::mutate(across(c(asv_id, contrast, hold, variable, model, test_type), ~factor(.x)))) %>%
  map(., ~.x %>% 
        dplyr::distinct(asv_id, contrast, hold, variable, pair1, pair2) %>%
        droplevels %>%
        split(., f = .$hold) %>%
        map(., ~.x %>%
              dplyr::mutate(sampling_time = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ hold,
                                                             .default = NA),
                            site= dplyr::case_when((hold %in% names(site_lookup)) ~ hold,
                                                   .default = NA)) %>%
              droplevels) %>%
        map(., ~.x %>%
              tidyr::pivot_longer(., cols = c(pair1, pair2),
                                  names_to = NULL,
                                  values_to = "value") %>%
              dplyr::mutate(sampling_time = dplyr::case_when(is.na(sampling_time) ~ value,
                                                             .default = sampling_time),
                            site = dplyr::case_when(is.na(site) ~ value,
                                                    .default = site)) %>%
              dplyr::select(-value) %>%
              dplyr::distinct(.)) %>%
        bind_rows(., .id = NULL) %>%
        dplyr::left_join(., usvi_sda_asvs_relabund.df, 
                         by = join_by(asv_id, site, sampling_time), relationship = "many-to-many", multiple = "all") %>%
        tidyr::drop_na(.) %>%
        dplyr::mutate(pair1 = gsub("([[:print:]]+)( - )(.*)$", "\\1", variable),
                      pair2 = gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)) %>%
        tidyr::pivot_longer(., cols = c(pair1, pair2),
                            names_to = NULL,
                            values_to = "pair") %>%
        dplyr::rowwise(.) %>%
        dplyr::filter((pair %in% sampling_time) | (pair %in% site)) %>%
        # dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
        dplyr::mutate(hold2 = dplyr::case_when((hold %in% names(site_lookup)) ~ "spatial", 
                                               .default = "temporal")) %>%
        dplyr::distinct(.) %>%
        droplevels) %>%
  map(., ~.x %>%
        dplyr::mutate(across(c(asv_id, contrast, hold, variable, sampling_time, site, pair, hold2, sample_id), ~factor(.x))) %>%
        dplyr::mutate(pair = factor(pair, levels = c(names(site_lookup), names(sampling_time_lookup))),
                      sampling_time = factor(sampling_time, levels = names(sampling_time_lookup)),
                      site = factor(site, levels = names(site_lookup))) %>%
        droplevels)

temp_list <- shared_sda_asvs_abund_cluster_list %>%
  map(., ~.x %>%
        dplyr::left_join(., (usvi_sda_asvs_compare_summary.df %>%
                               dplyr::select(!c("test_type", c(Domain:Species))) %>%
                               tidyr::pivot_longer(., cols = !c(asv_id),
                                                   names_to = "contrast",
                                                   values_to = "significance") %>%
                               dplyr::mutate(significance = dplyr::case_when(grepl("_", significance) ~ "B",
                                                                             grepl("^deseq", significance) ~ "D",
                                                                             grepl("^corncob", significance) ~ "C",
                                                                             .default = NA)) %>%
                               # dplyr::mutate(contrast = gsub("peak_photo", "afternoon", contrast)) %>%
                               tidyr::drop_na(.) %>%
                               droplevels),
                         by = join_by(asv_id, contrast)) %>%
        dplyr::mutate(label_y = ceiling(relabund)*1.1) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(desc(label_y)) %>%
        dplyr::select(asv_id, contrast, significance, label_y) %>% 
        dplyr::distinct(asv_id, contrast, significance, .keep_all = TRUE) %>% 
        tidyr::drop_na(.) %>%
        droplevels) %>%
  map(., ~.x %>%
        dplyr::mutate(across(c(asv_id, contrast, significance), ~factor(.x))) %>%
        droplevels)
shared_sda_asvs_abund_cluster_list <- map(names(shared_sda_asvs_abund_cluster_list),
                                ~temp_list[[.x]] %>%
                                  dplyr::right_join(., shared_sda_asvs_abund_cluster_list[[.x]] %>%
                                                      # dplyr::filter(asv_id %in% namevar) %>%
                                                      droplevels,
                                                    by = join_by(asv_id, contrast), relationship = "one-to-many", multiple = "first") %>%
                                  droplevels) %>%
  setNames(., names(shared_sda_asvs_abund_cluster_list))

temp_df <- shared_sda_asvs_abund_cluster_list %>%
  map(., ~.x %>%
        dplyr::mutate(across(!"relabund", ~factor(.x))) %>%
        dplyr::arrange(asv_id, contrast, significance) %>%
        dplyr::distinct(asv_id, contrast, significance)) %>%
  bind_rows(., .id = NULL) %>%
  dplyr::filter(grepl("B", significance)) %>%
  droplevels

#filter for only those 181 ASVs flagged as SDA in both methods:
shared_sda_asvs_abund_filtered_list <- shared_sda_asvs_abund_cluster_list %>%
  map(., ~.x %>%
        dplyr::select(asv_id, contrast, site, sampling_time, sample_id, relabund) %>%
        dplyr::mutate(across(!"relabund", ~factor(.x))) %>%
        ungroup %>%      
        dplyr::left_join(., temp_df, by = join_by(asv_id, contrast), relationship = "many-to-many", multiple = "all") %>%
        dplyr::filter(grepl("B", significance)) %>%
        droplevels) %>%
  # temp_list2 <- shared_sda_asvs_idx_filtered_list %>%
  bind_rows(., .id = NULL) %>%
  dplyr::distinct(asv_id, contrast, significance, .keep_all = FALSE) %>%
  droplevels

shared_sda_asvs_abund_filtered_list <- shared_sda_asvs_abund_filtered_list %>%
  dplyr::left_join(., bind_rows(shared_sda_asvs_abund_cluster_list, .id = "quantile_group") %>%
                     dplyr::select(-significance)) %>%
  dplyr::mutate(group_label = dplyr::case_when(grepl("dawn", sampling_time) ~ "_d",
                                               grepl("photo", sampling_time) ~ "_p",
                                               .default = NA) %>%
                  paste0(site, .)) %>%
  dplyr::mutate(group_label = factor(group_label, levels = names(group_labels_lookup))) %>%
  dplyr::mutate(quantile_group = factor(quantile_group, levels = c("top", "middle", "bottom", "rare"))) %>%
  split(., f = .$quantile_group) %>%
  map(., ~.x %>%
        droplevels)

rm(list = apropos("g3_sda_.*", mode = "list"))

# for(i in seq_along(shared_sda_asvs_abund_cluster_list)){
#   
#   namevar <- names(shared_sda_asvs_abund_cluster_list[i])
#   temp_df1 <- shared_sda_asvs_abund_cluster_list[[i]] %>%
#     droplevels
# 
#   temp_g1 <- (
#     ggplot(data = temp_df1 %>%
#              droplevels, aes(y = relabund, x = pair, group = pair))
#     + theme_bw() 
#     + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
#                   alpha = 1, show.legend = FALSE, color = "grey",
#                   position = position_dodge(0.8))
#     + geom_point(aes(fill = sampling_time, shape = site), 
#                  alpha = 0.8, color = "black", size = 1, show.legend = TRUE, 
#                  position = position_jitter(width = 0.2, seed = 48105))
#     + geom_text(aes(label = significance, x = 1.5, vjust = "outward", hjust = "center", group = contrast,
#                     y = label_y),
#                 colour = "black", fontface = "bold")
#     + scale_shape_manual(name = "Sampling site", values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#     + scale_fill_manual(name = "Sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#     + scale_x_discrete(labels = c(site_lookup, sampling_time_lookup), name = NULL)
#     + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.5)))
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical", 
#                                  override.aes = list(color = "black", shape = 21, size = 3)),
#              shape = guide_legend(order = 1, ncol = 1, title = "Sampling site",  direction = "vertical", 
#                                   override.aes = list(color = "black", size = 3)),
#              color = "none")
#     + theme(strip.text.y = element_text(size = rel(0.7), angle = 0), 
#             strip.text.x = element_text(size = rel(0.7), angle = 0), 
#             axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 3),
#             axis.title.y = element_blank())
#     + facet_grid(asv_id ~ contrast,
#                  drop = TRUE,
#                  labeller = global_labeller,
#                  scales = "free")
#     + ggtitle(paste0("Observed SDA ASVs in ", namevar, " relative abundance ranks"))
#   )
#   
#   assign(paste0("g3_sda_", i, "_", namevar), temp_g1, envir = .GlobalEnv, inherits = TRUE)
#   rm(temp_g1)
#   rm(namevar)
# }
# 
# 
# gpatch <- apropos("g3_sda_.*", mode = "list") %>%
#   lapply(., get) %>%
#   purrr::reduce(., `+`)
# g3_sda_asvs <- (gpatch & theme(legend.position = "none")) + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")
# g3_sda_asvs
# 
# if(!any(grepl("asvs_most_abund", list.files(projectpath, pattern = "usvi_sda_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_sda_asvs_most_abund-", Sys.Date(), ".png"),
#          g3_sda_asvs,
#          width = 24, height = 15, units = "in")
# }





# temp_df1 <- usvi_sda_asvs_both_list[[1]] %>%
#   droplevels %>%
#   dplyr::group_by(asv_id, contrast, hold, variable, sampling_time, site, pair, hold2, group_label) %>%
#   dplyr::summarise(mean = mean(relabund, na.rm = TRUE),
#                    sd = sd(relabund, na.rm = TRUE), .groups = "keep") %>%
#   dplyr::left_join(., shared_sda_asvs_idx_cluster %>%
#                      bind_rows(., .id = "rank") %>%
#                      dplyr::select(asv_id, rank) %>%
#                      dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#                      dplyr::mutate(rank = factor(rank, levels = c("top", "middle", "bottom", "rare"))), by = join_by(asv_id))
  



# How many SDA ASVs belonged to genera that were also SDA? ----------------
#not necessary, since we don't care about the genera
{
# usvi_sda_asvs_in_sda_genera.df <- dplyr::full_join((usvi_sda_genera_compare_summary.df %>%
#                                                       dplyr::ungroup(.) %>%
#                                                       dplyr::distinct(asv_id, test_type) %>%
#                                                       # dplyr::mutate(test_type = dplyr::case_when(grepl("all", test_type) ~ 2,
#                                                       #                                                grepl("_", test_type) ~ 1.5,
#                                                       #                                                .default = 1)) %>%
#                                                       dplyr::mutate(test_type = 1) %>%
#                                                       dplyr::rename(sda_agg_genus = "test_type")),
#                                                    (usvi_sda_asvs_compare_summary.df %>%
#                                                       dplyr::ungroup(.) %>%
#                                                       dplyr::distinct(asv_id, test_type) %>%
#                                                       dplyr::mutate(test_type = 1) %>%
#                                                       # dplyr::mutate(test_type = dplyr::case_when(grepl("all", test_type) ~ 2,
#                                                       #                                                .default = 1)) %>%
#                                                       dplyr::rename(sda_asv = "test_type"))) %>%
#   dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"), by = join_by(asv_id)) %>%
#   droplevels %>%
#   dplyr::mutate(taxonomy = gsub("^([[:alnum:]_]{9}\\:\\s)(.*)", "\\2", taxonomy))
# 
# usvi_sda_asvs_in_sda_genera_summary.df <- usvi_sda_asvs_in_sda_genera.df %>%
#   # tidyr::pivot_longer(., cols = !c("asv_id", "taxonomy"),
#   #                     names_to = "taxon_resolution",
#   #                     values_to = "test_type") %>%
#   # dplyr::distinct(asv_id, taxonomy, taxon_resolution, .keep_all = TRUE) %>%
#   # dplyr::group_by(taxonomy, taxon_resolution) %>%
#   dplyr::group_by(taxonomy) %>%
#   dplyr::summarise(num_sig_asvs = sum(sda_asv, na.rm = TRUE),
#                    num_sig_genera = sum(sda_agg_genus, na.rm = TRUE)) %>%
#   droplevels
# 
# usvi_sda_resolution_summary.df <- usvi_prok_asvs.taxonomy %>%
#   dplyr::group_by(taxonomy) %>%
#   dplyr::summarise(num_agglom_asvs = length(asv_id)) %>%
#   dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
#   dplyr::full_join(., usvi_sda_asvs_in_sda_genera_summary.df,
#                    by = join_by(taxonomy)) %>%
#   dplyr::full_join(., usvi_sda_asvs_in_sda_genera_summary.df %>% 
#                      dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
#                      dplyr::select(taxonomy, num_sig_asvs) %>%
#                      dplyr::rename(num_sig_asvs_not_genera = num_sig_asvs),
#                    by = join_by(taxonomy)) %>%
#   # dplyr::full_join(., usvi_sda_asvs_in_sda_genera_summary.df %>%
#   #                    dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
#   #                    dplyr::select(taxonomy, num_sig_genera) %>%
#   #                    dplyr::rename(num_sig_genera_not_asvs = num_sig_genera),
#   #                  by = join_by(taxonomy)) %>%
#   tidyr::drop_na(num_sig_asvs, num_sig_genera) %>%
#   # dplyr::mutate(num_not_sig_asvs = dplyr::case_when(is.na(num_sig_asvs_not_genera) ~ num_agglom_asvs - num_sig_asvs,
#   dplyr::mutate(num_not_sig_asvs = dplyr::case_when(!is.na(num_sig_asvs) ~ num_agglom_asvs - num_sig_asvs,
#                                                     .default = NA)) %>%
#   dplyr::select(-c(num_sig_asvs_not_genera)) %>%
#   droplevels
# 
# 
# #how many genera were found to be significant at either the ASV or genus-level?
# #170 genera:
# usvi_sda_resolution_summary.df %>% 
#   dplyr::distinct(taxonomy) %>%
#   dplyr::summarise(num_sig_genera_overall = length(taxonomy))
# 
# #how many total ASVs are in these agglomerated genera?
# #8061 ASVs:
# usvi_sda_resolution_summary.df %>%
#   # dplyr::filter(num_sig_genera > 0) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::summarise(num_agglom_asvs = sum(num_agglom_asvs))
# 
# 
# #how many taxa were found to be significant at the genus-level, but not individual ASVs in there? 
# #how many ASVs were agglomerated to genus level, and genus-level was SDA, but not SDA at an ASV level?
# #66 genera encompassing 4636 ASVs:
# usvi_sda_resolution_summary.df %>%
#   dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
#   dplyr::summarise(num_asvs_in_sig_genera = sum(num_agglom_asvs),
#                    num_genera = length(taxonomy))
# # dplyr::summarise(num_sig_genera_not_asvs = sum(num_sig_genera))
# 
# #how many genera were SDA and had SDA ASVs:
# #97 genera
# usvi_sda_resolution_summary.df %>%
#   dplyr::filter(num_sig_asvs > 0 & num_sig_genera > 0) %>%
#   dplyr::summarise(num_genera = length(taxonomy))
# 
# #how many ASVs are significant, and belong to SDA genera?
# #250 ASVs in 97 genera:
# usvi_sda_resolution_summary.df %>%
#   dplyr::filter(num_sig_genera > 0 & num_sig_asvs > 0) %>%
#   dplyr::summarise(num_asvs_in_sig_genera = sum(num_sig_asvs),
#                    num_not_sig_asvs = sum(num_agglom_asvs) - sum(num_sig_asvs),
#                    num_genera = length(taxonomy))
# 
# 
# 
# #how many asvs were found to be significant but not at the agglomerated genera? 
# #8 ASVs reprsenting 7 genera:
# usvi_sda_resolution_summary.df %>%
#   dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0)
# 
# #how many ASVS in those 7 genera were not SDA?
# #138 ASVs:
# usvi_sda_resolution_summary.df %>%
#   dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::summarise(num_agglom_asvs = sum(num_agglom_asvs) - sum(num_sig_asvs),
#                    num_genera = length(taxonomy))
# 
# #how many ASVs were in sig genera but not sig at ASV levels?
# #7803 ASVs
# usvi_sda_resolution_summary.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::summarise(num_not_sig_asvs = sum(num_agglom_asvs) - sum(num_sig_asvs),
#                    num_genera = length(taxonomy))
# 
# 
# 
# # Plot the relative abundances of the 8061 ASVs in 170 genera -------------
# 
# 
# usvi_sda_genera_relabund.df <- usvi_sda_asvs_in_sda_genera.df %>%
#   dplyr::filter(!is.na(sda_agg_genus)) %>%
#   droplevels %>%
#   dplyr::distinct(asv_id, taxonomy) %>%
#   dplyr::mutate(significant = 1) %>%
#   dplyr::bind_rows(., usvi_sda_asvs_in_sda_genera.df %>%
#                      dplyr::filter(is.na(sda_agg_genus) & !is.na(sda_asv)) %>%
#                      dplyr::distinct(taxonomy, .keep_all = FALSE) %>%
#                      droplevels) %>%
#   dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
#   dplyr::left_join(., usvi_sw_genus.taxa.df %>%
#                      dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
#                      dplyr::select(asv_id, taxonomy) %>%
#                      droplevels, by = join_by(taxonomy)) %>%
#   dplyr::mutate(genera_id = across(starts_with("asv_id")) %>% purrr::reduce(coalesce)) %>%
#   dplyr::select(genera_id, taxonomy, significant) %>%
#   dplyr::right_join(usvi_sw_genus.df %>%
#                      dplyr::filter(sample_id %in% usvi_selected_metadata[["sample_id"]]) %>%
#                      dplyr::group_by(sample_id) %>%
#                      dplyr::mutate(relabund = relabund(counts)) %>%
#                      droplevels, ., by = join_by("asv_id" == "genera_id")) %>%
#   dplyr::rename(genera_id = "asv_id") %>%
#   dplyr::mutate(significant = dplyr::case_when(is.na(significant) ~ 0,
#                                                .default = significant)) %>%
#   droplevels
# 
# usvi_sda_asvs_and_relatives_relabund.df <- usvi_sda_asvs_in_sda_genera.df %>%
#   dplyr::filter(!is.na(sda_asv)) %>%
#   dplyr::select(asv_id, taxonomy) %>%
#   dplyr::mutate(significant = 1) %>%
#   dplyr::full_join(., usvi_prok_asvs.taxonomy %>%
#                      dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
#                      dplyr::select(asv_id, taxonomy) %>%
#                      dplyr::semi_join(usvi_sda_asvs_in_sda_genera.df,
#                                       by = join_by(taxonomy))) %>%
#   split(., f = .$taxonomy) %>%
#   purrr::map(., ~.x %>%
#                dplyr::arrange(asv_id) %>%
#                tibble::rowid_to_column(var = "rank_order") %>%
#                droplevels) %>%
#   bind_rows(.) %>%
#   dplyr::right_join(usvi_prok_asvs.df %>%
#                      dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
#                       dplyr::group_by(sample_ID) %>%
#                       dplyr::mutate(relabund = relabund(counts)) %>%
#                       dplyr::rename(sample_id = "sample_ID") %>%
#                      droplevels, ., by = join_by(asv_id)) %>%
#   dplyr::mutate(significant = dplyr::case_when(is.na(significant) ~ 0,
#                                                .default = significant)) %>%
#   droplevels
# 
# #first, plot the relative abundances of the SDA genera and genera of ASVs found to be SDA at the individual level
# print(
#   ggplot(data = usvi_sda_genera_relabund.df)
#   + geom_boxplot(aes(x = genera_id, y = relabund, fill = taxonomy), shape = 21, show.legend = FALSE)
# )
# 
# 
# #what is the distribution of agglomerated ASVs in each of the 170 genera?
# temp_df <- usvi_sda_asvs_and_relatives_relabund.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::arrange(desc(rank_order)) %>%
#   dplyr::select(taxonomy, rank_order) %>%
#   dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy)))
# 
# 
# 
# #two version of histograms:
# {
#   print(ggplot(data = temp_df)
#         + geom_histogram(aes(x = log10(rank_order)),  
#                          # bins= 5, binwidth = 0.8,
#                          breaks = log10(c(1, 10, 50, 100, 500, 1500)),
#                          color = "black", fill = "grey")
#         + scale_x_continuous(name = "Abundance of ASVs in agglomerated genera",
#                              breaks = (log10(c(1, 10, 50, 100, 500)*(1.5))),
#                              labels = c("Between 1 and 10", "Between 50 and 11", "Between 100 and 51",
#                                         "Between 500 and 101", "More than 501")
#         )
#         + scale_y_continuous(name = "Number of genus-level taxa")
#         + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#                 panel.grid.minor = element_blank(),
#                 panel.grid.major.x = element_blank(),
#                 panel.grid.major.y = element_line(color = "grey"),
#                 strip.text.x = element_text(angle = 0),
#                 strip.text.y = element_text(angle = 0))
#   )
#   
#   print(ggplot(data = temp_df)
#         + geom_histogram(aes(x = rank_order), breaks = c(0, 10, 50, 100, 500, 1500), color = "black", fill = "grey")
#         + scale_x_continuous(transform = scales::as.transform(t_pseudolog10), breaks =c(3, 20, 75, 300, 750),
#                              labels = c("Between 1 and 10", "Between 50 and 11", "Between 100 and 51",
#                                         "Between 500 and 101", "More than 501"),
#                              name = "Abundance of ASVs in agglomerated genera")
#         + scale_y_continuous(name = "Number of genus-level taxa")
#         + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#                 panel.grid.minor = element_blank(),
#                 panel.grid.major.x = element_blank(),
#                 panel.grid.major.y = element_line(color = "grey"),
#                 strip.text.x = element_text(angle = 0),
#                 strip.text.y = element_text(angle = 0))
#   )  
#   }
# 
# 
# 
# #97 genera had individual ASVs that were also SDA. 
# #what is the distribution of the relative abundances of those 250 ASVs and the other ASVs agglomerated into those genera?
# sda_genera_with_sda_asvs_idx <- intersect((usvi_sda_genera_relabund.df %>%
#             dplyr::ungroup(.) %>%
#             dplyr::filter(significant > 0) %>%
#             dplyr::distinct(taxonomy) %>%
#             unlist %>% as.character), (usvi_sda_asvs_and_relatives_relabund.df %>%
#                                         dplyr::ungroup(.) %>%
#                                         dplyr::filter(significant > 0) %>%
#                                         dplyr::distinct(taxonomy) %>%
#                                         unlist %>% as.character))
# 
# 
# temp_df2 <- temp_df %>%
#   dplyr::filter(taxonomy %in% sda_genera_with_sda_asvs_idx) %>%
#   dplyr::left_join(., usvi_sda_asvs_and_relatives_relabund.df %>%
#                      dplyr::filter(taxonomy %in% sda_genera_with_sda_asvs_idx) %>%
#                      dplyr::ungroup(.) %>%
#                      dplyr::distinct(asv_id, taxonomy, significant) %>%
#                      dplyr::group_by(taxonomy) %>%
#                      dplyr::summarise(num_sig = sum(significant, na.rm = TRUE)),
#                    by = join_by(taxonomy)) %>%
#   dplyr::arrange(desc(num_sig), desc(rank_order)) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#   droplevels
# 
# 
# usvi_sda_genera_relabund.list <- usvi_sda_asvs_and_relatives_relabund.df %>%
#   dplyr::select(sample_id, asv_id, rank_order, relabund, taxonomy, significant) %>%
#   dplyr::rename(taxon = "asv_id") %>%
#   dplyr::mutate(res = "asv") %>%
#   bind_rows(., usvi_sda_genera_relabund.df %>%
#               dplyr::mutate(taxon = paste0("g_", genera_id)) %>%
#               dplyr::select(sample_id, taxon, relabund, taxonomy, significant) %>%
#               dplyr::mutate(res = "agg_genus") %>%
#               droplevels) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#   dplyr::mutate(rank_order = dplyr::case_when(is.na(rank_order) ~ 0,
#                                               .default = rank_order)) %>%
#   droplevels
# 
# temp_list <- temp_df2 %>%
#                          dplyr::slice_max(n = 10, order_by = num_sig) %>%
#                          dplyr::select(taxonomy) %>%
#   droplevels %>%
#   dplyr::left_join(., usvi_sda_genera_relabund.list,
#                    by = join_by(taxonomy), relationship = "many-to-many", multiple = "all") %>%
#   droplevels
# 
# #oneata time:
# {
#   
#   # 
#   # 
#   # namevar <- unique(usvi_sda_genera_relabund.list[["taxonomy"]])
#   # print(
#   #   ggplot(data = usvi_sda_genera_relabund.list)
#   #   + theme_bw()
#   #   + geom_boxplot(aes(x = rank_order, y = relabund, fill = factor(significant), group = rank_order, color =factor(significant)), 
#   #                  alpha = 0.8,
#   #                  outliers = TRUE, outlier.shape = NA)
#   #   + scale_y_continuous(name = "Relative abundance (%)",
#   #                        # transform = scales::as.transform(t_pseudolog10), 
#   #                        breaks = c(0, 0.1, 0.5, 1, 5, 10, 25),
#   #                        transform = "log1p", 
#   #                        expand = expansion(mult = c(0.01,0.1)))
#   #   + scale_x_continuous(name = "Rank order of taxon")
#   #   # + scale_fill_manual(values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white")
#   #   + scale_discrete_manual(aesthetics = c("color"), values = c("grey", "black"), breaks = c(0, 1),
#   #                           drop = TRUE)
#   #   + scale_discrete_manual(aesthetics = c("fill"), values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white",
#   #                           drop = TRUE)
#   #   + guides(color = "none")
#   #   + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#   #           panel.grid.minor = element_blank(),
#   #           panel.grid.major.x = element_blank(),
#   #           panel.grid.major.y = element_line(color = "grey"),
#   #           strip.text.x = element_text(angle = 0),
#   #           strip.text.y = element_text(angle = 0))
#   #   + ggtitle(paste0("Taxa of ", namevar))
#   # )
#   
#   
# }
# 
# y <- unique(temp_list[["taxonomy"]]) %>% unlist %>% as.character(.)
# 
# # for(i in seq_len(2)){
# for(i in seq_len(length(y))){
#   namevar <- y[i]
#   temp_df_to_plot <- temp_list %>%
#     dplyr::filter(taxonomy == namevar) %>%
#     droplevels
#   x_breaks <- temp_df_to_plot %>%
#     dplyr::distinct(rank_order) %>%
#     dplyr::arrange(rank_order) %>%
#     tibble::deframe(.)
#   x_labels <- c("Agglom", x_breaks[-1])
#   
#   temp_g <- print(
#     ggplot(data = temp_df_to_plot)
#     + theme_bw()
#     + geom_boxplot(aes(x = rank_order, y = relabund, fill = factor(significant), group = rank_order, color =factor(significant)), 
#                    alpha = 0.8,
#                    outliers = TRUE, outlier.shape = NA)
#     + scale_y_continuous(name = "Relative abundance (%)",
#                          transform = "log1p", 
#                          expand = expansion(mult = c(0.01,0.1)))
#     + scale_x_continuous(name = "Rank order of taxon", labels = x_labels, breaks = x_breaks)
#     + scale_discrete_manual(aesthetics = c("color"), values = c("grey", "black"), breaks = c(0, 1),
#                             drop = TRUE)
#     + scale_discrete_manual(aesthetics = c("fill"), values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white",
#                             drop = TRUE)
#     + guides(color = "none")
#     + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#             panel.grid.minor = element_blank(),
#             panel.grid.major.x = element_blank(),
#             panel.grid.major.y = element_line(color = "grey"),
#             strip.text.x = element_text(angle = 0),
#             strip.text.y = element_text(angle = 0))
#     + ggtitle(paste0(namevar))
#   )
#   assign(paste0("g_sda_asv_relabund_", i), temp_g, envir = .GlobalEnv)
# }
# 
# gpatch <- lapply(apropos("^g_sda_asv_relabund_.*$", mode = "list"),
#                  get) %>%
#   purrr::reduce(., `+`) + 
#   patchwork::plot_layout(guides = "collect")
#   # patchwork::plot_layout(guides = "collect") & theme(legend.position="none")
# 
# g_top_sda_asv_relabund <- gpatch + patchwork::plot_annotation(title = "Distribution of SDA ASVs and their SDA genera",
#                                                              tag_levels = "A")
# 
# 
# if(!any(grepl("top_sda_asv_relabund", list.files(projectpath, pattern = "usvi_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_top_sda_asv_relabund-", Sys.Date(), ".png"),
#          g_top_sda_asv_relabund, 
#          width = 16, height = 16, units = "in")
# }
# 
# #what about the 66 genera whose individual ASVs were not SDA?
# sda_genera_not_sda_asvs_idx <- setdiff((usvi_sda_genera_relabund.df %>%
#                                           dplyr::ungroup(.) %>%
#                                           dplyr::filter(significant > 0) %>%
#                                           dplyr::distinct(taxonomy) %>%
#                                           unlist %>% as.character), (usvi_sda_asvs_and_relatives_relabund.df %>%
#                                                                        dplyr::ungroup(.) %>%
#                                                                        dplyr::filter(significant > 0) %>%
#                                                                        dplyr::distinct(taxonomy) %>%
#                                                                        unlist %>% as.character))
# 
# #these are the 4636 ASVs belonging to the 66 SDA genera
# temp_df3 <- usvi_prok_asvs.taxonomy %>%
#   dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
#   dplyr::filter(taxonomy %in% sda_genera_not_sda_asvs_idx) %>%
#   dplyr::select(asv_id, taxonomy) %>%
#   split(., f = .$taxonomy) %>%
#   purrr::map(., ~.x %>%
#                dplyr::arrange(asv_id) %>%
#                tibble::rowid_to_column(var = "rank_order") %>%
#                droplevels) %>%
#   bind_rows(.) %>%
#   dplyr::left_join(., (usvi_prok_asvs.df %>%
#                          dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
#                          dplyr::group_by(sample_ID) %>%
#                          dplyr::mutate(relabund = relabund(counts)) %>%
#                          dplyr::rename(sample_id = "sample_ID") %>%
#                          droplevels),
#                    by = join_by(asv_id)) %>%
#   dplyr::arrange(desc(relabund)) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#   droplevels
# temp_df3 %>%
#   dplyr::mutate(relabund = dplyr::case_when(relabund > 0 ~ relabund,
#                                             .default = NA)) %>%
#   # dplyr::group_by(taxonomy) %>%
#   dplyr::summarise(relabund = sum(relabund, na.rm = TRUE), .by = taxonomy) %>% dplyr::arrange(desc(relabund)) %>%
#   # dplyr::ungroup(.) %>% dplyr::summarise(quantiles = quantile(relabund, probs = seq(0, 1, 0.25), na.rm = TRUE)) %>%
#   droplevels
# 
# #examinethe top 10 genera:
# temp_list2 <- usvi_sda_genera_relabund.list %>%
#   dplyr::filter(taxonomy %in% sda_genera_not_sda_asvs_idx) %>%
#   dplyr::filter(res == "agg_genus") %>%
#   dplyr::ungroup(.) %>%
#   dplyr::arrange(desc(relabund)) %>%
#   dplyr::distinct(taxon, taxonomy, .keep_all = TRUE) %>%
#   dplyr::slice_max(n = 10, order_by = relabund) %>%
#   dplyr::semi_join((usvi_sda_genera_relabund.list %>%
#                      dplyr::filter(taxonomy %in% sda_genera_not_sda_asvs_idx) %>%
#                      # dplyr::filter(res == "agg_genus") %>%
#                      droplevels), ., by = join_by(taxonomy, taxon)) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#   droplevels
# 
# #here are the ASVs belonging to those top 10 genera
# temp_list3 <- temp_list2 %>%
#   bind_rows(temp_df3 %>%
#               dplyr::filter(taxonomy %in% unique(temp_list2[["taxonomy"]])) %>%
#               droplevels %>%
#             dplyr::rename(taxon = "asv_id") %>%
#               dplyr::mutate(significant = 0)) %>%
#   dplyr::select(sample_id, taxon, rank_order, relabund, taxonomy, significant) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#   dplyr::mutate(rank_order = dplyr::case_when(is.na(rank_order) ~ 0,
#                                               .default = rank_order)) %>%
#   droplevels
# 
# y <- unique(temp_list3[["taxonomy"]]) %>% unlist %>% as.character(.)
# 
# # for(i in seq_len(2)){
# for(i in seq_len(length(y))){
#   namevar <- y[i]
#   temp_df_to_plot <- temp_list3 %>%
#     dplyr::filter(taxonomy == namevar) %>%
#     droplevels
#   x_breaks <- temp_df_to_plot %>%
#     ungroup(.) %>%
#     dplyr::distinct(rank_order) %>%
#     dplyr::arrange(rank_order) %>%
#     tibble::deframe(.)
#   x_labels <- c("Agglom", x_breaks[-1])
#   
#   temp_g <- print(
#     ggplot(data = temp_df_to_plot)
#     + theme_bw()
#     + geom_boxplot(aes(x = rank_order, y = relabund, fill = factor(significant), group = rank_order, color =factor(significant)), 
#                    alpha = 0.8,
#                    outliers = TRUE, outlier.shape = NA)
#     + scale_y_continuous(name = "Relative abundance (%)",
#                          transform = "log1p", 
#                          expand = expansion(mult = c(0.01,0.1)))
#     + scale_x_continuous(name = "Rank order of taxon", labels = x_labels, breaks = x_breaks)
#     + scale_discrete_manual(aesthetics = c("color"), values = c("grey", "black"), breaks = c(0, 1),
#                             drop = TRUE)
#     + scale_discrete_manual(aesthetics = c("fill"), values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white",
#                             drop = TRUE)
#     + guides(color = "none")
#     + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
#             panel.grid.minor = element_blank(),
#             panel.grid.major.x = element_blank(),
#             panel.grid.major.y = element_line(color = "grey"),
#             strip.text.x = element_text(angle = 0),
#             strip.text.y = element_text(angle = 0))
#     + ggtitle(paste0(namevar))
#   )
#   assign(paste0("g_sda_genus_relabund_", i), temp_g, envir = .GlobalEnv)
# }
# 
# gpatch2 <- lapply(apropos("^g_sda_genus_relabund_.*$", mode = "list"),
#                  get) %>%
#   purrr::reduce(., `+`) + 
#   patchwork::plot_layout(guides = "collect")
# # patchwork::plot_layout(guides = "collect") & theme(legend.position="none")
# 
# g_top_sda_genera_not_asvs_relabund <- gpatch2 + patchwork::plot_annotation(title = "Distribution of SDA genera and their non-SDA ASVs",
#                                                               tag_levels = "A")
# 
# 
# if(!any(grepl("top_sda_genera_not_asvs_relabund", list.files(projectpath, pattern = "usvi_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_top_sda_genera_not_asvs_relabund-", Sys.Date(), ".png"),
#          g_top_sda_genera_not_asvs_relabund, 
#          width = 16, height = 16, units = "in")
# }
}




