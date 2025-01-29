# 08_taxa_metabolites_corr.R

#evaluate realtionships to metabolites

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
# library(phyloseq)
library(ggvegan)
library(ggtext)
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)


# Custom functions --------------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_custom_functions.R"))){
  cat("Loading custom functions.")
  source(paste0(getwd(), "/", "00_custom_functions.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 


# Read in metadata --------------------------------------------------------

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


# if(!exists("ps_usvi", envir = .GlobalEnv)){
#   if(file.exists(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))){
#     # ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))
#     ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))
#   } else {
#     cli::cli_alert_warning("Please process the USVI data through Phyloseq.")
#     
#   }
# }
if(file.exists(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"))){
  metabolomics_sample_metadata <- readr::read_delim(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
} else {
  metabolomics_sample_metadata <- ps_usvi %>%
    phyloseq::sample_data(.) %>%
    tibble::as_tibble(., rownames = "sample_id") %>%
    tidyr::drop_na(metab_deriv_label) %>%
    tidyr::separate_longer_delim(., metab_deriv_label, delim = ", ") %>%
    dplyr::mutate(metab_deriv_label = stringr::str_remove_all(metab_deriv_label, "[[:punct:]]") %>%
                    stringr::str_replace_all(., "[[:alpha:]]{1,}", "CINAR_BC_")) %>%
    dplyr::mutate(metab_deriv_label = dplyr::case_when(grepl("Metab_219", sample_id) ~ "CINAR_BC_81B",
                                                       grepl("Metab_319", sample_id) ~ "CINAR_BC_81A",
                                                       .default = metab_deriv_label)) %>%
    droplevels
  
  readr::write_delim(metabolomics_sample_metadata, paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
  
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


#don't need to genus-level assessments:
{
  # usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
  #   dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus) %>%
  #   dplyr::distinct(Domain, Phylum, Class, Order, Family, Genus, .keep_all = TRUE) %>%
  #   droplevels
  # 
  # usvi_sw_genus.tbl <- usvi_prok_asvs.df %>%
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
  #   droplevels
  
}

# Prepare color palette ---------------------------------------------------

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

if(!exists("annotation_taxa_colors_list", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "annotation_taxa_colors_list.rds"))){
    annotation_taxa_colors_list <- readr::read_rds(paste0(projectpath, "/", "annotation_taxa_colors_list.rds"))
  } else {
    cli::cli_alert_warning("Please prepare a list of colors for taxonomy.")
  }
}
# Read in metabolites -----------------------------------------------------

if(file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
  temp_list <- readr::read_rds(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
  list2env(temp_list, envir = .GlobalEnv)
  rm(temp_list)
} else {
  cli::cli_alert_warning("Please tidy the metabolomics datasets.")
}



# Read in significant asvs/genera -----------------------------------------

#don't need to evaluate significant genera
{
  # if(any(grepl("rds$", list.files(projectpath, pattern = "usvi_sda_genera_compare.*")))){
  #   usvi_sda_genera_compare.df <- readr::read_rds(paste0(projectpath, "/", "usvi_sda_genera_compare.rds"))
  # } else {
  #   if(any(grepl("tsv$", list.files(projectpath, pattern = "usvi_sda_genera_compare.*")))){
  #     usvi_sda_genera_compare.df <- readr::read_delim(paste0(projectpath, "/", "usvi_sda_genera_compare.tsv"),
  #                                                     delim = "\t", col_names = TRUE, show_col_types = FALSE)
  #   } else {
  #     cli::cli_alert_info("Reading in rademu and deseq results to parse...")
  #     
  #     #read in the genera that were significant by site through rademu:
  #     if(any(grepl("genera_reef_highlow", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
  #       temp_file <- list.files(projectpath, pattern = "usvi_rademu_genera_reef_highlow.*.RData")[1]
  #       load(paste0(projectpath, "/", temp_file))
  #       rm(temp_file)
  #     }
  #     
  #     if(any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
  #       temp_file <- list.files(projectpath, pattern = "usvi_rademu_genera_site_time.*.RData")[1]
  #       load(paste0(projectpath, "/", temp_file))
  #       rm(temp_file)
  #     }
  #     
  #     rademu_res.df <- bind_rows((score_1_res_list %>%
  #                                   bind_rows(.) %>% #these were genera that were found to be higher in samples collected at dawn
  #                                   dplyr::mutate(group = "high_dawn") %>%
  #                                   dplyr::filter(pval < 0.05) %>%
  #                                   tidyr::drop_na(score_stat)),
  #                                (score_2_res_list %>%
  #                                   bind_rows(.) %>% # these were genera that were found to be higher in samples collected at peak photo
  #                                   dplyr::mutate(group = "low_dawn") %>%
  #                                   dplyr::filter(pval < 0.05) %>%
  #                                   tidyr::drop_na(score_stat)),
  #                                (score_4_res_list %>%
  #                                   bind_rows(.) %>% #these were genera that were found to be higher in samples collected in seagrass
  #                                   dplyr::mutate(group = "high_seagrass") %>%
  #                                   dplyr::filter(pval < 0.05) %>%
  #                                   tidyr::drop_na(score_stat)),
  #                                (score_5_res_list %>%
  #                                   bind_rows(.) %>% #these were genera that were found to be higher in samples collected in reef sites
  #                                   dplyr::mutate(group = "low_seagrass") %>%
  #                                   dplyr::filter(pval < 0.05) %>%
  #                                   tidyr::drop_na(score_stat)),
  #                                (score_lowreef_res_list %>%
  #                                   bind_rows(.) %>%
  #                                   dplyr::mutate(group = "low_reef") %>%
  #                                   dplyr::filter(pval < 0.05) %>%
  #                                   tidyr::drop_na(score_stat)),
  #                                (score_highreef_res_list %>%
  #                                   bind_rows(.) %>%
  #                                   dplyr::mutate(group = "high_reef") %>%
  #                                   dplyr::filter(pval < 0.05) %>%
  #                                   tidyr::drop_na(score_stat))) %>%
  #       dplyr::rename(asv_id = "category") %>%
  #       droplevels %>%
  #       dplyr::rename(contrast = "covariate") %>%
  #       dplyr::mutate(contrast = stringr::str_remove_all(contrast, "site") %>%
  #                       stringr::str_replace_all(., ":sampling_time", "_")) %>%
  #       dplyr::mutate(contrast = dplyr::case_when(grepl("reef", group) ~ paste0(contrast, " - LB_seagrass_all"),
  #                                                 grepl("dawn", group) ~ paste0(contrast, " - all_peak_photo"),
  #                                                 grepl("seagrass", group) ~ paste0(contrast, " - LB_seagrass_all"),
  #                                                 .default = contrast)) %>%
  #       dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
  #       dplyr::mutate(commonid = stringr::str_split_i(group, "_", 2)) %>%
  #       dplyr::arrange(commonid, asv_id) %>%
  #       dplyr::group_by(commonid) %>%
  #       dplyr::distinct(asv_id, group, contrast, .keep_all = TRUE) %>%
  #       droplevels
  #     
  #     ## read in the results from DESeq:
  #     
  #     if(any(grepl("asvs_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  #       temp_file <- list.files(projectpath, pattern = "usvi_deseq_asvs_site_time.*.RData")[1]
  #       load(paste0(projectpath, "/", temp_file))
  #       rm(temp_file)
  #     }
  #     if(any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  #       temp_file <- list.files(projectpath, pattern = "usvi_deseq_genera_site_time.*.RData")[1]
  #       load(paste0(projectpath, "/", temp_file))
  #       rm(temp_file)
  #     }
  #     usvi_sda_genera_compare.df <- rademu_res.df %>%
  #       dplyr::mutate(test_type = "rademu") %>%
  #       dplyr::rename(estimate = "score_stat") %>%
  #       dplyr::mutate(model = "auto_fit") %>%
  #       dplyr::mutate(sampling_time = dplyr::case_when(grepl("dawn|peak", contrast) ~ stringr::str_extract(contrast, "dawn|peak_photo"),
  #                                                      .default = "all")) %>%
  #       dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
  #       dplyr::mutate(Combo = dplyr::case_when(grepl("Yawzi_|Tektite_", contrast) ~ stringr::str_extract(contrast, "Yawzi_.|Tektite_."),
  #                                              .default = site)) %>%
  #       dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
  #       dplyr::mutate(baseline = dplyr::case_when(grepl("all_peak_photo", baseline) ~ "all_peak_photo",
  #                                                 grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
  #                                                 .default = "LB_seagrass_all")) %>%
  #       bind_rows(., (usvi_deseq_genera_abund.df %>%
  #                       # bind_rows(., (usvi_deseq_genera_abund_filtered.df %>%
  #                       dplyr::mutate(test_type = "deseq"))) %>%
  #       dplyr::mutate(pair1 = stringr::str_split_i(contrast, " - ", 1),
  #                     pair2 = stringr::str_split_i(contrast, " - ", 2))
  #   }
  # }
  # 
  # #which genera are flagged as significantly differentially abundant in both tests?
  # shared_sda_genera_idx <- usvi_sda_genera_compare.df %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  #   dplyr::group_by(asv_id) %>%
  #   dplyr::summarise(num_results = length(test_type)) %>%
  #   dplyr::filter(num_results > 1) %>%
  #   # dplyr::filter(num_results > 2) %>%
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
  # #which are in DESeq2 but not in rademu?
  # 
  # setdiff(unique(usvi_sda_genera_compare.df %>%
  #                  dplyr::ungroup(.) %>%
  #                  tidyr::drop_na(baseMean) %>%
  #                  dplyr::distinct(asv_id) %>%
  #                  unlist), unique(usvi_sda_genera_compare.df %>%
  #                                    dplyr::ungroup(.) %>%
  #                                    tidyr::drop_na(group) %>%
  #                                    dplyr::distinct(asv_id) %>%
  #                                    unlist))
  # #what are in radEmu not in DeSeq2?
  # setdiff(unique(usvi_sda_genera_compare.df %>%
  #                  dplyr::ungroup(.) %>%
  #                  tidyr::drop_na(group) %>%
  #                  dplyr::distinct(asv_id) %>%
  #                  unlist), unique(usvi_sda_genera_compare.df %>%
  #                                    dplyr::ungroup(.) %>%
  #                                    tidyr::drop_na(baseMean) %>%
  #                                    dplyr::distinct(asv_id) %>%
  #                                    unlist))
  # 
  # #how many genera were "low_seagrass" and "high_reef"?
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "high_reef") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "low_seagrass") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # # [1] "ASV_00002" "ASV_00007" "ASV_00017" "ASV_00023" "ASV_00031" "ASV_00072" "ASV_00129" "ASV_00138"
  # # [9] "ASV_00150" "ASV_00171"
  # 
  # #how many genera were "low_reef" and "high_seagrass"?
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "high_seagrass") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "low_reef") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # # [1] "ASV_00003" "ASV_00006" "ASV_00009" "ASV_00038" "ASV_00045" "ASV_00060" "ASV_00097" "ASV_00105"
  # # [9] "ASV_00111" "ASV_00124" "ASV_00131" "ASV_00251" "ASV_00271"
  # 
  # 
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "high_dawn") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "high_reef") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # # [1] "ASV_00002" "ASV_00007" "ASV_00017" "ASV_00023" "ASV_00031" "ASV_00072" "ASV_00129" "ASV_00138"
  # # [9] "ASV_00150" "ASV_00171"
  # 
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "low_dawn") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "low_reef") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # # [1] "ASV_00003" "ASV_00006" "ASV_00009" "ASV_00038" "ASV_00045" "ASV_00060" "ASV_00097" "ASV_00105"
  # # [9] "ASV_00111" "ASV_00124" "ASV_00131" "ASV_00215" "ASV_00227" "ASV_00243" "ASV_00251" "ASV_00282"
  # # [17] "ASV_01134"
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "low_dawn") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "high_reef") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # #0
  # 
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "high_dawn") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "low_reef") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # #0
  # 
  # intersect(unique(usvi_sda_genera_compare.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    tidyr::drop_na(group) %>%
  #                    dplyr::filter(group == "high_dawn") %>%
  #                    dplyr::distinct(asv_id) %>%
  #                    unlist), unique(usvi_sda_genera_compare.df %>%
  #                                      dplyr::ungroup(.) %>%
  #                                      tidyr::drop_na(group) %>%
  #                                      dplyr::filter(group == "high_seagrass") %>%
  #                                      dplyr::distinct(asv_id) %>%
  #                                      unlist))
  # #0
  # 
  # 
  # #how many are just significantly more abundant in seagrass?
  # usvi_sda_genera_compare.df %>%
  #   dplyr::ungroup(.) %>%
  #   tidyr::drop_na(group) %>%
  #   dplyr::distinct(asv_id, group) %>%
  #   dplyr::group_by(group) %>%
  #   dplyr::summarise(num_obs = length(asv_id))
  # 
  # usvi_sda_genera_compare.df %>%
  #   dplyr::ungroup(.) %>%
  #   tidyr::drop_na(baseMean) %>%
  #   dplyr::distinct(asv_id, baseline, .keep_all = TRUE) %>%
  #   dplyr::mutate(highlow = dplyr::case_when(log2FoldChange_MMSE < 0 ~ -1,
  #                                            log2FoldChange_MMSE > 0 ~ 1,
  #                                            .default = NA)) %>%
  #   dplyr::group_by(baseline, Combo, highlow) %>%
  #   dplyr::summarise(num_obs = length(asv_id))
  # 
  # 
}

# Correlations between DNA profile and metabolomics -----------------------


drop <- c("CINAR_BC_73")

#old method:
{
  # usvi_metab_df <- usvi_metabolomics.df %>%
  #   dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  #   dplyr::left_join(., (metabolomics_sample_metadata %>%
  #                          dplyr::filter(grepl("seawater", sample_type)) %>%
  #                          dplyr::select(sample_id, metab_deriv_label) %>%
  #                          droplevels),
  #                    by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
  #   dplyr::relocate(sample_id) %>%
  #   tidyr::pivot_longer(., cols = !c(sample_id, metab_deriv_label),
  #                       names_to = "simpleName",
  #                       values_to = "conc") %>%
  #   dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::group_by(sample_id, simpleName) %>%
  #   # dplyr::summarise(num = length(conc)) %>%
  #   dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE),  
  #                    num = length(conc),
  #                    # .by = c(sample_id, simpleName),
  #                    .groups = "keep",
  #                    sd = sd(conc, na.rm = TRUE)) %>%
  #   dplyr::rename(conc = "mean_conc") %>%
  #   dplyr::select(sample_id, simpleName, conc) %>%
  #   droplevels
  # 
  # usvi_metab_mat <- usvi_metab_df %>%
  #   dplyr::group_by(simpleName) %>%
  #   dplyr::mutate(conc = scales::rescale(conc)) %>%
  #   tidyr::pivot_wider(., id_cols = "sample_id",
  #                      names_from = "simpleName",
  #                      values_fill = 0,
  #                      values_from = "conc") %>%
  #   tibble::column_to_rownames(var = "sample_id") %>%
  #   as.matrix(.) %>%
  #     vegan::vegdist(., binary = FALSE, upper = TRUE,
  #                    # distance = "horn", 
  #                    distance = "bray",
  #                    autotransform = TRUE) %>%
  #   as.matrix(.)
  # usvi_asv.tbl <- usvi_prok_asvs.df %>%
  #   dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
  #   dplyr::filter(sample_ID %in% rownames(usvi_metab_mat)) %>%
  #   # dplyr::filter(sample_ID %in% (metadata %>% dplyr::filter(grepl("seawater", sample_type)) %>% dplyr::select(sample_id) %>% unlist)) %>%
  #   droplevels %>%
  #   tidyr::pivot_wider(., id_cols = "asv_id",
  #                      names_from = "sample_ID",
  #                      values_from = "counts",
  #                      values_fill = 0) %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   # dplyr::slice(which(rowSums(.) > 0)) %>%
  #   # usvi_asv.tbl <- ps_usvi %>%
  #   #   phyloseq::subset_samples(., sample_type == "seawater") %>%
  #   #   phyloseq::otu_table(.) %>%
  #   # as.data.frame %>%
  #   apply(., 2, relabund) %>% 
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   tidyr::pivot_longer(., cols = -c("asv_id"),
  #                       names_to = "sample",
  #                       values_to = "abundance") %>%
  #   dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
  #                                   log2(abundance+1), #log transform abundance (with +1 pseudocount)
  #                                   0)) %>%
  #   tidyr::pivot_wider(., id_cols = "sample",
  #                      values_from = "abundance",
  #                      # values_from = "logabund",
  #                      names_from = "asv_id") %>%
  #   tibble::column_to_rownames(var = "sample")
}

#because some 16S samples have technical replicates in the metabolite samples, 
#we need to either:
#A. average the metabolite profiles that are technical replicates, so there is 1 metabolome sample per microbiome sample
#B. keep the technical replicates, but duplicate the corresponding microbiome sample

#option A:
#with option A, dropping day 1 microbiome samples as well as the corresponding 16S sample to CINAR_BC_73,
#we lose 2076 ASVs representing between 0-8.2% of the community
#these 2076 ASVs contribute <1% to the samples from Days 2-5

if(file.exists(paste0(projectpath, "/", "usvi_dist_mat_list_optA", ".rds"))){
  usvi_dist_mat_list_optA <- readr::read_rds(paste0(projectpath, "/", "usvi_dist_mat_list_optA", ".rds"))
} else {
  usvi_metab_mat <- usvi_metabolomics.df %>%
    dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::relocate(sample_id) %>%
    tidyr::pivot_longer(., cols = !c(sample_id, metab_deriv_label),
                        names_to = "simpleName",
                        values_to = "conc") %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(sample_id, simpleName) %>%
    # dplyr::summarise(num = length(conc)) %>%
    dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE),  
                     num = length(conc),
                     # .by = c(sample_id, simpleName),
                     .groups = "keep",
                     sd = sd(conc, na.rm = TRUE)) %>%
    dplyr::rename(conc = "mean_conc") %>%
    dplyr::select(sample_id, simpleName, conc) %>%
    dplyr::mutate(conc = log2(conc + 1)) %>%
    tidyr::pivot_wider(., id_cols = "sample_id",
                       names_from = "simpleName",
                       values_fill = 0,
                       values_from = "conc") %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::vegdist(., binary = FALSE, upper = TRUE,
                   # distance = "horn",
                   distance = "bray",
                   autotransform = TRUE) %>%
    as.matrix(.)
  
  sample_relabel <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::select(metab_deriv_label, sample_id, site, sampling_day, sampling_time) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    dplyr::left_join(., metadata %>%
                       dplyr::select(sample_id, replicate),
                     by = join_by(sample_id)) %>%
    droplevels %>%
    dplyr::mutate(replA = LETTERS[replicate]) %>%
    dplyr::select(sample_id, site, sampling_day, sampling_time, replA) %>%
    dplyr::arrange(site, sampling_time, sampling_day, replA) %>%
    droplevels %>%
    tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time, replA), sep = "_", remove = FALSE)  %>%
    dplyr::distinct(sample_id, relabeled_sample) %>%
    tibble::deframe(.)
  
  # sample_relabel <- metabolomics_sample_metadata %>%
  #   dplyr::select(sample_id, site, sampling_day, sampling_time) %>%
  #   dplyr::distinct(., .keep_all = TRUE) %>%
  #   dplyr::arrange(site, sampling_time, sampling_day) %>%
  #   droplevels %>%
  #   tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time), sep = "_", remove = TRUE)  %>%
  #   tibble::deframe(.)
  
  
  # temp_mat <- usvi_prok_asvs.df %>%
  #   dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
  #   dplyr::right_join(., (metabolomics_sample_metadata %>%
  #                          dplyr::filter(grepl("seawater", sample_type)) %>%
  #                          dplyr::distinct(sample_id) %>%
  #                          droplevels),
  #                    by = join_by("sample_ID" == "sample_id"), multiple = "all", relationship = "many-to-many") %>%
  #   droplevels %>%
  #   tidyr::pivot_wider(., id_cols = "asv_id",
  #                      names_from = "sample_ID",
  #                      values_from = "counts",
  #                      values_fill = 0) %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   tidyr::drop_na(.) %>%
  #   as.matrix(.)
  # dropped_asv_idx <- names(which(rowSums(temp_mat[,colnames(usvi_metab_mat)]) == 0))
  # data.frame(`sample_id` = names(which(colSums(temp_mat[dropped_asv_idx,]) < 1))) %>%
  #   dplyr::left_join(., (metabolomics_sample_metadata %>%
  #                          dplyr::filter(grepl("seawater", sample_type)) %>%
  #                          dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
  #                          dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  #                          droplevels)) %>%
  #   dplyr::arrange(sampling_day, site, sampling_time)
  
  
  usvi_asv_mat <- usvi_prok_asvs.df %>%
    dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
    dplyr::right_join(., (metabolomics_sample_metadata %>%
                            dplyr::filter(grepl("seawater", sample_type)) %>%
                            dplyr::distinct(sample_id) %>%
                            droplevels),
                      by = join_by("sample_ID" == "sample_id"), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(sample_ID %in% colnames(usvi_metab_mat)) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "asv_id",
                       names_from = "sample_ID",
                       values_from = "counts",
                       values_fill = 0) %>%
    tibble::column_to_rownames(var = "asv_id") %>%
    apply(., 2, relabund) %>%
    as.data.frame(.) %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    tidyr::drop_na(.) %>%
    dplyr::select(colnames(usvi_metab_mat)) %>%
    t() %>%
    vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
                   autotransform = TRUE) %>%
    as.matrix(.)
  
  meta.seawater <- metabolomics_sample_metadata %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    dplyr::filter(sample_id %in% rownames(usvi_asv_mat)) %>%
    dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
    dplyr::select(!c(contains("label"), contains("dna_"))) %>%
    tibble::column_to_rownames(., var = "sample_id") %>%
    droplevels
  
  usvi_dist_mat_list_optA <- list(usvi_metab_mat, usvi_asv_mat, meta.seawater, sample_relabel) %>%
    setNames(., c("usvi_metab_mat", "usvi_asv_mat", "meta.seawater", "sample_relabel"))
  readr::write_rds(usvi_dist_mat_list_optA, paste0(projectpath, "/", "usvi_dist_mat_list_optA", ".rds"), compress = "gz")
}

#option B:
#in the microbiome distance matrix, CINAR_BC_39 and CINAR_BC_40 are duplicates
#CINAR_BC_65 and CINAR_BC_66 are duplicates

#in the metabolome, CINAR_BC_39 and CINAR_BC_40 have a Bray-Curtis dissimilarity of 0.04279126
#in the metabolome, CINAR_BC_65 and CINAR_BC_66 have a Bray-Curtis dissimilarity of 0.02544024
#among the non-technical replicate samples, the range of dissimilarities is 0.03306662 to 0.62614110, with a median dissimilarity of 0.24356958
if(file.exists(paste0(projectpath, "/", "usvi_dist_mat_list_optB", ".rds"))){
  usvi_dist_mat_list_optB <- readr::read_rds(paste0(projectpath, "/", "usvi_dist_mat_list_optB", ".rds"))
} else {
  sample_relabel <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::select(metab_deriv_label, sample_id, site, sampling_day, sampling_time) %>%
    dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
    dplyr::left_join(., metadata %>%
                       dplyr::select(sample_id, replicate),
                     by = join_by(sample_id)) %>%
    droplevels %>%
    dplyr::mutate(replA = LETTERS[replicate]) %>%
    dplyr::mutate(repl = seq_len(length(sample_id)), .by = c("sample_id")) %>% droplevels %>%
    dplyr::mutate(replA = dplyr::case_when((repl == 1) ~ replA,
                                           .default = NA)) %>%
    tidyr::fill(replA, .direction = "down") %>%
    dplyr::mutate(repl = dplyr::case_when(length(sample_id) > 1 ~ repl,
                                          .default = NA), .by = "sample_id") %>%
    droplevels %>%
    tidyr::unite("repl", c(replA, repl), sep = "", remove = TRUE, na.rm = TRUE) %>%
    dplyr::select(metab_deriv_label, site, sampling_day, sampling_time, repl) %>%
    dplyr::arrange(site, sampling_time, sampling_day, repl) %>%
    droplevels %>%
    tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time, repl), sep = "_", remove = FALSE)  %>%
    dplyr::distinct(metab_deriv_label, relabeled_sample) %>%
    tibble::deframe(.)
  
  temp_mat <- usvi_metabolomics.df %>%
    dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::relocate(sample_id) %>%
    tidyr::pivot_longer(., cols = !c(sample_id, metab_deriv_label),
                        names_to = "simpleName",
                        values_to = "conc") %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::mutate(conc = log2(conc + 1)) %>%
    tidyr::pivot_wider(., id_cols = "metab_deriv_label",
                       names_from = "simpleName",
                       values_fill = 0,
                       values_from = "conc") %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
    as.matrix(.) %>%
    vegan::vegdist(., binary = FALSE, upper = TRUE,
                   # distance = "horn",
                   distance = "bray",
                   autotransform = TRUE) %>%
    as.matrix(.)
  temp_mat <- temp_mat[intersect(names(sample_relabel), rownames(temp_mat)), intersect(names(sample_relabel), colnames(temp_mat))]
  
  
  
  # dropped_cinar_idx <- c("CINAR_BC_39", "CINAR_BC_40", "CINAR_BC_65", "CINAR_BC_66")
  # temp_mat2 <- temp_mat[c("CINAR_BC_39", "CINAR_BC_40"), c("CINAR_BC_39", "CINAR_BC_40")]
  # temp_mat2[lower.tri(temp_mat2, diag = FALSE)]
  # # [1] 0.04279126
  # 
  # temp_mat2 <- temp_mat[c("CINAR_BC_65", "CINAR_BC_66"), c("CINAR_BC_65", "CINAR_BC_66")]
  # temp_mat2[lower.tri(temp_mat2, diag = FALSE)]
  # # [1] 0.02544024
  # 
  # temp_mat2 <- temp_mat[!(rownames(temp_mat) %in% dropped_cinar_idx), !(colnames(temp_mat) %in% dropped_cinar_idx)]
  # 
  # quantile(temp_mat2[lower.tri(temp_mat2, diag = FALSE)], probs = seq(0, 1, 0.25), names = TRUE)
  # # 0%        25%        50%        75%       100% 
  # # 0.03306662 0.17085251 0.24356958 0.31371609 0.62614110 
  
  temp_mat2 <- usvi_prok_asvs.df %>%
    dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
    dplyr::right_join(., (metabolomics_sample_metadata %>%
                            dplyr::filter(grepl("seawater", sample_type)) %>%
                            dplyr::select(sample_id, metab_deriv_label) %>%
                            droplevels),
                      by = join_by("sample_ID" == "sample_id"), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(metab_deriv_label %in% colnames(temp_mat)) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "asv_id",
                       names_from = "metab_deriv_label",
                       values_from = "counts",
                       values_fill = 0) %>%
    tibble::column_to_rownames(var = "asv_id") %>%
    apply(., 2, relabund) %>%
    as.data.frame(.) %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    tidyr::drop_na(.) %>%
    dplyr::select(colnames(temp_mat)) %>%
    t() %>%
    vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
                   autotransform = TRUE) %>%
    as.matrix(.)
  
  temp_df <- metabolomics_sample_metadata %>%
    dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
    dplyr::filter(metab_deriv_label %in% rownames(temp_mat2)) %>%
    dplyr::select(metab_deriv_label, sample_id, sampling_time, sampling_day, site) %>%
    # dplyr::select(!c(contains("label"), contains("dna_"))) %>%
    tibble::column_to_rownames(., var = "metab_deriv_label") %>%
    droplevels
  
  usvi_dist_mat_list_optB <- list(temp_mat, temp_mat2, temp_df, sample_relabel) %>%
    setNames(., c("usvi_metab_mat", "usvi_asv_mat", "meta.seawater", "sample_relabel"))
  readr::write_rds(usvi_dist_mat_list_optB, paste0(projectpath, "/", "usvi_dist_mat_list_optB", ".rds"), compress = "gz")
}


# Plot heatmaps of the distance matrices ----------------------------------



#get the top 100 genera:
#don't need to do genus-level assessments:
{
  # usvi_top100_genus_mat <- usvi_sw_genus.tbl %>%
  #   dplyr::select(c(asv_id, rownames(usvi_metab_mat))) %>%
  #   # dplyr::filter(sample_ID %in% rownames(usvi_metab_mat)) %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   dplyr::mutate(TotAbund = rowSums(.)) %>%
  #   dplyr::arrange(desc(TotAbund)) %>%
  #   dplyr::slice_head(., n = 100) %>%
  #   dplyr::select(-TotAbund) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   tidyr::pivot_longer(., cols = -c("asv_id"),
  #                       names_to = "sample",
  #                       values_to = "abundance")  %>%
  #   droplevels %>%
  #   # dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
  #   #                                 log2(abundance+1), #log transform abundance (with +1 pseudocount)
  #   #                                 0)) %>%
  #   tidyr::pivot_wider(., id_cols = "sample",
  #                      values_from = "abundance",
  #                      # values_from = "logabund",
  #                      names_from = "asv_id") %>%
  #   dplyr::filter(sample %in% colnames(usvi_metab_mat)) %>%
  #   tibble::column_to_rownames(var = "sample") %>%
  #   tidyr::drop_na(.) %>%
  #   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
  #                  autotransform = TRUE) %>%
  #   as.matrix(.)
  # 
  # usvi_top25_genus_mat <- usvi_sw_genus.tbl %>%
  #   dplyr::select(c(asv_id, rownames(usvi_metab_mat))) %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   dplyr::mutate(TotAbund = rowSums(.)) %>%
  #   dplyr::arrange(desc(TotAbund)) %>%
  #   dplyr::slice_head(., n = 25) %>%
  #   dplyr::select(-TotAbund) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   tidyr::pivot_longer(., cols = -c("asv_id"),
  #                       names_to = "sample",
  #                       values_to = "abundance")  %>%
  #   droplevels %>%
  #   dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
  #                                   log2(abundance+1), #log transform abundance (with +1 pseudocount)
  #                                   0)) %>%
  #   tidyr::pivot_wider(., id_cols = "sample",
  #                      values_from = "abundance",
  #                      # values_from = "logabund",
  #                      names_from = "asv_id") %>%
  #   dplyr::filter(sample %in% colnames(usvi_metab_mat)) %>%
  #   tibble::column_to_rownames(var = "sample") %>%
  #   tidyr::drop_na(.) %>%
  #   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
  #                  autotransform = TRUE) %>%
  #   # as.matrix(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(var = "sample_id") %>% droplevels
  #   as.matrix(.)
  }



#plotting distance matrix of samples using metabolites 

try(f_run_chunk())
if(execchunk) {
  list2env(usvi_dist_mat_list_optB, envir = .GlobalEnv)
  temp_df <- usvi_metab_mat %>%
    as.matrix(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_1") %>%
    tidyr::pivot_longer(., cols = !c(sample_1),
                        names_to = "sample_2",
                        values_to = "dissimilarity") %>%
    dplyr::distinct(., .keep_all = FALSE) %>%
    dplyr::mutate(dissimilarity = dplyr::case_when(sample_1 == sample_2 ~ NA,
                                                   .default = dissimilarity)) %>%
    dplyr::mutate(sample_1 = factor(sample_1, levels = names(sample_relabel)),
                  sample_2 = factor(sample_2, levels = names(sample_relabel))) %>%
    droplevels
  
  # g1 <- print(
  #   ggplot(data = temp_df)
  #   + theme_bw() 
  #   + geom_tile(aes(x = sample_1, y = sample_2, 
  #                   # group = interaction(site, sampling_day, sampling_time),
  #                   fill = (1-dissimilarity)*100),
  #               show.legend = TRUE)
  #   + scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100), transform = "reverse",
  #                          aesthetics = "fill", expand = expansion(1.1,1.1), name = "Similarity (%)",
  #                          na.value = "white")
  #   + scale_x_discrete(labels = sample_relabel, name = "Sample identifier")
  #   + scale_y_discrete(labels = sample_relabel, name = "Sample identifier")
  #   + theme(strip.text.y = element_text(angle = 0),
  #           axis.text.x = element_text(angle = 90), 
  #           panel.grid.minor = element_blank(),
  #           panel.grid.major = element_blank())
  # )
  
  #if you want to plot the similarities in samples using distance matrices of metabolites or genera:
  
    reorder_mat <- intersect(names(sample_relabel), colnames(usvi_metab_mat))
    temp_mat <- usvi_metab_mat[reorder_mat, reorder_mat]
    
    temp_mat[upper.tri(temp_mat, diag = FALSE)] <- NA
    temp_df <- temp_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "sample_1") %>%
      tidyr::pivot_longer(., cols = !c(sample_1),
                          names_to = "sample_2",
                          values_to = "dissimilarity") %>%
      tidyr::drop_na(.) %>%
      dplyr::mutate(input = "metabolites") %>%
      dplyr::ungroup(.) %>%
      # dplyr::mutate(sample_1 = factor(sample_1, levels = names(sample_relabel)),
      #               sample_2 = factor(sample_2, levels = names(sample_relabel))) %>%
      dplyr::arrange(dissimilarity) %>%
      dplyr::mutate(sample_1 = factor(sample_1, levels = unique(.[["sample_1"]])),
                    sample_2 = factor(sample_2, levels = unique(.[["sample_2"]]))) %>%
      droplevels
    
    g1 <- (
      ggplot(data = temp_df %>%
               dplyr::mutate(dissimilarity = dplyr::case_when(sample_1 == sample_2 ~ NA,
                                                              .default = dissimilarity)))
      + theme_bw() 
      + geom_tile(aes(x = sample_1, y = sample_2, fill = (1-dissimilarity)*100, group = input),
                  show.legend = TRUE)
      + scale_fill_gradientn(colors = colorRampPalette(pals::parula(n = 3))(100), 
                             # transform = "reverse", limits = c(100, 0),
                             aesthetics = "fill", 
                             limits = c(0, 100),
                             # expand = expansion(1.1,1.1), 
                             name = "Similarity (%)",
                             na.value = "white")
      + scale_x_discrete(labels = sample_relabel, name = "Sample identifier")
      + scale_y_discrete(labels = sample_relabel, name = "Sample identifier")
      + theme(strip.text.y = element_text(angle = 0),
              axis.text.x = element_text(angle = 90), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    )
    
    
    temp_mat2 <- usvi_asv_mat[rownames(temp_mat), colnames(temp_mat)]
    temp_mat2[upper.tri(temp_mat2, diag = FALSE)] <- NA
    temp_df2 <- temp_mat2 %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "sample_2") %>%
      tidyr::pivot_longer(., cols = !c(sample_2),
                          names_to = "sample_1",
                          values_to = "dissimilarity") %>%
      tidyr::drop_na(.) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(input = "asv") %>%
      dplyr::arrange(dissimilarity) %>%
      dplyr::mutate(sample_1 = factor(sample_1, levels = unique(.[["sample_1"]])),
                    sample_2 = factor(sample_2, levels = unique(.[["sample_2"]]))) %>%
      dplyr::mutate(dissimilarity = dplyr::case_when((sample_2 %in% c("CINAR_BC_39", "CINAR_BC_40") & sample_1 %in% c("CINAR_BC_39", "CINAR_BC_40")) ~ NA,
                                                     (sample_2 %in% c("CINAR_BC_65", "CINAR_BC_66") & sample_1 %in% c("CINAR_BC_65", "CINAR_BC_66")) ~ NA,
                                                     # sample_1 == sample_2 ~ NA,
                                                     .default = dissimilarity)) %>%
      droplevels
    
    g2 <- (
      ggplot(data = temp_df2)
      + theme_bw() 
      + geom_tile(aes(x = sample_1, y = sample_2, fill = (1-dissimilarity)*100, group = input),
                  show.legend = TRUE)
      + scale_fill_gradientn(colors = colorRampPalette(pals::parula(n = 3))(100), 
                             # transform = "reverse", limits = c(100, 0),
                             aesthetics = "fill", 
                             limits = c(0, 100),
                             # expand = expansion(1.1,1.1), 
                             name = "Similarity (%)",
                             na.value = "white")
      + scale_x_discrete(labels = sample_relabel, name = "Sample identifier")
      + scale_y_discrete(labels = sample_relabel, name = "Sample identifier")
      + theme(strip.text.y = element_text(angle = 0),
              axis.text.x = element_text(angle = 90), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    )
    
    temp_df3 <- temp_df2 %>%
      bind_rows(., temp_df) %>%
      dplyr::mutate(input = factor(input, levels = c("asv", "metabolites"))) %>%
      dplyr::rowwise(.) %>%
      dplyr::mutate(dissimilarity = dplyr::case_when(sample_1 == sample_2 ~ NA,
                                                     .default = dissimilarity)) %>%
      dplyr::mutate(sample_1 = factor(sample_1, levels = unique(.[["sample_1"]])),
                    sample_2 = factor(sample_2, levels = unique(.[["sample_2"]]))) %>%
      # dplyr::mutate(sample_1 = factor(sample_1, levels = names(sample_relabel)),
      #               sample_2 = factor(sample_2, levels = names(sample_relabel))) %>%
      dplyr::arrange(input, sample_1, sample_2) %>%
      droplevels
    
    g3 <- (
      ggplot(data = temp_df3)
      + theme_bw() 
      + geom_tile(aes(x = sample_1, y = sample_2, fill = (1-dissimilarity)*100, group = input),
                  show.legend = TRUE)
      +  scale_fill_gradientn(colors = colorRampPalette(pals::parula(n = 3))(100), 
                              # transform = "reverse", limits = c(100, 0),
                              aesthetics = "fill", 
                              limits = c(0, 100),
                              # expand = expansion(1.1,1.1), 
                              name = "Similarity (%)",
                              na.value = "white")
      + scale_x_discrete(labels = sample_relabel,
                         # name = "Sample identifier")
                         name = "Metabolomics")
      + scale_y_discrete(labels = sample_relabel, 
                         # name = "Sample identifier")
                         name = "Dissimilarity by relative abundances of ASVs")
      + theme(strip.text.y = element_text(angle = 0),
              axis.text.x = element_text(angle = 90), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    )
  
  g4 <- g3 + facet_wrap(.~input) + scale_x_discrete(labels = sample_relabel, name = "Sample identifier") + scale_y_discrete(labels = sample_relabel, name = "Sample identifier")
  g3 <- g3 + patchwork::plot_annotation(title = "Similarity in profiles generated from ASV-level abundances, and metabolomics")
  g4 <- g4 + patchwork::plot_annotation(title = "Similarity in profiles generated from ASV-level abundances, and metabolomics")
  
  g5 <- (g1 + ggtitle("Metabolomics")) + (g2 + ggtitle("ASVs")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Similarity in profiles generated from ASV-level abundances, and metabolomics", tag_levels = "A")


  
  # if(!any(grepl("dissim_metab_asv-", list.files(projectpath, pattern = "usvi_.*.png")))){
  #   ggsave(paste0(projectpath, "/", "usvi_dissim_metab_asv-", Sys.Date(), ".png"),
  #          g5,
  #          width = 20,height = 10, units = "in")
  #   ggsave(paste0(projectpath, "/", "usvi_dissim_metab_asv_combo-", Sys.Date(), ".png"),
  #          g3,
  #          width = 10,height = 10, units = "in")
  # }
  
  
}


#trying out mrpp--in progress
{
  # temp_mat_metadata <- metabolomics_sample_metadata %>%
  #   dplyr::select(sample_id, site, sampling_day, sampling_time) %>%
  #   dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  #   droplevels %>%
  #   dplyr::ungroup(.) %>%
  #   # dplyr::mutate(across(c(site, sampling_day, sampling_time), ~as.character(.x))) %>%
  #   # dplyr::mutate(across(c(site, sampling_day, sampling_time), ~as.numeric(.x))) %>%
  #   tibble::column_to_rownames(var = "sample_id")
  #   # dplyr::select(sample_id, site) %>%
  #   # tibble::deframe(.)
  # 
  # metab_mrpp <- vegan::mrpp(usvi_metab_mat, temp_mat_metadata$site)
}

#mantel testing-- in progress
{
  
  # usvi_asv_df <- ps_usvi %>%
  #   phyloseq::subset_samples(., sample_type == "seawater") %>%
  #   phyloseq::otu_table(.) %>%
  #   t() %>%
  #   as.data.frame %>%
  #   tibble::rownames_to_column(var = "sample_id") %>%
  #   dplyr::left_join(., (metabolomics_sample_metadata %>%
  #                          dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  #                          dplyr::filter(grepl("seawater", sample_type)) %>%
  #                          dplyr::select(sample_id, metab_deriv_label) %>%
  #                          droplevels),
  #                    by = join_by(sample_id)) %>%
  #     dplyr::select(-sample_id) %>%
  #     tibble::column_to_rownames(var = "metab_deriv_label") %>%
  #   t() %>%
  #   apply(., 2, relabund) %>% 
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   tidyr::pivot_longer(., cols = -c("asv_id"),
  #                       names_to = "sample",
  #                       values_to = "abundance")
  # usvi_asv_mat <- usvi_asv_df %>%
  #   dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
  #                                   log2(abundance+1), #log transform abundance (with +1 pseudocount)
  #                                   0)) %>%
  #   tidyr::pivot_wider(., id_cols = "sample",
  #                      # values_from = "abundance",
  #                      values_from = "logabund",
  #                      names_from = "asv_id") %>%
  #   dplyr::filter(sample %in% colnames(usvi_metab_mat)) %>%
  #   tibble::column_to_rownames(var = "sample") %>%
  #   tidyr::drop_na(.) %>%
  #   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
  #                  autotransform = TRUE) %>%
  #   # as.matrix(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(var = "sample_id") %>% droplevels
  #   as.matrix(.)
  
  # usvi_mantel_res <- vegan::mantel(usvi_asv_mat, usvi_metab_mat, permutations = 999, parallel = nthreads,
  #                                  # method = "pearson") 
  #                                  method = "spearman")
  
  #at the ASV level,
  # Call:
  #   vegan::mantel(xdis = usvi_asv_df, ydis = usvi_metab_df, method = "spearman",      permutations = 999, parallel = nthreads) 
  # 
  # Mantel statistic r: -0.04485 
  # Significance: 0.786 
  
  
  #what about at the genus level?
  # usvi_genus_df <- ps_usvi %>%
  #   phyloseq::subset_samples(., sample_type == "seawater") %>%
  #   phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% # remove ASVs not present in any samples
  #   phyloseq::tax_glom(., taxrank = "Genus") %>%
  #   phyloseq::otu_table(.) %>%
  #   t() %>%
  #   as.data.frame %>%
  #   tibble::rownames_to_column(var = "sample_id") %>%
  #   # dplyr::filter(sample_ID %in% metabolomics_sample_metadata[["sample_id"]]) %>%
  #   dplyr::left_join(., (metabolomics_sample_metadata %>%
  #                          dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  #                          dplyr::filter(grepl("seawater", sample_type)) %>%
  #                          dplyr::select(sample_id, metab_deriv_label) %>%
  #                          droplevels),
  #                    by = join_by(sample_id)) %>%
  #   dplyr::select(-sample_id) %>%
  #   tibble::column_to_rownames(var = "metab_deriv_label") %>%
  #   t() %>%
  #   apply(., 2, relabund) %>% 
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   # tidyr::pivot_longer(., cols = -c("asv_id"),
  #   #                     names_to = "sample",
  #   #                     values_to = "abundance")  %>%
  #   droplevels
  
  # usvi_genus_mat <- usvi_genus_df %>%
  #   dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
  #                                   log2(abundance+1), #log transform abundance (with +1 pseudocount)
  #                                   0)) %>%
  #   tidyr::pivot_wider(., id_cols = "sample",
  #                      # values_from = "abundance",
  #                      values_from = "logabund",
  #                      names_from = "asv_id") %>%
  #   dplyr::filter(sample %in% colnames(usvi_metab_mat)) %>%
  #   tibble::column_to_rownames(var = "sample") %>%
  #   tidyr::drop_na(.) %>%
  #   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
  #                  autotransform = TRUE) %>%
  #   # as.matrix(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(var = "sample_id") %>% droplevels
  #   as.matrix(.)
  # 
  # usvi_mantel_genus_res <- vegan::mantel(usvi_genus_mat, usvi_metab_mat, permutations = 999, parallel = nthreads,
  #                                  method = "pearson")
  #                                  # method = "spearman")
  # Call:
  #   vegan::mantel(xdis = usvi_genus_df, ydis = usvi_metab_df, method = "pearson",      permutations = 999, parallel = nthreads) 
  # 
  # Mantel statistic r: -0.06304 
  # Significance: 0.846 
  
  
  # 
  # 
  # usvi_top_genus_mat <- usvi_top100_genus_mat
  # # usvi_top_genus_mat <- usvi_top25_genus_mat
  # usvi_mantel_genus_res <- vegan::mantel(usvi_top_genus_mat, usvi_metab_mat, permutations = 999, parallel = nthreads,
  #                                        # method = "pearson")
  #                                        method = "spearman")
  # usvi_mantel_genus_res
  
}



#16S profiles of Tektite samples ranged in simialrity between 62.1% to 90.8% similar
#between LB and Yawzi, 34.7% to 84.3% similar
#and between LB and Tektite, 25.9% to 72.3% similar
{
  # temp_mat <- usvi_asv_mat
  # # temp_mat <- usvi_top100_genus_mat
  # temp_mat <- temp_mat[rownames(temp_mat) %in% names(sample_relabel)[grepl("Tektite", sample_relabel)], colnames(temp_mat) %in% names(sample_relabel)[grepl("Tektite", sample_relabel)]]
  # 1-range(temp_mat[lower.tri(temp_mat, diag = FALSE)])
  # 
  # temp_mat <- usvi_asv_mat
  # # temp_mat <- usvi_top100_genus_mat
  # temp_mat <- temp_mat[rownames(temp_mat) %in% names(sample_relabel)[grepl("Tektite", sample_relabel)], colnames(temp_mat) %in% names(sample_relabel)[grepl("LB", sample_relabel)]]
  # 1-range(temp_mat[lower.tri(temp_mat, diag = FALSE)])
  # 
  # temp_mat <- usvi_asv_mat
  # # temp_mat <- usvi_top100_genus_mat
  # temp_mat <- temp_mat[rownames(temp_mat) %in% names(sample_relabel)[grepl("Yawzi", sample_relabel)], colnames(temp_mat) %in% names(sample_relabel)[grepl("LB", sample_relabel)]]
  # 1-range(temp_mat[lower.tri(temp_mat, diag = FALSE)])
}

#on the other hand, metabolomics profiles of Tektite samples ranged from 40.0% to 95.5% similar
#and between LB and Tektite, 52.1% to 86.2% similar
{
  # temp_mat <- usvi_metab_mat
  # temp_mat <- temp_mat[rownames(temp_mat) %in% names(sample_relabel)[grepl("Tektite", sample_relabel)], colnames(temp_mat) %in% names(sample_relabel)[grepl("Tektite", sample_relabel)]]
  # 1-range(temp_mat[lower.tri(temp_mat, diag = FALSE)])
  # 
  # temp_mat <- usvi_metab_mat
  # temp_mat <- temp_mat[rownames(temp_mat) %in% names(sample_relabel)[grepl("LB", sample_relabel)], colnames(temp_mat) %in% names(sample_relabel)[grepl("Tektite", sample_relabel)]]
  # 1-range(temp_mat[lower.tri(temp_mat, diag = FALSE)])
}





#recaling the concentrations of metabolits from 0 to 1 wrt metabolites
#using top 25 genera relative abundance:
# Call:
#   vegan::mantel(xdis = usvi_top_genus_mat, ydis = usvi_metab_mat,      method = "spearman", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: -0.006939 
# Significance: 0.551 


#calculating relative abundance wtihin sample of each metabolite, then rescaling from 0 to 1:
#using top 25 genera relative abundance:
# vegan::mantel(xdis = usvi_top_genus_mat, ydis = usvi_metab_mat,      method = "spearman", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: 0.04983 
# Significance: 0.115 

#calculating relative abundance wtihin sample of each metabolite, then rescaling from 0 to 1:
#using top 25 genera relative abundance log+1-transformed:
# vegan::mantel(xdis = usvi_top_genus_mat, ydis = usvi_metab_mat,      method = "spearman", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: 0.04983 
# Significance: 0.112 



#calculating relative abundance wtihin sample of each metabolite, then rescaling from 0 to 1:
#using top 100 genera rlative abundance:
# vegan::mantel(xdis = usvi_top_genus_mat, ydis = usvi_metab_mat,      method = "spearman", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: 0.05378 
# Significance: 0.095 


#calculating relative abundance wtihin sample of each metabolite, then rescaling from 0 to 1:
#using top 100 genera rlative abundance log+1 transformed:
# vegan::mantel(xdis = usvi_top_genus_mat, ydis = usvi_metab_mat,      method = "spearman", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: 0.05378 
# Significance: 0.105 


# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# vegan::mantel(xdis = usvi_top100_genus_mat, ydis = usvi_metab_mat,      method = "pearson", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: 0.0172 
#       Significance: 0.319 

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# vegan::mantel(xdis = usvi_top100_genus_mat, ydis = usvi_metab_mat,      method = "spearman", permutations = 999, parallel = nthreads) 
# 
# Mantel statistic r: 0.02933 
#       Significance: 0.275 
# 



# Spearman rank correlation coefficient calculations ----------------------


#option B:
if(file.exists(paste0(projectpath, "/", "spearman.test.optB.list", ".rds"))){
  spearman.test.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.optB.list", ".rds"))
} else { 
  list2env(usvi_dist_mat_list_optB, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics.df %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::relocate(sample_id) %>%
    tidyr::pivot_longer(., cols = !c(sample_id, metab_deriv_label),
                        names_to = "simpleName",
                        values_to = "conc") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(metab_deriv_label, simpleName) %>%
    # dplyr::summarise(num = length(conc)) %>%
    dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE), 
                     num = length(conc),
                     # .by = c(sample_id, simpleName),
                     .groups = "keep",
                     sd = sd(conc, na.rm = TRUE)) %>%
    dplyr::rename(conc = "mean_conc") %>%
    dplyr::mutate(log_conc = ifelse(!(is.na(conc) | (conc < 0)),
                                    log2(conc+1), #log transform abundance (with +1 pseudocount)
                                    0)) %>%
    dplyr::select(metab_deriv_label, simpleName, log_conc) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "metab_deriv_label",
                       values_from = "log_conc",
                       # values_from = "logabund",
                       names_from = "simpleName") %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
    droplevels
  
  usvi_asv.tbl <- usvi_prok_asvs.df %>%
    dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
    dplyr::right_join(., (metabolomics_sample_metadata %>%
                            dplyr::filter(grepl("seawater", sample_type)) %>%
                            dplyr::select(sample_id, metab_deriv_label) %>%
                            droplevels),
                      by = join_by("sample_ID" == "sample_id"), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(metab_deriv_label %in% rownames(usvi_metab.tbl)) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "asv_id",
                       names_from = "metab_deriv_label",
                       values_from = "counts",
                       values_fill = 0) %>%
    tibble::column_to_rownames(var = "asv_id") %>%
    apply(., 2, relabund) %>%
    as.data.frame(.) %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    tidyr::drop_na(.) %>%
    dplyr::select(rownames(usvi_metab.tbl)) %>%
    t()
  
  usvi_asv.tbl <- usvi_asv.tbl[rownames(usvi_metab.tbl),]
  
  spearman.test <- matrix(nrow = ncol(usvi_asv.tbl), ncol = ncol(usvi_metab.tbl))
  colnames(spearman.test) <- colnames(usvi_metab.tbl)
  rownames(spearman.test) <- colnames(usvi_asv.tbl)
  
  spearman.test.rho <- spearman.test
  
  y <- length(colnames(spearman.test))
  for(j in seq_len(y)){
    vector_metab <- usvi_metab.tbl[, j]
    # for(i in seq_len(2)){
    for(i in seq_len(nrow(spearman.test))){
      spearman.test[i, j] <- cor.test(vector_metab, usvi_asv.tbl[,i], method = "spearman", exact = FALSE) %>%
        purrr::pluck(., "p.value")
      spearman.test.rho[i, j] <- cor.test(vector_metab, usvi_asv.tbl[,i], method = "spearman", exact = FALSE) %>%
        purrr::pluck(., "estimate")
    }
  }
  
  
  #for a q-value of 0.05, the adjusted p threshold should be 0.0206
  #and we have X significant relatinships
  #For a q-value of 0.1, the adjusted p-threshold should be 0.0579
  #and we have Y significant relationships
  
  q_value <- 0.05
  padj_cutoff <- spearman.test %>%
    # apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
    apply(., 2, na.omit) %>%
    apply(., 2, function(x) ashr::qval.from.lfdr(x)) %>%
    unlist %>% as.matrix(.) %>%
    # quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)
    quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
  
  spearman.test.bh.corrected <- spearman.test %>%
    apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
    apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  
  # spearman.test.0.05.corrected <- spearman.test %>%
  #   apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
  #   # apply(., 2, function(x) ifelse(x < padj_cutoff, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  #   apply(., 2, function(x) ifelse(x < padj_cutoff[1], x, NA)) #drop the p.values > the adjusted p-value or did not compute
  # spearman.test.0.10.corrected <- spearman.test %>%
  #   apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
  #   # apply(., 2, function(x) ifelse(x < padj_cutoff, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  #   apply(., 2, function(x) ifelse(x < padj_cutoff[2], x, NA)) #drop the p.values > the adjusted p-value or did not compute
  # 
  # length(which(!is.na(spearman.test.0.10.corrected))) #7298
  # length(which(!is.na(spearman.test.0.05.corrected))) #5557
  # length(which(!is.na(spearman.test.bh.corrected))) #6503
  
  # temp_mat <- spearman.test.0.10.corrected[spearman.test.bh.corrected] %>%
  #   as.matrix(.)
  # 
  dend_asv <- spearman.test.rho %>%
    dist(t(.), method = "euclidean") %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram
  dend_metab <- spearman.test.rho %>%
    t() %>%
    dist(t(.), method = "euclidean") %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram
  
  spearman.test.df <- spearman.test.bh.corrected %>%
    tibble::as_tibble(., rownames = "asv_id") %>%
    tidyr::pivot_longer(., cols = !asv_id,
                        names_to = "simpleName",
                        values_to = "padj_bh") %>%
    droplevels %>%
    dplyr::left_join(., (spearman.test %>%
                           tibble::as_tibble(., rownames = "asv_id") %>%
                           tidyr::pivot_longer(., cols = !asv_id,
                                               names_to = "simpleName",
                                               values_to = "padj") %>%
                           tidyr::drop_na(.)),
                     by = join_by(asv_id, simpleName)) %>%
    dplyr::relocate(padj_bh, .after = "padj") %>%
    tidyr::drop_na(padj_bh) %>%
    dplyr::mutate(padj_05 = dplyr::case_when(padj_bh <= padj_cutoff[1] ~ padj_bh,
                                             .default = NA)) %>%
    dplyr::right_join(., (spearman.test.rho %>%
                            tibble::as_tibble(., rownames = "asv_id") %>%
                            tidyr::pivot_longer(., cols = !asv_id,
                                                names_to = "simpleName",
                                                values_to = "estimate")),
                      by = join_by(asv_id, simpleName)) %>%
    dplyr::mutate(sig = dplyr::case_when(!is.na(padj_bh) & !is.na(padj_05) ~ "sig",
                                         !is.na(padj_bh) & is.na(padj_05) ~ "maybe",
                                         .default = "not")) %>%
    dplyr::mutate(sig = factor(sig)) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
    dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
    dplyr::arrange(asv_id, simpleName) %>%
    dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
    # tidyr::drop_na(contains("padj")) %>%
    dplyr::mutate(label = signif(estimate, digits = 2)) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
    droplevels
  
  spearman.test.df %>%
    tidyr::drop_na(padj_05) %>%
    # dplyr::summarise(num_results = length(padj_bh), .by = "simpleName") %>%
    dplyr::summarise(num_results = length(padj_bh), .by = "asv_id") %>%
    dplyr::arrange(desc(num_results)) %>%
    dplyr::select(num_results) %>%
    tibble::deframe(.) %>%
    quantile(., probs = seq(0, 1, 0.25), names = TRUE)
  
  #passing the q-value cutoff of 0.05:
  #of the 49 metabolites, the maximum number of significantly correlated ASVs was 652, the median was 142
  #of the 1664 ASVs that had significant correlations with metabolite concentrations,
  #the maximum number of significantly correlated metabolites for an ASV was 24.  the median was 1.
  
  #using the BH adjusted p-values < 0.05:
  #of the 49 metabolites, the maximum number of significantly correlated ASVs was 658, the median was 154.
  #of the 1742 ASVs that had significant correlations with metabolite concentrations,
  #the maximum number of significantly correlated metabolites for an ASV was 27.  the median was 1.
  
  #using the q-value cutoff, we lose 78 ASVs that had an adjusted p-value < 0.05 but not below the 0.05 q-value cutoff.
  
  spearman.test.optB.list <- list(spearman.test.df, dend_asv, dend_metab, padj_cutoff) %>%
    setNames(., c("spearman.test.df", "dend_asv", "dend_metab", "padj_cutoff"))
  readr::write_rds(spearman.test.optB.list, paste0(projectpath, "/", "spearman.test.optB.list", ".rds"), compress = "gz")
}

#option A:
if(file.exists(paste0(projectpath, "/", "spearman.test.optA.list", ".rds"))){
  spearman.test.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.optA.list", ".rds"))
} else {
  list2env(usvi_dist_mat_list_optA, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics.df %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::relocate(sample_id) %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
    tidyr::pivot_longer(., cols = !c(sample_id, metab_deriv_label),
                        names_to = "simpleName",
                        values_to = "conc") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(sample_id, simpleName) %>%
    # dplyr::summarise(num = length(conc)) %>%
    dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE), 
                     num = length(conc),
                     # .by = c(sample_id, simpleName),
                     .groups = "keep",
                     sd = sd(conc, na.rm = TRUE)) %>%
    dplyr::rename(conc = "mean_conc") %>%
    dplyr::mutate(log_conc = ifelse(!(is.na(conc) | (conc < 0)),
                                    log2(conc+1), #log transform abundance (with +1 pseudocount)
                                    0)) %>%
    dplyr::select(sample_id, simpleName, log_conc) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "sample_id",
                       values_from = "log_conc",
                       # values_from = "logabund",
                       names_from = "simpleName") %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    droplevels
  
  usvi_asv.tbl <- usvi_prok_asvs.df %>%
    dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
    # dplyr::right_join(., (metabolomics_sample_metadata %>%
    #                         dplyr::filter(grepl("seawater", sample_type)) %>%
    #                         dplyr::select(sample_id, metab_deriv_label) %>%
    #                         droplevels),
    #                   by = join_by("sample_ID" == "sample_id"), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(sample_ID %in% rownames(usvi_metab.tbl)) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "asv_id",
                       names_from = "sample_ID",
                       values_from = "counts",
                       values_fill = 0) %>%
    tibble::column_to_rownames(var = "asv_id") %>%
    apply(., 2, relabund) %>%
    as.data.frame(.) %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    tidyr::drop_na(.) %>%
    dplyr::select(rownames(usvi_metab.tbl)) %>%
    t()
  
  usvi_asv.tbl <- usvi_asv.tbl[rownames(usvi_metab.tbl),]
  
  spearman.test <- matrix(nrow = ncol(usvi_asv.tbl), ncol = ncol(usvi_metab.tbl))
  colnames(spearman.test) <- colnames(usvi_metab.tbl)
  rownames(spearman.test) <- colnames(usvi_asv.tbl)
  
  spearman.test.rho <- spearman.test
  
  y <- length(colnames(spearman.test))
  for(j in seq_len(y)){
    vector_metab <- usvi_metab.tbl[, j]
    # for(i in seq_len(2)){
    for(i in seq_len(nrow(spearman.test))){
      spearman.test[i, j] <- cor.test(vector_metab, usvi_asv.tbl[,i], method = "spearman", exact = FALSE) %>%
        purrr::pluck(., "p.value")
      spearman.test.rho[i, j] <- cor.test(vector_metab, usvi_asv.tbl[,i], method = "spearman", exact = FALSE) %>%
        purrr::pluck(., "estimate")
    }
  }
  
  q_value <- 0.05
  padj_cutoff <- spearman.test %>%
    apply(., 2, na.omit) %>%
    apply(., 2, function(x) ashr::qval.from.lfdr(x)) %>%
    unlist %>% as.matrix(.) %>%
    quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
  
  spearman.test.bh.corrected <- spearman.test %>%
    apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
    apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute

  # spearman.test.0.05.corrected <- spearman.test %>%
  #   apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
  #   # apply(., 2, function(x) ifelse(x < padj_cutoff, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  #   apply(., 2, function(x) ifelse(x < padj_cutoff[1], x, NA)) #drop the p.values > the adjusted p-value or did not compute
  # spearman.test.0.10.corrected <- spearman.test %>%
  #   apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
  #   # apply(., 2, function(x) ifelse(x < padj_cutoff, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  #   apply(., 2, function(x) ifelse(x < padj_cutoff[2], x, NA)) #drop the p.values > the adjusted p-value or did not compute
  # 
  # length(which(!is.na(spearman.test.0.10.corrected))) #6645
  # length(which(!is.na(spearman.test.0.05.corrected))) #5558
  # length(which(!is.na(spearman.test.bh.corrected))) #6429
  
  dend_asv <- spearman.test.rho %>%
    dist(t(.), method = "euclidean") %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram
  dend_metab <- spearman.test.rho %>%
    t() %>%
    dist(t(.), method = "euclidean") %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram
  
  spearman.test.df <- spearman.test.bh.corrected %>%
    tibble::as_tibble(., rownames = "asv_id") %>%
    tidyr::pivot_longer(., cols = !asv_id,
                        names_to = "simpleName",
                        values_to = "padj_bh") %>%
    droplevels %>%
    dplyr::left_join(., (spearman.test %>%
                           tibble::as_tibble(., rownames = "asv_id") %>%
                           tidyr::pivot_longer(., cols = !asv_id,
                                               names_to = "simpleName",
                                               values_to = "padj") %>%
                           tidyr::drop_na(.)),
                     by = join_by(asv_id, simpleName)) %>%
    dplyr::relocate(padj_bh, .after = "padj") %>%
    tidyr::drop_na(padj_bh) %>%
    dplyr::mutate(padj_05 = dplyr::case_when(padj_bh <= padj_cutoff[1] ~ padj_bh,
                                             .default = NA)) %>%
    # dplyr::mutate(padj_10 = dplyr::case_when(padj_bh <= padj_cutoff[2] ~ padj_bh,
    #                                          .default = NA)) %>%
    # dplyr::full_join(., (spearman.test.0.10.corrected %>%
    #                        tibble::as_tibble(., rownames = "asv_id") %>%
    #                        tidyr::pivot_longer(., cols = !asv_id,
    #                                            names_to = "simpleName",
    #                                            values_to = "padj_10")),
    #                  by = join_by(asv_id, simpleName)) %>%
    # dplyr::full_join(., (spearman.test.0.05.corrected %>%
    #                        tibble::as_tibble(., rownames = "asv_id") %>%
    #                        tidyr::pivot_longer(., cols = !asv_id,
    #                                            names_to = "simpleName",
    #                                            values_to = "padj_05")),
    #                  by = join_by(asv_id, simpleName)) %>%
    dplyr::right_join(., (spearman.test.rho %>%
                            tibble::as_tibble(., rownames = "asv_id") %>%
                            tidyr::pivot_longer(., cols = !asv_id,
                                                names_to = "simpleName",
                                                values_to = "estimate")),
                      by = join_by(asv_id, simpleName)) %>%
    dplyr::mutate(sig = dplyr::case_when(!is.na(padj_bh) & !is.na(padj_05) ~ "sig",
                                         !is.na(padj_bh) & is.na(padj_05) ~ "maybe",
                                         .default = "not")) %>%
    dplyr::mutate(sig = factor(sig)) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
    dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
    dplyr::arrange(asv_id, simpleName) %>%
    dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
    # tidyr::drop_na(contains("padj")) %>%
    dplyr::mutate(label = signif(estimate, digits = 2)) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
    droplevels
  
  length(unique((spearman.test.df %>% tidyr::drop_na(padj_05))$asv_id))
  length(unique((spearman.test.df %>% tidyr::drop_na(padj_bh))$asv_id))
  
  spearman.test.df %>%
    # tidyr::drop_na(padj_05) %>%
    # dplyr::summarise(num_results = length(padj_bh), .by = "simpleName") %>%
    dplyr::summarise(num_results = length(padj_bh), .by = "asv_id") %>%
    dplyr::arrange(desc(num_results)) %>%
    dplyr::select(num_results) %>%
    tibble::deframe(.) %>%
    quantile(., probs = seq(0, 1, 0.25), names = TRUE)
  
  #passing the q-value cutoff of 0.05:
  #of the 49 metabolites, the maximum number of significantly correlated ASVs was 651, the median was 144
  #of the 1660 ASVs that had significant correlations with metabolite concentrations,
  #the maximum number of significantly correlated metabolites for an ASV was 24.  the median was 1.
  
  #using the BH adjusted p-values < 0.05:
  #of the 49 metabolites, the maximum number of significantly correlated ASVs was 658, the median was 172
  #of the 1734 ASVs that had significant correlations with metabolite concentrations,
  #the maximum number of significantly correlated metabolites for an ASV was 27.  the median was 1.
  
  #using the q-value cutoff, we lose 74 ASVs that had an adjusted p-value < 0.05 but not below the 0.05 q-value cutoff.
  
  spearman.test.optA.list <- list(spearman.test.df, dend_asv, dend_metab, padj_cutoff) %>%
    setNames(., c("spearman.test.df", "dend_asv", "dend_metab", "padj_cutoff"))
  readr::write_rds(spearman.test.optA.list, paste0(projectpath, "/", "spearman.test.optA.list", ".rds"), compress = "gz")
}

# Filter the spearman correlations ----------------------------------------


# setdiff(unique((spearman.test.optA.list[["spearman.test.df"]] %>% tidyr::drop_na(padj_bh))$asv_id),
#         unique((spearman.test.optB.list[["spearman.test.df"]] %>% tidyr::drop_na(padj_bh))$asv_id))
# setdiff(unique((spearman.test.optB.list[["spearman.test.df"]] %>% tidyr::drop_na(padj_bh))$asv_id),
#         unique((spearman.test.optA.list[["spearman.test.df"]] %>% tidyr::drop_na(padj_bh))$asv_id)) %in% unique((spearman.test.optB.list[["spearman.test.df"]] %>% tidyr::drop_na(padj_bh))$asv_id)
# 
# #12 ASVs are in spearman.test.optB.list that are not in spearman.test.optA.list
# #3 ASVs are in spearman.test.optA.list that are not in spearman.test.optB.list

spearman.test.df <- list(spearman.test.optA.list[["spearman.test.df"]],
                         spearman.test.optB.list[["spearman.test.df"]]) %>%
  setNames(., c("optA", "optB")) %>%
  bind_rows(., .id = "test_type")
# temp_df <- spearman.test.df %>%
#   tidyr::drop_na(padj_bh) %>%
#   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
#   dplyr::filter(consistent > 1) %>%
#   dplyr::select(test_type, asv_id, simpleName, padj_05) %>%
#   tidyr::pivot_wider(., id_cols = NULL,
#                      names_from = "test_type",
#                      values_from = "padj_05") %>%
#   droplevels
# 
# temp_df %>% tidyr::drop_na(optA) %>% droplevels %>% nrow(.)
# #5558 correlations
# temp_df %>% tidyr::drop_na(optB) %>% droplevels %>% nrow(.)
# #5557 correlations
# temp_df %>% tidyr::drop_na(.) %>% droplevels %>% nrow(.)
# #5468 correlations that are consistent between options A and B

#filter for the correlations between ASVs and metabolites, that are consistently significant, in the same direction, and observed in both options A and B
temp_df <- spearman.test.df %>%
  tidyr::drop_na(padj_bh) %>%
  dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
  dplyr::filter(consistent > 1) %>%
  dplyr::select(test_type, asv_id, simpleName, estimate) %>%
  tidyr::pivot_wider(., id_cols = NULL,
                     names_from = "test_type",
                     values_from = "estimate") %>%
  dplyr::mutate(consistent = dplyr::case_when((optA * optB) > 0 ~ 1,
                                              .default = NA)) %>%
  tidyr::drop_na(.) %>%
  dplyr::select(asv_id, simpleName) %>%
  dplyr::inner_join(., (spearman.test.df %>%
                         tidyr::drop_na(padj_bh) %>%
                         dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
                         dplyr::filter(consistent > 1) %>%
                         dplyr::select(test_type, asv_id, simpleName, padj_05) %>%
                         tidyr::pivot_wider(., id_cols = NULL,
                                            names_from = "test_type",
                                            values_from = "padj_05") %>%
                         tidyr::drop_na(.) %>%
                         droplevels %>%
                         dplyr::select(asv_id, simpleName)),
                   by = join_by(asv_id, simpleName)) %>%
  dplyr::left_join(., spearman.test.df %>%
                     dplyr::select(test_type, asv_id, simpleName, estimate, sig, label) %>%
                     droplevels,
                   by = join_by(asv_id, simpleName), relationship = "many-to-many", multiple = "all") %>%
  dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                     .default = NA)) %>%
  dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
  droplevels


# output list of SDA/sig ASVs ---------------------------------------------


#write out the list of ASVs that were significant, to make a tree:
#SDA ASVs:
temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_sda_asvs_compare_summary-.*.tsv"))
usvi_sda_asvs_compare_summary.df <- readr::read_delim(paste0(projectpath, "/", temp_file),
                                                      delim = "\t", col_names = TRUE, na = "", show_col_types = FALSE)
rm(temp_file)


usvi_sig_seqs_idx <- temp_df %>%  
  # dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
  # tidyr::drop_na(estimate) %>%
  # dplyr::summarise(num_results = length(estimate), .by = "asv_id") %>%
  # dplyr::filter(num_results > 1) %>%
  dplyr::distinct(asv_id, .keep_all = FALSE) %>%
  droplevels %>%
  dplyr::bind_rows(., usvi_sda_asvs_compare_summary.df %>%
                     dplyr::filter(test_type == "all") %>%
                     droplevels %>%
                     dplyr::distinct(asv_id, .keep_all = FALSE)) %>%
  dplyr::distinct(asv_id) %>%
  droplevels

usvi_sig_seqs_key.list <- usvi_prok_asvs.taxa %>%
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
  dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
  dplyr::filter(asv_id %in% usvi_sig_seqs_idx[["asv_id"]]) %>%
  dplyr::select(asv_id, taxonomy, Domain:Genus) %>%
  droplevels %>%
  # tidyr::nest(., .by = c("taxonomy", Domain:Genus), .key = "asvs") %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  dplyr::left_join(., usvi_prok_asvs.taxa, by = join_by(asv_id)) %>%
  droplevels %>%
  dplyr::select(asv_id, sequence) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::mutate(abundance = 1) %>%
  droplevels



library(dada2)
dada2::uniquesToFasta(usvi_sig_seqs_key.list, paste0(projectpath, "/", "usvi_prok_sig_asvs.fna"), ids = usvi_sig_seqs_key.list[["asv_id"]], mode = "w")

# temp_df %>% dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%

##look at only the strongest correlations
#the minimum abs(estimate) is 0.3929270

#645 ASVs have 2 or more significant correlations
#251 ASVs have 2 or more significant correlations where the abs(estimate) >= 0.5
#225 ASVs have 3 or more significant correlations above the estimate threshold
#4 or more: 210 ASVs, and there are 3-206 significant correlations with metabolites
#5 or more: 192 ASVs
# temp_df %>%
#   dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
#   tidyr::drop_na(filtered_estimate) %>%
#   dplyr::summarise(num_results = length(filtered_estimate), .by = "asv_id") %>%
#   dplyr::filter(num_results > 3) %>%
#   dplyr::select(num_results) %>%
#   tibble::deframe(.) %>% 
#   length(.)
# sig_corr_asvs_idx <- temp_df %>%  
#   dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
#   # tidyr::drop_na(estimate) %>%
#   # dplyr::summarise(num_results = length(estimate), .by = "asv_id") %>%
#   tidyr::drop_na(filtered_estimate) %>%
#   dplyr::summarise(num_results = length(filtered_estimate), .by = "asv_id") %>%
#   dplyr::filter(num_results > 4) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(asv_id) %>%
#   droplevels %>%
#   tibble::deframe(.)
# spearman.test.filtered.df <- temp_df %>%
#   dplyr::filter(asv_id %in% sig_corr_asvs_idx) %>%
#   # tidyr::drop_na(filtered_estimate) %>%
#   # dplyr::select(-filtered_estimate) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
#   # dplyr::arrange(asv_id, simpleName) %>%
#   droplevels
# spearman.test.filtered.df %>%
#   dplyr::summarise(num_results = length(estimate), .by = "simpleName") %>%
#   dplyr::arrange(desc(num_results)) %>%
#   droplevels

spearman.sig.metab.list <- temp_df %>%  
  dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
  tidyr::drop_na(estimate) %>%
  dplyr::select(simpleName, asv_id) %>%
  dplyr::mutate(simpleName = factor(simpleName, levels = (temp_df %>%  
                                                            dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
                                                            tidyr::drop_na(estimate) %>%
                                                            dplyr::summarise(num_results = length(estimate), .by = "simpleName") %>%
                                                            dplyr::arrange(desc(num_results)) %>%
                                                            dplyr::select(simpleName) %>%
                                                            tibble::deframe(.)))) %>%
  split(., f = .$simpleName) %>%
  map(., ~.x %>%
        dplyr::distinct(asv_id, .keep_all = FALSE) %>%
        unlist %>%
        as.character)

spearman.sig.strong.metab.list <- temp_df %>%  
  dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
  tidyr::drop_na(filtered_estimate) %>%
  dplyr::select(simpleName, asv_id) %>%
  dplyr::mutate(simpleName = factor(simpleName, levels = (temp_df %>%  
                                                            dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
                                                            tidyr::drop_na(filtered_estimate) %>%
                                                            dplyr::summarise(num_results = length(filtered_estimate), .by = "simpleName") %>%
                                                            dplyr::arrange(desc(num_results)) %>%
                                                            droplevels %>%
                                                            dplyr::select(simpleName) %>%
                                                            tibble::deframe(.)))) %>%
  split(., f = .$simpleName) %>%
  map(., ~.x %>%
        dplyr::distinct(asv_id, .keep_all = FALSE) %>%
        unlist %>%
        as.character)


spearman.sig_corr_asvs_idx_list <- NULL
spearman.sig_strong_corr_asvs_idx_list <- NULL
for(i in seq(2, 11, 1)){
  namevar <- paste0("sig_corr_asvs_", i, "_idx")
  temp_idx <- temp_df %>%  
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
    tidyr::drop_na(estimate) %>%
    dplyr::summarise(num_results = length(estimate), .by = "asv_id") %>%
    dplyr::filter(num_results > i) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id) %>%
    droplevels %>%
    tibble::deframe(.)
  temp_idx <- list(temp_idx) %>% 
    setNames(., namevar)
  spearman.sig_corr_asvs_idx_list <- append(spearman.sig_corr_asvs_idx_list, temp_idx)
  
  namevar <- paste0("sig_strong_corr_asvs_", i, "_idx")
  temp_idx <- temp_df %>%  
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
    tidyr::drop_na(filtered_estimate) %>%
    dplyr::summarise(num_results = length(filtered_estimate), .by = "asv_id") %>%
    dplyr::filter(num_results >= i) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id) %>%
    droplevels %>%
    tibble::deframe(.)
  temp_idx <- list(temp_idx) %>% 
    setNames(., namevar)
  spearman.sig_strong_corr_asvs_idx_list <- append(spearman.sig_strong_corr_asvs_idx_list, temp_idx)
  rm(temp_idx)
  rm(namevar)
}
# View(spearman.sig_strong_corr_asvs_idx_list)
# View(spearman.sig_corr_asvs_idx_list)





# spearman.test.filtered.df <- temp_df %>%
#   dplyr::filter(asv_id %in% sig_corr_asvs_idx) %>%
#   # tidyr::drop_na(filtered_estimate) %>%
#   # dplyr::select(-filtered_estimate) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
#   # dplyr::arrange(asv_id, simpleName) %>%
#   droplevels
spearman.test.strong.filtered.list <- spearman.sig_strong_corr_asvs_idx_list %>%
  map(., ~dplyr::filter(temp_df, asv_id %in% .x) %>%
        tidyr::drop_na(filtered_estimate) %>%
        dplyr::select(-filtered_estimate) %>%
        droplevels %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
        dplyr::arrange(asv_id) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
        dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
        # dplyr::arrange(asv_id, simpleName) %>%
        droplevels)


# Summary of Spearman rank correlation coefficient tests ------------------

#if we filter for only the strongest correlations (abs(estimate) >= 0.5) between ASVs and metabolites
#we resolved 82 ASVs that each had 11+ correlations to 28 metabolites
#10+ correlations: 110 ASVs correlated to 28 metabolites
#9+: 123 ASVs and 28 metabolites
#8+: 145 ASVs and 28 metabolites
#7+: 167 ASVs and 28 metabolites


spearman_summary.df <- spearman.test.strong.filtered.list %>%
  map(., ~.x %>%
        dplyr::distinct(asv_id, simpleName) %>%
        dplyr::reframe(num_asvs = length(unique(.[["asv_id"]])),
                       num_metabs = length(unique(.[["simpleName"]])))) %>%
  bind_rows(., .id = "num_corrs") %>%
  dplyr::mutate(num_corrs = gsub("([[:alpha:]_]+)([[:digit:]]{1,2})(_.*$)", "\\2", num_corrs)) %>%
  dplyr::mutate(num_metabs = factor(num_metabs)) %>%
  dplyr::mutate(num_corrs = as.numeric(num_corrs)) 

g6 <- (ggplot(data = spearman_summary.df)
       + theme_bw()
       + geom_point(aes(x = num_corrs, y = num_asvs, fill = num_metabs), color = "black", shape = 21, size = 3, alpha = 1.0)
       + scale_x_continuous(n.breaks = 10, name = "Minimum number of significant correlations for each ASV")
       + scale_y_continuous(name = "Total ASVs with significant correlations above threshold", 
                            expand = expansion(mult = c(0.1,0.1)))
       + scale_fill_viridis(option = "mako", discrete = TRUE, name = "Number of \nobserved metabolites \n in correlations")
       + ggtitle("Distribution of strong correlations for ASVs and metabolites")
       + theme(panel.spacing = unit(1, "lines"),
               panel.background = element_blank(),
               axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
               axis.text.y = element_text(vjust = 0.5, hjust = 1),
               panel.grid.major = element_blank(),
               # panel.grid.minor.y = element_blank(),
               # panel.grid.minor.x = element_blank(),
               panel.ontop = FALSE,
               strip.text.y = element_blank(),
               legend.position = "bottom")
)

if(!any(grepl("spearman_strong_corr_summar", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_spearman_strong_corr_summary-", Sys.Date(), ".png"),
         g6,
         width = 6, height = 6, units = "in")
}

# Plot the heatmaps of correlations ---------------------------------------

as.hclust(dend_metab) %>%
  cutree(., k = 2) %>%
  tibble::enframe(name = "metabolite", value = "grouping")

as.hclust(dend_asv) %>%
  cutree(., k = 2) %>%
  tibble::enframe(name = "asv_id", value = "grouping")

#try to plot the number of significant correlations across each row (ASV) and column (metabolite):
# temp_spearman.df <- spearman.test.strong.filtered.list[[9]] %>%
#   dplyr::select(-label)
# temp_spearman.df2 <- bind_rows(
#                                (temp_spearman.df %>%
#                                   dplyr::group_by(asv_id) %>%
#                                   dplyr::summarise(label = length(estimate)) %>%
#                                   dplyr::mutate(simpleName = "total correlations"))) %>%
#   dplyr::bind_rows(., temp_spearman.df %>%
#                      dplyr::summarise(label = length(estimate), .by = "simpleName") %>%
#                      dplyr::mutate(grouping = 1) %>%
#                      dplyr::mutate(asv_id = "total correlations")) %>%
#   # dplyr::bind_rows(., temp_spearman.df %>%
#   #                    dplyr::summarise(label = length(estimate), .by = "simpleName") %>%
#   #                    dplyr::mutate(grouping = 2) %>%
#   #                    dplyr::mutate(asv_id = "total correlations")) %>%
#   # dplyr::bind_rows(., data.frame(asv_id = "",
#   #                                label = NA,
#   #                                simpleName = "")) %>%
#   bind_rows(temp_spearman.df, .) %>%
#   dplyr::left_join(., as.hclust(dend_asv) %>%
#                      cutree(., k = 2) %>%
#                      tibble::enframe(name = "asv_id", value = "grouping"),
#                    by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
#   dplyr::mutate(grouping = across(starts_with("grouping")) %>% purrr::reduce(coalesce)) %>%
#   dplyr::select(-ends_with(c(".x", ".y"))) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = c(labels(dend_asv), "total correlations"))) %>%
#   # dplyr::mutate(asv_id = factor(asv_id, levels = c(labels(dend_asv), "", "total"))) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
#   # dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "", "total"))) %>%
#   # dplyr::group_by(simpleName) %>% tidyr::complete(grouping) %>%
#   droplevels
# # spearman.test.filtered.df <- spearman.test.strong.filtered.list[[1]]


# temp_spearman.df1 <- temp_spearman.df %>%
#   dplyr::left_join(., as.hclust(dend_asv) %>%
#                      cutree(., k = 2) %>%
#                      tibble::enframe(name = "asv_id", value = "grouping"),
#                    by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
#   split(., f = .$grouping) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         dplyr::group_by(asv_id, grouping) %>%
#         dplyr::summarise(label = length(estimate), .groups = "keep") %>%
#         dplyr::mutate(simpleName2 = "total correlations")) %>%
#   map(., ~.x %>%
#         dplyr::left_join(temp_spearman.df %>%
#                            dplyr::select(asv_id, simpleName, estimate), by = join_by(asv_id))) %>%
#   bind_rows(., .id = NULL) %>%
#   droplevels
# 
# temp_spearman.df2 <- temp_spearman.df1 %>%
#   split(., f = .$simpleName) %>%
#   map(., ~.x %>%
#         dplyr::group_by(grouping, simpleName) %>%
#         dplyr::summarise(label = length(estimate), .groups = "keep") %>%
#         dplyr::mutate(asv_id = "total correlations") %>%
#         # droplevels)
#         droplevels) %>%
#   bind_rows(., .id = NULL) %>%
#   bind_rows(., temp_spearman.df1 %>%
#               dplyr::select(asv_id, grouping, label, simpleName2) %>%
#               dplyr::rename(simpleName = "simpleName2") %>%
#               dplyr::distinct(asv_id, grouping, label, simpleName) %>%
#               droplevels) %>%
#   bind_rows(., (temp_spearman.df %>%
#                   dplyr::left_join(., as.hclust(dend_asv) %>%
#                                      cutree(., k = 2) %>%
#                                      tibble::enframe(name = "asv_id", value = "grouping"),
#                                    by = join_by(asv_id), relationship = "many-to-many", multiple = "all") )) %>%
#   dplyr::select(asv_id, grouping, simpleName, label, estimate, sig) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = c(labels(dend_asv), "total correlations"))) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
#   droplevels
# 
# print(
#   ggplot(data = temp_spearman.df2 %>%
#            droplevels, aes(x = asv_id, y = simpleName))
#   + theme_bw()
#   # + geom_raster(aes(fill = estimate), hjust = 0, alpha = 0.7, show.legend = TRUE)
#   + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
#   + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
#   +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100),
#                           transform = "reverse",
#                           aesthetics = "fill",
#                           limits = c(1, -1),
#                           # expand = expansion(1.1,1.1),
#                           # name = "Similarity (%)",
#                           na.value = "white")
#   + scale_color_manual(values = c("grey", "black"),
#                        breaks = c("not", "sig"),
#                        labels = c("not", "sig"))
#   + scale_x_discrete(labels = usvi_genera_relabel,
#                      expand = c(0,0),
#                      # name = "Sample identifier")
#                      name = "Taxon")
#   + scale_y_discrete(name = "Metabolite",
#                      expand = c(0,0))
#   + theme(panel.spacing = unit(1, "lines"),
#           # + theme(
#           panel.background = element_blank(),
#           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
#           axis.text.y = element_text(vjust = 0.5, hjust = 1),
#           # axis.ticks.x = element_blank(),
#           # axis.ticks.y = element_blank(),
#           # axis.minor.ticks.x.bottom = element_line(linewidth = 2),
#           panel.grid.major = element_blank(),
#           # panel.grid.minor = element_line(color = "black", linewidth = 2),
#           # panel.grid.major.y = element_line(color = "black"),
#           # panel.grid.major.x = element_line(color = "black"),
#           panel.grid.minor.y = element_blank(),
#           panel.grid.minor.x = element_blank(),
#           panel.ontop = FALSE,
#           strip.text.y = element_blank())
#   # + facet_wrap(grouping~., scales = "free_y", drop = TRUE)
#   + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
#                                override.aes = list(stroke = 1, color = "black")),
#            color = "none")
#   + coord_flip()
# )
try(f_run_chunk())
if(execchunk) {
for(i in seq_len(length(spearman.test.strong.filtered.list))){
# for(i in seq_len(2)){
  namevar <- names(spearman.test.strong.filtered.list)[i] %>%
    gsub("([[:alpha:]_]+)(\\d{1,2})(_idx)$", "\\2", .)
  title_plot <- paste0("ASVs with ", namevar, " or more significant correlations")
  temp_spearman.df <- spearman.test.strong.filtered.list[[i]]
  fig_width <- (round(nrow(temp_spearman.df)/100, digits = 0) + 2)
  
  temp_spearman.df2 <- temp_spearman.df %>%
    dplyr::summarise(label = length(estimate), .by = "asv_id") %>%
    dplyr::mutate(simpleName = "total correlations") %>%
    dplyr::bind_rows(., temp_spearman.df %>%
                       dplyr::summarise(label = length(estimate), .by = "simpleName") %>%
                       dplyr::mutate(asv_id = "total correlations")) %>%
    bind_rows(temp_spearman.df %>% dplyr::select(-label), .) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = c(labels(dend_asv), "total correlations"))) %>%
    dplyr::arrange(asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
    dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
    droplevels
  
  temp_g6 <- (
    ggplot(data = temp_spearman.df2 %>%
             droplevels, aes(x = asv_id, y = simpleName))
    + theme_bw() 
    + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
    + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
    +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100), 
                            transform = "reverse", aesthetics = "fill", 
                            limits = c(1, -1), na.value = "white")
    + scale_color_manual(values = c("grey", "black"), 
                         breaks = c("not", "sig"), 
                         labels = c("not", "sig"))
    + scale_x_discrete(labels = usvi_genera_relabel,
                       expand = c(0,0), name = "Taxon")
    + scale_y_discrete(name = "Metabolite", expand = c(0,0))
    + theme(panel.spacing = unit(1, "lines"),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
            axis.text.y = element_text(vjust = 0.5, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.ontop = FALSE,
            strip.text.y = element_blank())
    + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
                                 override.aes = list(stroke = 1, color = "black")),
             color = "none")
    # + coord_flip()
    + ggtitle(title_plot)
  )
  
  assign(paste0("g6_asvs_", namevar, "_results"), temp_g6, envir = .GlobalEnv)
  if(!any(grepl("spearman_.*_strong_corr", list.files(projectpath, pattern = "usvi_.*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_spearman_", namevar, "_strong_corr-", Sys.Date(), ".png"),
           temp_g6,
           width = fig_width, height = 12, units = "in")
  }
  rm(temp_g6)
  rm(title_plot)
  rm(namevar)
}

}



# 
# #filter out only rho abs(estimates) >= 0.5
# spearman.test.filtered.df <- spearman.test.corrected %>%
#   tibble::as_tibble(., rownames = "asv_id") %>%
#   tidyr::pivot_longer(., cols = !asv_id,
#                       names_to = "simpleName",
#                       values_to = "padj") %>%
#   dplyr::right_join(., (spearman.test.rho %>%
#                           tibble::as_tibble(., rownames = "asv_id") %>%
#                           tidyr::pivot_longer(., cols = !asv_id,
#                                               names_to = "simpleName",
#                                               values_to = "estimate")),
#                     by = join_by(asv_id, simpleName)) %>%
#   dplyr::mutate(rho_threshold = dplyr::case_when((abs(estimate) >= 0.5) ~ 1,
#                                                  .default = NA)) %>%
#   dplyr::left_join(., usvi_sda_genera_compare.df %>%
#                      dplyr::distinct(asv_id, test_type, .keep_all = FALSE) %>%
#                      droplevels,
#                    by = join_by(asv_id), multiple = "all", relationship = "many-to-many") %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
#   dplyr::arrange(asv_id) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
#   # dplyr::mutate(simpleName = factor(simpleName, levels = unique(.[["simpleName"]]))) %>%
#   dplyr::arrange(asv_id, simpleName) %>%
#   tidyr::drop_na(padj) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
#   dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "label")) %>%
#   droplevels


print(
  ggplot(data = spearman.test.filtered.df %>%
           tidyr::drop_na(filtered_estimate) %>%
           droplevels, aes(x = asv_id, y = simpleName))
  + theme_bw() 
  # + geom_raster(aes(fill = estimate), hjust = 0, alpha = 0.7, show.legend = TRUE)
  + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
  +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100), 
                          transform = "reverse",
                          aesthetics = "fill", 
                          limits = c(1, -1),
                          # expand = expansion(1.1,1.1), 
                          # name = "Similarity (%)",
                          na.value = "white")
  + scale_color_manual(values = c("grey", "black"), 
                       breaks = c("not", "sig"), 
                       labels = c("not", "sig"))
  + scale_x_discrete(labels = usvi_genera_relabel,
                     expand = c(0,0),
                     # name = "Sample identifier")
                     name = "Taxon")
  + scale_y_discrete(name = "Metabolite",
                     expand = c(0,0))
  + theme(panel.spacing = unit(1, "lines"),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
          axis.text.y = element_text(vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.ontop = FALSE,
          strip.text.y = element_blank())
  + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
                               override.aes = list(stroke = 1, color = "black")),
           color = "none")
  + coord_flip()
)


# Plot relative abundances of significant taxa ----------------------------


#from rademu:
# rademu_reef_res.df %>%
#   dplyr::distinct(asv_id, group) %>%
#   droplevels

#these genera flagged as significantly differentially abundant between reef sites and seagrass seawater, contributed between 14.0% - 45.6% of total sequences in each sample


#from spearman:
# spearman.test.filtered.df %>%
#   dplyr::select(asv_id) %>%
#   dplyr::distinct(.)



temp_df <- usvi_sda_genera_compare.df %>%
              dplyr::ungroup(.) %>%
                  # dplyr::distinct(asv_id, test_type, group, site, .keep_all = FALSE) %>%
  dplyr::mutate(group = dplyr::coalesce(group, sampling_time)) %>% dplyr::select(-sampling_time) %>% dplyr::distinct(asv_id, group, test_type) %>%
  # dplyr::distinct(asv_id, test_type, site, .keep_all = FALSE) %>%
  # dplyr::rename(group = "site") %>%
                  dplyr::rename(category = "test_type") %>%
  dplyr::distinct(asv_id, category, .keep_all = TRUE) %>%
                  # tidyr::drop_na(group) %>%
  dplyr::arrange(asv_id) %>%
  droplevels %>%
    dplyr::full_join(.,
              (spearman.test.filtered.df %>%
                 dplyr::select(asv_id) %>%
                 dplyr::distinct(.) %>%
                 dplyr::mutate(spearman = "rho"))) %>%
  dplyr::mutate(group = dplyr::case_when(is.na(group) ~ "not",
                                         .default = group)) %>%
  # tidyr::drop_na(.) %>%
  # tidyr::pivot_wider(., id_cols = c("asv_id", "spearman"),
  #                    names_from = "group",
  #                    values_from = "category") %>%
  # dplyr::select(-not) %>%
  droplevels
# temp_df <- (usvi_sda_genera_compare.df %>%
#               dplyr::ungroup(.) %>%
#               # dplyr::distinct(asv_id, sampling_time, test_type) %>%
#               # dplyr::rename(category = "sampling_time") %>%
#               # dplyr::rename(group = "test_type") %>%
#               # dplyr::mutate(category = "rademu") %>%
#               # temp_df <- (rademu_reef_res.df %>%
#               tidyr::drop_na(group) %>%
#               dplyr::distinct(asv_id, group) %>%
#               # dplyr::distinct(asv_id, Combo) %>%
#               # dplyr::rename(group = "Combo") %>%
#               dplyr::mutate(category = "rademu") %>%
#               droplevels) %>%
#   bind_rows(., (usvi_sda_genera_compare.df %>%
#                   dplyr::ungroup(.) %>%
#                   # dplyr::distinct(asv_id, sampling_time, test_type) %>%
#                   # dplyr::rename(category = "sampling_time") %>%
#                   # dplyr::rename(group = "test_type") %>%
#                   # dplyr::mutate(category = "rademu") %>%
#                   # temp_df <- (rademu_reef_res.df %>%
#                   tidyr::drop_na(baseMean) %>%
#                   dplyr::distinct(asv_id, test_type, Combo, .keep_all = FALSE) %>%
#                   dplyr::rename(group = "Combo") %>%
#                   dplyr::rename(category = "test_type") %>%
#                   tidyr::drop_na(group) %>%
#                   # dplyr::mutate(category = "deseq") %>%
#                   droplevels) ) %>%
#   dplyr::arrange(asv_id) %>%
#   # tidyr::pivot_wider(., id_cols = "asv_id",
#   #                    names_from = "category",
#   #                    values_from = "group") %>%
#   # bind_rows(., 
#   #           (spearman.test.filtered.df %>%
#   #              dplyr::select(asv_id) %>%
#   #              dplyr::distinct(.) %>%
#   #              # dplyr::mutate(category = "spearman") %>%
#   #              dplyr::mutate(group = "rho"))) %>%
#   droplevels

usvi_key_genera.df <- usvi_sw_genus.tbl %>%
  dplyr::select(c(asv_id, all_of(unique(metabolomics_sample_metadata[["sample_id"]])))) %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  dplyr::right_join(., usvi_sda_genera_compare.df %>%
                      dplyr::ungroup(.) %>%
                      dplyr::distinct(asv_id) %>%
                    dplyr::left_join(., selfish_sda_genera_idx %>%
                                       tibble::enframe(name = "test_type.x", value = "asv_id")) %>%
                      dplyr::left_join(., data.frame(asv_id = shared_sda_genera_idx,
                                                     test_type.y = "both")) %>%
                      tidyr::unite("category", c("test_type.x", "test_type.y"), remove = TRUE, na.rm = TRUE) %>%
                    # usvi_sda_genera_compare.df %>%
                    #   dplyr::ungroup(.) %>%
                    #   dplyr::mutate(group = dplyr::coalesce(group, sampling_time)) %>% dplyr::select(-sampling_time) %>% dplyr::distinct(asv_id, group, test_type) %>%
                    #   dplyr::rename(category = "test_type") %>%
                    #   dplyr::distinct(asv_id, category, .keep_all = FALSE) %>%
                    #   dplyr::arrange(asv_id) %>%
                    #   droplevels %>%
                      dplyr::full_join(.,
                                       (spearman.test.filtered.df %>%
                                          dplyr::select(asv_id) %>%
                                          dplyr::distinct(asv_id) %>%
                                          dplyr::mutate(spearman = "rho"))) %>%
                      droplevels,
                    by = join_by(asv_id)) %>%
  droplevels %>%
  tidyr::pivot_longer(., cols = starts_with("Metab_"),
                      names_to = "sample_id",
                      values_to = "relabund") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::distinct(sample_id, site,  sampling_day, sampling_time) %>%
                         droplevels),
                   by = join_by(sample_id)) %>%
  dplyr::distinct(asv_id, sample_id, .keep_all = TRUE) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  dplyr::mutate(across(c(asv_id, category, spearman, sample_id, site, sampling_day, sampling_time), ~factor(.x))) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
  droplevels


# usvi_key_genera.df %>%
#   dplyr::group_by(sample_id) %>%
#   dplyr::summarise(totalabund = sum(relabund)) %>%
#   dplyr::arrange((totalabund))
#   # dplyr::arrange(desc(totalabund))


#between 46.5% and 85.3% of sequences are reprsented in these taxa identified through either radEmu or spearmanr ank correlation, as significantly differentially abundant or r2>0.5 correlated with metabolite concentrations

temp_df <- usvi_sw_genus.taxa.df %>%
  dplyr::filter(asv_id %in% unique(usvi_key_genera.df[["asv_id"]])) %>%
  droplevels %>%
  dplyr::select(asv_id, Genus, Domain, Phylum, Class, Order, Family) %>%
  dplyr::distinct(asv_id, Genus, Domain, Phylum, Class, Order, Family, .keep_all = FALSE) %>%
  dplyr::mutate(across(c(Class, Order, Genus), ~ dplyr::case_when(is.na(.x) ~ "NA",
                                                                  .default = .x)))

usvi_key_taxonomy_colors.df <- dplyr::right_join(annotation_taxa_colors_list[["Genus"]] %>%
                                                   tibble::enframe(., name = "Genus", value = "color"), 
                                                 usvi_sw_genus.taxa.df %>%
                                                   dplyr::filter(asv_id %in% unique(usvi_key_genera.df[["asv_id"]])) %>%
                                                   droplevels %>%
                                                   dplyr::select(asv_id, Genus, Domain, Phylum, Class, Order, Family) %>%
                                                   dplyr::distinct(asv_id, Genus, Domain, Phylum, Class, Order, Family, .keep_all = FALSE) %>%
                                                   dplyr::mutate(across(c(Class, Order, Genus), ~ dplyr::case_when(is.na(.x) ~ "NA",
                                                                                                                   .default = .x))),
                                                 # dplyr::mutate(Class = dplyr::case_when(is.na(Class) ~ "NA",
                                                 #                                        .default = Class)),
                                                 by = join_by(Genus)) %>%
  droplevels %>%
  tidyr::unite("taxonomy_class", c(Domain, Phylum, Class), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_order", c(Domain, Phylum, Class, Order), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_family", c(Domain, Phylum, Class, Order, Family), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_genus", c(Domain, Phylum, Class, Order, Family, Genus), sep = ";", remove = TRUE, na.rm = TRUE) %>%
  # dplyr::distinct(color, .keep_all = TRUE) %>%
  dplyr::select(asv_id, color, contains("taxonomy")) %>%
  tidyr::pivot_longer(., cols = !c(color, asv_id),
                      # dplyr::select(color, contains("taxonomy")) %>%
                      # tidyr::pivot_longer(., cols = !c(color),
                      names_to = "taxonomic_level",
                      values_to = "taxonomy") %>%
  dplyr::mutate(taxonomy = dplyr::case_when(grepl("NA;NA", taxonomy) ~ "NA",
                                            .default = taxonomy)) %>%
  droplevels %>%
  dplyr::group_by(asv_id, color) %>%
  dplyr::mutate(color = dplyr::case_when(is.na(color) ~ sample(color_resample_list, size = 1, replace = FALSE),
                                         .default = color)) %>%
  dplyr::ungroup(.) %>%
  tidyr::drop_na(.) %>%
  dplyr::arrange(taxonomy) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  droplevels


asv_colors <- usvi_key_taxonomy_colors.df %>%
  dplyr::select(asv_id, color) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  tibble::deframe(.)
taxonomy_colors <- usvi_key_taxonomy_colors.df %>%
  dplyr::select(taxonomy, color) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  tibble::deframe(.)

taxonomy_colors_lookup <- data.frame(taxonomy = names(taxonomy_colors),
                                     first = stringr::str_split_i(names(taxonomy_colors), ";", 2),
                                     second = stringr::str_split_i(names(taxonomy_colors), ";", -1)) %>%
  tidyr::unite(label, c("first", "second"), sep = ";") %>%
  dplyr::select(label, taxonomy) %>%
  tibble::deframe(.)
asv_colors_lookup <- usvi_key_taxonomy_colors.df %>%
  dplyr::select(asv_id, taxonomy) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::mutate(first = stringr::str_split_i(taxonomy, ";", 2),
                second = stringr::str_split_i(taxonomy, ";", -1)) %>%
  tidyr::unite(label, c("first", "second"), sep = ";") %>%
  dplyr::select(label, asv_id) %>%
  tibble::deframe(.)

temp_df <- usvi_key_genera.df %>%
  dplyr::left_join(., usvi_key_taxonomy_colors.df %>%
                     dplyr::filter(grepl("genus", taxonomic_level)) %>%
                     dplyr::distinct(asv_id, taxonomy, .keep_all = FALSE) %>%
                     droplevels,
                   by = join_by(asv_id),
                   relationship = "many-to-many", multiple = "all") %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, taxonomy, sample_id, .keep_all = TRUE) %>%
  droplevels

temp_df %>%
  dplyr::distinct(asv_id, category, spearman) %>%
  # tidyr::drop_na(rademu) %>%
  # tidyr::drop_na(spearman) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::group_by(category, spearman) %>%
  dplyr::summarise(num_taxa = length(asv_id))
temp_df %>%
  dplyr::distinct(asv_id, category, spearman) %>%
  # tidyr::drop_na(rademu) %>%
  # tidyr::drop_na(spearman) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(num_taxa = length(asv_id))



g3 <- print(ggplot(data = temp_df %>%
                     dplyr::filter(!is.na(category) & !is.na(spearman)) %>%
                     droplevels, aes(y = relabund,
                                     x = sample_id,
                                     group = taxonomy,
                                     # color = taxonomy,
                                     # fill = taxonomy))
                                     fill = asv_id))
            + geom_bar(width = 0.90, show.legend = TRUE, position = "stack", stat = "identity")
            + geom_bar(color = "black", width = 0.90, show.legend = FALSE, position = "stack", stat = "identity")
            + theme_bw()
            + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
            + scale_x_discrete(labels = gsub("((Yawzi_)|(Tektite_)|(LB_seagrass_))", "", sample_relabel), name = "Sample identifier")
            + scale_fill_manual(values = asv_colors, 
                                breaks = names(asv_colors),
                                labels = names(asv_colors_lookup),
                                # + scale_fill_manual(values = taxonomy_colors, 
                                #                     breaks = names(taxonomy_colors),
                                #                     labels = names(taxonomy_colors_lookup),
                                drop = FALSE)
            + facet_grid(scales = "free_x", space = "free",
                         labeller = labeller(site = site_lookup),
                         category ~ site)
            + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    strip.text.x = element_text(angle = 0),
                    strip.text.y = element_text(angle = 0))
            + labs(x = "sample",
                   y = "Relative abundance (%)")
            # + coord_cartesian(ylim = c(0,60), expand = FALSE)
            + guides(fill = guide_legend(order = 2, ncol = 2, title = "Taxonomy", direction = "vertical",
                                         override.aes = list(stroke = 1, color = "black")),
                     color = "none")
)
g4 <- print(ggplot(data = temp_df %>%
                     dplyr::filter(is.na(category) & !is.na(spearman)) %>%
                     droplevels, aes(y = relabund,
                                     x = sample_id,
                                     group = taxonomy,
                                     # color = taxonomy,
                                     fill = taxonomy))
            + geom_bar(width = 0.90, show.legend = TRUE, position = "stack", stat = "identity")
            + geom_bar(color = "black", width = 0.90, show.legend = FALSE, position = "stack", stat = "identity")
            + theme_bw()
            + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
            + scale_x_discrete(labels = gsub("((Yawzi_)|(Tektite_)|(LB_seagrass_))", "", sample_relabel), name = "Sample identifier")
            + scale_fill_manual(values = taxonomy_colors, 
                                breaks = names(taxonomy_colors),
                                labels = names(taxonomy_colors_lookup),
                                drop = TRUE)
            + facet_grid(scales = "free_x", space = "free",
                         labeller = labeller(site = site_lookup),
                         . ~ site)
            + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    strip.text.x = element_text(angle = 0),
                    strip.text.y = element_text(angle = 0))
            + labs(x = "sample",
                   y = "Relative abundance (%)")
            # + coord_cartesian(ylim = c(0,100), expand = FALSE)
            + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
                                         override.aes = list(stroke = 1, color = "black")),
                     color = "none")
)


gpatch <- (g3 + ggtitle("SDA taxa, strong correlations with metabolites") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) / (g4 + ggtitle("Strong correlations with metabolites, not SDA taxa")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Relative abundance of key genera", tag_levels = "A")
gpatch

if(!any(grepl("key_genera_relabund", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_key_genera_relabund-", Sys.Date(), ".png"),
         gpatch,
         width = 20, height = 16, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_key_genera_relabund_sda-", Sys.Date(), ".png"),
         g3,
         width = 20, height = 10, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_key_genera_relabund_not_sda-", Sys.Date(), ".png"),
         g4,
         width = 16, height = 8, units = "in")
}



# make networkg raph ------------------------------------------------------


temp_graph <- spearman.test.filtered.df %>%
  dplyr::filter(!is.na(rho_threshold) | !is.na(group)) %>%
  dplyr::select(asv_id, simpleName, estimate) %>%
  droplevels %>%
  droplevels %>%
  ggsankey::make_long(., asv_id, simpleName, value = "estimate") %>%
  tidyr::drop_na(next_node) %>%
  dplyr::rename(from = "node",
                to = "next_node",
                estimate = "value") %>%
  dplyr::select(from, to, estimate) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  igraph::graph_from_data_frame(., directed = TRUE)

plot(temp_graph)
# igraph::degree_distribution(temp_graph)
# igraph::degree(temp_graph)
# igraph::component_distribution(temp_graph)
# igraph::count_components(temp_graph)
# which.max(igraph::degree(temp_graph))
# 
# edge_density(temp_graph)

#make a tidy-compatible dataframe to plot the graph:
temp_graph_coords <- igraph::layout_nicely(temp_graph, dim = 2) %>%
  tibble::as_tibble(.) %>%
  dplyr::mutate(nodes = vertex_attr(temp_graph, "name")) %>%
  dplyr::relocate(nodes) %>%
  dplyr::mutate(type = dplyr::case_when(grepl("ASV_", nodes) ~ "taxon",
                                        .default = "metabolite")) %>%
  dplyr::left_join(., spearman.test.filtered.df %>%
                     dplyr::select(asv_id, group) %>%
                     dplyr::distinct(asv_id, .keep_all = TRUE) %>%
                     droplevels,
                   by = join_by("nodes" =="asv_id"))

# print(
#   ggplot(data = temp_graph_coords, aes(x = V1, y = V2, fill = type))
#   + theme_bw()
#   + geom_point(shape = 21)
# )

temp_df <- temp_graph_coords %>%
  dplyr::select(nodes, V1, V2) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::rename(xstart = "V1",
                node_start = "nodes",
                ystart = "V2") %>%
  droplevels

temp_graph_edges <- spearman.test.filtered.df %>%
  dplyr::filter(!is.na(rho_threshold) | !is.na(group)) %>%
  dplyr::select(asv_id, simpleName, estimate) %>%
  droplevels %>%
  ggsankey::make_long(., asv_id, simpleName, value = "estimate") %>%
  tidyr::drop_na(next_node) %>%
  dplyr::rename(node_start = "node",
                node_end = "next_node",
                estimate = "value") %>%
  dplyr::select(node_start, node_end, estimate) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::left_join(., temp_df,
                   by = join_by(node_start)) %>%
  dplyr::left_join(., temp_df %>%
                     dplyr::rename(node_end = "node_start",
                                   xend = "xstart",
                                   yend = "ystart"),
                   by = join_by(node_end)) %>%
  # dplyr::left_join(., temp_graph_coords %>%
  #   dplyr::select(nodes, type, color, group) %>%
  #   droplevels,
  # by = join_by("node_start" == "nodes")) %>%
  dplyr::mutate(direction = dplyr::case_when((estimate < 0) ~ "negative",
                                             (estimate > 0) ~ "positive",
                                             .default = NA)) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  droplevels


#plot the taxa with higher abundance in reefs than in seagrass, and their metabolite correlations:

keep <- spearman.test.filtered.df %>%
  dplyr::filter(grepl("high", group)) %>%
  dplyr::filter(!is.na(rho_threshold)) %>%
  dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
  droplevels %>%
  tidyr::pivot_longer(., cols = everything(),
                      names_to = NULL,
                      values_to = "nodes") %>%
  dplyr::arrange(nodes) %>%
  distinct(.) %>%
  unlist(.)

temp_df2 <- temp_graph_coords %>%
  dplyr::filter(grepl(paste0(keep, collapse = "|"), nodes)) %>%
  dplyr::mutate(nodes = dplyr::case_when(grepl("ASV_", nodes) ~ dplyr::recode_factor(nodes, !!!usvi_genera_relabel),
                                         .default = nodes)) %>%
  dplyr::mutate(nodes = stringr::str_split_i(nodes, ";", -1)) %>%
  droplevels

temp_df1 <- temp_graph_edges %>%
  dplyr::right_join(., spearman.test.filtered.df %>%
                      dplyr::filter(grepl("high", group)) %>%
                      dplyr::filter(!is.na(rho_threshold)) %>%
                      dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
                      droplevels,
                    by = join_by("node_start" == "asv_id", "node_end" == "simpleName")) %>%
  dplyr::mutate(node_start = dplyr::case_when(grepl("ASV_", node_start) ~ dplyr::recode_factor(node_start, !!!usvi_genera_relabel),
                                              .default = node_start)) %>%
  dplyr::mutate(node_start = stringr::str_split_i(node_start, ";", -1)) %>%
  droplevels



#now do the lower-abudnace taxa and their metabolite correlations

keep <- spearman.test.filtered.df %>%
  dplyr::filter(grepl("low", group)) %>%
  dplyr::filter(!is.na(rho_threshold)) %>%
  dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
  droplevels %>%
  tidyr::pivot_longer(., cols = everything(),
                      names_to = NULL,
                      values_to = "nodes") %>%
  dplyr::arrange(nodes) %>%
  distinct(.) %>%
  unlist(.)

temp_df4 <- temp_graph_coords %>%
  dplyr::filter(grepl(paste0(keep, collapse = "|"), nodes)) %>%
  dplyr::mutate(nodes = dplyr::case_when(grepl("ASV_", nodes) ~ dplyr::recode_factor(nodes, !!!usvi_genera_relabel),
                                         .default = nodes)) %>%
  dplyr::mutate(nodes = stringr::str_split_i(nodes, ";", -1)) %>%
  droplevels

temp_df3 <- temp_graph_edges %>%
  dplyr::right_join(., spearman.test.filtered.df %>%
                      dplyr::filter(grepl("low", group)) %>%
                      dplyr::filter(!is.na(rho_threshold)) %>%
                      dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
                      droplevels,
                    by = join_by("node_start" == "asv_id", "node_end" == "simpleName")) %>%
  dplyr::mutate(node_start = dplyr::case_when(grepl("ASV_", node_start) ~ dplyr::recode_factor(node_start, !!!usvi_genera_relabel),
                                              .default = node_start)) %>%
  dplyr::mutate(node_start = stringr::str_split_i(node_start, ";", -1)) %>%
  droplevels

temp_list <- list(list(temp_df1, temp_df2) %>%
                    setNames(., c("connections", "nodes")), 
                  (list(temp_df3, temp_df4) %>%
                     setNames(., c("connections", "nodes")))) %>%
  setNames(., c("high_reef", "low_reef"))


g1 <- print(
  ggplot(data = temp_df1)
  + theme_bw()
  + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend,
                     linewidth = abs_estimate, color = direction),
                 arrow = arrow(length = unit(0.03, "npc"),
                               angle = 30, ends = "last", type = "closed"),
                 lineend = "butt", linejoin = "round",
                 alpha = 0.7, show.legend = TRUE)
  + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend),
                 lineend = "butt", linejoin = "round", linetype = 2, color = "black",
                 alpha = 0.7, show.legend = TRUE)
  + geom_point(data = temp_df2, aes(x = V1, y = V2, fill = type), shape = 21, size = 3, alpha = 1)
  + ggrepel::geom_label_repel(data = temp_df2, aes(x = V1, y = V2, label = nodes, fill = type),
                              direction = "both", segment.color = NA, seed = 123,
                              force = 1, force_pull = 1,
                              point.padding = unit(0.01, "npc"),
                              box.padding = unit(0.01, "npc"),
                              colour = "black", fontface = "bold")
  + scale_color_manual(values = c("#dec1aa", "#03c4da"),
                       breaks = c("negative", "positive"),
                       labels = c("negative", "positive"))
  + scale_fill_manual(values = c("#c5b8dc", "#b9d2b1"),
                      breaks = c("taxon", "metabolite"),
                      labels = c("taxon", "metabolite"))
  + scale_linewidth_continuous(name = "Strength of relationship")
  + scale_x_continuous(expand = expansion(0.1,0.1))
  + scale_y_continuous(expand = expansion(0.1,0.1))
  + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())
  + ggtitle("Taxa significantly more abundant in reef sites and strong correlations with metabolites")
)



for(i in 1:length(temp_list)){
  namevar <- eval(names(temp_list)[i]) %>%
    gsub("_", " in ", .)
  temp_df1 <- temp_list[[i]]$connections %>%
    as.data.frame(.)
  temp_df2 <- temp_list[[i]]$nodes %>%
    as.data.frame(.)
  
  g <- print(
    ggplot(data = temp_df1)
    # ggplot(data = temp_list[[i]]$connections)
    + theme_bw()
    + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend,
                       linewidth = abs_estimate, color = direction), 
                   arrow = arrow(length = unit(0.03, "npc"), 
                                 angle = 30, ends = "last", type = "closed"), 
                   lineend = "butt", linejoin = "round",
                   alpha = 0.7, show.legend = TRUE)
    + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend), 
                   lineend = "butt", linejoin = "round", linetype = 2, color = "black",
                   alpha = 0.7, show.legend = FALSE)
    + geom_point(data = temp_df2,
                 # + geom_point(data = temp_list[[i]]$node, 
                 aes(x = V1, y = V2, fill = type), shape = 21, size = 3, alpha = 1)
    + ggrepel::geom_label_repel(aes(x = V1, y = V2, label = nodes, fill = type),
                                # data = temp_list[[i]]$node, 
                                data = temp_df2, 
                                max.overlaps = 30,
                                direction = "both", segment.color = NA, seed = 123,
                                force = 1, force_pull = 1,
                                point.padding = unit(0.01, "npc"),
                                box.padding = unit(0.01, "npc"),
                                colour = "black", fontface = "bold")
    + scale_color_manual(values = c("#dec1aa", "#03c4da"),
                         breaks = c("negative", "positive"),
                         labels = c("negative", "positive"), drop = FALSE)
    + scale_fill_manual(values = c("#c5b8dc", "#b9d2b1"),
                        breaks = c("taxon", "metabolite"),
                        labels = c("taxon", "metabolite"), drop = FALSE)
    + scale_linewidth_continuous(name = "Strength of relationship", limits = c(0, 0.9))
    + scale_x_continuous(expand = expansion(0.1,0.1))
    + scale_y_continuous(expand = expansion(0.1,0.1))
    + theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
    + ggtitle(paste0("Taxa significantly ", namevar, " sites"))
  )
  
  assign(paste0("g", i), g, envir = .GlobalEnv, inherits = TRUE)
  rm(g)
  rm(namevar)
}

gpatch2  <- ((g1 + theme(legend.position = "none")) | g2) + patchwork::plot_layout(guides = "collect")
gpatch2 <- gpatch2 + patchwork::plot_annotation(title = "Network graphs of metabolite concentrations related to taxa abundances",
                                                tag_levels = "A")
gpatch2
if(!any(grepl("key_taxa_network", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_key_taxa_network-", Sys.Date(), ".png"),
         gpatch2,
         width = 18, height = 10, units = "in")
}

if(!any(grepl("spearman_df", list.files(projectpath, pattern = "usvi_.*.RData")))){
  save(spearman.test.corrected, spearman.test.filtered.df, file = paste0(projectpath, "/", "usvi_spearman_df-", Sys.Date(), ".RData"))
}


#pantothenic acid emerged as significant coral-produced metabolites (Weber et al., 2022)