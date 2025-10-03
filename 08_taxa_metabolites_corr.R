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


if(file.exists(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"))){
  metabolomics_sample_metadata <- readr::read_delim(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
} else {
  cli::cli_alert_warning("Please process the metabolomics sample metadata.")
}

potential_metab_outliers_idx <- c(73) %>%
  paste0("CINAR_BC_", .)
keep <- c("sample_id", "metab_deriv_label", "sample_type", "site", "sampling_time", "sampling_day", "sampling_date", "depth", "site_type",
          "fcm_Prochlorococcus", "fcm_Synechococcus", "fcm_Picoeukaryotes", "fcm_Unpigmented_cells",
          "replicate", "grouping", "PAR", "lumens", "lux", "temp")
usvi_selected_metadata <- metabolomics_sample_metadata %>%
  dplyr::filter(!(metab_deriv_label %in% potential_metab_outliers_idx)) %>%
  dplyr::filter(!(grepl("Day1", sampling_day))) %>% #drop day1 samples
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
  # tidyr::drop_na(.) %>% #removed sample Metab_306 due to lack of FCM measurements
  droplevels




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
  usvi_metab_mat <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(-LODflag) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id,  metab_deriv_label) %>%
                           droplevels),
                     by = c("metab_deriv_label")) %>%
    dplyr::mutate(across(c(sample_id, simpleName), ~factor(.x))) %>%
    dplyr::relocate(sample_id) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-metab_deriv_label) %>%
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
                       # values_fill = 0,
                       values_from = "conc") %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::vegdist(., binary = FALSE, upper = TRUE, na.rm = TRUE,
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
  
  temp_mat <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(-LODflag) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::relocate(sample_id) %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::mutate(conc = log2(conc + 1)) %>%
    tidyr::pivot_wider(., id_cols = "metab_deriv_label",
                       names_from = "simpleName",
                       # values_fill = 0,
                       values_from = "conc") %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
    as.matrix(.) %>%
    vegan::vegdist(., binary = FALSE, upper = TRUE, na.rm = TRUE,
                   # distance = "horn",
                   distance = "bray",
                   autotransform = TRUE) %>%
    as.matrix(.)
  temp_mat <- temp_mat[intersect(names(sample_relabel), rownames(temp_mat)), intersect(names(sample_relabel), colnames(temp_mat))]
  
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

#whatis the minimum number of samples in which a metabolite is observed?
#3: amMP
usvi_metabolomics_long.df %>%
  dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
  dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
  dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
  #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
  dplyr::filter(LODflag == 0) %>%
  dplyr::summarise(num_obs = length(conc), .by = c("simpleName")) %>%
  dplyr::arrange(num_obs)

# Plot heatmaps of the distance matrices ----------------------------------



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




{
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
  
}



# Spearman rank correlation coefficient calculations ----------------------


#option B:
if(file.exists(paste0(projectpath, "/", "spearman.test.optB.list", ".rds"))){
  spearman.test.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.optB.list", ".rds"))
} else { 
  list2env(usvi_dist_mat_list_optB, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, site, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::select(-LODflag) %>%
    dplyr::relocate(sample_id) %>%
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
                                    NA)) %>%
                                    # 0)) %>%
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
  spearman.test.n <- spearman.test
  
  y <- length(colnames(spearman.test))
  for(j in seq_len(y)){
    vector_metab <- usvi_metab.tbl[, j]
    for(i in seq_len(nrow(spearman.test))){
      vector_microb <- usvi_asv.tbl[,i]
      vector_microb <- vector_microb[!is.na(vector_metab)]
      vector_metab_na <- vector_metab[!is.na(vector_metab)]
      
      if(length(vector_metab_na) >= 3 & sum(vector_microb) > 0){ #if 3 of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
        spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
          purrr::pluck(., "p.value")
        spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
          purrr::pluck(., "estimate")
        spearman.test.n[i, j] <- length(vector_metab_na)
      } else {
        spearman.test[i, j] <- NA
        spearman.test.rho[i, j] <- NA
        spearman.test.n[i, j] <- length(vector_metab_na)
      }
    }
  }
  
  
  #for a q-value of 0.05, the adjusted p threshold should be 0.0206
  #and we have X significant relatinships
  #For a q-value of 0.1, the adjusted p-threshold should be 0.0579
  #and we have Y significant relationships
  
  padj_cutoff <- spearman.test %>%
    apply(., 2, na.omit) %>% unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
    quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
  
  spearman.test.bh.corrected <- spearman.test %>%
    apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
    apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  
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
    dplyr::mutate(label = signif(estimate, digits = 2)) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
    droplevels
  
  # spearman.test.df %>%
  #   tidyr::drop_na(padj_05) %>%
  #   dplyr::summarise(num_results = length(padj_bh), .by = "asv_id") %>%
  #   dplyr::arrange(desc(num_results)) %>%
  #   dplyr::select(num_results) %>%
  #   tibble::deframe(.) %>%
  #   quantile(., probs = seq(0, 1, 0.25), names = TRUE)
  
  #passing the q-value cutoff of 0.05:
  #of the 48 metabolites, the maximum number of significantly correlated ASVs was 652, the median was 142
  #of the 975 ASVs that had significant correlations with metabolite concentrations,
  #the maximum number of significantly correlated metabolites for an ASV was 23.  the median was 1.
  
  #using the BH adjusted p-values < 0.05:
  #of the 48 metabolites, the maximum number of significantly correlated ASVs was 658, the median was 154.
  #of the 1742 ASVs that had significant correlations with metabolite concentrations,
  #the maximum number of significantly correlated metabolites for an ASV was 27.  the median was 1.
  
  #using the q-value cutoff, we lose 78 ASVs that had an adjusted p-value < 0.05 but not below the 0.05 q-value cutoff.
  
  spearman.test.optB.list <- list(spearman.test.df, 
                                  spearman.test, spearman.test.rho, spearman.test.n,
                                  dend_asv, dend_metab, padj_cutoff) %>%
    setNames(., c("spearman.test.df", 
                  "spearman.test", "spearman.test.rho", "spearman.test.n",
                  "dend_asv", "dend_metab", "padj_cutoff"))
  
  # spearman.test.optB.list <- list(spearman.test.df, dend_asv, dend_metab, padj_cutoff) %>%
  #   setNames(., c("spearman.test.df", "dend_asv", "dend_metab", "padj_cutoff"))
  readr::write_rds(spearman.test.optB.list, paste0(projectpath, "/", "spearman.test.optB.list", ".rds"), compress = "gz")
}

#option A:
if(file.exists(paste0(projectpath, "/", "spearman.test.optA.list", ".rds"))){
  spearman.test.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.optA.list", ".rds"))
} else {
  list2env(usvi_dist_mat_list_optA, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
                           dplyr::select(sample_id, site, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::select(-LODflag) %>%
    dplyr::relocate(sample_id) %>%
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
                                    NA)) %>%
                                    # 0)) %>%
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
  spearman.test.n <- spearman.test  
  
  y <- length(colnames(spearman.test))
  for(j in seq_len(y)){
    vector_metab <- usvi_metab.tbl[, j]
    for(i in seq_len(nrow(spearman.test))){
      vector_microb <- usvi_asv.tbl[,i]
      vector_microb <- vector_microb[!is.na(vector_metab)]
      vector_metab_na <- vector_metab[!is.na(vector_metab)]
      
      if(length(vector_metab_na) >= 3 & sum(vector_microb) > 0){ #if 3 of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
        spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
          purrr::pluck(., "p.value")
        spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
          purrr::pluck(., "estimate")
        spearman.test.n[i, j] <- length(vector_metab_na)
      } else {
        spearman.test[i, j] <- NA
        spearman.test.rho[i, j] <- NA
        spearman.test.n[i, j] <- length(vector_metab_na)
      }
    }
  }
  
  padj_cutoff <- spearman.test %>%
    apply(., 2, na.omit) %>% unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
    quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
  
  spearman.test.bh.corrected <- spearman.test %>%
    apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
    apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
  
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
  
  spearman.test.optA.list <- list(spearman.test.df, 
                                  spearman.test, spearman.test.rho, spearman.test.n,
                                  dend_asv, dend_metab, padj_cutoff) %>%
    setNames(., c("spearman.test.df", 
                  "spearman.test", "spearman.test.rho", "spearman.test.n",
                  "dend_asv", "dend_metab", "padj_cutoff"))
  # spearman.test.optA.list <- list(spearman.test.df, dend_asv, dend_metab, padj_cutoff) %>%
  #   setNames(., c("spearman.test.df", "dend_asv", "dend_metab", "padj_cutoff"))
  readr::write_rds(spearman.test.optA.list, paste0(projectpath, "/", "spearman.test.optA.list", ".rds"), compress = "gz")
}


# Repeat Spearman on site-specific ----------------------------------------


#note: since some of the metabolites were measured/recorded above detection only a handful of times in each site
#the correlations can be junky
#throw out any cor.test result between a metabolite and ASV where there are fewer than 8 measurements of the metabolite in the site's samples (so <1/3 of samples had measurable metabolite)

#option A:
if(file.exists(paste0(projectpath, "/", "spearman.test.site.optA.list", ".rds"))){
  spearman.test.site.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.optA.list", ".rds"))
} else {
  list2env(usvi_dist_mat_list_optA, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(-LODflag) %>%
    dplyr::relocate(sample_id) %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(sample_id, simpleName) %>%
    dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE), 
                     num = length(conc),
                     .groups = "keep",
                     sd = sd(conc, na.rm = TRUE)) %>%
    dplyr::rename(conc = "mean_conc") %>%
    dplyr::mutate(log_conc = ifelse(!(is.na(conc) | (conc < 0)),
                                    log2(conc+1), #log transform abundance (with +1 pseudocount)
                                    # 0)) %>%
                                    NA)) %>%
    dplyr::select(sample_id, simpleName, log_conc) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "sample_id",
                       values_from = "log_conc",
                       names_from = "simpleName") %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    droplevels
  
  usvi_asv.tbl <- usvi_prok_asvs.df %>%
    dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
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
  
  usvi_metab_site.list <- metabolomics_sample_metadata %>%
      dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
      dplyr::distinct(sample_id, site) %>%
      droplevels %>%
    split(., f = .$site) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(sample_id) %>%
          dplyr::left_join(., usvi_metab.tbl %>%
                             tibble::as_tibble(rownames = "sample_id")) %>%
          # dplyr::select(which(colSums(.) > 0)) %>% #remove any that are insignificant.
          tibble::column_to_rownames(var = "sample_id"))
    
  usvi_asv_site.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
    dplyr::distinct(sample_id, site) %>%
    droplevels %>%
    split(., f = .$site) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(sample_id) %>%
          dplyr::inner_join(., usvi_asv.tbl %>%
                             tibble::as_tibble(rownames = "sample_id") %>%
                             droplevels) %>%
          tibble::column_to_rownames(var = "sample_id") %>%
          dplyr::select(which(colSums(.) > 0)) %>%
          droplevels)
  

  #for loop for each site:
  # for(site in names(usvi_asv_site.list)[3]){
  for(site in names(usvi_asv_site.list)){
    namevar <- site
    temp_asv.tbl <- usvi_asv_site.list[[site]]
    temp_metab.tbl <- usvi_metab_site.list[[site]]
    
    temp_spearman.test <- matrix(nrow = ncol(temp_asv.tbl), ncol = ncol(temp_metab.tbl))
    colnames(temp_spearman.test) <- colnames(temp_metab.tbl)
    rownames(temp_spearman.test) <- colnames(temp_asv.tbl)
    
    temp_spearman.test.rho <- temp_spearman.test
    temp_spearman.test.n <- temp_spearman.test
    
    y <- length(colnames(temp_spearman.test))
    
    # vector_metab <- temp_metab.tbl[, 38] #"spermidine 3" has only 4 non-NA entries in Yawzi. so the correlation is junky...
    # vector_microb <- temp_asv.tbl[, 1]
    #     vector_microb <- vector_microb[!is.na(vector_metab)]
    #     vector_metab_na <- vector_metab[!is.na(vector_metab)]
    #     cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
    #               purrr::pluck(., "p.value")
    for(j in seq_len(y)){
      vector_metab <- temp_metab.tbl[, j]
      for(i in seq_len(nrow(temp_spearman.test))){
        vector_microb <- temp_asv.tbl[,i]
        vector_microb <- vector_microb[!is.na(vector_metab)]
        vector_metab_na <- vector_metab[!is.na(vector_metab)]
        
        # if(length(vector_metab_na) > 0){ #if the metabolite has non-NA values for the samples in this site:
        # if(length(vector_metab_na) >= 8){ #if 1/3 of the samples have a non-NA value for that metabolite
        # if(length(vector_metab_na) >= 4){ #if 1/6 of the samples have a non-NA value for that metabolite
        if(length(vector_metab_na) >= 3 & sum(vector_microb) > 0){ #if 1/6 of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
          temp_spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
            purrr::pluck(., "p.value")
          temp_spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
            purrr::pluck(., "estimate")
          temp_spearman.test.n[i, j] <- length(vector_metab_na)
        } else {
          temp_spearman.test[i, j] <- NA
          temp_spearman.test.rho[i, j] <- NA
          temp_spearman.test.n[i, j] <- length(vector_metab_na)
        }
      }
    }
    
    ##for troubleshooting:
    # assign(paste0("temp_spearman.test.", namevar), temp_spearman.test, envir = .GlobalEnv)
    # assign(paste0("temp_spearman.test.rho.", namevar), temp_spearman.test.rho, envir = .GlobalEnv)
    
    padj_cutoff <- temp_spearman.test %>%
      apply(., 2, na.omit) %>%
      unlist %>%
      ashr::qval.from.lfdr(.) %>%
      as.matrix(.) %>%
      quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
    
    temp_spearman.test.bh.corrected <- temp_spearman.test %>%
      apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
      apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
    
    # temp_dend_asv <- temp_spearman.test.rho %>%
    #   dist(t(.), method = "euclidean") %>%
    #   hclust(method = "ward.D2") %>%
    #   as.dendrogram
    # temp_dend_metab <- temp_spearman.test.rho %>%
    #   t() %>%
    #   dist(t(.), method = "euclidean") %>%
    #   hclust(method = "ward.D2") %>%
    #   as.dendrogram
    
    temp_spearman.test.df <- temp_spearman.test.bh.corrected %>%
      tibble::as_tibble(., rownames = "asv_id") %>%
      tidyr::pivot_longer(., cols = !asv_id,
                          names_to = "simpleName",
                          values_to = "padj_bh") %>%
      droplevels %>%
      dplyr::left_join(., (temp_spearman.test %>%
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
      dplyr::right_join(., (temp_spearman.test.rho %>%
                              tibble::as_tibble(., rownames = "asv_id") %>%
                              tidyr::pivot_longer(., cols = !asv_id,
                                                  names_to = "simpleName",
                                                  values_to = "estimate")),
                        by = join_by(asv_id, simpleName)) %>%
      dplyr::mutate(sig = dplyr::case_when(!is.na(padj_bh) & !is.na(padj_05) ~ "sig",
                                           !is.na(padj_bh) & is.na(padj_05) ~ "maybe",
                                           .default = "not")) %>%
      dplyr::mutate(sig = factor(sig)) %>%
      # dplyr::mutate(asv_id = factor(asv_id, levels = labels(temp_dend_asv))) %>%
      # dplyr::mutate(simpleName = factor(simpleName, levels = labels(temp_dend_metab))) %>%
      dplyr::arrange(asv_id, simpleName) %>%
      dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      # tidyr::drop_na(contains("padj")) %>%
      dplyr::mutate(label = signif(estimate, digits = 2)) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
      droplevels
    
    spearman.test.site.list <- list(temp_spearman.test.df, 
                                    temp_spearman.test, temp_spearman.test.rho, temp_spearman.test.n,
                                    # temp_dend_asv, temp_dend_metab, 
                                    padj_cutoff) %>%
      setNames(., c("spearman.test.df", 
                    "spearman.test", "spearman.test.rho", "spearman.test.n",
                    # "dend_asv", "dend_metab", 
                    "padj_cutoff"))
    
    assign(paste0("spearman.test.siteA.", namevar, ".list"), spearman.test.site.list, envir = .GlobalEnv)
    rm(spearman.test.site.list)
    rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
  }
  
  
  spearman.test.site.optA.list <- lapply(apropos("^spearman.test.siteA.*$", mode = "list"),
                                         get) %>%
    setNames(., names(usvi_asv_site.list))
  readr::write_rds(spearman.test.site.optA.list, paste0(projectpath, "/", "spearman.test.site.optA.list", ".rds"), compress = "gz")

}
rm(list = apropos("^(spearman.test.siteA.)(.*)$", mode = "list"))


#option B:
if(file.exists(paste0(projectpath, "/", "spearman.test.site.optB.list", ".rds"))){
  spearman.test.site.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.optB.list", ".rds"))
} else {
  list2env(usvi_dist_mat_list_optB, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(-LODflag) %>%
    dplyr::relocate(sample_id) %>%
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
                                    # 0)) %>%
                                    NA)) %>%
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
  
  usvi_metab_site.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::distinct(metab_deriv_label, site) %>%
    droplevels %>%
    split(., f = .$site) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(metab_deriv_label) %>%
          dplyr::left_join(., usvi_metab.tbl %>%
                             tibble::as_tibble(rownames = "metab_deriv_label")) %>%
          tibble::column_to_rownames(var = "metab_deriv_label"))
  
  usvi_asv_site.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::distinct(metab_deriv_label, site) %>%
    droplevels %>%
    split(., f = .$site) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(metab_deriv_label) %>%
          dplyr::inner_join(., usvi_asv.tbl %>%
                              tibble::as_tibble(rownames = "metab_deriv_label") %>%
                              droplevels) %>%
          tibble::column_to_rownames(var = "metab_deriv_label") %>%
          dplyr::select(which(colSums(.) > 0)) %>%
          droplevels)
  
  
  #for loop for each site:
  for(site in names(usvi_asv_site.list)){
    namevar <- site
    temp_asv.tbl <- usvi_asv_site.list[[site]]
    temp_metab.tbl <- usvi_metab_site.list[[site]]
    
    temp_spearman.test <- matrix(nrow = ncol(temp_asv.tbl), ncol = ncol(temp_metab.tbl))
    colnames(temp_spearman.test) <- colnames(temp_metab.tbl)
    rownames(temp_spearman.test) <- colnames(temp_asv.tbl)
    
    temp_spearman.test.rho <- temp_spearman.test
    temp_spearman.test.n <- temp_spearman.test
    
    y <- length(colnames(temp_spearman.test))
    
    for(j in seq_len(y)){
      vector_metab <- temp_metab.tbl[, j]
      for(i in seq_len(nrow(temp_spearman.test))){
        vector_microb <- temp_asv.tbl[,i]
        vector_microb <- vector_microb[!is.na(vector_metab)]
        vector_metab_na <- vector_metab[!is.na(vector_metab)]
        
        if(length(vector_metab_na) >= 3 & sum(vector_microb) > 0){ #if 1/6 of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
          temp_spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
            purrr::pluck(., "p.value")
          temp_spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
            purrr::pluck(., "estimate")
          temp_spearman.test.n[i, j] <- length(vector_metab_na)
        } else {
          temp_spearman.test[i, j] <- NA
          temp_spearman.test.rho[i, j] <- NA
          temp_spearman.test.n[i, j] <- length(vector_metab_na)
        }
      }
    }
    
    padj_cutoff <- temp_spearman.test %>%
      apply(., 2, na.omit) %>%
      unlist %>%
      ashr::qval.from.lfdr(.) %>%
      as.matrix(.) %>%
      quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
    
    temp_spearman.test.bh.corrected <- temp_spearman.test %>%
      apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
      apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
    
    temp_spearman.test.df <- temp_spearman.test.bh.corrected %>%
      tibble::as_tibble(., rownames = "asv_id") %>%
      tidyr::pivot_longer(., cols = !asv_id,
                          names_to = "simpleName",
                          values_to = "padj_bh") %>%
      droplevels %>%
      dplyr::left_join(., (temp_spearman.test %>%
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
      dplyr::right_join(., (temp_spearman.test.rho %>%
                              tibble::as_tibble(., rownames = "asv_id") %>%
                              tidyr::pivot_longer(., cols = !asv_id,
                                                  names_to = "simpleName",
                                                  values_to = "estimate")),
                        by = join_by(asv_id, simpleName)) %>%
      dplyr::mutate(sig = dplyr::case_when(!is.na(padj_bh) & !is.na(padj_05) ~ "sig",
                                           !is.na(padj_bh) & is.na(padj_05) ~ "maybe",
                                           .default = "not")) %>%
      dplyr::mutate(sig = factor(sig)) %>%
      dplyr::arrange(asv_id, simpleName) %>%
      dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      dplyr::mutate(label = signif(estimate, digits = 2)) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
      droplevels
    
    spearman.test.site.list <- list(temp_spearman.test.df, 
                                    temp_spearman.test, temp_spearman.test.rho, temp_spearman.test.n,
                                    padj_cutoff) %>%
      setNames(., c("spearman.test.df", 
                    "spearman.test", "spearman.test.rho", "spearman.test.n",
                    "padj_cutoff"))
    
    assign(paste0("spearman.test.siteB.", namevar, ".list"), spearman.test.site.list, envir = .GlobalEnv)
    rm(spearman.test.site.list)
    rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
  }
  
  
  spearman.test.site.optB.list <- lapply(apropos("^spearman.test.siteB.*$", mode = "list"),
                                         get) %>%
    setNames(., names(usvi_asv_site.list))
  
  readr::write_rds(spearman.test.site.optB.list, paste0(projectpath, "/", "spearman.test.site.optB.list", ".rds"), compress = "gz")
  
}
rm(list = apropos("^(spearman.test.siteB.)(.*)$", mode = "list"))



# Spearman on site- and time-specific matrices ----------------------------


#note: since some of the metabolites were measured/recorded above detection only a handful of times in each site
#the correlations can be junky
#for each site x time, there will be at most 12 samples with ASVs and metabolomes.
#if we have fewer than 3 measurements of a metabolite across those 12 samples, do not attempt correlation calculation

#option A: 71 samples
if(file.exists(paste0(projectpath, "/", "spearman.test.site.time.optA.list", ".rds"))){
  spearman.test.site.time.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.time.optA.list", ".rds"))
} else {
  list2env(usvi_dist_mat_list_optA, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(-LODflag) %>%
    dplyr::relocate(sample_id) %>%
    dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(sample_id, simpleName) %>%
    dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE), 
                     num = length(conc),
                     .groups = "keep",
                     sd = sd(conc, na.rm = TRUE)) %>%
    dplyr::rename(conc = "mean_conc") %>%
    dplyr::mutate(log_conc = ifelse(!(is.na(conc) | (conc < 0)),
                                    log2(conc+1), #log transform abundance (with +1 pseudocount)
                                    # 0)) %>%
                                    NA)) %>%
    dplyr::select(sample_id, simpleName, log_conc) %>%
    droplevels %>%
    tidyr::pivot_wider(., id_cols = "sample_id",
                       values_from = "log_conc",
                       names_from = "simpleName") %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    droplevels
  
  usvi_asv.tbl <- usvi_prok_asvs.df %>%
    dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
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
  
  usvi_metab_site.time.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
    dplyr::mutate(grouping = paste0(site, ".", sampling_time)) %>%
    dplyr::distinct(sample_id, grouping) %>%
    droplevels %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(sample_id) %>%
          dplyr::left_join(., usvi_metab.tbl %>%
                             tibble::as_tibble(rownames = "sample_id")) %>%
          # dplyr::select(which(colSums(.) > 0)) %>% #remove any that are insignificant.
          tibble::column_to_rownames(var = "sample_id"))
  
  usvi_asv_site.time.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
    dplyr::mutate(grouping = paste0(site, ".", sampling_time)) %>%
    dplyr::distinct(sample_id, grouping) %>%
    droplevels %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(sample_id) %>%
          dplyr::inner_join(., usvi_asv.tbl %>%
                              tibble::as_tibble(rownames = "sample_id") %>%
                              droplevels) %>%
          tibble::column_to_rownames(var = "sample_id") %>%
          dplyr::select(which(colSums(.) > 0)) %>%
          droplevels)
  
  
  #for loop for each site and time:
  # for(grouping in names(usvi_asv_site.time.list)[3]){
  for(grouping in names(usvi_asv_site.time.list)){
    namevar <- grouping
    temp_asv.tbl <- usvi_asv_site.time.list[[grouping]]
    temp_metab.tbl <- usvi_metab_site.time.list[[grouping]]
    
    temp_spearman.test <- matrix(nrow = ncol(temp_asv.tbl), ncol = ncol(temp_metab.tbl))
    colnames(temp_spearman.test) <- colnames(temp_metab.tbl)
    rownames(temp_spearman.test) <- colnames(temp_asv.tbl)
    
    temp_spearman.test.rho <- temp_spearman.test
    temp_spearman.test.n <- temp_spearman.test
    
    y <- length(colnames(temp_spearman.test))
    
    # for(j in seq_len(y)){
    #   vector_metab <- temp_metab.tbl[, j]
    #   for(i in seq_len(nrow(temp_spearman.test))){
    #     vector_microb <- temp_asv.tbl[,i]
    #     vector_microb <- vector_microb[!is.na(vector_metab)]
    #     vector_metab_na <- vector_metab[!is.na(vector_metab)]
    #     
    #     if(length(vector_metab_na) >= 3 & sum(vector_microb) > 0){ #if 25% of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
    #       temp_spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
    #         purrr::pluck(., "p.value")
    #       temp_spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb, method = "spearman", exact = FALSE) %>%
    #         purrr::pluck(., "estimate")
    #       temp_spearman.test.n[i, j] <- length(vector_metab_na)
    #     } else {
    #       temp_spearman.test[i, j] <- NA
    #       temp_spearman.test.rho[i, j] <- NA
    #       if(length(vector_metab_na) < 3){
    #         temp_spearman.test.n[i, j] <- length(vector_metab_na) 
    #       } else {
    #         temp_spearman.test.n[i, j] <- length(which(vector_microb > 0))
    #       }
    #     }
    #   }
    # }
    
    for(j in seq_len(y)){
      vector_metab <- temp_metab.tbl[, j]
      vector_metab_na <- vector_metab[!is.na(vector_metab)]
      
      #if 25% of the samples have a NA value for that metabolite, stop the test for that metabolite against any ASV abundance.
      if(length(vector_metab_na) < 3){
        temp_spearman.test[, j] <- NA
        temp_spearman.test.rho[, j] <- NA
        temp_spearman.test.n[, j] <- length(vector_metab_na) 
      } else {
        for(i in seq_len(nrow(temp_spearman.test))){
          vector_microb <- temp_asv.tbl[,i]
          vector_microb_na <- vector_microb[!is.na(vector_metab)]  
          #if the total relative abudance of an ASV in samples with non-NA measurements of the metabolites is 0, stop the test and report the proportion of samples with non-zero relative abundance of that ASV
          if(sum(vector_microb_na) == 0){
            temp_spearman.test[i, j] <- NA
            temp_spearman.test.rho[i, j] <- NA
            temp_spearman.test.n[i, j] <- length(which(vector_microb > 0))/length(vector_microb)
          } else {
            if(length(vector_metab_na) >= 3 & sum(vector_microb_na) > 0){ #if 25% of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
              temp_spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb_na, method = "spearman", exact = FALSE) %>%
                purrr::pluck(., "p.value")
              temp_spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb_na, method = "spearman", exact = FALSE) %>%
                purrr::pluck(., "estimate")
              temp_spearman.test.n[i, j] <- length(vector_metab_na)
            }
          }
        }
      }
    }
    
    ##for troubleshooting:
    # assign(paste0("temp_spearman.test.", namevar), temp_spearman.test, envir = .GlobalEnv)
    # assign(paste0("temp_spearman.test.rho.", namevar), temp_spearman.test.rho, envir = .GlobalEnv)
    
    padj_cutoff <- temp_spearman.test %>%
      apply(., 2, na.omit) %>%
      unlist %>%
      ashr::qval.from.lfdr(.) %>%
      as.matrix(.) %>%
      quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
    
    temp_spearman.test.bh.corrected <- temp_spearman.test %>%
      apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
      apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
    
    temp_spearman.test.df <- temp_spearman.test.bh.corrected %>%
      tibble::as_tibble(., rownames = "asv_id") %>%
      tidyr::pivot_longer(., cols = !asv_id,
                          names_to = "simpleName",
                          values_to = "padj_bh") %>%
      droplevels %>%
      dplyr::left_join(., (temp_spearman.test %>%
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
      dplyr::right_join(., (temp_spearman.test.rho %>%
                              tibble::as_tibble(., rownames = "asv_id") %>%
                              tidyr::pivot_longer(., cols = !asv_id,
                                                  names_to = "simpleName",
                                                  values_to = "estimate")),
                        by = join_by(asv_id, simpleName)) %>%
      dplyr::mutate(sig = dplyr::case_when(!is.na(padj_bh) & !is.na(padj_05) ~ "sig",
                                           !is.na(padj_bh) & is.na(padj_05) ~ "maybe",
                                           .default = "not")) %>%
      dplyr::mutate(sig = factor(sig)) %>%
      dplyr::arrange(asv_id, simpleName) %>%
      dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      dplyr::mutate(label = signif(estimate, digits = 2)) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
      droplevels
    
    spearman.test.site.time.list <- list(temp_spearman.test.df, 
                                    temp_spearman.test, temp_spearman.test.rho, temp_spearman.test.n,
                                    # temp_dend_asv, temp_dend_metab, 
                                    padj_cutoff) %>%
      setNames(., c("spearman.test.df", 
                    "spearman.test", "spearman.test.rho", "spearman.test.n",
                    # "dend_asv", "dend_metab", 
                    "padj_cutoff"))
    
    assign(paste0("spearman.test.siteA.", namevar, ".list"), spearman.test.site.time.list, envir = .GlobalEnv)
    rm(spearman.test.site.time.list)
    rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
  }
  
  
  spearman.test.site.time.optA.list <- lapply(apropos("^spearman.test.siteA.*$", mode = "list"),
                                         get) %>%
    setNames(., names(usvi_asv_site.time.list))
  readr::write_rds(spearman.test.site.time.optA.list, paste0(projectpath, "/", "spearman.test.site.time.optA.list", ".rds"), compress = "gz")
  
}
rm(list = apropos("^(spearman.test.siteA.)(.*)$", mode = "list"))


#option B: 73 samples
if(file.exists(paste0(projectpath, "/", "spearman.test.site.time.optB.list", ".rds"))){
  spearman.test.site.time.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.time.optB.list", ".rds"))
} else {
  list2env(usvi_dist_mat_list_optB, envir = .GlobalEnv)
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, metab_deriv_label) %>%
                           droplevels),
                     by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::select(-LODflag) %>%
    dplyr::relocate(sample_id) %>%
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
                                    # 0)) %>%
                                    NA)) %>%
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
  
  usvi_metab_site.time.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::mutate(grouping = paste0(site, ".", sampling_time)) %>%
    dplyr::distinct(metab_deriv_label, grouping) %>%
    droplevels %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(metab_deriv_label) %>%
          dplyr::left_join(., usvi_metab.tbl %>%
                             tibble::as_tibble(rownames = "metab_deriv_label"),
                           by = join_by(metab_deriv_label)) %>%
          tibble::column_to_rownames(var = "metab_deriv_label"))
  
  usvi_asv_site.time.list <- metabolomics_sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::filter(metab_deriv_label %in% colnames(usvi_metab_mat)) %>%
    dplyr::mutate(grouping = paste0(site, ".", sampling_time)) %>%
    dplyr::distinct(metab_deriv_label, grouping) %>%
    droplevels %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(metab_deriv_label) %>%
          dplyr::inner_join(., usvi_asv.tbl %>%
                              tibble::as_tibble(rownames = "metab_deriv_label") %>%
                              droplevels, 
                            by = join_by(metab_deriv_label)) %>%
          tibble::column_to_rownames(var = "metab_deriv_label") %>%
          dplyr::select(which(colSums(.) > 0)) %>%
          droplevels)
  
  
  #for loop for each site x time:
  for(grouping in names(usvi_asv_site.time.list)){
    namevar <- grouping
    temp_asv.tbl <- usvi_asv_site.time.list[[grouping]]
    temp_metab.tbl <- usvi_metab_site.time.list[[grouping]]
    
    temp_spearman.test <- matrix(nrow = ncol(temp_asv.tbl), ncol = ncol(temp_metab.tbl))
    colnames(temp_spearman.test) <- colnames(temp_metab.tbl)
    rownames(temp_spearman.test) <- colnames(temp_asv.tbl)
    
    temp_spearman.test.rho <- temp_spearman.test
    temp_spearman.test.n <- temp_spearman.test
    
    y <- length(colnames(temp_spearman.test))
    
    for(j in seq_len(y)){
      vector_metab <- temp_metab.tbl[, j]
      vector_metab_na <- vector_metab[!is.na(vector_metab)]
      
      #if 25% of the samples have a NA value for that metabolite, stop the test for that metabolite against any ASV abundance.
      if(length(vector_metab_na) < 3){
        temp_spearman.test[, j] <- NA
        temp_spearman.test.rho[, j] <- NA
        temp_spearman.test.n[, j] <- length(vector_metab_na) 
      } else {
        for(i in seq_len(nrow(temp_spearman.test))){
          vector_microb <- temp_asv.tbl[,i]
          vector_microb_na <- vector_microb[!is.na(vector_metab)]  
          #if the total relative abudance of an ASV in samples with non-NA measurements of the metabolites is 0, stop the test and report the proportion of samples with non-zero relative abundance of that ASV
          if(sum(vector_microb_na) == 0){
            temp_spearman.test[i, j] <- NA
            temp_spearman.test.rho[i, j] <- NA
            temp_spearman.test.n[i, j] <- length(which(vector_microb > 0))/length(vector_microb)
          } else {
            if(length(vector_metab_na) >= 3 & sum(vector_microb_na) > 0){ #if 25% of the samples have a non-NA value for that metabolite, and the ASV was observed in those samples..
              temp_spearman.test[i, j] <- cor.test(vector_metab_na, vector_microb_na, method = "spearman", exact = FALSE) %>%
                purrr::pluck(., "p.value")
              temp_spearman.test.rho[i, j] <- cor.test(vector_metab_na, vector_microb_na, method = "spearman", exact = FALSE) %>%
                purrr::pluck(., "estimate")
              temp_spearman.test.n[i, j] <- length(vector_metab_na)
            }
          }
        }
      }
    }
    
    padj_cutoff <- temp_spearman.test %>%
      apply(., 2, na.omit) %>%
      unlist %>%
      ashr::qval.from.lfdr(.) %>%
      as.matrix(.) %>%
      quantile(., probs = seq(0.05, 0.1, 0.05), na.rm = TRUE, names = FALSE,type = 7) #get the possible p-adj cutoffs for different q-values
    
    temp_spearman.test.bh.corrected <- temp_spearman.test %>%
      apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
      apply(., 2, function(x) ifelse(x <= 0.05, x, NA)) #drop the p.values > the adjusted p-value or did not compute
    
    temp_spearman.test.df <- temp_spearman.test.bh.corrected %>%
      tibble::as_tibble(., rownames = "asv_id") %>%
      tidyr::pivot_longer(., cols = !asv_id,
                          names_to = "simpleName",
                          values_to = "padj_bh") %>%
      droplevels %>%
      dplyr::left_join(., (temp_spearman.test %>%
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
      dplyr::right_join(., (temp_spearman.test.rho %>%
                              tibble::as_tibble(., rownames = "asv_id") %>%
                              tidyr::pivot_longer(., cols = !asv_id,
                                                  names_to = "simpleName",
                                                  values_to = "estimate")),
                        by = join_by(asv_id, simpleName)) %>%
      dplyr::mutate(sig = dplyr::case_when(!is.na(padj_bh) & !is.na(padj_05) ~ "sig",
                                           !is.na(padj_bh) & is.na(padj_05) ~ "maybe",
                                           .default = "not")) %>%
      dplyr::mutate(sig = factor(sig)) %>%
      dplyr::arrange(asv_id, simpleName) %>%
      dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      dplyr::mutate(label = signif(estimate, digits = 2)) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
      droplevels
    
    spearman.test.site.time.list <- list(temp_spearman.test.df, 
                                    temp_spearman.test, temp_spearman.test.rho, temp_spearman.test.n,
                                    padj_cutoff) %>%
      setNames(., c("spearman.test.df", 
                    "spearman.test", "spearman.test.rho", "spearman.test.n",
                    "padj_cutoff"))
    
    assign(paste0("spearman.test.siteB.", namevar, ".list"), spearman.test.site.time.list, envir = .GlobalEnv)
    rm(spearman.test.site.time.list)
    rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
  }
  
  
  spearman.test.site.time.optB.list <- lapply(apropos("^spearman.test.siteB.*$", mode = "list"),
                                         get) %>%
    setNames(., names(usvi_asv_site.time.list))
  
  readr::write_rds(spearman.test.site.time.optB.list, paste0(projectpath, "/", "spearman.test.site.time.optB.list", ".rds"), compress = "gz")
  
}
rm(list = apropos("^(spearman.test.siteB.)(.*)$", mode = "list"))

rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
rm(list = apropos("^(temp_)(.*)(.tbl)$", mode = "list"))
rm(list = apropos("^(vector_m)(.*)(b)$", mode = "list"))
rm(list = apropos("^(vector_m)(.*)(b_na)$", mode = "list"))


# Compare what happens to sig relationships when we get granular ----------

#summary:
# using adjusted p-values for these granular datasets to filter for significance, may not be appropriate. 
#In all 7 comparisons, I observed that using a FDR q value = 0.10 would retain correlations where p-values < 0.1. 
#On the other hand, BH-adjusted p-values < 0.05 corresponded to a maximum raw p-value<0.0003213417 for the site-specific analyses, and maximum raw p.value<0.0001239889 for the site-and-time specific analyses, 
#which is really cutthroat.

#look at distribution of p-values
spearman.test.site.time.full.df <- list(spearman.test.site.time.optA.list, spearman.test.site.time.optB.list) %>%
  setNames(., c("optA", "optB")) %>%
  map(.,
      ~imap(.,
            ~.x[["spearman.test"]] %>%
              tibble::as_tibble(rownames = "asv_id") %>%
              tidyr::pivot_longer(., cols = -c("asv_id"),
                                  names_to = "simpleName",
                                  values_to = "p_value") %>%
              dplyr::mutate(grouping = .y)) %>%
        bind_rows(.)) %>%
  bind_rows(., .id = "test_type") %>%
  dplyr::left_join(., list(spearman.test.site.time.optA.list, spearman.test.site.time.optB.list) %>%
                     setNames(., c("optA", "optB")) %>%
                     map(.,
                         ~imap(.,
                               ~.x[["spearman.test.rho"]] %>%
                                 tibble::as_tibble(rownames = "asv_id") %>%
                                 tidyr::pivot_longer(., cols = -c("asv_id"),
                                                     names_to = "simpleName",
                                                     values_to = "estimate") %>%
                                 # dplyr::mutate(estimate = abs(estimate)) %>%
                                 dplyr::mutate(grouping = .y)) %>%
                           bind_rows(.)) %>%
                     bind_rows(., .id = "test_type"),
                   by = join_by(asv_id, test_type, grouping, simpleName)) %>%
  dplyr::mutate(padj_bh = p.adjust(p_value, "BH"))


# print(ggplot(data = spearman.test.site.time.full.df %>%
#                dplyr::filter(test_type == "optA" & grouping == "Yawzi.peak_photo") %>%
#                droplevels, aes(x = padj_bh, fill = test_type))
#       + scale_y_continuous(transform = "log10")
#       + geom_histogram( color = "black", binwidth = 0.01)
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(grouping~test_type, scales = "free_y", drop = TRUE)
#       )
# 
# print(ggplot(data = spearman.test.site.time.full.df %>%
#                dplyr::filter(test_type == "optA" & grouping == "Yawzi.peak_photo") %>%
#                tidyr::pivot_longer(., cols = c(p_value, padj_bh),
#                                    names_to = "metric",
#                                    values_to = "value") %>%
#                droplevels, 
#              aes(x = value, fill = metric))
#       + geom_histogram( color = "black", binwidth = 0.01, alpha = 0.5)
#       + scale_y_continuous(transform = "log10")
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(grouping~test_type, scales = "free_y", drop = TRUE)
# )
# print(ggplot(data = spearman.test.site.time.full.df %>%
#                dplyr::filter(test_type == "optA" & grouping == "Yawzi.peak_photo") %>%
#                tidyr::pivot_longer(., cols = c(p_value, padj_bh),
#                                    names_to = "metric",
#                                    values_to = "value") %>%
#                droplevels, 
#              aes(x = value, y= estimate, fill = metric))
#       + theme_bw()
#       + geom_point(shape = 21)
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(grouping~metric, scales = "free_y", drop = TRUE)
# )


#compare to the site-specific correlations
spearman.test.site.full.df <- list(spearman.test.site.optA.list, spearman.test.site.optB.list) %>%
  setNames(., c("optA", "optB")) %>%
  map(.,
      ~imap(.,
            ~.x[["spearman.test"]] %>%
              tibble::as_tibble(rownames = "asv_id") %>%
              tidyr::pivot_longer(., cols = -c("asv_id"),
                                  names_to = "simpleName",
                                  values_to = "p_value") %>%
              dplyr::mutate(grouping = .y)) %>%
        bind_rows(.)) %>%
  bind_rows(., .id = "test_type") %>%
  dplyr::left_join(., list(spearman.test.site.optA.list, spearman.test.site.optB.list) %>%
                     setNames(., c("optA", "optB")) %>%
                     map(.,
                         ~imap(.,
                               ~.x[["spearman.test.rho"]] %>%
                                 tibble::as_tibble(rownames = "asv_id") %>%
                                 tidyr::pivot_longer(., cols = -c("asv_id"),
                                                     names_to = "simpleName",
                                                     values_to = "estimate") %>%
                                 # dplyr::mutate(estimate = abs(estimate)) %>%
                                 dplyr::mutate(grouping = .y)) %>%
                           bind_rows(.)) %>%
                     bind_rows(., .id = "test_type"),
                   by = join_by(asv_id, test_type, grouping, simpleName)) %>%
  dplyr::mutate(padj_bh = p.adjust(p_value, "BH"))

#compare to the all-sample correlations
spearman.test.full.df <- list(spearman.test.optA.list, spearman.test.optB.list) %>%
  setNames(., c("optA", "optB")) %>%
  map(., ~.x[["spearman.test"]] %>%
        tibble::as_tibble(rownames = "asv_id") %>%
        tidyr::pivot_longer(., cols = -c("asv_id"),
                            names_to = "simpleName",
                            values_to = "p_value") %>%
        droplevels) %>%
  bind_rows(., .id = "test_type") %>%
  dplyr::left_join(., list(spearman.test.optA.list, spearman.test.optB.list) %>%
                     setNames(., c("optA", "optB")) %>%
                     map(., ~.x[["spearman.test.rho"]] %>%
                           tibble::as_tibble(rownames = "asv_id") %>%
                           tidyr::pivot_longer(., cols = -c("asv_id"),
                                               names_to = "simpleName",
                                               values_to = "estimate") %>%
                           # dplyr::mutate(estimate = abs(estimate)) %>%
                           droplevels) %>%
                     bind_rows(., .id = "test_type"),
                   by = join_by(asv_id, test_type, simpleName)) %>%
  dplyr::mutate(padj_bh = p.adjust(p_value, "BH")) %>%
  droplevels



# print(ggplot(data = spearman.test.full.df %>%
#                dplyr::filter(test_type == "optA") %>%
#                tidyr::pivot_longer(., cols = c(p_value, padj_bh),
#                                    names_to = "metric",
#                                    values_to = "value") %>%
#                droplevels, 
#              aes(x = value, y= estimate, fill = metric))
#       + theme_bw()
#       + geom_point(shape = 21)
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(.~metric, scales = "free_y", drop = TRUE)
# )
# 
# print(ggplot(data = spearman.test.full.df %>%
#                dplyr::filter(test_type == "optA") %>%
#                tidyr::pivot_longer(., cols = c(p_value, padj_bh),
#                                    names_to = "metric",
#                                    values_to = "value") %>%
#                droplevels, 
#              aes(x = value, fill = metric))
#       + geom_histogram( color = "black", binwidth = 0.01, alpha = 0.5)
#       + scale_y_continuous(transform = "log10")
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(.~test_type, scales = "free_y", drop = TRUE)
# )

#now combine themm
spearman.all.rho.p.df <- spearman.test.full.df %>%
  dplyr::filter(test_type == "optA") %>%
  dplyr::filter(p_value <= 0.1) %>%
  dplyr::distinct(asv_id, simpleName) %>%
  dplyr::mutate(grouping = "all") %>%
  droplevels %>%
  dplyr::inner_join(., spearman.test.site.time.full.df %>%
                      dplyr::filter(test_type == "optA") %>%
                      dplyr::filter(p_value <= 0.1) %>%
                      dplyr::distinct(asv_id, simpleName) %>%
                      droplevels, 
                    by = join_by(asv_id, simpleName)) %>%
  droplevels %>%
  dplyr::inner_join(., spearman.test.site.full.df %>%
                      dplyr::filter(test_type == "optA") %>%
                      dplyr::filter(p_value <= 0.1) %>%
                      dplyr::distinct(asv_id, simpleName) %>%
                      droplevels, 
                    by = join_by(asv_id, simpleName)) %>%
  droplevels %>%
  dplyr::left_join(., spearman.test.full.df %>%
                     dplyr::filter(test_type == "optA") %>%
                     # dplyr::mutate(grouping = "all") %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels, 
                   by = join_by(asv_id, simpleName)) %>%
  dplyr::bind_rows(., spearman.test.site.time.full.df %>%
                     dplyr::filter(test_type == "optA") %>%
                     dplyr::filter(p_value <= 0.1) %>%
                     dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels) %>%
  dplyr::bind_rows(., spearman.test.site.full.df %>%
                     dplyr::filter(test_type == "optA") %>%
                     dplyr::filter(p_value <= 0.1) %>%
                     dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels) %>%
  droplevels

spearman.all.rho.p.df <- spearman.all.rho.p.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, simpleName, grouping) %>%
  dplyr::filter(grepl("all", grouping), .by = c(asv_id, simpleName)) %>%
  dplyr::distinct(asv_id, simpleName) %>%
  dplyr::left_join(., spearman.all.rho.p.df,
                   by = join_by(asv_id, simpleName)) %>%
  tidyr::pivot_wider(., id_cols = c(asv_id, simpleName, grouping, test_type),
                     names_from = "metric",
                     values_from = "value") %>%
  dplyr::mutate(grouping = fct_relevel(grouping, "all", after = Inf)) %>%
  droplevels

# print(ggplot(data = spearman.all.rho.p.df %>%
#                droplevels, 
#              aes(x = p_value, y= estimate, fill = grouping))
#       + theme_bw()
#       + geom_point(shape = 21)
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(.~grouping, scales = "free_y", drop = TRUE)
# )


temp_g1 <- print(ggplot(data = spearman.all.rho.p.df %>%
                          dplyr::filter(grepl("all|LB", grouping)) %>%
                          dplyr::mutate(estimate = abs(estimate)) %>%
                          droplevels, 
                        aes(x = simpleName, y = estimate, fill = grouping))
                 + theme_bw()
                 + geom_point(shape = 21, show.legend = FALSE)
                 + geom_hline(yintercept = 0.5, color = "black")
                 + scale_y_continuous(name = "|Spearman's rho estimate|")
                 + facet_wrap(grouping~., nrow = 4, scales = "fixed", drop = FALSE, dir = "v")
                 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
temp_g2 <- print(ggplot(data = spearman.all.rho.p.df %>%
                          dplyr::filter(grepl("all|Yawzi", grouping)) %>%
                          dplyr::mutate(estimate = abs(estimate)) %>%
                          droplevels, 
                        aes(x = simpleName, y = estimate, fill = grouping))
                 + theme_bw()
                 + geom_point(shape = 21, show.legend = FALSE)
                 + geom_hline(yintercept = 0.5, color = "black")
                 + scale_y_continuous(name = "|Spearman's rho estimate|")
                 + facet_wrap(grouping~., nrow = 4, scales = "fixed", drop = FALSE, dir = "v")
                 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
temp_g3 <- print(ggplot(data = spearman.all.rho.p.df %>%
                          dplyr::filter(grepl("all|Tektite", grouping)) %>%
                          dplyr::mutate(estimate = abs(estimate)) %>%
                          droplevels, 
                        aes(x = simpleName, y = estimate, fill = grouping))
                 + theme_bw()
                 + geom_point(shape = 21, show.legend = FALSE)
                 + geom_hline(yintercept = 0.5, color = "black")
                 + scale_y_continuous(name = "|Spearman's rho estimate|")
                 + facet_wrap(grouping~., nrow = 4, scales = "fixed", drop = FALSE, dir = "v")
                 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)

gpatch <- temp_g1 | temp_g2 | temp_g3
gpatch <- gpatch + patchwork::plot_annotation(title = "Spearman's rho estimates for ASV-metabolite correlations where p < 0.10",
                                              subtitle = "Only those ASV-metabolite correlations consistently p < 0.10 are shown",
                                              tag_level = "A")
gpatch
ggsave(paste0(projectpath, "/", "spearman_sig_rho_all_vs_granular-", Sys.Date(), ".png"),
       gpatch,
       width = 12, height = 8, units = "in")

#in the all-sample spearman correlation analysis, for those correlations where p < 0.10, the range of Spearman's rho was 0-1
#of the ASV-metabolite correlations that were p < 0.10 in the all-sample, that were also in the site-specific correlations and p < 0.10, the range of Spearman's rho was smaller (0.35-1)
#finally, of the ASV-metabolite correlations that were p < 0.10 in the all-sample, that were also in the site- and time-specific correlations and p < 0.10, the range of Spearman's rho was even smaller (0.499-1)
# print(ggplot(data = spearman.test.full.df %>%
#                dplyr::filter(test_type == "optA") %>%
#                tidyr::pivot_longer(., cols = c(p_value, padj_bh),
#                                    names_to = "metric",
#                                    values_to = "value") %>%
#                droplevels, 
#              aes(x = value, fill = metric))
#       + geom_histogram( color = "black", binwidth = 0.01, alpha = 0.5)
#       + scale_y_continuous(transform = "log10")
#       + geom_vline(xintercept = 0.05, color = "black")
#       + facet_wrap(.~test_type, scales = "free_y", drop = TRUE)
# )
temp_g4 <- print(ggplot(data = spearman.test.full.df %>%
                          dplyr::filter(test_type == "optA") %>%
                                         tidyr::pivot_longer(., cols = c(p_value, padj_bh),
                                                             names_to = "metric",
                                                             values_to = "value") %>%
                                         droplevels,
                                       aes(x = value, fill = metric))
                 + theme_bw()
                 + geom_histogram( color = "black", alpha = 0.5)
                 + geom_vline(xintercept = 0.05, color = "black")
                 + scale_x_continuous(name = "p-value", 
                                      transform = scales::transform_boxcox(0.1),
                                      expand = c(0, 0))
                 + scale_y_continuous(transform = "log10", name = "Density", expand = c(0, 0), labels = scaleFUN0)
                 # + facet_wrap(~., nrow = 4, scales = "fixed", drop = FALSE, dir = "v")
                 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)),
                         legend.position = "bottom")
)
temp_g5 <- print(ggplot(data = spearman.test.site.full.df %>%
                          dplyr::filter(test_type == "optA") %>%
                          tidyr::pivot_longer(., cols = c(p_value, padj_bh),
                                              names_to = "metric",
                                              values_to = "value") %>%
                          droplevels,
                        aes(x = value, fill = metric))
                 + theme_bw()
                 + geom_histogram( color = "black", alpha = 0.5, show.legend = FALSE)
                 + geom_vline(xintercept = 0.05, color = "black")
                 + scale_x_continuous(name = "p-value", 
                                      transform = scales::transform_boxcox(0.1),
                                      expand = c(0, 0))
                 + scale_y_continuous(transform = "log10", name = "Density", expand = c(0, 0), labels = scaleFUN0)
                 + facet_wrap(grouping~., nrow = 4, scales = "fixed", drop = FALSE, dir = "v")
                 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
temp_g6 <- print(ggplot(data = spearman.test.site.time.full.df %>%
                          dplyr::filter(test_type == "optA") %>%
                          tidyr::pivot_longer(., cols = c(p_value, padj_bh),
                                              names_to = "metric",
                                              values_to = "value") %>%
                          droplevels,
                        aes(x = value, fill = metric))
                 + theme_bw()
                 + geom_histogram( color = "black", alpha = 0.5, show.legend = FALSE)
                 + geom_vline(xintercept = 0.05, color = "black")
                 + scale_x_continuous(name = "p-value", 
                                      transform = scales::transform_boxcox(0.1),
                                      expand = c(0, 0))
                 + scale_y_continuous(transform = "log10", name = "Density", expand = c(0, 0), labels = scaleFUN0)
                 + facet_wrap(grouping~., nrow = 3, scales = "fixed", drop = FALSE, dir = "h")
                 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
gpatch3_layout <- "
    AABBCCCC
    ##BBCCCC
    ##BBCCCC
  "

gpatch3 <- temp_g4 | temp_g5 | temp_g6
gpatch3 <- gpatch3 + patchwork::plot_layout(design =gpatch3_layout) + patchwork::plot_annotation(title = "Adjusted p-values for ASV-metabolite correlations",
                                                tag_levels = "A")
gpatch3
ggsave(paste0(projectpath, "/", "spearman_p_all_vs_granular-", Sys.Date(), ".png"),
       gpatch3,
       width = 16, height = 8, units = "in")


# SAR11-specific detour ---------------------------------------------------


#pull out SAR11-specific correlations
spearman_sar11.df <- spearman.test.site.time.full.df %>%
  dplyr::filter(asv_id %in% tibble::deframe(usvi_prok_asvs.taxa %>%
                                              dplyr::filter(if_any(everything(), ~grepl("SAR11", .x))) %>%
                                              dplyr::select(asv_id))) %>%
  bind_rows(., (spearman.test.site.full.df %>%
                  dplyr::filter(asv_id %in% tibble::deframe(usvi_prok_asvs.taxa %>%
                                                              dplyr::filter(if_any(everything(), ~grepl("SAR11", .x))) %>%
                                                              dplyr::select(asv_id))))) %>%
  bind_rows(., (spearman.test.full.df %>%
                  dplyr::filter(asv_id %in% tibble::deframe(usvi_prok_asvs.taxa %>%
                                                              dplyr::filter(if_any(everything(), ~grepl("SAR11", .x))) %>%
                                                              dplyr::select(asv_id))) %>%
                  dplyr::mutate(grouping = "all"))) %>%
  dplyr::filter(simpleName == "taurine") %>%
  tidyr::pivot_longer(., cols = c(p_value, padj_bh),
                      names_to = "metric",
                      values_to = "value") %>%
  dplyr::mutate(grouping = fct_relevel(grouping, "all", after = Inf)) %>%
  droplevels
#are there any in the more granular analyses that aren't in the all-sample analyses?
spearman_sar11_all_idx <- spearman_sar11.df %>%
  dplyr::filter(value <= 0.1) %>%
  dplyr::filter(grepl("all", grouping)) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id) %>%
  tibble::deframe(.)
spearman_sar11_granular_idx <- spearman_sar11.df %>%
  dplyr::filter(value <= 0.1) %>%
  dplyr::filter(!grepl("all", grouping)) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id) %>%
  tibble::deframe(.)
length(setdiff(spearman_sar11_granular_idx, spearman_sar11_all_idx)) #there are 9 SAR11 ASVs that have significant relationships to metabolites in the granular data, that are not in the all-sample results
length(setdiff(spearman_sar11_all_idx, spearman_sar11_granular_idx)) #there are 26 SAR11 ASVs that have significant relationships to metabolites in the all-sample results, that were not shared in the more granular results
length(unique(union(spearman_sar11_granular_idx, spearman_sar11_all_idx))) #74 SAR11 ASVs were significantly correlated with taurine in at least one grouping of data 


temp_g7 <- print(ggplot(data = spearman_sar11.df %>%
               dplyr::filter(test_type == "optA") %>%
                 dplyr::mutate(estimate = abs(estimate)) %>%
               droplevels,
             aes(x = value, y = estimate, fill = metric))
      + theme_bw()
      + geom_point(shape = 21, show.legend = TRUE)
      + geom_hline(yintercept = 0.5, color = "black")
      + geom_vline(xintercept = 0.1, color = "red")
      + scale_x_continuous(name = "p-value")
      + scale_y_continuous(name = "|Spearman's rho estimate|")
      # + facet_wrap(test_type~grouping, nrow = 4, scales = "fixed", drop = FALSE, dir = "h")
      + facet_wrap(.~grouping, nrow = 4, scales = "fixed", drop = FALSE, dir = "h")
      + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)

temp_g7 <- temp_g7 + patchwork::plot_annotation(title = "74 SAR11 ASVs with significant (p < 0.1) Spearman correlations to taurine concentrations",
                                                subtitle = "26 ASVs were significant in 'all' that were not in the granular subsets; 9 ASVs were significant in granular subsets but not in the 'all' dataset")
temp_g7
ggsave(paste0(projectpath, "/", "spearman_sar11_taurine_rho_all_vs_granular-", Sys.Date(), ".png"),
       temp_g7,
       width = 12, height = 8, units = "in")


temp_g8 <- print(ggplot(data = spearman.all.rho.p.df %>%
               droplevels)
      + theme_bw()
      + geom_histogram(aes(x = p_value), fill = "maroon", color = "black", binwidth = 0.01, alpha = 0.5, show.legend = TRUE)
      + geom_histogram(aes(x = padj_bh), fill = "seagreen", color = "black", binwidth = 0.01, alpha = 0.5, show.legend = TRUE)
      + geom_vline(xintercept = 0.05, color = "black")
      + scale_x_continuous(name = "p-values", 
                           transform = "identity")
      + scale_y_continuous(transform = "log10", name = "number of results", expand = c(0,0))
      + facet_wrap(.~grouping, nrow = 4, scales = "fixed", drop = FALSE, dir = "h")
)

temp_g9 <- print(ggplot(data = spearman.all.rho.p.df %>%
               droplevels,
             aes(x = p_value, y = padj_bh, fill = grouping))
      + theme_bw()
      + geom_point(shape = 21, show.legend = FALSE)
      + scale_x_continuous(name = "p-values", 
                           # transform = "log10", labels = scales::label_log(base = 10, digits = 5))
                           transform = scales::transform_boxcox(0.1))
      + scale_y_continuous(name = "adjusted p-values", 
                           # transform = "log10", labels = scales::label_log(base = 10, digits = 5))
                           transform = scales::transform_boxcox(0.1))
      # transform = "identity")
      + geom_hline(yintercept = 0.05, color = "red")
      + geom_vline(xintercept = 0.05, color = "black")
      + facet_wrap(.~grouping, nrow = 4, scales = "fixed", drop = FALSE, dir = "h")
)
gpatch2 <- temp_g8 + temp_g9 + patchwork::plot_annotation(title = "Adjusted p-values for ASV-metabolite correlations where p < 0.10",
                                                subtitle = "Only those ASV-metabolite correlations consistently p < 0.10 are shown",
                                                tag_levels = "A")
gpatch2
ggsave(paste0(projectpath, "/", "spearman_sig_p_all_vs_granular-", Sys.Date(), ".png"),
       gpatch2,
       width = 16, height = 8, units = "in")

# Filter the spearman correlations ----------------------------------------

#old method that used only BH-adjusted p-values to filter for significant correlations:
{
  # spearman.test.df <- list(spearman.test.optA.list[["spearman.test.df"]],
  #                          spearman.test.optB.list[["spearman.test.df"]]) %>%
  #   setNames(., c("optA", "optB")) %>%
  #   bind_rows(., .id = "test_type")
  # 
  # spearman.test.filtered.df <- spearman.test.df %>%
  #   tidyr::drop_na(padj_bh) %>%
  #   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
  #   dplyr::filter(consistent > 1) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::select(test_type, asv_id, simpleName, estimate) %>%
  #   tidyr::pivot_wider(., id_cols = NULL,
  #                      names_from = "test_type",
  #                      values_from = "estimate") %>%
  #   dplyr::mutate(consistent = dplyr::case_when((optA * optB) > 0 ~ 1,
  #                                               .default = NA)) %>%
  #   tidyr::drop_na(.) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::distinct(asv_id, simpleName) %>%
  #   dplyr::inner_join(., (spearman.test.df %>%
  #                          tidyr::drop_na(padj_bh) %>%
  #                          dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
  #                          dplyr::filter(consistent > 1) %>%
  #                          dplyr::select(test_type, asv_id, simpleName, padj_05) %>%
  #                          tidyr::pivot_wider(., id_cols = NULL,
  #                                             names_from = "test_type",
  #                                             values_from = "padj_05") %>%
  #                          tidyr::drop_na(.) %>%
  #                           dplyr::ungroup(.) %>%
  #                          dplyr::distinct(asv_id, simpleName) %>%
  #                           droplevels),
  #                    by = join_by(asv_id, simpleName)) %>%
  #   dplyr::left_join(., spearman.test.df %>%
  #                      dplyr::select(test_type, asv_id, simpleName, estimate, sig, label) %>%
  #                      dplyr::ungroup(.) %>%
  #                      droplevels,
  #                    by = join_by(asv_id, simpleName), relationship = "many-to-many", multiple = "first") %>%
  #   dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, .default = estimate)) %>%
  #   tidyr::drop_na(estimate) %>%
  #   dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
  #                                                      .default = NA)) %>%
  #   dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
  #   droplevels
  
  # (spearman.test.df %>%
  #     tidyr::drop_na(padj_bh) %>%
  #     dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
  #     dplyr::filter(consistent > 1) %>%
  #     dplyr::select(test_type, asv_id, simpleName, padj_05) %>%
  #     tidyr::pivot_wider(., id_cols = NULL,
  #                        names_from = "test_type",
  #                        values_from = "padj_05") %>%
  #     tidyr::drop_na(.) %>%
  #     dplyr::ungroup(.) %>%
  #     dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
  #     droplevels)
  
  
  #also site specific:
  
  # spearman.test.site.df <- list(spearman.test.site.optA.list, spearman.test.site.optB.list) %>%
  #   setNames(., c("optA", "optB")) %>%
  #   map(.,
  #       ~imap(.,
  #             ~.x[["spearman.test.df"]] %>%
  #               dplyr::mutate(site = .y)) %>%
  #         bind_rows(.)) %>%
  #   bind_rows(., .id = "test_type")
  # spearman.test.site.df <- spearman.test.site.df %>%
  #   dplyr::left_join(., (list(spearman.test.site.optA.list, spearman.test.site.optB.list) %>%
  #                          setNames(., c("optA", "optB")) %>%
  #                          map(.,
  #                              ~imap(.,
  #                                    ~.x[["spearman.test.n"]] %>%
  #                                      tibble::as_tibble(rownames = "asv_id") %>%
  #                                      dplyr::select(-asv_id) %>%
  #                                      dplyr::distinct(.) %>%
  #                                      dplyr::mutate(site = .y)) %>%
  #                                bind_rows(.)) %>%
  #                          bind_rows(., .id = "test_type")) %>%
  #                      tidyr::pivot_longer(., cols = !c(test_type, site),
  #                                          names_to = "simpleName",
  #                                          values_to = "num_obs_metab"),
  #                    by = join_by(test_type, site, simpleName))
  # 
  # #with site-specific, don't yet filter for whether the metabolite*ASV was found in both options A and B of handling metabolome replicate samples.
  # spearman.test.site.filtered.df <- spearman.test.site.df %>%
  #   tidyr::drop_na(padj_bh) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::distinct(asv_id, simpleName, site) %>%
  #   dplyr::inner_join(., (spearman.test.site.df %>%
  #                           tidyr::drop_na(padj_bh) %>%
  #                           dplyr::select(test_type, asv_id, simpleName, padj_05, site) %>%
  #                           tidyr::pivot_wider(., id_cols = NULL,
  #                                              names_from = "test_type",
  #                                              values_from = "padj_05") %>%
  #                           droplevels %>%
  #                           dplyr::distinct(asv_id, simpleName, site)),
  #                     by = join_by(asv_id, simpleName, site)) %>%
  #   dplyr::left_join(., spearman.test.site.df %>%
  #                      dplyr::select(test_type, asv_id, simpleName, estimate, sig, site, label, num_obs_metab) %>%
  #                      tidyr::pivot_wider(., id_cols = NULL,
  #                                         names_from = "test_type",
  #                                         values_from = "test_type") %>%
  #                      tidyr::unite("test_type", c(optA, optB), sep = "_", remove = TRUE, na.rm = TRUE) %>%
  #                      dplyr::arrange(abs(estimate)) %>%
  #                      dplyr::ungroup(.) %>%
  #                      dplyr::distinct(asv_id, simpleName, site, test_type, num_obs_metab, .keep_all = TRUE) %>%
  #                      droplevels,
  #                    by = join_by(asv_id, simpleName, site), relationship = "many-to-many", multiple = "all") %>%
  #   dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, .default = estimate)) %>%
  #   tidyr::drop_na(estimate) %>%
  #  dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
  #   droplevels
  
  # spearman.test.site.time.df <- list(spearman.test.site.time.optA.list, spearman.test.site.time.optB.list) %>%
  #   setNames(., c("optA", "optB")) %>%
  #   map(.,
  #       ~imap(.,
  #             ~.x[["spearman.test.df"]] %>%
  #               dplyr::mutate(grouping = .y)) %>%
  #         bind_rows(.)) %>%
  #   bind_rows(., .id = "test_type")
  # spearman.test.site.time.df <- list(spearman.test.site.time.optA.list, spearman.test.site.time.optB.list) %>%
  #   setNames(., c("optA", "optB")) %>%
  #   map(.,
  #       ~imap(.,
  #             ~.x[["spearman.test.df"]] %>%
  #               # tibble::as_tibble(rownames = "asv_id") %>%
  #               dplyr::mutate(grouping = .y)) %>%
  #         bind_rows(.)) %>%
  #   bind_rows(., .id = "test_type")
  # 
  # spearman.test.site.time.df <- spearman.test.site.time.df %>%
  #   dplyr::left_join(., (list(spearman.test.site.time.optA.list, spearman.test.site.time.optB.list) %>%
  #                          setNames(., c("optA", "optB")) %>%
  #                          map(.,
  #                              ~imap(.,
  #                                    ~.x[["spearman.test.n"]] %>%
  #                                      tibble::as_tibble(rownames = "asv_id") %>%
  #                                      # dplyr::select(-asv_id) %>%
  #                                      dplyr::distinct(.) %>%
  #                                      dplyr::mutate(grouping = .y)) %>%
  #                                bind_rows(.)) %>%
  #                          bind_rows(., .id = "test_type")) %>%
  #                      dplyr::distinct(test_type, asv_id, grouping, .keep_all = TRUE) %>%
  #                      tidyr::pivot_longer(., cols = !c(test_type, grouping, asv_id),
  #                                          names_to = "simpleName",
  #                                          values_to = "num_obs_metab"),
  #                    by = join_by(test_type, asv_id, grouping, simpleName))
  # 
  # spearman.test.site.time.filtered.df <- spearman.test.site.time.df %>%
  #   tidyr::drop_na(padj_bh) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::distinct(asv_id, simpleName, grouping) %>%
  #   dplyr::inner_join(., (spearman.test.site.time.df %>%
  #                           tidyr::drop_na(padj_bh) %>%
  #                           dplyr::select(test_type, asv_id, simpleName, padj_05, grouping) %>%
  #                           tidyr::pivot_wider(., id_cols = NULL,
  #                                              names_from = "test_type",
  #                                              values_from = "padj_05") %>%
  #                           droplevels %>%
  #                           dplyr::distinct(asv_id, simpleName, grouping)),
  #                     by = join_by(asv_id, simpleName, grouping)) %>%
  #   dplyr::left_join(., spearman.test.site.time.df %>%
  #                      dplyr::select(test_type, asv_id, simpleName, estimate, sig, grouping, label, num_obs_metab) %>%
  #                      tidyr::pivot_wider(., id_cols = NULL,
  #                                         names_from = "test_type",
  #                                         values_from = "test_type") %>%
  #                      tidyr::unite("test_type", c(optA, optB), sep = "_", remove = TRUE, na.rm = TRUE) %>%
  #                      dplyr::arrange(abs(estimate)) %>%
  #                      dplyr::ungroup(.) %>%
  #                      dplyr::distinct(asv_id, simpleName, grouping, test_type, num_obs_metab, .keep_all = TRUE) %>%
  #                      droplevels,
  #                    by = join_by(asv_id, simpleName, grouping), relationship = "many-to-many", multiple = "all") %>%
  #   dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, .default = estimate)) %>%
  #   tidyr::drop_na(estimate) %>%
  #   dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
  #   dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
  #                 sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
  #   droplevels
}
#so, instead of using BH adjusted p-values to filter for significant correlations, which seems only to
#be appropriate for the full dataset, use q.values

#filter for the correlations between ASVs and metabolites, that are consistently significant, in the same direction in options A and B if present in both

dend_metab <- list(spearman.test.optA.list[["dend_metab"]],
                                  spearman.test.optB.list[["dend_metab"]]) %>%
  setNames(., c("optA", "optB"))


#1. All samples used in correlation analyses:

# padj_cutoff <- list(spearman.test.optA.list[["padj_cutoff"]],
#      spearman.test.optB.list[["padj_cutoff"]]) %>%
#   setNames(., c("optA", "optB")) %>%
#   map(., ~.x %>% setNames(., c("q_05", "q_10")))

#if you want to recalculate the p-value threshold for other FDRs

padj_cutoff <- spearman.test.full.df %>%
  # dplyr::filter(abs(estimate) < 1) %>%
  # dplyr::rowwise(.) %>%
  split(., f = .$test_type) %>%
  map(., ~.x %>%
        dplyr::select(p_value) %>%
        tibble::deframe(.) %>% na.omit(.) %>%
        unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
        quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
        setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values

if(exists("padj_cutoff_list")){
  if(all(!grepl("spearman.test.full.df", names(padj_cutoff_list)))){
    padj_cutoff_list <- padj_cutoff_list %>%
      append(., (list(padj_cutoff) %>%
                   setNames(., c("spearman.test.full.df"))))
  } else {
    NULL
  }
} else {
  padj_cutoff_list <- list(padj_cutoff) %>%
    setNames(., c("spearman.test.full.df"))
}

# padj_cutoff
# # $optA
# # q_01        q_025         q_05         q_10 
# # 0.0000824359 0.0049043952 0.0317848477 0.0639543816 
# # 
# # $optB
# # q_01        q_025         q_05         q_10 
# # 0.0000772343 0.0047138053 0.0295009993 0.0621283345 
#there is a non-linear relationship between FDR cutoffs and p-value thresholds:
#the p-value threshold for a 2.5% FDR is approximately 1/10 of the 5% FDR
#the p-value for a 1% FDR is approximately 1/400 of the 5% FDR

#after discussion, we will stick to optA results only.
if(purrr::pluck_depth(padj_cutoff) > 2){
  cli::cli_alert_warning("Make sure your padj_cutoff list is appropriate for this step.")
} else {
  spearman.test.df <- spearman.test.full.df %>%
    tidyr::drop_na(p_value) %>%
    dplyr::filter(grepl("optA", test_type)) %>% #after discussion, we will stick to optA results only
    split(., f = .$test_type) %>%
    imap(., ~.x %>%
           dplyr::mutate(across(c(test_type, asv_id, simpleName), ~factor(.x))) %>%
           droplevels %>%
           dplyr::rowwise(.) %>%
           # dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                         # padj_10 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_10"] ~ p_value, .default = NA),
                         # padj_05 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_05"] ~ p_value, .default = NA),
                         # padj_025 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_025"] ~ p_value, .default = NA),
                         # padj_01 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_01"] ~ p_value, .default = NA)) %>%
           dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                         padj_01 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_01"] ~ p_value, .default = NA),
                         padj_025 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_025"] ~ p_value, .default = NA),
                         padj_05 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_05"] ~ p_value, .default = NA),
                         padj_10 = dplyr::case_when(p_value <= padj_cutoff[["optA"]]["q_10"] ~ p_value, .default = NA)) %>%
           tidyr::drop_na(padj_10) %>%
           dplyr::ungroup(.) %>%
           dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                     abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
           tidyr::drop_na(estimate) %>%
           dplyr::rowwise(.) %>%
           dplyr::mutate(sig = dplyr::case_when(
             !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
             !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
             !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
             !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
             !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
             .default = NA)) %>%
           dplyr::mutate(sig = factor(sig)) %>%
           dplyr::arrange(asv_id, simpleName) %>%
           dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
           dplyr::ungroup(.) %>%
           dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
           droplevels) %>%
    bind_rows(.)
  # spearman.test.df %>% dplyr::filter(sig %in% c("vsig", "sig_q01")) %>% dplyr::filter(abs(estimate) < 1) %>% dplyr::filter(abs(estimate) >= 0.5) %>% nrow(.)
  # spearman.test.df %>% dplyr::filter(sig %in% c("sig_q01")) %>% dplyr::filter(abs(estimate) < 1) %>% dplyr::filter(abs(estimate) >= 0.5) %>% nrow(.)
  #if you have optA and optB in spearman.test.site.df
  {
    # spearman.test.filtered.df <- spearman.test.df %>%
    #   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
    #   dplyr::filter(consistent > 1) %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::select(test_type, asv_id, simpleName, estimate) %>%
    #   tidyr::pivot_wider(., id_cols = NULL,
    #                      names_from = "test_type",
    #                      values_from = "estimate") %>%
    #   dplyr::mutate(consistent = dplyr::case_when((optA * optB) > 0 ~ 1,
    #                                               .default = NA)) %>%
    #   tidyr::drop_na(.) %>%
    #   dplyr::mutate(test_type = "optA_optB") %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::distinct(asv_id, simpleName, test_type) %>%
    #   dplyr::left_join(., spearman.test.df %>%
    #                      dplyr::ungroup(.) %>%
    #                      dplyr::filter(abs(estimate) < 1) %>%
    #                      dplyr::arrange(abs(estimate)) %>%
    #                      dplyr::select(asv_id, simpleName, estimate, sig),
    #                      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE),
    #                    by = join_by(asv_id, simpleName), relationship = "one-to-many", multiple = "first") %>%
    #   bind_rows(., (spearman.test.df %>%
    #                   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName")) %>%
    #                   dplyr::filter(consistent == 1) %>%
    #                   dplyr::ungroup(.) %>%
    #                   dplyr::distinct(test_type, asv_id, simpleName, estimate, sig))) %>%
    #   dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
    #                                                      .default = NA)) %>%
    #   dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
    #   dplyr::mutate(across(c(asv_id, simpleName, test_type, sig), ~factor(.x))) %>%
    #   droplevels  
    }
  
  
  spearman.test.filtered.df <- spearman.test.df %>%
    dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
    dplyr::ungroup(.) %>%
    dplyr::rowwise(.) %>%
    dplyr::filter(abs(estimate) < 1) %>%
    dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                       .default = NA)) %>%
    tidyr::drop_na(filtered_estimate) %>%
    # dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
    dplyr::mutate(across(c(asv_id, simpleName, test_type, sig), ~factor(.x))) %>%
    droplevels
}


#2. Site-specific samples used in correlation analyses:
#site-specific

# padj_cutoff <- list(spearman.test.site.optA.list,
#                     spearman.test.site.optB.list) %>%
#   setNames(., c("optA", "optB")) %>%
#   map_depth(., 2, ~.x %>% purrr::pluck("padj_cutoff") )

#if you want to recalculate the p-value threshold for other FDRs:

padj_cutoff <- spearman.test.site.full.df %>%
  dplyr::rowwise(.) %>%
  split(., f = .$test_type) %>%
  map(., ~.x %>%
        droplevels %>%
        split(., f = .$grouping) %>%
        map(., ~.x %>%
              droplevels %>%
              dplyr::select(p_value) %>%
              tibble::deframe(.) %>% na.omit(.) %>%
              unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
              quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
              setNames(., c("q_01", "q_025", "q_05", "q_10")))) #get the possible p-adj cutoffs for different q-values

if(exists("padj_cutoff_list")){
  if(all(!grepl("spearman.test.site.full.df", names(padj_cutoff_list)))){
    padj_cutoff_list <- padj_cutoff_list %>%
      append(., (list(padj_cutoff) %>%
                   setNames(., c("spearman.test.site.full.df"))))
  } else {
    NULL
  }
} else {
  padj_cutoff_list <- list(padj_cutoff) %>%
    setNames(., c("spearman.test.site.full.df"))
}
if(exists("spearman.test.site.filtered.df")){
  temp_df <- padj_cutoff %>%
    map(., ~.x %>%
          dplyr::bind_rows(., .id = "grouping")) %>%
    dplyr::bind_rows(., .id = "test_type") %>%
    # dplyr::left_join(., spearman.test.site.df %>%
    dplyr::left_join(., spearman.test.site.filtered.df %>%
                       dplyr::summarise(padj_01 = max(padj_01, na.rm = TRUE), .by = c("grouping", "test_type")),
                     by = join_by(test_type, grouping)) %>%
    tidyr::drop_na(padj_01) %>%
    droplevels
  if(!any(temp_df[["q_01"]] > temp_df[["padj_01"]])){
    cli::cli_alert_warning("Please reprocess the filtering step for correlation results.")
  }
  rm(temp_df)
} else {
  cli::cli_alert_warning("Please reprocess the filtering step for correlation results.")
  
  if(purrr::pluck_depth(padj_cutoff) < 2){
    cli::cli_alert_warning("Make sure your padj_cutoff list is appropriate for this step.")
  } else {
    spearman.test.site.df <- spearman.test.site.full.df %>%
      tidyr::drop_na(p_value) %>%
      dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels)
    spearman.test.site.df <- spearman.test.site.df %>%
      imap(., ~.x %>%
             dplyr::mutate(across(c(test_type, asv_id, simpleName, grouping), ~factor(.x))) %>%
             droplevels %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                           padj_01 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_01"] ~ p_value, .default = NA),
                           padj_025 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_025"] ~ p_value, .default = NA),
                           padj_05 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_05"] ~ p_value, .default = NA),
                           padj_10 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_10"] ~ p_value, .default = NA)) %>%
             tidyr::drop_na(padj_10) %>%
             dplyr::ungroup(.) %>%
             dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                       abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
             tidyr::drop_na(estimate) %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(sig = dplyr::case_when(
               !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
               !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
               .default = NA)) %>%
             dplyr::mutate(sig = factor(sig)) %>%
             dplyr::arrange(asv_id, simpleName) %>%
             dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
             dplyr::ungroup(.) %>%
             dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
             dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                           sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
             droplevels) %>%
      # setNames(., names(spearman.test.site.df)) %>%
      bind_rows(.)
    # spearman.test.site.df %>% dplyr::filter(sig %in% c("vsig", "sig_q01")) %>% dplyr::filter(abs(estimate) < 1) %>% dplyr::filter(abs(estimate) >= 0.5) %>% nrow(.)
    #if you have optA and optB in spearman.test.site.df
    {
      # spearman.test.site.df <- names(spearman.test.site.df) %>%
      #   # spearman.test.site.df <- c("optA", "optB") %>%
      #   imap(., ~spearman.test.site.df[[.x]] %>%
      #          dplyr::mutate(across(c(test_type, asv_id, simpleName, grouping), ~factor(.x))) %>%
      #          droplevels %>%
      #          dplyr::rowwise(.) %>%
      #          dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
      #                        padj_01 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_01"] ~ p_value, .default = NA),
      #                        padj_025 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_025"] ~ p_value, .default = NA),
      #                        padj_05 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_05"] ~ p_value, .default = NA),
      #                        padj_10 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_10"] ~ p_value, .default = NA)) %>%
      #          tidyr::drop_na(padj_10) %>%
      #          dplyr::ungroup(.) %>%
      #          dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
      #                                                    abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
      #          tidyr::drop_na(estimate) %>%
      #          dplyr::rowwise(.) %>%
      #          dplyr::mutate(sig = dplyr::case_when(
      #            !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
      #            !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
      #            !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
      #            !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
      #            !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
      #            .default = NA)) %>%
      #          dplyr::mutate(sig = factor(sig)) %>%
      #          dplyr::arrange(asv_id, simpleName) %>%
      #          dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      #          dplyr::ungroup(.) %>%
      #          dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
      #          dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
      #                        sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      #          droplevels) %>%
      #   setNames(., names(spearman.test.site.df)) %>%
      #   # setNames(., c("optA", "optB") ) %>%
      #   bind_rows(.)
      # spearman.test.site.filtered.df <- spearman.test.site.df %>%
      #   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName", "grouping")) %>%
      #   dplyr::filter(consistent > 1) %>%
      #   dplyr::ungroup(.) %>%
      #   dplyr::select(test_type, asv_id, simpleName, estimate, grouping) %>%
      #   tidyr::pivot_wider(., id_cols = NULL,
      #                      names_from = "test_type",
      #                      values_from = "estimate") %>%
      #   dplyr::mutate(consistent = dplyr::case_when((optA * optB) > 0 ~ 1,
      #                                               .default = NA)) %>%
      #   tidyr::drop_na(.) %>%
      #   dplyr::mutate(test_type = "optA_optB") %>%
      #   dplyr::ungroup(.) %>%
      #   dplyr::distinct(asv_id, simpleName, test_type, grouping) %>%
      #   dplyr::left_join(., spearman.test.site.df %>%
      #                      dplyr::ungroup(.) %>%
      #                      dplyr::filter(abs(estimate) < 1) %>%
      #                      dplyr::arrange(abs(estimate)) %>%
      #                      dplyr::select(asv_id, simpleName, estimate, grouping, sig),
      #                    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE),
      #                    by = join_by(asv_id, simpleName, grouping), relationship = "one-to-many", multiple = "first") %>%
      #   bind_rows(., (spearman.test.site.df %>%
      #                   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName" ,"grouping")) %>%
      #                   dplyr::filter(consistent == 1) %>%
      #                   dplyr::ungroup(.) %>%
      #                   dplyr::distinct(test_type, asv_id, simpleName, grouping, estimate, sig))) %>%
      #   dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
      #                                                      .default = NA)) %>%
      #   dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "grouping", "sig")) %>%
      #     dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
      #                   sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      #   dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      #   dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      #   dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      #   droplevels
      # 
      # spearman.test.site.filtered.df %>%
      #   # dplyr::filter(grepl("optA", test_type)) %>%
      #   tidyr::drop_na(filtered_estimate) %>%
      #   dplyr::filter(test_type != "optB") %>%
      #   droplevels %>%
      #   dplyr::group_by(site, sig) %>%
      #   dplyr::summarise(num_results = length(filtered_estimate))
      
      }
    
    
    spearman.test.site.filtered.df <- spearman.test.site.df %>%
      dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
      dplyr::ungroup(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
  }
  
}


spearman.test.site.df %>%
  dplyr::filter(grepl("optA", test_type)) %>%
  dplyr::filter(grepl("sig_q01|vsig", sig)) %>%
  dplyr::filter(abs(estimate) < 1) %>%
  droplevels %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(num_results = length(sig))

spearman.test.site.filtered.df %>%
  dplyr::filter(grepl("optA", test_type)) %>%
  dplyr::filter(grepl("sig_q01|vsig", sig)) %>%
  # tidyr::drop_na(padj_01) %>%
  droplevels %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(num_results = length(padj_01))

spearman.test.site.filtered.df %>%
  dplyr::filter(grepl("optA", test_type)) %>%
  # dplyr::filter(grepl("sig_q01|vsig", sig)) %>%
  tidyr::drop_na(padj_01) %>%
  droplevels %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(num_results = length(filtered_estimate))
# # A tibble: 3 × 2
# site        num_results
# <fct>             <int>
#   1 LB_seagrass         980
# 2 Yawzi               373
# 3 Tektite             442

#sig corrs: 980 Lameshur, 442 Tektite, 373 Yawzi
spearman_site_summary.df <- spearman.test.site.filtered.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, test_type, simpleName, site, sig, .keep_all = TRUE) %>%
  tidyr::drop_na(filtered_estimate) %>%
  dplyr::filter(sig %in% c("sig_q01", "vsig")) %>%
  dplyr::filter(grepl("optA", test_type)) %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        dplyr::ungroup(.) %>%
        # droplevels %>%
        # dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
        # dplyr::summarise(num_asvs = length(unique(.[["asv_id"]])),
        #                num_metabs = length(unique(.[["simpleName"]]))) %>%
        droplevels)
        # droplevels) %>% bind_rows(., .id = "site") %>% droplevels
# split(., f = .$asv_id) %>%
# map(., ~.x %>%
# dplyr::ungroup(.) %>%
# droplevels %>%
# dplyr::distinct(simpleName) %>% 
# dplyr::reframe(num_corrs = length(unique(.[["simpleName"]]))))
# dplyr::reframe(num_corrs = length(unique(.[["simpleName"]])))) %>% bind_rows(., .id = "asv_id")

spearman.test.site.filtered.df %>%
  dplyr::distinct(asv_id, test_type, simpleName, sig, grouping, .keep_all = TRUE) %>%
  dplyr::ungroup(.) %>%
  tidyr::drop_na(filtered_estimate) %>%
  dplyr::filter(grepl("sig_q01|vsig", sig)) %>%
  dplyr::filter(grepl("optA", test_type)) %>%
  dplyr::ungroup(.) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
  dplyr::distinct(asv_id, simpleName) %>%
  dplyr::reframe(num_asvs = length(unique(.[["asv_id"]])),
                 num_metabs = length(unique(.[["simpleName"]])))) %>%
  dplyr::bind_rows(., .id = "grouping") %>%
  dplyr::bind_rows(., (spearman.test.site.filtered.df %>%
                         dplyr::distinct(asv_id, test_type, simpleName, sig, grouping, .keep_all = TRUE) %>%
                         dplyr::ungroup(.) %>%
                         tidyr::drop_na(filtered_estimate) %>%
                         dplyr::filter(grepl("sig_q01|vsig", sig)) %>%
                         dplyr::filter(grepl("optA", test_type)) %>%
                         dplyr::ungroup(.) %>%
                         dplyr::distinct(asv_id, simpleName) %>%
                         dplyr::reframe(num_asvs = length(unique(.[["asv_id"]])),
                                        num_metabs = length(unique(.[["simpleName"]]))) %>%
                         dplyr::mutate(grouping = "total")))

# # A tibble: 4 × 3
# grouping    num_asvs num_metabs
# <chr>          <int>      <int>
#   1 LB_seagrass      324         42
# 2 Tektite          229         41
# 3 Yawzi            190         39
# 4 total            517         42

spearman_site_summary.df %>% dplyr::filter(num_corrs >= 10) %>% nrow(.)

# #we have this rds object saved: "usvi_spearman.sig.filtered.list-2025-04-07.rds"
# #it has 1239 significant correlations between ASVs and metabolites in Lameshur,
# #304 in Tektite, and 217 in Yawzi
# # temp_df <- dplyr::anti_join(`usvi_spearman.sig.filtered.list-2025-04-07`[["Tektite"]], spearman.test.site.filtered.df) %>%
# #   dplyr::bind_rows(., dplyr::anti_join(`usvi_spearman.sig.filtered.list-2025-04-07`[["Yawzi"]], spearman.test.site.filtered.df))
# #99 rows of data, which contains
# #43 entries for Yawzi, and
# #56 entries for Tektite
# # temp_df %>%
# #   tidyr::drop_na(padj_01) %>%
# #   dplyr::summarise(across(contains("padj_01"), list(min = min, max = max)), .by = c(grouping, site))
# # # A tibble: 2 × 4
# # grouping site    padj_01_min padj_01_max
# # <fct>    <fct>         <dbl>       <dbl>
# #   1 Tektite  Tektite      0.0105      0.0147
# # 2 Yawzi    Yawzi        0.0109      0.0144
# temp_df <- dplyr::anti_join(spearman.test.site.filtered.df, `usvi_spearman.sig.filtered.list-2025-04-07`[["Tektite"]]) %>%
#   dplyr::bind_rows(., dplyr::anti_join(spearman.test.site.filtered.df, `usvi_spearman.sig.filtered.list-2025-04-07`[["Yawzi"]]))
# #5435 rows of data:
# 
# temp_df %>%
#   tidyr::drop_na(padj_01) %>%
#   dplyr::summarise(across(contains("padj_01"), list(min = min, max = max)), .by = c(grouping, site))
# 
# # # A tibble: 3 × 4
# # grouping    site         padj_01_min padj_01_max
# # <fct>       <fct>              <dbl>       <dbl>
# #   1 LB_seagrass LB_seagrass 0.0000000198      0.0104
# # 2 Tektite     Tektite     0.0000314         0.0190
# # 3 Yawzi       Yawzi       0.0000556         0.0194
# 
# 
# #these minimum padj values are greater than the cutoff calculated in invididual sites:
# padj_cutoff[["optA"]][["LB_seagrass"]][["q_01"]]
# # 0.01036783
# padj_cutoff[["optA"]][["Yawzi"]][["q_01"]]
# #Yawzi: q_01 = 0.01941558
# #Tektite: q_01 = 0.01908171
# #it appears that I might have used the error cutoffs for Lameshur Bay, in Yawzi and Tektite correlations
# 

  



#3. Site- and time-specific samples used in correlation analyses:
#also site x time specific:
# padj_cutoff <- list(spearman.test.site.time.optA.list,
#                     spearman.test.site.time.optB.list) %>%
#   setNames(., c("optA", "optB")) %>%
#   map_depth(., 2, ~.x %>% purrr::pluck("padj_cutoff") )
padj_cutoff <- spearman.test.site.time.full.df %>%
  split(., f = .$test_type) %>%
  map(., ~.x %>%
        split(., f = .$grouping) %>%
        map(., ~.x %>%
              dplyr::select(p_value) %>%
              tibble::deframe(.) %>% na.omit(.) %>%
              unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
              quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
              setNames(., c("q_01", "q_025", "q_05", "q_10")))) #get the possible p-adj cutoffs for different q-values

if(exists("padj_cutoff_list")){
  if(all(!grepl("spearman.test.site.time.full.df", names(padj_cutoff_list)))){
    padj_cutoff_list <- padj_cutoff_list %>%
      append(., (list(padj_cutoff) %>%
                   setNames(., c("spearman.test.site.time.full.df"))))
  } else {
    NULL
  }
} else {
  padj_cutoff_list <- list(padj_cutoff) %>%
    setNames(., c("spearman.test.site.time.full.df"))
}

if(exists("spearman.test.site.time.filtered.df")){
  temp_df <- padj_cutoff %>%
    map(., ~.x %>%
          dplyr::bind_rows(., .id = "grouping")) %>%
    dplyr::bind_rows(., .id = "test_type") %>%
    # dplyr::left_join(., spearman.test.site.time.df %>%
    dplyr::left_join(., spearman.test.site.time.filtered.df %>%
                       dplyr::summarise(padj_01 = max(padj_01, na.rm = TRUE), .by = c("grouping", "test_type")),
                     by = join_by(test_type, grouping)) %>%
    tidyr::drop_na(padj_01) %>%
    droplevels
  if(!any(temp_df[["q_01"]] > temp_df[["padj_01"]])){
    cli::cli_alert_warning("Please reprocess the filtering step for correlation results.")
  }
  rm(temp_df)
} else {
  cli::cli_alert_info("Processing the filtering step for correlation results...")
  
  if(purrr::pluck_depth(padj_cutoff) < 2){
    cli::cli_alert_warning("Make sure your padj_cutoff list is appropriate for this step.")
  } else {
    spearman.test.site.time.df <- spearman.test.site.time.full.df %>%
      tidyr::drop_na(p_value) %>%
      dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels)
    spearman.test.site.time.df <- spearman.test.site.time.df %>%
      imap(., ~.x %>%
             dplyr::mutate(across(c(test_type, asv_id, simpleName, grouping), ~factor(.x))) %>%
             droplevels %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                           padj_01 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_01"] ~ p_value, .default = NA),
                           padj_025 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_025"] ~ p_value, .default = NA),
                           padj_05 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_05"] ~ p_value, .default = NA),
                           padj_10 = dplyr::case_when(p_value <= padj_cutoff[["optA"]][[.y]]["q_10"] ~ p_value, .default = NA)) %>%
             tidyr::drop_na(padj_10) %>%
             dplyr::ungroup(.) %>%
             dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                       abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
             tidyr::drop_na(estimate) %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(sig = dplyr::case_when(
               !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
               !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
               .default = NA)) %>%
             dplyr::mutate(sig = factor(sig)) %>%
             dplyr::arrange(asv_id, simpleName) %>%
             dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
             dplyr::ungroup(.) %>%
             dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
             dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                           sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
             droplevels) %>%
      # setNames(., names(spearman.test.site.df)) %>%
      bind_rows(.)
    
    
    spearman.test.site.time.filtered.df <- spearman.test.site.time.df %>%
      dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
      dplyr::ungroup(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
    
  }
  
  #if you have optA and optB test types still:
  {
    # spearman.test.site.time.df <- names(spearman.test.site.time.df) %>%
    #   # spearman.test.site.time.df <- c("optA", "optB") %>%
    #   imap(., ~spearman.test.site.time.df[[.x]] %>%
    #          dplyr::mutate(across(c(test_type, asv_id, simpleName, grouping), ~factor(.x))) %>%
    #          droplevels %>%
    #          dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
    #                        padj_01 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_01"] ~ p_value, .default = NA),
    #                        padj_025 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_025"] ~ p_value, .default = NA),
    #                        padj_05 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_05"] ~ p_value, .default = NA),
    #                        padj_10 = dplyr::case_when(p_value <= padj_cutoff[[.x]][[.y]]["q_10"] ~ p_value, .default = NA)) %>%
    #          tidyr::drop_na(padj_10) %>%
    #          dplyr::ungroup(.) %>%
    #          dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
    #                                                    abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
    #          tidyr::drop_na(estimate) %>%
    #          dplyr::mutate(sig = dplyr::case_when(
    #            !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
    #            !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
    #            !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
    #            !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
    #            !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
    #            .default = NA)) %>%
    #          dplyr::mutate(sig = factor(sig)) %>%
    #          dplyr::arrange(asv_id, simpleName) %>%
    #          dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
    #          dplyr::ungroup(.) %>%
    #          dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
    #          dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
    #                        sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
    #          droplevels) %>%
    #   setNames(., names(spearman.test.site.time.df)) %>%
    #   # setNames(., c("optA", "optB") ) %>%
    #   bind_rows(.)
    # spearman.test.site.time.filtered.df <- spearman.test.site.time.df %>%
    #   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName", "grouping")) %>%
    #   dplyr::filter(consistent > 1) %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::select(test_type, asv_id, simpleName, estimate, grouping) %>%
    #   tidyr::pivot_wider(., id_cols = NULL,
    #                      names_from = "test_type",
    #                      values_from = "estimate") %>%
    #   dplyr::mutate(consistent = dplyr::case_when((optA * optB) > 0 ~ 1,
    #                                               .default = NA)) %>%
    #   tidyr::drop_na(.) %>%
    #   dplyr::mutate(test_type = "optA_optB") %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::distinct(asv_id, simpleName, test_type, grouping) %>%
    #   dplyr::left_join(., spearman.test.site.time.df %>%
    #                      dplyr::ungroup(.) %>%
    #                      dplyr::filter(abs(estimate) < 1) %>%
    #                      dplyr::arrange(abs(estimate)) %>%
    #                      dplyr::select(asv_id, simpleName, estimate, grouping, sig),
    #                    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE),
    #                    by = join_by(asv_id, simpleName, grouping), relationship = "one-to-many", multiple = "first") %>%
    #   bind_rows(., (spearman.test.site.time.df %>%
    #                   dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName" ,"grouping")) %>%
    #                   dplyr::filter(consistent == 1) %>%
    #                   dplyr::ungroup(.) %>%
    #                   dplyr::distinct(test_type, asv_id, simpleName, grouping, estimate, sig))) %>%
    #   dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
    #                                                      .default = NA)) %>%
    #   dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "grouping", "sig")) %>%
    #   dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
    #                 sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
    #   dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
    #   dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
    #   dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
    #   droplevels
  }
}

if(!any(grepl("spearman_full", list.files(projectpath, pattern = "usvi_.*.RData")))){
  save(spearman.test.site.time.full.df, spearman.test.site.full.df, spearman.test.full.df,
       padj_cutoff_list,
       file = paste0(projectpath, "/", "usvi_spearman_full-", Sys.Date(), ".RData"))
}


if(!any(grepl("spearman.sig.filtered.list", list.files(projectpath, pattern = "usvi_.*.rds")))){
  spearman.test.filtered.lists <- list(spearman.test.site.time.filtered.df,
                    spearman.test.site.filtered.df,
                    spearman.test.filtered.df,
                    padj_cutoff_list) %>%
    setNames(., c("spearman.test.site.time.filtered.df",
                  "spearman.test.site.filtered.df",
                  "spearman.test.filtered.df",
                  "padj_cutoff_list"))
  readr::write_rds(spearman.test.filtered.lists, 
                   paste0(projectpath, "/", "usvi_spearman.sig.filtered.list-", Sys.Date(), ".rds"), 
                   compress = "gz")
  
  for(df in c("spearman.test.site.time.filtered.df",
              "spearman.test.site.filtered.df",
              "spearman.test.filtered.df")){
    temp_df <- get0(df, mode = "any")
    readr::write_delim(temp_df, paste0(projectpath, "/", "usvi_", df, "-", Sys.Date(), ".tsv"),
                       col_names = TRUE, delim = "\t")  
  }
  
  # save(spearman.test.site.time.strong.filtered.list,
  #      spearman.test.site.strong.filtered.list,
  #      spearman.test.strong.filtered.list,
  #      file = paste0(projectpath, "/", "usvi_spearman.sig.filtered.list-", Sys.Date(), ".RData"))
}


# output list of SDA/sig ASVs ---------------------------------------------

if(file.exists(paste0(projectpath, "/", "usvi_sig_seqs_phylogeny.df", ".tsv"))){
  usvi_sig_seqs_phylogeny.df <- readr::read_delim(paste0(projectpath, "/", "usvi_sig_seqs_phylogeny.df", ".tsv"), 
                                                  show_col_types = FALSE,
                                                  delim = "\t", col_names = TRUE)
  
} else {
  #write out the list of ASVs that were significant, to make a tree:
  #SDA ASVs:
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_sda_asvs_compare_summary-.*.tsv"))
  usvi_sda_asvs_compare_summary.df <- readr::read_delim(paste0(projectpath, "/", temp_file),
                                                        delim = "\t", col_names = TRUE, na = "", show_col_types = FALSE)
  rm(temp_file)
  
  
  # usvi_sig_seqs_idx <- spearman.test.filtered.df %>%  
  usvi_sig_seqs_idx <- bind_rows(spearman.test.filtered.df, spearman.test.site.filtered.df) %>%
    bind_rows(., spearman.test.site.time.filtered.df) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id, .keep_all = FALSE) %>%
    droplevels %>%
    dplyr::bind_rows(., usvi_sda_asvs_compare_summary.df %>%
                       dplyr::filter(test_type == "all") %>%
                       droplevels %>%
                       dplyr::distinct(asv_id, .keep_all = FALSE)) %>%
    dplyr::distinct(asv_id) %>%
    droplevels
  
  #2 additional ASVs found to be sig correlated site x time, that weren't in the SDA analyses, all Spearman, and site-specific Spearman
  # setdiff(unique(spearman.test.site.time.filtered.df[["asv_id"]]), usvi_sig_seqs_idx)
  # setdiff(unique(spearman.test.site.time.filtered.df[["asv_id"]]), unique(usvi_sda_asvs_compare_summary.df[["asv_id"]]))
  
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
  library(ape)
  library(stats)
  library(ggdendro)
  dada2::uniquesToFasta(usvi_sig_seqs_key.list, paste0(projectpath, "/", "usvi_prok_sig_asvs.fna"), ids = usvi_sig_seqs_key.list[["asv_id"]], mode = "w")
  
  #import the results from Silva
  keep_tax <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  
  usvi_sig_seqs_key.taxonomy <- usvi_prok_asvs.taxa %>%
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
    dplyr::distinct(asv_id, taxonomy, .keep_all = TRUE) %>%
    droplevels
  
  
  #option 1: use Silva tree to propagate taxonomy and relative arrangement of ASVs by taxonomy
  #look at the Silva tax tree and search for matching taxonomy:
  silva_ssu_tax.df <- readr::read_delim("~/projects/silva/tax_slv_ssu_138.2.txt", delim = "\t", col_names = FALSE, show_col_types = FALSE) %>%
    setNames(., c("taxonomy", "id", "tax_level", "X4", "version")) %>%
    dplyr::select(-X4) %>%
    tidyr::separate_wider_delim(taxonomy, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
    droplevels
  
  # usvi_sig_seqs_filtered_silva.df <- usvi_sig_seqs_key.taxonomy %>%
  #   dplyr::left_join(., silva_ssu_tax.df %>%
  #                      dplyr::select(-c(taxonomy)) %>%
  #                      droplevels, by = join_by(!!!keep_tax), relationship = "many-to-many", multiple = "all")
  # 
  usvi_sig_seqs_filtered_silva.df <- usvi_sig_seqs_key.taxonomy %>%
    dplyr::select(Domain:Genus) %>%
    dplyr::distinct(.) %>%
    dplyr::left_join(., silva_ssu_tax.df %>%
                       dplyr::select(id, tax_level, Domain:Genus) %>%
                       dplyr::distinct(.) %>%
                       droplevels, by = join_by(Domain, Phylum, Class, Order, Family, Genus), relationship = "many-to-many", multiple = "all") %>%
    tidyr::drop_na(.) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  
  silva_newick <- ape::read.tree("~/projects/silva/tax_slv_ssu_138.2.tre")
  # phyloseq::plot_tree(silva_newick, method = "sampledodge",
  #                     min.abundance = 5,
  #                     ladderize = TRUE)
  
  
  # usvi_sig_seqs_filtered_silva_idx <- usvi_sig_seqs_filtered_silva.df %>%
  #   dplyr::right_join(., usvi_sig_seqs_key.taxonomy %>%
  #                      dplyr::distinct(taxonomy, .keep_all = TRUE), relationship = "many-to-many", multiple = "all") %>%
  #   dplyr::select(asv_id, taxonomy, id) %>%
  #   dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  #   dplyr::select(asv_id, id) %>%
  #   tibble::deframe(.)
  #   # droplevels
  # 
  # silva_usvi_idx <- intersect(usvi_sig_seqs_filtered_silva_idx, silva_newick[["tip.label"]])
  
  # silva_newick_pruned <- ape::keep.tip(silva_newick, silva_usvi_idx)
  # phyloseq::plot_tree(silva_newick_pruned, method = "sampledodge",
  #                     min.abundance = 5,
  #                     ladderize = TRUE)
  # # usvi_sig_seqs_filtered_silva_idx %>%
  # #   tibble::enframe(name = "asv_id", value = "tip.label") %>%
  # #   dplyr::mutate(`tip.label` = as.character(`tip.label`)) %>%
  # #   dplyr::right_join(., tibble::tibble(`tip.label` = silva_usvi_idx))
  # # usvi_sig_seqs_filtered_silva_idx[match(usvi_sig_seqs_filtered_silva_idx, silva_usvi_idx)]
  # 
  # head(labels(silva_newick_pruned))
  # silva_usvi_renamed_idx <- usvi_sig_seqs_filtered_silva_idx[usvi_sig_seqs_filtered_silva_idx %in% silva_usvi_idx]
  # silva_usvi_renamed_idx <- names(silva_usvi_renamed_idx)
  # usvi_asvs_silva_newick_pruned <- silva_newick_pruned
  # labels(usvi_asvs_silva_newick_pruned) <- silva_usvi_renamed_idx
  # # phyloseq::plot_tree(usvi_asvs_silva_newick_pruned, method = "sampledodge",
  # #                     min.abundance = 5,
  #                     # ladderize = TRUE)
  # head(labels(usvi_asvs_silva_newick_pruned))
  # 
  # silva_newick_relabled_pruned <- silva_newick_pruned
  # labels(silva_newick_relabled_pruned)
  # silva_renamed_idx <- usvi_sig_seqs_filtered_silva.df %>%
  #   dplyr::right_join(., usvi_sig_seqs_key.taxonomy %>%
  #                       dplyr::distinct(taxonomy, .keep_all = TRUE), relationship = "many-to-many", multiple = "all") %>%
  #   dplyr::select(asv_id, taxonomy, id) %>%
  #   dplyr::mutate(taxonomy = paste0(asv_id, ": ", taxonomy)) %>%
  #   # dplyr::select(asv_id, id) %>%
  #   dplyr::select(taxonomy, id) %>%
  #   tibble::deframe(.)
  # silva_renamed_idx <- silva_renamed_idx[silva_renamed_idx %in% labels(silva_newick_relabled_pruned)]
  # silva_renamed_idx <- names(silva_renamed_idx)
  # labels(silva_newick_relabled_pruned) <- silva_renamed_idx
  # head(labels(silva_newick_relabled_pruned))
  # # phyloseq::plot_tree(usvi_asvs_silva_newick_pruned, method = "sampledodge",
  # #                     min.abundance = 5,
  # #                     ladderize = TRUE)
  
  # usvi_sig_seqs_phylogeny.df <- usvi_sig_seqs_filtered_silva.df %>%
  #   dplyr::mutate(id = factor(id, levels = labels(silva_newick))) %>%
  #   dplyr::arrange(id) %>%
  #     droplevels %>%
  #   dplyr::right_join(., usvi_sig_seqs_key.taxonomy, relationship = "many-to-many", multiple = "all") %>%
  #   dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, id) %>%
  #   tidyr::drop_na(Domain) %>%
  #   dplyr::mutate(across(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "tax_level", "asv_id", "taxonomy"), ~factor(.x))) %>%
  #   dplyr::mutate(arrangement = seq_len(nrow(.))) %>%
  #   droplevels
  
  #option 2:
  #made a tree with SINA (Arb web service) for 467 ASVs that, during analyses up to 1/30/2025, were found SDA
  #look at an ARB-generated tree
  
  usvi_arb.tree <- ape::read.tree("~/projects/apprill/usvi_temporal/arb-silva.de_2025-01-29_id1375389/arb-silva.de_2025-01-29_id1375389.tree")
  
  # phyloseq::plot_tree(usvi_arb.tree, method = "sampledodge",
  #                     min.abundance = 5,
  #                     ladderize = TRUE)
  # labels(usvi_arb.tree)
  # length(grep("ASV_", labels(usvi_arb.tree)))
  # usvi_arb.tree[["tip.labels"]]
  # plot(usvi_arb.tree)
  
  usvi_arb.taxonomy <- readr::read_delim("~/projects/apprill/usvi_temporal/arb-silva.de_2025-01-29_id1375389/arb-silva.de_align_resultlist_1375389.csv",
                                         delim = ";", col_names = TRUE, show_col_types = FALSE) %>%
    dplyr::select(sequence_identifier, identity, starts_with("lca_tax")) %>%
    tidyr::drop_na(.) %>%
    dplyr::mutate(across(starts_with("lca_tax"), ~dplyr::case_when(grepl("^Unclassified", .x) ~ NA,
                                                                   .default = .x))) %>%
    dplyr::mutate(across(starts_with("lca_tax"), ~stringr::str_remove_all(.x, '"'))) %>%
    dplyr::mutate(across(starts_with("lca_tax"), ~stringr::str_remove_all(.x, ';$')))
  
  #does the ARB taxonomy match our assignments?
  usvi_sig_seqs_arb.taxonomy <-  usvi_arb.taxonomy %>%
    dplyr::rename(asv_id = "sequence_identifier") %>%
    dplyr::select(asv_id, lca_tax_slv) %>%
    tidyr::separate_wider_delim(lca_tax_slv, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = TRUE) %>%
    droplevels %>%
    tidyr::pivot_longer(., cols = !c("asv_id"),
                        names_to = "level",
                        values_to = "arb") %>%
    dplyr::mutate(arb = gsub("ii", "i", arb)) %>%
    dplyr::left_join(., usvi_sig_seqs_key.taxonomy %>%
                       tidyr::pivot_longer(., cols = !c("asv_id"),
                                           names_to = "level",
                                           values_to = "dada2"),
                     by = join_by(asv_id, level)) %>%
    dplyr::group_by(asv_id, level) %>%
    dplyr::mutate(matched = dplyr::case_when((arb == dada2) ~ 1,
                                             # .default = NA)) %>% dplyr::group_by(asv_id) %>% tidyr::drop_na(.)
                                             .default = NA))
  # usvi_sig_seqs_arb.taxonomy %>%
  #   dplyr::group_by(asv_id) %>%
  #   dplyr::summarise(num_matched = sum(matched, na.rm = TRUE)) %>%
  #   dplyr::arrange(num_matched)
  
  usvi_arb_asvs_idx <- grep("ASV_", usvi_arb.tree[["tip.label"]])
  usvi_arb_pruned.tree <- ape::keep.tip(usvi_arb.tree, usvi_arb_asvs_idx)
  head(usvi_arb_pruned.tree[["tip.label"]])
  
  length(usvi_arb_pruned.tree[["tip.label"]])
  usvi_arb_renamed_idx <- data.frame(asv_id = usvi_arb_pruned.tree[["tip.label"]]) %>%
    dplyr::left_join(., usvi_prok_asvs.taxa %>%
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
                       dplyr::select(asv_id, taxonomy), 
                     by= join_by(asv_id)) %>%
    dplyr::mutate(taxonomy = across(starts_with("taxonomy")) %>% purrr::reduce(coalesce)) %>%
    dplyr::select(-ends_with(c(".x", ".y"))) %>%
    dplyr::filter(asv_id %in% usvi_arb_pruned.tree[["tip.label"]]) %>%
    dplyr::select(asv_id, taxonomy) %>%
    dplyr::mutate(taxonomy = paste0(asv_id, ": ", taxonomy)) %>%
    dplyr::select(taxonomy, asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = usvi_arb_pruned.tree[["tip.label"]])) %>%
    dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
    # droplevels %>% dplyr::arrange(asv_id)
    dplyr::arrange(asv_id) %>% tibble::deframe(.)
  
  usvi_arb_renamed_idx <- usvi_arb_renamed_idx[usvi_arb_renamed_idx %in% usvi_arb_pruned.tree[["tip.label"]]]
  head(names(usvi_arb_renamed_idx))
  length(names(usvi_arb_renamed_idx))
  usvi_arb_renamed_idx <- names(usvi_arb_renamed_idx)
  
  usvi_arb_pruned.tree[["tip.label"]] <- usvi_arb_renamed_idx
  
  # phyloseq::plot_tree(usvi_arb_pruned.tree, method = "sampledodge",
  #                     label.tips = "taxa_names",
  #                     ladderize = TRUE)
  
  usvi_arb_order_idx <- data.frame(tax_label = usvi_arb_pruned.tree[["tip.label"]]) %>%
    dplyr::mutate(asv_id = stringr::str_split_i(tax_label, ": ", 1)) %>%
    dplyr::mutate(arb_rangement = seq_len(nrow(.))) %>%
    droplevels
  
  
  #consolidate the ARB and Silva arrangements for the significant ASVs
  usvi_sig_seqs_phylogeny.df <- usvi_sig_seqs_filtered_silva.df %>%
    dplyr::mutate(id = factor(id, levels = silva_newick[["tip.label"]])) %>%
    dplyr::arrange(id) %>%
    droplevels %>%
    dplyr::right_join(., (usvi_sig_seqs_key.taxonomy %>%
                            dplyr::full_join(., data.frame(asv_id = grep("ASV_", usvi_arb.tree[["tip.label"]], value = TRUE)), by = join_by(asv_id))),
                      relationship = "many-to-many", multiple = "all") %>%
    #   dplyr::right_join(., usvi_sig_seqs_key.taxonomy, relationship = "many-to-many", multiple = "all") %>%
    dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, id) %>%
    tidyr::drop_na(Domain) %>%
    dplyr::mutate(across(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "tax_level", "asv_id", "taxonomy"), ~factor(.x))) %>%
    dplyr::mutate(arrangement = seq_len(nrow(.))) %>%
    dplyr::left_join(., usvi_arb_order_idx %>%
                       dplyr::select(asv_id, arb_rangement),
                     by = join_by(asv_id)) %>%
    dplyr::mutate(prop_arb_rangement = arb_rangement) %>%
    dplyr::arrange(arrangement, taxonomy) %>%
    dplyr::group_by(taxonomy) %>% tidyr::fill(prop_arb_rangement, .direction = "downup") %>%
    dplyr::ungroup(.) %>% tidyr::fill(prop_arb_rangement, .direction = "updown") %>%
    droplevels
  
  readr::write_delim(usvi_sig_seqs_phylogeny.df, paste0(projectpath, "/", "usvi_sig_seqs_phylogeny.df", ".tsv"), delim = "\t", col_names = TRUE)
  
}



# Combine Spearman correlation results from the 10 comparisons ------------

padj_cutoff_labels <- data.frame(sig = c("maybe", 
                                         # "sig", 
                                         "sig_q05",
                                         "sig_q025",
                                         "sig_q01",
                                         "vsig"),
                                   label = c("Q10", 
                                             # "Q05", 
                                             "Q05",
                                             "Q02.5",
                                             "Q01",
                                             "BH")) %>%
  tibble::deframe(.)

padj_cutoff_colors <- c("grey50", "tan", "gold", "lavender", "salmon") %>%
  setNames(., names(padj_cutoff_labels))

spearman.test.bh.dist.df <- bind_rows(spearman.test.df, spearman.test.site.df) %>%
  bind_rows(., spearman.test.site.time.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  tidyr::drop_na(padj_10) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  # tidyr::drop_na(padj_bh_05) %>%
  # dplyr::ungroup(.) %>%
  # dplyr::arrange(desc(padj_bh_05)) %>%
  dplyr::arrange(desc(padj_bh_05), .by_group = TRUE) %>%
  dplyr::select(grouping, site, sampling_time, p_value, padj_10, padj_05, padj_025, padj_01, padj_bh_05) %>%
  dplyr::distinct(grouping, site, sampling_time, .keep_all = TRUE) %>%
  droplevels
  
#plot the distribution of Spearman rho's when using all samples, compared to site-specific
spearman.test.rho.dist.df <- bind_rows(spearman.test.filtered.df, spearman.test.site.filtered.df) %>%
  bind_rows(., spearman.test.site.time.filtered.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  dplyr::distinct(asv_id, simpleName, grouping, site, sampling_time, sig, test_type, .keep_all = TRUE) %>%
  dplyr::select(asv_id, simpleName, grouping, site, sampling_time, estimate, sig, test_type) %>%
  dplyr::mutate(sig = factor(sig, levels = names(padj_cutoff_labels))) %>%
  dplyr::group_by(asv_id, simpleName,  grouping,site, sampling_time, sig, test_type) %>%
  # dplyr::distinct(estimate, .keep_all = TRUE) %>%
  droplevels

spearman.test.rho.dist.df %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(abs(estimate) >= 0.5) %>%
  # dplyr::distinct(asv_id,grouping, site, sampling_time, sig, test_type, .keep_all = TRUE) %>%
  dplyr::summarise(num_obs = length(asv_id), .by = c(grouping, site, sampling_time, sig)) %>%
  tidyr::complete(nesting(grouping, site, sampling_time), sig) %>%
  dplyr::mutate(num_obs = tidyr::replace_na(num_obs, 0)) %>%
  dplyr::arrange(site, sampling_time, sig) %>%
  dplyr::mutate(site = recode_factor(site, !!!site_lookup),
                sampling_time = recode_factor(sampling_time, !!!sampling_time_lookup)) %>%
  droplevels

spearman.test.rho.dist.summary.df <- spearman.test.rho.dist.df %>%
  dplyr::ungroup(.) %>%
  # dplyr::distinct(asv_id,grouping, site, sampling_time, sig, test_type, .keep_all = TRUE) %>%
  # dplyr::summarise(num_obs = length(asv_id), .by = c(grouping, site, sampling_time, sig)) %>%
  dplyr::summarise(num_obs = length(estimate), .by = c(grouping, site, sampling_time, sig)) %>%
  tidyr::complete(nesting(grouping, site, sampling_time), sig) %>%
  # dplyr::summarise(num_obs = length(estimate), .by = c(grouping, site, sampling_time, sig, test_type)) %>%
  # tidyr::complete(nesting(grouping, site, sampling_time, test_type), sig) %>%
  dplyr::mutate(num_obs = tidyr::replace_na(num_obs, 0)) %>%
  dplyr::arrange(site, sampling_time, sig) %>%
  dplyr::mutate(site = recode_factor(site, !!!site_lookup),
                sampling_time = recode_factor(sampling_time, !!!sampling_time_lookup)) %>%
  droplevels

print(spearman.test.rho.dist.summary.df, n = 50)
# # A tibble: 50 × 5
# grouping               site                  sampling_time sig      num_obs
# <chr>                  <fct>                 <fct>         <fct>      <int>
#   1 LB_seagrass            Lameshur Bay seagrass all           maybe        735
# 2 LB_seagrass            Lameshur Bay seagrass all           sig_q05     1625
# 3 LB_seagrass            Lameshur Bay seagrass all           sig_q025    1684
# 4 LB_seagrass            Lameshur Bay seagrass all           sig_q01      942
# 5 LB_seagrass            Lameshur Bay seagrass all           vsig          38
# 6 LB_seagrass.dawn       Lameshur Bay seagrass Dawn          maybe        451
# 7 LB_seagrass.dawn       Lameshur Bay seagrass Dawn          sig_q05      692
# 8 LB_seagrass.dawn       Lameshur Bay seagrass Dawn          sig_q025    1121
# 9 LB_seagrass.dawn       Lameshur Bay seagrass Dawn          sig_q01      655
# 10 LB_seagrass.dawn       Lameshur Bay seagrass Dawn          vsig           7
# 11 LB_seagrass.peak_photo Lameshur Bay seagrass Afternoon     maybe        531
# 12 LB_seagrass.peak_photo Lameshur Bay seagrass Afternoon     sig_q05      676
# 13 LB_seagrass.peak_photo Lameshur Bay seagrass Afternoon     sig_q025    1292
# 14 LB_seagrass.peak_photo Lameshur Bay seagrass Afternoon     sig_q01      823
# 15 LB_seagrass.peak_photo Lameshur Bay seagrass Afternoon     vsig           3
# 16 Tektite                Tektite Reef          all           maybe        945
# 17 Tektite                Tektite Reef          all           sig_q05     1440
# 18 Tektite                Tektite Reef          all           sig_q025     807
# 19 Tektite                Tektite Reef          all           sig_q01      410
# 20 Tektite                Tektite Reef          all           vsig           0
# 21 Tektite.dawn           Tektite Reef          Dawn          maybe        314
# 22 Tektite.dawn           Tektite Reef          Dawn          sig_q05      619
# 23 Tektite.dawn           Tektite Reef          Dawn          sig_q025     922
# 24 Tektite.dawn           Tektite Reef          Dawn          sig_q01      539
# 25 Tektite.dawn           Tektite Reef          Dawn          vsig           2
# 26 Tektite.peak_photo     Tektite Reef          Afternoon     maybe       1164
# 27 Tektite.peak_photo     Tektite Reef          Afternoon     sig_q05     1105
# 28 Tektite.peak_photo     Tektite Reef          Afternoon     sig_q025     579
# 29 Tektite.peak_photo     Tektite Reef          Afternoon     sig_q01      294
# 30 Tektite.peak_photo     Tektite Reef          Afternoon     vsig           3
# 31 Yawzi                  Yawzi Reef            all           maybe        762
# 32 Yawzi                  Yawzi Reef            all           sig_q05     1005
# 33 Yawzi                  Yawzi Reef            all           sig_q025    1044
# 34 Yawzi                  Yawzi Reef            all           sig_q01      406
# 35 Yawzi                  Yawzi Reef            all           vsig           0
# 36 Yawzi.dawn             Yawzi Reef            Dawn          maybe        642
# 37 Yawzi.dawn             Yawzi Reef            Dawn          sig_q05     1097
# 38 Yawzi.dawn             Yawzi Reef            Dawn          sig_q025     689
# 39 Yawzi.dawn             Yawzi Reef            Dawn          sig_q01      216
# 40 Yawzi.dawn             Yawzi Reef            Dawn          vsig           0
# 41 Yawzi.peak_photo       Yawzi Reef            Afternoon     maybe        348
# 42 Yawzi.peak_photo       Yawzi Reef            Afternoon     sig_q05      605
# 43 Yawzi.peak_photo       Yawzi Reef            Afternoon     sig_q025     888
# 44 Yawzi.peak_photo       Yawzi Reef            Afternoon     sig_q01      530
# 45 Yawzi.peak_photo       Yawzi Reef            Afternoon     vsig           5
# 46 all                    all                   all           maybe       5318
# 47 all                    all                   all           sig_q05     4646
# 48 all                    all                   all           sig_q025    2427
# 49 all                    all                   all           sig_q01        0
# 50 all                    all                   all           vsig        3727

#important to note, that the correlations tallied in "sig_q05", "sig_q025", and "sig_q01" also quality for "sig_q10" aka "maybe"
#similarly, "sig_q025" and "sig_q01" qualify for "sig_q05"
#and "sig_q025" qualifies for "sig_q01"
#sum them to their respective groups (e.g. we could therefore make a Venn diagram of these correlations based on their q thresholded p-values)

temp_df <- spearman.test.rho.dist.summary.df %>%
  dplyr::filter(!grepl("vsig", sig)) %>%
  dplyr::group_by(grouping, site, sampling_time) %>%
  dplyr::summarise(num_obs = sum(num_obs), .groups = "keep") %>%
  dplyr::mutate(sig = "maybe") %>%
  bind_rows(., (spearman.test.rho.dist.summary.df %>%
                  dplyr::filter(sig %in% c("sig_q05", "sig_q025", "sig_q01")) %>%
                  dplyr::group_by(grouping, site, sampling_time) %>%
                  dplyr::summarise(num_obs = sum(num_obs), .groups = "keep") %>%
                  dplyr::mutate(sig = "sig_q05"))) %>%
  bind_rows(., (spearman.test.rho.dist.summary.df %>%
                  dplyr::filter(sig %in% c("sig_q025", "sig_q01")) %>%
                  dplyr::group_by(grouping, site, sampling_time) %>%
                  dplyr::summarise(num_obs = sum(num_obs), .groups = "keep") %>%
                  dplyr::mutate(sig = "sig_q025"))) %>%
  bind_rows(., (spearman.test.rho.dist.summary.df %>%
                  dplyr::filter(grepl("sig_q01", sig)) %>%
                  droplevels)) %>%
  bind_rows(., (spearman.test.rho.dist.summary.df %>%
                  dplyr::filter(grepl("vsig", sig)) %>%
                  droplevels)) %>%
  dplyr::mutate(sig = factor(sig, levels = names(padj_cutoff_labels))) %>%
  droplevels

temp_g10 <- print(
  ggplot(data = temp_df, 
         aes(x = grouping, y = num_obs, fill = sig, group = grouping))
  + theme_bw()
  + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = TRUE, alpha = 0.7,
             position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
  + geom_text(aes(x = grouping, y = (num_obs+1)*1.2, label = num_obs, group = grouping),
              position = position_dodge2(width = 0.90, padding = 0.2, preserve = "total", reverse = TRUE))
  + scale_y_continuous(name = "Number of significant correlations", 
                       expand = expansion(mult = c(0.01,0.1)),
                       transform = scales::as.transform(t_pseudolog10))
  + scale_fill_manual(name = "Significance threshold", 
                      # values = c("grey50", "tan", "gold", "lavender", "salmon"), 
                      values = padj_cutoff_colors, breaks = names(padj_cutoff_colors),
                      labels = padj_cutoff_labels)
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)),
          legend.position = "bottom",
          axis.title.x = element_blank())
  )

temp_g10 <- temp_g10 + patchwork::plot_annotation(title = "Number of significant Spearman rank correlations between ASVs and metabolites",
                                      subtitle = "All samples, site-specific, and site- and time-specific datasets")
ggsave(paste0(projectpath, "/", "spearman_sig_threshold_dist-", Sys.Date(), ".png"),
       temp_g10,
       width = 10, height = 8, units = "in")

temp_g11 <- print(ggplot(data = spearman.test.rho.dist.df %>%
               dplyr::filter(grepl("LB", grouping)) %>%
               dplyr::mutate(estimate = abs(estimate)) %>%
               droplevels, 
             aes(x = simpleName, y = estimate, fill = sig))
      + theme_bw()
      + geom_point(shape = 21, show.legend = FALSE)
      + geom_hline(yintercept = 0.5, color = "black")
      + scale_y_continuous(name = "|Spearman's rho estimate|")
      + scale_fill_manual(name = "Significance threshold", 
                          # values = c("grey50", "tan", "gold", "lavender", "salmon"), 
                          values = padj_cutoff_colors,
                          breaks = names(padj_cutoff_colors),
                          labels = padj_cutoff_labels)
      + facet_grid(sig~grouping, scales = "fixed", drop = FALSE, 
                   labeller = labeller(sig = padj_cutoff_labels))
      + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
temp_g12 <- print(ggplot(data = spearman.test.rho.dist.df %>%
                           dplyr::filter(grepl("Yawzi", grouping)) %>%
                           dplyr::mutate(estimate = abs(estimate)) %>%
                           droplevels, 
                         aes(x = simpleName, y = estimate, fill = sig))
                  + theme_bw()
                  + geom_point(shape = 21, show.legend = FALSE)
                  + geom_hline(yintercept = 0.5, color = "black")
                  + scale_y_continuous(name = "|Spearman's rho estimate|")
                  + scale_fill_manual(name = "Significance threshold", 
                                      # values = c("grey50", "tan", "gold", "lavender", "salmon"), 
                                      values = padj_cutoff_colors,
                                      breaks = names(padj_cutoff_colors),
                                      labels = padj_cutoff_labels)
                  + facet_grid(sig~grouping, scales = "fixed", drop = FALSE, 
                               labeller = labeller(sig = padj_cutoff_labels))
                  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
temp_g13 <- print(ggplot(data = spearman.test.rho.dist.df %>%
                           dplyr::filter(grepl("Tektite", grouping)) %>%
                           dplyr::mutate(estimate = abs(estimate)) %>%
                           droplevels, 
                         aes(x = simpleName, y = estimate, fill = sig))
                  + theme_bw()
                  + geom_point(shape = 21, show.legend = FALSE, alpha = 0.5)
                  + geom_hline(yintercept = 0.5, color = "black")
                  + scale_y_continuous(name = "|Spearman's rho estimate|")
                  + scale_fill_manual(name = "Significance threshold", 
                                      # values = c("grey50", "tan", "gold", "lavender", "salmon"), 
                                      values = padj_cutoff_colors,
                                      breaks = names(padj_cutoff_colors),
                                      labels = padj_cutoff_labels)
                  + facet_grid(sig~grouping, scales = "fixed", drop = FALSE, 
                               labeller = labeller(sig = padj_cutoff_labels))
                  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
temp_g14 <- print(ggplot(data = spearman.test.rho.dist.df %>%
                           dplyr::filter(grepl("all", grouping)) %>%
                           dplyr::mutate(estimate = abs(estimate)), 
                         aes(x = simpleName, y = estimate, fill = sig))
                  + theme_bw()
                  + geom_point(shape = 21, show.legend = FALSE, alpha = 0.5)
                  + geom_hline(yintercept = 0.5, color = "black")
                  + scale_y_continuous(name = "|Spearman's rho estimate|")
                  + scale_fill_manual(name = "Significance threshold", 
                                      # values = c("grey50", "tan", "gold", "lavender", "salmon"), 
                                      values = padj_cutoff_colors,
                                      breaks = names(padj_cutoff_colors),
                                      labels = padj_cutoff_labels)
                  + facet_grid(sig~grouping, scales = "fixed", drop = FALSE, 
                               labeller = labeller(sig = padj_cutoff_labels))
                  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)))
)
gpatch4_layout <- "
    ABBBBCCCCDDDD
  "
gpatch4 <- temp_g14 | temp_g11 | temp_g12 | temp_g13 
gpatch4 <- gpatch4 + patchwork::plot_layout(design = gpatch4_layout) + patchwork::plot_annotation(title = "Spearman's rho estimates for ASV-metabolite correlations where p < 0.10",
                                                subtitle = "Only those ASV-metabolite correlations consistently p < 0.10 are shown",
                                                tag_level = "A")
gpatch4
ggsave(paste0(projectpath, "/", "spearman_p_fdr_dist-", Sys.Date(), ".png"),
       gpatch4,
       width = 28, height = 12, units = "in")



# Output results in tabular format ----------------------------------------


spearman.test.all.results.df <- bind_rows(spearman.test.df, spearman.test.site.df) %>%
  bind_rows(., spearman.test.site.time.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  tidyr::drop_na(padj_10) %>%
  dplyr::arrange(desc(padj_bh_05), .by_group = TRUE) %>%
  dplyr::select(test_type, grouping, site, sampling_time, asv_id, simpleName, estimate, sig, p_value, padj_10, padj_05, padj_025, padj_01, padj_bh_05) %>%
  droplevels
spearman.test.consistent.results.df <- spearman.test.all.results.df %>%
  dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName", "grouping", "site", "sampling_time")) %>%
  dplyr::filter(consistent > 1) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(test_type, asv_id, simpleName, estimate, grouping, site, sampling_time) %>%
  tidyr::pivot_wider(., id_cols = NULL,
                     names_from = "test_type",
                     values_from = "estimate") %>%
  dplyr::mutate(consistent = dplyr::case_when((optA * optB) > 0 ~ 1,
                                              .default = NA)) %>%
  tidyr::drop_na(.) %>%
  dplyr::mutate(test_type = "optA_optB") %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, asv_id, simpleName, grouping, site, sampling_time) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::left_join(., spearman.test.all.results.df %>%
                           # dplyr::select(test_type, asv_id, simpleName, estimate, grouping, site, sampling_time) %>%
                           droplevels,
                         by = join_by(asv_id, simpleName, grouping, site, sampling_time), relationship = "many-to-many", multiple = "all") %>%
        dplyr::summarise(mean_estimate = mean(estimate, na.rm = TRUE), 
                         across(c(p_value, padj_10, padj_05, padj_025, padj_01, padj_bh_05), list(max = max)),
                         .by = c(asv_id, simpleName, grouping, site, sampling_time)) %>%
        droplevels %>%
        dplyr::mutate(test_type = "optA_optB")) %>%
  bind_rows(., .id = NULL) %>%
  dplyr::relocate(test_type)
  
spearman.test.opt.results.df <- spearman.test.all.results.df %>%
  dplyr::mutate(consistent = length(test_type), .by = c("asv_id", "simpleName", "grouping", "site", "sampling_time")) %>%
  dplyr::filter(consistent < 2) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(-consistent) %>%
  droplevels

readr::write_delim(spearman.test.consistent.results.df, paste0(projectpath, "/", "spearman.test.consistent.results.df", ".tsv"),
                   delim = "\t", col_names = TRUE)
readr::write_delim(spearman.test.opt.results.df, paste0(projectpath, "/", "spearman.test.opt.results.df", ".tsv"),
                   delim = "\t", col_names = TRUE)
temp_df <- usvi_genera_relabel %>%
  tibble::enframe(name = "asv_id", value = "label")
readr::write_delim(temp_df, paste0(projectpath, "/", "usvi_genera_relabel", ".tsv"),
                   delim = "\t", col_names = TRUE)

temp_df <- spearman.test.opt.results.df %>%
# temp_df <- spearman.test.all.results.df %>%
  dplyr::ungroup(.) %>%
  # dplyr::filter(abs(estimate) >= 0.5) %>%
  dplyr::filter(grouping == "Tektite.dawn") %>%
  dplyr::filter( !is.na(padj_10) ) %>%
  dplyr::distinct(asv_id, simpleName, test_type, .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::summarise(across(starts_with("padj_"), ~list(length(na.omit(.x)))))
temp_df <- spearman.test.consistent.results.df %>%
  dplyr::ungroup(.) %>%
  # dplyr::filter(abs(mean_estimate) >= 0.5) %>%
  dplyr::filter(grouping == "Tektite.dawn") %>%
  dplyr::filter( !is.na(padj_10_max) ) %>%
  dplyr::distinct(asv_id, simpleName, test_type, .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::summarise(across(starts_with("padj_"), ~list(length(na.omit(.x)))))

# 
# print(ggplot(data = spearman.test.rho.dist.df, aes(x = estimate, fill = simpleName))
#       + geom_histogram( color = "black")
#       # + scale_discrete_manual(aesthetics = "fill", values = metab_colors, breaks = names(metab_colors), labels = names(metab_colors))
#       + facet_wrap(site~sampling_time, scales = "free_y", drop = TRUE)
#       )
# 
# print(ggplot(data = spearman.test.rho.dist.df %>%
#                dplyr::ungroup(.) %>%
#                # dplyr::distinct(asv_id, simpleName, site, .keep_all = TRUE) %>% dplyr::summarise(num_obs = length(simpleName), .by = c("asv_id", "site")))
#                dplyr::distinct(asv_id, site, .keep_all = TRUE) %>% dplyr::summarise(num_obs = length(asv_id), .by = c("site")))
#       + geom_col(aes(x = site, y= num_obs), color = "black")
# )


# Look at strongest correlations ------------------------------------------



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





dend_metab <- spearman.test.optA.list[["dend_metab"]]

spearman.sig.metab.list <- spearman.test.filtered.df %>%  
  dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
  tidyr::drop_na(estimate) %>%
  dplyr::select(simpleName, asv_id) %>%
  dplyr::mutate(simpleName = factor(simpleName, levels = (spearman.test.filtered.df %>%  
                                                            dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
                                                            tidyr::drop_na(estimate) %>%
                                                            dplyr::summarise(num_results = length(estimate), .by = "simpleName") %>%
                                                            dplyr::arrange(desc(num_results)) %>%
                                                            dplyr::select(simpleName) %>%
                                                            tibble::deframe(.)))) %>%
  split(., f = .$simpleName) %>%
  map(., ~.x %>%
        dplyr::distinct(asv_id, .keep_all = FALSE) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        dplyr::arrange(asv_id) %>%
        unlist %>%
        as.character)

spearman.sig.strong.metab.list <- spearman.test.filtered.df %>%  
  dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
  tidyr::drop_na(filtered_estimate) %>%
  dplyr::select(simpleName, asv_id) %>%
  dplyr::mutate(simpleName = factor(simpleName, levels = (spearman.test.filtered.df %>%  
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
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        dplyr::arrange(asv_id) %>%
        unlist %>%
        as.character)


spearman.sig_corr_asvs_idx_list <- NULL
spearman.sig_strong_corr_asvs_idx_list <- NULL
for(i in seq(2, 11, 1)){
  namevar <- paste0("sig_corr_asvs_", i, "_idx")
  temp_idx <- spearman.test.filtered.df %>%  
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
  temp_idx <- spearman.test.filtered.df %>%  
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


spearman.test.strong.filtered.list <- spearman.sig_strong_corr_asvs_idx_list %>%
  map(., ~dplyr::filter(spearman.test.filtered.df, asv_id %in% .x) %>%
        tidyr::drop_na(filtered_estimate) %>%
        dplyr::select(-filtered_estimate) %>%
        droplevels %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        # dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
        droplevels %>%
        dplyr::arrange(asv_id) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
        dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
        # dplyr::arrange(asv_id, simpleName) %>%
        droplevels)


# #look at site-specific correlations
# spearman.site.sig.metab.list <- spearman.test.site.filtered.df %>%  
#   dplyr::distinct(asv_id, simpleName, site, .keep_all = TRUE) %>% 
#   tidyr::drop_na(estimate) %>%
#   dplyr::select(simpleName, asv_id, site) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = (spearman.test.site.filtered.df %>%  
#                                                             dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
#                                                             tidyr::drop_na(estimate) %>%
#                                                             dplyr::summarise(num_results = length(estimate), .by = "simpleName") %>%
#                                                             dplyr::arrange(desc(num_results)) %>%
#                                                             dplyr::select(simpleName) %>%
#                                                             tibble::deframe(.)))) %>%
#   split(., f = .$simpleName) %>%
#   map(., ~.x %>%
#         dplyr::distinct(asv_id, .keep_all = FALSE) %>%
#         dplyr::mutate(asv_id = factor(asv_id, levels = usvi_sig_seqs_phylogeny.df[["asv_id"]])) %>%
#         dplyr::arrange(asv_id) %>%
#         unlist %>%
#         as.character)
# 
# #look at site- and time-specific correlations
# spearman.site.time.sig.metab.list <- spearman.test.site.time.filtered.df %>%  
#   dplyr::distinct(asv_id, simpleName, site, sampling_time, .keep_all = TRUE) %>% 
#   tidyr::drop_na(estimate) %>%
#   dplyr::select(simpleName, asv_id, site, sampling_time) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = (spearman.test.site.time.filtered.df %>%  
#                                                             dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
#                                                             tidyr::drop_na(estimate) %>%
#                                                             dplyr::summarise(num_results = length(estimate), .by = "simpleName") %>%
#                                                             dplyr::arrange(desc(num_results)) %>%
#                                                             dplyr::select(simpleName) %>%
#                                                             tibble::deframe(.)))) %>%
#   split(., f = .$simpleName) %>%
#   map(., ~.x %>%
#         dplyr::distinct(asv_id, .keep_all = FALSE) %>%
#         dplyr::mutate(asv_id = factor(asv_id, levels = usvi_sig_seqs_phylogeny.df[["asv_id"]])) %>%
#         dplyr::arrange(asv_id) %>%
#         unlist %>%
#         as.character)

# Summary of Spearman rank correlation coefficient tests ------------------

#if we filter for only the strongest correlations (abs(estimate) >= 0.5) between ASVs and metabolites
#we resolved 55 ASVs that each had 11+ correlations to 26 metabolites
#10+ correlations: 80 ASVs correlated
#9+: 107 ASVs 
#8+: 129 ASVs 
#7+: 153 ASVs 
#6+: 175 ASVs 
# # A tibble: 10 × 3
# num_corrs num_asvs num_metabs
# <dbl>    <int> <fct>     
#   1         2      253 31        
# 2         3      223 31        
# 3         4      202 30        
# 4         5      191 26        
# 5         6      175 26        
# 6         7      153 26        
# 7         8      129 26        
# 8         9      107 26        
# 9        10       80 26        
# 10        11       55 26        

spearman_summary.df <- spearman.test.strong.filtered.list %>%
  map(., ~.x %>%
        dplyr::distinct(asv_id, simpleName) %>%
        dplyr::reframe(num_asvs = length(unique(.[["asv_id"]])),
                       num_metabs = length(unique(.[["simpleName"]])))) %>%
  bind_rows(., .id = "num_corrs") %>%
  dplyr::mutate(num_corrs = gsub("([[:alpha:]_]+)([[:digit:]]{1,2})(_.*$)", "\\2", num_corrs)) %>%
  dplyr::mutate(num_metabs = factor(num_metabs)) %>%
  dplyr::mutate(num_corrs = as.numeric(num_corrs)) 

# g6 <- (ggplot(data = spearman_summary.df)
#        + theme_bw()
#        + geom_point(aes(x = num_corrs, y = num_asvs, fill = num_metabs), color = "black", shape = 21, size = 3, alpha = 1.0)
#        + scale_x_continuous(n.breaks = 10, name = "Minimum number of significant correlations for each ASV")
#        + scale_y_continuous(name = "Total ASVs with significant correlations above threshold", 
#                             expand = expansion(mult = c(0.1,0.1)))
#        + scale_fill_viridis(option = "mako", discrete = TRUE, name = "Number of \nobserved metabolites \n in correlations")
#        + ggtitle("Distribution of strong correlations for ASVs and metabolites")
#        + theme(panel.spacing = unit(1, "lines"),
#                panel.background = element_blank(),
#                axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
#                axis.text.y = element_text(vjust = 0.5, hjust = 1),
#                panel.grid.major = element_blank(),
#                # panel.grid.minor.y = element_blank(),
#                # panel.grid.minor.x = element_blank(),
#                panel.ontop = FALSE,
#                strip.text.y = element_blank(),
#                legend.position = "bottom")
# )
# g6
# 
# if(!any(grepl("spearman_strong_corr_summar", list.files(projectpath, pattern = "usvi_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_spearman_strong_corr_summary-", Sys.Date(), ".png"),
#          g6,
#          width = 6, height = 6, units = "in")
# }

if(!any(grepl("spearman_df", list.files(projectpath, pattern = "usvi_.*.RData")))){
  save(spearman.test.site.time.filtered.df, spearman.test.site.filtered.df, spearman.test.filtered.df, usvi_sig_seqs_phylogeny.df,
       spearman.site.time.sig.metab.list, spearman.site.sig.metab.list, spearman.sig.metab.list, spearman.sig.strong.metab.list, 
       file = paste0(projectpath, "/", "usvi_spearman_df-", Sys.Date(), ".RData"))
}

# Plot the heatmaps of correlations ---------------------------------------
# 

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
    # dplyr::mutate(asv_id = factor(asv_id, levels = c(labels(dend_asv), "total correlations"))) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = c(usvi_sig_seqs_phylogeny.df[["asv_id"]], "total correlations"))) %>%
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
  # if(!any(grepl("spearman_.*_strong_corr", list.files(projectpath, pattern = "usvi_.*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_spearman_", namevar, "_strong_corr-", Sys.Date(), ".png"),
           temp_g6,
           width = fig_width, height = 12, units = "in")
  # }
  rm(temp_g6)
  rm(title_plot)
  rm(namevar)
}

}



# Extract strong correlations for site- and sitextime specific ------------

#option A. average the metabolite profiles that are technical replicates, so there is 1 metabolome sample per microbiome sample (71 samples, LB: 23, Yawzi: 24, Tektite: 24)
#option B. keep the technical replicates, but duplicate the corresponding microbiome sample (73 samples, LB: 23, Yawzi: 25, Tektite: 25)



#plot the site-specific correlations

spearman.site.sig.filtered.list <- spearman.test.site.filtered.df %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        tidyr::drop_na(estimate) %>%
        droplevels %>%
        # dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        droplevels %>%
        dplyr::arrange(asv_id) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
        # dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
        droplevels)



rm(list = apropos(paste0("^g7_asvs_", ".*_.*_results$"), mode = "list"))

# for(i in c(1)){
for(i in seq_len(length(spearman.site.sig.filtered.list))){
  namevar <- names(spearman.site.sig.filtered.list)[i]
  title_plot <- recode(namevar, !!!site_lookup)
  temp_spearman.df <- spearman.site.sig.filtered.list[[i]] %>% split(., f = .$sig)
  
  temp_spearman.df2 <- temp_spearman.df %>%
    map(., ~.x %>%
          tibble::as_tibble(.) %>%
          dplyr::summarise(label = length(estimate), .by = "asv_id") %>%
          dplyr::mutate(simpleName = "total correlations")) 
  temp_spearman.df2 <- temp_spearman.df2 %>%
    imap(., ~.x %>% 
          dplyr::bind_rows(., temp_spearman.df[[.y]] %>%
                             dplyr::summarise(label = length(estimate), .by = "simpleName") %>%
                             dplyr::mutate(asv_id = "total correlations")))
  temp_spearman.df2 <- temp_spearman.df2 %>%
    imap(., ~.x %>%
           bind_rows(temp_spearman.df[[.y]], .) %>%
          # bind_rows(temp_spearman.df %>% dplyr::select(-label), .) %>%
          # dplyr::mutate(asv_id = factor(asv_id, levels = c(labels(dend_asv), "total correlations"))) %>%
          dplyr::mutate(asv_id = factor(asv_id, levels = c(usvi_sig_seqs_phylogeny.df[["asv_id"]], "total correlations"))) %>%
          dplyr::arrange(asv_id) %>%
          dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
          dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
          droplevels)
  
  for(j in seq_len(length(temp_spearman.df2))){
    sig_level <- names(temp_spearman.df2)[j]
    title_plot_sig <- paste0("ASVs with ", sig_level, " significant correlations in ", title_plot)
    fig_width <- max(10, (round(length(unique(temp_spearman.df2[[j]][["simpleName"]]))/10, digits = 0) + 8))
    fig_height <- max(8, (round(length(unique(temp_spearman.df2[[j]][["asv_id"]]))/100, digits = 0) + 4))
    temp_g7 <- (
      ggplot(data = temp_spearman.df2[[j]] %>%
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
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = TRUE,
              plot.title.position = "plot",
              title = element_text(hjust = 0),
              strip.text.y = element_blank())
      + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
                                   override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
      + ggtitle(title_plot_sig)
    )
    
    assign(paste0("g7_asvs_", namevar, "_", sig_level, "_results"), temp_g7, envir = .GlobalEnv)
    # if(!any(grepl(namevar, list.files(projectpath, pattern = "usvi_spearman_.*site_corr.*.png")))){
    # ggsave(paste0(projectpath, "/", "usvi_spearman_", namevar,"_", sig_level, "_corr-", Sys.Date(), ".png"),
    #        temp_g7,
    #        width = fig_width, height = fig_height,
    #        # height = 8, 
    #        units = "in")
    # }
    rm(temp_g7)
    rm(title_plot_sig)
  }
  # gpatch <- lapply(apropos(paste0("^g7_asvs_", namevar, "_.*_results$"), mode = "list"),
  #                  get) %>%
  #   purrr::reduce(., `/`) + 
  #   patchwork::plot_layout(guides = "collect") &
  #   theme(legend.position="none")
  # assign(paste0("g7_asvs_", namevar, "_all_results"), gpatch, envir = .GlobalEnv)
  
  rm(title_plot)
  rm(namevar)
}



#site x time specific
spearman.site.time.sig.filtered.list <- spearman.test.site.time.filtered.df %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        tidyr::drop_na(estimate) %>%
        droplevels %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        droplevels %>%
        dplyr::arrange(asv_id) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
        dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
        droplevels)


rm(list = apropos(paste0("^g8_asvs_", ".*_.*_results$"), mode = "list"))

# for(i in c(1)){
for(i in seq_len(length(spearman.site.time.sig.filtered.list))){
  namevar <- names(spearman.site.time.sig.filtered.list)[i]
  title_plot <- stringr::str_split_i(namevar, "\\.", 1) %>%
    recode(., !!!site_lookup)
  title_plot <- stringr::str_split_i(namevar, "\\.", 2) %>%
    recode(., !!!sampling_time_lookup) %>%
    paste0(title_plot, " at ", .)

  temp_spearman.df <- spearman.site.time.sig.filtered.list[[i]] %>% split(., f = .$sig)
  
  temp_spearman.df2 <- temp_spearman.df %>%
    map(., ~.x %>%
          tibble::as_tibble(.) %>%
          dplyr::summarise(label = length(estimate), .by = "asv_id") %>%
          dplyr::mutate(simpleName = "total correlations")) 
  temp_spearman.df2 <- temp_spearman.df2 %>%
    imap(., ~.x %>% 
           dplyr::bind_rows(., temp_spearman.df[[.y]] %>%
                              dplyr::summarise(label = length(estimate), .by = "simpleName") %>%
                              dplyr::mutate(asv_id = "total correlations")))
  temp_spearman.df2 <- temp_spearman.df2 %>%
    imap(., ~.x %>%
           bind_rows(temp_spearman.df[[.y]], .) %>%
           dplyr::mutate(asv_id = factor(asv_id, levels = c(usvi_sig_seqs_phylogeny.df[["asv_id"]], "total correlations"))) %>%
           dplyr::arrange(asv_id) %>%
           dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
           dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
           droplevels)
  
  for(j in seq_len(length(temp_spearman.df2))){
    sig_level <- names(temp_spearman.df2)[j]
    title_plot_sig <- paste0("ASVs with ", sig_level, " significant correlations in ", title_plot)
    fig_width <- max(10, (round(length(unique(temp_spearman.df2[[j]][["simpleName"]]))/10, digits = 0) + 8))
    fig_height <- max(8, (round(length(unique(temp_spearman.df2[[j]][["asv_id"]]))/100, digits = 0) + 4))
    temp_g8 <- (
      ggplot(data = temp_spearman.df2[[j]] %>%
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
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = TRUE,
              plot.title.position = "plot",
              title = element_text(hjust = 0),
              strip.text.y = element_blank())
      + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
                                   override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
      + ggtitle(title_plot_sig)
    )
    
    assign(paste0("g8_asvs_", namevar, "_", sig_level, "_results"), temp_g8, envir = .GlobalEnv)
    # if(!any(grepl(namevar, list.files(projectpath, pattern = "usvi_spearman_.*site_corr.*.png")))){
    # ggsave(paste0(projectpath, "/", "usvi_spearman_", namevar,"_", sig_level, "_corr-", Sys.Date(), ".png"),
    #        temp_g8,
    #        width = fig_width, height = fig_height,
    #        # height = 8,
    #        units = "in")
    # }
    rm(temp_g8)
    rm(title_plot_sig)
  }
  # gpatch <- lapply(apropos(paste0("^g7_asvs_", namevar, "_.*_results$"), mode = "list"),
  #                  get) %>%
  #   purrr::reduce(., `/`) + 
  #   patchwork::plot_layout(guides = "collect") &
  #   theme(legend.position="none")
  # assign(paste0("g7_asvs_", namevar, "_all_results"), gpatch, envir = .GlobalEnv)
  
  rm(title_plot)
  rm(namevar)
}

# Plot only the ASVs with the most correlations ---------------------------

#instead of plotting correlation estimates by significance level, summarize by density of correlations with ASVs

#first, filter for only those correlations where padj_05 is not NA
#then sumamrize number of ASVs with correlations to metabolites by site and sampling time


spearman.sig.site.time_corr_asvs_idx_list <- NULL #this represents ASVs with correlations to metabolites p < 0.05
spearman.sig.site.time_strong_corr_asvs_idx_list <- NULL #this represents ASVs with correlations to metabolites p < 0.01
for(i in seq(2, 11, 1)){
  namevar <- paste0("sig_corr_asvs_", i, "_idx")
  
  temp_idx <- spearman.test.site.time.filtered.df %>%  
    dplyr::filter(sig %in% c("sig_q05", "sig_q025", "sig_q01", "vsig")) %>%
    droplevels %>%
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
  spearman.sig.site.time_corr_asvs_idx_list <- append(spearman.sig.site.time_corr_asvs_idx_list, temp_idx)
  
  namevar <- paste0("sig_strong_corr_asvs_", i, "_idx")
  temp_idx <- spearman.test.site.time.filtered.df %>%  
    dplyr::filter(sig %in% c("sig_q01", "vsig")) %>%
    droplevels %>%
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
    tidyr::drop_na(estimate) %>%
    dplyr::summarise(num_results = length(estimate), .by = "asv_id") %>%
    # tidyr::drop_na(filtered_estimate) %>%
    # dplyr::summarise(num_results = length(filtered_estimate), .by = "asv_id") %>%
    dplyr::filter(num_results >= i) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id) %>%
    droplevels %>%
    tibble::deframe(.)
  temp_idx <- list(temp_idx) %>% 
    setNames(., namevar)
  spearman.sig.site.time_strong_corr_asvs_idx_list <- append(spearman.sig.site.time_strong_corr_asvs_idx_list, temp_idx)
  rm(temp_idx)
  rm(namevar)
}


spearman.test.site.time.strong.filtered.list <- spearman.sig.site.time_strong_corr_asvs_idx_list %>%
  map(., ~dplyr::filter(spearman.test.site.time.filtered.df, asv_id %in% .x) %>%
        dplyr::filter(sig %in% c("sig_q01", "vsig")) %>%
        droplevels %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                      sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                      grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
        droplevels %>%
        dplyr::arrange(asv_id) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
        droplevels)


spearman.sig.site_corr_asvs_idx_list <- NULL #this represents ASVs with correlations to metabolites p < 0.05
spearman.sig.site_strong_corr_asvs_idx_list <- NULL #this represents ASVs with correlations to metabolites p < 0.01
for(i in seq(2, 11, 1)){
  namevar <- paste0("sig_corr_asvs_", i, "_idx")
  
  temp_idx <- spearman.test.site.filtered.df %>%  
    dplyr::filter(sig %in% c("sig_q05", "sig_q025", "sig_q01", "vsig")) %>%
    droplevels %>%
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
  spearman.sig.site_corr_asvs_idx_list <- append(spearman.sig.site_corr_asvs_idx_list, temp_idx)
  
  namevar <- paste0("sig_strong_corr_asvs_", i, "_idx")
  temp_idx <- spearman.test.site.filtered.df %>%  
    dplyr::filter(sig %in% c("sig_q01", "vsig")) %>%
    droplevels %>%
    dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>% 
    tidyr::drop_na(estimate) %>%
    dplyr::summarise(num_results = length(estimate), .by = "asv_id") %>%
    # tidyr::drop_na(filtered_estimate) %>%
    # dplyr::summarise(num_results = length(filtered_estimate), .by = "asv_id") %>%
    dplyr::filter(num_results >= i) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(asv_id) %>%
    droplevels %>%
    tibble::deframe(.)
  temp_idx <- list(temp_idx) %>% 
    setNames(., namevar)
  spearman.sig.site_strong_corr_asvs_idx_list <- append(spearman.sig.site_strong_corr_asvs_idx_list, temp_idx)
  rm(temp_idx)
  rm(namevar)
}


spearman.test.site.strong.filtered.list <- spearman.sig.site_strong_corr_asvs_idx_list %>%
  map(., ~dplyr::filter(spearman.test.site.filtered.df, asv_id %in% .x) %>%
        dplyr::filter(sig %in% c("sig_q01", "vsig")) %>%
        droplevels %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
        dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                      sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                      grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
        droplevels %>%
        dplyr::arrange(asv_id) %>%
        dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
        droplevels)

if(!any(grepl("spearman.test.strong.filtered.list", list.files(projectpath, pattern = "usvi_.*.RData")))){
  save(spearman.test.site.time.strong.filtered.list,
       spearman.test.site.strong.filtered.list,
       spearman.test.strong.filtered.list,
       file = paste0(projectpath, "/", "usvi_spearman.test.strong.filtered.list-", Sys.Date(), ".RData"))
}



#what does it look liek foe the site-specific strongest correlatiosn with ASVs with 11 or more correlations?
{
#   temp_df <- spearman.test.site.strong.filtered.list[[10]]
# temp_df2 <- temp_df %>%
#         tibble::as_tibble(.) %>%
#   tidyr::complete(nesting(asv_id), grouping) %>%
#         dplyr::summarise(label = length(na.omit(estimate)), .by = c("asv_id", "grouping")) %>%
#         dplyr::mutate(simpleName = "total correlations") %>%
#   dplyr::bind_rows(., temp_df %>%
#                      tidyr::complete(nesting(simpleName), grouping) %>%
#                      dplyr::summarise(label = length(na.omit(estimate)), .by = c("simpleName", "grouping")) %>%
#                      dplyr::mutate(asv_id = "total correlations")) %>%
#          bind_rows(temp_df, .) %>%
#          dplyr::mutate(asv_id = factor(asv_id, levels = c(usvi_sig_seqs_phylogeny.df[["asv_id"]], "total correlations"))) %>%
#          dplyr::arrange(asv_id) %>%
#          dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#          dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
#   dplyr::ungroup(.) %>%
#          droplevels
# print(
#   ggplot(data = temp_df2 %>%
#            droplevels, aes(x = asv_id, y = simpleName))
#   + theme_bw() 
#   + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
#   + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
#   +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100), 
#                           transform = "reverse", aesthetics = "fill", 
#                           limits = c(1, -1), na.value = "white")
#   + scale_color_manual(values = c("grey", "black"), 
#                        breaks = c("not", "sig"), 
#                        labels = c("not", "sig"))
#   + scale_x_discrete(labels = usvi_genera_relabel,
#                      expand = c(0,0), name = "Taxon")
#   + scale_y_discrete(name = "Metabolite", expand = c(0,0))
#   + theme(panel.spacing = unit(1, "lines"),
#           panel.background = element_blank(),
#           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
#           axis.text.y = element_text(vjust = 0.5, hjust = 1),
#           panel.grid.major = element_blank(),
#           panel.grid.minor.y = element_blank(),
#           panel.grid.minor.x = element_blank(),
#           panel.ontop = TRUE,
#           # strip.text.y = element_blank(),
#           plot.title.position = "plot",
#           title = element_text(hjust = 0))
#   + facet_grid(grouping~., space = "fixed", scales = "free", drop = TRUE, 
#                labeller = labeller(grouping = site_lookup))
#   + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
#                                override.aes = list(stroke = 1, color = "black")),
#            color = "none")
#   + coord_flip()
# )
  }

#how much overlap between site-and-time-specific and site-specific and all-sample correlations?
strongest_idx <- "sig_strong_corr_asvs_11_idx"
spearman.test.sig_strong_corr_asvs_11.df <- bind_rows(spearman.test.site.time.strong.filtered.list[[strongest_idx]],
                                                      spearman.test.site.strong.filtered.list[[strongest_idx]],
                                                      spearman.test.strong.filtered.list[[strongest_idx]]) %>%
  dplyr::select(-filtered_estimate) %>%
  # dplyr::mutate(filtered_estimate = dplyr::case_when(is.na(filtered_estimate) & abs(estimate) >= 0.5 ~ estimate,
  #                                                    .default = filtered_estimate)) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  droplevels

{
  temp_df <- spearman.test.sig_strong_corr_asvs_11.df
  temp_df2 <- temp_df %>%
    tibble::as_tibble(.) %>%
    tidyr::complete(nesting(asv_id), grouping) %>%
    dplyr::summarise(label = length(na.omit(estimate)), .by = c("asv_id", "grouping")) %>%
    dplyr::mutate(simpleName = "total correlations") %>%
    dplyr::bind_rows(., temp_df %>%
                       tidyr::complete(nesting(simpleName), grouping) %>%
                       dplyr::summarise(label = length(na.omit(estimate)), .by = c("simpleName", "grouping")) %>%
                       dplyr::mutate(asv_id = "total correlations")) %>%
    bind_rows(temp_df, .) %>%
      dplyr::mutate(site = dplyr::case_when(is.na(site) & grouping == "all" ~ "all", 
                                            is.na(site) & grouping != "all" ~ stringr::str_split_i(grouping, "\\.", 1), 
                                            .default = site),
                    sampling_time = dplyr::case_when(is.na(sampling_time) & grepl("dawn|peak", grouping) ~ stringr::str_split_i(grouping, "\\.", 2),
                                                     is.na(sampling_time) & !grepl("dawn|peak", grouping) ~ "all",
                                                     .default = sampling_time)) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = c(usvi_sig_seqs_phylogeny.df[["asv_id"]], "total correlations"))) %>%
    dplyr::arrange(asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
    dplyr::mutate(simpleName = factor(simpleName, levels = c(labels(dend_metab), "total correlations"))) %>%
    dplyr::ungroup(.) %>%
    droplevels
  # temp_g9 <- print(
  #   ggplot(data = temp_df2 %>%
  #            dplyr::filter(grepl("all", grouping)) %>%
  #            droplevels, aes(x = asv_id, y = simpleName))
  #   + theme_bw()
  #   + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
  #   + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
  #   +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100),
  #                           transform = "reverse", aesthetics = "fill",
  #                           limits = c(1, -1), na.value = "white")
  #   + scale_color_manual(values = c("grey", "black"),
  #                        breaks = c("not", "sig"),
  #                        labels = c("not", "sig"))
  #   + scale_x_discrete(name = "Taxon",
  #                      # labels = usvi_genera_relabel,
  #                      expand = c(0,0))
  #   + scale_y_discrete(name = "Metabolite", expand = c(0,0))
  #   + theme(panel.spacing = unit(1, "lines"),
  #           panel.background = element_blank(),
  #           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
  #           axis.text.y = element_text(vjust = 0.5, hjust = 1),
  #           panel.grid.major = element_blank(),
  #           panel.grid.minor.y = element_blank(),
  #           panel.grid.minor.x = element_blank(),
  #           panel.ontop = TRUE,
  #           # strip.text.y = element_blank(),
  #           plot.title.position = "plot",
  #           title = element_text(hjust = 0))
  #   # + facet_wrap(grouping~., 
  #   + facet_grid(grouping~., space = "fixed",
  #                scales = "free", 
  #                # labeller = labeller(grouping = site_lookup)
  #                drop = TRUE
  #                )
  #   + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
  #                                override.aes = list(stroke = 1, color = "black")),
  #            color = "none")
  #   + coord_flip()
  # )
  # temp_g10 <- print(
  #   ggplot(data = temp_df2 %>%
  #            dplyr::filter(grepl("LB", grouping)) %>%
  #            droplevels, aes(x = asv_id, y = simpleName))
  #   + theme_bw()
  #   + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
  #   + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
  #   +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100),
  #                           transform = "reverse", aesthetics = "fill",
  #                           limits = c(1, -1), na.value = "white")
  #   + scale_color_manual(values = c("grey", "black"),
  #                        breaks = c("not", "sig"),
  #                        labels = c("not", "sig"))
  #   + scale_x_discrete(name = "Taxon",
  #                      # labels = usvi_genera_relabel,
  #                      expand = c(0,0))
  #   + scale_y_discrete(name = "Metabolite", expand = c(0,0))
  #   + theme(panel.spacing = unit(1, "lines"),
  #           panel.background = element_blank(),
  #           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
  #           axis.text.y = element_text(vjust = 0.5, hjust = 1),
  #           panel.grid.major = element_blank(),
  #           panel.grid.minor.y = element_blank(),
  #           panel.grid.minor.x = element_blank(),
  #           panel.ontop = TRUE,
  #           # strip.text.y = element_blank(),
  #           plot.title.position = "plot",
  #           title = element_text(hjust = 0))
  #   # + facet_wrap(grouping~., 
  #   + facet_grid(sampling_time~., space = "fixed",
  #                scales = "free", 
  #                # labeller = labeller(grouping = site_lookup)
  #                drop = TRUE
  #   )
  #   + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
  #                                override.aes = list(stroke = 1, color = "black")),
  #            color = "none")
  #   + coord_flip()
  # )
  # temp_g11 <- print(
  #   ggplot(data = temp_df2 %>%
  #            dplyr::filter(grepl("Yawzi", grouping)) %>%
  #            droplevels, aes(x = asv_id, y = simpleName))
  #   + theme_bw()
  #   + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
  #   + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
  #   +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100),
  #                           transform = "reverse", aesthetics = "fill",
  #                           limits = c(1, -1), na.value = "white")
  #   + scale_color_manual(values = c("grey", "black"),
  #                        breaks = c("not", "sig"),
  #                        labels = c("not", "sig"))
  #   + scale_x_discrete(name = "Taxon",
  #                      # labels = usvi_genera_relabel,
  #                      expand = c(0,0))
  #   + scale_y_discrete(name = "Metabolite", expand = c(0,0))
  #   + theme(panel.spacing = unit(1, "lines"),
  #           panel.background = element_blank(),
  #           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
  #           axis.text.y = element_text(vjust = 0.5, hjust = 1),
  #           panel.grid.major = element_blank(),
  #           panel.grid.minor.y = element_blank(),
  #           panel.grid.minor.x = element_blank(),
  #           panel.ontop = TRUE,
  #           # strip.text.y = element_blank(),
  #           plot.title.position = "plot",
  #           title = element_text(hjust = 0))
  #   # + facet_wrap(grouping~., 
  #   + facet_grid(grouping~., space = "fixed",
  #                scales = "free", 
  #                # labeller = labeller(grouping = site_lookup)
  #                drop = TRUE
  #   )
  #   + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
  #                                override.aes = list(stroke = 1, color = "black")),
  #            color = "none")
  #   + coord_flip()
  # )
  # temp_g12 <- print(
  #   ggplot(data = temp_df2 %>%
  #            dplyr::filter(grepl("Tektite", grouping)) %>%
  #            droplevels, aes(x = asv_id, y = simpleName))
  #   + theme_bw()
  #   + geom_tile(aes(fill = estimate), stat = "identity", color = "black", alpha = 0.7, show.legend = TRUE)
  #   + geom_text(aes(x = asv_id, y = simpleName, label = label), size = 3)
  #   +  scale_fill_gradientn(colors = colorRampPalette(pals::coolwarm(n = 3))(100),
  #                           transform = "reverse", aesthetics = "fill",
  #                           limits = c(1, -1), na.value = "white")
  #   + scale_color_manual(values = c("grey", "black"),
  #                        breaks = c("not", "sig"),
  #                        labels = c("not", "sig"))
  #   + scale_x_discrete(name = "Taxon",
  #                      # labels = usvi_genera_relabel,
  #                      expand = c(0,0))
  #   + scale_y_discrete(name = "Metabolite", expand = c(0,0))
  #   + theme(panel.spacing = unit(1, "lines"),
  #           panel.background = element_blank(),
  #           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = rel(0.8)),
  #           axis.text.y = element_text(vjust = 0.5, hjust = 1),
  #           panel.grid.major = element_blank(),
  #           panel.grid.minor.y = element_blank(),
  #           panel.grid.minor.x = element_blank(),
  #           panel.ontop = TRUE,
  #           # strip.text.y = element_blank(),
  #           plot.title.position = "plot",
  #           title = element_text(hjust = 0))
  #   # + facet_wrap(grouping~., 
  #   + facet_grid(grouping~., space = "fixed",
  #                scales = "free", 
  #                # labeller = labeller(grouping = site_lookup)
  #                drop = TRUE
  #   )
  #   + guides(fill = guide_legend(order = 2, ncol = 1, title = "Spearman estimate", direction = "vertical",
  #                                override.aes = list(stroke = 1, color = "black")),
  #            color = "none")
  #   + coord_flip()
  # )
}


# #if you wanted to de novo cluster the ASVs in this subset by taxonomic values...
# temp_dend <- usvi_sig_seqs_phylogeny.df %>%
#   tidyr::drop_na(id) %>%
#   dplyr::select(asv_id, Domain:Genus) %>%
#   tibble::column_to_rownames(var = "asv_id") %>%
#   dplyr::mutate(across(everything(), as.numeric)) %>%
#   dist(t(.), method = "euclidean") %>%
#   hclust(method = "ward.D2") %>%
#   as.dendrogram




# make networkg raph ------------------------------------------------------


temp_graph <- spearman.test.sig_strong_corr_asvs_11.df %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::filter(grepl("dawn", grouping)) %>%
  dplyr::select(asv_id, simpleName, estimate, site, sampling_time) %>%
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
# temp_graph_coords <- igraph::layout_nicely(temp_graph, dim = 2) %>%
#   tibble::as_tibble(.) %>%
#   dplyr::mutate(nodes = vertex_attr(temp_graph, "name")) %>%
#   dplyr::relocate(nodes) %>%
#   dplyr::mutate(type = dplyr::case_when(grepl("ASV_", nodes) ~ "taxon",
#                                         .default = "metabolite")) %>%
#   dplyr::left_join(., spearman.test.filtered.df %>%
#                      # dplyr::select(asv_id, group) %>%
#                      dplyr::distinct(asv_id, .keep_all = FALSE) %>%
#                      droplevels,
#                    by = join_by("nodes" =="asv_id"))
# 
# # print(
# #   ggplot(data = temp_graph_coords, aes(x = V1, y = V2, fill = type))
# #   + theme_bw()
# #   + geom_point(shape = 21)
# # )
# 
# temp_df <- temp_graph_coords %>%
#   dplyr::select(nodes, V1, V2) %>%
#   dplyr::distinct(., .keep_all = TRUE) %>%
#   dplyr::rename(xstart = "V1",
#                 node_start = "nodes",
#                 ystart = "V2") %>%
#   droplevels
# 
# temp_graph_edges <- spearman.test.filtered.df %>%
#   dplyr::filter(!is.na(rho_threshold) | !is.na(group)) %>%
#   dplyr::select(asv_id, simpleName, estimate) %>%
#   droplevels %>%
#   ggsankey::make_long(., asv_id, simpleName, value = "estimate") %>%
#   tidyr::drop_na(next_node) %>%
#   dplyr::rename(node_start = "node",
#                 node_end = "next_node",
#                 estimate = "value") %>%
#   dplyr::select(node_start, node_end, estimate) %>%
#   dplyr::distinct(., .keep_all = TRUE) %>%
#   droplevels %>%
#   dplyr::left_join(., temp_df,
#                    by = join_by(node_start)) %>%
#   dplyr::left_join(., temp_df %>%
#                      dplyr::rename(node_end = "node_start",
#                                    xend = "xstart",
#                                    yend = "ystart"),
#                    by = join_by(node_end)) %>%
#   # dplyr::left_join(., temp_graph_coords %>%
#   #   dplyr::select(nodes, type, color, group) %>%
#   #   droplevels,
#   # by = join_by("node_start" == "nodes")) %>%
#   dplyr::mutate(direction = dplyr::case_when((estimate < 0) ~ "negative",
#                                              (estimate > 0) ~ "positive",
#                                              .default = NA)) %>%
#   dplyr::mutate(abs_estimate = abs(estimate)) %>%
#   droplevels
# 
# 
# #plot the taxa with higher abundance in reefs than in seagrass, and their metabolite correlations:
# 
# keep <- spearman.test.filtered.df %>%
#   dplyr::filter(grepl("high", group)) %>%
#   dplyr::filter(!is.na(rho_threshold)) %>%
#   dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
#   droplevels %>%
#   tidyr::pivot_longer(., cols = everything(),
#                       names_to = NULL,
#                       values_to = "nodes") %>%
#   dplyr::arrange(nodes) %>%
#   distinct(.) %>%
#   unlist(.)
# 
# temp_df2 <- temp_graph_coords %>%
#   dplyr::filter(grepl(paste0(keep, collapse = "|"), nodes)) %>%
#   dplyr::mutate(nodes = dplyr::case_when(grepl("ASV_", nodes) ~ dplyr::recode_factor(nodes, !!!usvi_genera_relabel),
#                                          .default = nodes)) %>%
#   dplyr::mutate(nodes = stringr::str_split_i(nodes, ";", -1)) %>%
#   droplevels
# 
# temp_df1 <- temp_graph_edges %>%
#   dplyr::right_join(., spearman.test.filtered.df %>%
#                       dplyr::filter(grepl("high", group)) %>%
#                       dplyr::filter(!is.na(rho_threshold)) %>%
#                       dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
#                       droplevels,
#                     by = join_by("node_start" == "asv_id", "node_end" == "simpleName")) %>%
#   dplyr::mutate(node_start = dplyr::case_when(grepl("ASV_", node_start) ~ dplyr::recode_factor(node_start, !!!usvi_genera_relabel),
#                                               .default = node_start)) %>%
#   dplyr::mutate(node_start = stringr::str_split_i(node_start, ";", -1)) %>%
#   droplevels
# 
# 
# 
# #now do the lower-abudnace taxa and their metabolite correlations
# 
# keep <- spearman.test.filtered.df %>%
#   dplyr::filter(grepl("low", group)) %>%
#   dplyr::filter(!is.na(rho_threshold)) %>%
#   dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
#   droplevels %>%
#   tidyr::pivot_longer(., cols = everything(),
#                       names_to = NULL,
#                       values_to = "nodes") %>%
#   dplyr::arrange(nodes) %>%
#   distinct(.) %>%
#   unlist(.)
# 
# temp_df4 <- temp_graph_coords %>%
#   dplyr::filter(grepl(paste0(keep, collapse = "|"), nodes)) %>%
#   dplyr::mutate(nodes = dplyr::case_when(grepl("ASV_", nodes) ~ dplyr::recode_factor(nodes, !!!usvi_genera_relabel),
#                                          .default = nodes)) %>%
#   dplyr::mutate(nodes = stringr::str_split_i(nodes, ";", -1)) %>%
#   droplevels
# 
# temp_df3 <- temp_graph_edges %>%
#   dplyr::right_join(., spearman.test.filtered.df %>%
#                       dplyr::filter(grepl("low", group)) %>%
#                       dplyr::filter(!is.na(rho_threshold)) %>%
#                       dplyr::distinct(asv_id, simpleName, .keep_all = FALSE) %>%
#                       droplevels,
#                     by = join_by("node_start" == "asv_id", "node_end" == "simpleName")) %>%
#   dplyr::mutate(node_start = dplyr::case_when(grepl("ASV_", node_start) ~ dplyr::recode_factor(node_start, !!!usvi_genera_relabel),
#                                               .default = node_start)) %>%
#   dplyr::mutate(node_start = stringr::str_split_i(node_start, ";", -1)) %>%
#   droplevels
# 
# temp_list <- list(list(temp_df1, temp_df2) %>%
#                     setNames(., c("connections", "nodes")), 
#                   (list(temp_df3, temp_df4) %>%
#                      setNames(., c("connections", "nodes")))) %>%
#   setNames(., c("high_reef", "low_reef"))
# 
# 
# g1 <- print(
#   ggplot(data = temp_df1)
#   + theme_bw()
#   + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend,
#                      linewidth = abs_estimate, color = direction),
#                  arrow = arrow(length = unit(0.03, "npc"),
#                                angle = 30, ends = "last", type = "closed"),
#                  lineend = "butt", linejoin = "round",
#                  alpha = 0.7, show.legend = TRUE)
#   + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend),
#                  lineend = "butt", linejoin = "round", linetype = 2, color = "black",
#                  alpha = 0.7, show.legend = TRUE)
#   + geom_point(data = temp_df2, aes(x = V1, y = V2, fill = type), shape = 21, size = 3, alpha = 1)
#   + ggrepel::geom_label_repel(data = temp_df2, aes(x = V1, y = V2, label = nodes, fill = type),
#                               direction = "both", segment.color = NA, seed = 123,
#                               force = 1, force_pull = 1,
#                               point.padding = unit(0.01, "npc"),
#                               box.padding = unit(0.01, "npc"),
#                               colour = "black", fontface = "bold")
#   + scale_color_manual(values = c("#dec1aa", "#03c4da"),
#                        breaks = c("negative", "positive"),
#                        labels = c("negative", "positive"))
#   + scale_fill_manual(values = c("#c5b8dc", "#b9d2b1"),
#                       breaks = c("taxon", "metabolite"),
#                       labels = c("taxon", "metabolite"))
#   + scale_linewidth_continuous(name = "Strength of relationship")
#   + scale_x_continuous(expand = expansion(0.1,0.1))
#   + scale_y_continuous(expand = expansion(0.1,0.1))
#   + theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank())
#   + ggtitle("Taxa significantly more abundant in reef sites and strong correlations with metabolites")
# )
# 
# 
# 
# for(i in 1:length(temp_list)){
#   namevar <- eval(names(temp_list)[i]) %>%
#     gsub("_", " in ", .)
#   temp_df1 <- temp_list[[i]]$connections %>%
#     as.data.frame(.)
#   temp_df2 <- temp_list[[i]]$nodes %>%
#     as.data.frame(.)
#   
#   g <- print(
#     ggplot(data = temp_df1)
#     # ggplot(data = temp_list[[i]]$connections)
#     + theme_bw()
#     + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend,
#                        linewidth = abs_estimate, color = direction), 
#                    arrow = arrow(length = unit(0.03, "npc"), 
#                                  angle = 30, ends = "last", type = "closed"), 
#                    lineend = "butt", linejoin = "round",
#                    alpha = 0.7, show.legend = TRUE)
#     + geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend), 
#                    lineend = "butt", linejoin = "round", linetype = 2, color = "black",
#                    alpha = 0.7, show.legend = FALSE)
#     + geom_point(data = temp_df2,
#                  # + geom_point(data = temp_list[[i]]$node, 
#                  aes(x = V1, y = V2, fill = type), shape = 21, size = 3, alpha = 1)
#     + ggrepel::geom_label_repel(aes(x = V1, y = V2, label = nodes, fill = type),
#                                 # data = temp_list[[i]]$node, 
#                                 data = temp_df2, 
#                                 max.overlaps = 30,
#                                 direction = "both", segment.color = NA, seed = 123,
#                                 force = 1, force_pull = 1,
#                                 point.padding = unit(0.01, "npc"),
#                                 box.padding = unit(0.01, "npc"),
#                                 colour = "black", fontface = "bold")
#     + scale_color_manual(values = c("#dec1aa", "#03c4da"),
#                          breaks = c("negative", "positive"),
#                          labels = c("negative", "positive"), drop = FALSE)
#     + scale_fill_manual(values = c("#c5b8dc", "#b9d2b1"),
#                         breaks = c("taxon", "metabolite"),
#                         labels = c("taxon", "metabolite"), drop = FALSE)
#     + scale_linewidth_continuous(name = "Strength of relationship", limits = c(0, 0.9))
#     + scale_x_continuous(expand = expansion(0.1,0.1))
#     + scale_y_continuous(expand = expansion(0.1,0.1))
#     + theme(axis.text = element_blank(),
#             axis.ticks = element_blank(),
#             axis.title = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.grid.major = element_blank())
#     + ggtitle(paste0("Taxa significantly ", namevar, " sites"))
#   )
#   
#   assign(paste0("g", i), g, envir = .GlobalEnv, inherits = TRUE)
#   rm(g)
#   rm(namevar)
# }
# 
# gpatch2  <- ((g1 + theme(legend.position = "none")) | g2) + patchwork::plot_layout(guides = "collect")
# gpatch2 <- gpatch2 + patchwork::plot_annotation(title = "Network graphs of metabolite concentrations related to taxa abundances",
#                                                 tag_levels = "A")
# gpatch2
# if(!any(grepl("key_taxa_network", list.files(projectpath, pattern = "usvi_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_key_taxa_network-", Sys.Date(), ".png"),
#          gpatch2,
#          width = 18, height = 10, units = "in")
# }



#pantothenic acid emerged as significant coral-produced metabolites (Weber et al., 2022)