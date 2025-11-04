# 10_plot_corr.R

# plot correlations between ASVs and metabolites

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

if(file.exists(paste0(projectpath, "/", "usvi_sig_seqs_phylogeny.df", ".tsv"))){
  usvi_sig_seqs_phylogeny.df <- readr::read_delim(paste0(projectpath, "/", "usvi_sig_seqs_phylogeny.df", ".tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
} else {
  cli::cli_alert_warning("Please process the phylogeny of ASVs.")
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


# Read in correlation results ---------------------------------------------

#here were the alternatives evaluated given the 73 metabolome samples and 71 microbiome samples:
#A. average the metabolite profiles that are technical replicates, so there is 1 metabolome sample per microbiome sample
#B. keep the technical replicates, but duplicate the corresponding microbiome sample

#we want to use only option A (N = 71 samples)

to_import <- c("spearman.test.site.optA.list", #site-specific
               #"spearman.test.site.optB.list",  
               "spearman.test.optA.list", # #all sites together:
               #"spearman.test.optB.list", 
               "spearman.test.site.time.optA.list", #site- and time-specific
               #"spearman.test.site.time.optB.list",  
               NULL)

for(file in to_import){
  if(!exists(file, envir = .GlobalEnv)){
    namevar <- file
    if(file.exists(paste0(projectpath, "/", namevar, ".rds"))){
      cli::cli_alert_info("Importing this dataset: {namevar}")
      temp_df <- readr::read_rds(paste0(projectpath, "/", namevar,  ".rds"))
      assign(paste0(namevar), temp_df, envir = .GlobalEnv)
      rm(temp_df)
      to_import <- to_import[grep(namevar, to_import, value = FALSE, invert = TRUE)]
    } else {
      cli::cli_alert_warning("Please prepare this dataset: {namevar}")
    }
  }
}

# #all sites together:
# #option A:
# if(file.exists(paste0(projectpath, "/", "spearman.test.optA.list", ".rds"))){
#   spearman.test.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.optA.list", ".rds"))
# } 
# 
# #option B:
# if(file.exists(paste0(projectpath, "/", "spearman.test.optB.list", ".rds"))){
#   spearman.test.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.optB.list", ".rds"))
# }
# 
# 
# #site-specific
# #option A:
# if(file.exists(paste0(projectpath, "/", "spearman.test.site.optA.list", ".rds"))){
#   spearman.test.site.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.optA.list", ".rds"))
# }
# 
# #option B:
# if(file.exists(paste0(projectpath, "/", "spearman.test.site.optB.list", ".rds"))){
#   spearman.test.site.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.optB.list", ".rds"))
# }
# 
# #site- and time-specific
# 
# #option A: 71 samples
# if(file.exists(paste0(projectpath, "/", "spearman.test.site.time.optA.list", ".rds"))){
#   spearman.test.site.time.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.time.optA.list", ".rds"))
# } 
# #option B: 73 samples
# if(file.exists(paste0(projectpath, "/", "spearman.test.site.time.optB.list", ".rds"))){
#   spearman.test.site.time.optB.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.time.optB.list", ".rds"))
# }

#read in these 3 processed dfs: spearman.test.site.time.full.df, spearman.test.site.full.df, spearman.test.full.df
if(any(grepl("spearman_full", list.files(projectpath, pattern = "usvi_.*.RData")))){
    temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_spearman_full-.*.RData"))
    load(paste0(projectpath, "/", temp_file))
    rm(temp_file)
}


# #here is your categorization of metabolites by chemical category: 
dend_metab <- spearman.test.optA.list[["dend_metab"]]



# Filter for FDR ----------------------------------------------------------


#filter for only those passing FDR 1% in option A

if(any(grepl("spearman.sig.filtered.list", list.files(projectpath, pattern = "usvi_.*.rds")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_spearman.sig.filtered.list-*.rds"))
  temp_list <- readr::read_rds(paste0(projectpath, "/", temp_file))
  list2env(temp_list)
  rm(temp_list)
} else {
  for(dataset in c("spearman.test.site.time.full.df", 
                   "spearman.test.site.full.df", 
                   "spearman.test.full.df",
                   NULL)){
    temp_full.df <- get0(dataset, inherits = TRUE)
    namevar <- dataset %>%
      stringr::str_remove(pattern = ".full.df")
    
    if(!any(grepl("grouping", colnames(temp_full.df)))){
      temp_full.df <- temp_full.df %>%
        dplyr::mutate(grouping = "all.all")
    }
    temp_padj_cutoff <- temp_full.df %>%
      dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    
    temp_filtered.df <- temp_full.df %>%
      tidyr::drop_na(p_value) %>%
      dplyr::filter(grepl("optA", test_type)) %>% #after discussion, keep results from optA
      split(., f = .$grouping) %>%
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
      bind_rows(.)
    
    temp_filtered.df <- temp_filtered.df %>%
      dplyr::filter(!is.na(padj_01)) %>%
      dplyr::filter(grepl("optA", test_type)) %>%
      droplevels %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
      dplyr::mutate(across(c(asv_id, simpleName, test_type, sig, grouping), ~factor(.x))) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, test_type, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
    
    assign(paste0(namevar, ".filtered.df"), temp_filtered.df, envir = .GlobalEnv)
    rm(temp_full.df)
    rm(temp_padj_cutoff)
    rm(temp_filtered.df)
  }
  for(dataset in c("spearman.test.site.time.filtered.df", 
                   "spearman.test.site.filtered.df", 
                   "spearman.test.filtered.df",
                   NULL)){
    # if(!any(grepl("spearman.test.site.filtered.df", list.files(projectpath, pattern = "usvi_spearman.test.site.filtered.df.*.tsv")))){
    #   readr::write_delim(spearman.test.site.filtered.df, paste0(projectpath, "/", "spearman.test.site.filtered.df-", Sys.Date(), ".tsv"),
    #                      delim = "\t", col_names = TRUE)
    # }
    if(!any(grepl(dataset, list.files(projectpath, pattern = "usvi_spearman.test.*.tsv")))){
      temp_filtered.df <- get0(dataset, inherits = TRUE)
      readr::write_delim(temp_filtered.df, paste0(projectpath, "/", "usvi_", dataset, "-", Sys.Date(), ".tsv"),
                         delim = "\t", col_names = TRUE)
      rm(temp_filtered.df)
    }
    
  }
}

spearman.test.site.filtered.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(num_sda = length(test_type))
spearman.test.site.time.filtered.df %>%
  dplyr::group_by(grouping) %>%
  dplyr::summarise(num_sda = length(test_type))


# Prepare to plot ---------------------------------------------------------

if(any(grepl("spearman.sig.filtered.list", list.files(projectpath, pattern = "usvi_spearman.sig.filtered.list.*.rds")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_spearman.sig.filtered.list-.*.rds"))
  temp_list <- readr::read_rds(paste0(projectpath, "/", temp_file))
  spearman.sig.filtered.list <- temp_list[grep("site", names(temp_list), value = FALSE)] %>%
    bind_rows(., .id = NULL) %>%
    dplyr::ungroup(.) %>%
    # dplyr::filter(grepl("optA", test_type)) %>%
    dplyr::select(asv_id, simpleName, grouping, estimate, padj_bh_05, padj_01, filtered_estimate, site, sampling_time) %>%
    droplevels %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::rowwise(.) %>%
          tidyr::drop_na(filtered_estimate) %>%
          dplyr::filter(!is.na(padj_01) | !is.na(padj_bh_05)) %>%
          droplevels %>%
          dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
          droplevels %>%
          dplyr::arrange(asv_id) %>%
          dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
          # dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
          droplevels)
  
} else {
  spearman.sig.filtered.list <- spearman.test.site.filtered.df %>%
    bind_rows(., spearman.test.site.time.filtered.df) %>%
    dplyr::filter(grepl("optA", test_type)) %>%
    dplyr::select(asv_id, simpleName, grouping, estimate, padj_bh_05, padj_01, filtered_estimate, site, sampling_time) %>%
    droplevels %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          dplyr::rowwise(.) %>%
          tidyr::drop_na(filtered_estimate) %>%
          dplyr::filter(!is.na(padj_01) | !is.na(padj_bh_05)) %>%
          droplevels %>%
          dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_sig_seqs_phylogeny.df[["asv_id"]]))) %>%
          droplevels %>%
          dplyr::arrange(asv_id) %>%
          dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
          # dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
          droplevels)
  
  readr::write_rds(spearman.sig.filtered.list, paste0(projectpath, "/", "usvi_spearman.sig.filtered.list-", Sys.Date(), ".rds"), compress = "gz")
}

rm(list = apropos(paste0("^g8_asvs_", ".*_.*_results$"), mode = "list"))

for(i in c(1)){
# for(i in seq_len(length(spearman.sig.filtered.list))){
  namevar <- names(spearman.sig.filtered.list)[i]
  title_plot <- stringr::str_split_i(namevar, "\\.", 1) %>%
    recode(., !!!site_lookup)
  if(!is.na(stringr::str_split_i(namevar, "\\.", 2))){
    title_plot <- stringr::str_split_i(namevar, "\\.", 2) %>%
      recode(., !!!sampling_time_lookup) %>%
      paste0(title_plot, " at ", .)
  }
  title_plot_sig <- paste0("ASVs with significant correlations in ", title_plot)
  
  temp_spearman.df <- spearman.sig.filtered.list[[i]] %>%
    tibble::as_tibble(.)
  temp_dend <- temp_spearman.df %>%
    dplyr::select(asv_id, simpleName, filtered_estimate) %>%
    tidyr::pivot_wider(, id_cols = "asv_id", 
                       names_from = "simpleName",
                       values_fill = 0,
                       values_from = "filtered_estimate") %>%
    tibble::column_to_rownames(var = "asv_id") %>%
    t() %>%
    dist(t(.), method = "euclidean") %>%
    hclust(method = "ward.D2") %>%
    as.dendrogram
  
  temp_spearman.df2 <- temp_spearman.df %>%
    dplyr::summarise(label = length(estimate), .by = "asv_id") %>%
    dplyr::mutate(simpleName = "total correlations") %>%
    dplyr::arrange(desc(label)) %>%
    dplyr::mutate(simpleName = factor(simpleName)) %>%
    dplyr::bind_rows(., temp_spearman.df %>%
                       dplyr::summarise(label = length(estimate), .by = "simpleName") %>%
                       dplyr::mutate(asv_id = "total correlations") %>%
                       # dplyr::arrange(desc(label)) %>%
                       dplyr::mutate(asv_id = factor(asv_id)))
  
  temp_spearman.df2 <- temp_spearman.df2 %>%
    bind_rows(temp_spearman.df, .) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = c(usvi_sig_seqs_phylogeny.df[["asv_id"]], "total correlations"))) %>%
    dplyr::arrange(asv_id) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
    dplyr::mutate(asv_id = fct_relevel(asv_id, "total correlations", after = Inf)) %>%
    dplyr::mutate(simpleName = factor(simpleName, levels = unique(.[["simpleName"]]))) %>%
    dplyr::mutate(simpleName = fct_relevel(simpleName, "total correlations", after = Inf)) %>%
    droplevels
  
  
    fig_width <- max(10, (round(length(unique(temp_spearman.df2[["simpleName"]]))/10, digits = 0) + 8))
    fig_height <- max(8, (round(length(unique(temp_spearman.df2[["asv_id"]]))/100, digits = 0) + 4))
    temp_g8 <- (
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
    
    assign(paste0("g8_asvs_", namevar, "_results"), temp_g8, envir = .GlobalEnv)
   
    rm(temp_g8)
    rm(title_plot_sig)
  # gpatch <- lapply(apropos(paste0("^g7_asvs_", namevar, "_.*_results$"), mode = "list"),
  #                  get) %>%
  #   purrr::reduce(., `/`) + 
  #   patchwork::plot_layout(guides = "collect") &
  #   theme(legend.position="none")
  # assign(paste0("g7_asvs_", namevar, "_all_results"), gpatch, envir = .GlobalEnv)
  
  rm(title_plot)
  rm(namevar)
}

for(i in seq_len(length(spearman.sig.filtered.list))){
  namevar <- names(spearman.sig.filtered.list)[i]
  temp_g8 <- apropos(paste0("^g8_asvs_", namevar, "_results$"), mode = "any") %>%
    get0(.)
  fig_width <- max(10, (round(length(unique(temp_g8[["data"]][["simpleName"]]))/10, digits = 0) + 8))
  fig_height <- max(8, (round(length(unique(temp_g8[["data"]][["asv_id"]]))/100, digits = 0) + 4))
  # # if(!any(grepl(namevar, list.files(projectpath, pattern = "usvi_spearman_.*site_corr.*.png")))){
  # ggsave(paste0(projectpath, "/", "usvi_spearman_sig_strong_", namevar,"_corr-", Sys.Date(), ".png"),
  #        temp_g8,
  #        width = fig_width, height = fig_height,
  #        units = "in")
  # # }
  ggsave(paste0(projectpath, "/", "usvi_spearman_sig_strong_", namevar,"_corr-", Sys.Date(), ".svg"),
         temp_g8,
         width = fig_width, height = fig_height,
         units = "in")
}

# Graphnetwork ------------------------------------------------------------

#test it out:
{
  temp_spearman.df <- spearman.sig.filtered.list[[1]] %>%
    tibble::as_tibble(.)
  
  temp_spearman.tbl <- temp_spearman.df %>%
    dplyr::select(asv_id, simpleName, filtered_estimate) %>%
    ggsankey::make_long(., asv_id, simpleName, value = "filtered_estimate") %>%
    tidyr::drop_na(next_node) %>%
    dplyr::rename(from = "node",
                  to = "next_node",
                  estimate = "value") %>%
    dplyr::left_join(., temp_spearman.df %>%
                       dplyr::select(asv_id, simpleName, padj_01) %>%
                       droplevels,
                     by = join_by("from" == "asv_id", "to" == "simpleName")) %>%
    droplevels
  
  temp_spearman.graph.nodes <- temp_spearman.tbl %>%
    dplyr::select(from) %>%
    dplyr::rename(name = "from") %>%
    dplyr::mutate(omics_type = "ASV") %>%
    bind_rows(., (temp_spearman.tbl %>%
                    dplyr::select(to) %>%
                    dplyr::rename(name = "to") %>%
                    dplyr::mutate(omics_type = "metabolite"))) %>%
    dplyr::distinct(name, omics_type)
  
  spearman_network.nodes <- temp_spearman.graph.nodes %>%
    tibble::rowid_to_column(var = "Id") %>%
    dplyr::rename(Label = "name") %>%
    dplyr::relocate(Label)
  
  spearman_network.edges <- temp_spearman.tbl %>%
    dplyr::select(from, to, estimate) %>% 
    dplyr::inner_join(spearman_network.nodes, by = join_by("from" == "Label")) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(spearman_network.nodes, by = join_by("to" == "Label")) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::select(-contains("omics_type")) %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = abs(estimate)) %>%
    dplyr::select(Source, Target, Type, Weight, estimate) %>%
    droplevels
  
  # readr::write_delim(spearman_network.nodes, paste0(projectpath, "/", "spearman_network.nodes", ".csv"),
  #                    delim = ",", col_names = TRUE)
  # readr::write_delim(spearman_network.edges, paste0(projectpath, "/", "spearman_network.edges", ".csv"),
  #                    delim = ",", col_names = TRUE)
}
#it works, so iterate through the other lists and export the files for making networks.


for(i in seq_len(length(spearman.sig.filtered.list))){
  namevar <- names(spearman.sig.filtered.list)[i]
  temp_spearman.df <- spearman.sig.filtered.list[[i]] %>%
    tibble::as_tibble(.)
  temp_spearman.tbl <- temp_spearman.df %>%
    dplyr::select(asv_id, simpleName, filtered_estimate) %>%
    ggsankey::make_long(., asv_id, simpleName, value = "filtered_estimate") %>%
    tidyr::drop_na(next_node) %>%
    dplyr::rename(from = "node",
                  to = "next_node",
                  estimate = "value") %>%
    dplyr::left_join(., temp_spearman.df %>%
                       dplyr::select(asv_id, simpleName, padj_01) %>%
                       droplevels,
                     by = join_by("from" == "asv_id", "to" == "simpleName")) %>%
    droplevels
  
  
  temp_spearman.graph.nodes <- temp_spearman.tbl %>%
    dplyr::select(from) %>%
    dplyr::rename(name = "from") %>%
    dplyr::mutate(omics_type = "ASV") %>%
    bind_rows(., (temp_spearman.tbl %>%
                    dplyr::select(to) %>%
                    dplyr::rename(name = "to") %>%
                    dplyr::mutate(omics_type = "metabolite"))) %>%
    dplyr::distinct(name, omics_type) %>%
    dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"),
                     by = join_by("name" == "asv_id"))
  
  temp_spearman_network.nodes <- temp_spearman.graph.nodes %>%
    tibble::rowid_to_column(var = "Id") %>%
    dplyr::rename(Label = "name") %>%
    dplyr::relocate(Label) %>%
    dplyr::mutate(taxonomy = dplyr::case_when((is.na(taxonomy) & omics_type == "metabolite") ~ Label, 
                                              (is.na(taxonomy) & omics_type == "ASV") ~ Label,
                                              .default = taxonomy))
  
  temp_spearman_network.edges <- temp_spearman.tbl %>%
    dplyr::select(from, to, estimate, padj_01) %>% 
    dplyr::inner_join(temp_spearman_network.nodes, by = join_by("from" == "Label")) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(temp_spearman_network.nodes, by = join_by("to" == "Label")) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::select(-contains("omics_type")) %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = abs(estimate)) %>%
    dplyr::distinct(Source, Target, Type, Weight, estimate, padj_01) %>%
    droplevels
  
  readr::write_delim(temp_spearman_network.nodes, paste0(projectpath, "/", "spearman_network.", namevar, ".nodes", ".csv"),
                     delim = ",", col_names = TRUE)
  readr::write_delim(temp_spearman_network.edges, paste0(projectpath, "/", "spearman_network.", namevar, ".edges", ".csv"),
                     delim = ",", col_names = TRUE)
  rm(temp_spearman_network.nodes)
  rm(temp_spearman_network.edges)
}


#how to decide which nodes and edges to plot?

temp_spearman.df <- spearman.sig.filtered.list[["LB_seagrass"]] %>%
  tibble::as_tibble(.) %>%
  dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "label")) %>%
  dplyr::mutate(taxonomy = stringr::str_split_i(label, ": ", 2))


#what if we agglomerated by taxonomy and direction of correlation?
# temp_spearman_filtered.df <- temp_spearman.df %>%
#   dplyr::mutate(direction = round(estimate, digits = 0)) %>%
#   dplyr::mutate(direction = dplyr::case_when(direction < 0 ~ "down", direction > 0 ~ "up", .default = NA)) %>%
#   dplyr::mutate(direction = factor(direction)) %>%
#   dplyr::group_by(taxonomy, simpleName, direction) %>%
#   dplyr::summarise(num_obs = length(direction)) %>%
#   dplyr::mutate(estimate_direction = dplyr::case_when(direction == "up" ~ 0.5, direction == "down" ~ -0.5, .default = NA)) %>%
#   dplyr::arrange(taxonomy, simpleName, desc(num_obs))
# 
# temp_spearman_filtered.df %>%
#   dplyr::group_by(simpleName) %>%
#   dplyr::reframe(distribution = quantile(num_obs, probs = seq(1, 0, -0.25), na.rm = TRUE, names = TRUE)) %>%
#   dplyr::distinct(.) %>%
#   dplyr::arrange(desc(distribution)) %>%
#   dplyr::mutate(simpleName = factor(simpleName, levels = unique(.[["simpleName"]]))) %>%
#   dplyr::arrange(simpleName, desc(distribution))
# 
# readr::write_delim(temp_spearman_filtered.df, paste0(projectpath, "/", "temp_spearman_filtered.df", ".tsv"),
#                    delim = "\t", col_names = TRUE)
# 
# 
#simplify it for Gephi

{
  temp_spearman_filtered.graph.nodes <- temp_spearman_filtered.df %>%
    dplyr::ungroup(.) %>%
    dplyr::select(taxonomy, direction, num_obs) %>%
    dplyr::rename(name = "taxonomy") %>%
    dplyr::mutate(omics_type = "microbe") %>%
    dplyr::arrange(name, omics_type, direction) %>%
    bind_rows(., (temp_spearman_filtered.df %>%
                    dplyr::ungroup(.) %>%
                    dplyr::select(simpleName, direction, num_obs) %>%
                    dplyr::rename(name = "simpleName") %>%
                    dplyr::mutate(omics_type = "metabolite") %>%
                    dplyr::arrange(name, omics_type, direction))) %>%
    dplyr::distinct(name, direction, num_obs, omics_type)
  
  temp_spearman_filtered.graph.nodes <- temp_spearman_filtered.graph.nodes %>%
    tibble::rowid_to_column(var = "Id") %>%
    dplyr::rename(Label = "name") %>%
    dplyr::relocate(Label)
  
  temp_spearman_filtered.graph.edges <- temp_spearman_filtered.df %>%
    dplyr::select(taxonomy, simpleName, direction, num_obs, estimate_direction) %>% 
    dplyr::inner_join(temp_spearman_filtered.graph.nodes, by = join_by("taxonomy" == "Label", direction, num_obs)) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(temp_spearman_filtered.graph.nodes, by = join_by("simpleName" == "Label", direction, num_obs)) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::select(-contains("omics_type")) %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = num_obs) %>%
    dplyr::distinct(Source, Target, Type, Weight, estimate_direction) %>%
    droplevels
  
  temp_spearman.graph <- temp_spearman_filtered.graph.edges %>%
    dplyr::select(taxonomy, simpleName, estimate_direction) %>% dplyr::distinct(., .keep_all = TRUE) %>%
    igraph::graph_from_data_frame(., directed = TRUE, vertices = NULL)
  plot(temp_spearman.graph)

  # readr::write_delim(temp_spearman_filtered.graph.nodes, paste0(projectpath, "/", "temp_spearman_filtered.graph.nodes", ".tsv"),
  #                    delim = "\t", col_names = TRUE)
  # 
  # readr::write_delim(temp_spearman_filtered.graph.edges, paste0(projectpath, "/", "temp_spearman_filtered.graph.edges", ".tsv"),
  #                    delim = "\t", col_names = TRUE)

  }

#generalize the loop to export the results:

for(i in seq_len(length(spearman.sig.filtered.list))){
  namevar <- names(spearman.sig.filtered.list)[i]
  temp_spearman.df <- spearman.sig.filtered.list[[i]] %>%
    tibble::as_tibble(.) %>%
    dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "label")) %>%
    dplyr::mutate(taxonomy = stringr::str_split_i(label, ": ", 2))
  
  temp_spearman_filtered.df <- temp_spearman.df %>%
    dplyr::mutate(direction = round(estimate, digits = 0)) %>%
    dplyr::mutate(direction = dplyr::case_when(direction < 0 ~ "down", direction > 0 ~ "up", .default = NA)) %>%
    dplyr::mutate(direction = factor(direction)) %>%
    dplyr::group_by(taxonomy, simpleName, direction) %>%
    dplyr::summarise(num_obs = length(direction)) %>%
    dplyr::mutate(estimate_direction = dplyr::case_when(direction == "up" ~ 0.5, direction == "down" ~ -0.5, .default = NA)) %>%
    dplyr::arrange(taxonomy, simpleName, desc(num_obs))
  temp_spearman_filtered.graph.nodes <- temp_spearman_filtered.df %>%
    dplyr::ungroup(.) %>%
    dplyr::select(taxonomy, direction, num_obs) %>%
    dplyr::rename(name = "taxonomy") %>%
    dplyr::mutate(omics_type = "microbe") %>%
    dplyr::arrange(name, omics_type, direction) %>%
    bind_rows(., (temp_spearman_filtered.df %>%
                    dplyr::ungroup(.) %>%
                    dplyr::select(simpleName, direction, num_obs) %>%
                    dplyr::rename(name = "simpleName") %>%
                    dplyr::mutate(omics_type = "metabolite") %>%
                    dplyr::arrange(name, omics_type, direction))) %>%
    dplyr::distinct(name, direction, num_obs, omics_type)
  
  temp_spearman_filtered.graph.nodes <- temp_spearman_filtered.graph.nodes %>%
    tibble::rowid_to_column(var = "Id") %>%
    dplyr::rename(Label = "name") %>%
    dplyr::relocate(Label)
  
  temp_spearman_filtered.graph.edges <- temp_spearman_filtered.df %>%
    dplyr::select(taxonomy, simpleName, direction, num_obs, estimate_direction) %>% 
    dplyr::inner_join(temp_spearman_filtered.graph.nodes, by = join_by("taxonomy" == "Label", direction)) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(temp_spearman_filtered.graph.nodes, by = join_by("simpleName" == "Label", direction)) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::select(-contains("omics_type")) %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = num_obs) %>%
    dplyr::distinct(Source, Target, Type, Weight, estimate_direction) %>%
    droplevels
  
  readr::write_delim(temp_spearman_filtered.graph.nodes, paste0(projectpath, "/", "spearman_taxonomy_network.", namevar, ".nodes", ".tsv"),
                     delim = "\t", col_names = TRUE)
  
  readr::write_delim(temp_spearman_filtered.graph.edges, paste0(projectpath, "/", "spearman_taxonomy_network.", namevar, ".edges", ".tsv"),
                     delim = "\t", col_names = TRUE)
  rm(temp_spearman_filtered.graph.nodes)
  rm(temp_spearman_filtered.graph.edges)
}



#what if we retained only those ASVs with a lot of correlations?
{
  # temp_spearman_rankabund <- temp_spearman.df %>%
  #   dplyr::group_by(asv_id) %>%
  #   # dplyr::group_by(simpleName) %>%
  #   dplyr::summarise(num_obs = length(estimate)) %>%
  #   dplyr::arrange(desc(num_obs)) %>%
  #   dplyr::reframe(distribution = quantile(num_obs, probs = seq(1, 0, -0.25), na.rm = TRUE, names = TRUE)) %>%
  #   tibble::deframe(.)
  # 
  # 
  # temp_spearman_filtered.df <- temp_spearman.df %>%
  #   dplyr::right_join(., temp_spearman.df %>%
  #                       dplyr::group_by(asv_id) %>%
  #                       # dplyr::group_by(simpleName) %>%
  #                       dplyr::summarise(num_obs = length(estimate)) %>% dplyr::arrange(desc(num_obs)) %>% dplyr::filter(num_obs >= temp_spearman_rankabund[2]),
  #                     by = join_by(asv_id), relationship = "many-to-many", multiple = "all")
  # 
  # temp_spearman_filtered.df %>%
  #   dplyr::group_by(simpleName) %>%
  #   dplyr::summarise(num_obs = length(estimate)) %>%
  #   dplyr::arrange(desc(num_obs))
}



# temp_spearman.graph.nodes <- temp_spearman.tbl %>%
#   dplyr::select(from) %>%
#   dplyr::rename(name = "from") %>%
#   dplyr::mutate(omics_type = "ASV") %>%
#   bind_rows(., (temp_spearman.tbl %>%
#                   dplyr::select(to) %>%
#                   dplyr::rename(name = "to") %>%
#                   dplyr::mutate(omics_type = "metabolite"))) %>%
#   dplyr::distinct(name, omics_type)
# temp_spearman.graph <- temp_spearman.tbl %>%
#   dplyr::select(from, to, estimate) %>% dplyr::distinct(., .keep_all = TRUE) %>% 
#   igraph::graph_from_data_frame(., directed = TRUE, vertices = temp_spearman.graph.nodes)
# 
# plot(temp_spearman.graph)

# temp_spearman_graph.coords <- igraph::layout_nicely(temp_spearman.graph, dim = 2) %>%
#     tibble::as_tibble(.) %>%
#     dplyr::mutate(nodes = igraph::vertex_attr(temp_spearman.graph, "name")) %>%
#     dplyr::relocate(nodes) %>%
#     dplyr::mutate(type = dplyr::case_when(grepl("ASV_", nodes) ~ "taxon",
#                                           .default = "metabolite")) %>%
#   droplevels
# igraph::write_graph(temp_spearman.graph, paste0(projectpath, "/", "temp_spearman.graph", ".txt"),
#                     format = "graphml")
# 
# # igraph::edge_connectivity(temp_spearman.graph, source = NULL, target= NULL)
# temp_spearman.graph2 <- temp_spearman.tbl %>%
#   dplyr::select(from, to, estimate) %>% dplyr::distinct(., .keep_all = TRUE) %>% 
#   igraph::graph_from_data_frame(., directed = FALSE, vertices = temp_spearman.graph.vertex)
# 
# temp_spearman.graph2_blocks <- igraph::cohesive_blocks(temp_spearman.graph2, labels = TRUE)
# temp_spearman_blocks.graph2 <- igraph::graphs_from_cohesive_blocks(temp_spearman.graph2_blocks, temp_spearman.graph2)
# 
# igraph::max_cardinality(temp_spearman.graph2)
# 
# plot(temp_spearman_blocks.graph2[[2]])





