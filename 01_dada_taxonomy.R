## 01_dada_taxonomy.R

# import the amplicon fastqs for processing through DADA2

# Load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
library(BiocParallel)
library(dada2)
library(cli)
library(furrr)
library(progressr)

# Plan for resource allocation --------------------------------------------
if(grepl("arch64", Sys.getenv("R_PLATFORM"))){
  print("Using parallelly...")
  nthreads <- parallelly::availableCores() - 1
  future::plan(multisession, workers = nthreads)
  options(future.globals.maxSize = 10e9)
} else {
  print("Using data.table")
  nthreads <- data.table::getDTthreads()
  future::plan(sequential)
  options(future.globals.maxSize = 10e9)
}


# nthreads <- data.table::getDTthreads()
# cluster <- multidplyr::new_cluster(n = nthreads)
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

# Read in metadata --------------------------------------------------------



sample_metadata <- readr::read_delim(paste0(projectpath, "/", "usvi_metadata.txt"),
                                     delim = "\t",
                                     quote = "",
                                     col_names = TRUE,
                                     show_col_types = FALSE,
                                     num_threads = nthreads)



# find existing processed files -------------------------------------------

if(!exists("usvi_dada_merged", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_dada_merged", ".rds"))){
    cli::cli_alert_info("Reading in DADA2 processed dataset...")
    usvi_dada_merged <- readr::read_rds(paste0(projectpath, "/", "usvi_dada_merged", ".rds"))  
  }
}



if(file.exists(paste0(projectpath, "/", "usvi_dada_nochim", ".rds"))){
  cli::cli_alert_info("Reading in DADA2 non-chimeric dataset...")
  usvi_dada_nochim <- readr::read_rds(paste0(projectpath, "/", "usvi_dada_nochim", ".rds"))
  list2env(usvi_dada_nochim, envir = .GlobalEnv)
  rm(usvi_dada_nochim)
}

if(!exists("usvi_qc_summary_df", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_qc_summary_df", ".tsv"))){
    cli::cli_alert_info("Reading in DADA2 processing summary...")
    usvi_qc_summary_df <- readr::read_delim(paste0(projectpath, "/", "usvi_qc_summary_df", ".tsv"),
                                         show_col_types = FALSE, delim  = "\t", col_names = TRUE, num_threads = nthreads)
  } else {
    # if(!exists("usvi_dada_res", envir = .GlobalEnv)){
    #   if(any(grepl("usvi_dada_res", list.files(path = projectpath, pattern = ".rds", recursive = FALSE, include.dirs = TRUE, full.names = TRUE)))){
    #     cli::cli_alert_info("Loading DADA2-processed files...")
    #     varname <- sort(grep("usvi_dada_res", list.files(path = projectpath, pattern = ".rds", recursive = FALSE, include.dirs = TRUE, full.names = TRUE), value = TRUE), decreasing = TRUE)[1]
    #     usvi_dada_res <- readr::read_rds(varname)
    #     list2env(usvi_dada_res, envir = .GlobalEnv)
    #   } else {
    #     print("This output does not exist")
    #   }
    # }
    cli::cli_alert_warning("Please process the QC summary from DADA2 filterAndTrim.")
  }  
}


# Assign taxonomy ---------------------------------------------------------

if(!exists("usvi_nochim_asvs.taxa.df", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_nochim_asvs.taxa.tsv", ".gz"))){
    usvi_nochim_asvs.taxa <- readr::read_delim(paste0(projectpath, "/", "usvi_nochim_asvs.taxa.tsv", ".gz"), 
                                                delim = "\t", show_col_types = FALSE, 
                                                col_names = TRUE, num_threads = nthreads)    
  } else {
    cli::cli_alert_info("Assigning taxonomy to ASVs...")
    usvi_nochim_asvs.taxa <- dada2::assignTaxonomy(usvi_nochim_asvs, 
                                                    "/proj/omics/bioinfo/databases/dada2/silva_nr99_v138.1_wSpecies_train_set.fa.gz", 
                                                    minBoot = 50,
                                                    tryRC = TRUE,
                                                    multithread = nthreads)
    usvi_nochim_asvs.taxa <- usvi_nochim_asvs.taxa %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., var = "sequence")
    
  usvi_nochim_asvs.taxa.df <- usvi_nochim_asvs.taxa %>%
    left_join(., usvi_asvs_key) %>%
    dplyr::relocate(asv_id, .before = "sequence") %>%
    dplyr::relocate(sequence, .after = "Species")
  readr::write_delim(usvi_nochim_asvs.taxa.df, paste0(projectpath, "/", "usvi_nochim_asvs.taxa.tsv", ".gz"),
                     delim  = "\t", col_names = TRUE, num_threads = nthreads)
  
  }
}

if(!file.exists(paste0(projectpath, "/", "usvi_nochim_asvs.df", ".tsv", ".gz"))){
  
  usvi_nochim_asvs.df <- usvi_nochim_asvs %>%
    # t() %>%
    # as.data.frame() %>%
    # tibble::rownames_to_column(var = "sequence") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name") %>%
    dplyr::mutate(sample_ID = dplyr::case_when(grepl("_", sample_name) ~ sample_name,
                                               # grepl("bowhead", sample_name) ~ stringr::str_to_title(sample_name),
                                               .default = paste0("Metab_", sample_name))) %>%
    dplyr::relocate(sample_ID) %>%
    dplyr::select(-sample_name) %>%
    data.table::transpose(., keep.names = "sequence", make.names = "sample_ID") %>%
    right_join(., usvi_nochim_asvs.taxa.df %>%
                 dplyr::select(asv_id, sequence),
               by = "sequence") %>%
    droplevels %>%
    dplyr::relocate(asv_id, .before = "sequence") %>%
    dplyr::select(-sequence) %>%
    tidyr::pivot_longer(., cols = -c("asv_id"),
                        names_to = "sample_ID",
                        values_to = "counts")
  
  readr::write_delim(usvi_nochim_asvs.df,
                     file = paste0(projectpath, "/", "usvi_nochim_asvs.df", ".tsv", ".gz"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
}


#assign species if available to ASVs
# if(!exists("usvi_nochim_asvs.species.df", envir = .GlobalEnv)){
#   if(file.exists(paste0(projectpath, "/", "usvi_nochim_asvs.species.tsv", ".gz"))){
#     usvi_nochim_asvs.species <- readr::read_delim(paste0(projectpath, "/", "usvi_nochim_asvs.species.tsv", ".gz"),
#                                                    delim = "\t", show_col_types = FALSE, 
#                                                    col_names = TRUE, num_threads = nthreads)    
#   } else {
#     cli::cli_alert_info("Assigning species to ASVs...")
#     usvi_nochim_asvs.species <- dada2::assignSpecies(usvi_nochim_asvs, 
#                                                       "/proj/omics/bioinfo/databases/dada2/silva_species_assignment_v138.1.fa.gz", 
#                                                       tryRC = TRUE,
#                                                       n = 2000)
#     usvi_nochim_asvs.species <- usvi_nochim_asvs.species %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column(., var = "sequence")
#     
#   usvi_nochim_asvs.species.df <- usvi_nochim_asvs.species %>%
#     left_join(., usvi_asvs_key, by = "sequence") %>%
#     dplyr::relocate(asv_id, .before = "sequence") %>%
#     dplyr::relocate(sequence, .after = "Species")
#   readr::write_delim(usvi_nochim_asvs.species.df, paste0(projectpath, "/", "usvi_nochim_asvs.species.tsv", ".gz"), 
#                      delim  = "\t", col_names = TRUE, num_threads = nthreads)
#   }
# }



# Remove any non-bact and arch sequences ----------------------------------


#drop any chlorphyll, mitochondria, eukaryotic sequences
if(!exists("usvi_prok_asvs.taxa.df", envir = .GlobalEnv) & file.exists(paste0(projectpath, "/", "usvi_prok_asvs.taxa.tsv", ".gz"))){
  usvi_prok_asvs.taxa.df <- readr::read_delim(paste0(projectpath, "/", "usvi_prok_asvs.taxa.tsv", ".gz"),
                                               delim = "\t", col_names = TRUE, num_threads = nthreads, show_col_types = FALSE)
} else if(!file.exists(paste0(projectpath, "/", "usvi_prok_asvs.taxa.tsv", ".gz")) & exists("usvi_nochim_asvs.taxa.df", envir = .GlobalEnv)){
  drop <- c("Chloroplast", "mitochondria", "eukary")
  usvi_prok_asvs.taxa.df <- usvi_nochim_asvs.taxa.df %>%
    dplyr::filter(if_all(c("Kingdom":"Species"), ~!grepl(paste0(drop, collapse = "|"), .x, ignore.case = TRUE))) %>%
    droplevels %>%
    dplyr::mutate(across(c("Kingdom":"Species"), ~ifelse(grepl("Cyanobacteriia", .x),
                                                         "Cyanobacteria", .x))) %>%
    dplyr::mutate(across(c("Kingdom":"Species"), ~ifelse(grepl("iia", .x, ignore.case = FALSE),
                                                         gsub("iia", "ia", .),
                                                         .x)))
  readr::write_delim(usvi_prok_asvs.taxa.df,
                     file = paste0(projectpath, "/", "usvi_prok_asvs.taxa.tsv", ".gz"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
}

if(!file.exists(paste0(projectpath, "/", "usvi_prok_asvs.df", ".tsv", ".gz"))){
  usvi_prok_asvs.df <- usvi_nochim_asvs %>%
    # t() %>%
    # as.data.frame() %>%
    # tibble::rownames_to_column(var = "sequence") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name") %>%
    dplyr::mutate(sample_ID = dplyr::case_when(grepl("_", sample_name) ~ sample_name,
                                               # grepl("bowhead", sample_name) ~ stringr::str_to_title(sample_name),
                                               .default = paste0("Metab_", sample_name))) %>%
    dplyr::relocate(sample_ID) %>%
    dplyr::select(-sample_name) %>%
    data.table::transpose(., keep.names = "sequence", make.names = "sample_ID") %>%
    dplyr::right_join(., usvi_prok_asvs.taxa.df %>%
                        dplyr::select(asv_id, sequence),
                      by = "sequence") %>%
    dplyr::relocate(asv_id, .before = "sequence") %>%
    dplyr::select(-sequence) %>%
    tidyr::pivot_longer(., cols = -c("asv_id"),
                        names_to = "sample_ID",
                        values_to = "counts")
  
  readr::write_delim(usvi_prok_asvs.df,
                     file = paste0(projectpath, "/", "usvi_prok_asvs.df", ".tsv", ".gz"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
}


#add statistics to the summary file
usvi_prok_asvs_summary <- usvi_prok_asvs.df %>%
  tidyr::pivot_longer(., cols = !c(sequence, asv_id),
                      names_to = "sample_ID",
                      values_to = "counts") %>%
  dplyr::select(-sequence)

if(!any(grepl("reads.final", colnames(usvi_qc_summary_df)))){
  # cli::cli_alert_info("Adding reads post-taxonomy checking")
  #reads.final: the total reads after qc, chimera-checking, and taxonomy filtering
  usvi_qc_summary_df <- usvi_qc_summary_df %>%
    dplyr::left_join(., (usvi_prok_asvs_summary %>%
                    dplyr::group_by(sample_ID) %>%
                    dplyr::summarise(reads.final = sum(counts))),
              by = join_by(sample_ID))
  readr::write_delim(usvi_qc_summary_df, paste0(projectpath, "/", "usvi_qc_summary_df", ".tsv"),
                    delim  = "\t", col_names = TRUE, num_threads = nthreads)
}

