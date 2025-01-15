# 05c_deseq.R

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

no_progress <- function() {} #your placeholder function for progress reporting

idx_samples <- function(x) grep("samples", names(x), value = TRUE)
coalesce2 <- function(x, y, sep = ".") ifelse(x == y, coalesce(x, y, sep = sep), paste0(x, "_vs_", y))

set.alpha <- 0.05
set.seed(48105)

F_deseq_frac <- function(counts, mdata, mdata_filter, deseq_design, prefix,
                         p = no_progress){
  #here is a generalized function to conduct deseq2
  
  #counts : your counts table
  
  #mdata: your metadata
  
  #mdata_filter: optional selector of your data
  #provide the filter a list of named vectors, e.g.
  # metadata_filter <- list(`Fraction` = c("12L", "13H"))
  # metadata_filter <- list(`Type` = c("SIP"),
  #                         `Temperature` = c("30; 55"))
  
  #due to DESeq's design, all samples included via this mdata_filter will be used
  #to calculate a baseline expression for the genes
  #if you don't want that, you will need to modify your filter, or account for it post-hoc
  
  #deseq_design: the formula
  
  #`prefix`: name prefix for the output list
  
  if(!missing(prefix)){
    prefix <- rlang::parse_expr(prefix)
  } else {
    prefix <- "default"
  }
  
  if(exists(paste0(prefix, "_deseq_list"), envir = .GlobalEnv)){
    cli::cli_alert_warning("A DESeq list object already exists with this name: {prefix}_deseq_list", wrap = TRUE)
    cli::cli_abort("Ending DESeq attempt.")
  }
  
  #d
  deseq_design <- rlang::parse_expr(deseq_design)
  
  #m
  if(!(is.null(mdata_filter))){
    metadata <- map2(mdata_filter, names(mdata_filter),
                     ~mdata %>%
                       dplyr::filter(
                         grepl(paste0(.x, collapse = "|"), .data[[.y]], ignore.case = TRUE)
                       )) %>%
      purrr::reduce(inner_join) %>%
      bind_rows %>%
      dplyr::distinct(sample_name, .keep_all = TRUE) %>%
      droplevels
  } else {
    metadata <- mdata %>%
      dplyr::mutate(across(everything(.), ~factor(.x, levels = unique(.x)))) %>%
      droplevels
  }
  
  # metadata <- metadata %>%
  #     dplyr::mutate(across(everything(), ~factor(.x, levels = unique(.x)))) %>%
  #     droplevels
  
  #check the design formula against the available metadata to make sure we have enough samples
  mod_mat <- model.matrix(eval(deseq_design),
                          data = metadata)
  
  if(length(which(colSums(mod_mat) == 0)) > 0){ #none of these should be 0--if so, you need to fix your design/metadata coding
    idx <- names(colSums(mod_mat))[which(colSums(mod_mat) == 0)]
    idx <- cli::cli_vec(idx)
    cli::cli_alert_warning("One or more groups in the comparison does not have enough replicates:", wrap = TRUE)
    cli::cli_bullets_raw(idx)
    cli::cli_abort("Ending DESeq attempt.")
  }
  
  #data
  colnames(counts)[1] <- "gene_id"
  # if(!("sample_name" %in% colnames(counts))){
  # if(!("sample" %in% colnames(counts))){
  if(!any(grepl("sample", colnames(deseq_asvs.df)))){
    cli::cli_alert_info("Formatting the counts data to be a long table.")
    counts <- counts %>%
      tidyr::pivot_longer(., 
                          cols = !c(contains("gene"), contains("asv")),
                          names_to = "sample_name",
                          values_to = "counts")
    assign("counts", counts, envir = .GlobalEnv, inherits = FALSE)
  }
  matrix <- counts %>%
    dplyr::mutate(across(!contains("counts"), factor)) %>%
    dplyr::filter(sample_name %in% as.character(metadata[["sample_name"]])) %>%
    droplevels %>%
    tidyr::drop_na() %>%
    tidyr::pivot_wider(., id_cols = "gene_id",
                names_from = "sample_name",
                values_fill = 0,
                values_from = "counts") %>%
    tibble::column_to_rownames(., var = "gene_id") %>%
    dplyr::select(metadata[["sample_name"]]) %>%
    dplyr::slice(which(rowSums(.) > 0))
  
  dds <- DESeqDataSetFromMatrix(countData = matrix,
                                tidy = FALSE,
                                colData = metadata,
                                design = eval(deseq_design)) %>%
    DESeq(., parallel = TRUE,
          BPPARAM = bpparam_multi)
  
  # resultsNames(dds)
  # 
  # mod_mat <- model.matrix(design(dds), colData(dds))
  # 
  temp_list <- list(metadata, matrix, dds, mod_mat) %>%
    setNames(., c("metadata", "matrix", "dds", "mod_mat"))
  assign(paste0(prefix, "_deseq_list"), temp_list, envir = .GlobalEnv)
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

# long metabolomics dataset from Brianna confirms that CINAR_BC_81A is the BC sample associated with Tektite Metab_319
# and CINAR_BC_81B is the BC sample associated with LB_seagrass Metab_219


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
  variable = variable_labels_lookup,
  site = site_lookup,
  sampling_day = sampling_day_lookup,
  sampling_time = sampling_time_lookup,
  asv_id = usvi_genera_relabel,
  .multi_line = TRUE,
  .default = label_wrap_gen(25, multi_line = TRUE)
  # .default = label_value
)


# Prepare models for de ---------------------------------------------------

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



# DESeq2 on genera-level table --------------------------------------------


#don't need to do genera-level DESeq2
##let's try consolidating the ASV table to species- or genus-level first
{
  
  # usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
  #   # dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
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
  #   droplevels
  # 
  # usvi_sw_genus.tbl <- usvi_sw_genus.tbl %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   # apply(., 2, relabund) %>%
  #   # as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   as.data.frame(.)
  # usvi_sw_genus.mat <- usvi_sw_genus.tbl %>%
  #   t(.)
  # 
  # 
  # usvi_sw_genus.df <- usvi_sw_genus.mat %>%
  #   as.data.frame(.) %>%
  #   tibble::rownames_to_column(var = "sample_id") %>%
  #   tidyr::pivot_longer(., cols = !c(sample_id),
  #                       names_to = "asv_id",
  #                       values_to = "counts")
  #
  # 
  # deseq_design <- "~0 + site:sampling_time"
  # deseq_design <- rlang::parse_expr(deseq_design)
  # mod_mat <- model.matrix(eval(deseq_design),
  #                         data = usvi_selected_metadata)
  # 
  # if(length(grep("0", colSums(mod_mat))) > 0){ #none of these should be 0--if so, you need to fix your design/metadata coding
  #   idx <- names(colSums(mod_mat))[grep("0", colSums(mod_mat))]
  #   idx <- cli::cli_vec(idx)
  #   cli::cli_alert_warning("One or more groups in the comparison does not have enough replicates:", wrap = TRUE)
  #   cli::cli_bullets_raw(idx)
  #   cli::cli_abort("Ending DESeq attempt.")
  # }
  # 
  # 
  # deseq_df <- usvi_sw_genus.tbl %>%
  #   as.data.frame(.) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   tidyr::pivot_longer(., cols = !c(asv_id),
  #                       names_to = "sample_name",
  #                       values_to = "counts")
  # 
  # deseq_metadata <- usvi_selected_metadata %>%
  #   dplyr::rename(sample_name = "sample_id")
  # 
  # 
  # F_deseq_frac(mdata = deseq_metadata, 
  #              counts = deseq_df, 
  #              mdata_filter = NULL,
  #              deseq_design = "~0 + site:sampling_time", 
  #              prefix = "usvi_site_vs_time")
  # 
  # #examine the dispersions
  # plotDispEsts(usvi_site_vs_time_deseq_list[["dds"]])
  # 
  # list2env(usvi_site_vs_time_deseq_list, envir = .GlobalEnv)
  # dds <- DESeq2::estimateDispersions(dds, fitType="local")
  # dds <- DESeq2::DESeq(dds, fitType = "local")
  # 
  # plotDispEsts(dds)
  # usvi_site_vs_time_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  #   setNames(., c("metadata", "matrix", "dds", "mod_mat"))
  # 
  # 
  # 
  # 
}



# DESeq2 genus-level, time-specific ---------------------------------------

{
# deseq_df <- usvi_sw_genus.tbl %>%
#   as.data.frame(.) %>%
#   tibble::rownames_to_column(var = "asv_id") %>%
#   tidyr::pivot_longer(., cols = !c(asv_id),
#                       names_to = "sample_name",
#                       values_to = "counts")
# 
# deseq_metadata <- usvi_selected_metadata %>%
#   dplyr::rename(sample_name = "sample_id")
# metadata_filter <- list(`sampling_time` = c("dawn"))
# 
# F_deseq_frac(mdata = deseq_metadata, 
#              counts = deseq_df, 
#              mdata_filter = metadata_filter,
#              deseq_design = "~0 + site", 
#              prefix = "usvi_site_vs_dawn")
# 
# #examine the dispersions
# plotDispEsts(usvi_site_vs_dawn_deseq_list[["dds"]])
# 
# list2env(usvi_site_vs_dawn_deseq_list, envir = .GlobalEnv)
# dds <- DESeq2::estimateDispersions(dds, fitType="local")
# dds <- DESeq2::DESeq(dds, fitType = "local")
# 
# plotDispEsts(dds)
# usvi_site_vs_dawn_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
#   setNames(., c("metadata", "matrix", "dds", "mod_mat"))
# 
# 
# metadata_filter <- list(`sampling_time` = c("peak_photo"))
# 
# F_deseq_frac(mdata = deseq_metadata, 
#              counts = deseq_df, 
#              mdata_filter = metadata_filter,
#              deseq_design = "~0 + site", 
#              prefix = "usvi_site_vs_peak_photo")
# 
# #examine the dispersions
# plotDispEsts(usvi_site_vs_peak_photo_deseq_list[["dds"]])
# 
# list2env(usvi_site_vs_peak_photo_deseq_list, envir = .GlobalEnv)
# dds <- DESeq2::estimateDispersions(dds, fitType="local")
# dds <- DESeq2::DESeq(dds, fitType = "local")
# 
# plotDispEsts(dds)
# usvi_site_vs_peak_photo_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
#   setNames(., c("metadata", "matrix", "dds", "mod_mat"))
# 
# if(!file.exists(paste0(projectpath, "/", "usvi_site_vs_time_deseq.RData"))){
#   save(list = c("usvi_site_vs_time_manual_deseq_list", "usvi_site_vs_time_deseq_list", 
#                 "usvi_site_vs_peak_photo_manual_deseq_list", "usvi_site_vs_peak_photo_deseq_list",
#                 "usvi_site_vs_dawn_manual_deseq_list", "usvi_site_vs_dawn_deseq_list"),
#        file = paste0(projectpath, "/", "usvi_site_vs_time_deseq.RData"))
# }

}

# DESeq2 on all ASVs (not at genus-level) ---------------------------------

usvi_sw_asvs.df <- usvi_prok_asvs.df %>%
  dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = !c(asv_id),
                      names_to = "sample_name",
                      values_to = "counts") %>%
  droplevels

deseq_metadata <- usvi_selected_metadata %>%
  dplyr::rename(sample_name = "sample_id")


F_deseq_frac(mdata = deseq_metadata, 
             counts = usvi_sw_asvs.df, 
             mdata_filter = NULL,
             deseq_design = "~0 + site:sampling_time", 
             prefix = "usvi_asvs_site_vs_time")

#examine the dispersions
# plotDispEsts(usvi_asvs_site_vs_time_deseq_list[["dds"]])
list2env(usvi_asvs_site_vs_time_deseq_list, envir = .GlobalEnv)

# temp_dds <- DESeq2::estimateDispersions(dds, fitType="glmGamPoi")
# plotDispEsts(temp_dds)
# temp_dds2 <- DESeq2::estimateDispersions(dds, fitType="parametric")
# plotDispEsts(temp_dds2)
# temp_dds3 <- DESeq2::estimateDispersions(dds, fitType="local")
# plotDispEsts(temp_dds3)
# temp_dds4 <- DESeq2::estimateDispersions(dds, fitType="mean")
# plotDispEsts(temp_dds4)
# #local seems to have the best fitted 
# rm(temp_dds2)
# rm(temp_dds3)
# rm(temp_dds4)


dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_site_vs_time_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))




metadata_filter <- list(`sampling_time` = c("dawn"))
F_deseq_frac(mdata = deseq_metadata, 
             counts = usvi_sw_asvs.df,
             mdata_filter = metadata_filter,
             deseq_design = "~0 + site", 
             prefix = "usvi_asvs_site_vs_dawn")

#examine the dispersions
# plotDispEsts(usvi_asvs_site_vs_dawn_deseq_list[["dds"]])

list2env(usvi_asvs_site_vs_dawn_deseq_list, envir = .GlobalEnv)
dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_site_vs_dawn_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))




metadata_filter <- list(`sampling_time` = c("peak_photo"))
F_deseq_frac(mdata = deseq_metadata, 
             counts = usvi_sw_asvs.df, 
             mdata_filter = metadata_filter,
             deseq_design = "~0 + site", 
             prefix = "usvi_asvs_site_vs_peak_photo")


#examine the dispersions
# plotDispEsts(usvi_asvs_site_vs_peak_photo_deseq_list[["dds"]])

list2env(usvi_asvs_site_vs_peak_photo_deseq_list, envir = .GlobalEnv)
dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_site_vs_peak_photo_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))

if(!file.exists(paste0(projectpath, "/", "usvi_asvs_site_vs_time_deseq-", Sys.Date(), ".RData"))){
# if(!file.exists(paste0(projectpath, "/", "usvi_asvs_site_vs_time_deseq.RData"))){
  save(list = c("usvi_asvs_site_vs_time_manual_deseq_list", "usvi_asvs_site_vs_time_deseq_list", 
                "usvi_asvs_site_vs_peak_photo_manual_deseq_list", "usvi_asvs_site_vs_peak_photo_deseq_list",
                "usvi_asvs_site_vs_dawn_manual_deseq_list", "usvi_asvs_site_vs_dawn_deseq_list"),
       file = paste0(projectpath, "/", "usvi_asvs_site_vs_time_deseq-", Sys.Date(), ".RData"))
}


# DESeq for interreef comparisons at the same time ------------------------

deseq_metadata <- usvi_selected_metadata %>%
  dplyr::rename(sample_name = "sample_id") %>%
  droplevels

deseq_asvs.df <- usvi_prok_asvs.df %>%
  dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
  dplyr::filter(sample_ID %in% deseq_metadata[["sample_name"]]) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = !c(asv_id),
                      names_to = "sample_name",
                      values_to = "counts") %>%
  droplevels



metadata_filter <- list(`site` = c("Yawzi", "Tektite"))
F_deseq_frac(mdata = deseq_metadata, 
             counts = deseq_asvs.df, 
             mdata_filter = metadata_filter,
             deseq_design = "~0 + sampling_time:site", 
             prefix = "usvi_asvs_yawzi_vs_tektite")

#examine the dispersions
list2env(usvi_asvs_yawzi_vs_tektite_deseq_list, envir = .GlobalEnv)

dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_yawzi_vs_tektite_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))

# DESeq for site-specific dawn vs peak photo -----------------------------------


#repeat the DESeq process on using dawn as basline

metadata_filter <- list(`site` = c("LB_seagrass"))
F_deseq_frac(mdata = deseq_metadata, 
             counts = deseq_asvs.df, 
             mdata_filter = metadata_filter,
             deseq_design = "~0 + sampling_time", 
             prefix = "usvi_asvs_lameshur_dawn_vs_peak_photo")

#examine the dispersions
# plotDispEsts(usvi_asvs_lameshur_dawn_vs_peak_photo_deseq_list[["dds"]])
list2env(usvi_asvs_lameshur_dawn_vs_peak_photo_deseq_list, envir = .GlobalEnv)

dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_lameshur_dawn_vs_peak_photo_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))


#for tektite

metadata_filter <- list(`site` = c("Tektite"))
F_deseq_frac(mdata = deseq_metadata, 
             counts = deseq_asvs.df, 
             mdata_filter = metadata_filter,
             deseq_design = "~0 + sampling_time", 
             prefix = "usvi_asvs_tektite_dawn_vs_peak_photo")

#examine the dispersions
# plotDispEsts(usvi_asvs_tektite_dawn_vs_peak_photo_deseq_list[["dds"]])
list2env(usvi_asvs_tektite_dawn_vs_peak_photo_deseq_list, envir = .GlobalEnv)

dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_tektite_dawn_vs_peak_photo_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))



#for Yawzi

metadata_filter <- list(`site` = c("Yawzi"))
F_deseq_frac(mdata = deseq_metadata, 
             counts = deseq_asvs.df, 
             mdata_filter = metadata_filter,
             deseq_design = "~0 + sampling_time", 
             prefix = "usvi_asvs_yawzi_dawn_vs_peak_photo")

#examine the dispersions
list2env(usvi_asvs_yawzi_dawn_vs_peak_photo_deseq_list, envir = .GlobalEnv)

dds <- DESeq2::estimateDispersions(dds, fitType="local")
dds <- DESeq2::DESeq(dds, fitType = "local")

# plotDispEsts(dds)
usvi_asvs_yawzi_dawn_vs_peak_photo_manual_deseq_list <- list(metadata, matrix, dds, mod_mat) %>%
  setNames(., c("metadata", "matrix", "dds", "mod_mat"))


# Filter sig results ------------------------------------------------------

#for reef_time vs all LB_seagrass:

baseline_filter <- list(`site` = c("LB_seagrass"))

for(namevar in c(
# for(namevar in c("usvi_site_vs_time_manual_deseq_list", "usvi_site_vs_time_deseq_list",
                 "usvi_asvs_site_vs_time_manual_deseq_list", "usvi_asvs_site_vs_time_deseq_list")){
  # for(namevar in c("usvi_asvs_site_vs_time_manual_deseq_list", "usvi_asvs_site_vs_time_deseq_list")){
  prefix <- namevar %>% gsub("_deseq_list", "", .)
  deseq_list <- get0(namevar, )
  list2env(deseq_list, envir = .GlobalEnv)
  
  res <- NULL
  #identify what the baseline should consist of
  baseline_idx <- purrr::map2(baseline_filter, names(baseline_filter),
                              ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
    purrr::reduce(inner_join)
  
  baseline_vec <- colMeans(mod_mat[baseline_idx, ])
  #identify the conditions you are testing against the baseline, and determine interactions = contrasts
  condition_filter <- list(`site` = c("Yawzi", "Tektite"),
                           `sampling_time` = c("dawn", "peak_photo"))
  
  condition_idx_names <- purrr::map(condition_filter,
                                    ~unlist(., recursive = FALSE)) %>%
    purrr::reduce(interaction, sep = "_") %>%
    levels %>%
    as.character
  
  #brute force: list the model vectors
  for(i in c("Yawzi", "Tektite")){
    namevar <- i
    for(j in c("dawn", "peak_photo")){
      lightvar <- j
      if(is.null(dim(mod_mat[dds$site == namevar & dds$sampling_time == lightvar,]))){
        temp_vec <- mod_mat[dds$site == namevar & dds$sampling_time == lightvar, ]
      } else {
        temp_vec <- colMeans(mod_mat[dds$site == namevar & dds$sampling_time == lightvar, ])
      }
      assign(paste0(namevar, "_", lightvar), temp_vec, envir = .GlobalEnv)
      rm(temp_vec)
      rm(lightvar)
    }
    rm(namevar)
  }
  
  contrast_id <- condition_idx_names %>%
    as.data.frame() %>%
    dplyr::mutate(baseline = "baseline_vec") %>%
    dplyr::mutate(across(everything(), as.character))
  
  label_id <- paste0(condition_idx_names, " - LB_seagrass")
  
  y <- nrow(contrast_id)
  for(i in seq(y)){
    Sys.sleep(1/100)
    pair1 <- get0(contrast_id[i,1], inherits = TRUE)
    pair2 <- get0(contrast_id[i,2], inherits = TRUE)
    res1 <- DESeq2::results(dds,
                            contrast = pair1 - pair2,
                            parallel = TRUE,
                            BPPARAM = bpparam_multi,
                            alpha = set.alpha) %>%
      as.data.frame(.) %>%
      tibble::rownames_to_column("asv_id") %>%
      dplyr::left_join(., DESeq2::lfcShrink(dds,
                                            contrast = pair1 - pair2,
                                            parallel = TRUE,
                                            BPPARAM = bpparam_multi,
                                            type = "ashr",
                                            quiet = TRUE) %>%
                         as.data.frame(.) %>%
                         tibble::rownames_to_column("asv_id") %>%
                         dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                         dplyr::select(asv_id, log2FoldChange_MMSE),
                       by = "asv_id") %>%
      dplyr::mutate(contrast = label_id[i]) %>%
      dplyr::relocate(contrast, .before = "asv_id")
    if(exists("res")){
      res <- bind_rows(res, res1)
      rm(res1)}
    if(!exists("res")){
      res <- res1
      rm(res1)}
  }
  assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
  
  readr::write_rds(res, file = paste0(projectpath, "/", prefix, "_res_list", ".rds"))
  
}



# Extract results for time- and site-specific -----------------------------



baseline_filter <- list(`site` = c("LB_seagrass"))
for(delist in c(
# for(delist in c("usvi_site_vs_peak_photo_manual_deseq_list", "usvi_site_vs_peak_photo_deseq_list",
#                 "usvi_site_vs_dawn_manual_deseq_list", "usvi_site_vs_dawn_deseq_list",
                "usvi_asvs_site_vs_peak_photo_manual_deseq_list", "usvi_asvs_site_vs_peak_photo_deseq_list",
                "usvi_asvs_site_vs_dawn_manual_deseq_list", "usvi_asvs_site_vs_dawn_deseq_list")){
  # for(delist in c("usvi_asvs_site_vs_peak_photo_manual_deseq_list", "usvi_asvs_site_vs_peak_photo_deseq_list",
#                 "usvi_asvs_site_vs_dawn_manual_deseq_list", "usvi_asvs_site_vs_dawn_deseq_list")){
  
  prefix <- delist %>% gsub("_deseq_list", "", .)
  deseq_list <- get0(delist, )
  list2env(deseq_list, envir = .GlobalEnv)
  rm(deseq_list)
  # { 
    # prefix <- "usvi_site_vs_dawn_manual"
  # list2env(usvi_site_vs_dawn_manual_deseq_list, envir = .GlobalEnv)
  tempid <- stringr::str_extract(prefix, "dawn|peak_photo")
  res <- NULL
  #identify what the baseline should consist of
  baseline_idx <- purrr::map2(baseline_filter, names(baseline_filter),
                              ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
    purrr::reduce(inner_join)
  
  baseline_vec <- colMeans(mod_mat[baseline_idx, ])
  #identify the conditions you are testing against the baseline, and determine interactions = contrasts
  condition_filter <- list(`site` = c("Yawzi", "Tektite"),
                           `sampling_time` = tempid)
  
  condition_idx_names <- purrr::map(condition_filter,
                                    ~unlist(., recursive = FALSE)) %>%
    purrr::reduce(interaction, sep = "_") %>%
    levels %>%
    as.character
  
  #brute force: list the model vectors
  for(i in c("Yawzi", "Tektite")){
    namevar <- i
    for(j in tempid){
      lightvar <- j
      if(is.null(dim(mod_mat[dds$site == namevar & dds$sampling_time == lightvar,]))){
        temp_vec <- mod_mat[dds$site == namevar & dds$sampling_time == lightvar, ]
      } else {
        temp_vec <- colMeans(mod_mat[dds$site == namevar & dds$sampling_time == lightvar, ])
      }
      assign(paste0(namevar, "_", lightvar), temp_vec, envir = .GlobalEnv)
      rm(temp_vec)
      rm(lightvar)
    }
    rm(namevar)
  }
  
  contrast_id <- condition_idx_names %>%
    as.data.frame() %>%
    dplyr::mutate(baseline = "baseline_vec") %>%
    dplyr::mutate(across(everything(), as.character))
  
  label_id <- paste0(condition_idx_names, " - LB_seagrass_", tempid)
  
  y <- nrow(contrast_id)
  for(i in seq(y)){
    Sys.sleep(1/100)
    pair1 <- get0(contrast_id[i,1], inherits = TRUE)
    pair2 <- get0(contrast_id[i,2], inherits = TRUE)
    res1 <- DESeq2::results(dds,
                            contrast = pair1 - pair2,
                            parallel = TRUE,
                            BPPARAM = bpparam_multi,
                            alpha = set.alpha) %>%
      as.data.frame(.) %>%
      tibble::rownames_to_column("asv_id") %>%
      dplyr::left_join(., DESeq2::lfcShrink(dds,
                                            contrast = pair1 - pair2,
                                            parallel = TRUE,
                                            BPPARAM = bpparam_multi,
                                            type = "ashr",
                                            quiet = TRUE) %>%
                         as.data.frame(.) %>%
                         tibble::rownames_to_column("asv_id") %>%
                         dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                         dplyr::select(asv_id, log2FoldChange_MMSE),
                       by = "asv_id") %>%
      dplyr::mutate(contrast = label_id[i]) %>%
      dplyr::relocate(contrast, .before = "asv_id")
    if(exists("res")){
      res <- bind_rows(res, res1)
      rm(res1)}
    if(!exists("res")){
      res <- res1
      rm(res1)}
  }
  assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
  
  readr::write_rds(res, file = paste0(projectpath, "/", prefix, "_res_list", ".rds"))
  rm(prefix)
}



# For site-specific ---------------------------------------------------------

baseline_filter <- list(`sampling_time` = c("dawn"))
for(delist in c("usvi_asvs_lameshur_dawn_vs_peak_photo_manual_deseq_list", "usvi_asvs_lameshur_dawn_vs_peak_photo_deseq_list",
                "usvi_asvs_tektite_dawn_vs_peak_photo_manual_deseq_list", "usvi_asvs_tektite_dawn_vs_peak_photo_deseq_list",
                "usvi_asvs_yawzi_dawn_vs_peak_photo_manual_deseq_list", "usvi_asvs_yawzi_dawn_vs_peak_photo_deseq_list")){
# { delist <- "usvi_asvs_lameshur_dawn_vs_peak_photo_manual_deseq_list"
  prefix <- delist %>% gsub("_deseq_list", "", .)
  deseq_list <- get0(delist, )
  list2env(deseq_list, envir = .GlobalEnv)
  rm(deseq_list)
  tempid <- stringr::str_extract(prefix, "lameshur|yawzi|tektite")
  if("lameshur" %in% tempid){
    tempid <- "LB_seagrass"
  } else {
    tempid <- stringr::str_to_title(tempid)
  }
  
  res <- NULL
  #identify what the baseline should consist of
  baseline_idx <- purrr::map2(baseline_filter, names(baseline_filter),
                              ~grepl(paste0(.x, collapse = "|"), dds[[.y]], ignore.case = TRUE, perl = TRUE)) %>%
    purrr::reduce(inner_join)
  
  baseline_vec <- colMeans(mod_mat[baseline_idx, ])
  #identify the conditions you are testing against the baseline, and determine interactions = contrasts
  condition_filter <- list(`site` = tempid,
                           `sampling_time` = c("peak_photo"))
  
  condition_idx_names <- purrr::map(condition_filter,
                                    ~unlist(., recursive = FALSE)) %>%
    purrr::reduce(interaction, sep = "_") %>%
    levels %>%
    as.character
  
  #brute force: list the model vectors
  for(i in c("peak_photo")){
  # for(i in c("dawn", "peak_photo")){
    lightvar <- i
      sitevar <- tempid
      if(is.null(dim(mod_mat[dds$site == sitevar & dds$sampling_time == lightvar,]))){
        temp_vec <- mod_mat[dds$site == sitevar & dds$sampling_time == lightvar, ]
      } else {
        temp_vec <- colMeans(mod_mat[dds$site == sitevar & dds$sampling_time == lightvar, ])
      }
      assign(paste0(sitevar, "_", lightvar), temp_vec, envir = .GlobalEnv)
      rm(temp_vec)
      rm(lightvar)
    rm(sitevar)
  }
  
  contrast_id <- condition_idx_names %>%
    as.data.frame() %>%
    dplyr::mutate(baseline = "baseline_vec") %>%
    dplyr::mutate(across(everything(), as.character))
  
  label_id <- paste0(condition_idx_names, " - ", tempid, "_dawn")
  
  y <- nrow(contrast_id)
  for(i in seq(y)){
    Sys.sleep(1/100)
    pair1 <- get0(contrast_id[i,1], inherits = TRUE)
    pair2 <- get0(contrast_id[i,2], inherits = TRUE)
    res1 <- DESeq2::results(dds,
                            contrast = pair1 - pair2,
                            parallel = TRUE,
                            BPPARAM = bpparam_multi,
                            alpha = set.alpha) %>%
      as.data.frame(.) %>%
      tibble::rownames_to_column("asv_id") %>%
      dplyr::left_join(., DESeq2::lfcShrink(dds,
                                            contrast = pair1 - pair2,
                                            parallel = TRUE,
                                            BPPARAM = bpparam_multi,
                                            type = "ashr",
                                            quiet = TRUE) %>%
                         as.data.frame(.) %>%
                         tibble::rownames_to_column("asv_id") %>%
                         dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                         dplyr::select(asv_id, log2FoldChange_MMSE),
                       by = "asv_id") %>%
      dplyr::mutate(contrast = label_id[i]) %>%
      dplyr::relocate(contrast, .before = "asv_id")
    if(exists("res")){
      res <- bind_rows(res, res1)
      rm(res1)}
    if(!exists("res")){
      res <- res1
      rm(res1)}
  }
  assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
  
  readr::write_rds(res, file = paste0(projectpath, "/", prefix, "_res_list", ".rds"))
  rm(prefix)
}


# Inter-reef comparisons --------------------------------------------------

baseline_filter <- list(`site` = c("Tektite"))
for(delist in c("usvi_asvs_yawzi_vs_tektite_manual_deseq_list", "usvi_asvs_yawzi_vs_tektite_deseq_list")){
  # {   delist <- "usvi_asvs_yawzi_vs_tektite_manual_deseq_list"
  prefix <- delist %>% gsub("_deseq_list", "", .)
  deseq_list <- get0(delist, )
  list2env(deseq_list, envir = .GlobalEnv)
  rm(deseq_list)
  
  res <- NULL
  #identify what the baseline should consist of
  #baseline will be the corresponding time at Tektite.
  
  #identify the conditions you are testing against the baseline, and determine interactions = contrasts
  condition_filter <- list(`site` = c("Yawzi", "Tektite"),
                           `sampling_time` = c("dawn", "peak_photo"))
  
  condition_idx_names <- purrr::map(condition_filter,
                                    ~unlist(., recursive = FALSE)) %>%
    purrr::reduce(interaction, sep = "_") %>%
    levels %>%
    rev(.) %>%
    as.character
  
  #brute force: list the model vectors
  for(i in c("dawn", "peak_photo")){
    lightvar <- i
    for(j in c("Yawzi", "Tektite")){
      sitevar <- j
      if(is.null(dim(mod_mat[dds$site == sitevar & dds$sampling_time == lightvar,]))){
        temp_vec <- mod_mat[dds$site == sitevar & dds$sampling_time == lightvar, ]
      } else {
        temp_vec <- colMeans(mod_mat[dds$site == sitevar & dds$sampling_time == lightvar, ])
      }
      assign(paste0(sitevar, "_", lightvar), temp_vec, envir = .GlobalEnv)
      rm(temp_vec)  
      rm(sitevar)
    }
    rm(lightvar)
  }
  
  
  
  contrast_id <- condition_idx_names %>%
    matrix(., nrow = 2, byrow = TRUE) %>%
    as.data.frame() %>%
    # dplyr::mutate(baseline = "baseline_vec") %>%
    dplyr::mutate(across(everything(), as.character))
  
  label_id <- contrast_id %>%
    tidyr::unite(., label_id, sep = " - ", remove = TRUE) %>%
    unlist %>% as.character
  
  y <- nrow(contrast_id)
  for(i in seq(y)){
    Sys.sleep(1/100)
    pair1 <- get0(contrast_id[i,1], inherits = TRUE)
    pair2 <- get0(contrast_id[i,2], inherits = TRUE)
    res1 <- DESeq2::results(dds,
                            contrast = pair1 - pair2,
                            parallel = TRUE,
                            BPPARAM = bpparam_multi,
                            alpha = set.alpha) %>%
      as.data.frame(.) %>%
      tibble::rownames_to_column("asv_id") %>%
      dplyr::left_join(., DESeq2::lfcShrink(dds,
                                            contrast = pair1 - pair2,
                                            parallel = TRUE,
                                            BPPARAM = bpparam_multi,
                                            type = "ashr",
                                            quiet = TRUE) %>%
                         as.data.frame(.) %>%
                         tibble::rownames_to_column("asv_id") %>%
                         dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                         dplyr::select(asv_id, log2FoldChange_MMSE),
                       by = "asv_id") %>%
      dplyr::mutate(contrast = label_id[i]) %>%
      dplyr::relocate(contrast, .before = "asv_id")
    if(exists("res")){
      res <- bind_rows(res, res1)
      rm(res1)}
    if(!exists("res")){
      res <- res1
      rm(res1)}
  }
  assign(paste0(prefix, "_res_list"), res, envir = .GlobalEnv, inherits = TRUE)
  
  readr::write_rds(res, file = paste0(projectpath, "/", prefix, "_res_list", ".rds"))
  rm(prefix)
}

# Extract results ---------------------------------------------------------

rlist2 <- list(usvi_asvs_lameshur_dawn_vs_peak_photo_manual_res_list,
               usvi_asvs_yawzi_dawn_vs_peak_photo_manual_res_list,
               usvi_asvs_tektite_dawn_vs_peak_photo_manual_res_list,
               
               usvi_asvs_lameshur_dawn_vs_peak_photo_res_list,
               usvi_asvs_yawzi_dawn_vs_peak_photo_res_list,
               usvi_asvs_tektite_dawn_vs_peak_photo_res_list) %>%
  setNames(., c("manual_lameshur", "manual_yawzi", "manual_tektite", 
                
                "auto_lameshur", "auto_yawzi", "auto_tektite")) %>%
  bind_rows(., .id = "model")
# hist(rlist2[["pvalue"]])
rlist <- list(usvi_asvs_site_vs_time_manual_res_list,
              usvi_asvs_site_vs_dawn_manual_res_list,
              usvi_asvs_site_vs_peak_photo_manual_res_list,
              usvi_asvs_yawzi_vs_tektite_manual_res_list,
              usvi_asvs_site_vs_time_res_list,
              usvi_asvs_site_vs_dawn_res_list,
              usvi_asvs_site_vs_peak_photo_res_list, 
              usvi_asvs_yawzi_vs_tektite_res_list) %>%
  setNames(., c("manual_dispersion", "manual_dawn", "manual_photo", "manual_reef",
                "auto_dispersion", "auto_dawn", "auto_photo","auto_reef")) %>%
  bind_rows(., .id = "model")

rlist_pvalues <- c(rlist[["pvalue"]], rlist2[["pvalue"]])

hist(rlist_pvalues)

#although DESeq does multiple test corrections, their adjusted p-value cutoff is by default 0.1
#determine based on the number of tests we did, what the alpha should be for a q-value of 0.1 

q_value <- 0.1
padj_cutoff <- rlist_pvalues %>%
# padj_cutoff <- rlist2[["pvalue"]] %>%
  p.adjust(., method = "BH") %>% #multiple testing corrections, "BH" is an alias for "fdr" accordiong to p.adjust()
  na.omit(.) %>%
  ashr::qval.from.lfdr(.) %>%
  unlist %>%
  quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)

usvi_deseq_asvs_res.df <- rlist %>%
  dplyr::group_by(asv_id) %>%
  dplyr::filter(pvalue <= 0.05) %>%
  droplevels %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(padj_qcutoff = dplyr::case_when((padj < padj_cutoff) ~ padj,
                                                .default = NA)) %>%
  droplevels

usvi_deseq_asvs_intrasite_res.df <- rlist2 %>%
  dplyr::group_by(asv_id) %>%
  dplyr::filter(pvalue <= 0.05) %>%
  droplevels %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(padj_qcutoff = dplyr::case_when((padj < padj_cutoff) ~ padj,
                                                .default = NA)) %>%
  droplevels

contrast_labels2 <- usvi_deseq_asvs_intrasite_res.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(contrast, .keep_all = FALSE) %>%
  dplyr::mutate(contrast = factor(contrast, levels = unique(.[["contrast"]]))) %>%
  droplevels

contrast_labels <- usvi_deseq_asvs_res.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(contrast, .keep_all = FALSE) %>%
  dplyr::mutate(contrast = factor(contrast, levels = unique(.[["contrast"]]))) %>%
  bind_rows(., contrast_labels2) %>%
  droplevels




# Which genera were SDA? ---------------------------------------------

#don't need to examine this anymore, since we're not doing genera-level
{
# if(any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
#   temp_file <- list.files(projectpath, pattern = "usvi_deseq_genera_site_time-.*.RData")[1]
#   load(paste0(projectpath, "/", temp_file))
#   rm(temp_file)
# } else {
#   usvi_deseq_genera_res.taxa.df <- usvi_sw_genus.taxa.df %>%
#     dplyr::right_join(., usvi_deseq_genera_res.df %>%
#                         dplyr::distinct(asv_id, .keep_all = FALSE),
#                       by = join_by("asv_id"))
#   
#   usvi_deseq_genera_res.list <- usvi_deseq_genera_res.df %>%
#     dplyr::filter(grepl("manual", model)) %>%
#     tidyr::drop_na(padj_qcutoff) %>%
#     dplyr::select(model, contrast, asv_id, baseMean, padj, log2FoldChange_MMSE, padj_qcutoff) %>%
#     droplevels %>%
#     split(., f = .$asv_id) %>%
#     map(., ~.x %>%
#           dplyr::mutate(sampling_time = stringr::str_extract(contrast, "dawn|peak_photo")) %>%
#           dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi_|Tektite_") %>%
#                           stringr::str_remove(., "_")) %>%
#           dplyr::mutate(Combo = stringr::str_extract(contrast, "Yawzi_.|Tektite_.")) %>%
#           dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
#           # dplyr::mutate(baseline = dplyr::case_when(grepl("dawn|peak", baseline) ~ baseline,
#           dplyr::mutate(baseline = dplyr::case_when(grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
#                                                     .default = "LB_seagrass_all")) %>%
#           droplevels)
#   
#   
#   #so this analysis identified genera that were significantly differentially abundant in Tektite and Yawzi
#   #at dawn and peak photosynthesis sampling times
#   #but these taxa were also called "significant" if they were just more noticeably more (or less) abundant in seagrass
#   #and show up as enriched in both Yawzi and Tektite, or depleted in both. 
#   
#   #which genera were not just differentially abundant between the reef sites at different sampling times?
#   #remove those genera that were num_enriched = 4 or num_depleted = 4
#   #then look at those genera that had an ODD count for num_enriched or num_depleted
#   
#   usvi_deseq_genera_summary.df <- usvi_deseq_genera_res.list %>%
#     map(., ~.x %>%
#           dplyr::mutate(enriched = dplyr::case_when((log2FoldChange_MMSE > 0) ~ 1,
#                                                     .default = NA)) %>%
#           dplyr::mutate(depleted = dplyr::case_when((log2FoldChange_MMSE < 0) ~ 1,
#                                                     .default = NA)) %>%
#           dplyr::group_by(baseline) %>%
#           dplyr::summarise(num_enriched = sum(enriched, na.rm = TRUE),
#                            num_depleted = sum(depleted, na.rm = TRUE))) %>%
#     dplyr::bind_rows(., .id = "asv_id")
#   
#   usvi_deseq_genera_abund.df <- usvi_deseq_genera_res.list %>%
#     dplyr::bind_rows(., .id = "asv_id") %>%
#     dplyr::filter(grepl("manual", model)) %>%
#     droplevels %>%
#     dplyr::semi_join(., (usvi_deseq_genera_summary.df %>%
#                            dplyr::filter(if_all(c(num_enriched, num_depleted), ~.x < 4)) %>%
#                            dplyr::filter(if_any(c(num_enriched, num_depleted), ~gtools::odd(.x))) %>%
#                            dplyr::distinct(asv_id, baseline, .keep_all = FALSE)),
#                      by = join_by(asv_id, baseline)) %>%
#     # dplyr::filter(grepl("dawn|peak", baseline)) %>%
#     droplevels
#   
#   
#   # print(
#   #   ggplot(data = usvi_deseq_genera_abund.df %>%
#   #            droplevels)
#   #   + theme_bw() 
#   #   + geom_hline(aes(yintercept = 0), color = "black")
#   #   + geom_segment(aes(y = log2FoldChange_MMSE, yend = 0, x = asv_id, group = asv_id, 
#   #                      color = site), show.legend = FALSE)
#   #   + geom_point(aes(y = log2FoldChange_MMSE, x = asv_id, group = asv_id, fill = site),
#   #                alpha = 0.8, shape = 21, color = "black", size = 3, show.legend = TRUE, 
#   #                position = position_dodge(0.8))
#   #   # + geom_boxplot(aes(group = asv_id), alpha = 0.6, show.legend = FALSE,
#   #   #                position = position_dodge(0.8), fill = "grey", color = "black")
#   #   
#   #   + scale_fill_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
#   #   + scale_color_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
#   #   + scale_x_discrete(labels = usvi_genera_relabel)
#   #   + scale_y_continuous(name = "log2(fold change)")
#   #   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#   #                                override.aes = list(color = "black", shape = 21, size = 3)),
#   #            color = "none")
#   #   + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
#   #           strip.text.y = element_text(angle = 0))
#   #   + facet_wrap(baseline ~ sampling_time,
#   #                # ncol = 4, 
#   #                drop = TRUE,
#   #                # space = "fixed",
#   #                scales = "free_x", labeller = global_labeller)
#   #   + ggtitle("Significantly differentially abundant from LB seagrass samples")
#   # )
#   
#   #all the genera that were differentially abundant in seagrass regardless of time:
#   temp_df <- (usvi_deseq_genera_abund.df %>%
#                 dplyr::distinct(contrast, asv_id, site, sampling_time, .keep_all = TRUE) %>%
#                 dplyr::group_by(asv_id, site, model) %>%
#                 dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
#                 dplyr::arrange(desc(obs)) %>%
#                 dplyr::filter(gtools::even(obs)) %>%
#                 droplevels)
#   usvi_deseq_genera_abund_filtered.df <- usvi_deseq_genera_abund.df %>%
#     dplyr::anti_join(., (usvi_deseq_genera_abund.df %>%
#                            dplyr::distinct(contrast, asv_id, site, sampling_time, .keep_all = TRUE) %>%
#                            dplyr::group_by(asv_id, site, model) %>%
#                            dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
#                            dplyr::arrange(desc(obs)) %>%
#                            dplyr::filter(gtools::even(obs)) %>%
#                            droplevels), 
#                      # relationship = "many-to-many", multiple = "all",
#                      by = join_by(asv_id)) %>%
#     droplevels
# }
# 
# 
# # usvi_deseq_genera_abund_filtered.df <- usvi_deseq_genera_abund.df %>%
# #   # dplyr::inner_join(., (usvi_deseq_genera_abund.df %>%
# #                           dplyr::semi_join(., (usvi_deseq_genera_abund.df %>%
# #                           dplyr::group_by(asv_id, baseline, site) %>%
# #                           dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
# #                           dplyr::filter(gtools::odd(obs))), 
# #                     # relationship = "many-to-many", multiple = "all",
# #                     by = join_by(asv_id, baseline, site)) %>%
# #   droplevels
# 
# #these are the taxa that were flagged as SDA in a reef_time, but are just different from LB Seagrass as a whole
# 
# # temp_df3 <- usvi_deseq_genera_abund.df %>%
# #   dplyr::anti_join(., (usvi_deseq_genera_abund.df %>%
# #                          dplyr::group_by(asv_id, baseline, site) %>%
# #                          dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
# #                          dplyr::filter(gtools::odd(obs))), by = join_by(asv_id, baseline, site)) %>%
# #   droplevels
# g3 <- print(
#   ggplot(data = usvi_deseq_genera_abund_filtered.df %>%
#            droplevels)
#   + theme_bw() 
#   + geom_hline(aes(yintercept = 0), color = "black")
#   + geom_segment(aes(y = log2FoldChange_MMSE, yend = 0, x = asv_id, group = asv_id, 
#                      color = sampling_time), position = position_jitter(width = 0.2, seed = 48105), show.legend = FALSE)
#   + geom_point(aes(y = log2FoldChange_MMSE, x = asv_id, group = asv_id, fill = sampling_time, shape = sampling_time),
#                alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
#                position = position_jitter(width = 0.2, seed = 48105))
#   + scale_shape_manual(name = "Sampling site and time",
#                        values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
#   + scale_fill_manual(name = "Sampling site and time",
#                       values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_color_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_x_discrete(labels = usvi_genera_relabel)
#   + scale_y_continuous(name = "log2(fold change)")
#   + guides(color = "none")
#   + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
#           strip.text.y = element_text(angle = 0))
#   + facet_grid(baseline ~ site,
#                drop = TRUE,
#                scales = "fixed", labeller = global_labeller)
#   + ggtitle("Significantly differentially abundant from LB seagrass samples")
# )
# 
# # usvi_site_time_SDA_summary <- rlist %>%
# #   dplyr::rowwise(.) %>%
# #   dplyr::filter(padj < padj_cutoff) %>%
# #   dplyr::ungroup(.) %>%
# # dplyr::select(c(model, contrast, asv_id, padj)) %>%
# #   dplyr::mutate(contrast = factor(contrast, levels = unique(.[["contrast"]]))) %>%
# #   dplyr::mutate(model = factor(model)) %>%
# #   tidyr::pivot_longer(., cols = !c(model, "contrast", "asv_id"),
# #                       names_to = NULL,
# #                       values_to = "padj") %>%
# #   dplyr::ungroup(.) %>%
# #   tidyr::drop_na(.) %>%
# #   dplyr::distinct(model, contrast, asv_id, .keep_all = TRUE) %>%
# #   dplyr::group_by(contrast, model) %>%
# #   dplyr::summarise(Num_SDA = length(asv_id), .groups = "keep") %>%
# #   dplyr::full_join(., contrast_labels %>%
# #                      dplyr::select(contrast), by = "contrast") %>%
# #   dplyr::mutate(Num_SDA = tidyr::replace_na(Num_SDA, 0)) %>%
# #   dplyr::arrange(desc(Num_SDA)) %>%
# #   droplevels
# 
# # usvi_site_time_SDA_summary <- rlist %>%
# usvi_site_time_SDA_summary <- usvi_deseq_genera_res.df %>%
#   dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
#   dplyr::mutate(baseline = dplyr::case_when(grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
#                                             .default = "LB_seagrass_all")) %>%
#   dplyr::rowwise(.) %>%
#   dplyr::filter(padj < padj_cutoff) %>%
#   droplevels %>%
#   dplyr::semi_join(., (usvi_deseq_genera_summary.df %>%
#                          dplyr::filter(if_all(c(num_enriched, num_depleted), ~.x < 4)) %>%
#                          dplyr::filter(if_any(c(num_enriched, num_depleted), ~gtools::odd(.x))) %>%
#                          dplyr::distinct(asv_id, baseline, .keep_all = FALSE)),
#                    by = join_by(asv_id, baseline)) %>%
#   dplyr::select(c(model, contrast, asv_id, padj)) %>%
#   dplyr::mutate(contrast = factor(contrast, levels = unique(.[["contrast"]]))) %>%
#   # dplyr::mutate(model = factor(model)) %>%
#   dplyr::distinct(model, contrast, asv_id, .keep_all = TRUE) %>%
#   dplyr::group_by(contrast, model) %>%
#   dplyr::summarise(Num_SDA = length(asv_id), .groups = "keep") %>%
#   dplyr::full_join(., contrast_labels %>%
#                      dplyr::select(contrast), by = "contrast") %>%
#   dplyr::mutate(Num_SDA = tidyr::replace_na(Num_SDA, 0)) %>%
#   dplyr::arrange(desc(Num_SDA)) %>%
#   droplevels
# 
# 
# 
# #using all LB seagrass samples for baseline, yielded a higher number of SDA ASVs than using the time-specific LB seagrass samples for each comparison
# #using manual dispersion did not substantially change the number of SDA ASVs
# 
# 
# {
#   temp_df <- usvi_site_time_SDA_summary %>%
#     dplyr::mutate(sampling_time = stringr::str_extract(contrast, "dawn|peak_photo")) %>%
#     dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi_|Tektite_") %>%
#                     stringr::str_remove(., "_")) %>%
#     dplyr::mutate(Combo = stringr::str_extract(contrast, "Yawzi_.|Tektite_.")) %>%
#     dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
#     dplyr::mutate(baseline = dplyr::case_when(grepl("dawn|peak", baseline) ~ baseline,
#                                               .default = "LB_seagrass_all")) %>%
#     dplyr::mutate(label_y = dplyr::case_when( (Num_SDA < 10) ~ 10, 
#                                               .default = (Num_SDA+10)*1.1)) %>%
#     dplyr::mutate(model2 = recode(model, !!!model_dispersion_lookup)) %>%
#     droplevels
#   
#   g1 <- print(
#     ggplot(data = temp_df %>%
#              droplevels, aes(y = Num_SDA, 
#                              x = Combo, group = contrast, 
#                              fill = Combo))
#     + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
#     + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
#                position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
#     + theme_bw() 
#     + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#     + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
#             strip.text.y = element_text(angle = 0))
#     + scale_fill_manual(values = group_labels_colors, labels = group_labels_lookup, breaks = names(group_labels_lookup))
#     + scale_x_discrete(labels = group_labels_lookup)
#     + geom_text(aes(label = Num_SDA, x = Combo, vjust = "outward", hjust = "center", group = contrast,
#                     y = label_y),
#                 colour = "black", fontface = "bold")
#     + labs(x = "Sample type", y = "Number of SDA genera")
#     # + facet_grid(rows = vars(model), space = "fixed", scales = "fixed", labeller = labeller(model = model_dispersion_lookup, .multi_line = TRUE))
#     + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                  override.aes = list(color = "black", shape = 21, size = 3)),
#              color = "none")
#     + ggtitle("Baseline abundance from LB seagrass samples")
#     + facet_wrap(model ~ .,
#                  # ncol = 4, 
#                  drop = TRUE,
#                  # space = "fixed",
#                  scales = "fixed", labeller = global_labeller)
#   )
#   
#   }
# 
# if(!any(grepl("genera_dawn_photo_sda", list.files(projectpath, pattern = "usvi_deseq_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_deseq_genera_dawn_photo_sda-", Sys.Date(), ".png"),
#          g1,
#          width = 10,  height = 8, units = "in")
# }
# 
# 
# 
# 
# if(!any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
#   save(usvi_deseq_genera_res.list, usvi_deseq_genera_abund_filtered.df, usvi_deseq_genera_abund.df, usvi_deseq_genera_res.df, file = paste0(projectpath, "/", "usvi_deseq_genera_site_time-", Sys.Date(), ".RData"))
# }
}

# Determine whether to use manual or auto dispersions ---------------------

#determine whether to proceed with manual or automatica dispersions in the models


usvi_asvs_site_time_SDA_summary <- bind_rows(usvi_deseq_asvs_res.df, usvi_deseq_asvs_intrasite_res.df) %>%
  tidyr::drop_na(padj_qcutoff) %>%
  dplyr::filter(!(padj > 0.05)) %>%
  # usvi_asvs_site_time_SDA_summary <- usvi_deseq_asvs_res.df %>%
  dplyr::select(c(model, contrast, asv_id, padj)) %>%
  dplyr::mutate(contrast = factor(contrast, levels = unique(.[["contrast"]]))) %>%
  tidyr::pivot_longer(., cols = !c(model, "contrast", "asv_id"),
                      names_to = "metric",
                      values_to = "padj") %>%
  tidyr::drop_na(.) %>%
  dplyr::filter(padj < padj_cutoff) %>%
  # dplyr::filter(padj < set.alpha) %>%
  droplevels %>%
  dplyr::group_by(contrast, model) %>%
  dplyr::summarise(Num_SDA = length(asv_id), .groups = "keep") %>%
  dplyr::bind_rows(data.frame(contrast = c("Yawzi_peak_photo - Yawzi_dawn"), model = c("auto_yawzi"))) %>% #since many of the identified ASVs were not SDA
  # dplyr::full_join(., contrast_labels %>%
  #                    dplyr::select(contrast), by = "contrast") %>%
  dplyr::mutate(Num_SDA = tidyr::replace_na(Num_SDA, 0)) %>%
  dplyr::arrange(desc(Num_SDA)) %>%
  droplevels

#plot the number of SDA ASVs when using manually-calculated dispersions, or model-automatically calculated
{
  temp_df <- usvi_asvs_site_time_SDA_summary %>%
    dplyr::filter(!grepl("dispersion", model)) %>%
    dplyr::mutate(sampling_time = stringr::str_extract(contrast, "dawn|peak_photo")) %>%
    dplyr::mutate(site = stringr::str_extract(contrast, "^LB_seagrass_|Yawzi_|Tektite_") %>%
                    stringr::str_remove(., "_$")) %>%
    dplyr::mutate(Combo = stringr::str_extract(contrast, "LB_seagrass_.|Yawzi_.|Tektite_.")) %>%
    dplyr::mutate(dispersion_model = stringr::str_extract_all(model, "(manual|auto)") %>%
                    stringr::str_trim(., "both")) %>%
    dplyr::mutate(hold = stringr::str_remove_all(model, "manual_|auto_")) %>%
    dplyr::mutate(hold = dplyr::case_when(grepl("lameshur", hold) ~ "LB_seagrass",
                                          grepl("photo", hold) ~ "peak_photo",
                                          grepl("yawzi|tektite|lb", hold) ~ stringr::str_to_sentence(hold),
                                          grepl("reef", model) ~ sampling_time,
                                          .default = hold)) %>%
    dplyr::mutate(variable = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ stringr::str_remove_all(contrast, paste0("_", hold)),
                                              (hold %in% names(site_lookup)) ~ stringr::str_remove_all(contrast, paste0(site, "_")),
                                              .default = NA)) %>%
    dplyr::mutate(hold2 = dplyr::case_when(grepl("reef", model) ~ "interreef",
                                           grepl("dawn|photo", hold) ~ "sampling time",
                                           .default = "site")) %>%
    dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
    dplyr::mutate(baseline2 = stringr::str_remove_all(baseline, "_dawn|_peak_photo")) %>%
    dplyr::mutate(label_y = dplyr::case_when( (Num_SDA < 10) ~ 20, 
                                              .default = (Num_SDA+10)*1.1)) %>%
    dplyr::mutate(model2 = recode(model, !!!model_dispersion_lookup)) %>%
    droplevels
  
  # g2_temporal <- print(
  #   ggplot(data = temp_df %>%
  #            dplyr::filter(grepl("time", hold2)) %>%
  #            droplevels, aes(y = Num_SDA, 
  #                            group = site, x = site, fill = site))
  #                            # group = contrast, x = Combo, fill = Combo))
  #   + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
  #   + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
  #              position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
  #   + theme_bw() 
  #   + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  #   + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
  #           strip.text.y = element_text(angle = 0))
  #   + scale_fill_manual(values = group_labels_colors, labels = group_labels_lookup, breaks = names(group_labels_lookup))
  #   # + scale_x_discrete(labels = group_labels_lookup)
  #   + geom_text(aes(label = Num_SDA, vjust = "outward", hjust = "center", 
  #                   x = site, group = contrast,
  #                   # x = Combo, group = contrast,
  #                   y = label_y),
  #               colour = "black", fontface = "bold")
  #   + labs(x = "Sample type", y = "Number of SDA ASVs")
  #   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
  #                                override.aes = list(color = "black", shape = 21, size = 3)),
  #            color = "none")
  #   + ggtitle("Baseline: Lameshur Bay seagrass at same sampling time")
  #   
  #   + facet_grid(sampling_time ~ dispersion_model, space = "free_x",
  #                # + facet_wrap(sampling_time~dispersion_model, strip.position = "top",
  #                scales = "free", drop = TRUE, labeller = global_labeller)
  # )
  g2_temporal <- print(
    ggplot(data = temp_df %>%
             dplyr::filter(!grepl("site", hold2)) %>%
             droplevels, aes(y = Num_SDA, 
                             group = variable, x = variable, fill = variable))
    # group = contrast, x = Combo, fill = Combo))
    + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
               position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + theme_bw() 
    + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
            strip.text.y = element_text(angle = 0))
    + scale_fill_manual(values = variable_labels_colors, labels = variable_labels_lookup, breaks = names(variable_labels_lookup))
    # + scale_fill_manual(values = group_labels_colors, labels = group_labels_lookup, breaks = names(group_labels_lookup))
    # + scale_x_discrete(labels = group_labels_lookup)
    + geom_text(aes(label = Num_SDA, vjust = "outward", hjust = "center", 
                    group = variable, x = variable,
                    # x = Combo, group = contrast,
                    y = label_y),
                colour = "black", fontface = "bold")
    + labs(x = "Sample type", y = "Number of SDA ASVs")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             color = "none")
    + ggtitle("Baseline: Same sampling time")
    
    + facet_grid(dispersion_model ~ sampling_time, space = "free_x",
                 # + facet_wrap(sampling_time~dispersion_model, strip.position = "top",
                 scales = "free", drop = TRUE, labeller = global_labeller)
  )
  # g2_temporal 
  # g2_reef <- print(
  #   ggplot(data = temp_df %>%
  #            dplyr::filter(grepl("reef", hold2)) %>%
  #            droplevels, aes(y = Num_SDA, 
  #                            x = sampling_time, group = variable, 
  #                            fill = Combo))
  #   + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
  #   + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
  #              position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
  #   + theme_bw() 
  #   + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  #   + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
  #           strip.text.y = element_text(angle = 0))
  #   + scale_fill_manual(values = group_labels_colors, labels = group_labels_lookup, breaks = names(group_labels_lookup))
  #   + scale_x_discrete(labels = sampling_time_lookup)
  #   + geom_text(aes(label = Num_SDA, x = sampling_time, vjust = "outward", hjust = "center", group = contrast,
  #                   y = label_y),
  #               colour = "black", fontface = "bold")
  #   + labs(x = "Sample type", y = "Number of SDA ASVs")
  #   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
  #                                override.aes = list(color = "black", shape = 21, size = 3)),
  #            color = "none")
  #   + ggtitle("Between Yawzi and Tektite Reefs, same sampling time")
  #   + facet_grid(. ~ dispersion_model, space = "free_x",
  #                # + facet_wrap(sampling_time~dispersion_model, strip.position = "top",
  #                scales = "free", drop = TRUE, labeller = global_labeller)
  # )
  
  g2_spatial <- print(
    ggplot(data = temp_df %>%
             dplyr::filter(grepl("site", hold2)) %>%
             droplevels, aes(y = Num_SDA, 
                             x = site, group = contrast, 
                             fill = Combo))
    + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE,
               position = position_dodge2(padding = 0.1, preserve = "total", reverse = TRUE))
    + theme_bw() 
    + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
            strip.text.y = element_text(angle = 0))
    + scale_fill_manual(values = group_labels_colors, labels = group_labels_lookup, breaks = names(group_labels_lookup))
    + scale_x_discrete(labels = site_lookup)
    + geom_text(aes(label = Num_SDA, x = site, vjust = "outward", hjust = "center", group = contrast,
                    y = label_y),
                colour = "black", fontface = "bold")
    + labs(x = "Sample type", y = "Number of SDA ASVs")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             color = "none")
    + ggtitle("Baseline: site-specific dawn sampling time")
    # + facet_wrap(site~.,
    # + facet_grid(.~site, space = "free",
    #              scales = "free", drop = TRUE, labeller = global_labeller)
    + facet_grid(. ~ dispersion_model, space = "free_x",
                 # + facet_wrap(site~dispersion_model, strip.position = "top",
                 scales = "free", drop = TRUE, labeller = global_labeller)
  )
}
g2 <- ((g2_temporal + g2_spatial) & guides(fill = "none", color = "none") )+ patchwork::plot_annotation(title = "Number of SDA ASVs between sites at the same time, or between times in the same site",
                                                                                                        subtitle = "Used either manual dispersion estimates or automatically calculated dispersion estimates for baseline ASV abundance.")

if(!any(grepl("asvs_spatiotemporal_sda", list.files(projectpath, pattern = "usvi_deseq_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_deseq_spatiotemporal_photo_sda-", Sys.Date(), ".png"),
         g2,
         width = 10,  height = 8, units = "in")
}

#across the board, more ASVs were called significantly differentially abundant when using manual dispersions in the model
#using sampling-time specific LB Seagrass samples as baseline, also increased the number of SDA ASVs per comparison



# ASV-level deseq results -------------------------------------------------

if(any(grepl("asvs_intrasite_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_asvs_intrasite_time-.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
  usvi_deseq_asvs_intrasite_res.list <- usvi_deseq_asvs_intrasite_res.df %>%
    # dplyr::filter(grepl("manual", model)) %>%
    tidyr::drop_na(padj_qcutoff) %>%
    dplyr::filter(!(padj > 0.05)) %>%
    dplyr::select(model, contrast, asv_id, baseMean, padj, log2FoldChange_MMSE, lfcSE, padj_qcutoff) %>%
    droplevels %>%
    split(., f = .$asv_id) %>%
    map(., ~.x %>%
          dplyr::mutate(site = stringr::str_extract(contrast, "LB_seagrass_|Yawzi_|Tektite_") %>%
                          stringr::str_remove(., "_$")) %>%
          dplyr::mutate(hold = site,
                        variable = stringr::str_remove_all(contrast, paste0(hold, "_")) %>%
                          paste0("(", ., ")")) %>%
          dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
          droplevels)
  usvi_deseq_asvs_intrasite_summary.df <- usvi_deseq_asvs_intrasite_res.list %>%
    map(., ~.x %>%
          dplyr::mutate(enriched = dplyr::case_when((log2FoldChange_MMSE > 0) ~ 1,
                                                    .default = NA)) %>%
          dplyr::mutate(depleted = dplyr::case_when((log2FoldChange_MMSE < 0) ~ 1,
                                                    .default = NA)) %>%
          dplyr::group_by(baseline) %>%
          dplyr::summarise(num_enriched = sum(enriched, na.rm = TRUE),
                           num_depleted = sum(depleted, na.rm = TRUE))) %>%
    dplyr::bind_rows(., .id = "asv_id")
  
  usvi_deseq_asvs_intrasite_abund.df <- usvi_deseq_asvs_intrasite_res.list %>%
    dplyr::bind_rows(., .id = "asv_id") %>%
    dplyr::filter(grepl("manual", model)) %>%
    # droplevels %>%
    # dplyr::semi_join(., (usvi_deseq_asvs_intrasite_summary.df %>%
    #                        # dplyr::filter(if_all(c(num_enriched, num_depleted), ~.x < 2)) %>%
    #                        # dplyr::filter(if_any(c(num_enriched, num_depleted), ~gtools::odd(.x))) %>%
    #                        dplyr::filter(if_any(c(num_enriched, num_depleted), ~gtools::even(.x))) %>%
    #                        dplyr::distinct(asv_id, baseline, .keep_all = FALSE)),
    #                  by = join_by(asv_id, baseline)) %>%
    # dplyr::filter(grepl("dawn|peak", baseline)) %>%
    droplevels
  
  
  if(!any(grepl("asvs_intrasite_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
    save(usvi_deseq_asvs_intrasite_res.list, usvi_deseq_asvs_intrasite_abund.df, usvi_deseq_asvs_intrasite_res.df, 
         file = paste0(projectpath, "/", "usvi_deseq_asvs_intrasite_time-", Sys.Date(), ".RData"))
  }
}


#now do it at the ASV level between sites
if(any(grepl("asvs_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  # temp_file <- list.files(projectpath, pattern = "usvi_deseq_asvs_site_time-.*.RData")[1]
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_asvs_site_time-.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
  
  # (unique(usvi_deseq_asvs_res.df[["asv_id"]]) %in% (usvi_prok_decontam_idx %>% dplyr::filter(keep == FALSE) %>% dplyr::select(asv_id) %>% tibble::deframe(.)))
  
  usvi_deseq_asvs_res.list <- usvi_deseq_asvs_res.df %>%
    # dplyr::filter(grepl("manual", model)) %>%
    tidyr::drop_na(padj_qcutoff) %>%
    dplyr::filter(!(padj > 0.05)) %>%
    dplyr::select(model, contrast, asv_id, baseMean, padj, log2FoldChange_MMSE, lfcSE, padj_qcutoff) %>%
    droplevels %>%
    split(., f = .$asv_id) %>%
    map(., ~.x %>%
          dplyr::mutate(sampling_time = stringr::str_extract(contrast, "dawn|peak_photo")) %>%
          dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi_|Tektite_") %>%
                          stringr::str_remove(., "_")) %>%
          dplyr::mutate(Combo = stringr::str_extract(contrast, "Yawzi_.|Tektite_.")) %>%
          dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
          # dplyr::mutate(baseline = dplyr::case_when(grepl("dawn|peak", baseline) ~ baseline,
          dplyr::mutate(baseline = dplyr::case_when(grepl("LB_seagrass_", baseline) ~ "LB_seagrass_time",
                                                    grepl("Tektite", baseline) ~ "Tektite_time",
                                                    .default = "LB_seagrass_all")) %>%
          droplevels)
  
  
  #so this analysis identified genera that were significantly differentially abundant in Tektite and Yawzi
  #at dawn and peak photosynthesis sampling times
  #but these taxa were also called "significant" if they were just more noticeably more (or less) abundant in seagrass
  #and show up as enriched in both Yawzi and Tektite, or depleted in both. 
  
  #which genera were not just differentially abundant between the reef sites at different sampling times?
  #remove those genera that were num_enriched = 4 or num_depleted = 4
  #then look at those genera that had an ODD count for num_enriched or num_depleted
  
  usvi_deseq_asvs_summary.df <- usvi_deseq_asvs_res.list %>%
    map(., ~.x %>%
          dplyr::mutate(enriched = dplyr::case_when((log2FoldChange_MMSE > 0) ~ 1,
                                                    .default = NA)) %>%
          dplyr::mutate(depleted = dplyr::case_when((log2FoldChange_MMSE < 0) ~ 1,
                                                    .default = NA)) %>%
          dplyr::group_by(baseline) %>%
          dplyr::summarise(num_enriched = sum(enriched, na.rm = TRUE),
                           num_depleted = sum(depleted, na.rm = TRUE))) %>%
    dplyr::bind_rows(., .id = "asv_id")
  
  usvi_deseq_asvs_abund.df <- usvi_deseq_asvs_res.list %>%
    dplyr::bind_rows(., .id = "asv_id") %>%
    dplyr::filter(grepl("manual", model)) %>%
    droplevels %>%
    dplyr::semi_join(., (usvi_deseq_asvs_summary.df %>%
                           dplyr::filter(if_all(c(num_enriched, num_depleted), ~.x < 4)) %>%
                           dplyr::filter(if_any(c(num_enriched, num_depleted), ~gtools::odd(.x))) %>%
                           dplyr::distinct(asv_id, baseline, .keep_all = FALSE)),
                     by = join_by(asv_id, baseline)) %>%
    # dplyr::filter(grepl("dawn|peak", baseline)) %>%
    droplevels
  (usvi_deseq_asvs_abund.df %>%
      dplyr::group_by(asv_id, baseline, site) %>%
      dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
      dplyr::filter(gtools::odd(obs)))
  # 
  # usvi_deseq_asvs_abund_filtered.df <- usvi_deseq_asvs_abund.df %>%
  #   dplyr::semi_join(., (usvi_deseq_asvs_abund.df %>%
  #                          dplyr::group_by(asv_id, baseline, site) %>%
  #                          dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
  #                          dplyr::filter(gtools::odd(obs))), by = join_by(asv_id, baseline, site)) %>%
  #   droplevels
  usvi_deseq_asvs_abund_filtered.df <- usvi_deseq_asvs_abund.df %>%
    dplyr::anti_join(., (usvi_deseq_asvs_abund.df %>%
                           dplyr::distinct(contrast, asv_id, site, sampling_time, .keep_all = TRUE) %>%
                           dplyr::group_by(asv_id, site, model) %>%
                           dplyr::summarise(obs = length(log2FoldChange_MMSE)) %>%
                           dplyr::arrange(desc(obs)) %>%
                           dplyr::filter(gtools::even(obs)) %>%
                           droplevels), 
                     # relationship = "many-to-many", multiple = "all",
                     by = join_by(asv_id)) %>%
    droplevels
  
  if(!any(grepl("asvs_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
    save(usvi_deseq_asvs_res.list, usvi_deseq_asvs_abund_filtered.df, usvi_deseq_asvs_res.df, file = paste0(projectpath, "/", "usvi_deseq_asvs_site_time-", Sys.Date(), ".RData"))
  }
  
}




# Plot the SDA ASVs using manual dispersions ------------------------------


#49 ASVs in total were found to be SDA between dawn and afternoon, within sites
length(usvi_deseq_asvs_intrasite_summary.df[["asv_id"]][duplicated(usvi_deseq_asvs_intrasite_summary.df[["asv_id"]])])
#3 ASVs were temporally SDA in more than one site.


#total of 181 ASVs in usvi_deseq_asvs_abund.df:
#87 ASVs: 10 or more
#36 ASVs: 100 or more
#12 ASVs have a baseMean of 500 or more

usvi_deseq_asvs_abund.df %>%
# usvi_deseq_asvs_abund_filtered.df %>%
  dplyr::filter(!grepl("dispersion", model)) %>%
  droplevels %>%
  dplyr::distinct(asv_id, baseMean) %>%
  dplyr::ungroup(.) %>%
  dplyr::filter((baseMean >= 500)) %>%
  dplyr::distinct(asv_id) %>%
  droplevels

#total of 148 ASVs in usvi_deseq_asvs_abund_filtered.df

#so group the ASVs by rank abundance.
temp_rankabund_idx <- usvi_deseq_asvs_intrasite_abund.df %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::arrange(desc(baseMean), .by_group = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, baseMean) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  tibble::deframe(.) %>%
  quantile(., probs = seq(0.75, 0.25, -0.25), names = FALSE) %>% 
  setNames(., c("top", "middle", "bottom"))

usvi_deseq_asvs_intrasite_rankabund.df <- usvi_deseq_asvs_intrasite_abund.df %>%
  # dplyr::filter(!grepl("dispersion", model)) %>%
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
  dplyr::right_join(., usvi_deseq_asvs_intrasite_abund.df, by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
  droplevels


g3 <- print(
  ggplot(data = usvi_deseq_asvs_intrasite_rankabund.df %>%
           droplevels)
  + theme_bw() 
  + geom_hline(aes(yintercept = 0), color = "black")
  + geom_segment(aes(y = log2FoldChange_MMSE, yend = 0, x = asv_id, group = asv_id, 
                     color = site), 
                 position = position_dodge2(width = 0.2, preserve = "total"), 
                 # position = position_jitter(width = 0.2, seed = 48105), 
                 show.legend = FALSE)
  + geom_errorbar(aes(x = asv_id, ymin = (log2FoldChange_MMSE - lfcSE), ymax = (log2FoldChange_MMSE + lfcSE),
                      group = asv_id, color = site), 
                  position = position_dodge2(width = 0.2, preserve = "total"),
                  show.legend = FALSE)
  + geom_point(aes(y = log2FoldChange_MMSE, x = asv_id, group = asv_id, fill = site, shape = site),
               position = position_dodge2(width = 0.2, preserve = "total"),
               # position = position_jitter(width = 0.2, seed = 48105),
               alpha = 0.8, color = "black", size = 3, show.legend = TRUE)
  + scale_shape_manual(name = "Sampling site and time",
                       values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  + scale_fill_manual(name = "Sampling site and time",
                      values = site_colors, labels = site_lookup, breaks = names(site_lookup))
  + scale_color_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
  + scale_x_discrete(labels = usvi_genera_relabel)
  + scale_y_continuous(name = "log2(fold change) from dawn")
  + guides(color = "none",
           fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
                               override.aes = list(color = "black", size = 3)),
           shape = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
                                override.aes = list(color = "black", size = 3)))
  + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
          legend.position = "right", 
          strip.text.y = element_text(angle = 0))
  + facet_grid(rows = NULL, cols = vars(quantile_group),  space = "free",
               drop = TRUE, scales = "free", labeller = global_labeller)
  + ggtitle("Significantly differentially abundant from dawn sampling time")
)


temp_rankabund_idx <- usvi_deseq_asvs_abund.df %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::arrange(desc(baseMean), .by_group = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(asv_id, baseMean) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  tibble::deframe(.) %>%
  quantile(., probs = seq(0.75, 0.25, -0.25), names = FALSE) %>% 
  setNames(., c("top", "middle", "bottom"))

usvi_deseq_asvs_abund_rankabund.df <- usvi_deseq_asvs_abund.df %>%
  # dplyr::filter(!grepl("dispersion", model)) %>%
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
  dplyr::right_join(., usvi_deseq_asvs_abund_filtered.df, by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
  dplyr::filter(!grepl("dispersion", model)) %>%
  dplyr::mutate(hold = stringr::str_remove_all(model, "manual_|auto_")) %>%
  dplyr::mutate(hold = dplyr::case_when(grepl("lameshur", hold) ~ "LB_seagrass",
                                        grepl("photo", hold) ~ "peak_photo",
                                        grepl("yawzi|tektite|lb", hold) ~ stringr::str_to_sentence(hold),
                                        grepl("reef", model) ~ sampling_time,
                                        .default = hold)) %>%
  dplyr::mutate(variable = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ stringr::str_remove_all(contrast, paste0("_", hold)),
                                            (hold %in% names(site_lookup)) ~ stringr::str_remove_all(contrast, paste0(site, "_")),
                                            .default = NA)) %>%
  dplyr::mutate(hold2 = dplyr::case_when(grepl("reef", model) ~ "interreef",
                                         grepl("dawn|photo", hold) ~ "sampling time",
                                         .default = "site")) %>%
  dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
  dplyr::mutate(baseline2 = stringr::str_remove_all(baseline, "_dawn|_peak_photo")) %>%
  droplevels

g4 <- print(
  ggplot(data = usvi_deseq_asvs_abund_rankabund.df %>%
           dplyr::filter(!grepl("dispersion", model)) %>%
           droplevels)
  + theme_bw() 
  + geom_hline(aes(yintercept = 0), color = "black")
  + geom_segment(aes(y = log2FoldChange_MMSE, yend = 0, x = asv_id, group = asv_id, 
                     color = sampling_time), 
                 position = position_dodge2(width = 0.2, preserve = "total"),
                 # position = position_jitter(width = 0.2, seed = 48105), 
                 show.legend = FALSE)
  + geom_errorbar(aes(x = asv_id, ymin = (log2FoldChange_MMSE - lfcSE), ymax = (log2FoldChange_MMSE + lfcSE),
                      group = asv_id, color = variable), 
                  position = position_dodge2(width = 0.2, preserve = "total"),
                  show.legend = FALSE)
  + geom_point(aes(y = log2FoldChange_MMSE, x = asv_id, group = asv_id, fill = sampling_time, shape = sampling_time),
               alpha = 0.8, color = "black", size = 3,
               # position = position_jitter(width = 0.2, seed = 48105)
               position = position_dodge2(width = 0.2, preserve = "total"), 
               show.legend = TRUE)
  + scale_shape_manual(name = "Sampling time",
                       values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
  + scale_fill_manual(name = "Sampling time",
                      values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_color_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_x_discrete(labels = usvi_genera_relabel)
  + scale_y_continuous(name = "log2(fold change), hold sampling time")
  + guides(color = "none")
  + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
          strip.text.y = element_text(angle = 0))
  # + facet_wrap(. ~ quantile_group, nrow = 1,
  + facet_grid(rows = vars(variable), cols = vars(quantile_group),  space = "fixed",
               drop = TRUE, scales = "free_x", labeller = global_labeller)
  + ggtitle("Significantly differentially abundant between reefs and seagrass")
)
# g4 <- print(
#   ggplot(data = usvi_deseq_asvs_abund_rankabund.df %>%
#            dplyr::filter(!grepl("dispersion", model)) %>%
#            droplevels)
#   + theme_bw() 
#   + geom_hline(aes(yintercept = 0), color = "black")
#   + geom_segment(aes(y = log2FoldChange_MMSE, yend = 0, x = asv_id, group = asv_id, 
#                      color = sampling_time), 
#                  position = position_dodge2(width = 0.2, preserve = "total"),
#                  # position = position_jitter(width = 0.2, seed = 48105), 
#                  show.legend = FALSE)
#   + geom_errorbar(aes(x = asv_id, ymin = (log2FoldChange_MMSE - lfcSE), ymax = (log2FoldChange_MMSE + lfcSE),
#                       group = asv_id, color = sampling_time), 
#                   position = position_dodge2(width = 0.2, preserve = "total"),
#                   show.legend = FALSE)
#   + geom_point(aes(y = log2FoldChange_MMSE, x = asv_id, group = asv_id, fill = sampling_time, shape = sampling_time),
#                alpha = 0.8, color = "black", size = 3,
#                # position = position_jitter(width = 0.2, seed = 48105)
#                position = position_dodge2(width = 0.2, preserve = "total"), 
#                show.legend = TRUE)
#   + scale_shape_manual(name = "Sampling time",
#                        values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
#   + scale_fill_manual(name = "Sampling time",
#                       values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_color_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_x_discrete(labels = usvi_genera_relabel)
#   + scale_y_continuous(name = "log2(fold change) from contemporary LB seagrass")
#   + guides(color = "none")
#   + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
#           strip.text.y = element_text(angle = 0))
#   # + facet_wrap(. ~ quantile_group, nrow = 1,
#   + facet_grid(rows = vars(site), cols = vars(quantile_group),  space = "free",
#                drop = TRUE, scales = "free", labeller = global_labeller)
#   + ggtitle("Significantly differentially abundant between reefs and seagrass")
# )




if(!any(grepl("asvs_intrasite_temporal_sda", list.files(projectpath, pattern = "usvi_deseq_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_deseq_asvs_intrasite_temporal_sda-", Sys.Date(), ".png"),
         g3,
         width = 10,  height = 8, units = "in")
}

if(!any(grepl("asvs_intratemporal_sda", list.files(projectpath, pattern = "usvi_deseq_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_deseq_asvs_intratemporal_sda-", Sys.Date(), ".png"),
         g4,
         width = 10,  height = 8, units = "in")
}


