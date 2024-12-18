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

no_progress <- function() {} #your placeholder function for progress reporting

idx_samples <- function(x) grep("samples", names(x), value = TRUE)
coalesce2 <- function(x, y, sep = ".") ifelse(x == y, coalesce(x, y, sep = sep), paste0(x, "_vs_", y))

set.alpha <- 0.05
set.seed(48105)

t_pseudolog10 <- scales::new_transform("pseudolog10", 
                                       function(x)( log10(x+1)), 
                                       function(x)( 10^(x) - 1), 
                                       domain = c(0, Inf))

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


#let's try consolidating the ASV table to species- or genus-level first
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
  droplevels

usvi_sw_genus.tbl <- usvi_sw_genus.tbl %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  # apply(., 2, relabund) %>%
  # as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.)

# usvi_sw_genus.mat <- usvi_sw_genus.tbl %>%
#   t(.)
# 
# usvi_sw_genus.df <- usvi_sw_genus.mat %>%
usvi_sw_genus.df <- usvi_sw_genus.tbl %>%
  t(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  tidyr::pivot_longer(., cols = !c(sample_id),
                      names_to = "asv_id",
                      values_to = "counts")


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




# Read in results from DESeq2 ---------------------------------------------

if(!any(grepl("asvs_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_asvs_site_time.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
}
#usvi_deseq_asvs_res.list, usvi_deseq_asvs_abund_filtered.df, usvi_deseq_asvs_res.df

if(!any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_deseq_.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_deseq_genera_site_time.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
}
#usvi_deseq_genera_res.list, usvi_deseq_genera_abund_filtered.df, usvi_deseq_genera_abund.df, usvi_deseq_genera_res.df


# Read in results from radEmu ---------------------------------------------



#read in the genera that were significant by site through rademu:
if(any(grepl("genera_reef_highlow", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
  # temp_file <- list.files(projectpath, pattern = "usvi_rademu_genera_reef_highlow.*.RData")[1]
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_rademu_genera_reef_highlow.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
}

# rademu_reef_res.df <- score_lowreef_res_list %>%
#   bind_rows(.) %>%
#   dplyr::mutate(group = "low_reef") %>%
#   bind_rows(., score_highreef_res_list %>%
#               bind_rows(.) %>%
#               dplyr::mutate(group = "high_reef")) %>%
#   dplyr::rename(asv_id = "category") %>%
#   droplevels

if(any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
  # temp_file <- list.files(projectpath, pattern = "usvi_rademu_genera_site_time.*.RData")[1]
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_rademu_genera_site_time.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
}

rademu_res.df <- bind_rows((score_1_res_list %>%
                              bind_rows(.) %>% #these were genera that were found to be higher in samples collected at dawn
                              dplyr::mutate(group = "high_dawn") %>%
                              dplyr::filter(pval < 0.05) %>%
                              tidyr::drop_na(score_stat)),
                           (score_2_res_list %>%
                              bind_rows(.) %>% # these were genera that were found to be higher in samples collected at peak photo
                              dplyr::mutate(group = "low_dawn") %>%
                              dplyr::filter(pval < 0.05) %>%
                              tidyr::drop_na(score_stat)),
                           (score_4_res_list %>%
                              bind_rows(.) %>% #these were genera that were found to be higher in samples collected in seagrass
                              dplyr::mutate(group = "high_seagrass") %>%
                              dplyr::filter(pval < 0.05) %>%
                              tidyr::drop_na(score_stat)),
                           (score_5_res_list %>%
                              bind_rows(.) %>% #these were genera that were found to be higher in samples collected in reef sites
                              dplyr::mutate(group = "low_seagrass") %>%
                              dplyr::filter(pval < 0.05) %>%
                              tidyr::drop_na(score_stat)),
                           (score_lowreef_res_list %>%
                              bind_rows(.) %>%
                              dplyr::mutate(group = "low_reef") %>%
                              dplyr::filter(pval < 0.05) %>%
                              tidyr::drop_na(score_stat)),
                           (score_highreef_res_list %>%
                              bind_rows(.) %>%
                              dplyr::mutate(group = "high_reef") %>%
                              dplyr::filter(pval < 0.05) %>%
                              tidyr::drop_na(score_stat))) %>%
  dplyr::rename(asv_id = "category") %>%
  droplevels %>%
  dplyr::rename(contrast = "covariate") %>%
  dplyr::mutate(contrast = stringr::str_remove_all(contrast, "site") %>%
                  stringr::str_replace_all(., ":sampling_time", "_")) %>%
  dplyr::mutate(contrast = dplyr::case_when(grepl("reef", group) ~ paste0(contrast, " - LB_seagrass_all"),
                                            grepl("dawn", group) ~ paste0(contrast, " - all_peak_photo"),
                                            grepl("seagrass", group) ~ paste0(contrast, " - LB_seagrass_all"),
                                            .default = contrast)) %>%
  dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
  dplyr::mutate(commonid = stringr::str_split_i(group, "_", 2)) %>%
  dplyr::arrange(commonid, asv_id) %>%
  dplyr::group_by(commonid) %>%
  dplyr::distinct(asv_id, group, contrast, .keep_all = TRUE) %>%
  droplevels

#for the temporal nuances:
#subtract all genera that were just significantly higher/lower in seagrass regardless of time
# temp_list <- rademu_res.df %>%
#   split(., f = .$commonid) %>%
#   map(., ~.x %>%
#         droplevels)
rademu_seagrass_res.df <- rademu_res.df %>%
  dplyr::filter(commonid == "seagrass") %>%
  dplyr::distinct(asv_id, group, .keep_all = TRUE) %>%
  dplyr::arrange(asv_id) %>%
  droplevels

rademu_reef_res.df <- rademu_res.df %>%
  dplyr::filter(commonid == "reef") %>%
  dplyr::distinct(asv_id, group, .keep_all = TRUE) %>%
  dplyr::arrange(asv_id) %>%
  droplevels

rademu_time_res.df <- rademu_res.df %>%
  dplyr::filter(commonid == "dawn") %>%
  dplyr::distinct(asv_id, group, .keep_all = TRUE) %>%
  dplyr::arrange(asv_id) %>%
  droplevels

#these are the genera that were shared between high/low reef and high/low seagrass (23):
rademu_site_genera_idx <- intersect(rademu_reef_res.df[["asv_id"]], rademu_seagrass_res.df[["asv_id"]])
# union(rademu_reef_res.df[["asv_id"]], rademu_seagrass_res.df[["asv_id"]])

# #these are the genera that were not shared between high/low reef and high/low seagrass (13):
# setdiff(rademu_reef_res.df[["asv_id"]], rademu_seagrass_res.df[["asv_id"]])



#now find the difference between the genera flagged as spatiotemporally different, and the site different:
# setdiff(rademu_reef_res.df[["asv_id"]], rademu_time_res.df[["asv_id"]])
setdiff(rademu_time_res.df[["asv_id"]], rademu_site_genera_idx)
# [1] "ASV_00008" "ASV_00012" "ASV_00215" "ASV_00227" "ASV_00243" "ASV_00282" "ASV_01134"



# Read in results from corncob --------------------------------------------

#read in the genera that were significant via corncob:
# if(any(grepl("cc_dt_sda", list.files(projectpath, pattern = ".*ps_usvi_filtered-.*.RData")))){
#   temp_file <- data.table::last(list.files(projectpath, pattern = "cc_dt_sda_ps_usvi_filtered*.RData"))
#   load(paste0(projectpath, "/", temp_file))
#   rm(temp_file)
# }

if(file.exists(paste0(projectpath, "/", "cc_dt_usvi_summary.df", ".tsv"))){
  cc_dt_usvi_summary.df <- readr::read_delim(paste0(projectpath, "/", "cc_dt_usvi_summary.df", ".tsv"),
                                             delim = "\t", col_names = TRUE)
}


#since rademu couldn't proceed at the ASV level, use just corncob results at the agglomerated genera-level 

corncob_res.df <- cc_dt_usvi_summary.df %>%
  dplyr::filter(grepl("agg_genus", taxon_resolution)) %>%
  dplyr::mutate(test_type = "corncob") %>%
  dplyr::rename(padj = "p_adj") %>%
  dplyr::rename(model = "test") %>%
  dplyr::select(contrast, hold, variable, asv_id, p_fdr, padj, taxonomy, test_type, model) %>%
  dplyr::distinct(asv_id, hold, variable, .keep_all = TRUE) %>%
  # tidyr::separate_longer_delim(variable, delim = " - ") %>%
  dplyr::mutate(pair1 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\1", variable), "_", hold),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\1", variable)))) %>%
  dplyr::mutate(pair2 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\3", variable), "_", hold),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)))) %>%
  # dplyr::mutate(contrast = paste0(pair1, " - ", pair2)) %>%
  dplyr::select(contrast, hold, variable, asv_id, padj, p_fdr, test_type, model) %>%
  # dplyr::mutate(sampling_time = dplyr::case_when(grepl("dawn|peak", contrast) ~ stringr::str_extract(contrast, "dawn|peak_photo"),
  #                                                .default = "all")) %>%
  # dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
  # dplyr::mutate(Combo = dplyr::case_when(grepl("Yawzi_|Tektite_", contrast) ~ stringr::str_extract(contrast, "Yawzi_.|Tektite_."),
  #                                        .default = site)) %>%
  # dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
  # dplyr::mutate(baseline = dplyr::case_when(grepl("all_peak_photo", baseline) ~ "all_peak_photo",
  #                                           grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
  #                                           .default = "LB_seagrass_all")) %>%
  droplevels


usvi_sda_genera_compare.df <- rademu_res.df %>%
  dplyr::mutate(test_type = "rademu") %>%
  dplyr::rename(estimate = "score_stat") %>%
  dplyr::mutate(model = "auto_fit") %>%
  dplyr::mutate(sampling_time = dplyr::case_when(grepl("dawn|peak", contrast) ~ stringr::str_extract(contrast, "dawn|peak_photo"),
                                                 .default = "all")) %>%
  dplyr::mutate(site = stringr::str_extract(contrast, "Yawzi|Tektite")) %>%
  dplyr::mutate(Combo = dplyr::case_when(grepl("Yawzi_|Tektite_", contrast) ~ stringr::str_extract(contrast, "Yawzi_.|Tektite_."),
                                         .default = site)) %>%
  dplyr::mutate(baseline = gsub("([[:print:]]+)( - )(.*)$", "\\3", contrast)) %>%
  dplyr::mutate(baseline = dplyr::case_when(grepl("all_peak_photo", baseline) ~ "all_peak_photo",
                                            grepl("dawn|peak", baseline) ~ "LB_seagrass_time",
                                            .default = "LB_seagrass_all")) %>%
  bind_rows(., (usvi_deseq_genera_abund.df %>%
                  # bind_rows(., (usvi_deseq_genera_abund_filtered.df %>%
                  dplyr::select(-padj_qcutoff) %>%
                  dplyr::mutate(test_type = "deseq"))) %>%
  dplyr::mutate(pair1 = stringr::str_split_i(contrast, " - ", 1),
                pair2 = stringr::str_split_i(contrast, " - ", 2)) %>%
  bind_rows(., corncob_res.df)



#which genera are flagged as significantly differentially abundant in all tests?
shared_sda_genera_idx <- usvi_sda_genera_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(test_type)) %>%
  # dplyr::filter(num_results > 1) %>%
  dplyr::filter(num_results > 2) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id) %>%
  unlist %>%
  as.character
# droplevels


#which are not?
selfish_sda_genera_idx <- usvi_sda_genera_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(test_type)) %>%
  dplyr::filter(num_results <= 1) %>%
  # dplyr::filter(num_results <= 2) %>%
  dplyr::distinct(asv_id) %>%
  dplyr::left_join(., usvi_sda_genera_compare.df %>%
                     dplyr::distinct(asv_id, test_type, .keep_all = FALSE),
                   by = join_by(asv_id)) %>%
  dplyr::select(test_type, asv_id) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  tibble::deframe(.)
# unlist %>%
# as.character

usvi_sda_genera_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, test_type, pval, padj, p_fdr, .keep_all = FALSE) %>%
  dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(asv_id, test_type, p_value) %>%
  dplyr::distinct(asv_id, test_type, .keep_all = TRUE) %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "test_type",
                     values_from = "p_value") %>%
  dplyr::mutate(across(c("rademu", "deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
                                                                          .default = NA))) %>%
  tidyr::unite(test_type, c("rademu", "deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
  droplevels

semiself_sda_genera_idx <- usvi_sda_genera_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(test_type)) %>%
  # dplyr::filter(num_results <= 1) %>%
  dplyr::filter(num_results <= 2) %>%
  dplyr::distinct(asv_id) %>%
  dplyr::anti_join(., tibble::enframe(selfish_sda_genera_idx, name = "test_type", value = "asv_id")) %>%
  dplyr::left_join(., (usvi_sda_genera_compare.df %>%
                         dplyr::ungroup(.) %>%
                         dplyr::distinct(asv_id, test_type, pval, padj, p_fdr, .keep_all = FALSE) %>%
                         dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
                         dplyr::select(asv_id, test_type, p_value) %>%
                         dplyr::distinct(asv_id, test_type, .keep_all = TRUE) %>%
                         tidyr::pivot_wider(., id_cols = "asv_id",
                                            names_from = "test_type",
                                            values_from = "p_value") %>%
                         dplyr::mutate(across(c("rademu", "deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
                                                                                                 .default = NA))) %>%
                         tidyr::unite(test_type, c("rademu", "deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
                         droplevels),
                   # dplyr::left_join(., usvi_sda_genera_compare.df %>%
                   #                    dplyr::select(asv_id, test_type) %>%
                   #                    tidyr::pivot_wider(., id_cols = "asv_id",
                   #                                       names_from = "metric",
                   #                                       values_from = "test_type") %>%
                   #                    dplyr::distinct(asv_id, test_type, .keep_all = FALSE),
                   by = join_by(asv_id)) %>%
  dplyr::select(test_type, asv_id) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  tibble::deframe(.)

#save the results as a table listing which genera are significant in which comparison/test

usvi_sda_genera_compare_summary.df <- usvi_sda_genera_compare.df %>%
  dplyr::distinct(asv_id, contrast, test_type) %>%
  dplyr::arrange(asv_id) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = c(contrast),
                     values_from = c(test_type),
                     values_fill = NA) %>%
  droplevels %>%
  dplyr::left_join(., selfish_sda_genera_idx %>%
                     tibble::enframe(name = "test_type.x", value = "asv_id")) %>%
  dplyr::left_join(., semiself_sda_genera_idx %>%
                     tibble::enframe(name = "test_type.z", value = "asv_id")) %>%
  dplyr::left_join(., data.frame(asv_id = shared_sda_genera_idx,
                                 test_type.y = "all")) %>%
  tidyr::unite("test_type", c("test_type.x", "test_type.z", "test_type.y"), remove = TRUE, na.rm = TRUE) %>%
  dplyr::left_join(., usvi_sw_genus.taxa.df, by = join_by(asv_id)) %>%
  droplevels

readr::write_delim(usvi_sda_genera_compare_summary.df, paste0(projectpath, "/", "usvi_sda_genera_compare_summary.tsv"),
                   delim = "\t", col_names = TRUE, na = "")
if(!any(grepl("rds$", list.files(projectpath, pattern = "usvi_sda_genera_compare.*.rds")))){
  readr::write_rds(usvi_sda_genera_compare.df, 
                   paste0(projectpath, "/", "usvi_sda_genera_compare.rds"))
  #also save the results with the pvals and estimates
  readr::write_delim(usvi_sda_genera_compare.df, 
                     paste0(projectpath, "/", "usvi_sda_genera_compare.tsv"),
                     delim = "\t", col_names = TRUE, na = "")
}



# Plot the genera that were SDA across all 3 methods ----------------------

shared_sda_genera_idx_cluster <- usvi_sda_genera_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_genera_idx) %>%
  dplyr::distinct(contrast, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(contrast)) %>%
  dplyr::arrange(num_results) %>%
  dplyr::mutate(grouping = c(rep(1:3, each = 5, length.out = length(shared_sda_genera_idx)))) %>%
  dplyr::arrange(asv_id) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        dplyr::select(asv_id, num_results) %>%
        droplevels)

rm(list = apropos("g_sda_.*", mode = "list"))

for(i in seq_along(shared_sda_genera_idx_cluster)){
  namevar <- shared_sda_genera_idx_cluster[[i]] %>%
    purrr::pluck("asv_id")
  temp_df2 <- usvi_sda_genera_compare.df %>%
    dplyr::filter(asv_id %in% namevar) %>%
    droplevels
  g <- print(
    ggplot(data = temp_df2)
    + theme_bw() 
    + geom_boxplot(aes(y = relabund, x = sampling_time, group = sampling_time), 
                   alpha = 0.6, show.legend = FALSE, color = "black",
                   position = position_dodge(0.8))
    + geom_point(aes(y = relabund, x = sampling_time, group = sampling_time, fill = sampling_time, shape = site),
                 alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
                 position = position_jitter(width = 0.2, seed = 48105))
    + scale_shape_manual(name = "Sampling site and time",
                         values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
    + scale_fill_manual(name = "Sampling time",
                        values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
    + scale_x_discrete(labels = sampling_time_lookup, name = "Sampling time")
    + scale_y_continuous(name = "Relative abundance (%)")
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             color = "none")
    # + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
    + theme(axis.text.x = element_blank(),
            strip.text.y = element_text(size = rel(0.7)),
            axis.title.y = element_blank())
    + facet_grid(asv_id ~ site,
                 drop = TRUE,
                 # labeller = labeller(site = site_lookup, asv_id = stringr::str_wrap(usvi_genera_relabel)),
                 labeller = global_labeller,
                 scales = "free")
    + ggtitle(paste0("Group ", i))
  )
  assign(paste0("g_sda_", i), g, envir = .GlobalEnv, inherits = TRUE)
  rm(g)
  rm(namevar)
}

gpatch <- apropos("g_sda_.*", mode = "list") %>%
  lapply(., get) %>%
  # purrr::reduce(., `+ theme(legend.position = "none")`) %>%
  purrr::reduce(., `+`)
g2_sda_genera <- gpatch + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")
g2_sda_genera

if(!any(grepl("genera_relabund", list.files(projectpath, pattern = "usvi_sda_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_sda_genera_relabund-", Sys.Date(), ".png"),
         g2_sda_genera,
         width = 15, height = 10, units = "in")
}


# Plot the relative abundances of ASVs found sig --------------------------


corncob_asv_res.df <- cc_dt_usvi_summary.df %>%
  dplyr::filter(grepl("asv", taxon_resolution)) %>%
  dplyr::mutate(test_type = "corncob") %>%
  dplyr::rename(padj = "p_adj") %>%
  dplyr::rename(model = "test") %>%
  dplyr::select(contrast, hold, variable, asv_id, p_fdr, padj, taxonomy, test_type, model) %>%
  dplyr::distinct(asv_id, hold, variable, .keep_all = TRUE) %>%
  dplyr::mutate(pair1 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\1", variable), "_", hold),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\1", variable)))) %>%
  dplyr::mutate(pair2 = dplyr::case_when((hold %in% names(sampling_time_lookup)) ~ paste0(gsub("([[:print:]]+)( - )(.*)$", "\\3", variable), "_", hold),
                                         .default = paste0(hold, "_", gsub("([[:print:]]+)( - )(.*)$", "\\3", variable)))) %>%
  dplyr::select(contrast, hold, variable, asv_id, padj, p_fdr, test_type, model) %>%
  droplevels


usvi_sda_asvs_compare.df <- bind_rows(
  (usvi_deseq_asvs_abund_filtered.df %>%
     # (usvi_deseq_asvs_res.df %>%
     dplyr::select(-padj_qcutoff) %>%
     dplyr::mutate(test_type = "deseq"))) %>%
  dplyr::mutate(pair1 = stringr::str_split_i(contrast, " - ", 1),
                pair2 = stringr::str_split_i(contrast, " - ", 2)) %>%
  bind_rows(., corncob_asv_res.df)



#which genera are flagged as significantly differentially abundant in all tests?
shared_sda_asvs_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(test_type)) %>%
  dplyr::filter(num_results > 1) %>%
  # dplyr::filter(num_results > 2) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id) %>%
  unlist %>%
  as.character
# droplevels


#which are not?
selfish_sda_asvs_idx <- usvi_sda_asvs_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(test_type, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(test_type)) %>%
  dplyr::filter(num_results <= 1) %>%
  # dplyr::filter(num_results <= 2) %>%
  dplyr::distinct(asv_id) %>%
  dplyr::left_join(., usvi_sda_asv_compare.df %>%
                     dplyr::distinct(asv_id, test_type, .keep_all = FALSE),
                   by = join_by(asv_id)) %>%
  dplyr::select(test_type, asv_id) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  tibble::deframe(.)

usvi_sda_asvs_compare.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, test_type, padj, p_fdr, .keep_all = FALSE) %>%
  dplyr::mutate(p_value = across(starts_with("p")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(asv_id, test_type, p_value) %>%
  dplyr::distinct(asv_id, test_type, .keep_all = TRUE) %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "test_type",
                     values_from = "p_value") %>%
  dplyr::mutate(across(c("deseq", "corncob"), ~dplyr::case_when(!is.na(.x) ~ deparse(substitute(.x)),
                                                                .default = NA))) %>%
  tidyr::unite(test_type, c("deseq", "corncob"), sep = "_", na.rm = TRUE, remove = TRUE) %>%
  droplevels

#save the results as a table listing which genera are significant in which comparison/test

usvi_sda_asvs_compare_summary.df <- usvi_sda_asvs_compare.df %>%
  dplyr::distinct(asv_id, contrast, test_type) %>%
  dplyr::arrange(asv_id) %>%
  droplevels %>%
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
  dplyr::left_join(., usvi_prok_asvs.taxonomy, by = join_by(asv_id)) %>%
  droplevels

readr::write_delim(usvi_sda_asvs_compare_summary.df, paste0(projectpath, "/", "usvi_sda_asvs_compare_summary.tsv"),
                   delim = "\t", col_names = TRUE, na = "")




usvi_sda_asvs_relabund.df <- usvi_prok_asvs.df %>%
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

shared_sda_asvs_idx_cluster <- usvi_sda_asvs_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_idx) %>%
  dplyr::distinct(contrast, asv_id, .keep_all = FALSE) %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(num_results = length(contrast)) %>%
  dplyr::arrange(num_results) %>%
  dplyr::mutate(grouping = c(rep(1:5, each = 13, length.out = length(shared_sda_asvs_idx)))) %>%
  dplyr::arrange(asv_id) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        dplyr::select(asv_id, num_results) %>%
        droplevels)

rm(list = apropos("g_sda_.*", mode = "list"))

for(i in seq_along(shared_sda_asvs_idx_cluster)){
  namevar <- shared_sda_asvs_idx_cluster[[i]] %>%
    purrr::pluck("asv_id")
  temp_df2 <- usvi_sda_asvs_relabund.df %>%
    dplyr::filter(asv_id %in% namevar) %>%
    droplevels
  g <- print(
    ggplot(data = temp_df2)
    + theme_bw() 
    # + geom_boxplot(aes(y = relabund, x = sampling_time, group = sampling_time), 
    #                alpha = 0.6, show.legend = FALSE, color = "black", outliers = TRUE, outlier.shape = NA,
    #                position = position_dodge(0.8))
    + geom_violin(draw_quantiles = c(0.5), trim = TRUE, scale = "area",
                  aes(y = relabund, x = sampling_time, group = sampling_time), 
                  alpha = 0.6, show.legend = FALSE, color = "grey", 
                  position = position_dodge(0.8))
    + geom_point(aes(y = relabund, x = sampling_time, group = sampling_time, fill = sampling_time, shape = site),
                 alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
                 position = position_jitter(width = 0.2, seed = 48105))
    + scale_shape_manual(name = "Sampling site and time",
                         values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
    + scale_fill_manual(name = "Sampling time",
                        values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
    + scale_x_discrete(labels = sampling_time_lookup, name = "Sampling time")
    + scale_y_continuous(name = "Relative abundance (%)", expand = expansion(mult = c(0.1,0.1)))
    + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
                                 override.aes = list(color = "black", shape = 21, size = 3)),
             color = "none")
    # + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
    + theme(axis.text.x = element_blank(),
            strip.text.y = element_text(size = rel(0.7)),
            axis.title.y = element_blank())
    + facet_grid(asv_id ~ site,
                 drop = TRUE,
                 # labeller = labeller(site = site_lookup, asv_id = stringr::str_wrap(usvi_genera_relabel)),
                 labeller = global_labeller,
                 scales = "free")
    + ggtitle(paste0("Group ", i))
  )
  assign(paste0("g_sda_", i), g, envir = .GlobalEnv, inherits = TRUE)
  rm(g)
  rm(namevar)
}

gpatch <- apropos("g_sda_.*", mode = "list") %>%
  lapply(., get) %>%
  # purrr::reduce(., `+ theme(legend.position = "none")`) %>%
  purrr::reduce(., `+`)
g2_sda_asvs <- gpatch + patchwork::plot_layout(guides = "collect")  & theme(legend.position = "none")
# g2_sda_asvs

if(!any(grepl("asvs_relabund", list.files(projectpath, pattern = "usvi_sda_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund-", Sys.Date(), ".png"),
         g2_sda_asvs,
         width = 15, height = 24, units = "in")
}


for(i in seq_along(shared_sda_asvs_idx_cluster)){
  namevar <- paste0("g_sda_", i)
  temp_g2_sda_asvs <- get(namevar, inherits = TRUE)
  temp_g2_sda_asvs <- temp_g2_sda_asvs + theme(legend.position = "none")
  ggsave(paste0(projectpath, "/", "usvi_sda_asvs_relabund_group", i, "-", Sys.Date(), ".png"),
         temp_g2_sda_asvs,
         width = 8, height = 16, units = "in")
}



# How many SDA ASVs belonged to genera that were also SDA? ----------------

usvi_sda_asvs_in_sda_genera.df <- dplyr::full_join((usvi_sda_genera_compare_summary.df %>%
                                                      dplyr::ungroup(.) %>%
                                                      dplyr::distinct(asv_id, test_type) %>%
                                                      # dplyr::mutate(test_type = dplyr::case_when(grepl("all", test_type) ~ 2,
                                                      #                                                grepl("_", test_type) ~ 1.5,
                                                      #                                                .default = 1)) %>%
                                                      dplyr::mutate(test_type = 1) %>%
                                                      dplyr::rename(sda_agg_genus = "test_type")),
                                                   (usvi_sda_asvs_compare_summary.df %>%
                                                      dplyr::ungroup(.) %>%
                                                      dplyr::distinct(asv_id, test_type) %>%
                                                      dplyr::mutate(test_type = 1) %>%
                                                      # dplyr::mutate(test_type = dplyr::case_when(grepl("all", test_type) ~ 2,
                                                      #                                                .default = 1)) %>%
                                                      dplyr::rename(sda_asv = "test_type"))) %>%
  dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"), by = join_by(asv_id)) %>%
  droplevels %>%
  dplyr::mutate(taxonomy = gsub("^([[:alnum:]_]{9}\\:\\s)(.*)", "\\2", taxonomy))

usvi_sda_asvs_in_sda_genera_summary.df <- usvi_sda_asvs_in_sda_genera.df %>%
  # tidyr::pivot_longer(., cols = !c("asv_id", "taxonomy"),
  #                     names_to = "taxon_resolution",
  #                     values_to = "test_type") %>%
  # dplyr::distinct(asv_id, taxonomy, taxon_resolution, .keep_all = TRUE) %>%
  # dplyr::group_by(taxonomy, taxon_resolution) %>%
  dplyr::group_by(taxonomy) %>%
  dplyr::summarise(num_sig_asvs = sum(sda_asv, na.rm = TRUE),
                   num_sig_genera = sum(sda_agg_genus, na.rm = TRUE)) %>%
  droplevels

usvi_sda_resolution_summary.df <- usvi_prok_asvs.taxonomy %>%
  dplyr::group_by(taxonomy) %>%
  dplyr::summarise(num_agglom_asvs = length(asv_id)) %>%
  dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
  dplyr::full_join(., usvi_sda_asvs_in_sda_genera_summary.df,
                   by = join_by(taxonomy)) %>%
  dplyr::full_join(., usvi_sda_asvs_in_sda_genera_summary.df %>% 
                     dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
                     dplyr::select(taxonomy, num_sig_asvs) %>%
                     dplyr::rename(num_sig_asvs_not_genera = num_sig_asvs),
                   by = join_by(taxonomy)) %>%
  # dplyr::full_join(., usvi_sda_asvs_in_sda_genera_summary.df %>%
  #                    dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
  #                    dplyr::select(taxonomy, num_sig_genera) %>%
  #                    dplyr::rename(num_sig_genera_not_asvs = num_sig_genera),
  #                  by = join_by(taxonomy)) %>%
  tidyr::drop_na(num_sig_asvs, num_sig_genera) %>%
  # dplyr::mutate(num_not_sig_asvs = dplyr::case_when(is.na(num_sig_asvs_not_genera) ~ num_agglom_asvs - num_sig_asvs,
  dplyr::mutate(num_not_sig_asvs = dplyr::case_when(!is.na(num_sig_asvs) ~ num_agglom_asvs - num_sig_asvs,
                                                    .default = NA)) %>%
  dplyr::select(-c(num_sig_asvs_not_genera)) %>%
  droplevels


#how many genera were found to be significant at either the ASV or genus-level?
#170 genera:
usvi_sda_resolution_summary.df %>% 
  dplyr::distinct(taxonomy) %>%
  dplyr::summarise(num_sig_genera_overall = length(taxonomy))

#how many total ASVs are in these agglomerated genera?
#8061 ASVs:
usvi_sda_resolution_summary.df %>%
  # dplyr::filter(num_sig_genera > 0) %>%
  dplyr::ungroup(.) %>%
  dplyr::summarise(num_agglom_asvs = sum(num_agglom_asvs))


#how many taxa were found to be significant at the genus-level, but not individual ASVs in there? 
#how many ASVs were agglomerated to genus level, and genus-level was SDA, but not SDA at an ASV level?
#66 genera encompassing 4636 ASVs:
usvi_sda_resolution_summary.df %>%
  dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
  dplyr::summarise(num_asvs_in_sig_genera = sum(num_agglom_asvs),
                   num_genera = length(taxonomy))
# dplyr::summarise(num_sig_genera_not_asvs = sum(num_sig_genera))

#how many genera were SDA and had SDA ASVs:
#97 genera
usvi_sda_resolution_summary.df %>%
  dplyr::filter(num_sig_asvs > 0 & num_sig_genera > 0) %>%
  dplyr::summarise(num_genera = length(taxonomy))

#how many ASVs are significant, and belong to SDA genera?
#250 ASVs in 97 genera:
usvi_sda_resolution_summary.df %>%
  dplyr::filter(num_sig_genera > 0 & num_sig_asvs > 0) %>%
  dplyr::summarise(num_asvs_in_sig_genera = sum(num_sig_asvs),
                   num_not_sig_asvs = sum(num_agglom_asvs) - sum(num_sig_asvs),
                   num_genera = length(taxonomy))



#how many asvs were found to be significant but not at the agglomerated genera? 
#8 ASVs reprsenting 7 genera:
usvi_sda_resolution_summary.df %>%
  dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0)

#how many ASVS in those 7 genera were not SDA?
#138 ASVs:
usvi_sda_resolution_summary.df %>%
  dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
  dplyr::ungroup(.) %>%
  dplyr::summarise(num_agglom_asvs = sum(num_agglom_asvs) - sum(num_sig_asvs),
                   num_genera = length(taxonomy))

#how many ASVs were in sig genera but not sig at ASV levels?
#7803 ASVs
usvi_sda_resolution_summary.df %>%
  dplyr::ungroup(.) %>%
  dplyr::summarise(num_not_sig_asvs = sum(num_agglom_asvs) - sum(num_sig_asvs),
                   num_genera = length(taxonomy))



# Plot the relative abundances of the 8061 ASVs in 170 genera -------------


usvi_sda_genera_relabund.df <- usvi_sda_asvs_in_sda_genera.df %>%
  dplyr::filter(!is.na(sda_agg_genus)) %>%
  droplevels %>%
  dplyr::distinct(asv_id, taxonomy) %>%
  dplyr::mutate(significant = 1) %>%
  dplyr::bind_rows(., usvi_sda_asvs_in_sda_genera.df %>%
                     dplyr::filter(is.na(sda_agg_genus) & !is.na(sda_asv)) %>%
                     dplyr::distinct(taxonomy, .keep_all = FALSE) %>%
                     droplevels) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  dplyr::left_join(., usvi_sw_genus.taxa.df %>%
                     dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
                     dplyr::select(asv_id, taxonomy) %>%
                     droplevels, by = join_by(taxonomy)) %>%
  dplyr::mutate(genera_id = across(starts_with("asv_id")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(genera_id, taxonomy, significant) %>%
  dplyr::right_join(usvi_sw_genus.df %>%
                     dplyr::filter(sample_id %in% usvi_selected_metadata[["sample_id"]]) %>%
                     dplyr::group_by(sample_id) %>%
                     dplyr::mutate(relabund = relabund(counts)) %>%
                     droplevels, ., by = join_by("asv_id" == "genera_id")) %>%
  dplyr::rename(genera_id = "asv_id") %>%
  dplyr::mutate(significant = dplyr::case_when(is.na(significant) ~ 0,
                                               .default = significant)) %>%
  droplevels

usvi_sda_asvs_and_relatives_relabund.df <- usvi_sda_asvs_in_sda_genera.df %>%
  dplyr::filter(!is.na(sda_asv)) %>%
  dplyr::select(asv_id, taxonomy) %>%
  dplyr::mutate(significant = 1) %>%
  dplyr::full_join(., usvi_prok_asvs.taxonomy %>%
                     dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
                     dplyr::select(asv_id, taxonomy) %>%
                     dplyr::semi_join(usvi_sda_asvs_in_sda_genera.df,
                                      by = join_by(taxonomy))) %>%
  split(., f = .$taxonomy) %>%
  purrr::map(., ~.x %>%
               dplyr::arrange(asv_id) %>%
               tibble::rowid_to_column(var = "rank_order") %>%
               droplevels) %>%
  bind_rows(.) %>%
  dplyr::right_join(usvi_prok_asvs.df %>%
                     dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
                      dplyr::group_by(sample_ID) %>%
                      dplyr::mutate(relabund = relabund(counts)) %>%
                      dplyr::rename(sample_id = "sample_ID") %>%
                     droplevels, ., by = join_by(asv_id)) %>%
  dplyr::mutate(significant = dplyr::case_when(is.na(significant) ~ 0,
                                               .default = significant)) %>%
  droplevels

#first, plot the relative abundances of the SDA genera and genera of ASVs found to be SDA at the individual level
print(
  ggplot(data = usvi_sda_genera_relabund.df)
  + geom_boxplot(aes(x = genera_id, y = relabund, fill = taxonomy), shape = 21, show.legend = FALSE)
)


#what is the distribution of agglomerated ASVs in each of the 170 genera?
temp_df <- usvi_sda_asvs_and_relatives_relabund.df %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(desc(rank_order)) %>%
  dplyr::select(taxonomy, rank_order) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy)))



#two version of histograms:
{
  print(ggplot(data = temp_df)
        + geom_histogram(aes(x = log10(rank_order)),  
                         # bins= 5, binwidth = 0.8,
                         breaks = log10(c(1, 10, 50, 100, 500, 1500)),
                         color = "black", fill = "grey")
        + scale_x_continuous(name = "Abundance of ASVs in agglomerated genera",
                             breaks = (log10(c(1, 10, 50, 100, 500)*(1.5))),
                             labels = c("Between 1 and 10", "Between 50 and 11", "Between 100 and 51",
                                        "Between 500 and 101", "More than 501")
        )
        + scale_y_continuous(name = "Number of genus-level taxa")
        + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "grey"),
                strip.text.x = element_text(angle = 0),
                strip.text.y = element_text(angle = 0))
  )
  
  print(ggplot(data = temp_df)
        + geom_histogram(aes(x = rank_order), breaks = c(0, 10, 50, 100, 500, 1500), color = "black", fill = "grey")
        + scale_x_continuous(transform = scales::as.transform(t_pseudolog10), breaks =c(3, 20, 75, 300, 750),
                             labels = c("Between 1 and 10", "Between 50 and 11", "Between 100 and 51",
                                        "Between 500 and 101", "More than 501"),
                             name = "Abundance of ASVs in agglomerated genera")
        + scale_y_continuous(name = "Number of genus-level taxa")
        + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "grey"),
                strip.text.x = element_text(angle = 0),
                strip.text.y = element_text(angle = 0))
  )  
  }



#97 genera had individual ASVs that were also SDA. 
#what is the distribution of the relative abundances of those 250 ASVs and the other ASVs agglomerated into those genera?
sda_genera_with_sda_asvs_idx <- intersect((usvi_sda_genera_relabund.df %>%
            dplyr::ungroup(.) %>%
            dplyr::filter(significant > 0) %>%
            dplyr::distinct(taxonomy) %>%
            unlist %>% as.character), (usvi_sda_asvs_and_relatives_relabund.df %>%
                                        dplyr::ungroup(.) %>%
                                        dplyr::filter(significant > 0) %>%
                                        dplyr::distinct(taxonomy) %>%
                                        unlist %>% as.character))


temp_df2 <- temp_df %>%
  dplyr::filter(taxonomy %in% sda_genera_with_sda_asvs_idx) %>%
  dplyr::left_join(., usvi_sda_asvs_and_relatives_relabund.df %>%
                     dplyr::filter(taxonomy %in% sda_genera_with_sda_asvs_idx) %>%
                     dplyr::ungroup(.) %>%
                     dplyr::distinct(asv_id, taxonomy, significant) %>%
                     dplyr::group_by(taxonomy) %>%
                     dplyr::summarise(num_sig = sum(significant, na.rm = TRUE)),
                   by = join_by(taxonomy)) %>%
  dplyr::arrange(desc(num_sig), desc(rank_order)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  droplevels


usvi_sda_genera_relabund.list <- usvi_sda_asvs_and_relatives_relabund.df %>%
  dplyr::select(sample_id, asv_id, rank_order, relabund, taxonomy, significant) %>%
  dplyr::rename(taxon = "asv_id") %>%
  dplyr::mutate(res = "asv") %>%
  bind_rows(., usvi_sda_genera_relabund.df %>%
              dplyr::mutate(taxon = paste0("g_", genera_id)) %>%
              dplyr::select(sample_id, taxon, relabund, taxonomy, significant) %>%
              dplyr::mutate(res = "agg_genus") %>%
              droplevels) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  dplyr::mutate(rank_order = dplyr::case_when(is.na(rank_order) ~ 0,
                                              .default = rank_order)) %>%
  droplevels

temp_list <- temp_df2 %>%
                         dplyr::slice_max(n = 10, order_by = num_sig) %>%
                         dplyr::select(taxonomy) %>%
  droplevels %>%
  dplyr::left_join(., usvi_sda_genera_relabund.list,
                   by = join_by(taxonomy), relationship = "many-to-many", multiple = "all") %>%
  droplevels

#oneata time:
{
  
  # 
  # 
  # namevar <- unique(usvi_sda_genera_relabund.list[["taxonomy"]])
  # print(
  #   ggplot(data = usvi_sda_genera_relabund.list)
  #   + theme_bw()
  #   + geom_boxplot(aes(x = rank_order, y = relabund, fill = factor(significant), group = rank_order, color =factor(significant)), 
  #                  alpha = 0.8,
  #                  outliers = TRUE, outlier.shape = NA)
  #   + scale_y_continuous(name = "Relative abundance (%)",
  #                        # transform = scales::as.transform(t_pseudolog10), 
  #                        breaks = c(0, 0.1, 0.5, 1, 5, 10, 25),
  #                        transform = "log1p", 
  #                        expand = expansion(mult = c(0.01,0.1)))
  #   + scale_x_continuous(name = "Rank order of taxon")
  #   # + scale_fill_manual(values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white")
  #   + scale_discrete_manual(aesthetics = c("color"), values = c("grey", "black"), breaks = c(0, 1),
  #                           drop = TRUE)
  #   + scale_discrete_manual(aesthetics = c("fill"), values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white",
  #                           drop = TRUE)
  #   + guides(color = "none")
  #   + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
  #           panel.grid.minor = element_blank(),
  #           panel.grid.major.x = element_blank(),
  #           panel.grid.major.y = element_line(color = "grey"),
  #           strip.text.x = element_text(angle = 0),
  #           strip.text.y = element_text(angle = 0))
  #   + ggtitle(paste0("Taxa of ", namevar))
  # )
  
  
}

y <- unique(temp_list[["taxonomy"]]) %>% unlist %>% as.character(.)

# for(i in seq_len(2)){
for(i in seq_len(length(y))){
  namevar <- y[i]
  temp_df_to_plot <- temp_list %>%
    dplyr::filter(taxonomy == namevar) %>%
    droplevels
  x_breaks <- temp_df_to_plot %>%
    dplyr::distinct(rank_order) %>%
    dplyr::arrange(rank_order) %>%
    tibble::deframe(.)
  x_labels <- c("Agglom", x_breaks[-1])
  
  temp_g <- print(
    ggplot(data = temp_df_to_plot)
    + theme_bw()
    + geom_boxplot(aes(x = rank_order, y = relabund, fill = factor(significant), group = rank_order, color =factor(significant)), 
                   alpha = 0.8,
                   outliers = TRUE, outlier.shape = NA)
    + scale_y_continuous(name = "Relative abundance (%)",
                         transform = "log1p", 
                         expand = expansion(mult = c(0.01,0.1)))
    + scale_x_continuous(name = "Rank order of taxon", labels = x_labels, breaks = x_breaks)
    + scale_discrete_manual(aesthetics = c("color"), values = c("grey", "black"), breaks = c(0, 1),
                            drop = TRUE)
    + scale_discrete_manual(aesthetics = c("fill"), values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white",
                            drop = TRUE)
    + guides(color = "none")
    + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey"),
            strip.text.x = element_text(angle = 0),
            strip.text.y = element_text(angle = 0))
    + ggtitle(paste0(namevar))
  )
  assign(paste0("g_sda_asv_relabund_", i), temp_g, envir = .GlobalEnv)
}

gpatch <- lapply(apropos("^g_sda_asv_relabund_.*$", mode = "list"),
                 get) %>%
  purrr::reduce(., `+`) + 
  patchwork::plot_layout(guides = "collect")
  # patchwork::plot_layout(guides = "collect") & theme(legend.position="none")

g_top_sda_asv_relabund <- gpatch + patchwork::plot_annotation(title = "Distribution of SDA ASVs and their SDA genera",
                                                             tag_levels = "A")


if(!any(grepl("top_sda_asv_relabund", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_top_sda_asv_relabund-", Sys.Date(), ".png"),
         g_top_sda_asv_relabund, 
         width = 16, height = 16, units = "in")
}

#what about the 66 genera whose individual ASVs were not SDA?
sda_genera_not_sda_asvs_idx <- setdiff((usvi_sda_genera_relabund.df %>%
                                          dplyr::ungroup(.) %>%
                                          dplyr::filter(significant > 0) %>%
                                          dplyr::distinct(taxonomy) %>%
                                          unlist %>% as.character), (usvi_sda_asvs_and_relatives_relabund.df %>%
                                                                       dplyr::ungroup(.) %>%
                                                                       dplyr::filter(significant > 0) %>%
                                                                       dplyr::distinct(taxonomy) %>%
                                                                       unlist %>% as.character))

#these are the 4636 ASVs belonging to the 66 SDA genera
temp_df3 <- usvi_prok_asvs.taxonomy %>%
  dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
  dplyr::filter(taxonomy %in% sda_genera_not_sda_asvs_idx) %>%
  dplyr::select(asv_id, taxonomy) %>%
  split(., f = .$taxonomy) %>%
  purrr::map(., ~.x %>%
               dplyr::arrange(asv_id) %>%
               tibble::rowid_to_column(var = "rank_order") %>%
               droplevels) %>%
  bind_rows(.) %>%
  dplyr::left_join(., (usvi_prok_asvs.df %>%
                         dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
                         dplyr::group_by(sample_ID) %>%
                         dplyr::mutate(relabund = relabund(counts)) %>%
                         dplyr::rename(sample_id = "sample_ID") %>%
                         droplevels),
                   by = join_by(asv_id)) %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  droplevels
temp_df3 %>%
  dplyr::mutate(relabund = dplyr::case_when(relabund > 0 ~ relabund,
                                            .default = NA)) %>%
  # dplyr::group_by(taxonomy) %>%
  dplyr::summarise(relabund = sum(relabund, na.rm = TRUE), .by = taxonomy) %>% dplyr::arrange(desc(relabund)) %>%
  # dplyr::ungroup(.) %>% dplyr::summarise(quantiles = quantile(relabund, probs = seq(0, 1, 0.25), na.rm = TRUE)) %>%
  droplevels

#examinethe top 10 genera:
temp_list2 <- usvi_sda_genera_relabund.list %>%
  dplyr::filter(taxonomy %in% sda_genera_not_sda_asvs_idx) %>%
  dplyr::filter(res == "agg_genus") %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::distinct(taxon, taxonomy, .keep_all = TRUE) %>%
  dplyr::slice_max(n = 10, order_by = relabund) %>%
  dplyr::semi_join((usvi_sda_genera_relabund.list %>%
                     dplyr::filter(taxonomy %in% sda_genera_not_sda_asvs_idx) %>%
                     # dplyr::filter(res == "agg_genus") %>%
                     droplevels), ., by = join_by(taxonomy, taxon)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  droplevels

#here are the ASVs belonging to those top 10 genera
temp_list3 <- temp_list2 %>%
  bind_rows(temp_df3 %>%
              dplyr::filter(taxonomy %in% unique(temp_list2[["taxonomy"]])) %>%
              droplevels %>%
            dplyr::rename(taxon = "asv_id") %>%
              dplyr::mutate(significant = 0)) %>%
  dplyr::select(sample_id, taxon, rank_order, relabund, taxonomy, significant) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  dplyr::mutate(rank_order = dplyr::case_when(is.na(rank_order) ~ 0,
                                              .default = rank_order)) %>%
  droplevels

y <- unique(temp_list3[["taxonomy"]]) %>% unlist %>% as.character(.)

# for(i in seq_len(2)){
for(i in seq_len(length(y))){
  namevar <- y[i]
  temp_df_to_plot <- temp_list3 %>%
    dplyr::filter(taxonomy == namevar) %>%
    droplevels
  x_breaks <- temp_df_to_plot %>%
    ungroup(.) %>%
    dplyr::distinct(rank_order) %>%
    dplyr::arrange(rank_order) %>%
    tibble::deframe(.)
  x_labels <- c("Agglom", x_breaks[-1])
  
  temp_g <- print(
    ggplot(data = temp_df_to_plot)
    + theme_bw()
    + geom_boxplot(aes(x = rank_order, y = relabund, fill = factor(significant), group = rank_order, color =factor(significant)), 
                   alpha = 0.8,
                   outliers = TRUE, outlier.shape = NA)
    + scale_y_continuous(name = "Relative abundance (%)",
                         transform = "log1p", 
                         expand = expansion(mult = c(0.01,0.1)))
    + scale_x_continuous(name = "Rank order of taxon", labels = x_labels, breaks = x_breaks)
    + scale_discrete_manual(aesthetics = c("color"), values = c("grey", "black"), breaks = c(0, 1),
                            drop = TRUE)
    + scale_discrete_manual(aesthetics = c("fill"), values = c("white", "navy"), breaks = c(0, 1), labels = c("not", "significant"), name = "Significance", na.value = "white",
                            drop = TRUE)
    + guides(color = "none")
    + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey"),
            strip.text.x = element_text(angle = 0),
            strip.text.y = element_text(angle = 0))
    + ggtitle(paste0(namevar))
  )
  assign(paste0("g_sda_genus_relabund_", i), temp_g, envir = .GlobalEnv)
}

gpatch2 <- lapply(apropos("^g_sda_genus_relabund_.*$", mode = "list"),
                 get) %>%
  purrr::reduce(., `+`) + 
  patchwork::plot_layout(guides = "collect")
# patchwork::plot_layout(guides = "collect") & theme(legend.position="none")

g_top_sda_genera_not_asvs_relabund <- gpatch2 + patchwork::plot_annotation(title = "Distribution of SDA genera and their non-SDA ASVs",
                                                              tag_levels = "A")


if(!any(grepl("top_sda_genera_not_asvs_relabund", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_top_sda_genera_not_asvs_relabund-", Sys.Date(), ".png"),
         g_top_sda_genera_not_asvs_relabund, 
         width = 16, height = 16, units = "in")
}

{
# #make a bubble plot?
# {
#   temp_df <- usvi_sda_resolution_summary.df %>%
#     dplyr::arrange(desc(num_sig_asvs)) %>%
#     # dplyr::arrange(taxonomy) %>%
#     tibble::rowid_to_column(var = "order") %>%
#     tidyr::pivot_longer(., cols = !c("taxonomy", "order"),
#                         names_to = "metric",
#                         values_to = "count")  %>%
#     dplyr::mutate(metric = factor(metric, levels = c("num_agglom_asvs", "num_not_sig_asvs", "num_sig_asvs", "num_sig_genera"))) %>%
#     dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#     dplyr::filter(!grepl("num_sig_genera", metric)) %>%
#     droplevels
#   
#   print(
#     ggplot(data = temp_df)
#     + geom_point(aes(x = metric, y = taxonomy, fill = metric, size = count), shape = 21)
#     # + scale_size_continuous(transform = "log10", reverse = TRUE)
#     + scale_size_area(transform = "pseudo_log", max_size = 10, breaks = c(1, 10, 100, 1000, 1500))
#     + coord_flip()
#   )
#   
# }
# 
# # tempsankey_df <- mtcars %>%
# #   make_long(cyl, vs, am, gear, carb)
# 
# #make a sankey graph?
# 
# # usvi_sda_sankey_maybe.df <- usvi_prok_asvs.taxonomy %>%
# #   dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
# #   # dplyr::mutate(genera_id = paste0("g_", asv_id)) %>%
# #   dplyr::select(asv_id, taxonomy) %>%
# #   dplyr::semi_join(., usvi_sda_asvs_in_sda_genera.df, by = join_by(taxonomy)) %>%
# #   droplevels %>%
# #   dplyr::left_join(., usvi_sw_genus.taxa.df %>%
# #                      dplyr::mutate(genera_id = paste0("g_", asv_id)) %>%
# #                      dplyr::select(genera_id, asv_id),
# #                    by = join_by(asv_id)) %>%
# #   dplyr::arrange(taxonomy, asv_id) %>%
# #   dplyr::group_by(taxonomy) %>%
# #   tidyr::fill(genera_id, .direction = "down") %>%
# #   # dplyr::left_join(., usvi_sda_asvs_in_sda_genera.df %>%
# #   #                    dplyr::filter(!is.na(sda_agg_genus) & is.na(sda_asv)) %>%
# #   #                    # dplyr::filter(!is.na(sda_asv) & !is.na(sda_agg_genus)) %>%
# #   #                    # dplyr::mutate(genera_id = paste0("g_", asv_id)) %>% 
# #   #                    dplyr::select(asv_id) %>% dplyr::mutate(significant = 1), 
# #   #                  by = join_by(asv_id), relationship = "one-to-many", multiple = "first") %>%
# #   droplevels
# # 
# # usvi_sda_sankey_maybe.df <- usvi_sda_asvs_in_sda_genera.df %>%
# #   dplyr::filter(!is.na(sda_agg_genus) & is.na(sda_asv)) %>%
# #   dplyr::mutate(resolution = "genus") %>%
# #   dplyr::mutate(label = paste0("g_", asv_id)) %>%
# #   dplyr::select(resolution, label) %>% dplyr::mutate(significant = 1) %>%
# #   bind_rows(., usvi_sda_asvs_in_sda_genera.df %>%
# #               dplyr::filter(is.na(sda_agg_genus) & !is.na(sda_asv)) %>%
# #               dplyr::mutate(resolution = "genus") %>%
# #               dplyr::mutate(label = paste0("g_", asv_id)) %>%
# #               dplyr::select(resolution, label) %>% dplyr::mutate(significant = 0)) %>%
# #   bind_rows(., usvi_sda_asvs_in_sda_genera.df %>%
# #               dplyr::filter(!is.na(sda_agg_genus) & !is.na(sda_asv)) %>%
# #               dplyr::mutate(resolution = "asv") %>%
# #               dplyr::mutate(label = asv_id) %>%
# #               dplyr::select(resolution, label) %>% dplyr::mutate(significant = 1))
# # 
# # 
# # usvi_sda_asvs_in_sda_genera.df %>%
# #   dplyr::filter(asv_id %in% usvi_sw_genus.taxa.df[["asv_id"]]) %>%
# #   dplyr::mutate(genera_id = paste0("g_", asv_id)) %>%
# #   dplyr::select(genera_id, taxonomy) %>%
# #   dplyr::full_join(., usvi_sda_asvs_in_sda_genera.df %>%
# #                      dplyr::filter(!is.na(sda_agg_genus) & is.na(sda_asv)) %>%
# #                      # dplyr::filter(!is.na(sda_asv) & !is.na(sda_agg_genus)) %>%
# #                      dplyr::mutate(genera_id = paste0("g_", asv_id)) %>% dplyr::select(genera_id, taxonomy) %>% dplyr::mutate(significant = 1), 
# #                    by = join_by(taxonomy, genera_id)) %>%
# #                      # dplyr::select(asv_id, taxonomy),
# #                    # by = join_by(taxonomy)) %>%
# #   dplyr::bind_rows(., usvi_sda_asvs_in_sda_genera.df %>%
# #                      dplyr::filter(!is.na(sda_asv)) %>%
# #                      dplyr::select(asv_id, taxonomy) %>%
# #                      dplyr::mutate(significant = 1))
# 
# 
# usvi_sda_sankey_maybe.df <- usvi_sda_resolution_summary.df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::summarise(count = sum(num_agglom_asvs)) %>%
#   dplyr::mutate(node = "all_asvs") %>%
#   dplyr::mutate(degree = 1) %>%
#   dplyr::select(degree, node, count) %>%
#   dplyr::mutate(sub_of = "all_asvs") %>%
#   dplyr::full_join(., usvi_sda_resolution_summary.df %>%
#                      dplyr::filter(num_sig_asvs > 0 & num_sig_genera > 0) %>%
#                      dplyr::summarise(count = length(taxonomy)) %>%
#                      dplyr::mutate(degree = 2) %>%
#                      dplyr::mutate(node = "sig_genera_w_sig_asvs") %>%
#                      dplyr::mutate(sub_of = "all_asvs") %>%
#                      dplyr::mutate(significant = 1)) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_genera > 0 & num_sig_asvs > 0) %>%
#                          dplyr::summarise(count = sum(num_sig_asvs)) %>%
#                          dplyr::mutate(degree = 3) %>%
#                          dplyr::mutate(node = "sig_asv") %>%
#                          dplyr::mutate(sub_of = "sig_genera_w_sig_asvs") %>%
#                          dplyr::mutate(significant = 1))) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_genera > 0 & num_sig_asvs > 0) %>%
#                          dplyr::summarise(count = sum(num_agglom_asvs) - sum(num_sig_asvs)) %>%
#                          dplyr::mutate(degree = 3) %>%
#                          dplyr::mutate(node = "not_sig_asv") %>%
#                          dplyr::mutate(sub_of = "sig_genera_w_sig_asvs") %>%
#                          dplyr::mutate(significant = 0)
#   )) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
#                          dplyr::summarise(count = length(taxonomy)) %>%
#                          dplyr::mutate(node = "sig_genera_not_sig_asvs") %>%
#                          dplyr::mutate(degree = 2) %>%
#                          dplyr::mutate(sub_of = "all_asvs") %>%
#                          dplyr::mutate(significant = 1))) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
#                          dplyr::summarise(count = sum(num_agglom_asvs)) %>%
#                          dplyr::mutate(node = "not_sig_asv") %>%
#                          dplyr::mutate(degree = 3) %>%
#                          dplyr::mutate(sub_of = "sig_genera_not_sig_asvs") %>%
#                          dplyr::mutate(significant = 0))) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_asvs == 0 & num_sig_genera > 0) %>%
#                          dplyr::summarise(count = sum(num_agglom_asvs) - sum(num_agglom_asvs)) %>%
#                          dplyr::mutate(node = "sig_asv") %>%
#                          dplyr::mutate(degree = 3) %>%
#                          dplyr::mutate(sub_of = "sig_genera_not_sig_asvs") %>%
#                          dplyr::mutate(significant = 1))) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
#                          dplyr::summarise(count = length(taxonomy)) %>%
#                          dplyr::mutate(node = "not_sig_genera_sig_asvs") %>%
#                          dplyr::mutate(degree = 2) %>%
#                          dplyr::mutate(sub_of = "all_asvs") %>%
#                          dplyr::mutate(significant = 0))) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
#                          dplyr::summarise(count = sum(num_sig_asvs)) %>%
#                          dplyr::mutate(node = "sig_asv") %>%
#                          dplyr::mutate(degree = 3) %>%
#                          dplyr::mutate(sub_of = "not_sig_genera_sig_asvs") %>%
#                          dplyr::mutate(significant = 1))) %>%
#   dplyr::full_join(., (usvi_sda_resolution_summary.df %>%
#                          dplyr::filter(num_sig_asvs > 0 & num_sig_genera == 0) %>%
#                          dplyr::summarise(count = sum(num_agglom_asvs) - sum(num_sig_asvs)) %>%
#                          dplyr::mutate(node = "not_sig_asv") %>%
#                          dplyr::mutate(degree = 3) %>%
#                          dplyr::mutate(sub_of = "not_sig_genera_sig_asvs") %>%
#                          dplyr::mutate(significant = 0))) %>%
#   dplyr::rename(next_node = "node",
#                 node = "sub_of") %>%
#   dplyr::relocate(node, .after = degree) %>%
#   dplyr::mutate(degree = dplyr::case_when(degree == 1 ~ "all",
#                                           degree == 2 ~ "genus",
#                                           degree == 3 ~ "asv",
#                                           .default = NA)) %>%
#   droplevels
# 
# temp_df <- usvi_sda_sankey_maybe.df %>%
#   dplyr::filter(grepl("all", node)) %>%
#   dplyr::select(degree, next_node, count, significant) %>%
#   dplyr::rename(grouping = "next_node",
#                 resolution = "degree") %>%
#   bind_rows(., usvi_sda_sankey_maybe.df %>%
#               dplyr::filter(!grepl("all", node)) %>%
#               dplyr::select(degree, node, count, significant) %>%
#               dplyr::rename(grouping = "node",
#                             resolution = "degree") %>%
#               droplevels) %>%
#   dplyr::mutate(resolution = factor(resolution, levels = c("all", "genus", "asv")))
# 
# 
# print(
#   ggplot(data = temp_df)
#   + geom_col(aes(x = resolution, y = count, fill = factor(significant)))
#   + scale_fill_discrete()
# )
# 
# 
# # temp_df2 <- temp_df %>%
# # dplyr::filter(grepl("all", resolution)) %>%
# #   dplyr::rename(x = "resolution", node= "grouping", value = "count") %>%
# #   dplyr::mutate(next_x = "genus", next_node = NA) %>%
# #   bind_rows(., (temp_df %>%
# #                   dplyr::filter(grepl("genus", resolution)) %>%
# #                   dplyr::rename(x = "resolution", node= "grouping", value = "count") %>%
# #                   dplyr::mutate(next_x = "asv", next_node = NA))) %>%
# #   bind_rows(., (temp_df %>%
# #                   dplyr::filter(grepl("asv", resolution)) %>%
# #                   dplyr::rename(x = "resolution", node= "grouping", value = "count") %>%
# #                   dplyr::mutate(next_x =NA, next_node = NA)))
# # 
# # temp_df3 <- temp_df %>%
# #   ggsankey::make_long(resolution, grouping, value = count)
# 
# #spending too much time on this.
# {
#   
#   
#   temp_df4 <- temp_df %>%
#     dplyr::mutate(degree = dplyr::case_when(degree == 1 ~ "all",
#                                             degree == 2 ~ "genus",
#                                             degree == 3 ~ "asv",
#                                             .default = NA)) %>%
#     dplyr::filter(grepl("all", node)) %>%
#     dplyr::rename(genus = "count") %>%
#     dplyr::rename(grouping = "next_node") %>%
#     dplyr::select(degree, grouping, genus, significant) %>%
#     droplevels
#   temp_df3 <- temp_df %>%
#     dplyr::mutate(degree = dplyr::case_when(degree == 1 ~ "all",
#                                             degree == 2 ~ "genus",
#                                             degree == 3 ~ "asv",
#                                             .default = NA)) %>%
#     dplyr::filter(!grepl("all", node)) %>%
#     dplyr::rename(asv = "count") %>%
#     dplyr::rename(grouping = "node") %>%
#     dplyr::select(degree, grouping, asv, significant) %>%
#     droplevels
#   
#   temp_df2 <- tidyr::pivot_wider(temp_df3, id_cols = c("degree", "significant"),
#                                  names_from = "grouping",
#                                  values_from = "asv") %>%
#     bind_rows(., (temp_df4 %>%
#                     tidyr::pivot_wider(., id_cols = c("degree", "significant"),
#                                        names_from = "grouping",
#                                        values_from = "genus"))) %>%
#     tidyr::drop_na(significant) %>%
#     dplyr::select(-all_asvs) %>%
#     # ggsankey::make_long(sig_genera_w_sig_asvs, sig_genera_not_sig_asvs, not_sig_genera_sig_asvs) %>%
#     droplevels
#   
#   # temp_df2 <- temp_df %>%
#   #   dplyr::mutate(degree = dplyr::case_when(degree == 1 ~ "all",
#   #                                           degree == 2 ~ "genus",
#   #                                           degree == 3 ~ "asv",
#   #                                           .default = NA)) %>%
#   #   tidyr::pivot_wider(., id_cols = c("node", "next_node"),
#   #                      names_from = c("degree"),
#   #                      values_from = "count") %>%
#   #   dplyr::mutate(asv = dplyr::coalesce(asv, all)) %>%
#   #   dplyr::select(-all) %>%
#   #   # tidyr::pivot_longer(., cols = c("genus", "asv"),
#   #   #                     names_to = NULL,
#   #   #                     values_to = "count") %>%
#   #   #  # tidyr::drop_na(count) %>%
#   #   droplevels
#   
#   temp_df2 <- usvi_sda_resolution_summary.df %>%
#     dplyr::select(taxonomy, contains("asv")) %>%
#     ggsankey::make_long(num_agglom_asvs, num_sig_asvs, num_not_sig_asvs) %>%
#     droplevels
#   {
#     temp_df <- usvi_prok_asvs.taxonomy %>%
#       dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
#       dplyr::select(asv_id, taxonomy) %>%
#       droplevels %>%
#       dplyr::semi_join(., usvi_sda_asvs_in_sda_genera.df %>% 
#                          dplyr::filter(!is.na(sda_asv)) %>%
#                          dplyr::distinct(taxonomy), 
#                        by = join_by(taxonomy)) %>%
#       dplyr::anti_join(., usvi_sda_asvs_in_sda_genera.df %>% 
#                          dplyr::filter(!is.na(sda_asv)) %>%
#                          dplyr::distinct(asv_id, taxonomy)) %>%
#       dplyr::mutate(taxonomy_level = "asv",
#                     significance = 0) %>%
#       # dplyr::mutate(not_sig_asv = 0) %>%
#       bind_rows((usvi_sda_asvs_in_sda_genera.df %>%
#                    dplyr::filter(!is.na(sda_asv)) %>%
#                    dplyr::distinct(asv_id, taxonomy) %>%
#                    dplyr::mutate(taxonomy_level = "asv",
#                                  significance = 1) %>%
#                    # dplyr::mutate(is_sig_asv = 2) %>%
#                    droplevels), .) %>%
#       bind_rows(., (usvi_sda_asvs_in_sda_genera.df %>%
#                       dplyr::filter(!is.na(sda_agg_genus)) %>%
#                       dplyr::distinct(asv_id, taxonomy) %>%
#                       dplyr::mutate(taxonomy_level = "genus",
#                                     significance = 1) %>%
#                       # dplyr::mutate(asv_id = paste0("g_", asv_id)) %>%
#                       droplevels)) %>%
#       bind_rows(., (usvi_sda_asvs_in_sda_genera.df %>%
#                       dplyr::filter(!is.na(sda_asv) & is.na(sda_agg_genus)) %>%
#                       dplyr::distinct(asv_id, taxonomy) %>%
#                       dplyr::mutate(taxonomy_level = "asv",
#                                     significance = 1) %>%
#                       droplevels)) %>%
#       dplyr::distinct(asv_id, taxonomy, taxonomy_level, significance, .keep_all = TRUE) %>%
#       # dplyr::mutate(base_asv = 1) %>%
#       tidyr::pivot_wider(., id_cols = NULL,
#                          names_from = "taxonomy_level",names_prefix = "sig_",
#                          # values_fill = 0,
#                          values_from = "significance") %>%
#       droplevels
#     
#     
#     #node, x, next_node, next_x
#     #x options: sig_asv, asv,  not_sig_asv
#     #node options: 0, 1, 2
#     library(ggsankey)
#     
#     temp_df2 <- temp_df %>%
#       dplyr::filter(grepl("Synechococcus", taxonomy)) %>%
#       ggsankey::make_long(asv_id, genus, asv) %>%
#       # ggsankey::make_long(genus, asv, asv_id) %>%
#       # dplyr::filter(!(x == "asv_id")) %>%
#       dplyr::left_join(., usvi_prok_asvs.taxonomy %>%
#                          dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
#                          dplyr::select(asv_id, taxonomy) %>%
#                          dplyr::rename(label = "taxonomy") %>%
#                          droplevels,
#                        by = join_by("node" == "asv_id"), multiple = "first", relationship = "many-to-many") %>%
#       droplevels
#     
#   }
#   
# }
# 
# {
#   (
#     ggplot(data = temp_df2, 
#            aes(x = x, next_x = next_x,
#                node = node, next_node = next_node,        
#                # label = label,
#                color = factor(node),
#                fill = factor(node)))
#     + ggsankey::geom_sankey()
#     # + ggsankey::geom_alluvial(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)
#     # # + ggsankey::geom_alluvial_label(size = 3, color = "white", fill = "gray40")
#     # + ggsankey::theme_alluvial(base_size = 18)
#     # + ggsankey::geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)
#     # # + ggsankey::geom_sankey_label(size = 3, color = "white", fill = "gray40")
#     # + ggsankey::theme_sankey(base_size = 18)
#     # + coord_flip()
#     + theme_dark()
#     + theme(
#       axis.title = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks = element_blank(),
#       panel.grid = element_blank())
#     # + scale_fill_manual(values = annotation_taxa_colors,
#     #                     labels = names(annotation_taxa_colors),
#     #                     breaks = names(annotation_taxa_colors),
#     #                     drop = TRUE)
#     # + scale_color_manual(values = annotation_taxa_colors,
#     #                      labels = names(annotation_taxa_colors),
#     #                      breaks = names(annotation_taxa_colors),
#     #                      drop = TRUE)
#   )
#   
# }
# 
# usvi_sda_resolution_summary.df %>%
#   # tidyr::pivot_longer(cols = !c("taxonomy"),
#   #                     names_to = "metric",
#   #                     values_to = "num_asvs") %>%
#   # dplyr::mutate(num_asvs = tidyr::replace_na(num_asvs, 0)) %>%
#   dplyr::summarise(num_sig_genera = sum(num_sig_genera, na.rm = TRUE))
# 
# 
# # usvi_sda_resolution_venn.list <- usvi_sda_resolution_summary.df %>%
# #   tidyr::pivot_longer(., cols = !c("taxonomy"),
# #                       names_to = "grouping",
# #                       values_to = "count") %>%
# #   dplyr::mutate(grouping = factor(grouping, levels = c("num_agglom_asvs", "num_sig_genera", "num_sig_asvs", "num_sig_genera_not_asvs", "num_sig_asvs_not_genera"))) %>%
# #   # tidyr::complete(taxonomy, grouping, fill = list(count = 0)) %>%
# #   dplyr::filter(!grepl("not", grouping)) %>%
# #   droplevels %>%
# #   # split(., f = .$grouping) %>%
# #   split(., f = .$taxonomy) %>%
# #   map(., ~.x %>%
# #         dplyr::select(grouping, count) %>%
# #         droplevels) %>%
# #   map(., ~.x %>%
# #         tibble::deframe(.)) %>%
# # map(., ~.x %>%
# #         as.list)
# 
# #this makes a list for each taxonomy, always containing two character vectors: sda_asv and sda_agg_genus
# #the membership of those vectors may have NA, meaning that the taxon was observed in the other vector but not that one
# 
# # usvi_sda_resolution_venn.list <- usvi_sda_asvs_in_sda_genera.df %>%
# #   tidyr::pivot_longer(., cols = !c("asv_id", "taxonomy"),
# #                       names_to = "taxon_resolution",
# #                       values_to = "test_type") %>%
# #   dplyr::distinct(asv_id, taxonomy, taxon_resolution, .keep_all = TRUE) %>%
# #   split(., f = .$taxonomy) %>%
# #   map(., ~.x %>%
# #         dplyr::mutate(asv_id = dplyr::case_when(!is.na(test_type) ~ asv_id,
# #                                                 .default = NA)) %>%
# #         split(., f = .$taxon_resolution) %>%
# #         map(., ~.x %>%
# #               dplyr::select(asv_id) %>%
# #               tibble::deframe(.) %>%
# #               Filter(Negate(function(x) is.null(unlist(x))), .)) #remove the Null lists
# #   )
# 
# #this makes a list for each taxonomy, each of which contains up to two character vectors: sda_asv and sda_agg_genus
# #the membership of those vectors is only the sig.diff.abund ASV or genus
# 
# usvi_sda_resolution_venn.list <- usvi_sda_asvs_in_sda_genera.df %>%
#   dplyr::mutate(across(starts_with("sda_"), ~dplyr::case_when(!is.na(.x) ~ asv_id,
#                                                               .default = NA))) %>%
#   tidyr::pivot_longer(., cols = !c("taxonomy"),
#                       names_to = "taxon_resolution",
#                       values_to = "asv_id") %>%
#   dplyr::filter(!(taxon_resolution == "asv_id")) %>%
#   split(., f = .$taxonomy) %>%
#   map(., ~.x %>%
#         dplyr::select(starts_with("taxon_resolution"), asv_id) %>%
#         tidyr::drop_na(asv_id) %>%
#         split(., f = .$taxon_resolution) %>%
#         map(., ~.x %>% dplyr::select(asv_id))) %>%
#   map(., ~.x %>%
#         purrr::list_flatten(., name_spec = "{outer}")) %>%
#   map_depth(., 2, ~.x %>% simplify) %>%
#   map(., ~.x %>%
#         Filter(Negate(function(x) is.null(unlist(x))), .))
# 
# # temp_list <- usvi_prok_asvs.taxonomy %>%
# #                      dplyr::mutate(taxonomy = gsub(";", "; ", taxonomy)) %>%
# #                      dplyr::select(taxonomy, asv_id) %>%
# #   dplyr::semi_join(., usvi_sda_resolution_summary.df,
# #                    by = join_by(taxonomy))
# 
# # p1 <- ggVennDiagram(temp_list[["Alphaproteobacteria; AEGEAN-169 marine group"]])
# # p1 <- ggVennDiagram(temp_list[["Alphaproteobacteria; SAR11 Clade Ia"]])
# p1 <- ggVennDiagram(usvi_sda_resolution_venn.list[["Alphaproteobacteria; SAR11 Clade Ia"]])
# p1
}


