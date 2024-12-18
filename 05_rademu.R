# 05_rademu.R

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

# library(tidyverse)
# # library(phyloseq)
# library(ggvegan)
# library(ggtext)
# library(viridis)
# library(patchwork)
# library(viridisLite)
# library(pals)
# library(radEmu)

library(magrittr)
library(dplyr)
library(ggplot2)
library(stringr)
library(radEmu)

# Custom functions --------------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_custom_functions.R"))){
  cat("Loading custom functions.")
  source(paste0(getwd(), "/", "00_custom_functions.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 

F_parallel_rademu <- function(category, mdata, rademu_design, df, modmat, covariate){
  rademu_design <- rlang::parse_expr(rademu_design)
  temp_res <- emuFit(formula = eval(rademu_design),
                     data = mdata,
                     Y = df, 
                     fitted_model = modmat,
                     refit = FALSE,
                     # test_kj = data.frame(k = covariate_to_test, 
                     test_kj = data.frame(k = covariate, 
                                          j = category))
  # this was for just one covariate to test, but if we have multiple provided in the function call:
  # temp_res <- tibble::tribble(
  #   ~covariate, ~category, ~score_stat, ~pval,
  #   temp_res$coef$covariate[category], temp_res$coef$category[category], temp_res$coef$score_stat[category], temp_res$coef$pval[category]
  # )
  # return(temp_res)
  temp_res2 <- temp_res$coef %>%
    dplyr::filter(category_num == {{category}}) %>%
    dplyr::select(covariate, category, score_stat, pval)
  return(temp_res2)
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
  droplevels
temp_df <- usvi_prok_asvs.taxa %>%
  dplyr::filter(keep == FALSE) %>%
  droplevels
  
#replace NA in taxonomy with last known level

usvi_prok_filled.taxa.df <- usvi_prok_asvs.taxa %>%
  dplyr::filter(keep == TRUE) %>%
  dplyr::select(-keep) %>%
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
  cli::cli_alert_warning("Please process the metabolomics sample data previously.")
}


#replace NA in taxonomy with last known level

usvi_prok_filled.taxa.df <- usvi_prok_asvs.taxa %>%
  dplyr::left_join(., usvi_prok_decontam_idx, by = join_by(asv_id)) %>%
  dplyr::filter(keep == TRUE) %>%
  dplyr::select(-keep) %>%
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





# Prepare color palette ---------------------------------------------------

# site_lookup <- data.frame(site = c("LB_seagrass", "Tektite", "Yawzi", "control_extraction", "control_pcr", "control_seq"),
#                           label = c("Lameshur Bay seagrass", "Tektite Reef", "Yawzi Reef",
#                                     "Control (DNA Extraction)", "Control (PCR)", "Control (Sequencing)")) %>%
#   tibble::deframe(.)
# site_colors <- pals::kelly(22)[6:(5+length(site_lookup))] %>%
#   # site_colors <- viridisLite::cividis(n = length(site_lookup), direction = 1) %>%
#   setNames(., names(site_lookup))
# sampling_time_lookup <- data.frame(sampling_time = c("dawn", "peak_photo"),
#                                    label = c("Dawn", "Peak photosynthesis")) %>%
#   tibble::deframe(.)
# sampling_time_colors <- pals::ocean.haline(n = length(sampling_time_lookup)) %>%
#   setNames(., names(sampling_time_lookup))
# sampling_day_lookup <- data.frame(sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5"),
#                                   label = c("20210122", "20210123", "20210124", "20210125", "20210126")) %>%
#   tibble::deframe(.)
# sampling_day_colors <- pals::ocean.thermal(n = length(sampling_day_lookup)) %>%
#   setNames(., names(sampling_day_lookup))


# Read in metabolites -----------------------------------------------------


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




# Conduct rademu on genera in the sites -----------------------------------


# drop <- c("CINAR_BC_73", "CINAR_BC_43", "CINAR_BC_105", "CINAR_BC_69")
drop <- c("CINAR_BC_73")

#remember with radEmu we have to unload phyloseq :)
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


usvi_sw_genus.mat <- usvi_sw_genus.tbl %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.) %>%
  t() 


#using just the four taxa:

#genera with relative abundance related to sample site:
#Cryomorphaceae
#Flavobacteriaceae
#Prochlorococcus
#SAR86

{
  # keep <- c("SAR86 clade", "Prochlorococcus MIT9313", "Flavobacteriaceae", "Cryomorphaceae")
  # keep_idx <- usvi_sw_genus.taxa.df %>%
  #   dplyr::filter(grepl(paste0(keep, collapse = "|"), Genus)) %>%
  #   droplevels %>%
  #   dplyr::select(asv_id) %>%
  #   unlist %>%
  #   as.character()
  # 
  # usvi_selected_genus.mat <- usvi_sw_genus.mat %>%
  #   as.data.frame() %>%
  #   dplyr::select(all_of(keep_idx)) %>%
  #   droplevels
  # 
  # usvi_selected_fit <- radEmu::emuFit(formula = ~ site,
  #                                     data = usvi_selected_metadata,
  #                                     Y = usvi_selected_genus.mat,
  #                                     run_score_tests = FALSE)
  # 
  # usvi_selected_rademu.df <- usvi_selected_fit$coef %>%
  #   tibble::as_tibble(.) %>%
  #   dplyr::rename(asv_id = "category") %>%
  #   dplyr::left_join(., usvi_prok_filled.taxa.df %>%
  #                      dplyr::select(asv_id, Genus) %>%
  #                      droplevels) %>%
  #   dplyr::filter(abs(estimate) > 1) %>%
  #   dplyr::arrange(Genus, desc(abs(estimate))) %>%
  #   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  #   droplevels
  # 
  # ggplot(usvi_selected_rademu.df) + 
  #   geom_point(aes(x = asv_id, y = estimate, color = Genus), size = .5) +
  #   geom_errorbar(aes(x = asv_id, ymin = lower, ymax = upper, color = Genus), width = .25) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   facet_wrap(.~covariate) +
  #   labs(x = "Category",
  #        y = "Estimate") + 
  #   coord_cartesian(ylim = c(-5,10))
  # 
  # covariate_to_test <- which("siteTektite" == usvi_selected_fit$B %>% rownames)
  # temp_res <-radEmu::emuFit(formula = ~ site, 
  #                           data = usvi_selected_metadata,
  #                           Y = usvi_selected_genus.mat,
  #                           # return_both_score_pvals = TRUE,
  #                           fitted_model = usvi_selected_fit,
  #                           refit = FALSE,
  #                           test_kj = data.frame(k = covariate_to_test,
  #                                                j = c(1:5)),
  #                           run_score_tests = FALSE) 
  # 
  
  # F_rademu <- function(category){
  #   temp_res <- emuFit(formula = ~ site,
  #                      data = usvi_selected_metadata,
  #                      Y = usvi_selected_genus.mat,
  #                      fitted_model = usvi_selected_fit,
  #                      refit = FALSE,
  #                      test_kj = data.frame(k = covariate_to_test,
  #                                           j = category))
  #   temp_res <- data.frame(covariate = temp_res$coef$covariate[category],
  #                          category = temp_res$coef$category[category],
  #                          score_stat = temp_res$coef$score_stat[category],
  #                          pval = temp_res$coef$pval[category])
  #   return(temp_res)
  # }
  
  # score_res <- parallel::mclapply(asvs_to_test,
  #                       F_rademu,
  #                       mc.cores = nthreads)
  # if (!is.null(score_res)) {
  #   full_score <- sapply(1:length(score_res), 
  #                        function(x) score_res[[x]]$coef$score_stat[x])
  #   full_pval <- sapply(1:length(score_res), 
  #                       function(x) score_res[[x]]$coef$pval[x])
  #   full_coef <- usvi_selected_fit$coef %>%
  #     dplyr::select(-score_stat, -pval) %>%
  #     dplyr::filter(category_num %in% asvs_to_test) %>%
  #     dplyr::filter(covariate == "siteTektite") %>%
  #     dplyr::mutate(score_stat = full_score,
  #                   pval = full_pval) %>%
  #   dplyr::arrange(full_pval)
  #   full_coef
  # }
}

#do it all:

# Rademu on just site -----------------------------------------------------


usvi_sw_fit <- radEmu::emuFit(formula = ~ site, 
                                    data = usvi_selected_metadata,
                                    Y = usvi_sw_genus.mat,
                                    run_score_tests = FALSE) 

usvi_sw_rademu.df <- usvi_sw_fit$coef %>%
  tibble::as_tibble(.) %>%
  dplyr::rename(asv_id = "category") %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(asv_id, Genus) %>%
                     droplevels) %>%
  dplyr::filter(abs(estimate) > 1) %>%
  dplyr::arrange(desc(abs(estimate)), Genus) %>%
  dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  droplevels

ggplot(usvi_sw_rademu.df) + 
  geom_point(aes(x = asv_id, y = estimate, color = Genus), size = .5, show.legend = FALSE) +
  geom_errorbar(aes(x = asv_id, ymin = lower, ymax = upper, color = Genus), width = .25, show.legend = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(.~covariate) +
  labs(x = "Category",
       y = "Estimate") + 
  coord_cartesian(ylim = c(-5,10))




#conduct robust test scores, first on the genera that have estimates > 0 in both Yawzi and Tektite compared to LB seagrass
# covariate_to_test <- which("siteTektite" == usvi_sw_fit$B %>% rownames)
covariate_to_test <- which(grepl("site", usvi_sw_fit$B %>% rownames))
# this small test works:
# # temp_res <- F_parallel_rademu(category = 2, mdata = usvi_selected_metadata, rademu_design = " ~ site", covariate = covariate_to_test,
# #                               df = usvi_selected_genus.mat, modmat = usvi_selected_fit)

# #this test works as well:
# temp_res <- F_parallel_rademu(category = 2, mdata = usvi_selected_metadata, rademu_design = " ~ site", covariate = covariate_to_test,
#                   df = usvi_sw_genus.mat, modmat = usvi_sw_fit)
# # 
# # temp_res$coef %>%
# #   dplyr::filter(category_num == 2) %>%
# #   dplyr::select(covariate, category, score_stat, pval)
# # 
# # temp_res2 <- tibble::tribble(
# #   ~covariate, ~category, ~score_stat, ~pval,
# #   temp_res$coef$covariate[category], temp_res$coef$category[category], temp_res$coef$score_stat[category], temp_res$coef$pval[category]
# # )

# asvs_to_test <- seq_len(length(keep_idx))

taxa_to_test.df <- data.frame(asv_id = colnames(usvi_sw_genus.mat)) %>%
  dplyr::left_join(., bind_rows((usvi_sw_rademu.df %>%
                     dplyr::group_by(asv_id) %>%
                     dplyr::filter(estimate > 0) %>%
                     dplyr::summarise(obs = length(covariate)) %>%
                     dplyr::ungroup(.) %>%
                     dplyr::filter(obs == 2) %>%
                     dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                     dplyr::mutate(group = "high_reefs")),
                     (usvi_sw_rademu.df %>%
                        dplyr::group_by(asv_id) %>%
                        dplyr::filter(estimate < 0) %>%
                        dplyr::summarise(obs = length(covariate)) %>%
                        dplyr::ungroup(.) %>%
                        dplyr::filter(obs == 2) %>%
                        dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                        dplyr::mutate(group = "low_reefs"))),
                   by = join_by(asv_id)) %>%
  tibble::rowid_to_column(var = "rowid") %>%
  droplevels

# score_res_list2 <- furrr::future_map(asvs_to_test[1:5],
#                                     F_parallel_rademu,
#                                     mdata = usvi_selected_metadata, rademu_design = " ~ site",
#                                     df = usvi_sw_genus.mat, modmat = usvi_sw_fit,
#                                     .options = furrr::furrr_options(seed = 48105, globals = TRUE),
#                                     .progress = TRUE)
# 
# score_res_list <- parallel::mclapply(asvs_to_test[1:5],
#                                      F_parallel_rademu,
#                                      mdata = usvi_selected_metadata, rademu_design = " ~ site",
#                                     df = usvi_sw_genus.mat, modmat = usvi_sw_fit,
#                                     mc.cores = nthreads)

#go big or go home.

asvs_to_test <- taxa_to_test.df %>%
  dplyr::filter(grepl("high_reefs", group)) %>%
  dplyr::select(rowid) %>%
  tibble::deframe(.)

score_highreef_res_list <- parallel::mclapply(asvs_to_test,
                                     F_parallel_rademu,
                                     mdata = usvi_selected_metadata, rademu_design = " ~ site", covariate = covariate_to_test,
                                     df = usvi_sw_genus.mat, modmat = usvi_sw_fit,
                                     mc.cores = nthreads)
score_highreef_res.df <- score_highreef_res_list %>%
  bind_rows(., .id = "asv_id")

asvs_to_test <- taxa_to_test.df %>%
  dplyr::filter(grepl("low", group)) %>%
  dplyr::select(rowid) %>%
  tibble::deframe(.)

score_lowreef_res_list <- parallel::mclapply(asvs_to_test,
                                              F_parallel_rademu,
                                              mdata = usvi_selected_metadata, rademu_design = " ~ site", covariate = covariate_to_test,
                                              df = usvi_sw_genus.mat, modmat = usvi_sw_fit,
                                              mc.cores = nthreads)
score_lowreef_res.df <- score_lowreef_res_list %>%
  bind_rows(., .id = "asv_id")

if(!any(grepl("genera_reef_highlow", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
  save(score_highreef_res_list, score_lowreef_res_list, file = paste0(projectpath, "/", "usvi_rademu_genera_reef_highlow-", Sys.Date(), ".RData"))
}


# Rademu on site + time ---------------------------------------------------


#what if we do site and time of day?


#try first with the four taxa to make sure radEmu can take a multicomponent formula with interactions:

{
  keep <- c("SAR86 clade", "Prochlorococcus MIT9313", "Flavobacteriaceae", "Cryomorphaceae")
  keep_idx <- usvi_sw_genus.taxa.df %>%
    dplyr::filter(grepl(paste0(keep, collapse = "|"), Genus)) %>%
    droplevels %>%
    dplyr::select(asv_id) %>%
    unlist %>%
    as.character()
  
  # metadata <- usvi_selected_metadata
  usvi_selected_genus.mat <- usvi_sw_genus.mat %>%
    as.data.frame() %>%
    dplyr::select(all_of(keep_idx)) %>%
    droplevels
  
  # temp_rad_selected_no_ave.tbl <- radEmu::emuFit(data = usvi_selected_metadata,
  #                                     # formula = ~ 0 + site:sampling_time,
  #                                     # formula = ~ 0 + site + site:sampling_time + sampling_time,
  #                                     formula = ~ site + sampling_time:site, #this groups per-site, all samples as a covariate irrespective of sampling time (i.e. "Yawzi dawn + peakPhoto" and "Tektite dawn + peakphoto"), compared to LB_Seagrass_dawn as baseline/intercept
  #                                     Y = usvi_selected_genus.mat,
  #                                     # return_nullB = TRUE,
  #                                     run_score_tests = FALSE) %>%
  #   purrr::pluck("B")
  temp_rad_selected_no_ave2.tbl <- radEmu::emuFit(data = usvi_selected_metadata,
                                                 formula = ~ 0 + sampling_time:site, #this is the formula we need to use
                                                 Y = usvi_selected_genus.mat,
                                                 # return_nullB = TRUE,
                                                 run_score_tests = FALSE) %>%
    purrr::pluck("B")
  
  #just to make sure, calculate the B matrix using just Yawzi samples.
  #if the estimates for these four taxa for covariate "peak_photo" are the same as "Yawzi:peak_photo", and "Intercept" or "dawn" is the same as "Yawzi:dawn", then we're okay
  {
    # small_metadata <- usvi_selected_metadata %>%
    #   dplyr::filter(grepl("Yawzi", site)) %>%
    #   droplevels
    # small_y <- usvi_selected_genus.mat[rownames(small_metadata),]
    # 
    # temp_rad_selected_no_ave3.tbl <- radEmu::emuFit(data = small_metadata,
    #                                                 formula = ~ sampling_time, #this is the formula we need to use
    #                                                 Y = small_y,
    #                                                 # return_nullB = TRUE,
    #                                                 run_score_tests = FALSE) %>%
    #   purrr::pluck("B")
    # 
    }
  #it's okay :)
  
  # 
  # 

  usvi_selected_fit <- radEmu::emuFit(data = usvi_selected_metadata,
                                      formula = ~ 0 + site:sampling_time,
                                      Y = usvi_selected_genus.mat,
                                      # return_nullB = TRUE,
                                      run_score_tests = FALSE)

  # #the means for each ASV across each group:
  # usvi_selected_rademu.tbl <- usvi_selected_fit[["B"]] %>%
  #   tibble::as_tibble(., rownames = "covariate") %>%
  #   tidyr::pivot_longer(., cols = !c(covariate),
  #                       names_to = "asv_id",
  #                       values_to = "estimate") %>%
  #   droplevels
  # #radEmu uses the first combination of each variables' first factors, as the Intercept
  # 
  # # temp_df <- dplyr::left_join(usvi_selected_rademu.tbl,
  # #                      (usvi_selected_rademu.tbl %>%
  # #   dplyr::filter(grepl("LB", covariate) & grepl("dawn", covariate)) %>%
  # #     dplyr::select(asv_id, estimate) %>%
  # #     dplyr::rename(control = "estimate") %>%
  # #     # dplyr::mutate(covariate = "control") %>%
  # #   droplevels)) %>%
  # #   dplyr::mutate()
  # #   
  # # 
  # usvi_selected_rademu.df <- usvi_selected_fit$coef %>%
  #   tibble::as_tibble(.) %>%
  #   dplyr::rename(asv_id = "category") %>%
  #   dplyr::left_join(., usvi_prok_filled.taxa.df %>%
  #                      dplyr::select(asv_id, Genus) %>%
  #                      droplevels) %>%
  #   # dplyr::filter(abs(estimate) > 1) %>%
  #   dplyr::arrange(Genus, desc(abs(estimate))) %>%
  #   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  #   droplevels
  # 
  # ggplot(usvi_selected_rademu.df) +
  #   geom_point(aes(x = asv_id, y = estimate, color = Genus), size = .5) +
  #   geom_errorbar(aes(x = asv_id, ymin = lower, ymax = upper, color = Genus), width = .25) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   facet_wrap(.~covariate) +
  #   labs(x = "Category",
  #        y = "Estimate") +
  #   coord_cartesian(ylim = c(-5,10))
  # 
  # covariate_to_test <- which("siteTektite" == usvi_selected_fit$B %>% rownames)
  # temp_res <-radEmu::emuFit(formula = ~ site, 
  #                           data = usvi_selected_metadata,
  #                           Y = usvi_selected_genus.mat,
  #                           # return_both_score_pvals = TRUE,
  #                           fitted_model = usvi_selected_fit,
  #                           refit = FALSE,
  #                           test_kj = data.frame(k = covariate_to_test,
  #                                                j = c(1:5)),
  #                           run_score_tests = FALSE) 
  # 
  
  # F_rademu <- function(category){
  #   temp_res <- emuFit(formula = ~ site,
  #                      data = usvi_selected_metadata,
  #                      Y = usvi_selected_genus.mat,
  #                      fitted_model = usvi_selected_fit,
  #                      refit = FALSE,
  #                      test_kj = data.frame(k = covariate_to_test,
  #                                           j = category))
  #   temp_res <- data.frame(covariate = temp_res$coef$covariate[category],
  #                          category = temp_res$coef$category[category],
  #                          score_stat = temp_res$coef$score_stat[category],
  #                          pval = temp_res$coef$pval[category])
  #   return(temp_res)
  # }
  
  # score_res <- parallel::mclapply(asvs_to_test,
  #                       F_rademu,
  #                       mc.cores = nthreads)
  # if (!is.null(score_res)) {
  #   full_score <- sapply(1:length(score_res), 
  #                        function(x) score_res[[x]]$coef$score_stat[x])
  #   full_pval <- sapply(1:length(score_res), 
  #                       function(x) score_res[[x]]$coef$pval[x])
  #   full_coef <- usvi_selected_fit$coef %>%
  #     dplyr::select(-score_stat, -pval) %>%
  #     dplyr::filter(category_num %in% asvs_to_test) %>%
  #     dplyr::filter(covariate == "siteTektite") %>%
  #     dplyr::mutate(score_stat = full_score,
  #                   pval = full_pval) %>%
  #   dplyr::arrange(full_pval)
  #   full_coef
  # }
}

#would it help to make a "dummy" sample(s) that is the average of genera across all seawater samples, for a baselien of comparison?
#no:
{
  # usvi_sw_genus_ave.mat <- usvi_sw_genus.tbl %>%
  #   bind_cols(., 
  #             "average1" = c(usvi_sw_genus.tbl %>%
  #                              tibble::column_to_rownames(var = "asv_id") %>%
  #                              apply(., 1, mean) ),
  #             "average2" = c(usvi_sw_genus.tbl %>%
  #                              tibble::column_to_rownames(var = "asv_id") %>%
  #                              apply(., 1, mean) )) %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame(.) %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   as.data.frame(.) %>%
  #   t() 
  # 
  # usvi_selected_ave_metadata <- metabolomics_sample_metadata %>%
  #   dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  #   dplyr::filter(grepl("seawater", sample_type)) %>%
  #   dplyr::select(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
  #   dplyr::mutate(across(c(metab_deriv_label, sample_id, sample_type, sampling_date, sampling_time, sampling_day, site), ~factor(.x))) %>%
  #   # dplyr::select(sample_id, metab_deriv_label) %>%
  #   dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  #   # dplyr::mutate(rownames = metab_deriv_label) %>% tibble::column_to_rownames(var = "metab_deriv_label") %>%
  #   bind_rows(., tibble::tribble(
  #     ~sample_id, ~sample_type, ~site, ~metab_deriv_label, ~sampling_time, ~sampling_date, ~sampling_day,
  #     "average1", "control", "control", "average1", "dawn", "control", "control",
  #     # "average2", "control", "control", "average2", "peak_photo", "control", "control"
  #     "average2", "control", "control", "average2", "dawn", "control", "control"
  #   )) %>%
  #   dplyr::mutate(across(c(site, sample_type, 
  #                          # sampling_time, 
  #                          sampling_date, sampling_day), ~forcats::fct_relevel(.x, "control"))) %>%
  #   dplyr::select(sample_id, metab_deriv_label, sample_type, site, sampling_time, sampling_day) %>%
  #   dplyr::mutate(rownames = sample_id) %>% tibble::column_to_rownames(var = "rownames") %>%
  #   droplevels
  # 
  # usvi_selected_genus_ave.mat <- usvi_sw_genus_ave.mat %>%
  #   as.data.frame() %>%
  #   dplyr::select(all_of(keep_idx)) %>%
  #   droplevels
  # 
  # 
  # temp_rad_selected_ave.tbl <- radEmu::emuFit(data = usvi_selected_ave_metadata,
  #                                  formula = ~ 0 + site:sampling_time, #this is the formula we need to use
  #                                  Y = usvi_selected_genus_ave.mat,
  #                                  # return_nullB = TRUE,
  #                                  run_score_tests = FALSE) %>%
  #   purrr::pluck("B")
  # 
  
}



usvi_sw_site_time_fit <- radEmu::emuFit(data = usvi_selected_metadata,
                                        # formula = ~ site + sampling_time + site:sampling_time, 
                                        # formula = ~ 0 + site:sampling_time,
                                        # formula = ~ sample_type,
                                        formula = ~ 0 + site:sampling_time,
                              Y = usvi_sw_genus.mat,
                              run_score_tests = FALSE) 

# #what does the "baseline" (LB at dawn) look like for asvs?
temp_df <- usvi_sw_site_time_fit[["B"]] %>%
  tibble::as_tibble(., rownames = "covariate") %>%
  dplyr::filter(grepl("LB", covariate)) %>%
  dplyr::filter(grepl("dawn", covariate)) %>%
  tidyr::pivot_longer(., cols = !c(covariate),
                      names_to = "asv_id",
                      values_to = "estimate") %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(asv_id, Genus) %>%
                     droplevels) %>%
  # dplyr::filter(abs(estimate) > 1) %>%
  dplyr::arrange(desc(abs(estimate)), Genus) %>%
  dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  droplevels

# temp_df2 <- usvi_sw_site_time_fit[["B"]] %>%
#   tibble::as_tibble(., rownames = "covariate") %>%
#   dplyr::filter(!grepl("LB", covariate) | !grepl("dawn",  covariate)) %>%
#   # dplyr::filter(grepl("dawn", covariate)) %>%
#   tidyr::pivot_longer(., cols = !c(covariate),
#                       names_to = "asv_id",
#                       values_to = "estimate") %>%
#   dplyr::left_join(., usvi_prok_filled.taxa.df %>%
#                      dplyr::select(asv_id, Genus) %>%
#                      droplevels) %>%
#   # dplyr::filter(abs(estimate) > 1) %>%
#   dplyr::arrange(desc(abs(estimate)), Genus) %>%
#   dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
#   droplevels %>%
#   dplyr::left_join(., temp_df %>%
#                      dplyr::select(asv_id, estimate) %>%
#                      dplyr::rename(control_estimate = "estimate") %>%
#                      droplevels,
#                    by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
#   dplyr::mutate(adj_estimate = estimate - control_estimate)

# # #radEmu uses the first combination of each variables' first factors, as the Intercept
# # 
# ggplot(temp_df %>%
#          # dplyr::filter(abs(estimate) > 1) %>%
#          droplevels) +
#   geom_point(aes(x = asv_id, y = estimate, color = Genus), size = .5, show.legend = FALSE) +
#   # geom_errorbar(aes(x = asv_id, ymin = lower, ymax = upper, color = Genus), width = .25, show.legend = FALSE) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_wrap(.~covariate) +
#   labs(x = "Category",
#        y = "Estimate") +
#   coord_cartesian(ylim = c(-5,10))


usvi_sw_site_time_rademu.df <- usvi_sw_site_time_fit[["B"]] %>%
  tibble::as_tibble(., rownames = "covariate") %>%
  dplyr::filter(grepl("LB", covariate)) %>%
  dplyr::filter(grepl("dawn", covariate)) %>%
  tidyr::pivot_longer(., cols = !c(covariate),
                      names_to = "asv_id",
                      values_to = "estimate") %>%
  dplyr::right_join(., usvi_sw_site_time_fit$coef %>%
                      tibble::as_tibble(.) %>%
                      dplyr::rename(asv_id = "category") %>%
                      droplevels, by = join_by(covariate, asv_id), relationship = "many-to-many", multiple = "all") %>%
  dplyr::mutate(estimate = across(contains("estimate")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(-ends_with(c(".x", ".y"))) %>%
    dplyr::arrange(asv_id, desc(covariate)) %>%
  dplyr::ungroup(.) %>%
  dplyr::group_by(asv_id) %>%
  tidyr::fill(., category_num, .direction = "down") %>%
  dplyr::relocate(estimate, .after = category_num) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(asv_id, Genus) %>%
                     droplevels, by = join_by(asv_id)) %>%
  dplyr::arrange(desc(abs(estimate)), Genus) %>%
  dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
  droplevels %>%
  dplyr::left_join(., usvi_sw_site_time_fit[["B"]] %>%
                     tibble::as_tibble(., rownames = "covariate") %>%
                     dplyr::filter(!grepl("LB", covariate) | !grepl("dawn",  covariate)) %>%
                     tidyr::pivot_longer(., cols = !c(covariate),
                                         names_to = "asv_id",
                                         values_to = "estimate") %>%
                     dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                                        dplyr::select(asv_id, Genus) %>%
                                        droplevels) %>%
                     dplyr::arrange(desc(abs(estimate)), Genus) %>%
                     dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
                     droplevels %>%
                     dplyr::left_join(., temp_df %>%
                                        dplyr::select(asv_id, estimate) %>%
                                        dplyr::rename(control_estimate = "estimate") %>%
                                        droplevels,
                                      by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
                     dplyr::mutate(adj_estimate = estimate - control_estimate) %>%
                     dplyr::select(covariate, asv_id, adj_estimate),
                    by = join_by(covariate, asv_id)) %>%
  dplyr::relocate(adj_estimate, .after = asv_id) %>%
  droplevels

# ggplot(usvi_sw_site_time_rademu.df %>%
#          # dplyr::filter(abs(estimate) > 1) %>%
#          dplyr::filter(abs(adj_estimate) > 1) %>%
#          droplevels) +
#   geom_point(aes(x = asv_id, y = estimate, color = Genus), size = .5, show.legend = FALSE) +
#   geom_errorbar(aes(x = asv_id, ymin = lower, ymax = upper, color = Genus), width = .25, show.legend = FALSE) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_wrap(.~covariate) +
#   labs(x = "Category",
#        y = "Estimate") +
#   coord_cartesian(ylim = c(-5,10))




#conduct robust test scores, first on the genera that have estimates > 0 in both Yawzi and Tektite compared to LB seagrass
# covariate_to_test <- which("siteTektite" == usvi_sw_fit$B %>% rownames)
# covariate_to_test <- which(grepl("site", usvi_sw_site_time_fit$B %>% rownames))

# this small test works:
# # temp_res <- F_parallel_rademu(category = 2, mdata = usvi_selected_metadata, rademu_design = " ~ site", covariate = covariate_to_test,
# #                               df = usvi_selected_genus.mat, modmat = usvi_selected_fit)

# #this test works as well:
# temp_res <- F_parallel_rademu(category = 2, mdata = usvi_selected_metadata, rademu_design = " ~ site", covariate = covariate_to_test,
#                   df = usvi_sw_genus.mat, modmat = usvi_sw_fit)
# # 
# # temp_res$coef %>%
# #   dplyr::filter(category_num == 2) %>%
# #   dplyr::select(covariate, category, score_stat, pval)
# # 
# # temp_res2 <- tibble::tribble(
# #   ~covariate, ~category, ~score_stat, ~pval,
# #   temp_res$coef$covariate[category], temp_res$coef$category[category], temp_res$coef$score_stat[category], temp_res$coef$pval[category]
# # )

# asvs_to_test <- seq_len(length(keep_idx))

# temp_df <- usvi_sw_site_time_rademu.df %>%
#   # dplyr::filter(grepl("dawn", covariate)) %>%
#   # dplyr::filter(!grepl("LB", covariate)) %>%
#   dplyr::filter(abs(adj_estimate) > 1) %>%
#   dplyr::filter(adj_estimate > 0) %>%
#   # dplyr::filter(abs(estimate) > 1) %>%
#   # dplyr::filter(estimate > 0) %>%
#   dplyr::group_by(asv_id) %>%
#   # dplyr::summarise(obs = length(covariate)) %>%
#   dplyr::ungroup(.) %>%
#   # dplyr::filter(obs >= 2) %>%
#   # dplyr::distinct(asv_id, .keep_all = TRUE) %>%
#   droplevels

# taxa_to_test.df <- data.frame(asv_id = colnames(usvi_sw_genus.mat)) %>%
#   dplyr::left_join(., bind_rows((usvi_sw_site_time_rademu.df %>%
#                                    dplyr::filter(grepl("dawn", covariate)) %>%
#                                    # dplyr::filter(abs(estimate) > 1) %>%
#                                    # dplyr::filter(estimate > 0) %>%
#                                    dplyr::filter(abs(adj_estimate) > 1) %>%
#                                    dplyr::filter(adj_estimate > 0) %>%
#                                    dplyr::group_by(asv_id) %>%
#                                    dplyr::summarise(obs = length(covariate)) %>%
#                                    dplyr::ungroup(.) %>%
#                                    dplyr::filter(obs >= 2) %>%
#                                    dplyr::distinct(asv_id, .keep_all = FALSE) %>%
#                                    dplyr::mutate(group1 = "high_dawn")))) %>%
#   dplyr::left_join(., bind_rows((usvi_sw_site_time_rademu.df %>%
#                                    dplyr::filter(grepl("photo", covariate)) %>%
#                                    # dplyr::filter(abs(estimate) > 1) %>%
#                                    # dplyr::filter(estimate > 0) %>%
#                                    dplyr::filter(abs(adj_estimate) > 1) %>%
#                                    dplyr::filter(adj_estimate < 0) %>%
#                                    # dplyr::filter(adj_estimate > 0) %>%
#                                    dplyr::group_by(asv_id) %>%
#                                    dplyr::summarise(obs = length(covariate)) %>%
#                                    dplyr::ungroup(.) %>%
#                                    dplyr::filter(obs >= 2) %>%
#                                    dplyr::distinct(asv_id, .keep_all = FALSE) %>%
#                                    dplyr::mutate(group2 = "low_dawn"))),
#                                    # dplyr::mutate(group2 = "high_photo"))),
#                    by = join_by(asv_id)) %>%
#   tibble::rowid_to_column(var = "rowid") %>%
#   # tidyr::drop_na(group1, group2) %>%
#   droplevels

taxa_to_test2.df <- data.frame(asv_id = colnames(usvi_sw_genus.mat)) %>%
  dplyr::left_join(., bind_rows((usvi_sw_site_time_rademu.df %>%
                                   dplyr::filter(grepl("LB", covariate)) %>%
                                   dplyr::filter(abs(adj_estimate) > 1) %>%
                                   dplyr::filter(adj_estimate > 0) %>%
                                   droplevels %>%
                                   dplyr::group_by(asv_id) %>%
                                   dplyr::summarise(obs = length(covariate)) %>%
                                   dplyr::ungroup(.) %>%
                                   dplyr::filter(obs >= 1) %>%
                                   dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                                   dplyr::mutate(group = "high_seagrass")),
                                (usvi_sw_site_time_rademu.df %>%
                                   dplyr::filter(!grepl("LB", covariate)) %>%
                                   dplyr::filter(abs(adj_estimate) > 1) %>%
                                   dplyr::filter(adj_estimate < 0) %>%
                                   droplevels %>%
                                   dplyr::group_by(asv_id) %>%
                                   dplyr::summarise(obs = length(covariate)) %>%
                                   dplyr::ungroup(.) %>%
                                   dplyr::filter(obs >= 4) %>%
                                   dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                                   dplyr::mutate(group = "high_seagrass")),
                                (usvi_sw_site_time_rademu.df %>%
                                   dplyr::filter(grepl("LB", covariate)) %>%
                                   dplyr::filter(abs(adj_estimate) > 1) %>%
                                   dplyr::filter(adj_estimate < 0) %>%
                                   droplevels %>%
                                   dplyr::group_by(asv_id) %>%
                                   dplyr::summarise(obs = length(covariate)) %>%
                                   dplyr::ungroup(.) %>%
                                   dplyr::filter(obs >= 1) %>%
                                   dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                                   dplyr::mutate(group = "low_seagrass")),
                                (usvi_sw_site_time_rademu.df %>%
                                   dplyr::filter(!grepl("LB", covariate)) %>%
                                   dplyr::filter(abs(adj_estimate) > 1) %>%
                                   dplyr::filter(adj_estimate > 0) %>%
                                   droplevels %>%
                                   dplyr::group_by(asv_id) %>%
                                   dplyr::summarise(obs = length(covariate)) %>%
                                   dplyr::ungroup(.) %>%
                                   dplyr::filter(obs >= 4) %>%
                                   dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                                   dplyr::mutate(group = "low_seagrass"))),
                   by = join_by(asv_id)) %>%
  tibble::rowid_to_column(var = "rowid") %>%
  tidyr::drop_na(group) %>%
  droplevels
                                

taxa_to_test.df <- data.frame(asv_id = colnames(usvi_sw_genus.mat)) %>%
  dplyr::left_join(., bind_rows((usvi_sw_site_time_rademu.df %>%
                                   dplyr::filter(grepl("dawn", covariate)) %>%
                                   dplyr::filter(!grepl("LB", covariate)) %>%
                                   dplyr::filter(abs(adj_estimate) > 1) %>%
                                   dplyr::filter(adj_estimate > 0) %>%
                                   # dplyr::filter(abs(estimate) > 1) %>%
                                   # dplyr::filter(estimate > 0) %>%
                                   droplevels %>%
                                   dplyr::group_by(asv_id) %>%
                                   dplyr::summarise(obs = length(covariate)) %>%
                                   dplyr::ungroup(.) %>%
                                   dplyr::filter(obs >= 2) %>%
                                   dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                                   dplyr::mutate(group = "high_dawn")),
                   (usvi_sw_site_time_rademu.df %>%
                      dplyr::filter(grepl("dawn", covariate)) %>%
                      dplyr::filter(!grepl("LB", covariate)) %>%
                      dplyr::filter(abs(adj_estimate) > 1) %>%
                      dplyr::filter(adj_estimate < 0) %>%
                                   # dplyr::filter(abs(estimate) > 1) %>%
                                   # dplyr::filter(estimate > 0) %>%
                                   droplevels %>%
                      dplyr::group_by(asv_id) %>%
                                   dplyr::summarise(obs = length(covariate)) %>%
                                   dplyr::ungroup(.) %>%
                                   dplyr::filter(obs >= 2) %>%
                                   dplyr::distinct(asv_id, .keep_all = FALSE) %>%
                     dplyr::mutate(group = "low_dawn"))),
                   by = join_by(asv_id)) %>%
  tibble::rowid_to_column(var = "rowid") %>%
  tidyr::drop_na(group) %>%
  droplevels

# score_res_list2 <- furrr::future_map(asvs_to_test[1:5],
#                                     F_parallel_rademu,
#                                     mdata = usvi_selected_metadata, rademu_design = " ~ site",
#                                     df = usvi_sw_genus.mat, modmat = usvi_sw_fit,
#                                     .options = furrr::furrr_options(seed = 48105, globals = TRUE),
#                                     .progress = TRUE)
# 
# score_res_list <- parallel::mclapply(asvs_to_test[1:5],
#                                      F_parallel_rademu,
#                                      mdata = usvi_selected_metadata, rademu_design = " ~ site",
#                                     df = usvi_sw_genus.mat, modmat = usvi_sw_fit,
#                                     mc.cores = nthreads)

#go big or go home.
covariate_to_test <- which(grepl("dawn", usvi_sw_site_time_fit$B %>% rownames))
asvs_to_test <- taxa_to_test.df %>%
  dplyr::filter(grepl("high", group)) %>%
  dplyr::filter(grepl("dawn", group)) %>%
  dplyr::select(rowid) %>%
  tibble::deframe(.)

score_1_res_list <- parallel::mclapply(asvs_to_test,
                                              F_parallel_rademu,
                                              mdata = usvi_selected_metadata, 
                                              # rademu_design = " ~ site * sampling_time", 
                                       rademu_design = " ~ 0 + site:sampling_time", 
                                              covariate = covariate_to_test,
                                              df = usvi_sw_genus.mat, 
                                              modmat = usvi_sw_site_time_fit,
                                              mc.cores = nthreads)
score_1_res.df <- score_1_res_list %>%
  bind_rows(., .id = "asv_id") %>%
  tidyr::drop_na(score_stat)

covariate_to_test <- which(grepl("dawn", usvi_sw_site_time_fit$B %>% rownames))
asvs_to_test <- taxa_to_test.df %>%
  dplyr::filter(grepl("low", group)) %>%
  dplyr::filter(grepl("dawn", group)) %>%
  dplyr::select(rowid) %>%
  tibble::deframe(.)

score_2_res_list <- parallel::mclapply(asvs_to_test,
                                       F_parallel_rademu,
                                       mdata = usvi_selected_metadata, 
                                       rademu_design = " ~ 0 + site:sampling_time", 
                                       covariate = covariate_to_test,
                                       df = usvi_sw_genus.mat, 
                                       modmat = usvi_sw_site_time_fit,
                                             mc.cores = nthreads)
score_2_res.df <- score_2_res_list %>%
  bind_rows(., .id = "asv_id") %>%
  tidyr::drop_na(score_stat)


# #can't do site vs site here...
# covariate_to_test <- which(grepl("Yawzi", usvi_sw_site_time_fit$B %>% rownames))
# asvs_to_test <- taxa_to_test2.df %>%
#   dplyr::select(rowid) %>%
#   tibble::deframe(.)
# 
# score_3_res_list <- parallel::mclapply(asvs_to_test,
#                                        F_parallel_rademu,
#                                        mdata = usvi_selected_metadata,
#                                        rademu_design = " ~ 0 + site:sampling_time",
#                                        covariate = covariate_to_test,
#                                        df = usvi_sw_genus.mat,
#                                        modmat = usvi_sw_site_time_fit,
#                                        mc.cores = nthreads)



covariate_to_test <- which(!grepl("LB", usvi_sw_site_time_fit$B %>% rownames))
asvs_to_test <- taxa_to_test2.df %>%
  dplyr::filter(grepl("high", group)) %>%
  dplyr::filter(grepl("seagrass", group)) %>%
  dplyr::select(rowid) %>%
  tibble::deframe(.)

score_4_res_list <- parallel::mclapply(asvs_to_test,
                                       F_parallel_rademu,
                                       mdata = usvi_selected_metadata, 
                                       rademu_design = " ~ 0 + site:sampling_time", 
                                       covariate = covariate_to_test,
                                       df = usvi_sw_genus.mat, 
                                       modmat = usvi_sw_site_time_fit,
                                       mc.cores = nthreads)

# covariate_to_test <- which(grepl("Tektite", usvi_sw_site_time_fit$B %>% rownames))
asvs_to_test <- taxa_to_test2.df %>%
  dplyr::filter(grepl("low", group)) %>%
  dplyr::filter(grepl("seagrass", group)) %>%
  dplyr::select(rowid) %>%
  tibble::deframe(.)

score_5_res_list <- parallel::mclapply(asvs_to_test,
                                       F_parallel_rademu,
                                       mdata = usvi_selected_metadata, 
                                       rademu_design = " ~ 0 + site:sampling_time", 
                                       covariate = covariate_to_test,
                                       df = usvi_sw_genus.mat, 
                                       modmat = usvi_sw_site_time_fit,
                                       mc.cores = nthreads)

if(!any(grepl("genera_site_time", list.files(projectpath, pattern = "usvi_rademu_.*.RData")))){
  save(score_1_res_list, score_2_res_list, score_4_res_list, score_5_res_list,
       file = paste0(projectpath, "/", "usvi_rademu_genera_site_time-", Sys.Date(), ".RData"))
}

rademu_reef_res.df <- bind_rows((score_1_res_list %>%
                                   bind_rows(.) %>% #these were genera that were found to be higher in samples collected at dawn
                                   dplyr::mutate(group = "high_dawn") %>%
                                   tidyr::drop_na(score_stat)),
                                (score_2_res_list %>%
                                   bind_rows(.) %>% # these were genera that were found to be higher in samples collected at peak photo
                                   dplyr::mutate(group = "low_dawn") %>%
                                   tidyr::drop_na(score_stat)),
                                (score_4_res_list %>%
                                   bind_rows(.) %>% #these were genera that were found to be higher in samples collected in seagrass
                                   dplyr::mutate(group = "high_seagrass") %>%
                                   tidyr::drop_na(score_stat)),
                                (score_5_res_list %>%
                                   bind_rows(.) %>% #these were genera that were found to be higher in samples collected in reef sites
                                   dplyr::mutate(group = "low_seagrass") %>%
                                   tidyr::drop_na(score_stat))) %>%
  dplyr::rename(asv_id = "category") %>%
  droplevels
