# 03_asv_summary.R

# Load packages -----------------------------------------------------------
library(tidyverse)

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library(BiocManager)
library(BiocParallel)
library(cli)
library(future)
library(furrr)
library(progressr)
library(parallel)
library(data.table)
# library(phyloseq) #there is something weird about using radEmu with phyloseq in the environment. best not to use it at all.
# library(randomcoloR)
# library(cowplot)
# library(plyr)
library(vegan)
library(corncob)
# library(writexl)
# library(readxl)
# library(ggrepel)
if (!require("remotes", quietly = TRUE)){
  install.packages("remotes")
  remotes::install_github("statdivlab/radEmu")
}
library(magrittr)
library(ggplot2)
library(stringr)
library(radEmu)


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


nthreads <- data.table::getDTthreads()
cluster <- multidplyr::new_cluster(n = nthreads)
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


# Read in prepared metadata and phyloseq obejct ---------------------------


temporal_lookup <- list(sampling_time = c("dawn", "peak_photo", NA),
                        sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5", NA),
                        sampling_date = c(20210122, 20210123, 20210124, 20210125, 20210126, NA),
                        site = c("Tektite", "Yawzi", "LB_seagrass", NA),
                        sample_type = c("seawater", "control_extraction", "control_pcr", "control_seq"))



if(file.exists(paste0(projectpath, "/", "usvi_metadata_tidy.txt"))){
  metadata <- readr::read_delim(paste0(projectpath, "/", "usvi_metadata_tidy.txt"),
                                delim = "\t",
                                quote = "",
                                col_names = TRUE,
                                # num_threads = nthreads,
                                show_col_types = FALSE)
} else {
  cli::cli_alert_warning("Please tidy the metadata.")
}


to_import <- c("usvi_prok_asvs.df", "usvi_prok_asvs.taxa")

for(file in to_import){
  if(!exists(file, envir = .GlobalEnv)){
    namevar <- file
    if(file.exists(paste0(projectpath, "/", namevar, ".tsv", ".gz"))){
      cli::cli_alert_info("Importing this dataset: {namevar}")
      temp_df <- readr::read_delim(paste0(projectpath, "/", namevar, ".tsv", ".gz"),
                                   # num_threads = nthreads,
                                   col_names = TRUE, show_col_types = FALSE, delim = "\t")
      assign(paste0(namevar), temp_df, envir = .GlobalEnv)
      rm(temp_df)
      to_import <- to_import[grep(namevar, to_import, value = FALSE, invert = TRUE)]
    } else {
      cli::cli_alert_warning("Please prepare this dataset: {namevar}")
    }
  }
}



# #replace NA in taxonomy with last known level
# 
# usvi_prok_filled.taxa.df <- usvi_prok_asvs.taxa %>%
## usvi_prok_filled.taxa.df <- phyloseq::tax_table(ps_usvi) %>%
#   as.data.frame(.) %>%
#   tibble::as_tibble(., rownames = "asv_id") %>%
#   dplyr::mutate(Phylum = coalesce(Phylum, Domain)) %>%
#   dplyr::mutate(Class = coalesce(Class, Phylum)) %>%
#   dplyr::mutate(Order = coalesce(Order, Class)) %>%
#   dplyr::mutate(Family = coalesce(Family, Order)) %>%
#   dplyr::mutate(Genus = coalesce(Genus, Family)) %>%
#   dplyr::mutate(Species = coalesce(Species, Genus)) %>%
#   dplyr::relocate(asv_id) %>%
#   droplevels


chosen_genera <- c("Endozoicomonas", "Synechococcus CC9902", "Prochlorococcus MIT9313")
chosen_asvs <- usvi_prok_asvs.taxa %>%
  dplyr::filter(Genus %in% chosen_genera) %>%
  droplevels %>%
  dplyr::select(asv_id) %>%
  unlist %>%
  as.character
chosen_samples <- (metadata %>%
                     dplyr::filter(sample_type == "seawater") %>%
                     dplyr::filter(site == "LB_seagrass") %>%
                     droplevels %>%
                     dplyr::select(sample_ID) %>%
                     unlist %>%
                     as.character)

temp_asv_tbl <- usvi_prok_asvs.df %>%
  dplyr::filter(sample_ID %in% chosen_samples) %>%
  dplyr::filter(asv_id %in% chosen_asvs) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "sample_ID",
                     names_from = "asv_id",
                     values_from = "counts") %>%
  tibble::column_to_rownames(var = "sample_ID") %>%
  droplevels
sum(rowSums(temp_asv_tbl) == 0)
sum(colSums(temp_asv_tbl) == 0)
temp_asv_tbl <- temp_asv_tbl %>%
  dplyr::select(!which(colSums(temp_asv_tbl) == 0)) %>%
  droplevels


temp_meta_tbl <- metadata %>%
  dplyr::filter(sample_ID %in% chosen_samples) %>%
  dplyr::mutate(sample_type = factor(sample_type, levels = temporal_lookup[["sample_type"]]),
                sampling_time = factor(sampling_time, levels = temporal_lookup[["sampling_time"]]),
                sampling_date = factor(sampling_date, levels = temporal_lookup[["sampling_date"]]),
                sampling_day = factor(sampling_day, levels = temporal_lookup[["sampling_day"]])) %>%
  dplyr::arrange(match(sample_ID, rownames(temp_asv_tbl))) %>%
  dplyr::select(contains("sampl"), "site") %>%
  tibble::column_to_rownames(var = "sample_ID") %>%
  droplevels

#if we have phyloseq loaded, detach it:
if(any(grepl("phyloseq", search()))){
  cli::cli_alert_danger("You seem to have phyloseq in the environment, which does not play nicely with radEmu!")
  unloadNamespace("phyloseq")
}

# Run radEmu --------------------------------------------------------------



ch_fit <- emuFit(formula = ~ sampling_time,
                 data = temp_meta_tbl,
                 Y = temp_asv_tbl,
                 run_score_tests = FALSE)


ch_df <- ch_fit$coef %>%
  dplyr::left_join(., usvi_prok_asvs.taxa %>%
                     dplyr::select(asv_id, Genus) %>%
                     droplevels,
                   by = join_by("category" == "asv_id")) %>%
  dplyr::mutate(cat_small = paste0(Genus, "_", category)) %>%
  dplyr::mutate(cat_small = factor(cat_small, levels = cat_small[order(Genus)])) %>%
  droplevels


ggplot(ch_df) + 
  geom_point(aes(x = cat_small, y = estimate,color = Genus), size = .5) +
  geom_errorbar(aes(x = cat_small, ymin = lower, ymax = upper, color = Genus), width = .25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Category",
       y = "Estimate") + 
  coord_cartesian(ylim = c(-5,10))

#two genera of Prochlorococcus has different estimates:
#ASV_00007 and ASV_00091

#p-value test for two taxa:

taxa_to_test <- c(grep(paste0(c("ASV_00007", "ASV_00091"), collapse = "|"), colnames(temp_asv_tbl)))
covariate_to_test <- grep("peak_photo", ch_fit$B %>% rownames)
two_robust_score_tests <- emuFit(formula = ~ sampling_time,
                                 data = temp_meta_tbl,
                                 B = ch_fit,
                                 test_kj = data.frame(k = covariate_to_test, 
                                                      j = taxa_to_test), 
                                 Y = temp_asv_tbl)

two_robust_score_tests$coef[taxa_to_test, c("covariate", "category", "estimate", "pval")]

data.frame(counts = temp_asv_tbl[, "ASV_00091"],
           group = temp_meta_tbl$sampling_time) %>%
  mutate(ASV_00091 = counts > 0) %>%
  group_by(group, ASV_00091) %>%
  count()
#this ASV_00091 was found in 3 samples collected during peak_photo, and not at all in Dawn samples

data.frame(counts = temp_asv_tbl[, "ASV_00007"],
           group = temp_meta_tbl$sampling_time) %>%
  mutate(ASV_00007 = counts > 0) %>%
  group_by(group, ASV_00007) %>%
  count()
#whereas ASV_00007 was found in all samples regardless of collection time


# temp_rademu_fit <- radEmu::emuFit(formula = ~ sampling_time, 
#                                   data = temp_meta_tbl,
#                                   Y = temp_asv_tbl,
#                                   run_score_tests = FALSE)
