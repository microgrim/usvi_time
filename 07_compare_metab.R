#07_compare_metab.R


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
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)
library(scales)

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
  cli::cli_alert_warning("Please process the metabolomics sample data previously.")
}
if(file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
  temp_list <- readr::read_rds(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
  list2env(temp_list, envir = .GlobalEnv)
  rm(temp_list)
} else {
  cli::cli_alert_warning("Please process the metabolomics data previously.")
}

drop <- c("CINAR_BC_73")

sample_relabel <- metabolomics_sample_metadata %>%
  dplyr::select(sample_id, site, sampling_day, sampling_time) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  droplevels %>%
  tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time), sep = "_", remove = TRUE)  %>%
  tibble::deframe(.)


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


if(!exists("usvi_sus_metabolites_idx", envir = .GlobalEnv)){
  #here are metabolites where we don't have LODs reported:
  temp_df <- usvi_metabolomics_old_long.df %>%
    dplyr::arrange(LOD) %>%
    dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
    dplyr::arrange(metabolites) %>%
    dplyr::filter(is.na(LOD)) %>%
    droplevels
  
  #we are dropping 4-aminobenozic acid and potentiall glutamic acid
  
  # usvi_metabolomics.unique.summary.df <- usvi_metabolomics_long.df %>%
  #   dplyr::select(metabolites, adaptedDervLabel, concentration, LOD, LOQ) %>%
  #   dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
  #   dplyr::mutate(dummy = dplyr::case_when((10*conc == LOD) ~ 0, 
  #                                          (conc < LOD) ~ 0,
  #                                          (is.na(conc)) ~ NA,
  #                                          .default = 1)) %>%
  #   dplyr::mutate(conc = dplyr::case_when((conc >= LOD) ~ conc,
  #                                         .default = NA)) %>%
  #   dplyr::group_by(simpleName) %>%
  #   dplyr::summarise(unique_obs = length(unique(na.omit(conc))),
  #                    total_obs = sum(dummy, na.rm = TRUE)) %>%
  #   dplyr::arrange((unique_obs)) %>%
  #   dplyr::mutate(simpleName = factor(simpleName, levels = unique(c(temp_df[["metabolites"]], .[["simpleName"]])))) %>%
  #   tidyr::pivot_longer(., cols = !c("simpleName"),
  #                       names_to = "metric",
  #                       values_to = "obs") %>%
  #   droplevels
  # 
  # g_sus <- print(ggplot(data = usvi_metabolomics.unique.summary.df)
  #       + geom_col(aes(x = simpleName, y = obs, fill = metric, group = interaction(simpleName, metric)),
  #                  position = position_dodge2(width = 0.5, preserve = "total", padding = 0.1), color = "grey20")
  #       + theme_bw()
  #       + scale_discrete_manual(aesthetics = "fill", values = c("navy", "limegreen"), labels = c("Total observations", "Unique concentrations"), breaks = c("total_obs", "unique_obs"))
  #       + scale_x_discrete(name = "Metabolite")
  #       + scale_y_continuous(name = "Number of observations", expand = expansion(mult = c(0, 0.1)))
  #       + theme(panel.grid= element_blank(),
  #               strip.text.y = element_text(size = rel(0.7), angle = 0),
  #               strip.text.x = element_text(size = rel(0.7), angle = 0),
  #               legend.position = "bottom",
  #               axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = rel(1)))
  #       # + facet_grid(metric~., scales = "fixed", space = "fixed", drop = TRUE)
  #       )
  # g_sus <- g_sus + patchwork::plot_annotation(title = "Number of observations per metabolite across samples",
  #                                             subtitle = "Metabolites missing reported LODs are to the left")
  # 
  # ggsave(paste0(projectpath, "/", "usvi_sus_metabolite_dist-", Sys.Date(), ".png"),
  #        g_sus,
  #        width = 8, height = 6, units = "in")
  
  
  # usvi_sus_metabolites_idx <- data.frame(metabolites = c("2'deoxyguanosine", "HMP", "adenosine", "inosine", "pyridoxine", "4-aminobenzoic acid"))
  drop_metab <- c("4-aminobenzoic acid")
  usvi_sus_metabolites_idx <- usvi_metabolomics_old_long.df %>%
    dplyr::arrange(LOD) %>%
    dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
    dplyr::arrange(metabolites) %>%
    dplyr::filter(is.na(LOD) | metabolites %in% drop_metab) %>%
    # dplyr::filter(metabolites %in% sus_metab_threshold[["simpleName"]]) %>%
    droplevels
  
  
  # usvi_metab_outlier_summary.df <- usvi_metab_summary.df %>%
  #   dplyr::left_join(., (usvi_metab_summary.df %>%
  #                          dplyr::group_by(site, metabolite) %>%
  #                          dplyr::summarise(site_median = median(concentration, na.rm = TRUE))), by = join_by(site, metabolite)) %>%
  #   dplyr::mutate(rescaled_conc = dplyr::case_when((concentration != 0) ~ log2(concentration/site_median), .default = NA), .by = c(site, metabolite)) %>%
  #   dplyr::arrange(desc(site_median)) %>%
  #   dplyr::mutate(metabolite = factor(metabolite, levels = unique(.[["metabolite"]]))) %>%
  #   droplevels
  # 
  # 
  # usvi_metab_sus_summary.df <- usvi_metabolomics_long.df %>%
  #   dplyr::select(metabolites, adaptedDervLabel, concentration, LOD, LOQ) %>%
  #   dplyr::rename(metabolite = "metabolites", metab_deriv_label = "adaptedDervLabel") %>%
  #   dplyr::filter((metabolite %in% usvi_sus_metabolites_idx[["metabolites"]]) | metabolite == "taurine") %>%
  #   dplyr::left_join(., (metabolomics_sample_metadata %>%
  #                          dplyr::filter(grepl("seawater", sample_type)) %>%
  #                          dplyr::select(sample_id, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
  #                          droplevels),
  #                    by = c("metab_deriv_label")) %>%
  #   dplyr::left_join(., usvi_sus_metabolites_idx %>%
  #                      dplyr::mutate(sus_nonunique = dplyr::case_when((is.na(LODflag) | LODflag == 1) ~ "LOD",
  #                                                                     # (LODflag == TRUE) ~ "questionable bio",
  #                                                                     .default = "questionable bio")) %>%
  #                      dplyr::select(metabolites, sus_nonunique) %>%
  #                      droplevels,
  #                    by = join_by("metabolite" == "metabolites")) %>%
  #   # dplyr::mutate(sus_nonunique = factor(sus_nonunique, levels = unique(.[["sus_nonunique"]]))) %>%
  #   dplyr::mutate(sus_nonunique = dplyr::case_when(is.na(sus_nonunique) ~ "real", .default = sus_nonunique)) %>%
  #   dplyr::arrange(sus_nonunique) %>%
  #   dplyr::mutate(across(c(sample_id, metabolite, sampling_time, sampling_day, site, metab_deriv_label), ~factor(.x))) %>%
  #   # dplyr::left_join(., (usvi_metabolomics_long.df %>%
  #   #                        dplyr::arrange(LOD) %>%
  #   #                        dplyr::ungroup(.) %>%
  #   #                        dplyr::distinct(metabolites, LOD, LOQ, .keep_all = FALSE) %>%
  #   #                        dplyr::rename(metabolite = "metabolites") %>%
  #   #                        # droplevels %>%
  #   #                        # tidyr::pivot_longer(., cols = c("LOD", "LOQ"),
  #   #                        #                     names_to = "limit",
  #   #                        #                     values_to = "threshold") %>%
  #   #                        droplevels),
  #   #                  by = join_by(metabolite), multiple = "all", relationship = "many-to-many") %>%
  #   droplevels
  # 
  # g_usvi_metab_sus <- print(
  #   ggplot(data = usvi_metab_sus_summary.df)
  #   + geom_hline(aes(yintercept = LOD), color = "red", linetype = "dotted")
  #   + geom_hline(aes(yintercept = LOQ), color = "grey", linetype = "dashed")
  #   + geom_point(aes(x = metabolite, y = concentration, fill = sampling_time, shape = site, group = interaction(site, sampling_time)),
  #                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
  #                alpha = 1.0, size = 3)
  #   + theme_bw() 
  #   + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Concentration")
  #   + scale_discrete_manual(aesthetics = c("shape"), values = c(22, 21, 23), labels = site_lookup, breaks = names(site_lookup),
  #                             drop = TRUE)
  #   + scale_discrete_manual(aesthetics = c("fill"), values = sampling_time_colors, labels = names(sampling_time_lookup), 
  #                           breaks = sampling_time_lookup, 
  #                           drop = TRUE)
  #   + facet_wrap(metabolite ~ sus_nonunique, shrink = TRUE, scales = "free")
  # )
  # 
  # ggsave(paste0(projectpath, "/", "g_usvi_metab_sus-", Sys.Date(), ".png"),
  #        g_usvi_metab_sus,
  #        width = 10, height = 10, units = "in")
  
  #not recommended to drop outlier values, but here's a way to determine them:
  {
    #these values of metabolite concentrations are potentially outliers, based on all measured concentrations across sites and times:
    # usvi_sus_metabolites_conc.df <- usvi_metab_summary.df %>%
    #   dplyr::group_by(metabolite) %>%
    #   dplyr::summarise(mean_conc = mean(concentration, na.rm = TRUE),
    #                    conc_5 = quantile(concentration, probs = 0.05, na.rm = TRUE, type = 7),
    #                    conc_95 = quantile(concentration, probs = 0.95, na.rm = TRUE, type = 7),
    #                    outlier1 = outliers::outlier(concentration),
    #                    outlier2 = outliers::outlier(concentration, opposite = TRUE)) %>%
    #   dplyr::mutate(potential_high_outlier = dplyr::case_when((outlier1 > conc_95) ~ outlier1,
    #                                                           (outlier2 > conc_95) ~ outlier2,
    #                                                           .default = NA),
    #                 potential_low_outlier = dplyr::case_when((outlier1 < conc_5) ~ outlier1,
    #                                                          (outlier2 < conc_5) ~ outlier2,
    #                                                          .default = NA)) %>%
    #   dplyr::select(-c(outlier1, outlier2)) %>%
    #   dplyr::rowwise(.) %>%
    #   dplyr::mutate(high_FC = potential_high_outlier/conc_95,
    #                 low_FC = conc_5/potential_low_outlier) %>%
    #   dplyr::mutate(log2high_FC = dplyr::case_when(log2(high_FC) > 1.05 ~ log2(high_FC),
    #                                                .default = NA),
    #                 log2low_FC = dplyr::case_when(log2(low_FC) > 1.05 ~ log2(low_FC),
    #                                               .default = NA)) %>%
    #   droplevels %>%
    #   dplyr::mutate(drop_high = dplyr::case_when(!is.na(log2high_FC) ~ potential_high_outlier,
    #                                              .default = NA),
    #                 drop_low = dplyr::case_when(!is.na(log2low_FC) ~ potential_low_outlier,
    #                                             .default = NA)) %>%
    #   dplyr::select(metabolite, starts_with("drop")) %>%
    #   tidyr::pivot_longer(., cols = starts_with("drop"),
    #                       names_to = "direction",
    #                       values_to = "concentration") %>%
    #   droplevels
    # #drop 19 values:
    # usvi_metab_filtered.df <- dplyr::anti_join(usvi_metab_summary.df, usvi_sus_metabolites_conc.df,
    #                                            by = join_by(metabolite, concentration)) %>%
    #   droplevels
  }
}


# Plot heatmap of metab concentrations ------------------------------------

#also, there are some entries that are dummy values: they are reported as 1/10 the LOD
# the measurement should be removed
usvi_metab_median_list <- usvi_metabolomics_long.df %>%
  dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
  dplyr::rename(metabolite = "metabolites", metab_deriv_label = "adaptedDervLabel") %>%
  dplyr::filter(LODflag == 0) %>%
  dplyr::select(-LODflag) %>%
  dplyr::filter(!(metab_deriv_label %in% drop)) %>%
  dplyr::filter(!(metabolite %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
  dplyr::left_join(., metabolomics_sample_metadata %>%
                     dplyr::distinct(metab_deriv_label, sample_id, sampling_day, site, sampling_time) %>%
                     droplevels, by = join_by(metab_deriv_label)) %>%
  dplyr::filter(!grepl("Day1", sampling_day)) %>%
  droplevels %>%
  dplyr::mutate(across(c(sample_id, metab_deriv_label, metabolite, sampling_day, sampling_time, site), ~factor(.x))) %>%
  dplyr::group_by(sample_id, metab_deriv_label, sampling_day, sampling_time, site) %>%
  tidyr::complete(metabolite) %>%
  dplyr::mutate(concentration = tidyr::replace_na(concentration, 0)) %>%
  dplyr::group_by(metabolite, site, sampling_day) %>%
  dplyr::mutate(sum_conc = sum(concentration, na.rm = TRUE)) %>%
  dplyr::mutate(norm_conc = dplyr::case_when(sum_conc > 0 ~ (concentration/sum_conc),
                                                 .default = 0)) %>%
                                             # .default = NA)) %>%
  dplyr::group_by(metabolite, sampling_day, sampling_time, site) %>%
  dplyr::summarise(mean = mean(concentration, na.rm = TRUE),
                   mean_norm = mean(sum(norm_conc, na.rm = TRUE)),
                   median = median(concentration, na.rm = TRUE),
                   median_norm = median(sum(norm_conc, na.rm = TRUE)),
                   sd = sd(concentration, na.rm = TRUE),
                   .groups = "keep") %>%
  dplyr::mutate(median_label = signif(median, digits = 2),
                median_norm_label = signif(median_norm, digits = 2)) %>%
  dplyr::mutate(across(c(median_label, median_norm_label), ~dplyr::case_when(.x > 0 ~ .x,
                                                                             .default = NA))) %>%
  dplyr::mutate(across(c(median_norm, median, mean, mean_norm), ~dplyr::case_when(.x > 0 ~ .x,
                                                                             .default = NA))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = names(sampling_time_lookup)),
                sampling_day = factor(sampling_day, levels = names(sampling_day_lookup))) %>%
  droplevels %>%
  split(., f = .$site) %>%
  map(., ~.x %>%
        dplyr::group_by(site, metabolite) %>%
        # dplyr::ungroup(.) %>%
        dplyr::mutate(present = sum(median, na.rm = TRUE)) %>%
        tidyr::drop_na(present) %>%
        dplyr::select(-present) %>%
        dplyr::mutate(group_label = dplyr::case_when(grepl("dawn", sampling_time) ~ "_d",
                                                     grepl("photo", sampling_time) ~ "_p",
                                                     .default = NA) %>%
                        paste0(site, .)) %>%
        dplyr::mutate(sampling_time = recode(sampling_time, !!!sampling_time_lookup)) %>%
        dplyr::arrange(sampling_time, sampling_day) %>%
        droplevels)

# # usvi_metab_median_list[[1]] %>%
# #         dplyr::ungroup(.) %>%
# #         dplyr::summarise(present = sum(median, na.rm = FALSE), .by = c(site, metabolite)) %>%
# #         tidyr::drop_na(present)
# #       
temp_df1 <- usvi_metab_median_list[[1]] %>%
  droplevels
# 
# temp_df2 <- temp_df1 %>%
#   # dplyr::group_by(metabolite, group_label) %>%
#   # # dplyr::group_by(metabolite, sampling_day, group_label) %>%
#   # dplyr::mutate(median = median/sum(median, na.rm = TRUE)) %>%
#   dplyr::mutate(median_norm = tidyr::replace_na(median_norm, 0)) %>%
#   droplevels
# 
# temp_breaks2 <- temp_df2 %>% dplyr::ungroup(.) %>% dplyr::select(median_norm) %>% tibble::deframe(.) %>% 
#   quantile(., probs = seq(0, 1, 0.1), names = FALSE, na.rm = TRUE) %>% round(., digits = 1) %>% unique(.)
# temp_breaks2 <- c(temp_breaks2, 0, 1)%>% unique(.) %>% sort(.)
# # temp_pal2 <- colorRampPalette(pals::gnuplot(n = 10)[-1])(length(temp_breaks2))
# temp_pal2 <- colorRampPalette(pals::gnuplot(100))(100)
# temp_nbreaks2 <- (temp_breaks2)/2 %>% unique(.) %>% sort(.)
# 
# print(
#   ggplot(data = temp_df2 %>%
#            dplyr::filter(metabolite == "amMP") %>%
#            droplevels, aes(x = metabolite, y = interaction(sampling_time, sampling_day), group = interaction(group_label)))
#   + theme_bw() 
#   + geom_tile(aes(fill = median_norm), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
#   # + geom_text(aes(label = median_norm_label), size = rel(2))
#   +  scale_fill_gradientn(aesthetics = "fill",
#                           # na.value = NA,
#                           na.value = "white",
#                           colours = temp_pal2, 
#                           breaks = temp_breaks2,
#                           limits = c(NA, NA),
#                           # limits = range(temp_breaks2),
#                           values = temp_breaks2)
#   + scale_x_discrete(name = "Metabolite", 
#                      # labels = usvi_genera_relabel, 
#                      expand = c(0,0), position = "bottom")
#   + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y)
#   + theme(panel.spacing = unit(1, "lines"),
#           panel.background = element_blank(),
#           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
#           axis.text.y = element_text(vjust = 0.5, hjust = 1),
#           panel.grid.major = element_blank(),
#           panel.grid.minor.y = element_blank(),
#           panel.grid.minor.x = element_blank(),
#           panel.ontop = FALSE,
#           strip.text.y = element_blank())
#   + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE, draw.llim = TRUE,
#                                  # + guides(fill = guide_colorbar(nbin = 10, show.limits = TRUE, display = "rectangles", draw.ulim = TRUE, draw.llim = TRUE,
#                                  title = "Normalized \nconcentration \ndensity", direction = "vertical",
#                                  theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5)),
#                                  override.aes = list(stroke = 1, color = "black")),
#            color = "none")
#   + coord_flip()
#   + ggtitle(paste0("Normalized concentration of metabolites in ", namevar2))
# )


# temp_df1 %>%
#   dplyr::group_by(metabolite, sampling_day) %>%
#   dplyr::summarise(sum_norm = sum(median_norm, na.rm = TRUE))

rm(list = apropos("g3_metab_.*", mode = "list"))
for(i in seq_along(usvi_metab_median_list)){
# {i <- 1
  namevar1 <- names(usvi_metab_median_list[i])
  namevar2 <- namevar1 %>%
    gsub("LB_seagrass", "Lameshur Bay", .) %>%
    gsub("peak_photo", "afternoon", .)
  temp_df1 <- usvi_metab_median_list[[i]] %>%
    droplevels
  if(nrow(temp_df1) > 0){
    temp_breaks1 <- temp_df1 %>% dplyr::ungroup(.) %>% dplyr::select(median) %>% tibble::deframe(.)
    temp_pal1 <- colorRampPalette(viridisLite::viridis(100))(100)
    # temp_pal1 <- colorRampPalette(pals::cubicyf(100))(100)
    label_y <- temp_df1 %>% dplyr::ungroup(.) %>% 
      dplyr::distinct(sampling_time, sampling_day) %>%
      dplyr::mutate(label = sampling_day,
                    relabeled = paste0(sampling_time, ".", sampling_day)) %>%
       dplyr::select(label, relabeled) %>%
      tibble::deframe(.)
    
    temp_g3a <- (
      ggplot(data = temp_df1 %>%
               droplevels, aes(x = metabolite, y = interaction(sampling_time, sampling_day), group = group_label))
      + theme_bw()
      + geom_tile(aes(fill = median), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      +  scale_fill_gradientn(aesthetics = "fill",
                              # na.value = "white",
                              na.value = "grey",
                              colours = temp_pal1,
                              transform = T_log10p1(),
                              limits = range(log10p1_br(temp_breaks1)),
                              labels = log10p1_lab_na(temp_breaks1),
                              breaks = log10p1_br(temp_breaks1))
      + scale_x_discrete(name = "Metabolite", 
                         # labels = usvi_genera_relabel, 
                         expand = c(0,0))
      + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y)
      + theme(panel.spacing = unit(1, "lines"),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              strip.text.y = element_blank())
      + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE,  draw.llim = TRUE,
                                     title = paste0("Median \nconcentration \n", "(\U00B5M)"), direction = "vertical",
                                     theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5),
                                                   legend.text.position = "right"),
                                     override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
      + ggtitle(paste0("Concentration of metabolites in ", namevar2))
    )
    
    temp_df2 <- temp_df1 %>%
      # dplyr::mutate(median_norm = tidyr::replace_na(median_norm, 0)) %>%
      droplevels
    
    temp_breaks2 <- temp_df2 %>% dplyr::ungroup(.) %>% dplyr::select(median_norm) %>% tibble::deframe(.) %>% 
      quantile(., probs = seq(0, 1, 0.1), names = FALSE, na.rm = TRUE) %>% round(., digits = 1) %>% unique(.)
    temp_breaks2 <- c(temp_breaks2, 0, 1)%>% unique(.) %>% sort(.)
    temp_pal2 <- colorRampPalette(pals::gnuplot(100))(100)
    temp_nbreaks2 <- (temp_breaks2)/2 %>% unique(.) %>% sort(.)
    
    temp_g3b <- (
      ggplot(data = temp_df2 %>%
               droplevels, aes(x = metabolite, y = interaction(sampling_time, sampling_day), group = interaction(group_label)))
      + theme_bw() 
      + geom_tile(aes(fill = median_norm), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      +  scale_fill_gradientn(aesthetics = "fill",
                              # na.value = "black",
                              na.value = "grey20",
                              limits = range(temp_breaks2),
                              # limits = c(NA, NA),
                              colours = temp_pal2, 
                              breaks = temp_breaks2,
                              values = temp_breaks2)
      + scale_x_discrete(name = "Metabolite", 
                         # labels = usvi_genera_relabel, 
                         expand = c(0,0), position = "bottom")
      + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y)
      + theme(panel.spacing = unit(1, "lines"),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              strip.text.y = element_blank())
      + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE, draw.llim = TRUE,
                                     title = "Daily \nnormalized \nconcentration", direction = "vertical",
                                     theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5)),
                                     override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
      + ggtitle(paste0("Normalized concentration of metabolites in ", namevar2))
    )
    label_y2 <- temp_df2 %>% dplyr::ungroup(.) %>% 
      dplyr::distinct(sampling_time, sampling_day) %>%
      dplyr::mutate(label = sampling_day,
                    relabeled = paste0(sampling_day, ".", sampling_time)) %>%
      dplyr::select(label, relabeled) %>%
      tibble::deframe(.)
    temp_g3c <- (
      ggplot(data = temp_df2 %>%
               droplevels, aes(x = metabolite, y = interaction(sampling_day, sampling_time), group = interaction(group_label)))
      + theme_bw() 
      + geom_tile(aes(fill = median_norm), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      +  scale_fill_gradientn(aesthetics = "fill",
                              # na.value = "black",
                              na.value = "grey20",
                              limits = range(temp_breaks2),
                              # limits = c(NA, NA),
                              colours = temp_pal2, 
                              breaks = temp_breaks2,
                              values = temp_breaks2)
      + scale_x_discrete(name = "Metabolite", 
                         # labels = usvi_genera_relabel, 
                         expand = c(0,0), position = "bottom")
      + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y2)
      + theme(panel.spacing = unit(1, "lines"),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              strip.text.y = element_blank())
      + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE, draw.llim = TRUE,
                                     title = "Daily \nnormalized \nconcentration", direction = "vertical",
                                     theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5)),
                                     override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
    )
    temp_g3d <- (
      ggplot(data = temp_df1 %>%
               droplevels, aes(x = metabolite, y = interaction(sampling_day, sampling_time), group = interaction(group_label)))
      + theme_bw() 
      + geom_tile(aes(fill = median), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      +  scale_fill_gradientn(aesthetics = "fill",
                              # na.value = "white",
                              na.value = "grey",
                              colours = temp_pal1,
                              transform = T_log10p1(),
                              limits = range(log10p1_br(temp_breaks1)),
                              labels = log10p1_lab_na(temp_breaks1),
                              breaks = log10p1_br(temp_breaks1))
      + scale_x_discrete(name = "Metabolite", 
                         # labels = usvi_genera_relabel, 
                         expand = c(0,0))
      + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y2)
      + theme(panel.spacing = unit(1, "lines"),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              strip.text.y = element_blank())
      + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE,  draw.llim = TRUE,
                                     title = paste0("Median \nconcentration \n", "(\U00B5M)"), direction = "vertical",
                                     theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5),
                                                   legend.text.position = "right"),
                                     override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
    )
    temp_g3b <- (temp_g3b + (temp_g3c + theme(axis.text.y.left = element_blank()))) + patchwork::plot_layout(guides = "collect")
    temp_g3a <- (temp_g3a + (temp_g3d + theme(axis.text.y.left = element_blank()))) + patchwork::plot_layout(guides = "collect")
    temp_g3 <- (temp_g3a ) | (temp_g3b)
    assign(paste0("g3_metab_", i, "_", namevar1), temp_g3, envir = .GlobalEnv, inherits = TRUE)
    assign(paste0("g3_metab_", i, "_", namevar1, "_med_relabund"), temp_g3a, envir = .GlobalEnv, inherits = TRUE)
    assign(paste0("g3_metab_", i, "_", namevar1, "_norm_relabund"), temp_g3b, envir = .GlobalEnv, inherits = TRUE)
    
  #   # if(!any(grepl(paste0(namevar1, "relconc", Sys.Date(), collapse = "&"), list.files(projectpath, pattern = "usvi_metab_g3_relconc_.*.png")))){
    # ggsave(paste0(projectpath, "/", "usvi_metab_g3_relconc_", i, "_", namevar1 , "-", Sys.Date(), ".png"),
    #        temp_g3,
    #        width = 20, height = 10, units = "in")
    ggsave(paste0(projectpath, "/", "usvi_metab_g3_relconc_", i, "_", namevar1 , "-", Sys.Date(), ".svg"),
           temp_g3,
           width = 20, height = 10, units = "in")
  #   # }
    
    rm(list = apropos("temp_df.*", mode = "list"))
    rm(list = apropos("namevar.*", mode = "list"))
    rm(list = apropos("temp_breaks.*", mode = "list"))
    rm(list = apropos("temp_labels.*", mode = "list"))
    rm(list = apropos("temp_pal.*", mode = "list"))

  }
}


#these are the SDA metabolites in each site by sampling time:
usvi_sda_metabs <- list(c("homoserine betaine", "cysteine 2", "cysteate"),
                        c("tryptophan", "pantothenic acid", "kynurenine", "homoserine betaine", "histidine", "cysteine 2"),
                        c("uridine", "tryptophan", "tryptamine", "taurine", "homoserine betaine", "guanosine", "glutamine", "cysteine 2", "asparagine", "5'UMP", "5'AMP")) %>%
  setNames(., c("Tektite", "Yawzi", "LB_seagrass"))
usvi_sig_metab_median_list <- imap(usvi_metab_median_list, 
                                   ~.x %>%
                                     dplyr::filter(metabolite %in% usvi_sda_metabs[[.y]]) %>%
                                     droplevels)
usvi_sig_metab_figh <- map(usvi_sda_metabs, ~.x %>% length(.))
    

rm(list = apropos("g4_metab_.*", mode = "list"))
for(i in seq_along(usvi_sig_metab_median_list)){
  # {i <- 3
  namevar1 <- names(usvi_sig_metab_median_list[i])
  namevar2 <- namevar1 %>%
    gsub("LB_seagrass", "Lameshur Bay", .) %>%
    gsub("peak_photo", "afternoon", .)
  temp_df1 <- usvi_sig_metab_median_list[[namevar1]] %>%
    droplevels
  fig_height <- round(usvi_sig_metab_figh[[namevar1]]/2 + 2)
  
  if(nrow(temp_df1) > 0){ #if this dataframe exists and has entries:
    temp_breaks1 <- temp_df1 %>% dplyr::ungroup(.) %>% dplyr::select(median) %>% tibble::deframe(.)
    temp_pal1 <- colorRampPalette(viridisLite::viridis(100))(100)
    label_y <- temp_df1 %>% dplyr::ungroup(.) %>% 
      dplyr::distinct(sampling_time, sampling_day) %>%
      dplyr::mutate(label = sampling_day,
                    relabeled = paste0(sampling_time, ".", sampling_day)) %>%
      dplyr::select(label, relabeled) %>%
      tibble::deframe(.)
    
    temp_df2 <- temp_df1 %>%
      # dplyr::mutate(median_norm = tidyr::replace_na(median_norm, 0)) %>% #this converts any entry for a metabolite where it wasn't measured in that sample, to 0 
      ## but instead, let's just use the na.value option in scale_fill_gradientn
      droplevels
    
    temp_breaks2 <- temp_df2 %>% dplyr::ungroup(.) %>% dplyr::select(median_norm) %>% tibble::deframe(.) %>% 
      quantile(., probs = seq(0, 1, 0.1), names = FALSE, na.rm = TRUE) %>% round(., digits = 1) %>% unique(.)
    temp_breaks2 <- c(temp_breaks2, 0, 1)%>% unique(.) %>% sort(.)
    temp_pal2 <- colorRampPalette(pals::gnuplot(100))(100)
    temp_nbreaks2 <- (temp_breaks2)/2 %>% unique(.) %>% sort(.)
    
    label_y2 <- temp_df2 %>% dplyr::ungroup(.) %>% 
      dplyr::distinct(sampling_time, sampling_day) %>%
      dplyr::mutate(label = sampling_day,
                    relabeled = paste0(sampling_day, ".", sampling_time)) %>%
      dplyr::select(label, relabeled) %>%
      tibble::deframe(.)
    
    temp_g3a <- (
      ggplot(data = temp_df1 %>%
               droplevels, aes(x = metabolite, y = interaction(sampling_time, sampling_day), group = group_label))
      + theme_bw()
      + geom_tile(aes(fill = median), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      +  scale_fill_gradientn(aesthetics = "fill",
                              na.value = "grey",
                              colours = temp_pal1,
                              transform = T_log10p1(),
                              limits = range(log10p1_br(temp_breaks1)),
                              labels = log10p1_lab_na(temp_breaks1),
                              breaks = log10p1_br(temp_breaks1))
      + scale_x_discrete(name = "Metabolite", 
                         # labels = usvi_genera_relabel, 
                         expand = c(0,0))
      + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y)
      + theme(panel.spacing = unit(1, "lines"),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              strip.text.y = element_blank())
      + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE, draw.llim = TRUE,
                                     title = paste0("Median \nconcentration \n", "(\U00B5M)"), direction = "vertical",
                                     theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5)),
                                     override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
      + ggtitle(paste0("Concentration of significant metabolites in ", namevar2))
    )
    
    temp_g3d <- temp_g3a + (aes(x = metabolite, y = interaction(sampling_day, sampling_time), group = interaction(group_label), fill = median)) + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y2)
    
    { 
    # temp_g3d <- (
    #   ggplot(data = temp_df1 %>%
    #            droplevels, aes(x = metabolite, y = interaction(sampling_day, sampling_time), group = interaction(group_label)))
    #   + theme_bw() 
    #   + geom_tile(aes(fill = median), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
    #   +  scale_fill_gradientn(aesthetics = "fill",
    #                           # na.value = "white",
    #                           na.value = "grey",
    #                           colours = temp_pal1,
    #                           transform = T_log10p1(),
    #                           limits = range(log10p1_br(temp_breaks1)),
    #                           labels = log10p1_lab_na(temp_breaks1),
    #                           breaks = log10p1_br(temp_breaks1))
    #   + scale_x_discrete(name = "Metabolite", 
    #                      # labels = usvi_genera_relabel, 
    #                      expand = c(0,0))
    #   + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y2)
    #   + theme(panel.spacing = unit(1, "lines"),
    #           panel.background = element_blank(),
    #           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
    #           axis.text.y = element_text(vjust = 0.5, hjust = 1),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor.y = element_blank(),
    #           panel.grid.minor.x = element_blank(),
    #           panel.ontop = FALSE,
    #           strip.text.y = element_blank())
    #   + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE,  draw.llim = TRUE,
    #                                  title = paste0("Median \nconcentration \n", "(\U00B5M)"), direction = "vertical",
    #                                  theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5),
    #                                                legend.text.position = "right"),
    #                                  override.aes = list(stroke = 1, color = "black")),
    #            color = "none")
    #   + coord_flip()
    # )
    }
    
    
    temp_g3b <- (
      ggplot(data = temp_df2 %>%
               droplevels, aes(x = metabolite, y = interaction(sampling_time, sampling_day), group = interaction(group_label)))
      + theme_bw() 
      + geom_tile(aes(fill = median_norm), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      +  scale_fill_gradientn(aesthetics = "fill",
                              na.value = "black",
                              # na.value = "grey20",
                              limits = range(temp_breaks2),
                              # limits = c(NA, NA),
                              colours = temp_pal2, 
                              breaks = temp_breaks2,
                              values = temp_breaks2)
      + scale_x_discrete(name = "Metabolite", 
                         # labels = usvi_genera_relabel, 
                         expand = c(0,0), position = "bottom")
      + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y)
      + theme(panel.spacing = unit(1, "lines"),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
              axis.text.y = element_text(vjust = 0.5, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              strip.text.y = element_blank())
      + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE, draw.llim = TRUE,
                                     title = "Daily \nnormalized \nconcentration", direction = "vertical",
                                     theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5)),
                                     override.aes = list(stroke = 1, color = "black")),
               color = "none")
      + coord_flip()
      + ggtitle(paste0("Normalized concentration of significant metabolites in ", namevar2))
    )

    temp_g3c <- temp_g3b + (aes(x = metabolite, y = interaction(sampling_day, sampling_time), group = interaction(group_label), fill = median_norm)) + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y2)
    
    {
      # temp_g3c <- (
      #   ggplot(data = temp_df2 %>%
      #            droplevels, aes(x = metabolite, y = interaction(sampling_day, sampling_time), group = interaction(group_label)))
      #   + theme_bw() 
      #   + geom_tile(aes(fill = median_norm), stat = "identity", color = "black", alpha = 1.0, show.legend = TRUE)
      #   +  scale_fill_gradientn(aesthetics = "fill",
      #                           na.value = "black",
      #                           # na.value = "grey20",
      #                           limits = range(temp_breaks2),
      #                           # limits = c(NA, NA),
      #                           colours = temp_pal2, 
      #                           breaks = temp_breaks2,
      #                           values = temp_breaks2)
      #   + scale_x_discrete(name = "Metabolite", 
      #                      # labels = usvi_genera_relabel, 
      #                      expand = c(0,0), position = "bottom")
      #   + scale_y_discrete(name = "Sampling time", expand = c(0,0), labels = label_y2)
      #   + theme(panel.spacing = unit(1, "lines"),
      #           panel.background = element_blank(),
      #           axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
      #           axis.text.y = element_text(vjust = 0.5, hjust = 1),
      #           panel.grid.major = element_blank(),
      #           panel.grid.minor.y = element_blank(),
      #           panel.grid.minor.x = element_blank(),
      #           panel.ontop = FALSE,
      #           strip.text.y = element_blank())
      #   + guides(fill = guide_colorbar(ncol = 1, draw.ulim = TRUE, draw.llim = TRUE,
      #                                  title = "Daily \nnormalized \nconcentration", direction = "vertical",
      #                                  theme = theme(legend.ticks = element_line(color = "black", linewidth = 0.5)),
      #                                  override.aes = list(stroke = 1, color = "black")),
      #            color = "none")
      #   + coord_flip()
      # )
      }

    
    temp_g3b <- (temp_g3b + (temp_g3c + theme(axis.text.y.left = element_blank()))) + patchwork::plot_layout(guides = "collect")
    temp_g3a <- (temp_g3a + (temp_g3d + theme(axis.text.y.left = element_blank()))) + patchwork::plot_layout(guides = "collect")
    temp_g3 <- (temp_g3a ) | (temp_g3b)
    assign(paste0("g4_metab_", i, "_", namevar1), temp_g3, envir = .GlobalEnv, inherits = TRUE)
    assign(paste0("g4_metab_", i, "_", namevar1, "_med_relabund"), temp_g3a, envir = .GlobalEnv, inherits = TRUE)
    assign(paste0("g4_metab_", i, "_", namevar1, "_norm_relabund"), temp_g3b, envir = .GlobalEnv, inherits = TRUE)
    
    if(!any(grepl(paste0(namevar1, "relconc", Sys.Date(), collapse = "&"), list.files(projectpath, pattern = "usvi_metab_g4_sig_relconc_.*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_metab_g4_sig_relconc_", i, "_", namevar1 , "-", Sys.Date(), ".png"),
           temp_g3,
           scale = 1,
           width = 20, height = fig_height, units = "in")
      ggsave(paste0(projectpath, "/", "usvi_metab_g4_sig_relconc_", i, "_", namevar1 , "-", Sys.Date(), ".svg"),
             temp_g3,
             scale = 1,
             width = 20, height = fig_height, units = "in")
    }
    
    rm(list = apropos("temp_df.*", mode = "list"))
    rm(list = apropos("namevar.*", mode = "list"))
    rm(list = apropos("temp_breaks.*", mode = "list"))
    rm(list = apropos("temp_labels.*", mode = "list"))
    rm(list = apropos("temp_pal.*", mode = "list"))
    rm(list = apropos("temp_g3.*", mode = "list"))
  }
    # rm(fig_height)
}



# Other metabolite exploration --------------------------------------------

# usvi_metab_summary.df <- usvi_metabolomics_long.df %>%
#   dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
#   dplyr::bind_rows(., usvi_metab_cinar_bc_73.df %>%
#                      dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
#                      droplevels) %>%
#   dplyr::rename(metabolite = "metabolites", metab_deriv_label = "adaptedDervLabel") %>%
#   dplyr::filter(!(metabolite %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
#   #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
#   dplyr::filter(LODflag == 0) %>%
#   dplyr::left_join(., (metabolomics_sample_metadata %>%
#                          dplyr::filter(grepl("seawater", sample_type)) %>%
#                          dplyr::select(sample_id, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
#                          droplevels),
#                    by = c("metab_deriv_label")) %>%
#   dplyr::select(-LODflag) %>%
#   dplyr::mutate(across(c(sample_id, metabolite, sampling_time, sampling_day, site, metab_deriv_label), ~factor(.x))) %>%
#   droplevels
usvi_metab_filtered.df <- usvi_metabolomics_long.df %>%
  dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
  dplyr::bind_rows(., usvi_metab_cinar_bc_73.df %>%
                     dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
                     droplevels) %>%
  dplyr::rename(metabolite = "metabolites", metab_deriv_label = "adaptedDervLabel") %>%
  dplyr::filter(!(metabolite %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
  #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
  dplyr::filter(LODflag == 0) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                         droplevels),
                   by = c("metab_deriv_label")) %>%
  dplyr::select(-LODflag) %>%
  dplyr::mutate(across(c(sample_id, metabolite, sampling_time, sampling_day, site, metab_deriv_label), ~factor(.x))) %>%
  droplevels

usvi_metab_filtered.df <- usvi_metab_filtered.df %>%
  dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
  dplyr::left_join(., (usvi_metab_filtered.df %>%
                         dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
                         dplyr::group_by(metabolite) %>%
                         dplyr::summarise(med_conc = median(concentration, na.rm = TRUE),
                                          conc_5 = quantile(concentration, probs = 0.05, na.rm = TRUE, type = 7),
                                          conc_95 = quantile(concentration, probs = 0.95, na.rm = TRUE, type = 7)) %>%
                         dplyr::arrange(desc(conc_95)) %>%
                         dplyr::mutate(grouping = c(rep(1:5, each = 10, length.out = length(unique(.[["metabolite"]]))))) %>%
                         dplyr::select(metabolite, grouping)),
                   by = join_by(metabolite), relationship = "many-to-many", multiple = "all")

if(!file.exists(paste0(projectpath, "/", "usvi_metab_colors", ".rds"))){
  metab_colors <- c("#d699cf", "#43bb4e", "#a44ec9", "#74bf3e", "#6161d3", "#adb92f", "#bf3b9e", "#47c581", "#d63b83", "#4f8a29", "#d27ddf", "#9eb553", "#9384e2", "#dcae3d", "#4e62a7", "#e48c26", "#6a94db", "#e95a36", "#3abec8", "#de334d", "#59c9ab", "#be3220", "#4aa5d4", "#be5d27", "#3d9c7e", "#e86cbc", "#82bf71", "#854a98", "#96841e", "#a776b6", "#5f6913", "#cf4160", "#4e965b", "#b6507d", "#2f7038", "#e2828c", "#277257", "#e67c63", "#56642b", "#8e4f73", "#c0af5b", "#9a4555", "#7d8c42", "#a64837", "#a6b174", "#e39155", "#827233", "#d3a06f", "#af7821") %>%
    setNames(., c(usvi_metab_filtered.df %>%
                    dplyr::distinct(metabolite, grouping) %>%
                    dplyr::arrange(grouping) %>%
                    dplyr::select(metabolite) %>%
                    tibble::deframe(.)))
  
  readr::write_rds(metab_colors, paste0(projectpath, "/", "usvi_metab_colors", ".rds"))  
}

#if you want to plot the range of concentrations of the 49 metabolites:
{
  metab_bw_legend <- print(ggplot(data = data.frame(
    (usvi_metab_filtered.df %>%
       dplyr::distinct(metabolite, grouping, sampling_time) %>%
       droplevels)),
    aes(x = 1, y = 1, fill = metabolite, shape = sampling_time)) 
    + geom_blank() 
    + geom_point() 
    + theme_void()
    + guides(fill = guide_legend(order = 1, ncol = 5, title = "Metabolite",  direction = "vertical", override.aes = list(color = "black", shape = 23, size = 4, stroke = 1)),
             shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                  override.aes = list(color = "black", stroke = 1, size = 4)),
             color = "none")
    + theme(legend.background = element_rect(fill = NA, colour = "grey30"))
    + scale_discrete_manual(aesthetics = c("shape"), values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA),
                            drop = TRUE)
    + scale_discrete_manual(aesthetics = c("fill"), values = metab_colors, labels =names(metab_colors), 
                            breaks = names(metab_colors), 
                            drop = TRUE)) %>%
    g_legend()
  
  
  for(i in seq_len(5)){
    metab_grouping <- i
    temp_df <- usvi_metab_filtered.df %>%
      dplyr::filter(grouping == metab_grouping) %>%
      droplevels
    
    temp_g <- (
      ggplot(data = temp_df)
      + geom_boxplot(aes(x = metabolite, y = concentration, group = interaction(metabolite, site, sampling_time)), 
                     color = "black", position = position_dodge2(padding = 0.2, preserve = "single"),
                     outliers = FALSE, show.legend = FALSE)
      + geom_point(aes(x = metabolite, y = concentration, fill = metabolite, group = interaction(metabolite, site, sampling_time), shape = sampling_time), 
                   position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
                   alpha = 1.0, size = 3)
      + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Concentration")
      + scale_fill_manual(values = metab_colors, labels = names(metab_colors), name = "Metabolite")
      + scale_x_discrete(name = "Metabolite")
      + scale_shape_manual(values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
      + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
              panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
              panel.grid = element_blank(),
              axis.title.x = element_blank(),
              legend.position = "bottom",
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
              legend.text = element_text(size = 12, colour = "grey30"))
      + facet_grid(space = "fixed", rows = NULL, cols = vars(site),
                   drop = TRUE, scales = "free", shrink = TRUE,
                   labeller = labeller(site = site_lookup))
      + guides(color = "none", 
               # fill = "none",
               fill = guide_legend(order = 2, ncol = 2, title = "Metabolite", direction = "vertical",
                                   override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
               shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                    override.aes = list(color = "black", stroke = 1, size = 2)))
    )
    assign(paste0("g2_metab_", i), temp_g, envir = .GlobalEnv, inherits = TRUE)
    rm(temp_g)
  }
  
  gpatch <- lapply(apropos("^g2_metab_.$", mode = "list"),
                   get) %>%
    purrr::reduce(., `+`) + 
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position="none")
  
  gpatch_layout <- "
  AADD
  BBEE
  CCFF
"
  
  g2_metab_filtered <- (gpatch + metab_bw_legend) + patchwork::plot_layout(design = gpatch_layout) + patchwork::plot_annotation(title = "Concentrations of metabolites across sites and times",
                                                                                                                                # subtitle = "Removed 19 outlier values and CINAR_BC_73 prior to plotting", 
                                                                                                                                tag_levels = "A")
  
  g2_metab_filtered
}
# if(!any(grepl("metab_filtered", list.files(projectpath, pattern = "usvi_.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_metab_filtered-", Sys.Date(), ".png"),
#          g2_metab_filtered, 
#          width = 16, height = 12, units = "in")
# }

#simple, just one:
{
  # 
  # 
  # g2_mtab_conc <- print(
  #   # ggplot(data = usvi_metab_summary.df)
  #   ggplot(data = usvi_metab_filtered.df %>%
  #            dplyr::filter(grouping == 1) %>%
  #            droplevels)
  #   + geom_boxplot(aes(x = metabolite, y = concentration, group = interaction(metabolite, site, sampling_time)), 
  #                  color = "black", position = position_dodge2(padding = 0.2, preserve = "single"),
  #                  outliers = FALSE, show.legend = FALSE)
  #   + geom_point(aes(x = metabolite, y = concentration, fill = metabolite, group = interaction(metabolite, site, sampling_time), shape = sampling_time), 
  #                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
  #                alpha = 1.0, size = 3)
  #   + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Concentration")
  #   # + scale_fill_manual(values = pals::kelly(10), name = "Metabolite")
  #   + scale_fill_manual(values = metab_colors, labels = names(metab_colors), name = "Metabolite")
  #   + scale_shape_manual(values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
  #   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #           panel.grid = element_blank(),
  #           legend.position = "right",
  #           legend.key = element_blank(),
  #           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  #           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           legend.text = element_text(size = 12, colour = "grey30"))
  #   + facet_grid(space = "fixed", rows = NULL, cols = vars(site),
  #                # + facet_wrap(grouping~site, 
  #                drop = TRUE, scales = "free", shrink = TRUE,
  #                labeller = labeller(site = site_lookup))
  #   # + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #   #              labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  #   + guides(color = "none", 
  #            fill = guide_legend(order = 2, ncol = 2, title = "Metabolite", direction = "vertical",
  #                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #            shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                 override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
}




# Check glutamic acid concentrations --------------------------------------

#tangent, not necessary:
#Melissa determined that glutamic acid and o-acetyl serine derivatized fragments show up as the same peak
#so concentrations of glutamic acid might incorporate measured o-acetyl serine if it is not independently identified in Skyline
#we don't have o-acetyl-serine in this dataset, but...
#do glutamic acid, cysteine, glutathione,and serine trend in any way?

{
  # temp_sus_metab_idx <- c("glutamic acid", "cysteine", "glutathione", "serine")
  # temp_df <- usvi_metab_summary.df %>%
  #   dplyr::filter(grepl(paste0("^", temp_sus_metab_idx, collapse = "|"), metabolite)) %>%
  #   dplyr::distinct(metab_deriv_label, metabolite, .keep_all = TRUE) %>%
  #   dplyr::mutate(across(c(sample_id, metabolite, sampling_time, sampling_day, site, metab_deriv_label), ~factor(.x))) %>%
  #   droplevels
  # temp_tbl <- temp_df %>%
  #   tidyr::pivot_wider(., id_cols = c("metab_deriv_label", "sampling_time", "sampling_day", "site"),
  #                      names_from = "metabolite", 
  #                      values_from = "concentration")
  # temp_lm <- lm(serine ~ `glutamic acid` - 1, temp_tbl)
  # anova(temp_lm)
  # summary(temp_lm) #equivalent to temp_lm$coefficients
  # 
  # lm(`glutamic acid` ~ `serine` - 1, temp_tbl) %>%
  #   summary(.)
  # lm(`cysteine 2` ~ serine - 1, temp_tbl) %>%
  #   summary(.)
  # lm(`glutathione 2` ~ `cysteine 2` - 1, temp_tbl) %>%
  #   summary(.)
  # lm(`glutathione 2` ~ `serine` - 1, temp_tbl) %>%
  #   summary(.)
  # 
  # ggplot(temp_tbl) + geom_point(aes(x = `glutamic acid`, y = serine))
  # 
  # temp_g <- print(
  #   ggplot(data = temp_df)
  #   + geom_point(aes(x = interaction(site,sampling_time, sampling_day), y = concentration, fill = metabolite, group = interaction(site,sampling_time, sampling_day), shape = sampling_time),
  #                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
  #                alpha = 1.0, size = 3)
  #   + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Concentration", transform = "log10")
  #   + scale_shape_manual(values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
  #   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #           panel.grid = element_blank(),
  #           legend.position = "right",
  #           legend.key = element_blank(),
  #           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  #           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           legend.text = element_text(size = 12, colour = "grey30"))
  #   + facet_grid(.~sampling_time, drop = TRUE, scales = "free", space = "fixed",
  #                labeller = labeller(sampling_time = sampling_time_lookup))
  #   + guides(color = "none",
  #            fill = guide_legend(order = 2, ncol = 2, title = "Metabolite", direction = "vertical",
  #                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #            shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                 override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
  # 
}



# NMDS on NPOC/TN profiles ------------------------------------------------

#unnecessary at this time:
{
  # #Brianna determined these were potential outliers in the metabolomes:
  # # LB seagrass peak: 73, 117
  # # LB seagrass dawn: 109,110
  # # Yawzi peak: 67, 69
  # # Yawzi dawn: none
  # # Tektite peak: 98
  # # Tektite dawn: 79
  # 
  # potential_metab_outliers_idx <- c(73, 117, 109, 110, 67, 69, 98, 79) %>%
  #   paste0("CINAR_BC_", .)
  # sample_relabel2 <- metabolomics_sample_metadata %>%
  #   dplyr::select(sample_id, metab_deriv_label, site, sampling_day, sampling_time) %>%
  #   dplyr::distinct(., .keep_all = TRUE) %>%
  #   dplyr::arrange(site, sampling_time, sampling_day) %>%
  #   droplevels %>%
  #   dplyr::select(sample_id, metab_deriv_label) %>%
  #   dplyr::mutate(metab_deriv_label = dplyr::case_when(metab_deriv_label %in% potential_metab_outliers_idx ~ paste0(metab_deriv_label, "*"),
  #                                                      .default = metab_deriv_label)) %>%
  #   tibble::deframe(.)
  # 
  # cinar_npoc_tn.df <- readxl::read_excel(paste0(dirname(projectpath), "/CINAR_NPOCandTNdata.2022.04.06.xlsx")) %>%
  #   dplyr::mutate(doc_label = dplyr::case_when(grepl("CINAR", filename) ~ stringr::str_remove_all(filename, "CINAR_"),
  #                                              .default = filename)) %>%
  #   droplevels
  # 
  # usvi_npoc_tn.df <- ps_usvi %>%
  #   phyloseq::sample_data(.) %>%
  #   tibble::as_tibble(rownames = "sample_id") %>%
  #   dplyr::filter(grepl("seawater", sample_type)) %>%
  #   droplevels %>%
  #   dplyr::select(sample_id, sample_type, site, sampling_date, sampling_time, sampling_day, doc_label) %>%
  #   dplyr::mutate(doc_label = factor(doc_label)) %>%
  #   dplyr::left_join(., metadata %>%
  #                      dplyr::select(sample_id, sample_order, sample_order_all),
  #                    by = join_by(sample_id)) %>%
  #   dplyr::left_join(., cinar_npoc_tn.df %>%
  #                      dplyr::filter(grepl("CINAR", filename)) %>%
  #                      dplyr::select(filename, `NPOC(uM)`, `TN(uM)`, doc_label) %>%
  #                      dplyr::mutate(doc_label = factor(doc_label)) %>%
  #                      droplevels,
  #                    by = join_by(doc_label)) %>%
  #   dplyr::rename(NPOC = "NPOC(uM)",
  #                 TN = "TN(uM)") %>%
  #   tidyr::pivot_longer(., cols = !c(contains("sampl"), site, contains("metab"), filename, doc_label),
  #                       names_to = "nutrient",
  #                       values_to = "concentration") %>%
  #   dplyr::left_join(., metabolomics_sample_metadata %>%
  #                      dplyr::select(sample_id, metab_deriv_label),
  #                    by = join_by(sample_id), relationship = "many-to-many", multiple = "all") %>%
  #   dplyr::distinct(., .keep_all = TRUE) %>%
  #   dplyr::arrange(site, sampling_time, sampling_day) %>%
  #   dplyr::mutate(sample_id = factor(sample_id, levels = unique(.$sample_id)),
  #                 metab_deriv_label = factor(metab_deriv_label, levels = unique(.$metab_deriv_label))) %>%
  #   dplyr::mutate(potential_outlier = dplyr::case_when(metab_deriv_label %in% potential_metab_outliers_idx ~ "*",
  #                                                      .default = NA)) %>%
  #   droplevels
  # 
  # usvi_npoc_tn.tbl <- usvi_npoc_tn.df %>%
  #   dplyr::distinct(sample_id, doc_label, nutrient, .keep_all = TRUE) %>%
  #   tidyr::pivot_wider(., id_cols = c("sample_id", "doc_label"),
  #                      names_from = "nutrient",
  #                      values_from = "concentration") %>%
  #   tidyr::drop_na(.) %>%
  #   dplyr::right_join(., metabolomics_sample_metadata %>%
  #                       dplyr::select(sample_id, metab_deriv_label),
  #                     by = join_by(sample_id), multiple = "all", relationship = "many-to-many") %>%
  #   tidyr::drop_na(.)
  # outliers::outlier(usvi_npoc_tn.tbl[,"NPOC"])
  # outliers::outlier(usvi_npoc_tn.tbl[,"TN"])
  # 
  # min(usvi_npoc_tn.tbl[,"NPOC"])
  # min(usvi_npoc_tn.tbl[,"TN"])
  # 
  # # temp_df <- usvi_npoc_tn.df %>%
  # usvi_npoc_tn.df <- usvi_npoc_tn.df %>%
  #   dplyr::mutate(potential_outlier_y = dplyr::case_when((nutrient == "NPOC" & !is.na(potential_outlier)) ~ signif(min(usvi_npoc_tn.tbl[,"NPOC"]), digits = 0),
  #                                                        (nutrient == "TN" & !is.na(potential_outlier)) ~ signif(min(usvi_npoc_tn.tbl[,"TN"]), digits = 0),
  #                                                        .default = NA)) %>%
  #   droplevels
  # 
  # g_nuts1 <- print(ggplot(data = usvi_npoc_tn.df %>%
  #                           dplyr::filter(!grepl("73|43", metab_deriv_label)) %>%
  #                           droplevels)
  #                  + geom_boxplot(aes(x = interaction(sampling_day, sampling_time), y = concentration, fill = nutrient, group = interaction(sampling_day, sampling_time)))
  #                  + geom_point(aes(x = interaction(sampling_day, sampling_time),
  #                                   y = concentration, fill = nutrient, 
  #                                   group = interaction(sampling_day, sampling_time)), 
  #                               position = position_jitterdodge(jitter.width = 0.1),
  #                               shape = 21)
  #                  + theme_bw()
  #                  + facet_grid(nutrient ~ site, drop = TRUE, scales = "free", space = "free_x",
  #                               labeller = labeller(site = site_lookup))
  #                  + scale_x_discrete(labels = sample_relabel, name = "Sample time")
  #                  + scale_y_continuous(name = "Concentration (uM)")
  #                  + theme(strip.text.y = element_text(angle = 0),
  #                          axis.text.x = element_text(angle = 90), 
  #                          panel.grid.minor = element_blank(),
  #                          panel.grid.major = element_blank())
  # )
  # 
  # g_nuts2 <- print(ggplot(data = usvi_npoc_tn.df %>%
  #                           dplyr::distinct(sample_id, nutrient, .keep_all = TRUE) %>%
  #                           droplevels)
  #                  + geom_boxplot(aes(x = sample_id, y = concentration, fill = sampling_time, group = interaction(sampling_day, sampling_time)))
  #                  + geom_point(aes(x = sample_id, y = concentration, fill = sampling_time, group = interaction(sampling_day, sampling_time)), 
  #                               position = position_jitterdodge(jitter.width = 0.1, seed= 48105),
  #                               shape = 21)
  #                  + theme_bw()
  #                  + geom_text(aes(label = potential_outlier, x = sample_id, y = potential_outlier_y, group = interaction(sampling_day, sampling_time)),
  #                              stat = "identity", position = position_jitter(width = 0.1, height = 0, seed = 48105), size =5,
  #                              colour = "black", fontface = "bold")
  #                  + scale_fill_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                  + facet_grid(nutrient ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                               labeller = labeller(site = site_lookup))
  #                  + scale_x_discrete(labels = sample_relabel2, name = "Sample")
  #                  + scale_y_continuous(name = "Concentration (uM)")
  #                  + theme(strip.text.y = element_text(angle = 0),
  #                          axis.text.x = element_text(angle = 90), 
  #                          panel.grid.minor = element_blank(),
  #                          panel.grid.major = element_blank())
  # )
  # 
  # if(!any(grepl("nutrients", list.files(projectpath, pattern = "usvi_.*.png")))){
  #   ggsave(paste0(projectpath, "/", "usvi_nutrients-", Sys.Date(), ".png"),
  #          g_nuts1, 
  #          width = 10, height = 6, units = "in")
  # }
  # 
}

