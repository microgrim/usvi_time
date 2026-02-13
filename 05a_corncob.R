# 05a_corncob.R

#use corncob instead of radEmu to analyze differentially abundant ASVs and genera

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
library(phyloseq)
library(corncob)
library(splines2)
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

no_progress <- function() {} #your placeholder function for progress reporting

idx_samples <- function(x) grep("samples", names(x), value = TRUE)
coalesce2 <- function(x, y, sep = ".") ifelse(x == y, coalesce(x, y, sep = sep), paste0(x, "_vs_", y))

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

usvi_prok_asvs.taxonomy <- usvi_prok_asvs.taxa %>%
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
  dplyr::arrange(taxonomy) %>%
  dplyr::mutate(across(contains("_id"), ~factor(.x, levels = unique(.x)))) %>%
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


# Global labellers/renamers -----------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_global_labellers.R"))){
  cat("Loading project-wide labellers and lookups.")
  source(paste0(getwd(), "/", "00_global_labellers.R"), local = FALSE,
         echo = FALSE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
} 

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
  contrast = contrast_labels_lookup,
  .multi_line = TRUE,
  .default = label_wrap_gen(25, multi_line = TRUE)
  # .default = label_value
)




# Import phyloseq object --------------------------------------------------

drop <- c("CINAR_BC_73")

if(!exists("ps_usvi", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))){
    # ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))
    ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))
  } else {
    cli::cli_alert_warning("Please process the USVI data through Phyloseq.")
    
  }
}

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

ps_usvi_filtered <- ps_usvi %>%
  phyloseq::prune_samples(rownames(phyloseq::sample_data(ps_usvi)) %in% rownames(usvi_selected_metadata), .) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove ASVs not present in any samples

# #let's try consolidating the ASV table to species- or genus-level first
# usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
#   dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus) %>%
#   dplyr::distinct(Domain, Phylum, Class, Order, Family, Genus, .keep_all = TRUE) %>%
#   dplyr::left_join(., usvi_prok_asvs.taxonomy %>%
#                      dplyr::select(asv_id, taxonomy) %>%
#                      droplevels,
#                    by = join_by(asv_id)) %>%
#   droplevels
# 
# 
# usvi_sw_genus.tbl <- phyloseq::otu_table(ps_usvi_filtered) %>%
#   as.data.frame(.) %>%
#   tibble::as_tibble(rownames = "asv_id") %>%
#   otu_to_taxonomy(., usvi_prok_filled.taxa.df, level = "Genus") %>%
#   dplyr::left_join(., usvi_sw_genus.taxa.df,
#                    by = join_by(Domain, Phylum, Class, Order, Family, Genus)) %>%
#   dplyr::relocate(asv_id) %>%
#   dplyr::select(-c(Domain, Phylum, Class, Order, Family, Genus, taxonomy)) %>%
#   droplevels  %>%
#   tibble::column_to_rownames(var = "asv_id") %>%
#   dplyr::slice(which(rowSums(.) > 0)) %>%
#   as.data.frame(.)
# 
# ps_usvi_filtered_agglom <- phyloseq::phyloseq(phyloseq::otu_table(usvi_sw_genus.tbl,
#                                                   taxa_are_rows=TRUE),
#                               phyloseq::sample_data(usvi_selected_metadata),
#                               phyloseq::tax_table(usvi_sw_genus.taxa.df %>%
#                                                     dplyr::filter(asv_id %in% rownames(usvi_sw_genus.tbl)) %>%
#                                                     droplevels %>%
#                                                     tibble::column_to_rownames(var = "asv_id") %>%
#                                                     as.matrix))


# #not a phyloseq object:
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





# Prepare for corncob on phyloseq objects ---------------------------------


##summary of setting up corncob:
#corncob implements FDR/Benjamini-Hochberg adjustment to p-values
#however, many of those taxa that are "significant" have one or both of:
##A. really large standard errors in their estimates
##B. estimates centered around 0 (so not very substantial)
#for example, corncob found 114 significantly differentially abundant genera btw Yawzi, Tektite, and LB at peak photo
#but over 80 of them had abs(Estimate) < 2, and several had StdErrors ~ abs(Estimate)
#so it's more rigorous to use those p-values from corncob, perform a q-value adjustment with ashr::qval.from.lfdr, and filter the p-values that pass
#in this example, we go from 114 "significant" genera to 21

#try out testing for differential abundance and dispersion across sites, 
#controlling for the effect of sampling_day on dispersion and of site on abundance:
#first, compare using samples from just peak_photo in Yawzi, tektite, and LB.
{
  temp_ps <- phyloseq::prune_samples(rownames(phyloseq::sample_data(ps_usvi_filtered_agglom)) %in% 
                                       (usvi_selected_metadata %>%
                                          dplyr::filter(grepl("peak", sampling_time)) %>%
                                          rownames()), ps_usvi_filtered_agglom)
  
  temp_res_cc <- try(corncob::differentialTest(formula = ~ site + sampling_day,
                                               phi.formula = ~ site + sampling_day,
                                               formula_null = ~ site,
                                               phi.formula_null = ~ sampling_day,
                                               test = "Wald",
                                               data = temp_ps,
                                               fdr = "fdr",
                                               fdr_cutoff = set.alpha), silent = TRUE)
  # plot(temp_res_cc, level = "Genus")
  # summary(temp_res_cc)
  #pull out the coefficients for each significant taxon in each site or sampling day:
  
  
  temp_res_cc_list <- NULL
  # for(i in c(1:2)){
  for(i in seq_len(length(temp_res_cc[["significant_taxa"]]))){
    temp_df <- as.data.frame(temp_res_cc[["significant_models"]][[i]][["coefficients"]]) %>%
      tibble::as_tibble(rownames = "condition") %>%
      dplyr::mutate(asv_id = temp_res_cc[["significant_taxa"]][i])
    temp_res_cc_list <- bind_rows(temp_res_cc_list, temp_df)
  }
  
  
  
  
  #so the corncob test implements a fdr_cutoff, but is it fair?
  q_value <- 0.10
  padj_cutoff <- temp_res_cc[["p"]] %>%
    p.adjust(., method = "BH") %>% #multiple testing corrections, "BH" is an alias for "fdr" accordiong to p.adjust()
    na.omit(.) %>%
    ashr::qval.from.lfdr(.) %>%
    unlist %>%
    quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)
  
  temp_df <- tibble::enframe(temp_res_cc[["p"]], name = "asv_id", value = "p") %>%
    dplyr::mutate(p_bh_adj = dplyr::case_when(p <= padj_cutoff ~ p,
                                              .default = NA)) %>%
    dplyr::left_join(., tibble::enframe(temp_res_cc[["p_fdr"]], name = "asv_id", value = "p_fdr")) %>%
    dplyr::mutate(p_cc = dplyr::case_when(p_fdr <= set.alpha ~ p_fdr,
                                          .default = NA))
  
  
  
  
  #using corncob's implementation of multiple test correction and p < 0.05, we have 114 significant genera
  
  #using a BH post-hoc test correction
  #for q = 0.1 on p, 21 significant genera
  #for q = 0.05 on p, 12 significant genera
  #for q = 0.1 on the p_fdr, we have 17 significant genera
  #for q = 0.05 on the p_fdr, we have 9 significant genera
  
  #these taxa are included in the list of 114 identified in corncob as significant; there aren't any excluded
  
  temp_res_cc_df <- temp_res_cc_list %>%
    dplyr::rename(StdError = `Std. Error`) %>%
    dplyr::left_join(., tibble::enframe(temp_res_cc[["p"]], name = "asv_id", value = "p") %>%
                       dplyr::mutate(p_bh_adj = dplyr::case_when(p <= padj_cutoff ~ p,
                                                                 .default = NA)) %>%
                       dplyr::left_join(., tibble::enframe(temp_res_cc[["p_fdr"]], name = "asv_id", value = "p_fdr")) %>%
                       dplyr::mutate(p_cc = dplyr::case_when(p_fdr <= set.alpha ~ p_fdr,
                                                             .default = NA)) %>%
                       dplyr::select(asv_id, p_bh_adj, p_cc),
                     by = join_by(asv_id)) %>%
    dplyr::mutate(metric = dplyr::case_when(grepl("^mu", condition) ~ "mu",
                                            grepl("^phi", condition) ~ "phi",
                                            .default = NA)) %>%
    dplyr::mutate(variable = stringr::str_remove_all(condition, paste0(metric, "."))) %>%
    dplyr::mutate(variable = dplyr::case_when(grepl("Intercept", variable) ~ "LB_seagrass",
                                              grepl("site", variable) ~ stringr::str_remove_all(variable, "site"),
                                              grepl("sampling_day", variable) ~ stringr::str_remove_all(variable, "sampling_day"),
                                              .default = variable))
  
  #plot the coefficient estimates and std.errors for these 114 taxa
  temp_df_to_plot <- temp_res_cc_df %>%
    dplyr::left_join(., (temp_res_cc_df %>%
                           dplyr::filter(grepl("LB", variable)) %>%
                           dplyr::rename(Estimate_LB = Estimate) %>%
                           dplyr::select(Estimate_LB, asv_id, metric) %>%
                           droplevels),
                     by = join_by(asv_id, metric), relationship = "many-to-many", multiple = "all") %>%
    dplyr::filter(metric == "mu") %>%
    # dplyr::mutate(Estimate = Estimate - Estimate_LB) %>%
    # dplyr::filter(!grepl("LB", variable)) %>%
    dplyr::arrange(desc(Estimate)) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.$asv_id))) %>%
    droplevels
  g <- print(
    ggplot(data = temp_df_to_plot)
    + geom_errorbar(aes(x = asv_id, ymin = Estimate-StdError, ymax = Estimate+StdError), color = "black", width = .3, position=position_dodge(.9)) 
    + geom_point(aes(x = asv_id, y = Estimate, group = interaction(asv_id, variable)), fill = "grey", alpha = 0.5,
                 size = 1, pch = 21) 
    + geom_point(data = temp_df_to_plot %>%
                   tidyr::drop_na(p_bh_adj) %>%
                   droplevels,
                 aes(x = asv_id, y = Estimate, fill = variable, group = interaction(asv_id, variable)),
                 size = 4, pch = 21) 
    + scale_x_discrete(limits = rev(levels(temp_df_to_plot$asv_id))) 
    # + scale_fill_manual(values = palette) 
    # + geom_line() 
    + coord_flip()
    + theme_bw()
    + labs(x = "Taxa", y = "Coefficient", title = "Peak photosynthesis")
    + geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")
    + facet_grid(.~variable)
  )
  
  #for q = 0.1 on p-values, 21 significant genera using ashr
  temp_df2 <- temp_res_cc_df %>%
    dplyr::filter(condition == unique(temp_res_cc_list$condition)[2]) %>%
    dplyr::select(asv_id, Estimate, StdError, p_bh_adj, p_cc) %>%
    dplyr::arrange(Estimate) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
    tidyr::pivot_longer(., cols = starts_with("p_"),
                        names_to = "parameter",
                        values_to = "value") %>%
    dplyr::mutate(Estimate = dplyr::case_when(!is.na(value) ~ Estimate,
                                              .default = NA),
                  StdError = dplyr::case_when(!is.na(value) ~ StdError,
                                              .default = NA)) %>%
    dplyr::mutate(parameter = factor(parameter, levels = c("p_cc", "p_bh_adj"), labels = c("corncob's adjusted p-value", "ashr adjusted p-value"))) %>%
    droplevels 
  
  temp_list <- temp_df2 %>%
    split(., f = .$parameter) %>%
    map(., ~.x %>%
          tidyr::drop_na(value) %>%
          dplyr::select(Estimate) %>%
          tibble::deframe(.))
  
  t.test(temp_list[[1]], temp_list[[2]])
  kruskal.test(temp_list)
  
  
  temp_g1 <- print(ggplot(data = temp_df2)
                   + theme_bw()
                   + geom_histogram(aes(x = Estimate, fill = parameter), color = "black", bins = 10)
                   + scale_y_continuous(name = "Number of coefficients")
                   + facet_wrap(.~parameter)
                   + guides(fill = "none")
  )
  temp_g2 <- print(ggplot(data = temp_df2)
                   + theme_bw()
                   + geom_errorbar(aes(x = asv_id, ymin = Estimate-StdError, ymax = Estimate+StdError), color = "black", width = .3, position=position_dodge(.9)) 
                   + geom_point(aes(x = asv_id, y = Estimate, fill = parameter), shape = 21)
                   + coord_flip()
                   + guides(fill = "none")
                   # + scale_y_continuous(name = "Number of coefficients")
                   + facet_wrap(.~parameter)
  )
  
  temp_g3 <- print(ggplot(data = temp_df2)
                   + theme_bw()
                   + geom_histogram(aes(x = StdError, fill = parameter), color = "black", bins = 10)
                   + scale_y_continuous(name = "Number of coefficients")
                   + facet_wrap(.~parameter)
                   + guides(fill = "none")
  )
  
  temp_g <- temp_g1 / temp_g2
  if(!file.exists(paste0(projectpath, "/", "modeled_estimates_from_corncob", ".png"))){
    ggsave(paste0(projectpath, "/", "modeled_estimates_from_corncob", ".png"),
           temp_g,
           width = 10, height = 10, units = "in")  
  }
  
}





# Full scale corncob ------------------------------------------------------


#loop to handle agglomerated and not-agglomerated table
ps_objects <- c("ps_usvi_filtered")
# ps_objects <- c("ps_usvi_filtered", "ps_usvi_filtered_agglom")

model_design <- rlang::parse_expr("~ 0 + site:sampling_time")

extant_cc_dt_files <- list.files(path = projectpath, pattern = "^cc_dt_res_")

for(i in ps_objects){
  contrast_vars <- stringr::str_remove_all(deparse(model_design), "(\\d)") %>% 
    stringr::str_remove_all(., "(\\s)|(\\~)") %>% #remove trailing tilde and spaces
    stringr::str_remove_all(., "^((\\:)|(\\+)|(\\*))*") %>% #remove any beginning punctuation as part of the formula 
    stringr::str_split(., pattern = "(\\:)|(\\+)|(\\*)") %>%
    simplify(.)
  namevar <- paste0(c(i, contrast_vars), collapse = "-")
  if(!exists(paste0("cc_dt_res_", namevar, "_list"), inherits = TRUE)){
    
    dt_list <- grep(namevar, extant_cc_dt_files, value = TRUE) %>%
      grep("list", ., value = TRUE)
    
    if(!(purrr::is_empty(dt_list)) & file.exists(paste0(projectpath, "/", dt_list))){
      cli::cli_alert_success(paste0("Reading in extant results list for corncob::differentialTest for this PhyloSeq object: ", namevar), wrap = TRUE)
      progressr::with_progress({
        cc_dt_results_list <- readr::read_rds(file = paste0(gprojectpath, "/", dt_list))
      })
      assign(paste0("cc_dt_res_", namevar, "_list"), cc_dt_results_list, envir = .GlobalEnv, inherits = TRUE)
      rm(cc_dt_results_list)
    } else {
      cli::cli_alert_info(paste0("Running corncob::differentialTest on this PhyloSeq object: ", namevar), wrap = TRUE)
      ps_list <- get0(i, inherits = TRUE) #get the phyloseq object
      
      #format metadata to make the model matrix
      metadata <- ps_list %>%
        phyloseq::sample_data(.) %>%
        t() %>%
        as.data.frame() %>%
        tidyr::drop_na() %>%
        t() %>%
        as.data.frame() %>%
        droplevels %>%
        dplyr::mutate(across(any_of(contrast_vars), factor)) %>%
        dplyr::mutate(site = factor(site, levels = names(site_lookup)),
                      sampling_time = factor(sampling_time, levels = names(sampling_time_lookup))) %>%
        droplevels %>%
        dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                      sampling_time = fct_relevel(sampling_time, "dawn")) %>%
        dplyr::arrange(across(contrast_vars))
      
      #check the model matrix is full
      {
        mod_mat <- stats::model.matrix(eval(model_design),
                                       data = metadata)
        
        if(any(colSums(mod_mat) == 0)){ #none of these should be 0--if so, you need to fix your design/metadata coding
          idx <- names(which(colSums(mod_mat) == 0))
          idx <- cli::cli_vec(idx)
          cli::cli_alert_warning("One or more groups in the comparison does not have enough replicates:", wrap = TRUE)
          cli::cli_bullets_raw(idx)
          cli::cli_abort("Ending modeling attempt.")
        }
        
        if(exists("contrast_id")){
          rm(contrast_id)
        }
        }
        contrast_names <- colnames(mod_mat)
        
        #provide covariates to build pairs of comparisons
        {
          covariates <- metadata %>%
                dplyr::select(all_of(contrast_vars)) %>%
            dplyr::arrange(across(contrast_vars)) %>%
                dplyr::distinct(.) %>%
                simplify %>%
                as.character %>%
            unique(.)
          for(item in covariates){
            temp_contrast <- grep(item, contrast_names, value = TRUE, ignore.case = TRUE)  %>% as_tibble_row(., .name_repair = "unique")
            if(exists("contrast_id", inherits = TRUE)){
              contrast_id <- bind_rows(contrast_id, temp_contrast)
            } else {
              contrast_id <- temp_contrast
            }
          }

          contrast_id_df <- apply(contrast_id, 1, function(x) {
            x %>%
              t() %>%
              rev() %>% #this set up the contrast of "Tektite - Yawzi" when really I wanted "Yawzi - Tektite"
              combn(., 2) %>%
              t() %>%
              as.data.frame() %>%
              setNames(., c("pair1", "pair2")) %>%
              tidyr::unite(contrast, sep = " - ", remove = FALSE) %>%
              droplevels
          }) %>%
            setNames(., covariates) %>%
            bind_rows(., .id = "hold1") %>%
            dplyr::mutate(label_id = contrast) %>%
            dplyr::mutate(label_id = stringr::str_remove_all(label_id, paste0("[[:alpha:]_]{1,}", .$hold1)) %>%
                            stringr::str_remove_all(., paste0(c(contrast_vars, ":"), collapse = "|"))) %>%
            dplyr::mutate(label_id = paste0(.$hold1, " (", label_id) %>%
                            paste0(., ")")) %>%
            dplyr::mutate(var1 = contrast) %>%
            dplyr::mutate(var1 = stringr::str_remove_all(var1, paste0("[[:alpha:]_]{1,}", .$hold1)) %>%
                            stringr::str_remove_all(., paste0(c(":", covariates), collapse = "|")) %>%
                            stringr::str_remove_all(., " - [[:alpha:]_]{1,}$")) %>%
            dplyr::relocate(label_id) %>%
            tidyr::drop_na(pair1, pair2)
            # tidyr::drop_na(pair1)
        }
        
      
      ##iterate through the length of the contrast_id_df dataframe, going through pairs of sample groups
        
      rm(list = apropos("(^cc_dt_res_)([0-9])+$", mode = "list"))
      for(k in 1:length(contrast_id_df$label_id)){
        contrast_k <- contrast_id_df[k, "label_id"]
        cli::cli_progress_message("Processing contrast for differentialTest.")
        cli::cli_progress_step("Working on dataset {k}: {contrast_k}", spinner = TRUE)
        # #errored on 4th dataset: peak_photo (Yawzi - Tektite)
        # for(k in 4:5){ #for troubleshooting:
        ##remove the samples from the dataset that aren't in this iteration

        idx <- mod_mat[, (contrast_id_df[k,] %>% dplyr::select(starts_with("pair")) %>%
                            simplify %>%
                            stringr::str_remove_all(., "(\\(|\\)|\\-)") %>%
                            stringr::str_split(., "(\\s)+") %>%
                            simplify %>%
                            unique)]
        idx <- idx[which(rowSums(idx) > 0),]
        idx <- rownames(idx)
        temp_ps <- phyloseq::prune_samples(idx, ps_list)
        temp_list <- list()

        if(any(grepl("0", model_design))){ #tidy up the formula to remove 0
          # model_design_ready <- rlang::parse_expr(stringr::str_remove_all(deparse(model_design), "(\\d)|(\\+)|(\\w*):") %>% stringr::str_squish(.)) #use only the final covariate in the differentialTest
          model_design_ready <- contrast_id_df[k,] %>% dplyr::select(starts_with("var1")) %>% simplify %>% unique(.) %>% #use only the final covariate in the differentialTest
            paste0("~", .)
          model_design_ready <- rlang::parse_expr(model_design_ready)

        }

        tryCatch({
          temp_res_da <- try(corncob::differentialTest(formula = eval(model_design_ready),
                                                       phi.formula = eval(model_design_ready),
                                                       formula_null = ~ 1,
                                                       phi.formula_null = eval(model_design_ready),
                                                       test = "Wald",
                                                       data = temp_ps,
                                                       fdr = "fdr",
                                                       fdr_cutoff = set.alpha), silent = TRUE)
          temp_res_dv <- try(corncob::differentialTest(formula = eval(model_design_ready),
                                                       phi.formula = eval(model_design_ready),
                                                       formula_null = eval(model_design_ready),
                                                       phi.formula_null = ~ 1,
                                                       test = "Wald",
                                                       data = temp_ps,
                                                       fdr = "fdr",
                                                       fdr_cutoff = set.alpha), silent = TRUE)

          if ( !inherits(temp_res_da, "try-error") & !inherits(temp_res_dv, "try-error") ) {
            temp_list <- list(temp_res_da, temp_res_dv) %>%
              setNames(., c("cc_da", "cc_dv"))
          } else if (!inherits(temp_res_da, "try-error")){
            temp_list <- list(temp_res_da) %>%
              setNames(., c("cc_da"))
          } else if ( !inherits(temp_res_dv, "try-error")){
            temp_list <- list(temp_res_dv) %>%
              setNames(., c("cc_dv"))
          }
        }, error = function(e){
          warning("differentialTest did not have enough variables", call. = FALSE)
        })

        assign(paste0("cc_dt_res_", k), temp_list, envir = .GlobalEnv, inherits = TRUE)
        rm(temp_list)
        rm(temp_res_da)
        rm(temp_res_dv)
        cli::cli_progress_update()
      }
      cli::cli_progress_done()

      cc_dt_results_list <- apropos("(^cc_dt_res_)([0-9])+$", mode = "list") %>%
        lapply(., get0) %>%
        setNames(., contrast_id_df[["label_id"]])
      if(!exists(paste0("cc_dt_res_", namevar, "_list"), inherits = TRUE)){
        assign(paste0("cc_dt_res_", namevar, "_list"), cc_dt_results_list, envir = .GlobalEnv, inherits = TRUE)
        if(!any(grepl(namevar, list.files(projectpath, pattern = "^cc_dt_res_.*_list-.*.rds")))){
          readr::write_rds(cc_dt_results_list, file = paste0(projectpath, "/", "cc_dt_res_", namevar, "_list-", Sys.Date(), ".rds"))
        }
        rm(cc_dt_results_list)
      }
      rm(list = apropos("(^cc_dt_res_)([0-9])+$", mode = "list"))
    }
  } else {
    cli::cli_alert_danger("differentialTest results list already exists for this PhyloSeq object:", wrap = TRUE)
    cli::cli_bullets_raw(namevar)
  }
}

#want to see one plotted?

# plot(`cc_dt_res_ps_usvi_filtered_agglom-site-sampling_time_list`[["Yawzi (peak_photo - dawn)"]][["cc_da"]], level = "Genus")


#let's plot the p-values and their adjusted
{
  temp_res_cc <- `cc_dt_res_ps_usvi_filtered-site-sampling_time_list`[["Yawzi (peak_photo - dawn)"]][["cc_da"]]
  
  
  # temp_df1 <- `cc_dt_res_ps_usvi_filtered-site-sampling_time_list`[["Yawzi (peak_photo - dawn)"]][["cc_da"]][["p"]] %>%
  #   tibble::enframe(name = "asv_id", value = "p_val") %>%
  #   dplyr::mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
  #   tidyr::drop_na(.) %>%
  #   dplyr::mutate(q_val = ashr::qval.from.lfdr((p_val)),
  #                 q_adj = ashr::qval.from.lfdr((p_adj))) %>%
  #   droplevels
  # temp_df <- temp_df1 %>%
  #   dplyr::select(- starts_with("q_")) %>%
  #   tidyr::pivot_longer(., cols = !c("asv_id"),
  #                       names_to = "parameter",
  #                       values_to = "value") %>%
  #   dplyr::mutate(type = "p_values") %>%
  #   bind_rows(., (temp_df1 %>%
  #                   dplyr::select(asv_id,starts_with("q_")) %>%
  #                   tidyr::pivot_longer(., cols = !c("asv_id"),
  #                                       names_to = "parameter",
  #                                       values_to = "value") %>%
  #                   dplyr::mutate(type = "q_values") %>%
  #                   droplevels)) %>%
  #   dplyr::mutate(parameter = factor(parameter, levels = c("p_val", "p_adj", "q_val", "q_adj"))) %>%
  #   droplevels
  # temp_df_labels <- temp_df %>%
  #   dplyr::group_by(parameter) %>%
  #   dplyr::reframe(quant_05 = quantile(value, probs = 0.05, na.rm = TRUE),
  #                  quant_10 = quantile(value, probs = 0.10, na.rm = TRUE)) %>%
  #   tidyr::pivot_longer(., cols = !c("parameter"),
  #                       names_to = "quantile",
  #                       values_to = "value")
  
  temp_res_cc_list <- NULL
  # for(i in c(1:2)){
  for(i in seq_len(length(temp_res_cc[["significant_taxa"]]))){
    temp_df <- as.data.frame(temp_res_cc[["significant_models"]][[i]][["coefficients"]]) %>%
      tibble::as_tibble(rownames = "condition") %>%
      dplyr::mutate(asv_id = temp_res_cc[["significant_taxa"]][i])
    temp_res_cc_list <- bind_rows(temp_res_cc_list, temp_df)
  }
  
  temp_df1 <- temp_res_cc[["p"]] %>%
    tibble::enframe(name = "asv_id", value = "p_val")
  
  padj_cutoff <- temp_df1[["p_val"]] %>%
    # p.adjust(., method = "BH") %>% #multiple testing corrections, "BH" is an alias for "fdr" accordiong to p.adjust()
    na.omit(.) %>%
    ashr::qval.from.lfdr(.) %>%
    unlist %>%
    quantile(., probs = c(0.05, 0.10, 0.01), na.rm = TRUE, names = FALSE,type = 7)
  #for q = 0.10, we have a threshold of p-value = 0.00733
  #for q = 0.05, we have threshold of p-value = 0.000753
  #for q = 0.01, we have a threshold of p-value = 0.00000111
  
  
  temp_df2 <- temp_res_cc[["p_fdr"]] %>%
    tibble::enframe(name = "asv_id", value = "p_fdr")
  temp_df2 <- temp_df1 %>%
    dplyr::mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    dplyr::mutate(p_bh_05 = dplyr::case_when(p_val <= padj_cutoff[1] ~ p_val,
                                             .default = NA)) %>%
    dplyr::mutate(p_bh_10 = dplyr::case_when(p_val <= padj_cutoff[2] ~ p_val,
                                             .default = NA)) %>%
    dplyr::mutate(p_bh_01 = dplyr::case_when(p_val <= padj_cutoff[3] ~ p_val,
                                             .default = NA)) %>%
    dplyr::left_join(., temp_df2, by = join_by(asv_id)) %>%
    dplyr::mutate(p_cc = dplyr::case_when(p_fdr <= set.alpha ~ p_fdr,
                                          .default = NA)) %>%
    dplyr::mutate(p_val = dplyr::case_when(p_val <= set.alpha ~ p_val,
                                           .default = NA)) %>%
    dplyr::select(asv_id, p_val, p_adj, p_cc, starts_with("p_bh"))
  
  temp_df <- temp_df2 %>%
    tidyr::pivot_longer(., cols = starts_with("p_"),
                        names_to = "parameter",
                        values_to = "value")
  # temp_df %>%
  #   dplyr::group_by(parameter) %>%
  #   dplyr::summarise(num_SDA = length(na.omit(value)))
  #669 results from corncob::differentialTest on Yawzi temporal
  #ASVs where p-value < 0.05: 73
  #corncob-adjusted p-value < 0.05: 30
  #using q-value = 0.10 on raw p-values to determine the p-value threshold for significance: 47
  #using q-value = 0.05 on raw p-values to determine the p-value threshold for significance: 21
  #using q-value = 0.01 on raw p-values to determine the p-value threshold for significance: 5
  
  
  temp_res_cc_df <- temp_res_cc_list %>%
    dplyr::rename(StdError = `Std. Error`) %>%
    dplyr::left_join(., temp_df2 %>%
                       dplyr::select(asv_id, p_cc, starts_with("p_bh")),
                     by = join_by(asv_id)) %>%
    dplyr::mutate(metric = dplyr::case_when(grepl("^mu", condition) ~ "mu",
                                            grepl("^phi", condition) ~ "phi",
                                            .default = NA)) %>%
    dplyr::mutate(variable = stringr::str_remove_all(condition, paste0(metric, "."))) %>%
    dplyr::mutate(variable = dplyr::case_when(grepl("Intercept", variable) ~ "dawn",
                                              grepl("site", variable) ~ stringr::str_remove_all(variable, "site"),
                                              grepl("sampling_day", variable) ~ stringr::str_remove_all(variable, "sampling_day"),
                                              grepl("sampling_time", variable) ~ stringr::str_remove_all(variable, "sampling_time"),
                                              .default = variable)) %>%
    # dplyr::mutate(variable = paste0("Yawzi_", variable)) %>%
    droplevels
  
  #plot the coefficient estimates and std.errors for these 114 taxa
  temp_df_to_plot <- temp_res_cc_df %>%
    dplyr::filter(metric == "mu") %>%
    dplyr::arrange(desc(Estimate)) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.$asv_id))) %>%
    droplevels
  
  temp_res_cc_summary.df <- temp_df_to_plot %>%
    dplyr::select(-p_bh_10) %>%
    dplyr::mutate(ratio = abs(StdError/Estimate)) %>%
    dplyr::select(ratio, starts_with("p_")) %>%
    tidyr::pivot_longer(., cols = starts_with("p_"),
                        names_to = "parameter", 
                        values_to = "value") %>%
    split(., f = .$parameter) %>%
    map(., ~.x %>%
          tidyr::drop_na(.) %>%
          droplevels %>%
          dplyr::select(ratio) %>%
          tibble::deframe(.) %>%
          quantile(., probs = seq(0, 1, 0.25), na.rm = TRUE)) %>%
    bind_rows(., .id = "p_method")
  # temp_res_cc_summary.df
  # #with the corncob-adjusted p-values, the 30 sig ASVs ratios of StdError:Estimate ranged between 0.007-0.324
  # #with q-value=0.05, the 21 sig ASVs with p-values passing that threshold had ratios of 0.007-0.277
  # #with q-value=0.01, the 5 sig ASVs with p-values passing that threshold had ratios of 0.008-0.198
  
  g <- print(
    ggplot(data = temp_df_to_plot)
    + geom_errorbar(aes(x = asv_id, ymin = Estimate-StdError, ymax = Estimate+StdError), color = "black", width = .3, position=position_dodge(.9)) 
    + geom_point(aes(x = asv_id, y = Estimate, group = interaction(asv_id, variable)), fill = "grey", alpha = 0.5,
                 size = 1, pch = 21)
    + geom_point(data = temp_df_to_plot %>%
                   tidyr::drop_na(p_bh_05) %>%
                   droplevels,
                 aes(x = asv_id, y = Estimate, fill = variable, group = interaction(asv_id, variable)),
                 size = 4, pch = 21) 
    + scale_x_discrete(limits = rev(levels(temp_df_to_plot$asv_id))) 
    # + scale_fill_manual(values = palette) 
    # + geom_line() 
    + coord_flip()
    + theme_bw()
    + labs(x = "Taxa", y = "Coefficient", title = "Yawzi")
    + geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")
    + facet_grid(.~variable, scales = "free")
  )
  
  temp_df2 <- temp_res_cc_df %>%
    dplyr::filter(condition == unique(temp_res_cc_list$condition)[2]) %>%
    dplyr::select(asv_id, Estimate, StdError, p_cc, p_bh_05, p_bh_01) %>%
    dplyr::arrange(Estimate) %>%
    dplyr::mutate(asv_id = factor(asv_id, levels = unique(.[["asv_id"]]))) %>%
    tidyr::pivot_longer(., cols = starts_with("p_"),
                        names_to = "parameter",
                        values_to = "value") %>%
    dplyr::mutate(Estimate = dplyr::case_when(!is.na(value) ~ Estimate,
                                              .default = NA),
                  StdError = dplyr::case_when(!is.na(value) ~ StdError,
                                              .default = NA)) %>%
    dplyr::mutate(parameter = factor(parameter, levels = c("p_cc", "p_bh_05", "p_bh_01"), labels = c("corncob's adjusted p-value", "ashr adjusted p-value (q =0.05)", "ashr adjusted p-value (q =0.01)"))) %>%
    droplevels 
  
  temp_g1 <- print(ggplot(data = temp_df %>%
                            dplyr::filter(grepl("p_cc|p_bh_0", parameter)) %>%
                            tidyr::drop_na(.) %>%
                            dplyr::mutate(parameter = factor(parameter, levels = c("p_cc", "p_bh_05", "p_bh_01"), labels = c("corncob's adjusted p-value", "ashr adjusted p-value (q =0.05)", "ashr adjusted p-value (q =0.01)"))) %>%
                            droplevels)
                   # + geom_boxplot(aes(x = parameter, y = value, fill = parameter), outliers = TRUE, outlier.shape = NA)
                   + geom_violin(aes(x = parameter, y = value, fill = parameter), draw_quantiles = c(0.05, 0.1))
                   + geom_hline(aes(yintercept = 0.05), color = "black")
                   # + geom_hline(aes(yintercept = 0.10), color = "maroon")
                   + scale_y_continuous(transform = "log10", name = "p-values")
                   + theme_bw()
                   + theme(axis.text.x = element_blank(),
                           legend.position = "right")
  )
  temp_g2 <- print(ggplot(data = temp_df2 %>%
                            tidyr::drop_na(.) %>%
                            droplevels)
                   + theme_bw()
                   + geom_histogram(aes(x = Estimate, fill = parameter), color = "black", bins = 10)
                   + facet_wrap(.~parameter)
                   + theme(legend.position = "none")
  )
  temp_g3 <- print(ggplot(data = temp_df2 %>%
                            tidyr::drop_na(.) %>%
                            droplevels)
                   + theme_bw()
                   + geom_violin(aes(x = parameter, y = Estimate, fill = parameter))
                   + geom_point(aes(x = parameter, y = Estimate, group = interaction(parameter, asv_id)), 
                                position = position_jitter(width = 0.1, height = 0),
                                fill = "grey", alpha = 0.5, size = 1, pch = 21)
                   + theme(axis.text.x = element_blank(),
                           legend.position = "none")
  )
  # temp_g2 <- print(ggplot(data = temp_df %>%
  #                           dplyr::filter(type == "q_values"))
  #                  # + geom_boxplot(aes(x = parameter, y = value, fill = parameter), outliers = TRUE, outlier.shape = NA)
  #                  + geom_violin(aes(x = parameter, y = value, fill = parameter), draw_quantiles = c(0.05, 0.1))
  #                  + geom_crossbar(data = temp_df_labels %>% dplyr::filter(grepl("q", parameter)), 
  #                                  aes(x = parameter, y = value, ymin = value, ymax = value), color = "black")
  #                  + geom_text(data = temp_df_labels %>% dplyr::filter(grepl("q", parameter)),
  #                              aes(x = parameter, y = value, label = quantile), nudge_x = 0.05)
  #                  + scale_y_continuous(transform = "log10", name = "Values", labels = scales::label_number_auto())
  #                  # + geom_hline(aes(yintercept = 0.05), color = "black")
  #                  # + geom_hline(aes(yintercept = 0.10), color = "maroon")
  # )
  
  temp_g <- temp_g2 / (temp_g3 + temp_g2) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Yawzi afternoon vs dawn SDA ASVs", tag_levels = "A")
  if(!file.exists(paste0(projectpath, "/", "modeled_p_values_from_corncob", ".png"))){
    ggsave(paste0(projectpath, "/", "modeled_p_values_from_corncob", ".png"),
           temp_g,
           width = 10, height = 10, units = "in")  
  }
  
}


# Summarize corncob results -----------------------------------------------

#use q-value = 0.05 on raw p-values to filter significant ASVS
q_value <- 0.05
# ds_objects <- c("ps_usvi_filtered_agglom-site-sampling_time", "ps_usvi_filtered-site-sampling_time")
ds_objects <- c("ps_usvi_filtered-site-sampling_time")


# dt_results_df <- `cc_dt_res_ps_usvi_filtered-site-sampling_time_list`[1:2] %>%
#   Filter(Negate(function(x) is.null(unlist(x))), .) %>% #remove the Null lists
#   map(., ~imap_dfr(., ~.x %>%
#                        purrr::pluck("p") %>%
#                        tibble::enframe(name = "asv_id", value = "p_val") %>%
#                        dplyr::mutate(p_fdr = p.adjust(p_val, "BH")) %>%
#                        dplyr::mutate(test = .y)) %>% split(., f = .$test))
# 
# padj_cutoff <- dt_results_df %>%
#   map_depth(., 2, ~.x %>%
#                      dplyr::select(p_val) %>%
#                      tibble::deframe(.) %>%
#                      na.omit(.) %>%
#                      ashr::qval.from.lfdr(.) %>%
#                      unlist %>%
#                      quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)) %>%
#   map(., ~imap(., ~.x %>% tibble::enframe(., name = NULL, value = "q_cutoff") %>%
#                  dplyr::mutate(test = .y))) %>%
#   map(., ~.x %>% bind_rows(., .id = NULL))
# 
# temp_list <- dt_results_df %>%
#   map(., ~.x %>%
#         bind_rows(., .id = NULL)) %>%
#   map2(., padj_cutoff, ~.x %>% dplyr::left_join(.y, by = join_by(test))) %>%
#   map(., ~.x %>%
#           dplyr::mutate(p_adj = dplyr::case_when(p_val <= q_cutoff ~ p_val,
#                                                  .default = NA)) %>%
#           dplyr::filter(!is.na(p_fdr)) %>%
#           droplevels) %>%
#   map(., ~.x %>%
#         split(., f = .$test) %>%
#         map(., ~.x %>%
#               droplevels))

for(i in ds_objects){
  dt_list <- paste0("cc_dt_res_", i, "_list")
  
  if(!exists(paste0("cc_dt_res_", i, "_df"), envir = .GlobalEnv, inherits = TRUE)){
    cli::cli_alert_info(paste0("Summarizing differentialTest results for this list: ", dt_list), wrap = TRUE)
    dt_list <- get0(dt_list, inherits = TRUE) #get the results list object
    
    #future_map takes a really long time on my MBP, which stinks for parallel potential.
    {
      # dt_results_df <- dt_list %>%
      #   Filter(Negate(function(x) is.null(unlist(x))), .) %>% #remove the Null lists
      #   furrr::future_map(., ~.x %>%
      #                       imap_dfr(., ~.x %>%
      #                                  purrr::pluck("p") %>%
      #                                  # tibble::as_tibble_col(column_name = "p_adj") %>% #this drops the ASV_id in later map
      #                                  tibble::as_tibble(rownames = "asv_id") %>%
      #                                  dplyr::rename(p_val = "value") %>%
      #                                  dplyr::mutate(p_fdr = p.adjust(p_val, "BH")) %>%
      #                                  dplyr::mutate(test = .y))) %>%
      #   furrr::future_map(., (~.x %>%
      #                           dplyr::mutate(across(c(asv_id, test), ~as.factor(.x))) %>%
      #                           # dplyr::group_by(asv_id) %>%
      #                           # dplyr::filter(p_fdr <= 0.05) %>%
      #                           # dplyr::filter(!is.na(p_fdr)) %>%
      #                         droplevels), .progress = TRUE) %>%
      # bind_rows(., .id = "contrast") %>%
      # tidyr::pivot_longer(., cols = !c(contrast, asv_id, test),
      #                     names_to = "metric",
      #                     values_to = "value") %>%
      # # dplyr::left_join(., usvi_prok_asvs.taxonomy %>%
      # #                    dplyr::select(asv_id, taxonomy),
      # #                  by = join_by(asv_id)) %>%
      # droplevels
      
      #post-hoc p adjustment:
      {
        # padj_cutoff <- dt_results_df %>%
        #   dplyr::select(p_val) %>%
        #   tibble::deframe(.) %>%
        #   na.omit(.) %>%
        # # padj_cutoff <- dt_results_df[["p_fdr"]] %>%
        # #   p.adjust(., method = "BH") %>% #multiple testing corrections, "BH" is an alias for "fdr" accordiong to p.adjust()
        #   ashr::qval.from.lfdr(.) %>%
        #   unlist %>%
        #   quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)
        # 
        # dt_results_df <- dt_results_df %>%
        #   dplyr::mutate(p_adj = dplyr::case_when(p_val <= padj_cutoff ~ p_val,
        #   # dplyr::mutate(p_adj = dplyr::case_when(p_fdr <= padj_cutoff ~ p_fdr,
        #                                          .default = NA)) %>%
        #   dplyr::filter(!is.na(p_fdr)) %>%
        #   droplevels  
      }
    }
    
    dt_results_df <- dt_list %>%
      Filter(Negate(function(x) is.null(unlist(x))), .) %>% #remove the Null lists
      map(., ~imap_dfr(., ~.x %>%
                         purrr::pluck("p") %>%
                         tibble::enframe(name = "asv_id", value = "p_val") %>%
                         dplyr::mutate(p_fdr = p.adjust(p_val, "BH")) %>%
                         dplyr::mutate(test = .y)) %>% split(., f = .$test))
    
    padj_cutoff <- dt_results_df %>%
      map_depth(., 2, ~.x %>%
                  dplyr::select(p_val) %>%
                  tibble::deframe(.) %>%
                  na.omit(.) %>%
                  ashr::qval.from.lfdr(.) %>%
                  unlist %>%
                  quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)) %>%
      map(., ~imap(., ~.x %>% tibble::enframe(., name = NULL, value = "q_cutoff") %>%
                     dplyr::mutate(test = .y))) %>%
      map(., ~.x %>% bind_rows(., .id = NULL))
    
    dt_results_df <- dt_results_df %>%
      map(., ~.x %>%
            bind_rows(., .id = NULL)) %>%
      map2(., padj_cutoff, ~.x %>% dplyr::left_join(.y, by = join_by(test))) %>%
      map(., ~.x %>%
            dplyr::mutate(p_adj = dplyr::case_when(p_val <= q_cutoff ~ p_val,
                                                   .default = NA)) %>%
            dplyr::filter(!is.na(p_fdr)) %>%
            droplevels)
    
    assign(paste0("cc_dt_res_", i, "_df"), dt_results_df, envir = .GlobalEnv, inherits = TRUE)
    readr::write_rds(dt_results_df, file = paste0(projectpath, "/", "cc_dt_res_", i, "_df", ".rds"))
    rm(dt_results_df)
    rm(dt_list)
  } else {
    cli::cli_alert_danger("Summary of differentialTest results list already exists for this object:", wrap = TRUE)
    cli::cli_bullets_raw(dt_list)
  }
}

# Evaluate the summarized results -----------------------------------------

cc_dt_usvi_summary.df <- `cc_dt_res_ps_usvi_filtered-site-sampling_time_df` %>%
  bind_rows(., .id = "contrast") %>% 
  dplyr::mutate(taxon_resolution = "asv") %>%
# cc_dt_usvi_summary.df <- list(`cc_dt_res_ps_usvi_filtered_agglom-site-sampling_time_df`, 
#                 `cc_dt_res_ps_usvi_filtered-site-sampling_time_df`) %>%
#   setNames(., c("agg_genus", "asv")) %>%
  # bind_rows(., .id = "taxon_resolution") %>%
  dplyr::mutate(hold = gsub("([[:alpha:]_]{1,})([[:print:]]{1,})$", "\\1", contrast) %>%
                  stringr::str_squish(.)) %>%
  dplyr::mutate(variable = gsub("([[:alpha:]_]{1,})([[:print:]]{1,})$", "\\2", contrast) %>%
                  stringr::str_remove_all(., "(\\(|\\))") %>%
                  stringr::str_squish(.)) %>%
  dplyr::mutate(hold = factor(hold, levels = c(names(sampling_time_lookup), names(site_lookup)))) %>%
  dplyr::arrange(hold, asv_id, taxon_resolution) %>%
  dplyr::relocate(hold, variable, .after = "contrast") %>%
  droplevels

if(!any(grepl(Sys.Date(), list.files(projectpath, pattern = "$cc_dt_usvi_summary.*.tsv")))){
# if(!file.exists(paste0(projectpath, "/", "cc_dt_usvi_summary.df", ".tsv"))){
  readr::write_delim(cc_dt_usvi_summary.df, paste0(projectpath, "/", "cc_dt_usvi_summary.df-", Sys.Date(), ".tsv"),
                     delim = "\t", col_names = TRUE)
}

#to check the numbers of ASVs identified as SDA in corncob:
#in the spreadsheet, where the columsn look like:
#contrast	hold	variable	asv_id	p_val	p_fdr	test	q_cutoff	p_adj	taxon_resolution
#dawn (Tektite - LB_seagrass)	dawn	Tektite - LB_seagrass	ASV_00644	0.010849634	0.045206808	cc_da	1.12E-14	NA	asv
#
#sort first by p_fdr =< 0.05 to get 810 observations
#then distinct the ASVs to get 207 unique ASVs



# Summarise SDA taxa ------------------------------------------------------


# #build lists of the SDA taxa in each comparison obtained through corncob:
# temp_df1 <- phyloseq::otu_table(ps_usvi_filtered_agglom) %>%
#   t() %>%
#   as.data.frame() %>%
#   apply(., 1, relabund) %>%
#   tibble::as_tibble(., rownames = "asv_id") %>%
#   tidyr::pivot_longer(., cols = !c(asv_id),
#                       names_to = "sample_id",
#                       values_to = "relabund") %>%
#   dplyr::mutate(across(c(asv_id, sample_id), ~factor(.x))) %>%
#   droplevels
# 
# cc_dt_sda_ps_usvi_filtered_agglom_list <- cc_dt_usvi_summary.df %>%
#   dplyr::filter(grepl("agg_genus", taxon_resolution)) %>%
#   split(., f = .$contrast) %>%
#   map(., ~.x %>%
#         dplyr::distinct(asv_id, hold, variable, .keep_all = FALSE) %>%
#         tidyr::separate_longer_delim(variable, delim = " - ") %>%
#         droplevels %>%
#         dplyr::left_join(., phyloseq::sample_data(ps_usvi_filtered_agglom) %>%
#                            tibble::as_tibble(.) %>%
#                            dplyr::select(sample_id, sampling_time, site),
#                          relationship = "many-to-many", multiple = "all",
#                          by = join_by("hold" == "sampling_time", "variable" == "site")) %>%
#         dplyr::left_join(., phyloseq::sample_data(ps_usvi_filtered_agglom) %>%
#                            tibble::as_tibble(.) %>%
#                            dplyr::select(sample_id, sampling_time, site),
#                          relationship = "many-to-many", multiple = "all",
#                          by = join_by("hold" == "site", "variable" == "sampling_time")) %>%
#         dplyr::mutate(sample_id = across(contains("sample_id")) %>% purrr::reduce(coalesce)) %>%
#         dplyr::select(-ends_with(c(".x", ".y"))) %>%
#         droplevels) %>%
#   map(., ~.x %>%
#         dplyr::mutate(asv_id = factor(asv_id),
#                       sample_id = factor(sample_id)) %>%
#         dplyr::left_join(., temp_df1,
#                          relationship = "many-to-many", multiple = "all",
#                          by = join_by(sample_id, asv_id)) %>%
#         droplevels)
# 
# temp_df2 <- phyloseq::otu_table(ps_usvi_filtered) %>%
#   t() %>%
#   as.data.frame() %>%
#   apply(., 1, relabund) %>%
#   tibble::as_tibble(., rownames = "asv_id") %>%
#   tidyr::pivot_longer(., cols = !c(asv_id),
#                       names_to = "sample_id",
#                       values_to = "relabund") %>%
#   dplyr::mutate(across(c(asv_id, sample_id), ~factor(.x))) %>%
#   droplevels
# 
# cc_dt_sda_ps_usvi_filtered_list <- cc_dt_usvi_summary.df %>%
#   dplyr::filter(grepl("asv", taxon_resolution)) %>%
#   split(., f = .$contrast) %>%
#   map(., ~.x %>%
#         dplyr::distinct(asv_id, hold, variable, .keep_all = FALSE) %>%
#         tidyr::separate_longer_delim(variable, delim = " - ") %>%
#         droplevels %>%
#         dplyr::left_join(., phyloseq::sample_data(ps_usvi_filtered) %>%
#                            tibble::as_tibble(., rownames = "sample_id") %>%
#                            dplyr::select(sample_id, sampling_time, site),
#                          relationship = "many-to-many", multiple = "all",
#                          by = join_by("hold" == "sampling_time", "variable" == "site")) %>%
#         dplyr::left_join(., phyloseq::sample_data(ps_usvi_filtered) %>%
#                            tibble::as_tibble(., rownames = "sample_id") %>%
#                            dplyr::select(sample_id, sampling_time, site),
#                          relationship = "many-to-many", multiple = "all",
#                          by = join_by("hold" == "site", "variable" == "sampling_time")) %>%
#         dplyr::mutate(sample_id = across(contains("sample_id")) %>% purrr::reduce(coalesce)) %>%
#         dplyr::select(-ends_with(c(".x", ".y"))) %>%
#         droplevels) %>%
#   map(., ~.x %>%
#         dplyr::mutate(asv_id = factor(asv_id),
#                       sample_id = factor(sample_id)) %>%
#         dplyr::left_join(., temp_df2,
#                          relationship = "many-to-many", multiple = "all",
#                          by = join_by(sample_id, asv_id)) %>%
#         droplevels)
# 
# rm(temp_df1)
# rm(temp_df2)
# 
# if(!any(grepl("cc_dt_sda", list.files(projectpath, pattern = ".*ps_usvi_filtered-.*.RData")))){
#   save(cc_dt_sda_ps_usvi_filtered_list, cc_dt_sda_ps_usvi_filtered_agglom_list, file = paste0(projectpath, "/", "cc_dt_sda_ps_usvi_filtered-", Sys.Date(), ".RData"))
#   # readr::write_rds(cc_dt_sda_ps_usvi_filtered_list, paste0(projectpath, "/", "cc_dt_sda_ps_usvi_filtered_list", ".rds", compress = "gz"))
#   # readr::write_rds(cc_dt_sda_ps_usvi_filtered_agglom_list, paste0(projectpath, "/", "cc_dt_sda_ps_usvi_filtered_agglom_list", ".rds", compress = "gz"))
# }
# 
# 
# #summary of results:
# cc_dt_usvi_summary_table <- cc_dt_usvi_summary.df %>%
#   dplyr::group_by(contrast, taxon_resolution, test) %>%
#   dplyr::summarise(num_sda_taxa_cc = length(asv_id)) %>%
#   # tidyr::complete(test, fill = list(num_sda_taxa_cc = 0)) %>%
#   dplyr::left_join(., cc_dt_usvi_summary.df %>%
#                      tidyr::drop_na(p_adj) %>%
#                      droplevels %>%
#                      dplyr::group_by(contrast, taxon_resolution, test) %>%
#                      dplyr::summarise(num_sda_taxa_padj = length(asv_id)) %>%
#                      droplevels,
#                    by = join_by(contrast, taxon_resolution, test)) %>%
#   tidyr::complete(test, fill = list(num_sda_taxa_cc = 0, num_sda_taxa_padj = 0))
# 
# #using corncob's multiple test correction criteria:
# #site-specific temporal (dawn vs peak_photo) differences in abundances (SDA) and variability (SDV):
# #LB_seagrass had 7 genera that were SDA between dawn and peak_photo, and 0 were SDV; 5 ASVs were SDA and 1 ASV was SDV
# #Tektite had 3 genera that were SDA and SDV. 9 ASVs were SDV but 0 were SDA
# #Yawzi had 24 genera that were SDA and 0 were SDV. 43 ASVs were SDA and 13 ASVs were SDV
# 
# #dawn- specific (Yawzi vs seagrass, Tektite vs seagrass, Tektite vs Yawzi)
# #between Tektite and LB, 124 genera were SDA, 6 were SDV. 149 ASVs were SDA, 11 were SDV
# #between Yawzi and LB, 138 genera were SDA, 32 were SDV. 144 ASVs were SDA, 83 were SDV
# #between Yawzi and Tektite, 23 genera were SDA, 15 were SDV. 21 ASVs were SDA, 57 were SDV
# 
# #peak_photo specific (Yawzi vs seagrass, Tektite vs seagrass, Tektite vs Yawzi)
# #between Tektite and LB, 95 genera were SDA, 13 were SDV. 92 ASVs were SDA, 38 were SDV
# #between Yawzi and LB, 103 genera were SDA, 43 were SDV. 104 ASVs were SDA, 112 wereSDV
# #between Yawzi and Tektite, 4 genera were SDA, 6 were SDV. 2 ASVs were SDA, 65 were SDV
# 
# 
# # Plot the relabunds of genera found SDA/SDV through cc -------------------
# 
# #corncob says 83 genera were SDA between 	peak_photo (Tektite - LB_seagrass)
# #but post-test correction narrows it to 5 genera
# 
# temp_df <- cc_dt_sda_ps_usvi_filtered_agglom_list[["peak_photo (Yawzi - Tektite)"]] %>%
#   dplyr::left_join(., cc_dt_usvi_summary.df %>%
#                      dplyr::filter(grepl("agg_genus", taxon_resolution)) %>%
#                      dplyr::filter(contrast == "peak_photo (Yawzi - Tektite)") %>%
#                      dplyr::filter(test == "cc_da") %>%
#                      dplyr::select(asv_id, hold, test, p_fdr, p_adj, taxonomy) %>%
#                      droplevels, by = join_by(hold, asv_id), relationship = "many-to-many", multiple = "all") %>%
#   dplyr::rename(sampling_time = "hold", site = "variable") %>%
#   droplevels
# 
# g1_sigboth <- print(
#   ggplot(data = temp_df %>%
#            tidyr::drop_na(p_adj) %>%
#            droplevels)
#   + theme_bw() 
#   + geom_boxplot(aes(y = relabund, x = site, group = interaction(site, sampling_time)), 
#                  alpha = 0.6, show.legend = FALSE, color = "black",
#                  position = position_dodge(0.8))
#   + geom_point(aes(y = relabund, group = sampling_time, x = site, fill = sampling_time, shape = site),
#                alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
#                position = position_jitter(width = 0.2, seed = 48105))
#   + scale_shape_manual(name = "Sampling site and time",
#                        values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#   + scale_fill_manual(name = "Sampling time",
#                       values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_x_discrete(labels = site_lookup, name = "Site")
#   + scale_y_continuous(name = "Relative abundance (%)")
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                override.aes = list(color = "black", shape = 21, size = 3)),
#            color = "none")
#   + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
#   # # + theme(axis.text.x = element_blank(),
#   #         strip.text.y = element_text(size = rel(0.7)),
#   #         axis.title.y = element_blank())
#   # + facet_grid(space = "fixed",
#   # + facet_wrap(
#   #              asv_id ~ ., 
#   #              drop = TRUE,
#   #              # labeller = global_labeller,
#   #              scales = "free_y")
# )
# 
# g1_sigcc <- print(
#   ggplot(data = temp_df %>%
#            dplyr::filter(is.na(p_adj)) %>%
#            droplevels)
#   + theme_bw() 
#   + geom_boxplot(aes(y = relabund, x = site, group = interaction(site, sampling_time)), 
#                  alpha = 0.6, show.legend = FALSE, color = "black",
#                  position = position_dodge(0.8))
#   + geom_point(aes(y = relabund, group = sampling_time, x = site, fill = sampling_time, shape = site),
#                alpha = 0.8, color = "black", size = 3, show.legend = TRUE, 
#                position = position_jitter(width = 0.2, seed = 48105))
#   + scale_shape_manual(name = "Sampling site and time",
#                        values = c(22, 21, 23), labels = c(site_lookup, sampling_time_lookup), breaks = c(names(site_lookup), sampling_time_lookup))
#   + scale_fill_manual(name = "Sampling time",
#                       values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_x_discrete(labels = site_lookup, name = "Site")
#   + scale_y_continuous(name = "Relative abundance (%)")
#   + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling site and time",  direction = "vertical", 
#                                override.aes = list(color = "black", shape = 21, size = 3)),
#            color = "none")
#   + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
#   + facet_wrap(
#     asv_id ~ ., 
#     drop = TRUE,
#     scales = "free_y")
# )
