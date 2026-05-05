# 11_asv_transformations.R

# do alternative transformations ont he ASV data affect NMDS, pemranova, etc?


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
  
  
}

# Load additional packages ------------------------------------------------

library(tidyverse)
library(phyloseq)
library(ggvegan)
library(ggtext)
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)
library(UpSetR)


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


if(ncol(usvi_prok_decontam_idx) == 1){
  usvi_prok_decontam_idx <- usvi_prok_decontam_idx %>%
    tidyr::separate_wider_delim(1, names = c("asv_id", "keep"), delim = " ", cols_remove = TRUE, too_few = "align_start")
}


usvi_prok_asvs.taxa <- usvi_prok_asvs.taxa %>%
  dplyr::left_join(., usvi_prok_decontam_idx, by = join_by(asv_id)) %>%
  dplyr::filter(keep == TRUE) %>%
  dplyr::select(-keep) %>%
  droplevels
if(any(grepl("Proteobacteria", usvi_prok_asvs.taxa[["Phylum"]], ignore.case = TRUE))){
  
  #update taxonomy:
  update_taxonomy.df <- readr::read_delim("~/projects/silva/update_taxonomy_sorted.tsv", comment = "#",
                                          delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
    # dplyr::select(-c(level)) %>%
    tidyr::separate_wider_delim(original, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
    droplevels
  update_taxonomy.lookup <- update_taxonomy.df %>%
    dplyr::select(original, new) %>%
    # dplyr::select(new, original) %>%
    tibble::deframe(.)
  
  
  temp_df <- usvi_prok_asvs.taxa %>%
    dplyr::select(-sequence) %>%
    # dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Proteobacteria",
    #                                         grepl("Alphaproteobacteria", Class) ~ "Proteobacteria",
    # .default = Phylum)) %>%
    tidyr::unite("taxonomy", c(all_of(keep_tax), "Species"), sep = ";", remove = FALSE) %>%
    # dplyr::rename(original = "taxonomy") %>%
    dplyr::mutate(Phylum = paste(Domain, Phylum, sep = ";")) %>%
    dplyr::mutate(Class = paste(Phylum, Class, sep = ";")) %>%
    dplyr::mutate(Order = paste(Class, Order, sep = ";")) %>%
    dplyr::mutate(Family = paste(Order, Family, sep = ";")) %>%
    dplyr::mutate(Genus = paste(Family, Genus, sep = ";")) %>%
    tidyr::pivot_longer(., cols = c(all_of(keep_tax), "Species"),
                        names_to = "level",
                        values_to = "value") %>%
    dplyr::mutate(level = stringr::str_to_lower(level)) %>%
    dplyr::mutate(value = dplyr::case_when(grepl("Synechococcus", value) ~ stringr::str_remove_all(value, " .*$"),
                                           .default = value)) %>%
    dplyr::rename(original = "value") %>%
    droplevels
  
  temp_df2 <- temp_df %>%
    dplyr::left_join(., update_taxonomy.df %>%
                       dplyr::select(level, original, new) %>%
                       droplevels) %>%
    tidyr::drop_na(new) %>%
    # dplyr::mutate(`value` = `new`) %>%
    dplyr::select(taxonomy, new) %>%
    tidyr::separate_wider_delim(new, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
    tidyr::pivot_longer(., cols = all_of(keep_tax),
                        names_to = "level",
                        values_to = "value") %>%
    dplyr::mutate(level = stringr::str_to_lower(level)) %>%
    tidyr::drop_na(value) %>%
    droplevels 
  
  temp_df3 <- temp_df %>%
    dplyr::distinct(taxonomy, asv_id, level) %>%
    droplevels %>%
    dplyr::left_join(., temp_df2 %>%
                       dplyr::distinct(taxonomy, level, value) %>%
                       droplevels,
                     relationship = "many-to-many", multiple = "all") %>%
    # dplyr::select(-new) %>%
    dplyr::mutate(`level` = stringr::str_to_sentence(`level`)) %>%
    tidyr::pivot_wider(., id_cols = c("asv_id", "taxonomy"),
                       names_from = "level",
                       values_from = "value") %>%
    droplevels %>%
    dplyr::select(asv_id, taxonomy, all_of(keep_tax), Species) %>%
    # # # dplyr::select(-level) %>%
    dplyr::right_join(., (usvi_prok_asvs.taxa %>%
                            dplyr::select(-sequence) %>%
                            # tidyr::separate_wider_delim(taxonomy, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
                            # dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Proteobacteria",
                            #                                         grepl("Alphaproteobacteria", Class) ~ "Proteobacteria",
                            #                                         .default = Phylum)) %>%
                            tidyr::unite("taxonomy", c(all_of(keep_tax), "Species"), sep = ";", remove = FALSE) %>%
                            # dplyr::distinct(taxonomy,asv_id, .keep_all = FALSE) %>%
                            droplevels),
                      by = join_by("taxonomy", "asv_id"), relationship = "many-to-many", multiple = "all") %>%
    dplyr::mutate(Domain = across(starts_with("Domain")) %>% purrr::reduce(coalesce)) %>%
    dplyr::mutate(Phylum = across(starts_with("Phylum")) %>% purrr::reduce(coalesce)) %>%
    dplyr::mutate(Class = across(starts_with("Class")) %>% purrr::reduce(coalesce)) %>%
    dplyr::mutate(Order = across(starts_with("Order")) %>% purrr::reduce(coalesce)) %>%
    dplyr::mutate(Family = across(starts_with("Family")) %>% purrr::reduce(coalesce)) %>%
    dplyr::mutate(Genus = across(starts_with("Genus")) %>% purrr::reduce(coalesce)) %>%
    dplyr::mutate(Species = across(starts_with("Species")) %>% purrr::reduce(coalesce)) %>%
    dplyr::select(-ends_with(c(".x", ".y"))) %>%
    tidyr::unite("new_taxonomy", c(all_of(keep_tax)), sep = ";", remove = FALSE) %>%
    dplyr::mutate(new_taxonomy = dplyr::case_when(grepl(";NA", new_taxonomy) ~ stringr::str_remove_all(new_taxonomy, ";NA"),
                                                  .default = new_taxonomy)) %>%
    droplevels   
  
  usvi_prok_asvs_renamed.taxa <- usvi_prok_asvs.taxa %>%
    dplyr::select(asv_id, sequence) %>%
    dplyr::left_join(., temp_df3) %>%
    dplyr::select(-taxonomy) %>%
    dplyr::rename(taxonomy = "new_taxonomy") %>%
    dplyr::select(all_of(colnames(usvi_prok_asvs.taxa))) %>%
    droplevels
  
  usvi_prok_filled.taxa.df <- usvi_prok_asvs_renamed.taxa %>%
    dplyr::mutate(Phylum = coalesce(Phylum, Domain)) %>%
    dplyr::mutate(Class = coalesce(Class, Phylum)) %>%
    dplyr::mutate(Order = coalesce(Order, Class)) %>%
    dplyr::mutate(Family = coalesce(Family, Order)) %>%
    dplyr::mutate(Genus = coalesce(Genus, Family)) %>%
    dplyr::mutate(Species = coalesce(Species, Genus)) %>%
    dplyr::mutate(across(everything(), ~factor(.x))) %>%
    dplyr::relocate(asv_id) %>%
    droplevels
  
  usvi_prok_asvs.taxa <- usvi_prok_asvs_renamed.taxa
  
} else {
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
    # dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
    #                                         grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
    #                                         .default = Phylum)) %>%
    
    droplevels
  
}


if(file.exists(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"))){
  metabolomics_sample_metadata <- readr::read_delim(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
} else {
  if(!exists("ps_usvi", envir = .GlobalEnv)){
    if(!file.exists(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))){
      cli::cli_alert_warning("Please process the USVI data through Phyloseq.")
    } else {
      ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))
      
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
  }
}


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





# Explore whether CLR transformation matters for ASVs ---------------------

set.seed(48105)

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
  dplyr::mutate(site = factor(site, levels = names(site_lookup))) %>%
  droplevels %>%
  dplyr::mutate(rownames = sample_id) %>% tibble::column_to_rownames(var = "rownames") %>%
  dplyr::mutate(grouping = interaction(site, sampling_time)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_id, replicate)) %>%
  # tidyr::drop_na(.) %>% #removed sample Metab_306 due to lack of FCM measurements
  droplevels

usvi_asv.tbl <- usvi_prok_asvs.df %>%
  dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  dplyr::filter(!grepl("Metab_280", sample_ID)) %>% #drop the sample corresponding to CINAR_BC_73: Metab_280
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  # usvi_asv.tbl <- ps_usvi %>%
  #   phyloseq::subset_samples(., sample_type == "seawater") %>%
  #   phyloseq::otu_table(.) %>%
  #   as.data.frame %>%
  # apply(., 2, relabund) %>% 
  # as.data.frame(.) %>%
  # dplyr::slice(which(rowSums(.) > 0)) %>%
  # tibble::rownames_to_column(var = "asv_id") %>%
  # tidyr::pivot_longer(., cols = -c("asv_id"),
  #                     names_to = "sample",
  #                     values_to = "abundance") %>%
  # dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
  #                                 log2(abundance+1), #log transform abundance (with +1 pseudocount)
  #                                 0)) %>%
  droplevels

transformations <- c("counts", "relabund", "clr", "relabund_clr", "rclr", "relabund_rclr")
transformed_usvi_asv_counts.tbl <- usvi_asv.tbl %>%
  t()
transformed_usvi_asv_relabund.tbl <- usvi_asv.tbl %>%
  apply(., 2, relabund)  %>%
  t()
transformed_usvi_asv_clr.tbl <- usvi_asv.tbl %>%
  vegan::decostand(., method = "clr", 2, na.rm = TRUE, pseudocount = 1) %>%
  t()
transformed_usvi_asv_relabund_clr.tbl <- usvi_asv.tbl %>%
  apply(., 2, relabund)  %>%
  vegan::decostand(., method = "clr", 2, na.rm = TRUE, pseudocount = 0.1) %>%
  t()
transformed_usvi_asv_rclr.tbl <- usvi_asv.tbl %>%
  vegan::decostand(., method = "rclr", 2, na.rm = TRUE, impute = TRUE) %>%
  t()
transformed_usvi_asv_relabund_rclr.tbl <- usvi_asv.tbl %>%
  apply(., 2, relabund)  %>%
  vegan::decostand(., method = "rclr", 2, na.rm = TRUE, impute = TRUE) %>%
  t()

transformations_usvi_asv.list <- list(transformed_usvi_asv_counts.tbl,
                                      transformed_usvi_asv_relabund.tbl,
                                      transformed_usvi_asv_clr.tbl,
                                      transformed_usvi_asv_relabund_clr.tbl,
                                      transformed_usvi_asv_rclr.tbl,
                                      transformed_usvi_asv_relabund_rclr.tbl) %>%
  setNames(., c("counts", "relabund", "clr", "relabund_clr", "rclr", "relabund_rclr"))


nmds_keep <- data.frame(transformation = c("relabund", "relabund", "relabund", "relabund_clr", "clr", "clr", "relabund_clr", "relabund_clr",
                                           "relabund_clr", "clr",
                                           "relabund_rclr", "rclr", "relabund_rclr", "rclr"),
                        method = c("horn", "bray", "euclidean", "chord", "chord", "bray", "bray", "horn",
                                   "euclidean", "euclidean",
                                   "horn", "horn", "euclidean", "euclidean"))
temp_nmds <- transformations_usvi_asv.list[["clr"]] %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "euclidean",trymax = 100, weakties = FALSE, autotransform = FALSE, 
                 tidy = TRUE,
                 wascores = TRUE)
#...
# Run 99 stress 0.1196607 
# Run 100 stress 0.1209625 
# *** Best solution was not repeated -- monoMDS stopping criteria:
#   16: no. of iterations >= maxit
# 83: stress ratio > sratmax
# 1: scale factor of the gradient < sfgrmin

transformed_usvi_asv_nmds_keep.list <- NULL
transformed_usvi_asv_nmds.list <- NULL

# for(i in seq(1)){
for(i in seq(nrow(nmds_keep))){
  meth <- nmds_keep[i, "method"]
  transformation <- nmds_keep[i, "transformation"]
  namevar <- paste0(transformation, "_", meth)
  
  temp_transformed.tbl <- transformations_usvi_asv.list %>%
    purrr::pluck( {transformation} ) %>%
    as.matrix(.)
  
  temp_nmds.list <- temp_transformed.tbl %>%
    # vegan::vegdist(., method = {meth}, binary = FALSE, na.rm = TRUE) %>% 
    vegan::metaMDS(., distance = {meth},trymax = 100, weakties = FALSE, autotransform = FALSE, 
                   # tidy = FALSE,
                   tidy = TRUE,
                   wascores = TRUE)
  
  temp_nmds.df <- c("stress", "tries", "points", "species") %>%
    map(., ~purrr::pluck(temp_nmds.list, .x)) %>%
    setNames(., c("stress", "tries", "points", "species"))
  temp_nmds.df <- temp_nmds.df %>%
    purrr::modify_at("points", ~.x %>%
                       tibble::as_tibble(rownames = "label") %>%
                       dplyr::rename(NMDS1 = "MDS1", NMDS2 = "MDS2") %>% 
                       dplyr::mutate(type = "sample") %>%
                       dplyr::left_join(., (metadata %>%
                                              dplyr::filter(grepl("seawater", sample_type)) %>%
                                              dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, replicate) %>%
                                              droplevels),
                                        by = c("label" = "sample_id")) %>%
                       dplyr::mutate(site = factor(site, levels = names(site_lookup)),
                                     sampling_time = factor(sampling_time, levels = names(sampling_time_lookup))))
  if(length(temp_nmds.df[["species"]]) > 1){
    temp_nmds.df <- temp_nmds.df %>%
      purrr::modify_at("species", ~.x %>%
                         tibble::as_tibble(rownames = "label") %>%
                         dplyr::rename(NMDS1 = "MDS1", NMDS2 = "MDS2") %>%
                         dplyr::mutate(type = "asv_id") %>%
                         dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"),
                                          by = c("label" = "asv_id")) %>%
                         droplevels)
  } else {
    temp_nmds.df <- temp_nmds.df %>%
      purrr::discard_at("species")
    # map(., ~.x %>% Filter(Negate(function(x) length(unlist(x)) == 0), .)) %>%
    # # Filter(Negate(function(x) is.null(unlist(x))), .)
    # Filter(Negate(function(x) is.na(unlist(x))), .)
  }
  
  
  if(exists("transformed_usvi_asv_nmds_keep.list", envir = .GlobalEnv)){
    transformed_usvi_asv_nmds_keep.list <- append(transformed_usvi_asv_nmds_keep.list, list(temp_nmds.df)) %>%
      setNames(., c(names(transformed_usvi_asv_nmds_keep.list), namevar))
    transformed_usvi_asv_nmds.list <- append(transformed_usvi_asv_nmds.list, list(temp_nmds.list)) %>%
      setNames(., c(names(transformed_usvi_asv_nmds.list), namevar))
  } else {
    transformed_usvi_asv_nmds_keep.list <- list(temp_nmds.df) %>%
      setNames(., namevar)
    transformed_usvi_asv_nmds.list <- list(temp_nmds.list) %>%
      setNames(., namevar)
  }
  rm(temp_nmds.df)
  rm(temp_nmds.list)
  rm(temp_transformed.tbl)
  rm(namevar)
}

View(transformed_usvi_asv_nmds_keep.list)
View(transformed_usvi_asv_nmds.list)
if(!file.exists(paste0(projectpath, "/", "transformed_usvi_asv_nmds_keep.list", ".rds"))){
  readr::write_rds(transformed_usvi_asv_nmds_keep.list, 
                   paste0(projectpath, "/", "transformed_usvi_asv_nmds_keep.list", ".rds"), 
                   compress = "gz")
  readr::write_rds(transformed_usvi_asv_nmds.list, 
                   paste0(projectpath, "/", "transformed_usvi_asv_nmds.list", ".rds"), 
                   compress = "gz")
}

# Warning messages:
#   1: In distfun(comm, method = distance, na.rm = TRUE, ...) :
#   results may be meaningless because data have negative entries
# in method “bray”
# 2: In metaMDSdist(comm, distance = distance, autotransform = autotransform,  :
#                     some dissimilarities exceed expected maximum 1
#                   3: In distfun(comm, method = distance, na.rm = TRUE, ...) :
#                     results may be meaningless because data have negative entries
#                   in method “bray”
#                   4: In metaMDSdist(comm, distance = distance, autotransform = autotransform,  :
#                                       some dissimilarities exceed expected maximum 1
#                                     5: In distfun(comm, method = distance, na.rm = TRUE, ...) :
#                                       results may be meaningless because data have negative entries
#                                     in method “horn”
which(!map(transformed_usvi_asv_nmds_keep.list, ~.x[["species"]] %>% length(.)) > 1)
# relabund_clr_chord               clr_chord                clr_bray       relabund_clr_bray 
# 4                       5                       6                       7 
# relabund_clr_horn  relabund_clr_euclidean           clr_euclidean      relabund_rclr_horn 
# 8                       9                      10                      11 
# rclr_horn relabund_rclr_euclidean          rclr_euclidean 
# 12                      13                      14 
#so any distance matrix with CLR or RCLR transformation applied to it

metrics <- c("horn",
             # "chisq",  #chi-squared does not work for negative values, such as those obtained after CLR transformation
             # "aitchison", #this accepts pseudocount parameter from decostand
             "chord",
             "bray",
             "euclidean",
             "hellinger"
)


usvi_asv_nmds_results.df <- transformed_usvi_asv_nmds_keep.list %>%
  map(., purrr::pluck("stress")) %>%
  list(., map(transformed_usvi_asv_nmds_keep.list, purrr::pluck("tries"))) %>%
  setNames(., c("stress", "tries")) %>%
  bind_rows(., .id = "type") %>%
  data.table::transpose(keep.names= "listlabel", make.names = "type") %>%
  dplyr::mutate(transformation = gsub("^(\\w+)(_)(\\w+)$", "\\1", listlabel)) %>%
  dplyr::mutate(method = stringr::str_extract_all(listlabel, paste0(metrics, collapse = "|"))) %>%
  dplyr::mutate(method = stringr::str_trim(method, "both")) %>%
  dplyr::select(transformation, method, stress, tries, listlabel) %>%
  droplevels


{
  #   for(meth in metrics){
  #     #counts: can do all metrics
  #     temp_nmds.list <- transformed_usvi_asv_counts.tbl %>%
  #       vegan::vegdist(., method = {meth}, binary = FALSE, na.rm = TRUE, pseudocount = 0.1) %>% 
  #       vegan::metaMDS(., distance = {meth},trymax = 100, autotransform = FALSE, tidy = TRUE)
  #     if(exists("usvi_asv_nmds_results.df", envir = .GlobalEnv)){
  #       temp_df <- c("stress", "tries") %>%
  #         map(., ~purrr::pluck(temp_nmds.list, .x)) %>%
  #         setNames(., c("stress", "tries")) %>%
  #         append(c("transformation" = "counts", "method" = {meth} ), .)  %>%
  #         dplyr::bind_rows(.)
  #       usvi_asv_nmds_results.df <- bind_rows(usvi_asv_nmds_results.df, temp_df)
  #       rm(temp_df)
  #     } else {
  #       usvi_asv_nmds_results.df <- temp_df
  #       rm(temp_df)
  #     }
  #     rm(temp_nmds.list)
  #     
  #     #relabund: can do all metrics
  #     temp_nmds.list <- transformed_usvi_asv_relabund.tbl %>%
  #       vegan::vegdist(., method = {meth}, binary = FALSE, na.rm = TRUE) %>% 
  #       vegan::metaMDS(., distance = {meth},trymax = 100, autotransform = FALSE, tidy = TRUE)
  #     if(exists("usvi_asv_nmds_results.df", envir = .GlobalEnv)){
  #       temp_df <- c("stress", "tries") %>%
  #         map(., ~purrr::pluck(temp_nmds.list, .x)) %>%
  #         setNames(., c("stress", "tries")) %>%
  #         append(c("transformation" = "relabund", "method" = {meth} ), .)  %>%
  #         dplyr::bind_rows(.)
  #       usvi_asv_nmds_results.df <- bind_rows(usvi_asv_nmds_results.df, temp_df)
  #       rm(temp_df)
  #     } else {
  #       usvi_asv_nmds_results.df <- temp_df
  #       rm(temp_df)
  #     }
  #     rm(temp_nmds.list)
  #     
  #     #clr: can't do Hellinger or Bray
  #     if(!(meth %in% c("bray", "hellinger"))){
  #       temp_nmds.list <- transformed_usvi_asv_clr.tbl %>%
  #         vegan::vegdist(., method = {meth}, binary = FALSE, na.rm = TRUE) %>% 
  #         vegan::metaMDS(., distance = {meth},trymax = 100, autotransform = FALSE, tidy = TRUE)
  #       if(exists("usvi_asv_nmds_results.df", envir = .GlobalEnv)){
  #         temp_df <- c("stress", "tries") %>%
  #           map(., ~purrr::pluck(temp_nmds.list, .x)) %>%
  #           setNames(., c("stress", "tries")) %>%
  #           append(c("transformation" = "clr", "method" = {meth} ), .)  %>%
  #           dplyr::bind_rows(.)
  #         usvi_asv_nmds_results.df <- bind_rows(usvi_asv_nmds_results.df, temp_df)
  #         rm(temp_df)
  #       } else {
  #         usvi_asv_nmds_results.df <- temp_df
  #         rm(temp_df)
  #       }
  #       rm(temp_nmds.list)
  #       
  #       
  #       
  #       #relabund_clr: can't do Hellinger
  #       temp_nmds.list <- transformed_usvi_asv_relabund_clr.tbl %>%
  #         vegan::vegdist(., method = {meth}, binary = FALSE, na.rm = TRUE) %>%
  #         vegan::metaMDS(., distance = {meth},trymax = 100, autotransform = FALSE, tidy = TRUE)
  #       if(exists("usvi_asv_nmds_results.df", envir = .GlobalEnv)){
  #         temp_df <- c("stress", "tries") %>%
  #           map(., ~purrr::pluck(temp_nmds.list, .x)) %>%
  #           setNames(., c("stress", "tries")) %>%
  #           append(c("transformation" = "relabund_clr", "method" = {meth} ), .)  %>%
  #           dplyr::bind_rows(.)
  #         usvi_asv_nmds_results.df <- bind_rows(usvi_asv_nmds_results.df, temp_df)
  #         rm(temp_df)
  #       } else {
  #         usvi_asv_nmds_results.df <- temp_df
  #         rm(temp_df)
  #       }
  #       rm(temp_nmds.list)
  #       # } else {
  #       #   
  #     }
  #   }
  #   readr::write_delim(usvi_asv_nmds_results.df, 
  #                      paste0(projectpath, "/", "usvi_asv_nmds_results.df", "-", Sys.Date(), ".tsv"),
  #                      delim = "\t")
}

usvi_asv_nmds_results.df %>%
  dplyr::arrange(stress) %>%
  dplyr::select(transformation, method, stress, tries) %>%
  # dplyr::slice_min(stress, n =1, by = "method") %>%
  droplevels
# transformation    method     stress tries
# 1            rclr euclidean 0.05522947    20
# 2   relabund_rclr euclidean 0.05522947    20
# 3        relabund      bray 0.07511278    20
# 4    relabund_clr     chord 0.07520584    35
# 5        relabund      horn 0.07573668    20
# 6    relabund_clr      horn 0.07579326    60
# 7    relabund_clr euclidean 0.07630130   100
# 8        relabund euclidean 0.07790935    20
# 9    relabund_clr      bray 0.08699991    42
# 10            clr     chord 0.10154309   100
# 11            clr euclidean 0.11085590   100
# 12           rclr      horn 0.16253916    20
# 13  relabund_rclr      horn 0.17425883    20
# 14            clr      bray 0.35380128   100

#originally in submitted manuscript: MH on relabund (stress = 0.0757)
#best using CLR: Chord distance on relabund+CLR (stress = 0.07520584)
#Just CLR + any distance, best stress is obtained with Chord distance (stress = 0.10154309)
#rclr + euclidean (stress = 0.05522947) is better than any other RCLR + distance (next best stress: rclr+MH, stress = 0.16253916)

to_plot <- usvi_asv_nmds_results.df %>%
  dplyr::arrange(stress) %>%
  dplyr::mutate(transformation = factor(transformation, levels = transformations)) %>%
  droplevels %>%
  # dplyr::distinct(stress, .keep_all = TRUE) %>%
  dplyr::filter(transformation %in% c("relabund", "clr", "rclr")) %>%
  # dplyr::filter(stress < 0.2 ) %>%
  dplyr::arrange(transformation) %>%
  dplyr::filter((stress < 0.2) | (method %in% c("bray", "euclidean", "horn"))) %>%
  # droplevels
  droplevels %>% dplyr::select(listlabel) %>% tibble::deframe(.)
# transformation    method     stress tries              listlabel
# 1            clr euclidean 0.11085590   100          clr_euclidean
# 2           rclr euclidean 0.05522947    20         rclr_euclidean
# 3           rclr      horn 0.16253916    20              rclr_horn
# 4       relabund      bray 0.07511278    20          relabund_bray
# 5       relabund      horn 0.07573668    20          relabund_horn
# 6       relabund euclidean 0.07790935    20     relabund_euclidean
# 7   relabund_clr      horn 0.07579326    60      relabund_clr_horn
# 8   relabund_clr euclidean 0.07630130   100 relabund_clr_euclidean
# 9   relabund_clr      bray 0.08699991    42      relabund_clr_bray

rm(list = apropos("g_transformed_usvi_asv_nmds_*", mode = "any"))
#plot the results:
for(i in to_plot){
  # for(i in seq(length(transformed_usvi_asv_nmds_keep.list))){
  # namevar <- names(transformed_usvi_asv_nmds_keep.list)[i]
  namevar <- i
  temp_order <- which(to_plot == {i})
  temp_df <- transformed_usvi_asv_nmds_keep.list[[i]][["points"]]
  temp_stress <-  transformed_usvi_asv_nmds_keep.list[[i]][["stress"]] %>%
    round(., digits = 3)
  
  temp_trans <- gsub("^(\\w+)(_)(\\w+)$", "\\1", namevar)
  temp_meth <- stringr::str_extract_all(namevar, paste0(metrics, collapse = "|"))
  if(grepl("_clr", temp_trans)){
    temp_title <- gsub("_clr", " + CLR", temp_trans)
  } else {
    temp_title <- temp_trans
  }
  if(grepl("relabund", temp_title)){
    temp_title <- gsub("relabund", "Relative abundance", temp_title)
  }
  
  
  temp_title <- paste0(temp_meth, " on ", temp_title) %>%
    # temp_title <- stringr::str_replace_all(temp_title, "_", " ") %>%
    stringr::str_to_sentence(.)
  
  if(grepl("clr", temp_title, ignore.case = TRUE)){
    temp_title <- gsub("clr", "CLR", temp_title)
  }
  
  
  temp_nmds1 <- (ggplot(data = temp_df %>%
                          dplyr::filter(type == "sample") %>%
                          droplevels,
                        aes(x = NMDS1, y = NMDS2))
                 + theme_bw()
                 + geom_point(aes(shape = sampling_time, fill = site), color = "black",
                              size = 3, stroke = 1, alpha = 1)
                 + scale_discrete_manual(aesthetics = c("shape"), values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA), drop = TRUE)
                 + scale_discrete_manual(aesthetics = c("fill"), values = site_colors, labels = site_lookup, breaks = names(site_lookup),drop = TRUE)
                 + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                         panel.background = element_blank(), 
                         panel.border = element_rect(fill = "NA", colour = "grey30"),
                         panel.grid = element_blank(),
                         plot.title = element_text(size = rel(0.9)),
                         plot.subtitle = element_text(size = rel(0.7)),
                         legend.position = "right",
                         legend.key = element_blank(),
                         legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                         legend.text = element_text(size = 12, colour = "grey30"))
                 # + facet_wrap(. ~ setome, drop = TRUE, scales = "free", axes = "all", shrink = FALSE,
                 #              labeller = global_labeller)
                 + scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS2")
                 + scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS1")
                 + guides(color = "none",
                          fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                              override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                          shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                               override.aes = list(color = "black", stroke = 1, size = 2)))
                 + ggtitle(paste0("NMDS using ", temp_title), 
                           subtitle = paste0("(stress: ", temp_stress, ")"))
  )
  assign(paste0("g_transformed_usvi_asv_nmds_", temp_order, "_", namevar), temp_nmds1, envir = .GlobalEnv)
  rm(temp_nmds1)
  rm(namevar)
  rm(temp_title)
  rm(temp_stress)
  rm(temp_df)
}


gpatch <- apropos("g_transformed_usvi_asv_nmds_*", mode = "any") %>%
  lapply(., get0) %>%
  purrr::reduce(`+`)
gpatch <- gpatch + patchwork::plot_layout(nrow = 3, guides = "collect") + patchwork::plot_annotation(tag_levels = "A")
gpatch

ggsave(paste0(projectpath, "/", "transformed_usvi_asv_nmds-", Sys.Date(), ".png"),
       gpatch,
       width = 16, height = 10, units = "in")




#how to test how different the communities are?

# temp_gof <- vegan::goodness(transformed_usvi_asv_nmds.list[[1]])
# temp_gof <- vegan::stressplot(transformed_usvi_asv_nmds.list[[4]]) 
#plots a graph with a non-metric fit and R2, and a linear model with R2

# temp_res <- vegan::protest(transformed_usvi_asv_nmds.list[[1]], transformed_usvi_asv_nmds.list[[2]], 
#                            scores = "sites")
# 
# temp_res[["ss"]] # Procrustes Sum of Squares (m12 squared)
# temp_res[["t0"]] # Correlation in a symmetric Procrustes rotation
# temp_res[["signif"]] #Significance
# vegan::protest(transformed_usvi_asv_nmds.list[[1]], transformed_usvi_asv_nmds.list[[2]], 
#                scores = "sites")
# # Procrustes Sum of Squares (m12 squared):        0.121 
# # Correlation in a symmetric Procrustes rotation: 0.9376 
# # Significance:  0.001 
# # 
# # Permutation: free
# # Number of permutations: 999

transformations_to_test <- combn(to_plot, 2) %>%
  setNames(., c("V1", "V2")) %>%
  tibble::as_tibble(.) %>%
  data.table::transpose(.)
transformations_procrustes.df <- NULL
# transformations_procrustes.df <- matrix(nrow = nrow(transformations_to_test), ncol = 5)

for(i in seq(nrow(transformations_to_test))){
  var1 <- transformations_to_test[i, "V1"]
  var2 <- transformations_to_test[i, "V2"]
  temp_res <- vegan::protest(transformed_usvi_asv_nmds.list[[{var1}]], transformed_usvi_asv_nmds.list[[{var2}]], 
                             scores = "sites")
  temp_df <- c("ss", "t0", "signif") %>%
    map(., ~purrr::pluck(temp_res, .x)) %>%
    setNames(., c("ss", "t0", "signif")) %>%
    bind_rows(.) %>%
    dplyr::mutate(V1 = {var1},
                  V2 = {var2} ) %>%
    dplyr::relocate(V1, V2) %>%
    droplevels
  if(exists("transformations_procrustes.df", envir = .GlobalEnv)){
    transformations_procrustes.df <- bind_rows(transformations_procrustes.df, temp_df)  
  } else {
    transformations_procrustes.df <- temp_df
  }
  
}

# temp_res[["ss"]] # Procrustes Sum of Squares (m12 squared)
# temp_res[["t0"]] # Correlation in a symmetric Procrustes rotation
# temp_res[["signif"]] #Significance

#which NMDS has the best correlation to relabund + MH?
transformations_procrustes.df %>%
  dplyr::arrange(desc(`t0`)) %>%
  dplyr::filter(if_any(c(V1, V2), ~grepl("relabund_horn", .x)))
# # A tibble: 6 × 5
# V1             V2                     ss    t0 signif
# <chr>          <chr>               <dbl> <dbl>  <dbl>
#   1 relabund_horn  relabund_euclidean 0.0706 0.964  0.001
# 2 relabund_bray  relabund_horn      0.121  0.938  0.001
# 3 clr_chord      relabund_horn      0.505  0.703  0.001
# 4 clr_euclidean  relabund_horn      0.518  0.694  0.001
# 5 rclr_euclidean relabund_horn      0.584  0.645  0.001
# 6 rclr_horn      relabund_horn      0.974  0.162  0.277
transformations_procrustes.df %>%
  dplyr::arrange(desc(`t0`)) %>%
dplyr::filter(if_any(c(V1, V2), ~grepl("rclr_euclidean", .x)))
# A tibble: 7 × 5
# V1                 V2                ss     t0 signif
# <chr>              <chr>          <dbl>  <dbl>  <dbl>
#   1 clr_chord          rclr_euclidean 0.247 0.868   0.001
# 2 clr_euclidean      rclr_euclidean 0.307 0.832   0.001
# 3 relabund_bray      rclr_euclidean 0.458 0.736   0.001
# 4 relabund_horn      rclr_euclidean 0.584 0.645   0.001
# 5 relabund_euclidean rclr_euclidean 0.643 0.598   0.001
# 6 rclr_euclidean     rclr_horn      0.984 0.125   0.589
# 7 clr_bray           rclr_euclidean 0.991 0.0939  0.78 

# vegan::protest(transformed_usvi_asv_nmds.list[["relabund_horn"]], 
#                # transformed_usvi_asv_nmds.list[["rclr_horn"]], 
#                # transformed_usvi_asv_nmds.list[["relabund_euclidean"]], 
#                transformed_usvi_asv_nmds.list[["relabund_bray"]], 
#                # scores = "species")
#                scores = "sites")
# 
# 
# transformation_procrustes_1 <- vegan::protest(transformed_usvi_asv_nmds.list[["relabund_horn"]], 
#                                               # transformed_usvi_asv_nmds.list[["rclr_horn"]], 
#                                               # transformed_usvi_asv_nmds.list[["relabund_euclidean"]], 
#                                               transformed_usvi_asv_nmds.list[["relabund_bray"]], 
#                                               # scores = "species")
#                            scores = "sites")
# plot(transformation_procrustes_1)
# 
# #a plot with the two least significant/weakest correlated NMDS procrustes results looks like a pom pom ball
# #a plot with the two most significant/strongest correlated NMDS procrustes results looks like a calm wind direction map
# 
# 
# 
# #the plot with Morisita Horn on relative abundance and Euclidean on rCLR looks like the wind direction map from a windier day
# transformation_procrustes_2 <- vegan::protest(transformed_usvi_asv_nmds.list[["relabund_horn"]], 
#                                               transformed_usvi_asv_nmds.list[["rclr_euclidean"]],
#                                               # transformed_usvi_asv_nmds.list[["relabund_euclidean"]], 
#                                               scores = "sites")
# plot(transformation_procrustes_2)
# 
# 
# transformations_procrustes.df %>%
#   dplyr::arrange(desc(`t0`)) %>%
#   dplyr::filter(if_any(c(V1, V2), ~grepl("rclr", .x)))
#   # dplyr::filter(`signif` < 0.05)
# # # A tibble: 11 × 5
# # V1             V2                    ss    t0 signif
# # <chr>          <chr>              <dbl> <dbl>  <dbl>
# #   1 clr_chord      rclr_euclidean     0.245 0.869  0.001
# # 2 clr_euclidean  rclr_euclidean     0.307 0.832  0.001
# # 3 rclr_euclidean relabund_bray      0.458 0.736  0.001
# # 4 rclr_euclidean relabund_horn      0.584 0.645  0.001
# # 5 rclr_euclidean relabund_euclidean 0.643 0.598  0.001
# # 6 rclr_horn      relabund_horn      0.974 0.162  0.277
# # 7 rclr_horn      relabund_bray      0.976 0.156  0.285
# # 8 rclr_horn      relabund_euclidean 0.977 0.152  0.36 
# # 9 clr_chord      rclr_horn          0.977 0.151  0.336
# # 10 clr_euclidean  rclr_horn          0.983 0.130  0.462
# # 11 rclr_euclidean rclr_horn          0.984 0.125  0.579
# 
# #if we wanted "species" scores... need to use either relabund + bray, horn, or euclidean
# #Euclidean on rCLR was best correlated to Bray-Curtis on Relative abundance
# 
# 
# #the plot with Morisita Horn on relative abundance and Euclidean on rCLR looks like the wind direction map from a windier day
# transformation_procrustes_3 <- vegan::protest(transformed_usvi_asv_nmds.list[["rclr_euclidean"]],
#                                               transformed_usvi_asv_nmds.list[["rclr_horn"]], 
#                                               # transformed_usvi_asv_nmds.list[["relabund_bray"]], 
#                                               # transformed_usvi_asv_nmds.list[["relabund_euclidean"]], 
#                                               scores = "sites")
# plot(transformation_procrustes_3)
# 
# 
#       
# gpatch2 <- (plot(transformation_procrustes_1) + ggtitle("BC vs MH on relative abundance")) + (plot(transformation_procrustes_2) + ggtitle("rCLR + Euclidean vs relative abundance + MH")) + (plot(transformation_procrustes_3) + ggtitle("Euclidean vs BC on rCLR"))
# plot(gpatch2)
# 
# #if we wanted to brute force "species" in the NMDS... transpose the input community matrix.
# # # transformed_usvi_asv_nmds.list[["rclr_euclidean"]]
# # 
# # temp_transformed.tbl <- transformations_usvi_asv.list %>%
# #   purrr::pluck( "rclr" ) %>%
# #   as.matrix(.) %>%
# #   t(.)
# # #this is a long, process-intense function... don't do it yet.
# # # temp_nmds.list <- temp_transformed.tbl %>%
# # #   vegan::metaMDS(., distance = "euclidean",trymax = 100, weakties = FALSE, autotransform = FALSE, 
# # #                  # tidy = FALSE,
# # #                  tidy = TRUE,
# # #                  wascores = TRUE)
# # # temp_nmds.df <- c("stress", "tries", "points", "species") %>%
# # #   map(., ~purrr::pluck(temp_nmds.list, .x)) %>%
# # #   setNames(., c("stress", "tries", "points", "species"))
# # # temp_nmds.df <- temp_nmds.df %>%
# # #   purrr::modify_at("points", ~.x %>%
# # #                      tibble::as_tibble(rownames = "label") %>%
# # #                      dplyr::rename(NMDS1 = "MDS1", NMDS2 = "MDS2") %>% 
# # #                      dplyr::mutate(type = "sample") %>%
# # #                      dplyr::left_join(., (metadata %>%
# # #                                             dplyr::filter(grepl("seawater", sample_type)) %>%
# # #                                             dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, replicate) %>%
# # #                                             droplevels),
# # #                                       by = c("label" = "sample_id")) %>%
# # #                      dplyr::mutate(site = factor(site, levels = names(site_lookup)),
# # #                                    sampling_time = factor(sampling_time, levels = names(sampling_time_lookup))))
# # # if(length(temp_nmds.df[["species"]]) > 1){
# # #   temp_nmds.df <- temp_nmds.df %>%
# # #     purrr::modify_at("species", ~.x %>%
# # #                        tibble::as_tibble(rownames = "label") %>%
# # #                        dplyr::rename(NMDS1 = "MDS1", NMDS2 = "MDS2") %>%
# # #                        dplyr::mutate(type = "asv_id") %>%
# # #                        dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"),
# # #                                         by = c("label" = "asv_id")) %>%
# # #                        droplevels)
# # # } else {
# # #   temp_nmds.df <- temp_nmds.df %>%
# # #     purrr::discard_at("species")
# # # }
# 
# #nah. probably better to just use the relabund-provided scores




#take these:
plot_procrustes.df <- data.frame(V1 = c("relabund_horn", "relabund_horn", "relabund_horn", "rclr_euclidean"),
                                 V2 = c("relabund_bray", "clr_euclidean", "rclr_euclidean", "rclr_horn"))
plot_procrustes.list <- NULL
for(i in seq(nrow(plot_procrustes.df))){
  var1 <- plot_procrustes.df[i, "V1"]
  var2 <- plot_procrustes.df[i, "V2"]
  namevar <- paste0({var2}, ":", {var1})
  temp_res <- vegan::protest(transformed_usvi_asv_nmds.list[[{var1}]], 
                             transformed_usvi_asv_nmds.list[[{var2}]], 
                             scores = "sites")
  temp_df <- c("Yrot", "X") %>%
    map(., ~purrr::pluck(temp_res, .x)) %>%
    setNames(., c("Yrot", "X"))
  
  temp_df <- temp_df %>%
    purrr::modify_at("Yrot", ~.x %>%
                       tibble::as_tibble(rownames = "sample_name") %>%
                       dplyr::rename(Y_NMDS1 = "V1", Y_NMDS2 = "V2") %>%
                       dplyr::mutate(`procrustes` = namevar) %>%
                       # dplyr::mutate(`ordination` = {var2}) %>%
                       # dplyr::mutate(`type` = "rotated_Y") %>%
                       droplevels)
  temp_df <- temp_df %>%
    purrr::modify_at("X", ~.x %>%
                       tibble::as_tibble(rownames = "sample_name") %>%
                       dplyr::rename(X_NMDS1 = "NMDS1", X_NMDS2 = "NMDS2") %>%
                       dplyr::mutate(`procrustes` = namevar) %>%
                       # dplyr::mutate(`ordination` = {var1} ) %>%
                       # dplyr::mutate(`type` = "target_X") %>%
                       droplevels)
  temp_df <- dplyr::full_join(temp_df[["X"]], temp_df[["Yrot"]]) %>%
    # temp_df <- temp_df %>%
    # dplyr::bind_rows(.) %>%
    dplyr::mutate(`target_X` = {var1}) %>%
    dplyr::mutate(`rotated_Y` = {var2}) %>%
    dplyr::relocate(`procrustes`) %>%
    # dplyr::relocate(ordination, type) %>%
    droplevels
  if(exists("plot_procrustes.list", envir = .GlobalEnv)){
    plot_procrustes.list <- append(plot_procrustes.list, list(temp_df))   %>%
      setNames(., c(names(plot_procrustes.list), namevar))
  } else {
    plot_procrustes.list <- list(temp_df) %>%
      setNames(., namevar)
  }
}
{
  temp_df <- plot_procrustes.list[[2]]
  namevar <- temp_df %>% dplyr::select(procrustes) %>% dplyr::distinct(.) %>% as.character(.)
  # temp_title <- gsub(":", " against ", namevar)
  temp_title <- gsub(":", ") against (", namevar)
  temp_title <- paste0("(", temp_title, ")")
  if(grepl("relabund", temp_title)){
    temp_title <- gsub("relabund", "Relative abundance", temp_title, ignore.case = TRUE)
  }
  if(grepl("clr_", temp_title)){
    temp_title <- gsub("clr", "CLR", temp_title)
  }
  temp_title <- temp_title %>%
    stringr::str_replace_all(., "bray", "Bray-Curtis") %>%
    stringr::str_replace_all(., "horn", "Morisita-Horn") %>%
    stringr::str_replace_all(., "euclidean", "Euclidean") %>%
    gsub("_", " + ", .)
  var1 <- temp_df[["target_X"]] %>% unique(.)
  var2 <- temp_df[["rotated_Y"]] %>% unique(.)
  temp_stress <- transformations_procrustes.df %>%
    dplyr::filter(if_all(c(V1, V2), ~grepl(paste0(c({var1}, {var2}), collapse = "|"), .x))) %>%
    dplyr::filter(grepl(paste0("^", {var1}), V1) &
                    grepl(paste0("^", {var2}), V2)) %>%
    # dplyr::select(`t0`, `signif`) %>%
    droplevels
  temp_sig <- temp_stress[["signif"]] %>% unique(.)
  temp_stress <- temp_stress[["t0"]] %>% unique(.) %>% round(., digits = 3)
  
  print(
    ggplot(data = temp_df)
    + geom_segment(aes(x = X_NMDS1, y = X_NMDS2, group = sample_name,
                       xend = Y_NMDS1, yend = Y_NMDS2),
                   linewidth = 1, color = "limegreen",
                   # arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                   lineend = "butt", linejoin = "round",
                   alpha = 0.7, show.legend = FALSE)
    + geom_point(aes(x = X_NMDS1, y = X_NMDS2, group = sample_name), color = "black", shape = 21, alpha = 0.8, size = 2, fill = "blue")
    + geom_point(aes(x = Y_NMDS1, y = Y_NMDS2, group = sample_name), color = "grey", shape = 25, alpha = 0.5, size = 2, fill = "limegreen")
    + scale_y_continuous(name = "NMDS1")
    + scale_x_continuous(name = "NMDS2")
    + theme_bw()
    + theme(panel.grid.major = element_blank(),
            # axis.text.x = element_text(angle = 90, hjust = 0.5),
            panel.grid.minor = element_blank(),
            strip.text.y = element_text(angle = 90))
    + ggtitle({temp_title},
              subtitle = paste0("stress: ", {temp_stress}, " (p = ", {temp_sig}, ")"))
    # subtitle = "(Original in blue) compared to (Rotation in green)")
  )
}

# i <- 2
for(i in seq(length(plot_procrustes.list))){
  
  temp_df <- plot_procrustes.list[[i]]
  namevar <- temp_df %>% dplyr::select(procrustes) %>% dplyr::distinct(.) %>% as.character(.)
  gname <- gsub(":", "_vs_", namevar)
  temp_title <- gsub(":", ") against (", namevar)
  temp_title <- paste0("(", temp_title, ")")
  if(grepl("relabund", temp_title)){
    temp_title <- gsub("relabund", "Relative abundance", temp_title, ignore.case = TRUE)
  }
  if(grepl("clr_", temp_title)){
    temp_title <- gsub("clr", "CLR", temp_title)
  }
  temp_title <- temp_title %>%
    stringr::str_replace_all(., "bray", "Bray-Curtis") %>%
    stringr::str_replace_all(., "horn", "Morisita-Horn") %>%
    stringr::str_replace_all(., "euclidean", "Euclidean") %>%
    gsub("_", " + ", .)
  var1 <- temp_df[["target_X"]] %>% unique(.)
  var2 <- temp_df[["rotated_Y"]] %>% unique(.)
  temp_stress <- transformations_procrustes.df %>%
    dplyr::filter(if_all(c(V1, V2), ~grepl(paste0(c({var1}, {var2}), collapse = "|"), .x))) %>%
    # dplyr::select(`t0`, `signif`)
    droplevels
  if(nrow(temp_stress) > 1){
    temp_stress <- temp_stress %>%
      dplyr::filter(grepl(paste0("^", "(", {var1}, "|", {var2}, ")"), V1) &
                      grepl(paste0("^", "(", {var1}, "|", {var2}, ")"), V2)) %>%
      droplevels
  }
  temp_sig <- temp_stress[["signif"]] %>% unique(.)
  temp_stress <- temp_stress[["t0"]] %>% unique(.) %>% round(., digits = 3)
  
  temp_g <- print(
    ggplot(data = temp_df)
    + geom_segment(aes(x = X_NMDS1, y = X_NMDS2, group = sample_name,
                       xend = Y_NMDS1, yend = Y_NMDS2),
                   linewidth = 1, color = "limegreen",
                   lineend = "butt", linejoin = "round",
                   alpha = 0.7, show.legend = FALSE)
    + geom_point(aes(x = X_NMDS1, y = X_NMDS2, group = sample_name), color = "black", shape = 21, alpha = 0.8, size = 2, fill = "blue")
    + geom_point(aes(x = Y_NMDS1, y = Y_NMDS2, group = sample_name), color = "grey", shape = 25, alpha = 0.5, size = 2, fill = "limegreen")
    + scale_y_continuous(name = "NMDS1")
    + scale_x_continuous(name = "NMDS2")
    + theme_bw()
    + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = rel(0.9)),
            plot.subtitle = element_text(size = rel(0.7)),
            strip.text.y = element_text(angle = 90))
    + ggtitle({temp_title},
              subtitle = paste0("stress: ", {temp_stress}, " (p = ", {temp_sig}, ")"))
    # subtitle = "(Original in blue) compared to (Rotation in green)")
  )
  assign(paste0("g_procrustes_", i, "_", gname), temp_g, envir = .GlobalEnv)
  rm(temp_df)
  rm(temp_title)
  rm(temp_stress)
  rm(temp_sig)
  rm(namevar) 
  rm(gname)
}

gpatch2 <- apropos("g_procrustes_*", mode = "any") %>%
  lapply(., get0) %>%
  purrr::reduce(`+`)
gpatch2 <- gpatch2 + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(tag_levels = "A")
gpatch2
ggsave(paste0(projectpath, "/", "transformed_usvi_asvs_procrustes-", Sys.Date(), ".png"),
       gpatch2,
       width = 10, height = 8, units = "in")


#listen. why use MH instead of BC on relative abundances?
#from https://academic.oup.com/biometrics/article/62/2/361/7321785?login=false
#
# both BC and mH "ignore the effect of missing species"
# "Bray–Curtis index performs well when sampling fractions are equal, as the theory predicts in Section 4. However, this index exhibits unreasonably large positive biases when sampling fractions are not equal. It is suggested that this index be used with caution."
# "Morisita–Horn index systematically underestimates similarity" and "is negatively biased by undersampling"
#but...
# "Setting aside our adjusted estimators, we see that [Morisita–Horn] index is the least sensitive to sample size and performs well in terms of relative biases. Its general performance is superior to our unadjusted indices (at the expense of insensitivity to rarer species), but inferior to the adjusted ones."

#why use MH:
#1. Sequence counts range 10x across samples
#2. Data are sparse: while there are over 10,000 ASVs in this dataset, between 169-827 ASVs were actually observed in nonzero abundances in each sample.

usvi_asv.tbl %>%
  colSums(.) %>%
  range(.)
# [1]  13915 114423
usvi_asv.tbl %>%
  apply(., 2, relabund) %>%
  quantile(., probs = c(0, 0.5, 1), na.rm = TRUE)
# 0%      50%     100% 
# 0.00000  0.00000 37.30608 
temp_tbl <- usvi_asv.tbl %>% 
  apply(., 2, na_if, 0) %>%
  apply(., 2, function(x) ifelse(is.na(x), 0, 1))
temp_tbl %>%
  apply(., 2, sum, na.rm = TRUE) %>%
  range(.)
# [1] 169 827

# Going back to PERMANOVA -------------------------------------------------


#does the results from adonis2 change when using a different ordination?

drop <- c("CINAR_BC_73")
keep <- c("sample_id", "metab_deriv_label", "sample_type", "site", "sampling_time", "sampling_day", "sampling_date", "depth", "site_type",
          "fcm_Prochlorococcus", "fcm_Synechococcus", "fcm_Picoeukaryotes", "fcm_Unpigmented_cells",
          "replicate", "grouping", "PAR", "lumens", "lux", "temp")

usvi_sw_asv.tbl <- usvi_prok_asvs.df %>%
  dplyr::filter(sample_ID %in% usvi_selected_metadata[["sample_id"]]) %>%
  dplyr::filter(!grepl("Metab_280", sample_ID)) %>% #drop the sample corresponding to CINAR_BC_73: Metab_280
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  droplevels %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  as.data.frame(.) %>% 
  # dplyr::select(colnames(usvi_sw_genus.tbl)) %>%
  droplevels


meta.microb <- usvi_selected_metadata %>%
  dplyr::select(intersect(colnames(usvi_selected_metadata), keep)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = colnames(usvi_sw_asv.tbl))) %>%
  dplyr::arrange(sample_id) %>%
  # tidyr::drop_na(.) %>%
  dplyr::mutate(grouping2 = paste0(grouping, ".", sampling_day)) %>%
  dplyr::mutate(across(c(sample_id, sample_type, site, sampling_time, sampling_day, metab_deriv_label, site_type, grouping, grouping2), ~factor(.x))) %>%
  dplyr::mutate(site = fct_relevel(site, "LB_seagrass"),
                grouping = fct_relevel(grouping, "LB_seagrass.dawn"),
                grouping2 = fct_relevel(grouping2, "LB_seagrass.dawn.Day2")) %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  droplevels

dist_mh_usvi_asv.d <- transformed_usvi_asv_relabund.tbl %>%
# dist_mh_usvi_asv.d <- usvi_sw_asv.tbl %>% 
#   dplyr::select(rownames(meta.microb)) %>%
#   t() %>%
  vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE)

dist_rclr_usvi_asv.d <- transformed_usvi_asv_rclr.tbl %>% 
  # dplyr::select(rownames(meta.microb)) %>%
  # t() %>%
  vegan::vegdist(., method = "euclidean", binary = FALSE, na.rm = TRUE)


with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 999,
                                 formula = dist_mh_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
                                 parallel = nthreads, by = "terms"))
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# vegan::adonis2(formula = dist_mh_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time, data = meta.microb, permutations = 999, method = "horn", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1  0.07471 0.06034 15.6914  0.001 ***
#   sampling_day                     3  0.10842 0.08757  7.5906  0.001 ***
#   site                             2  0.74508 0.60182 78.2471  0.001 ***
#   sampling_day:site                6  0.06811 0.05502  2.3843  0.038 *  
#   sampling_time:sampling_day:site 11  0.01796 0.01451  0.3430  0.940    
# Residual                        47  0.22377 0.18074                   
# Total                           70  1.23805 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

with(meta.microb, vegan::adonis2(data = meta.microb, method = "euclidean", permutations = 999,
                                 formula = dist_rclr_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
                                 parallel = nthreads, by = "terms"))
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# vegan::adonis2(formula = dist_rclr_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time, data = meta.microb, permutations = 999, method = "euclidean", by = "terms", parallel = nthreads)
# Df SumOfSqs      R2       F Pr(>F)    
# sampling_time                    1    466.2 0.02861  5.9322  0.003 ** 
#   sampling_day                     3   2800.8 0.17188 11.8791  0.001 ***
#   site                             2   7908.1 0.48531 50.3112  0.001 ***
#   sampling_day:site                6    757.0 0.04646  1.6054  0.106    
# sampling_time:sampling_day:site 11    668.9 0.04105  0.7737  0.749    
# Residual                        47   3693.8 0.22669                   
# Total                           70  16294.9 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis_rclr_asv.res <- vegan::adonis2(data = meta.microb, 
                                 method = "euclidean", permutations = 999,
                                 formula = dist_rclr_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
                                 parallel = nthreads, by = "terms")
adonis_mh_asv.res <- vegan::adonis2(data = meta.microb, 
                                      method = "horn", permutations = 999,
                                      formula = dist_mh_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
                                      parallel = nthreads, by = "terms")

vegan::mantel(dist_rclr_usvi_asv.d, dist_mh_usvi_asv.d, method = "spearman")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
#   vegan::mantel(xdis = dist_rclr_usvi_asv.d, ydis = dist_mh_usvi_asv.d,      method = "spearman") 
# 
# Mantel statistic r: 0.6922 
# Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0698 0.0909 0.1118 0.1299 
# Permutation: free
# Number of permutations: 999

vegan::mantel(dist_rclr_usvi_asv.d, dist_mh_usvi_asv.d, method = "pearson")
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# vegan::mantel(xdis = dist_rclr_usvi_asv.d, ydis = dist_mh_usvi_asv.d,      method = "pearson") 
# 
# Mantel statistic r: 0.7014 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0753 0.1026 0.1189 0.1355 
# Permutation: free
# Number of permutations: 999

#Using a Mantel test, the rCLR-distance matrix is significantly correlated 
#with the MH-distance matrix (p < 0.005, R ~ 0.70)


# Estimating % variation from adonis2 model -------------------------------

#previous calculation:
{
# permanova_asv_estimate_variation.df <- with(meta.microb, vegan::adonis2(data = meta.microb, method = "horn", permutations = 9999,
#                                                                         # strata = site,
#                                                                         formula = dist_usvi_asv.d ~ sampling_time + sampling_day + site + (sampling_day:site)/sampling_time,
#                                                                         parallel = nthreads, by = "terms"))
# permanova_asv_estimate_variation.df <- permanova_asv_estimate_variation.df %>%
#   tibble::rownames_to_column(var = "term") %>%
#   dplyr::mutate(term = as.character(term)) %>%
#   dplyr::mutate(term = dplyr::case_when(term == "sampling_time:sampling_day:site" ~ "pooled",
#                                         .default = term)) %>%
#   dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
#                                       .default = NA))
# permanova_asv_estimate_variation.df <- permanova_asv_estimate_variation.df %>%
#   dplyr::mutate(permanova_asv_estimate_variation.df %>%
#                   dplyr::filter(!grepl("Total", term)) %>%
#                   dplyr::summarise(TotalEMS = sum(MS))) %>%
#   dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/TotalEMS,
#                                           .default = NA)) %>%
#   dplyr::select(-c(TotalEMS, R2)) %>%
#   dplyr::mutate(sq_root = MS^(1/2)) %>%
#   dplyr::mutate(perc_variation = sq_root/(sum(sq_root, na.rm = TRUE))*100)
}


permanova_asv_estimate_variation.df <- list((adonis_rclr_asv.res %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::mutate(term = as.character(term)) %>%
  dplyr::mutate(term = dplyr::case_when(term == "sampling_time:sampling_day:site" ~ "pooled",
                                        .default = term)) %>%
  dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                      .default = NA))),
  (adonis_mh_asv.res %>%
     tibble::rownames_to_column(var = "term") %>%
     dplyr::mutate(term = as.character(term)) %>%
     dplyr::mutate(term = dplyr::case_when(term == "sampling_time:sampling_day:site" ~ "pooled",
                                           .default = term)) %>%
     dplyr::mutate(MS = dplyr::case_when(!grepl("Total", term) ~ SumOfSqs/Df,
                                         .default = NA)))) %>%
  setNames(., c("rCLR", "MH"))

temp_df <- permanova_asv_estimate_variation.df %>%
  map(., ~dplyr::mutate(.x %>%
                          dplyr::filter(!grepl("Total", term)) %>%
                          dplyr::summarise(TotalEMS = sum(MS), .groups = "drop")))
permanova_asv_estimate_variation.df <- names(permanova_asv_estimate_variation.df) %>%
  map(., ~permanova_asv_estimate_variation.df[[.x]] %>%
        # dplyr::full_join(., temp_df[[.x]]) %>%
        dplyr::rowwise(.) %>%
        dplyr::mutate(V_term = dplyr::case_when(!grepl("Total", term) ~ 100*MS/temp_df[[.x]]["TotalEMS"],
                                                .default = NA)) %>%
        # dplyr::bind_rows(.) %>%
        dplyr::select(-c(R2)) %>%
        # dplyr::select(-c(TotalEMS, R2)) %>%
        dplyr::mutate(sq_root = MS^(1/2)) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(perc_variation = sq_root/(sum(sq_root, na.rm = TRUE))*100) %>%
        tidyr::unnest(., "V_term") %>%
        dplyr::bind_rows(.) %>%
        droplevels) %>%
  setNames(., c("rCLR", "MH")) %>%
  dplyr::bind_rows(., .id = "ordination")

permanova_asv_estimate_variation.df
# # A tibble: 14 × 10
# ordination term                 Df   SumOfSqs      F `Pr(>F)`        MS TotalEMS sq_root perc_variation
# <chr>      <chr>             <dbl>      <dbl>  <dbl>    <dbl>     <dbl>    <dbl>   <dbl>          <dbl>
#   1 rCLR       sampling_time         1   466.      5.93     0.007   4.66e+2    8.30  21.6             15.1 
# 2 rCLR       sampling_day          3  2801.     11.9      0.001   9.34e+2   16.6   30.6             21.4 
# 3 rCLR       site                  2  7908.     50.3      0.001   3.95e+3   70.4   62.9             44.0 
# 4 rCLR       sampling_day:site     6   757.      1.61     0.099   1.26e+2    2.25  11.2              7.86
# 5 rCLR       pooled               11   669.      0.774    0.736   6.08e+1    1.08   7.80             5.46
# 6 rCLR       Residual             47  3694.     NA       NA       7.86e+1    1.40   8.87             6.20
# 7 rCLR       Total                70 16295.     NA       NA      NA         NA     NA               NA   
# 8 MH         sampling_time         1     0.0747 15.7      0.001   7.47e-2   14.9    0.273           21.2 
# 9 MH         sampling_day          3     0.108   7.59     0.001   3.61e-2    7.21   0.190           14.7 
# 10 MH         site                  2     0.745  78.2      0.001   3.73e-1   74.3    0.610           47.3 
# 11 MH         sampling_day:site     6     0.0681  2.38     0.049   1.14e-2    2.27   0.107            8.26
# 12 MH         pooled               11     0.0180  0.343    0.953   1.63e-3    0.326  0.0404           3.13
# 13 MH         Residual             47     0.224  NA       NA       4.76e-3    0.950  0.0690           5.35
# 14 MH         Total                70     1.24   NA       NA      NA         NA     NA               NA   


permanova_asv_estimate_variation.df %>%
  dplyr::arrange(term) %>%
  dplyr::select(term, ordination, perc_variation, `F`, `Pr(>F)`) %>%
  droplevels
# # A tibble: 14 × 5
# term              ordination perc_variation      F `Pr(>F)`
# <chr>             <chr>               <dbl>  <dbl>    <dbl>
#   1 Residual          rCLR                 6.20 NA       NA    
# 2 Residual          MH                   5.35 NA       NA    
# 3 Total             rCLR                NA    NA       NA    
# 4 Total             MH                  NA    NA       NA    
# 5 pooled            rCLR                 5.46  0.774    0.736
# 6 pooled            MH                   3.13  0.343    0.953
# 7 sampling_day      rCLR                21.4  11.9      0.001
# 8 sampling_day      MH                  14.7   7.59     0.001
# 9 sampling_day:site rCLR                 7.86  1.61     0.099
# 10 sampling_day:site MH                   8.26  2.38     0.049
# 11 sampling_time     rCLR                15.1   5.93     0.007
# 12 sampling_time     MH                  21.2  15.7      0.001
# 13 site              rCLR                44.0  50.3      0.001
# 14 site              MH                  47.3  78.2      0.001


# "Bi plot" ASVs on NMDS --------------------------------------------------

#what would the top SDA ASVs look like on the NMDS?
#grab the IDs of curated SDA ASVs:
# shared_sda_asvs_abund_idx
# shared_sda_asvs_idx

if(any(grepl("sda_asvs_compare", list.files(projectpath, pattern = "usvi_sda_asvs_compare-.*.RData")))){
  temp_file <- data.table::last(list.files(projectpath, pattern = "usvi_sda_asvs_compare-.*.RData"))
  load(paste0(projectpath, "/", temp_file))
  rm(temp_file)
} else {
  cli::cli_alert_info("Please read in results from DESeq2 and corncob to compare.")
}


#remember, only the Relative Abundance + (MH or Bray or Euclidean) NMDS generated coordinates for ASVs

usvi_sda_asvs_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_abund_idx) %>%
  droplevels %>%
  dplyr::summarise(numObs = length(asv_id), .by = "hold")
# # A tibble: 5 × 2
# hold        numObs
# <chr>        <int>
#   1 dawn           175
# 2 peak_photo     144
# 3 LB_seagrass     24
# 4 Yawzi           34
# 5 Tektite          2

usvi_sda_asvs_compare.df %>%
  dplyr::filter(asv_id %in% shared_sda_asvs_abund_idx) %>%
  dplyr::distinct(hold, asv_id, .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::summarise(numObs = length(asv_id), .by = "hold")
# # A tibble: 5 × 2
# hold        numObs
# <chr>        <int>
#   1 dawn            36
# 2 peak_photo      36
# 3 LB_seagrass     19
# 4 Yawzi           22
# 5 Tektite          2

temp_df <- transformed_usvi_asv_nmds_keep.list[["relabund_horn"]][["species"]] %>%
  # dplyr::filter(label %in% shared_sda_asvs_idx) %>%
  dplyr::filter(label %in% shared_sda_asvs_abund_idx) %>%
  dplyr::rename(asv_id = `label`) %>%
  dplyr::select(asv_id, NMDS1, NMDS2, taxonomy) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(Phylum, Order, asv_id) %>%
                     droplevels,
                   by = join_by(asv_id)) %>%
  droplevels

g_usvi_asv_nmds_relabund_horn_sda <- (g_transformed_usvi_asv_nmds_2_relabund_horn
           + geom_segment(data = temp_df,
                          aes(x = 0, y = 0, group = asv_id,
                              color = Order,
                              xend = NMDS1, yend = NMDS2),
                          linewidth = 0.5, 
                          # color = "limegreen",
                          arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                          lineend = "butt", linejoin = "round",
                          alpha = 0.7, show.legend = FALSE)
           + ggrepel::geom_text_repel(data = temp_df,
                                      aes(label = str_wrap(taxonomy, 25), x = NMDS1, y = NMDS2, color = Order),
                                      size = rel(2), min.segment.length = 0,
                                      # point.size = NA, segment.color = NA,
                                      arrow = arrow(type = "open", length = unit(0.02, "npc")), segment.color = "grey", 
                                      # colour = "black",
                                      direction = "both", seed = 48105, max.overlaps = 15,  show.legend =TRUE,
                                      force = 2, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"),
                                      fontface = "bold")
           # + geom_text(data = temp_df,
           #             aes(label = str_wrap(taxonomy, 25), x = NMDS1, y = NMDS2, color = Order),
           #             size = rel(2), check_overlap = FALSE, 
           #             position = "nudge",
           #             # colour = "black", 
           #             fontface = "bold")
           + guides(color = guide_legend(order = 3, ncol = 1, title = "Bacterial Order", direction = "vertical",
                                         override.aes = list(size = 3, label = "->")),
                    fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                        override.aes = list(color = "black", label = "", stroke = 1, shape = 21, size = 2)),
                    shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                         override.aes = list(color = "black", label = "", stroke = 1, size = 2)),
                    size = "none")
           + ggtitle(label = "NMDS using Morisita-Horn on relative abundance",
                     subtitle = "labels: 36 most abundant ASVs that were SDA between sites, and temporally")
           )
g_usvi_asv_nmds_relabund_horn_sda

ggsave(paste0(projectpath, "/", "g_usvi_asv_nmds_relabund_horn_sda-", Sys.Date(), ".png"),
       g_usvi_asv_nmds_relabund_horn_sda,
       width = 10, height = 8, units = "in")




# Read in metabolites -----------------------------------------------------

if(file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
  temp_list <- readr::read_rds(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
  list2env(temp_list, envir = .GlobalEnv)
  rm(temp_list)
} else {
  cli::cli_alert_warning("Please tidy the metabolomics datasets.")
}

# Rerun Correlations ------------------------------------------------------


#recalculate using rCLR transformed ASV table
to_import <- c("spearman.test.rclr.optA.list",
               "spearman.test.site.rclr.optA.list", 
               "spearman.test.site.time.rclr.optA.list")

for(file in to_import){
  if(!exists(file, envir = .GlobalEnv)){
    namevar <- file
    if(file.exists(paste0(projectpath, "/", namevar, ".rds"))){
      cli::cli_alert_info("Importing this dataset: {namevar}")
      temp_df <- readr::read_rds(paste0(projectpath, "/", namevar, ".rds"))
      assign(paste0(namevar), temp_df, envir = .GlobalEnv)
      rm(temp_df)
      to_import <- to_import[grep(namevar, to_import, value = FALSE, invert = TRUE)]
    } else {
      cli::cli_alert_warning("Please prepare this dataset: {namevar}")
    }
  }
}
if(length(to_import) > 0){
  
  usvi_metab.tbl <- usvi_metabolomics_long.df %>%
    dplyr::select(metabolites, adaptedDervLabel, concentration, LODflag) %>%
    dplyr::rename(simpleName = "metabolites", metab_deriv_label = "adaptedDervLabel", conc = "concentration") %>%
    dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
    #   dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
    dplyr::filter(LODflag == 0) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           # dplyr::filter(sample_id %in% colnames(usvi_metab_mat)) %>%
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
  
  usvi_asv_rclr.tbl <- usvi_asv.tbl %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    dplyr::select(rownames(usvi_metab.tbl)) %>%
    vegan::decostand(., method = "rclr", 2, na.rm = TRUE, impute = TRUE) %>%
    t()
  
  usvi_asv_rclr.tbl <- usvi_asv_rclr.tbl[rownames(usvi_metab.tbl),]
  
  # if(!exists("spearman.test.rclr.optA.list", envir = .GlobalEnv)){
  #   if(file.exists(paste0(projectpath, "/", "spearman.test.rclr.optA.list", ".rds"))){
  #     spearman.test.rclr.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.rclr.optA.list", ".rds"))
  #   }
  # }
  if(any(grepl("spearman.test.rclr.optA.list", to_import))){
    #all samples, all times:
    {
      spearman.test <- matrix(nrow = ncol(usvi_asv_rclr.tbl), ncol = ncol(usvi_metab.tbl))
      colnames(spearman.test) <- colnames(usvi_metab.tbl)
      rownames(spearman.test) <- colnames(usvi_asv_rclr.tbl)
      
      spearman.test.rho <- spearman.test
      spearman.test.n <- spearman.test  
      
      y <- length(colnames(spearman.test))
      for(j in seq_len(y)){
        vector_metab <- usvi_metab.tbl[, j]
        for(i in seq_len(nrow(spearman.test))){
          vector_microb <- usvi_asv_rclr.tbl[,i]
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
        # dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
        dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
        dplyr::arrange(asv_id, simpleName) %>%
        dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
        # tidyr::drop_na(contains("padj")) %>%
        dplyr::mutate(label = signif(estimate, digits = 2)) %>%
        dplyr::ungroup(.) %>%
        dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
        droplevels
      
      spearman.test.rclr.optA.list <- list(spearman.test.df,
                                           spearman.test, spearman.test.rho, spearman.test.n,
                                           # dend_asv, dend_metab,
                                           padj_cutoff) %>%
        setNames(., c("spearman.test.df",
                      "spearman.test", "spearman.test.rho", "spearman.test.n",
                      # "dend_asv", "dend_metab", 
                      "padj_cutoff"))
    }
    readr::write_rds(spearman.test.rclr.optA.list, paste0(projectpath, "/", "spearman.test.rclr.optA.list", ".rds"), compress = "gz")
    
  }
  # 
  # if(!exists("spearman.test.site.rclr.optA.list", envir = .GlobalEnv)){
  #   if(file.exists(paste0(projectpath, "/", "spearman.test.site.rclr.optA.list", ".rds"))){
  #     spearman.test.site.rclr.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.rclr.optA.list", ".rds"))
  #   }
  # }
  if(any(grepl("spearman.test.site.rclr.optA.list", to_import))){
    #by site, option A:
    {
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
              dplyr::inner_join(., usvi_asv_rclr.tbl %>%
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
      
      
      spearman.test.site.rclr.optA.list <- lapply(apropos("^spearman.test.siteA.*$", mode = "list"),
                                                  get) %>%
        setNames(., names(usvi_asv_site.list))
      
    }
    readr::write_rds(spearman.test.site.rclr.optA.list, paste0(projectpath, "/", "spearman.test.site.rclr.optA.list", ".rds"), compress = "gz")
    
  }
  # 
  # if(!exists("spearman.test.site.time.rclr.optA.list", envir = .GlobalEnv)){
  #   if(file.exists(paste0(projectpath, "/", "spearman.test.site.time.rclr.optA.list", ".rds"))){
  #     spearman.test.site.time.rclr.optA.list <- readr::read_rds(paste0(projectpath, "/", "spearman.test.site.time.rclr.optA.list", ".rds"))
  #   }
  # }
  if(any(grepl("spearman.test.site.time.rclr.optA.list", to_import))){
    
    #by site and time, option A:
    {
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
              dplyr::inner_join(., usvi_asv_rclr.tbl %>%
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
        
        assign(paste0("spearman.test.siteA.rclr.", namevar, ".list"), spearman.test.site.time.list, envir = .GlobalEnv)
        rm(spearman.test.site.time.list)
        rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
      }
      
      
      spearman.test.site.time.rclr.optA.list <- lapply(apropos("^spearman.test.siteA.rclr.*$", mode = "list"),
                                                       get) %>%
        setNames(., names(usvi_asv_site.time.list))
    }
    
    readr::write_rds(spearman.test.site.time.rclr.optA.list, paste0(projectpath, "/", "spearman.test.site.time.rclr.optA.list", ".rds"), compress = "gz")
    
  }
  
  rm(list = apropos("^(spearman.test.siteA.)(.*)$", mode = "list"))
  rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
  rm(list = apropos("^(temp_)(.*)(.tbl)$", mode = "list"))
  rm(list = apropos("^(vector_m)(.*)(b)$", mode = "list"))
  rm(list = apropos("^(vector_m)(.*)(b_na)$", mode = "list"))
  
}



# Analyze results ---------------------------------------------------------

if(!exists("spearman.test.site.time.rclr.full.df", envir = .GlobalEnv)){
  if(any(grepl("spearman_rclr_full", list.files(projectpath, pattern = "usvi_.*.RData")))){
    load(list.files(projectpath, "usvi_spearman_rclr_full-.*.RData")[1])
    #spearman.test.site.time.rclr.full.df, spearman.test.site.rclr.full.df, spearman.test.rclr.full.df
    
  } else {
    spearman.test.site.time.rclr.full.df <- spearman.test.site.time.rclr.optA.list %>%
      map(.,
          ~.x[["spearman.test"]] %>%
            tibble::as_tibble(rownames = "asv_id") %>%
            tidyr::pivot_longer(., cols = -c("asv_id"),
                                names_to = "simpleName",
                                values_to = "p_value") %>%
            bind_rows(.)) %>%
      bind_rows(., .id = "grouping") %>%
      dplyr::left_join(., spearman.test.site.time.rclr.optA.list %>%
                         map(.,
                             ~.x[["spearman.test.rho"]] %>%
                               tibble::as_tibble(rownames = "asv_id") %>%
                               tidyr::pivot_longer(., cols = -c("asv_id"),
                                                   names_to = "simpleName",
                                                   values_to = "estimate") %>%
                               bind_rows(.)) %>%
                         bind_rows(., .id = "grouping"),
                       by = join_by(asv_id, grouping, simpleName)) %>%
      dplyr::mutate(padj_bh = p.adjust(p_value, "BH")) %>%
      droplevels
    
    View(spearman.test.site.time.rclr.full.df)
    
    #compare to the site-specific correlations
    spearman.test.site.rclr.full.df <- spearman.test.site.rclr.optA.list %>%
      map(.,
          ~.x[["spearman.test"]] %>%
            tibble::as_tibble(rownames = "asv_id") %>%
            tidyr::pivot_longer(., cols = -c("asv_id"),
                                names_to = "simpleName",
                                values_to = "p_value") %>%
            bind_rows(.)) %>%
      bind_rows(., .id = "grouping") %>%
      dplyr::left_join(., spearman.test.site.rclr.optA.list %>%
                         map(.,
                             ~.x[["spearman.test.rho"]] %>%
                               tibble::as_tibble(rownames = "asv_id") %>%
                               tidyr::pivot_longer(., cols = -c("asv_id"),
                                                   names_to = "simpleName",
                                                   values_to = "estimate") %>%
                               bind_rows(.)) %>%
                         bind_rows(., .id = "grouping"),
                       by = join_by(asv_id, grouping, simpleName)) %>%
      dplyr::mutate(padj_bh = p.adjust(p_value, "BH"))
    
    #compare to the all-sample correlations
    spearman.test.rclr.full.df <- spearman.test.rclr.optA.list %>%
      purrr::pluck(., "spearman.test") %>%
      tibble::as_tibble(rownames = "asv_id") %>%
      tidyr::pivot_longer(., cols = -c("asv_id"),
                          names_to = "simpleName",
                          values_to = "p_value") %>%
      droplevels %>%
      dplyr::left_join(., (spearman.test.rclr.optA.list %>%
                             purrr::pluck(., "spearman.test.rho") %>%
                             tibble::as_tibble(rownames = "asv_id") %>%
                             tidyr::pivot_longer(., cols = -c("asv_id"),
                                                 names_to = "simpleName",
                                                 values_to = "estimate") %>%
                             droplevels),
                       by = join_by(asv_id, simpleName)) %>%
      dplyr::mutate(padj_bh = p.adjust(p_value, "BH")) %>%
      dplyr::mutate(grouping = "all") %>%
      droplevels
    
    
    
    if(!any(grepl("spearman_rclr_full", list.files(projectpath, pattern = "usvi_.*.RData")))){
      save(spearman.test.site.time.rclr.full.df, spearman.test.site.rclr.full.df, spearman.test.rclr.full.df,
           file = paste0(projectpath, "/", "usvi_spearman_rclr_full-", Sys.Date(), ".RData"))
    }
    
    
  }
}



if(!exists("spearman.test.site.time.rclr.filtered.df", envir = .GlobalEnv)){
  
  if(any(grepl("spearman.rclr.sig.filtered.list", list.files(projectpath, pattern = "usvi_.*.rds")))){
    temp_list <- readr::read_rds(list.files(projectpath, pattern = "usvi_spearman.rclr.sig.filtered.list-.*.rds", full.names = TRUE)[1])
    list2env(temp_list, envir = .GlobalEnv)
    rm(temp_list)
    # spearman.test.rclr.filtered.lists <- list(spearman.test.site.time.rclr.filtered.df,
    #                                           spearman.test.site.rclr.filtered.df,
    #                                           spearman.test.rclr.filtered.df)
  } else {
    
    #recalculate FDRs
    padj_cutoff <- spearman.test.rclr.full.df %>%
      dplyr::select(p_value) %>%
      tibble::deframe(.) %>% na.omit(.) %>%
      unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
      quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
      setNames(., c("q_01", "q_025", "q_05", "q_10")) #get the possible p-adj cutoffs for different q-values
    
    
    spearman.test.rclr.df <- spearman.test.rclr.full.df %>%
      tidyr::drop_na(p_value) %>%
      dplyr::mutate(across(c(asv_id, simpleName), ~factor(.x))) %>%
      droplevels %>%
      dplyr::rowwise(.) %>%
      # dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
      #               padj_01 = dplyr::case_when(p_value <= padj_cutoff["q_01"] ~ p_value, .default = NA),
      #               padj_025 = dplyr::case_when(p_value <= padj_cutoff["q_025"] ~ p_value, .default = NA),
      #               padj_05 = dplyr::case_when(p_value <= padj_cutoff["q_05"] ~ p_value, .default = NA),
      #               padj_10 = dplyr::case_when(p_value <= padj_cutoff["q_10"] ~ p_value, .default = NA)) %>%
      dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                    padj_01 = dplyr::case_when(p_value <= padj_cutoff["q_01"] ~ p_value, .default = NA),
                    padj_025 = dplyr::case_when(p_value <= padj_cutoff["q_025"] ~ p_value, .default = NA),
                    padj_05 = dplyr::case_when(p_value <= padj_cutoff["q_05"] ~ p_value, .default = NA),
                    padj_10 = dplyr::case_when(p_value <= padj_cutoff["q_10"] ~ p_value, .default = NA)) %>%
      tidyr::drop_na(padj_10) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
      tidyr::drop_na(estimate)
    
    # temp_df <- spearman.test.rclr.df %>%
    spearman.test.rclr.df <- spearman.test.rclr.df %>%
      dplyr::rowwise(.) %>%
      dplyr::mutate(sig = dplyr::case_when(
        # !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
        !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
        !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
        !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
        !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
        .default = NA)) %>%
      dplyr::mutate(sig = dplyr::case_when(
        !is.na(padj_bh_05) & is.na(sig) ~ "vsig", #meaning that the adjusted p-value is below 0.05
        .default = sig)) %>%
      dplyr::mutate(sig = factor(sig)) %>%
      dplyr::arrange(asv_id, simpleName) %>%
      dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
      droplevels %>%
      bind_rows(.)
    
    
    
    #Now site-specific
    padj_cutoff <- spearman.test.site.rclr.full.df %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    
    spearman.test.site.rclr.df <- spearman.test.site.rclr.full.df %>%
      tidyr::drop_na(p_value) %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels)
    
    spearman.test.site.rclr.df <- spearman.test.site.rclr.df %>%
      imap(., ~.x %>%
             dplyr::mutate(across(c(asv_id, simpleName, grouping), ~factor(.x))) %>%
             droplevels %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                           padj_01 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_01"] ~ p_value, .default = NA),
                           padj_025 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_025"] ~ p_value, .default = NA),
                           padj_05 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_05"] ~ p_value, .default = NA),
                           padj_10 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_10"] ~ p_value, .default = NA)) %>%
             tidyr::drop_na(padj_10) %>%
             dplyr::ungroup(.) %>%
             dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                       abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
             tidyr::drop_na(estimate) %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(sig = dplyr::case_when(
               # !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
               !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
               .default = NA)) %>%
             dplyr::mutate(sig = dplyr::case_when(
               !is.na(padj_bh_05) & is.na(sig) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               .default = sig)) %>%
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
    
    #Now site- and time-specific
    padj_cutoff <- spearman.test.site.time.rclr.full.df %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    
    spearman.test.site.time.rclr.df <- spearman.test.site.time.rclr.full.df %>%
      tidyr::drop_na(p_value) %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels)
    spearman.test.site.time.rclr.df <- spearman.test.site.time.rclr.df %>%
      imap(., ~.x %>%
             dplyr::mutate(across(c(asv_id, simpleName, grouping), ~factor(.x))) %>%
             droplevels %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                           padj_01 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_01"] ~ p_value, .default = NA),
                           padj_025 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_025"] ~ p_value, .default = NA),
                           padj_05 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_05"] ~ p_value, .default = NA),
                           padj_10 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_10"] ~ p_value, .default = NA)) %>%
             tidyr::drop_na(padj_10) %>%
             dplyr::ungroup(.) %>%
             dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                       abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
             tidyr::drop_na(estimate) %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(sig = dplyr::case_when(
               # !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
               !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
               .default = NA)) %>%
             dplyr::mutate(sig = dplyr::case_when(
               !is.na(padj_bh_05) & is.na(sig) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               .default = sig)) %>%
             dplyr::mutate(sig = factor(sig)) %>%
             dplyr::arrange(asv_id, simpleName) %>%
             dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
             dplyr::ungroup(.) %>%
             dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
             dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                           sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
             droplevels) %>%
      bind_rows(.)
    
    spearman.test.rclr.filtered.df <- spearman.test.rclr.df %>%
      dplyr::rowwise(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      # dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
      dplyr::mutate(across(c(asv_id, simpleName, sig), ~factor(.x))) %>%
      droplevels
    
    spearman.test.site.rclr.filtered.df <- spearman.test.site.rclr.df %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
    
    spearman.test.site.time.rclr.filtered.df <- spearman.test.site.time.rclr.df %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
    
    # if(!any(grepl("spearman.rclr.sig.filtered.list", list.files(projectpath, pattern = "usvi_.*.rds")))){
    spearman.test.rclr.filtered.lists <- list(spearman.test.site.time.rclr.filtered.df,
                                              spearman.test.site.rclr.filtered.df,
                                              spearman.test.rclr.filtered.df) %>%
      setNames(., c("spearman.test.site.time.rclr.filtered.df",
                    "spearman.test.site.rclr.filtered.df",
                    "spearman.test.rclr.filtered.df"))
    readr::write_rds(spearman.test.rclr.filtered.lists, 
                     paste0(projectpath, "/", "usvi_spearman.rclr.sig.filtered.list-", Sys.Date(), ".rds"), 
                     compress = "gz")
    
    for(df in c("spearman.test.site.time.rclr.filtered.df",
                "spearman.test.site.rclr.filtered.df",
                "spearman.test.rclr.filtered.df")){
      temp_df <- get0(df, mode = "any")
      readr::write_delim(temp_df, paste0(projectpath, "/", "usvi_", df, "-", Sys.Date(), ".tsv"),
                         col_names = TRUE, delim = "\t")  
    }
    # }
  }
  
}

if(!exists("padj_cutoff_list", envir = .GlobalEnv)){
  if(any(grepl("spearman.test.rclr.padj_cutoff", list.files(projectpath, pattern = "usvi_.*.rds")))){
    padj_cutoff_list <- readr::read_rds(list.files(projectpath, "usvi_spearman.test.rclr.padj_cutoff-.*.rds")[1])
  } else {
  
  #save the calculated FDRs as well
  
  padj_cutoff <- spearman.test.rclr.full.df %>%
    dplyr::select(p_value) %>%
    tibble::deframe(.) %>% na.omit(.) %>%
    unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
    quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
    setNames(., c("q_01", "q_025", "q_05", "q_10")) #get the possible p-adj cutoffs for different q-values
  temp_df <- spearman.test.site.rclr.full.df %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          droplevels %>%
          dplyr::select(p_value) %>%
          tibble::deframe(.) %>% na.omit(.) %>%
          unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
          quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
          setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
  temp_df2 <- spearman.test.site.time.rclr.full.df %>%
    split(., f = .$grouping) %>%
    map(., ~.x %>%
          dplyr::select(p_value) %>%
          tibble::deframe(.) %>% na.omit(.) %>%
          unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
          quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
          setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
  
  padj_cutoff_list <- padj_cutoff %>%
    list(.) %>% 
    setNames(., c("all")) %>%
    append(., temp_df) %>%
    append(., temp_df2)
  readr::write_rds(padj_cutoff_list, paste0(projectpath, "/", "usvi_spearman.test.rclr.padj_cutoff-", Sys.Date(), ".rds"))
  
  padj_cutoff <- padj_cutoff_list %>%
    dplyr::bind_rows(., .id = "model")
  
  readr::write_delim(padj_cutoff, paste0(projectpath, "/", "usvi_spearman.test.rclr.padj_cutoff-", Sys.Date(), ".tsv"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
  
  }
  
}


if(!exists("spearman.test.site.time.rclr.df", envir = .GlobalEnv)){
  
    # temp_idx <- c("asv_id","simpleName", "p_value","estimate","padj_bh","grouping","padj_bh_05", 
    #               "padj_01","padj_025","padj_05","padj_10","sig")
    temp_drop <- c("filtered_estimate")
    
    spearman.test.rclr.df <- spearman.test.rclr.filtered.df %>%
      dplyr::select(-c(all_of(temp_drop))) %>%
      # dplyr::select(all_of(temp_idx)) %>%
      droplevels
    spearman.test.site.rclr.df <- spearman.test.site.rclr.filtered.df %>%
      dplyr::select(-c(all_of(temp_drop))) %>%
      # dplyr::select(all_of(temp_idx), "site") %>%
      droplevels
    spearman.test.site.time.rclr.df <- spearman.test.site.time.rclr.filtered.df %>%
      dplyr::select(-c(all_of(temp_drop))) %>%
      # dplyr::select(all_of(temp_idx), "site", "sampling_time") %>%
      droplevels
    
    
    # spearman.test.rclr.sig.lists <- list(spearman.test.site.time.rclr.df,
    #                                      spearman.test.site.rclr.df,
    #                                      spearman.test.rclr.df) %>%
    #   setNames(., c("spearman.test.site.time.rclr.df",
    #                 "spearman.test.site.rclr.df",
    #                 "spearman.test.rclr.df"))
    # readr::write_rds(spearman.test.rclr.sig.lists, 
    #                  paste0(projectpath, "/", "usvi_spearman.rclr.sig.list-", Sys.Date(), ".rds"), 
    #                  compress = "gz")
}


# Parse through Spearman results ------------------------------------------

#now combine themm
spearman.rclr.all.rho.p.df <- spearman.test.rclr.full.df %>%
  dplyr::filter(p_value <= 0.1) %>%
  dplyr::distinct(asv_id, simpleName) %>%
  droplevels %>%
  dplyr::inner_join(., spearman.test.site.time.rclr.full.df %>%
                      dplyr::filter(p_value <= 0.1) %>%
                      dplyr::distinct(asv_id, simpleName) %>%
                      droplevels, 
                    by = join_by(asv_id, simpleName)) %>%
  droplevels %>%
  dplyr::inner_join(., spearman.test.site.rclr.full.df %>%
                      dplyr::filter(p_value <= 0.1) %>%
                      dplyr::distinct(asv_id, simpleName) %>%
                      droplevels, 
                    by = join_by(asv_id, simpleName)) %>%
  droplevels %>%
  dplyr::left_join(., spearman.test.rclr.full.df %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels, 
                   by = join_by(asv_id, simpleName)) %>%
  dplyr::bind_rows(., spearman.test.site.time.rclr.full.df %>%
                     dplyr::filter(p_value <= 0.1) %>%
                     dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels) %>%
  dplyr::bind_rows(., spearman.test.site.rclr.full.df %>%
                     dplyr::filter(p_value <= 0.1) %>%
                     dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels) %>%
  droplevels

spearman.rclr.all.rho.p.df <- spearman.rclr.all.rho.p.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, simpleName, grouping) %>%
  dplyr::filter(grepl("all", grouping), .by = c(asv_id, simpleName)) %>%
  dplyr::distinct(asv_id, simpleName) %>%
  dplyr::left_join(., spearman.rclr.all.rho.p.df,
                   by = join_by(asv_id, simpleName)) %>%
  tidyr::pivot_wider(., id_cols = c(asv_id, simpleName, grouping),
                     names_from = "metric",
                     values_from = "value") %>%
  dplyr::mutate(grouping = fct_relevel(grouping, "all", after = Inf)) %>%
  droplevels

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

spearman.test.rclr.bh.dist.df <- bind_rows(spearman.test.rclr.df, spearman.test.site.rclr.df) %>%
  bind_rows(., spearman.test.site.time.rclr.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  tidyr::drop_na(padj_10) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  # tidyr::drop_na(padj_bh_05) %>%
  # dplyr::ungroup(.) %>%
  # dplyr::arrange(desc(padj_bh_05)) %>%
  dplyr::arrange(desc(padj_bh_05), .by_group = TRUE) %>%
  dplyr::select(grouping, site, sampling_time, p_value, padj_10, padj_05, padj_025, padj_01, padj_bh_05) %>%
  dplyr::distinct(grouping, site, sampling_time, .keep_all = TRUE) %>%
  droplevels

#plot the distribution of Spearman rho's when using all samples, compared to site-specific
spearman.test.rclr.rho.dist.df <- bind_rows(spearman.test.rclr.filtered.df, spearman.test.site.rclr.filtered.df) %>%
  bind_rows(., spearman.test.site.time.rclr.filtered.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  dplyr::distinct(asv_id, simpleName, grouping, site, sampling_time, sig, .keep_all = TRUE) %>%
  dplyr::select(asv_id, simpleName, grouping, site, sampling_time, estimate, sig) %>%
  # dplyr::mutate(sig = factor(sig, levels = names(padj_cutoff_labels))) %>%
  dplyr::group_by(asv_id, simpleName,  grouping,site, sampling_time, sig) %>%
  # dplyr::distinct(estimate, .keep_all = TRUE) %>%
  droplevels

spearman.test.rclr.rho.dist.df %>%
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






spearman.test.rclr.rho.dist.summary.df <- spearman.test.rclr.rho.dist.df %>%
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
spearman.test.rclr.rho.dist.summary.df %>%
  dplyr::filter(grepl("all", sampling_time))
# # A tibble: 16 × 5
# grouping    site         sampling_time sig      num_obs
# <fct>       <fct>        <fct>         <fct>      <int>
#   1 LB_seagrass Seagrass     all           sig_q05     2938
# 2 LB_seagrass Seagrass     all           sig_q025     822
# 3 LB_seagrass Seagrass     all           sig_q01        0
# 4 LB_seagrass Seagrass     all           maybe       9042
# 5 Yawzi       Yawzi Reef   all           sig_q05      673
# 6 Yawzi       Yawzi Reef   all           sig_q025    2056
# 7 Yawzi       Yawzi Reef   all           sig_q01     1066
# 8 Yawzi       Yawzi Reef   all           maybe       8016
# 9 Tektite     Tektite Reef all           sig_q05     1926
# 10 Tektite     Tektite Reef all           sig_q025    2229
# 11 Tektite     Tektite Reef all           sig_q01      893
# 12 Tektite     Tektite Reef all           maybe       5507
# 13 all         all          all           sig_q05     4633
# 14 all         all          all           sig_q025    3286
# 15 all         all          all           sig_q01     2146
# 16 all         all          all           maybe      12213
print(spearman.test.rclr.rho.dist.summary.df, n = 50)





# Compare rCLR+Spearman to relabund+Spearman ------------------------------


for(df in c("spearman.test.site.time.filtered.df",
            "spearman.test.site.filtered.df",
            "spearman.test.filtered.df")){
  namevar <- rev(list.files(projectpath, {df}))[1]
  temp_df <- readr::read_delim(paste0(projectpath, "/", namevar),
                     col_names = TRUE, delim = "\t", show_col_types = FALSE)
  assign(df,temp_df, envir = .GlobalEnv)
  rm(temp_df)
  rm(namevar)
}


corr_fdr01_MH_site_sig <- spearman.test.site.filtered.df %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig), ~factor(.x))) %>%
  filter(sig %in% c("sig_q01","vsig")) %>%
  droplevels

# spearman.test.site.filtered.df %>%
#   dplyr::mutate(across(c(asv_id, simpleName, grouping, sig), ~factor(.x))) %>%
#   filter(sig %in% c("sig_q01","vsig")) %>% # [1] 1795
#   filter(sig %in% c("sig_q01", "sig_q025")) %>% # [1] 1757
#   nrow(.)


#ubiquitous ASVs: the unique group of ASVs with 10 or more correlations with metabolites
#if we see that ASV-metabolite correlation in multiple sites:
corr_fdr01_MH_site_sig %>%
  distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
  count(asv_id, name = "num_correlations") %>%   # Count occurrences
  arrange(desc(num_correlations))%>% 
  distinct()%>%
  dplyr::filter(num_correlations >= 10) %>%
  nrow(.)
# [1] 21

corr_fdr01_MH_site_sig %>%
  count(asv_id, name = "num_correlations") %>%   # Count occurrences
  arrange(desc(num_correlations))%>% 
  distinct()%>%
  dplyr::filter(num_correlations >= 10) %>%
  nrow(.)
# [1] 24

#there are 3 ASVs that are correlated with 10 or more metabolites, in 2 or more sites:
corr_fdr01_MH_site_freq_sig_df <- dplyr::full_join((corr_fdr01_MH_site_sig %>%
                                                   distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
                                                   count(asv_id, name = "num_unique_correlations") %>%   # Count occurrences
                                                   arrange(desc(num_unique_correlations))),
                                                (corr_fdr01_MH_site_sig %>%
                                                   count(asv_id, name = "num_correlations") %>%   # Count occurrences
                                                   arrange(desc(num_correlations)))) %>%
  dplyr::filter(num_correlations >= 10) %>%
  droplevels
corr_fdr01_MH_site_freq_sig_df %>%
  dplyr::filter(num_unique_correlations < 10 ) %>%
  droplevels %>%
  dplyr::distinct(asv_id) %>%
  tibble::deframe(.)
# [1] ASV_00024 ASV_00026 ASV_00028

corr_fdr01_MH_site_sig %>%
  distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
  count(simpleName,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
  arrange(desc(num_correlations))%>%
  head()
# # A tibble: 6 × 2
# simpleName         num_correlations
# <chr>                         <int>
#   1 pantothenic acid                216
# 2 ciliatine                       141
# 3 cysteate                        100
# 4 homoserine betaine               96
# 5 kynurenine                       81
# 6 GABA                             80


#This is our overall count of significant and strong (abs(rho) >= 0.5) ASV-Metabolite correlations:
corr_fdr01_MH_site_sig_unique_venn_list <- corr_fdr01_MH_site_sig %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::distinct(corrlink) %>%
        tibble::deframe(.))


#now to count how many UNIQUE ASV-metab associations there are across all sites
#first, see how many metabolite-ASV correlations are observed in 2 or more sites ("cosmopolitan")
corr_fdr01_MH_site_sig_unique_df <- corr_fdr01_MH_site_sig %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  split(., f = .$corrlink) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::summarise(num_sites = length(unique(.[["grouping"]]))) %>%
        droplevels
  ) %>%
  bind_rows(., .id = "corrlink") %>%
  dplyr::arrange(desc(num_sites))
#if the num_sites > 1, then the metabolite-ASV was significantly correlated in 2 or more sites

#so which ASVs are observed with multiple correlations to metabolites?
corr_fdr01_MH_site_sig_freq_unique_df <- corr_fdr01_MH_site_sig %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  dplyr::filter(asv_id %in% corr_fdr01_MH_site_freq_sig_df[["asv_id"]]) %>%
  # dplyr::left_join(., corr_fdr01_site_freq_sig_df, by = join_by(asv_id), multiple = "all", relationship = "many-to-many") %>% 
  dplyr::distinct(corrlink, asv_id, simpleName, grouping, .keep_all = TRUE) %>%
  droplevels 

#list the cosmopolitan 
corr_fdr01_MH_site_sig_freq_unique_df %>%
  count(corrlink,name = "num_occ") %>%
  dplyr::arrange(desc(num_occ)) %>%
  dplyr::filter(num_occ > 1)
# # A tibble: 16 × 2
# corrlink                          num_occ
# <fct>                               <int>
#   1 ASV_00026:pantothenic acid              3
# 2 ASV_00018:aspartate                     2
# 3 ASV_00018:glutamic acid                 2
# 4 ASV_00018:glycine                       2
# 5 ASV_00018:lysine 2                      2
# 6 ASV_00018:ornithine 2                   2
# 7 ASV_00018:threonine                     2
# 8 ASV_00018:valine                        2
# 9 ASV_00024:pantothenic acid              2
# 10 ASV_00028:sn-glycerol 3-phosphate       2
# 11 ASV_00054:sn-glycerol 3-phosphate       2
# 12 ASV_00109:cysteate                      2
# 13 ASV_00134:lysine 2                      2
# 14 ASV_00218:isethionate                   2
# 15 ASV_00360:alanine                       2
# 16 ASV_00360:threonine                     2



#so filter for only those ubiquitous ASVs (10+ correlations) correlations:
corr_fdr01_MH_site_sig_freq_unique_venn_list <- corr_fdr01_MH_site_sig_freq_unique_df %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::distinct(corrlink, grouping) %>%
        dplyr::select(corrlink) %>%
        tibble::deframe(.))



corr_fdr01_MH_site_strong_asv_df <- corr_fdr01_MH_site_sig %>%
  dplyr::ungroup(.) %>%
  dplyr::rowwise(.) %>%
  dplyr::filter(abs(filtered_estimate) >= 0.9) %>%
  droplevels

corr_fdr01_MH_site_sig_asv_venn_list <- corr_fdr01_MH_site_sig %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        count(asv_id,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
        arrange(desc(num_correlations))%>%
        droplevels %>%
        dplyr::select(asv_id) %>%
        tibble::deframe(.)
  )

corr_fdr01_MH_site_sig_strong_venn_list <- corr_fdr01_MH_site_sig %>%
  dplyr::filter(asv_id %in% corr_fdr01_MH_site_strong_asv_df[["asv_id"]]) %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::filter(abs(filtered_estimate) >= 0.9) %>%
        # count(asv_id,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
        # arrange(desc(num_correlations))%>%
        # filter(num_correlations >= 10)%>%
        droplevels %>%
        # dplyr::distinct(asv_id) %>%
        dplyr::distinct(corrlink) %>% #how many times
        tibble::deframe(.)
  )
# corr_fdr01_MH_site_venn_list <- corr_fdr01_MH_site_sig %>%
#   # dplyr::filter(asv_id %in% corr_fdr01_MH_site_strong_asv_df[["asv_id"]]) %>%
#   dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
#   dplyr::mutate(corrlink = factor(corrlink)) %>%
#   split(., f = .$grouping) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         dplyr::filter(abs(filtered_estimate) >= 0.5) %>%
#         # count(asv_id,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
#         # arrange(desc(num_correlations))%>%
#         # filter(num_correlations >= 10)%>%
#         droplevels %>%
#         # dplyr::distinct(asv_id) %>%
#         dplyr::distinct(corrlink) %>% #how many times
#         tibble::deframe(.)
#   )

corr_fdr01_MH_site_venn_list <- corr_fdr01_MH_site_sig %>%
  dplyr::mutate(pairing = paste(asv_id, simpleName, sep = ":")) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        # droplevels %>%
        split(., f = .$sig) %>%
        map(., ~.x %>%
              dplyr::distinct(pairing) %>%
              # dplyr::distinct(asv_id, simpleName) %>%
              # dplyr::distinct(asv_id) %>%
              # Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
              tibble::deframe(.) %>%
              unlist) %>%
        purrr::list_flatten(.)
  ) %>%
  purrr::list_flatten(.)
corr_fdr01_MH_site_venn_list <- corr_fdr01_MH_site_venn_list %>%
  setNames(., (names(corr_fdr01_MH_site_venn_list) %>%
                 stringr::str_replace_all(., site_lookup) %>%
                 stringr::str_replace_all(., padj_cutoff_labels) %>%
                 stringr::str_replace_all(., "_Q0", " FDR ") %>%
                 stringr::str_replace_all(., "_BH", " BH 5") %>%
                 stringr::str_c(., "%")))

temp_g1 <- UpSetR::upset(UpSetR::fromList(corr_fdr01_MH_site_venn_list),
                         # main.bar.color = "black",
                         mainbar.y.label = "Number of strong, significant\nASV-metabolite correlations",
                         set.metadata = NULL,
                         sets.x.label = "Set size",
                         sets = c(names(corr_fdr01_MH_site_venn_list)),
                         empty.intersections = NULL,
                         order.by = "degree", decreasing = F,
                         main.bar.color = viridis::turbo(n = 9),
                         scale.intersections = "identity",
                         keep.order = T,
                         query.legend = "bottom")
print(temp_g1)

#now for the rCLR transformed + Spearman correlations

#note that the adjusted p-value for FDR q = 5% is now smaller than the BH-adjusted p-value from the Spearman test.
#so filtering for (!is.na(padj_bh_05) & !is.na(sig_q01)) is much more stringent

#union(maybe, sig_q01, sig_q025, sig_05)
# [1] 7496
corr_fdr01_rCLR_site_sig <- spearman.test.site.rclr.filtered.df %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig), ~factor(.x))) %>%
  # filter(sig %in% c("sig_q01","vsig")) %>% # [1] 1959 rows
  # filter(sig %in% c("sig_q01","sig_q05")) %>% # [1] 7496 rows
  # filter(sig %in% c("sig_q01","sig_q025")) %>% # [1] 7066 rows
  filter(sig %in% c("sig_q01","vsig", "sig_q025")) %>% #7066
  droplevels
nrow(corr_fdr01_rCLR_site_sig)

corr_fdr01_rCLR_site_sig_strong_venn_list <- corr_fdr01_rCLR_site_sig %>%
  dplyr::filter(asv_id %in% corr_fdr01_rCLR_site_strong_asv_df[["asv_id"]]) %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::filter(abs(filtered_estimate) >= 0.9) %>%
        # count(asv_id,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
        # arrange(desc(num_correlations))%>%
        # filter(num_correlations >= 10)%>%
        droplevels %>%
        # dplyr::distinct(asv_id) %>%
        dplyr::distinct(corrlink) %>% #how many times
        tibble::deframe(.)
  )

#upset plot:
compare_spearman.test.site.df <- corr_fdr01_MH_site_sig %>% #this df has only "sig_q01" and "vsig" correlations. "vsig" means that the adjusted p-value is below 0.05, "sig_q01" means that p-value is below the 1% FDR threshold for p-values
  dplyr::select(asv_id, simpleName, grouping, site, sig) %>%
  dplyr::mutate(model = "MH") %>%
  dplyr::bind_rows(., corr_fdr01_rCLR_site_sig %>% #this df has sig_q01, sig_q025, sig_05, and maybe. be more stringent; retain only "sig_q01"
                     filter(sig %in% c("sig_q01","vsig", "sig_q025")) %>%
                     # dplyr::left_join(., spearman.test.site.rclr.filtered.df %>% 
                     dplyr::select(asv_id, simpleName, grouping, site, sig) %>%
                     dplyr::mutate(model = "rCLR") %>%
                     droplevels) %>%
  droplevels

# usvi_SDE_MH_rCLR_venn.list <- compare_spearman.test.site.df %>%
#   dplyr::mutate(pairing = paste(asv_id, simpleName, sep = ":")) %>%
#   split(., f = .$grouping) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         split(., f = .$model) %>%
#         map(., ~.x %>%
#               droplevels %>%
#               dplyr::ungroup(.) %>%
#               dplyr::distinct(pairing) %>%
#               # dplyr::distinct(asv_id, simpleName) %>%
#               # dplyr::distinct(asv_id) %>%
#               Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
#               unlist) %>%
#         purrr::list_flatten(.)
#   )
# View(usvi_SDE_MH_rCLR_venn.list)
# intersect(usvi_SDE_MH_rCLR_venn.list[[1]][["MH"]], usvi_SDE_MH_rCLR_venn.list[[1]][["rCLR"]])

# labels_vennlist <- map(usvi_SDE_MH_rCLR_venn.list,
#                        ~.x %>%
#                          # purrr::list_flatten(., name_spec = "{inner}_{outer}") %>%
#                          purrr::list_flatten(., name_spec = "{outer}_{inner}"))
# params_list <- list("MH" = grep(paste0(c("MH"), collapse = "|"), names(labels_vennlist[[1]]), value = TRUE, ignore.case = FALSE),
#                     "rCLR" = grep(paste0(c("rCLR"), collapse = "|"), names(labels_vennlist[[1]]), value = TRUE, ignore.case = FALSE),
#                     "MH \u2229 rCLR" = grep(paste0(c("MH", "rCLR"), collapse = "|"), names(labels_vennlist[[1]]), value = TRUE, ignore.case = FALSE))
# 
# temp_venn.list <- usvi_SDE_MH_rCLR_venn.list[[1]]
# 
# temp_g1 <- UpSetR::upset(UpSetR::fromList(temp_venn.list),
#                          main.bar.color = "black",
#                          mainbar.y.label = "Intersection size",
#                          set.metadata = NULL,
#                          sets.x.label = "Set size",
#                          sets = c(names(temp_venn.list)),
#                          empty.intersections = NULL,
#                          order.by = "degree", decreasing = F,
#                          queries = list(list(query = intersects, params = params_list[1], color = "red", active = T, query.name = paste(names(params_list[1]))),
#                                         list(query = intersects, params = params_list[2], color = "blue", active = T, query.name = paste(names(params_list[2]))),
#                                         list(query = intersects, params = params_list[3], color = "darkgreen", active = T, query.name = paste(names(params_list[3])))),
#                          keep.order = T,
#                          query.legend = "bottom")
# print(temp_g1)

usvi_SDE_MH_rCLR_venn.list <- compare_spearman.test.site.df %>%
  dplyr::mutate(pairing = paste(asv_id, simpleName, sep = ":")) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        split(., f = .$model) %>%
        map(., ~.x %>%
              droplevels %>%
              dplyr::ungroup(.) %>%
              dplyr::distinct(pairing) %>%
              # dplyr::distinct(asv_id, simpleName) %>%
              # dplyr::distinct(asv_id) %>%
              Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
              unlist) %>%
        purrr::list_flatten(.)
  )

usvi_SDE_MH_rCLR_venn.list <- compare_spearman.test.site.df %>%
  dplyr::mutate(pairing = paste(asv_id, simpleName, sep = ":")) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        # droplevels %>%
        split(., f = .$model) %>%
        map(., ~.x %>%
              dplyr::distinct(pairing) %>%
              # dplyr::distinct(asv_id, simpleName) %>%
              # dplyr::distinct(asv_id) %>%
              # Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
              tibble::deframe(.) %>%
              unlist) %>%
        purrr::list_flatten(.)
  ) %>%
  purrr::list_flatten(.)

usvi_SDE_MH_rCLR_venn.list <- usvi_SDE_MH_rCLR_venn.list %>%
  setNames(., (names(usvi_SDE_MH_rCLR_venn.list) %>%
                 stringr::str_replace_all(., site_lookup) %>%
                 stringr::str_replace_all(., "_", " ") %>%
                 stringr::str_replace_all(., "MH", "Morisita-Horn")))


summary(usvi_SDE_MH_rCLR_venn.list)
# Length Class  Mode     
# Seagrass Morisita-Horn      980   -none- character
# Seagrass rCLR               822   -none- character
# Tektite Reef Morisita-Horn  442   -none- character
# Tektite Reef rCLR          3122   -none- character
# Yawzi Reef Morisita-Horn    373   -none- character
# Yawzi Reef rCLR            3122   -none- character


upset_length <- data.frame("pair1" = c(names(usvi_SDE_MH_rCLR_venn.list)),
                           "pair2" = NA) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(n_intersect = length(usvi_SDE_MH_rCLR_venn.list[[ {pair1} ]])) %>%
  dplyr::mutate( pair1 = factor(pair1)) %>%
  dplyr::arrange(pair1) %>%
  dplyr::bind_rows(., (names(usvi_SDE_MH_rCLR_venn.list) %>%
                         combn(., 2) %>%
                         t() %>% as.data.frame(.) %>%
                         setNames(., c("pair1", "pair2")) %>%
                         droplevels %>%
                         dplyr::rowwise(.) %>%
                         dplyr::mutate(n_intersect = length(intersect(usvi_SDE_MH_rCLR_venn.list[[{pair1}]], usvi_SDE_MH_rCLR_venn.list[[{pair2}]]))) %>%
                         # dplyr::filter(n_intersect > 0) %>%
                         dplyr::mutate(site = pair1) %>%
                         dplyr::mutate(site = stringr::str_remove_all(pair1, "( )(rCLR|M.*H.*$)")) %>%
                         dplyr::filter(grepl(site, pair2)) %>%
                         droplevels)) %>%
  dplyr::arrange(pair1) %>%
  dplyr::mutate(label = dplyr::case_when(is.na(pair2) ~ pair1,
                                         !is.na(pair2) ~ paste0(site, " MH \u2229 rCLR"),
                                         .default = NA)) %>%
  droplevels

params_list <- upset_length %>%
  dplyr::select(label, pair1, pair2) %>%
  split(., f = .$label) %>%
  purrr::map(., ~.x %>%
               dplyr::select(pair1, pair2) %>%
               simplify %>%
               na.omit(.) %>%
               unlist %>%
               as.character
  ) %>%
  purrr::list_flatten(.)
               
# temp_venn2.list <- usvi_SDE_MH_rCLR_venn.list
# rlang::parse_expr(deparse(as.list(params_list[[1]])))
# # list("Seagrass Morisita-Horn", "Seagrass rCLR")
# rlang::parse_expr(deparse(as.list(params_list[[7]])))
# # list("Yawzi Reef Morisita-Horn", "Yawzi Reef rCLR")
# # temp_idx <- deparse(as.list(params_list[[1]])) 
# temp_idx <- c(deparse(as.list(params_list[[1]])), deparse(as.list(params_list[[7]])))

temp_g2 <- UpSetR::upset(UpSetR::fromList(usvi_SDE_MH_rCLR_venn.list),
                         # main.bar.color = "black",
                         mainbar.y.label = "Intersection of strong, significant\nASV-metabolite correlations",
                         # set.metadata = NULL,
                         sets.x.label = "Set size",
                         sets = c(names(usvi_SDE_MH_rCLR_venn.list)),
                         empty.intersections = NULL,
                         order.by = "degree", decreasing = F,
                         # intersections = list(list("Seagrass Morisita-Horn"),
                         #                      list("Yawzi Reef rCLR")),
                         # intersections = list(list(eval(params_list[[1]])),
                         #                      list(eval(params_list[[4]]))),
                         queries = list(list(query = intersects, params = params_list[1], active = T, query.name = paste(names(params_list[1]))),
                                        list(query = intersects, params = params_list[2], active = T, query.name = paste(names(params_list[2]))),
                                        list(query = intersects, params = params_list[3], active = T, query.name = paste(names(params_list[3]))),
                                        list(query = intersects, params = params_list[4],active = T, query.name = paste(names(params_list[4]))),
                                        list(query = intersects, params = params_list[5], active = T, query.name = paste(names(params_list[5]))),
                                        list(query = intersects, params = params_list[6],  active = T, query.name = paste(names(params_list[6]))),
                                        # list(query = intersects, params = params_list[7], active = T, query.name = paste(names(params_list[7]))),
                                        list(query = intersects, params = params_list[8],  active = T, query.name = paste(names(params_list[8]))),
                                        list(query = intersects, params = params_list[9],  active = T, query.name = paste(names(params_list[9])))),
                         # main.bar.color = viridis::turbo(n = nrow(upset_length)),
                         # main.bar.color = viridis::turbo(n = 23),
                         # scale.intersections = "identity",
                         keep.order = T,   set_size.show = FALSE, 
                         # group.by = "sets",
                         query.legend = "bottom")
print(temp_g2)


#how many intersections between rCLR correlations?
usvi_SDE_asv_rCLR_venn.list <- corr_fdr01_rCLR_site_sig %>%
  dplyr::filter(abs(estimate) >= 0.5) %>%  
  filter(sig %in% c("sig_q01","sig_q025", "sig_q05")) %>% # [1] 1571
  dplyr::mutate(pairing = paste(asv_id, simpleName, sep = ":")) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        droplevels %>%
        split(., f = .$sig) %>%
        map(., ~.x %>%
              dplyr::distinct(pairing) %>%
              # dplyr::distinct(asv_id, simpleName) %>%
              # dplyr::distinct(asv_id) %>%
              # Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
              tibble::deframe(.) %>%
              unlist) %>%
        purrr::list_flatten(.)
  ) %>%
  purrr::list_flatten(.)


# [1] "Seagrass FDR 5%"       "Seagrass FDR 2.5%"     "Seagrass FDR 1%"      
# [4] "Tektite Reef FDR 5%"   "Tektite Reef FDR 2.5%" "Tektite Reef FDR 1%"  
# [7] "Yawzi Reef FDR 5%"     "Yawzi Reef FDR 2.5%"   "Yawzi Reef FDR 1%" 
usvi_SDE_asv_rCLR_venn.list <- usvi_SDE_asv_rCLR_venn.list %>%
  setNames(., (names(usvi_SDE_asv_rCLR_venn.list) %>%
                 stringr::str_replace_all(., site_lookup) %>%
                 stringr::str_replace_all(., padj_cutoff_labels) %>%
                 stringr::str_replace_all(., "_Q0", " FDR ") %>%
                 stringr::str_c(., "%")))

summary(usvi_SDE_asv_rCLR_venn.list)
# Length Class  Mode     
# Seagrass FDR 5%       2938   -none- character
# Seagrass FDR 2.5%      822   -none- character
# Tektite Reef FDR 5%   1926   -none- character
# Tektite Reef FDR 2.5% 2229   -none- character
# Tektite Reef FDR 1%    893   -none- character
# Yawzi Reef FDR 5%      673   -none- character
# Yawzi Reef FDR 2.5%   2056   -none- character
# Yawzi Reef FDR 1%     1066   -none- character



upset_length <- data.frame("pair1" = c(names(usvi_SDE_asv_rCLR_venn.list)),
                           "pair2" = NA) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(n_intersect = length(usvi_SDE_asv_rCLR_venn.list[[ {pair1} ]])) %>%
  dplyr::mutate( pair1 = factor(pair1)) %>%
  dplyr::arrange(pair1) %>%
  dplyr::bind_rows(., (names(usvi_SDE_asv_rCLR_venn.list) %>%
                         combn(., 2) %>%
                         t() %>% as.data.frame(.) %>%
                         setNames(., c("pair1", "pair2")) %>%
                         droplevels %>%
                         dplyr::rowwise(.) %>%
                         dplyr::mutate(n_intersect = length(intersect(usvi_SDE_asv_rCLR_venn.list[[{pair1}]], usvi_SDE_asv_rCLR_venn.list[[{pair2}]]))) %>%
                         dplyr::filter(n_intersect > 0) %>%
                         droplevels)) %>%
  dplyr::arrange(pair1) %>%
  dplyr::filter(n_intersect > 0) %>%
  droplevels


temp_g3 <- UpSetR::upset(UpSetR::fromList(usvi_SDE_asv_rCLR_venn.list),
                         # main.bar.color = "black",
                         mainbar.y.label = "Number of strong, significant\nASV-metabolite correlations",
                         set.metadata = NULL,
                         sets.x.label = "Set size",
                         sets = c(names(usvi_SDE_asv_rCLR_venn.list)),
                         empty.intersections = NULL,
                         order.by = "degree", decreasing = F,
                         main.bar.color = viridis::turbo(n = nrow(upset_length)),
                         scale.intersections = "identity",
                         keep.order = T,
                         query.legend = "bottom")
print(temp_g3)

# #Sankey graph of what happened to the correlations:
# temp_df <- compare_spearman.test.site.df %>%
#   dplyr::slice_sample(n = 100, by = "grouping") %>%
#   dplyr::mutate(recoded_sig = dplyr::case_when(
#     (sig_1 == "n.s.") & (sig_2 != "n.s.") ~ -1,
#     (sig_1 != "n.s.") & (sig_2 != "n.s.") ~ 0,
#     (sig_1 != "n.s.") & (sig_2 == "n.s.") ~ 1,
#     .default = NA)) 
# temp_df %>%
#   dplyr::count(recoded_sig, grouping, name = "obs") %>%
#   dplyr::mutate(label = dplyr::case_when(
#     (recoded_sig == -1) ~ "MH",
#     (recoded_sig ==  0) ~ "both",
#     (recoded_sig == 1) ~ "rCLR",
#     .default = NA)) %>%
#   dplyr::select(recoded_sig, grouping, label) %>%
#   droplevels
# temp_sankey <- temp_df %>%
#   # dplyr::arrange(taxonomy_string) %>%
#   split(., f = .$grouping) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         dplyr::mutate(recoded_sig = dplyr::case_when(
#           (sig_1 == "n.s.") & (sig_2 != "n.s.") ~ -1,
#           (sig_1 != "n.s.") & (sig_2 != "n.s.") ~ 0,
#           (sig_1 != "n.s.") & (sig_2 == "n.s.") ~ 1,
#           .default = NA)) %>%
#         ggsankey::make_long(., model_1, model_2, value = "recoded_sig") %>%
#         # ggsankey::make_long(., asv_id, simpleName, value = "recoded_sig") %>%
#         # ggsankey::make_long(., asv_id, simpleName, sig_1, sig_2) %>%
#         # ggsankey::make_long(., c(asv_id, simpleName)) %>%
#         # ggsankey::make_long(., c(sig_1, sig_2)) %>%
#         # dplyr::left_join(., (full_coassembly_taxonomy_skani_best.df %>%
#         #                        dplyr::select(-c(contains("taxonomy"))) %>%
#         #                        tidyr::pivot_longer(., cols = c(all_of(keep_tax), contains("taxonomy")),
#         #                                            names_to = NULL, values_to = "value") %>%
#         #                        tidyr::drop_na(value) %>%
#         #                        dplyr::mutate(value = factor(value, levels = unique(.[["value"]]))) %>%
#         #                        dplyr::count(value, name = "obs") %>%
#         #                        dplyr::filter(!grepl("Unclassified", value)) %>%
#         #                        dplyr::mutate(label = dplyr::case_when((obs < 5) ~ NA,
#         #                                                               (obs >= 5) ~ value)) %>%
#         #                        dplyr::select(value, label) %>%
#         #                        droplevels),
#         #                  by = join_by("node" == "value"), multiple = "first", relationship = "many-to-many") %>%
#         dplyr::rename(recoded_sig = "value") %>%
#         droplevels)
# temp_sankey <- temp_sankey %>%
#   imap(., ~.x %>%
#         dplyr::left_join(., temp_df %>%
#                            dplyr::filter(grouping == .y) %>%
#                            droplevels %>%
#                            dplyr::count(recoded_sig, grouping, name = "obs") %>%
#                            dplyr::mutate(label = dplyr::case_when(
#                              (recoded_sig == -1) ~ "MH",
#                              (recoded_sig ==  0) ~ "both",
#                              (recoded_sig == 1) ~ "rCLR",
#                              .default = NA)) %>%
#                            dplyr::select(recoded_sig, label) %>%
#                            droplevels,
#                          by = join_by("recoded_sig"), multiple = "first", relationship = "many-to-many") %>%
#         # dplyr::mutate(label = dplyr::case_when(num_taxa > 1 ~ label,
#         #                                        .default = NA)) %>%
#         # dplyr::arrange(desc(num_taxa)) %>%
#         dplyr::mutate(across(c(node, next_node), ~factor(.x, ordered = FALSE, levels = unique(.x)))))
# 
# g8 <- (
#   # ggplot(data = temp_df[[1]],
#   ggplot(data = temp_sankey[[1]],
#          aes(x = x, next_x = next_x,
#              node = node, next_node = next_node,        
#              label = label,
#              color = factor(node),
#              fill = factor(node)))
#   + ggsankey::geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)
#   + ggsankey::geom_sankey_label(size = 3, color = "white", fill = "gray40")
#   + ggsankey::theme_sankey(base_size = 18)
#   # + ggsankey::geom_alluvial(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)
#   # + ggsankey::geom_alluvial_label(size = 3, color = "white", fill = "gray40")
#   # + ggsankey::theme_alluvial(base_size = 18)
#   + theme_dark()
#   + theme(
#     axis.title = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank())
#   # + scale_fill_manual(values = annotation_taxa_colors,
#   #                     labels = names(annotation_taxa_colors),
#   #                     breaks = names(annotation_taxa_colors),
#   #                     drop = TRUE)
#   # + scale_color_manual(values = annotation_taxa_colors,
#   #                      labels = names(annotation_taxa_colors),
#   #                      breaks = names(annotation_taxa_colors),
#   #                      drop = TRUE)
# )

{
  #ubiquitous ASVs: the unique group of ASVs with 10 or more correlations with metabolites
  #if we see that ASV-metabolite correlation in multiple sites:
  corr_fdr01_rCLR_site_sig %>%
    distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
    count(asv_id, name = "num_correlations") %>%   # Count occurrences
    arrange(desc(num_correlations))%>% 
    distinct()%>%
    # dplyr::filter(num_correlations >= 2) %>% #1536 with the filtering step
    nrow(.)
  # [1] 4073
  
  # corr_fdr01_rCLR_site_sig %>%
  #   count(asv_id, name = "num_correlations") %>%   # Count occurrences
  #   arrange(desc(num_correlations))%>% 
  #   distinct()%>%
  #   # dplyr::filter(num_correlations >= 10) %>%
  #   nrow(.)
  # # [1] 24
  # 
  # # #there are 3 ASVs that are correlated with 10 or more metabolites, in 2 or more sites:
  # corr_fdr01_rCLR_site_freq_sig_df <- dplyr::full_join((corr_fdr01_rCLR_site_sig %>%
  #                                                       distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
  #                                                       count(asv_id, name = "num_unique_correlations") %>%   # Count occurrences
  #                                                       arrange(desc(num_unique_correlations))),
  #                                                    (corr_fdr01_rCLR_site_sig %>%
  #                                                       count(asv_id, name = "num_correlations") %>%   # Count occurrences
  #                                                       arrange(desc(num_correlations)))) %>%
  #   dplyr::filter(num_correlations >= 10) %>%
  #   droplevels
  # corr_fdr01_rCLR_site_freq_sig_df %>%
  #   dplyr::filter(num_unique_correlations < 10 ) %>%
  #   droplevels %>%
  #   dplyr::distinct(asv_id) %>%
  #   tibble::deframe(.)
  # # # [1] ASV_00024 ASV_00026 ASV_00028
  # # 
  # corr_fdr01_rCLR_site_sig %>%
  #   distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
  #   count(simpleName,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
  #   arrange(desc(num_correlations))%>%
  #   head()
  # # # A tibble: 6 × 2
  # # simpleName         num_correlations
  # # <chr>                         <int>
  # #   1 pantothenic acid                216
  # # 2 ciliatine                       141
  # # 3 cysteate                        100
  # # 4 homoserine betaine               96
  # # 5 kynurenine                       81
  # # 6 GABA                             80
  # 
  # 
  # #This is our overall count of significant and strong (abs(rho) >= 0.5) ASV-Metabolite correlations:
  # corr_fdr01_rCLR_site_sig_unique_venn_list <- corr_fdr01_rCLR_site_sig %>%
  #   dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  #   dplyr::mutate(corrlink = factor(corrlink)) %>%
  #   split(., f = .$grouping) %>%
  #   map(., ~.x %>%
  #         droplevels %>%
  #         dplyr::distinct(corrlink) %>%
  #         tibble::deframe(.))
  # 
  # 
  # #now to count how many UNIQUE ASV-metab associations there are across all sites
  # #first, see how many metabolite-ASV correlations are observed in 2 or more sites ("cosmopolitan")
  # corr_fdr01_rCLR_site_sig_unique_df <- corr_fdr01_rCLR_site_sig %>%
  #   dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  #   dplyr::mutate(corrlink = factor(corrlink)) %>%
  #   split(., f = .$corrlink) %>%
  #   map(., ~.x %>%
  #         droplevels %>%
  #         dplyr::summarise(num_sites = length(unique(.[["grouping"]]))) %>%
  #         droplevels
  #   ) %>%
  #   bind_rows(., .id = "corrlink") %>%
  #   dplyr::arrange(desc(num_sites))
  # #if the num_sites > 1, then the metabolite-ASV was significantly correlated in 2 or more sites
  # 
  # #so which ASVs are observed with multiple correlations to metabolites?
  # corr_fdr01_rCLR_site_sig_freq_unique_df <- corr_fdr01_rCLR_site_sig %>%
  #   dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  #   dplyr::mutate(corrlink = factor(corrlink)) %>%
  #   dplyr::filter(asv_id %in% corr_fdr01_rCLR_site_freq_sig_df[["asv_id"]]) %>%
  #   # dplyr::left_join(., corr_fdr01_site_freq_sig_df, by = join_by(asv_id), multiple = "all", relationship = "many-to-many") %>% 
  #   dplyr::distinct(corrlink, asv_id, simpleName, grouping, .keep_all = TRUE) %>%
  #   droplevels 
  # 
  # #list the cosmopolitan 
  # corr_fdr01_rCLR_site_sig_freq_unique_df %>%
  #   count(corrlink,name = "num_occ") %>%
  #   dplyr::arrange(desc(num_occ)) %>%
  #   dplyr::filter(num_occ > 1)
  # # # A tibble: 16 × 2
  # # corrlink                          num_occ
  # # <fct>                               <int>
  # #   1 ASV_00026:pantothenic acid              3
  # # 2 ASV_00018:aspartate                     2
  # # 3 ASV_00018:glutamic acid                 2
  # # 4 ASV_00018:glycine                       2
  # # 5 ASV_00018:lysine 2                      2
  # # 6 ASV_00018:ornithine 2                   2
  # # 7 ASV_00018:threonine                     2
  # # 8 ASV_00018:valine                        2
  # # 9 ASV_00024:pantothenic acid              2
  # # 10 ASV_00028:sn-glycerol 3-phosphate       2
  # # 11 ASV_00054:sn-glycerol 3-phosphate       2
  # # 12 ASV_00109:cysteate                      2
  # # 13 ASV_00134:lysine 2                      2
  # # 14 ASV_00218:isethionate                   2
  # # 15 ASV_00360:alanine                       2
  # # 16 ASV_00360:threonine                     2
  # 
  # 
  # 
  # #so filter for only those ubiquitous ASVs (10+ correlations) correlations:
  # corr_fdr01_rCLR_site_sig_freq_unique_venn_list <- corr_fdr01_rCLR_site_sig_freq_unique_df %>%
  #   split(., f = .$grouping) %>%
  #   map(., ~.x %>%
  #         droplevels %>%
  #         dplyr::distinct(corrlink, grouping) %>%
  #         dplyr::select(corrlink) %>%
  #         tibble::deframe(.))
  # 
  # 
  # 
  # temp_g3 <- ggVennDiagram::ggVennDiagram(corr_fdr01_rCLR_site_sig_unique_venn_list, 
  #                                         set_size = NA,
  #                                         label = "count",
  #                                         label_color = "white",
  #                                         label_geom = "label") +
  #   scale_fill_viridis(discrete = FALSE, option = "turbo",
  #                      # limits = c(0, max(myBreaks)),
  #                      guide = "none",
  #                      name = "ASV count") +
  #   scale_x_continuous(expand = expansion(mult = 0.2)) +
  #   ggtitle("Total") +
  #   theme_void()
  # 
  # temp_g4 <- ggVennDiagram::ggVennDiagram(corr_fdr01_rCLR_site_sig_freq_unique_venn_list, 
  #                                         set_size = NA,
  #                                         label = "count",
  #                                         label_color = "white",
  #                                         label_geom = "label") +
  #   scale_fill_viridis(discrete = FALSE, option = "turbo",
  #                      # limits = c(0, max(myBreaks)),
  #                      guide = "none",
  #                      name = "ASV count") +
  #   scale_x_continuous(expand = expansion(mult = 0.2)) +
  #   ggtitle("Occurences of correlations between \nmetabolites and ASVs with multiple \n(10+) correlations") +
  #   theme_void()
  # 
  # 
  # gvenn2 <- temp_g3 + temp_g4 + patchwork::plot_layout(guides = "collect") + 
  #   patchwork::plot_annotation(title = "How many ASVs are strongly (>0.5) correlated with metabolites in each site?", tag_levels = "A")
  # gvenn2
  # 
  corr_fdr01_rCLR_site_strong_asv_df <- corr_fdr01_rCLR_site_sig %>%
    dplyr::ungroup(.) %>%
    dplyr::rowwise(.) %>%
    # dplyr::filter(abs(filtered_estimate) >= 0.9) %>% #only 4 observations
    # dplyr::filter(abs(filtered_estimate) >= 0.8) %>% #249 observations
    dplyr::filter(abs(filtered_estimate) >= 0.85) %>% #66 observations
    droplevels
  # 
  # corr_fdr01_rCLR_site_sig_asv_venn_list <- corr_fdr01_rCLR_site_sig %>%
  #   split(., f = .$grouping) %>%
  #   map(., ~.x %>%
  #         droplevels %>%
  #         count(asv_id,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
  #         arrange(desc(num_correlations))%>%
  #         droplevels %>%
  #         dplyr::select(asv_id) %>%
  #         tibble::deframe(.)
  #   )
  # 
  # corr_fdr01_rCLR_site_sig_strong_venn_list <- corr_fdr01_rCLR_site_sig %>%
  #   dplyr::filter(asv_id %in% corr_fdr01_rCLR_site_strong_asv_df[["asv_id"]]) %>%
  #   dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  #   dplyr::mutate(corrlink = factor(corrlink)) %>%
  #   split(., f = .$grouping) %>%
  #   map(., ~.x %>%
  #         droplevels %>%
  #         dplyr::filter(abs(filtered_estimate) >= 0.9) %>%
  #         # count(asv_id,name = "num_correlations") %>%  # Count occurrences of each metabolite-ASV pair
  #         # arrange(desc(num_correlations))%>%
  #         # filter(num_correlations >= 10)%>%
  #         droplevels %>%
  #         # dplyr::distinct(asv_id) %>%
  #         dplyr::distinct(corrlink) %>% #how many times
  #         tibble::deframe(.)
  #   )
  # 
  # 
  # temp_g1 <- ggVennDiagram::ggVennDiagram(corr_fdr01_rCLR_site_sig_asv_venn_list, 
  #                                         set_size = NA,
  #                                         label = "count",
  #                                         label_color = "white",
  #                                         label_geom = "label") +
  #   scale_fill_viridis(discrete = FALSE, option = "turbo",
  #                      # limits = c(0, max(myBreaks)),
  #                      guide = "none",
  #                      name = "ASV count") +
  #   scale_x_continuous(expand = expansion(mult = 0.2)) +
  #   ggtitle("Number of ASVs with significant correlations to metabolites") +
  #   theme_void()
  # 
  # temp_g2 <- ggVennDiagram::ggVennDiagram(corr_fdr01_rCLR_site_sig_strong_venn_list,
  #                                         set_size = NA,
  #                                         label = "count",
  #                                         label_color = "white",
  #                                         label_geom = "label") +
  #   scale_fill_viridis(discrete = FALSE, option = "turbo",
  #                      # limits = c(0, max(myBreaks)),
  #                      guide = "none",
  #                      name = "ASV count") +
  #   scale_x_continuous(expand = expansion(mult = 0.2)) +
  #   ggtitle("ASVs with very strong correlations to metabolites") +
  #   theme_void()
  # 
  # gvenn1 <- temp_g1 + temp_g2 + patchwork::plot_layout(guides = "collect") + 
  #   patchwork::plot_annotation(title = "How many ASVs are significantly correlated with metabolites in each site?", tag_levels = "A")
  # gvenn1
  
}

corr_fdr01_rCLR_site_sig_unique_df <- corr_fdr01_rCLR_site_sig %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  split(., f = .$corrlink) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::summarise(num_sites = length(unique(.[["grouping"]]))) %>%
        droplevels
  ) %>%
  bind_rows(., .id = "corrlink") %>%
  dplyr::arrange(desc(num_sites))
# # A tibble: 1,126 × 2
# corrlink                   num_occ
# <fct>                        <int>
#   1 ASV_00130:pantothenic acid       2
# 2 ASV_00138:arginine               2

#no ASV-metabolite correlations are observed across all 3 sites
#1126 ASV-metabolite correlations are significant and strong at 2 of 3 sites

corr_fdr01_rCLR_site_freq_sig_df <- dplyr::full_join((corr_fdr01_rCLR_site_sig %>%
                                                        dplyr::ungroup(.) %>%
                                                              distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
                                                              count(asv_id, name = "num_unique_correlations") %>%   # Count occurrences
                                                              arrange(desc(num_unique_correlations))),
                                                           (corr_fdr01_rCLR_site_sig %>%
                                                              dplyr::ungroup(.) %>%
                                                              count(asv_id, name = "num_correlations") %>%   # Count occurrences
                                                              arrange(desc(num_correlations)))) %>%
  # dplyr::filter(num_correlations >= 5) %>% #was 10 for the MH-transformed ASV-metabolite
  # dplyr::filter(num_correlations >= 10) %>% #was 10 for the MH-transformed ASV-metabolite
  dplyr::filter(num_correlations >= 2) %>% #2 is the threshold for the rCLR-transformed genus-metabolite correlations
  # dplyr::filter(num_unique_correlations >= 3) %>% #2 is the threshold for the rCLR-transformed genus-metabolite correlations
  droplevels
# corr_fdr01_rCLR_site_freq_sig_df <- dplyr::full_join((corr_fdr01_rCLR_site_sig %>%
#                                                               distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
#                                                               count(asv_id, name = "num_unique_correlations") %>%   # Count occurrences
#                                                               arrange(desc(num_unique_correlations))),
#                                                            (corr_fdr01_rCLR_site_sig %>%
#                                                               count(asv_id, name = "num_correlations") %>%   # Count occurrences
#                                                               arrange(desc(num_correlations)))) %>%
#   # dplyr::filter(num_correlations >= 10) %>% #was 10 for the MH-transformed ASV-metabolite 
#   dplyr::filter(num_correlations >= 2) %>%
#   droplevels
#1809 ASVs have 2+ correlations

#so which ASVs are observed with multiple correlations to metabolites?
corr_fdr01_rCLR_site_sig_freq_unique_df <- corr_fdr01_rCLR_site_sig %>%
  dplyr::distinct(asv_id, simpleName, grouping, site, .keep_all = TRUE) %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  dplyr::filter(asv_id %in% corr_fdr01_rCLR_site_freq_sig_df[["asv_id"]]) %>%
  # dplyr::left_join(., corr_fdr01_site_freq_sig_df, by = join_by(asv_id), multiple = "all", relationship = "many-to-many") %>% 
  dplyr::distinct(corrlink, asv_id, simpleName, grouping, .keep_all = TRUE) %>%
  droplevels 

corr_fdr01_rCLR_site_sig_freq_unique_df %>%
  count(corrlink,name = "num_occ") %>%
  dplyr::arrange(desc(num_occ)) %>%
  dplyr::filter(num_occ > 1)
# # A tibble: 273 × 2
# corrlink                   num_occ
# <fct>                        <int>
#   1 ASV_00725:pantothenic acid       2
# 2 ASV_01074:pantothenic acid       2
# 3 ASV_01210:pantothenic acid       2
#none of these correlations between ASV-metabolite are in more than 2 of 3 sites


# Repeat Spearman with taxglom --------------------------------------------

#recalculate using rCLR transformed ASV table
to_import <- c("spearman.test.genus.rclr.optA.list",
               "spearman.test.genus.site.rclr.optA.list", 
               "spearman.test.genus.site.time.rclr.optA.list")

for(file in to_import){
  if(!exists(file, envir = .GlobalEnv)){
    namevar <- file
    if(file.exists(paste0(projectpath, "/", namevar, ".rds"))){
      cli::cli_alert_info("Importing this dataset: {namevar}")
      temp_df <- readr::read_rds(paste0(projectpath, "/", namevar, ".rds"))
      assign(paste0(namevar), temp_df, envir = .GlobalEnv)
      rm(temp_df)
      to_import <- to_import[grep(namevar, to_import, value = FALSE, invert = TRUE)]
    } else {
      cli::cli_alert_warning("Please prepare this dataset: {namevar}")
    }
  }
}
if(length(to_import) > 0){
  
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
  
  keep_tax <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  
  usvi_genus_rclr.tbl <- usvi_asv.tbl %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    dplyr::select(rownames(usvi_metab.tbl)) %>%
    tibble::rownames_to_column(var = "asv_id") %>%
    otu_to_taxonomy(., usvi_prok_asvs.taxa, level = "Genus") %>% droplevels
  # otu_to_taxonomy(., usvi_prok_asvs.taxa, level = "Genus") %>%
  # vegan::decostand(., method = "rclr", 2, na.rm = TRUE, impute = TRUE) %>%
  # t()
  
  
  usvi_prok_genus.taxa <- usvi_genus_rclr.tbl %>%
    dplyr::select(all_of(keep_tax)) %>%
    dplyr::left_join(., usvi_prok_asvs.taxa %>%
                       dplyr::select(all_of(keep_tax), asv_id) %>%
                       droplevels, 
                     relationship = "one-to-many", multiple = "first") %>%
    dplyr::mutate(p_asv_id = paste0("p", asv_id)) %>%
    dplyr::select(-asv_id) %>%
    dplyr::left_join(., usvi_prok_asvs.taxa %>%
                       dplyr::select(all_of(keep_tax), asv_id) %>%
                       droplevels, 
                     relationship = "one-to-many", multiple = "all")
  
  usvi_genus_rclr.tbl <- usvi_genus_rclr.tbl %>%
    # temp_df <- usvi_genus_rclr.tbl %>%
    dplyr::left_join(., usvi_prok_genus.taxa,
                     relationship = "many-to-many", multiple = "first") %>%
    dplyr::select(p_asv_id, rownames(usvi_metab.tbl)) %>%
    tibble::column_to_rownames(var = "p_asv_id") %>%
    # tibble::column_to_rownames(var = "p_asv_id") %>%      droplevels
    vegan::decostand(., method = "rclr", 2, na.rm = TRUE, impute = TRUE) %>%
    t()
  
  # temp_df <- temp_df[rownames(usvi_metab.tbl),]
  usvi_genus_rclr.tbl <- usvi_genus_rclr.tbl[rownames(usvi_metab.tbl),]
  
  if(any(grepl("spearman.test.genus.rclr.optA.list", to_import))){
    #all samples, all times:
    {
      spearman.test <- matrix(nrow = ncol(usvi_genus_rclr.tbl), ncol = ncol(usvi_metab.tbl))
      colnames(spearman.test) <- colnames(usvi_metab.tbl)
      rownames(spearman.test) <- colnames(usvi_genus_rclr.tbl)
      
      spearman.test.rho <- spearman.test
      spearman.test.n <- spearman.test  
      
      y <- length(colnames(spearman.test))
      for(j in seq_len(y)){
        vector_metab <- usvi_metab.tbl[, j]
        for(i in seq_len(nrow(spearman.test))){
          vector_microb <- usvi_genus_rclr.tbl[,i]
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
        # dplyr::mutate(asv_id = factor(asv_id, levels = labels(dend_asv))) %>%
        dplyr::mutate(simpleName = factor(simpleName, levels = labels(dend_metab))) %>%
        dplyr::arrange(asv_id, simpleName) %>%
        dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
        # tidyr::drop_na(contains("padj")) %>%
        dplyr::mutate(label = signif(estimate, digits = 2)) %>%
        dplyr::ungroup(.) %>%
        dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
        droplevels
      
      spearman.test.genus.rclr.optA.list <- list(spearman.test.df,
                                                 spearman.test, spearman.test.rho, spearman.test.n,
                                                 # dend_asv, dend_metab,
                                                 padj_cutoff) %>%
        setNames(., c("spearman.test.df",
                      "spearman.test", "spearman.test.rho", "spearman.test.n",
                      # "dend_asv", "dend_metab", 
                      "padj_cutoff"))
    }
    readr::write_rds(spearman.test.genus.rclr.optA.list, paste0(projectpath, "/", "spearman.test.genus.rclr.optA.list", ".rds"), compress = "gz")
    
  }
  
  if(any(grepl("spearman.test.genus.site.rclr.optA.list", to_import))){
    #by site, option A:
    {
      temp_df <- usvi_metab.tbl %>%
        # tibble::as_tibble(rownames = "sample_id") %>%
        # data.table::transpose(., keep.names = "simpleName", make.names = "sample_id") %>%
        t() %>%
        as.data.frame(.) %>%
        # dplyr::select(which(colSums(.) > 0)) #remove any that are insignificant.
        tibble::as_tibble(rownames = "simpleName") %>% droplevels
      
      
      usvi_metab_site.list <- metabolomics_sample_metadata %>%
        dplyr::filter(grepl("seawater", sample_type)) %>%
        dplyr::filter(sample_id %in% rownames(usvi_genus_rclr.tbl)) %>%
        dplyr::distinct(sample_id, site) %>%
        droplevels %>%
        split(., f = .$site) %>%
        map(., ~.x %>%
              droplevels %>%
              dplyr::select(sample_id) %>%
              dplyr::left_join(., usvi_metab.tbl %>%
                                 tibble::rownames_to_column(var = "sample_id")) %>%
              # dplyr::select(which(colSums(.) > 0)) %>% #remove any that are insignificant.
              tibble::column_to_rownames(var = "sample_id") %>%
              droplevels)
      
      # temp_df <- usvi_metab_site.list %>%
      #   map(., ~.x %>%
      #         # data.table::transpose(., keep.names = "simpleName", make.names = "sample_id") %>%
      #         t() %>%
      #         as.data.frame(.) %>%
      #         # dplyr::select(which(colSums(.) > 0)) #remove any that are insignificant.
      #         # tibble::as_tibble(rownames = "simpleName") %>%
      #         droplevels
      #   )
      
      usvi_asv_site.list <- metabolomics_sample_metadata %>%
        dplyr::filter(grepl("seawater", sample_type)) %>%
        dplyr::filter(sample_id %in% rownames(usvi_genus_rclr.tbl)) %>%
        dplyr::distinct(sample_id, site) %>%
        droplevels %>%
        split(., f = .$site) %>%
        map(., ~.x %>%
              droplevels %>%
              dplyr::select(sample_id) %>%
              dplyr::inner_join(., usvi_genus_rclr.tbl %>%
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
        
        spearman.test.genus.site.list <- list(temp_spearman.test.df, 
                                              temp_spearman.test, temp_spearman.test.rho, temp_spearman.test.n,
                                              # temp_dend_asv, temp_dend_metab, 
                                              padj_cutoff) %>%
          setNames(., c("spearman.test.df", 
                        "spearman.test", "spearman.test.rho", "spearman.test.n",
                        # "dend_asv", "dend_metab", 
                        "padj_cutoff"))
        
        assign(paste0("spearman.test.genus.siteA.", namevar, ".list"), spearman.test.genus.site.list, envir = .GlobalEnv)
        rm(spearman.test.genus.site.list)
        rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
      }
      
      
      spearman.test.genus.site.rclr.optA.list <- lapply(apropos("^spearman.test.genus.siteA.*$", mode = "list"),
                                                        get) %>%
        setNames(., names(usvi_asv_site.list))
      
    }
    readr::write_rds(spearman.test.genus.site.rclr.optA.list, paste0(projectpath, "/", "spearman.test.genus.site.rclr.optA.list", ".rds"), compress = "gz")
    
  }
  
  if(any(grepl("spearman.test.genus.site.time.rclr.optA.list", to_import))){
    
    #by site and time, option A:
    {
      usvi_metab_site.time.list <- metabolomics_sample_metadata %>%
        dplyr::filter(grepl("seawater", sample_type)) %>%
        dplyr::filter(sample_id %in% rownames(usvi_genus_rclr.tbl)) %>%
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
        dplyr::filter(sample_id %in% rownames(usvi_genus_rclr.tbl)) %>%
        dplyr::mutate(grouping = paste0(site, ".", sampling_time)) %>%
        dplyr::distinct(sample_id, grouping) %>%
        droplevels %>%
        split(., f = .$grouping) %>%
        map(., ~.x %>%
              droplevels %>%
              dplyr::select(sample_id) %>%
              dplyr::inner_join(., usvi_genus_rclr.tbl %>%
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
        
        spearman.test.genus.site.time.list <- list(temp_spearman.test.df, 
                                                   temp_spearman.test, temp_spearman.test.rho, temp_spearman.test.n,
                                                   # temp_dend_asv, temp_dend_metab, 
                                                   padj_cutoff) %>%
          setNames(., c("spearman.test.df", 
                        "spearman.test", "spearman.test.rho", "spearman.test.n",
                        # "dend_asv", "dend_metab", 
                        "padj_cutoff"))
        
        assign(paste0("spearman.test.genus.siteA.rclr.", namevar, ".list"), spearman.test.genus.site.time.list, envir = .GlobalEnv)
        rm(spearman.test.genus.site.time.list)
        rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
      }
      
      
      spearman.test.genus.site.time.rclr.optA.list <- lapply(apropos("^spearman.test.genus.siteA.rclr.*$", mode = "list"),
                                                             get) %>%
        setNames(., names(usvi_asv_site.time.list))
    }
    
    readr::write_rds(spearman.test.genus.site.time.rclr.optA.list, paste0(projectpath, "/", "spearman.test.genus.site.time.rclr.optA.list", ".rds"), compress = "gz")
    
  }
  
  rm(list = apropos("^(spearman.test.genus.siteA.)(.*)$", mode = "list"))
  rm(list = apropos("^(temp_spearman.test)(.*)$", mode = "list"))
  rm(list = apropos("^(temp_)(.*)(.tbl)$", mode = "list"))
  rm(list = apropos("^(vector_m)(.*)(b)$", mode = "list"))
  rm(list = apropos("^(vector_m)(.*)(b_na)$", mode = "list"))
  
}



# Analyze genus-level Spearman correlations -------------------------------


if(!exists("spearman.test.genus.site.time.rclr.full.df", envir = .GlobalEnv)){
  if(any(grepl("spearman_rclr_full", list.files(projectpath, pattern = "usvi_.*.RData")))){
    load(list.files(projectpath, "usvi_spearman_rclr_full-.*.RData")[1])
    #spearman.test.genus.site.time.rclr.full.df, spearman.test.genus.site.rclr.full.df, spearman.test.rclr.full.df
    
  } else {
    spearman.test.genus.site.time.rclr.full.df <- spearman.test.genus.site.time.rclr.optA.list %>%
      map(.,
          ~.x[["spearman.test"]] %>%
            tibble::as_tibble(rownames = "asv_id") %>%
            tidyr::pivot_longer(., cols = -c("asv_id"),
                                names_to = "simpleName",
                                values_to = "p_value") %>%
            bind_rows(.)) %>%
      bind_rows(., .id = "grouping") %>%
      dplyr::left_join(., spearman.test.genus.site.time.rclr.optA.list %>%
                         map(.,
                             ~.x[["spearman.test.rho"]] %>%
                               tibble::as_tibble(rownames = "asv_id") %>%
                               tidyr::pivot_longer(., cols = -c("asv_id"),
                                                   names_to = "simpleName",
                                                   values_to = "estimate") %>%
                               bind_rows(.)) %>%
                         bind_rows(., .id = "grouping"),
                       by = join_by(asv_id, grouping, simpleName)) %>%
      dplyr::mutate(padj_bh = p.adjust(p_value, "BH")) %>%
      droplevels
    
    View(spearman.test.genus.site.time.rclr.full.df)
    
    #compare to the site-specific correlations
    spearman.test.genus.site.rclr.full.df <- spearman.test.genus.site.rclr.optA.list %>%
      map(.,
          ~.x[["spearman.test"]] %>%
            tibble::as_tibble(rownames = "asv_id") %>%
            tidyr::pivot_longer(., cols = -c("asv_id"),
                                names_to = "simpleName",
                                values_to = "p_value") %>%
            bind_rows(.)) %>%
      bind_rows(., .id = "grouping") %>%
      dplyr::left_join(., spearman.test.genus.site.rclr.optA.list %>%
                         map(.,
                             ~.x[["spearman.test.rho"]] %>%
                               tibble::as_tibble(rownames = "asv_id") %>%
                               tidyr::pivot_longer(., cols = -c("asv_id"),
                                                   names_to = "simpleName",
                                                   values_to = "estimate") %>%
                               bind_rows(.)) %>%
                         bind_rows(., .id = "grouping"),
                       by = join_by(asv_id, grouping, simpleName)) %>%
      dplyr::mutate(padj_bh = p.adjust(p_value, "BH"))
    
    #compare to the all-sample correlations
    spearman.test.genus.rclr.full.df <- spearman.test.genus.rclr.optA.list %>%
      purrr::pluck(., "spearman.test") %>%
      tibble::as_tibble(rownames = "asv_id") %>%
      tidyr::pivot_longer(., cols = -c("asv_id"),
                          names_to = "simpleName",
                          values_to = "p_value") %>%
      droplevels %>%
      dplyr::left_join(., (spearman.test.genus.rclr.optA.list %>%
                             purrr::pluck(., "spearman.test.rho") %>%
                             tibble::as_tibble(rownames = "asv_id") %>%
                             tidyr::pivot_longer(., cols = -c("asv_id"),
                                                 names_to = "simpleName",
                                                 values_to = "estimate") %>%
                             droplevels),
                       by = join_by(asv_id, simpleName)) %>%
      dplyr::mutate(padj_bh = p.adjust(p_value, "BH")) %>%
      dplyr::mutate(grouping = "all") %>%
      droplevels
    
    
    
    if(!any(grepl("spearman_genus_rclr_full", list.files(projectpath, pattern = "usvi_.*.RData")))){
      save(spearman.test.genus.site.time.rclr.full.df, spearman.test.genus.site.rclr.full.df, spearman.test.genus.rclr.full.df,
           file = paste0(projectpath, "/", "usvi_spearman_genus_rclr_full-", Sys.Date(), ".RData"))
    }
    
    
  }
}



if(!exists("spearman.test.genus.site.time.rclr.filtered.df", envir = .GlobalEnv)){
  
  if(any(grepl("spearman.genus.rclr.sig.filtered.list", list.files(projectpath, pattern = "usvi_.*.rds")))){
    temp_list <- readr::read_rds(list.files(projectpath, pattern = "usvi_spearman.genus.rclr.sig.filtered.list-.*.rds", full.names = TRUE)[1])
    list2env(temp_list, envir = .GlobalEnv)
    rm(temp_list)
    # spearman.test.rclr.filtered.lists <- list(spearman.test.genus.site.time.rclr.filtered.df,
    #                                           spearman.test.genus.site.rclr.filtered.df,
    #                                           spearman.test.rclr.filtered.df)
  } else {
    
    #recalculate FDRs
    padj_cutoff <- spearman.test.genus.rclr.full.df %>%
      dplyr::select(p_value) %>%
      tibble::deframe(.) %>% na.omit(.) %>%
      unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
      quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
      setNames(., c("q_01", "q_025", "q_05", "q_10")) #get the possible p-adj cutoffs for different q-values
    
    
    spearman.test.genus.rclr.df <- spearman.test.genus.rclr.full.df %>%
      tidyr::drop_na(p_value) %>%
      dplyr::mutate(across(c(asv_id, simpleName), ~factor(.x))) %>%
      droplevels %>%
      dplyr::rowwise(.) %>%
      # dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
      #               padj_01 = dplyr::case_when(p_value <= padj_cutoff["q_01"] ~ p_value, .default = NA),
      #               padj_025 = dplyr::case_when(p_value <= padj_cutoff["q_025"] ~ p_value, .default = NA),
      #               padj_05 = dplyr::case_when(p_value <= padj_cutoff["q_05"] ~ p_value, .default = NA),
      #               padj_10 = dplyr::case_when(p_value <= padj_cutoff["q_10"] ~ p_value, .default = NA)) %>%
      dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                    padj_01 = dplyr::case_when(p_value <= padj_cutoff["q_01"] ~ p_value, .default = NA),
                    padj_025 = dplyr::case_when(p_value <= padj_cutoff["q_025"] ~ p_value, .default = NA),
                    padj_05 = dplyr::case_when(p_value <= padj_cutoff["q_05"] ~ p_value, .default = NA),
                    padj_10 = dplyr::case_when(p_value <= padj_cutoff["q_10"] ~ p_value, .default = NA)) %>%
      tidyr::drop_na(padj_10) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
      tidyr::drop_na(estimate)
    
    # temp_df <- spearman.test.genus.rclr.df %>%
    spearman.test.genus.rclr.df <- spearman.test.genus.rclr.df %>%
      dplyr::rowwise(.) %>%
      dplyr::mutate(sig = dplyr::case_when(
        # !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
        !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
        !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
        !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
        !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
        .default = NA)) %>%
      dplyr::mutate(sig = dplyr::case_when(
        !is.na(padj_bh_05) & is.na(sig) ~ "vsig", #meaning that the adjusted p-value is below 0.05
        .default = sig)) %>%
      dplyr::mutate(sig = factor(sig)) %>%
      dplyr::arrange(asv_id, simpleName) %>%
      dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(asv_id, simpleName, .keep_all = TRUE) %>%
      droplevels %>%
      bind_rows(.)
    
    
    
    #Now site-specific
    padj_cutoff <- spearman.test.genus.site.rclr.full.df %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    
    spearman.test.genus.site.rclr.df <- spearman.test.genus.site.rclr.full.df %>%
      tidyr::drop_na(p_value) %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels)
    
    spearman.test.genus.site.rclr.df <- spearman.test.genus.site.rclr.df %>%
      imap(., ~.x %>%
             dplyr::mutate(across(c(asv_id, simpleName, grouping), ~factor(.x))) %>%
             droplevels %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                           padj_01 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_01"] ~ p_value, .default = NA),
                           padj_025 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_025"] ~ p_value, .default = NA),
                           padj_05 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_05"] ~ p_value, .default = NA),
                           padj_10 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_10"] ~ p_value, .default = NA)) %>%
             tidyr::drop_na(padj_10) %>%
             dplyr::ungroup(.) %>%
             dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                       abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
             tidyr::drop_na(estimate) %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(sig = dplyr::case_when(
               # !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
               !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
               .default = NA)) %>%
             dplyr::mutate(sig = dplyr::case_when(
               !is.na(padj_bh_05) & is.na(sig) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               .default = sig)) %>%
             dplyr::mutate(sig = factor(sig)) %>%
             dplyr::arrange(asv_id, simpleName) %>%
             dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
             dplyr::ungroup(.) %>%
             dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
             dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                           sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
             droplevels) %>%
      # setNames(., names(spearman.test.genus.site.df)) %>%
      bind_rows(.)
    
    #Now site- and time-specific
    padj_cutoff <- spearman.test.genus.site.time.rclr.full.df %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    
    spearman.test.genus.site.time.rclr.df <- spearman.test.genus.site.time.rclr.full.df %>%
      tidyr::drop_na(p_value) %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels)
    spearman.test.genus.site.time.rclr.df <- spearman.test.genus.site.time.rclr.df %>%
      imap(., ~.x %>%
             dplyr::mutate(across(c(asv_id, simpleName, grouping), ~factor(.x))) %>%
             droplevels %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(padj_bh_05 = dplyr::case_when(padj_bh <= 0.05 ~ padj_bh, .default = NA),
                           padj_01 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_01"] ~ p_value, .default = NA),
                           padj_025 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_025"] ~ p_value, .default = NA),
                           padj_05 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_05"] ~ p_value, .default = NA),
                           padj_10 = dplyr::case_when(p_value <= padj_cutoff[[.y]]["q_10"] ~ p_value, .default = NA)) %>%
             tidyr::drop_na(padj_10) %>%
             dplyr::ungroup(.) %>%
             dplyr::mutate(estimate = dplyr::case_when(abs(estimate) == 1 ~ NA, 
                                                       abs(round(1/estimate, digits = 3)) > 1 ~ estimate, .default = NA)) %>%
             tidyr::drop_na(estimate) %>%
             dplyr::rowwise(.) %>%
             dplyr::mutate(sig = dplyr::case_when(
               # !is.na(padj_bh_05) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               !is.na(padj_01) ~ "sig_q01", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_025) ~ "sig_q025", #meaning that q-tested p-value are below their respective thresholds
               !is.na(padj_05) ~ "sig_q05", #meaning that the q-tested p-value is below the 5% FDR
               !is.na(padj_10) ~ "maybe", #meaning that the q-tested p-value is below the 10% FDR
               .default = NA)) %>%
             dplyr::mutate(sig = dplyr::case_when(
               !is.na(padj_bh_05) & is.na(sig) ~ "vsig", #meaning that the adjusted p-value is below 0.05
               .default = sig)) %>%
             dplyr::mutate(sig = factor(sig)) %>%
             dplyr::arrange(asv_id, simpleName) %>%
             dplyr::filter(if_any(contains("padj"), ~!is.na(.x))) %>%
             dplyr::ungroup(.) %>%
             dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
             dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                           sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
             droplevels) %>%
      bind_rows(.)
    
    spearman.test.genus.rclr.filtered.df <- spearman.test.genus.rclr.df %>%
      dplyr::rowwise(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      # dplyr::slice_max(abs(estimate), by = c("asv_id", "simpleName", "sig")) %>%
      dplyr::mutate(across(c(asv_id, simpleName, sig), ~factor(.x))) %>%
      droplevels
    
    spearman.test.genus.site.rclr.filtered.df <- spearman.test.genus.site.rclr.df %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
    
    spearman.test.genus.site.time.rclr.filtered.df <- spearman.test.genus.site.time.rclr.df %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(abs(estimate) < 1) %>%
      dplyr::mutate(filtered_estimate = dplyr::case_when(abs(estimate) >= 0.5 ~ estimate,
                                                         .default = NA)) %>%
      tidyr::drop_na(filtered_estimate) %>%
      dplyr::mutate(site = stringr::str_split_i(grouping, "\\.", 1),
                    sampling_time = stringr::str_split_i(grouping, "\\.", 2)) %>%
      dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
      dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
      dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
      droplevels
    
    # if(!any(grepl("spearman.rclr.sig.filtered.list", list.files(projectpath, pattern = "usvi_.*.rds")))){
    spearman.test.genus.rclr.filtered.lists <- list(spearman.test.genus.site.time.rclr.filtered.df,
                                              spearman.test.genus.site.rclr.filtered.df,
                                              spearman.test.genus.rclr.filtered.df) %>%
      setNames(., c("spearman.test.genus.site.time.rclr.filtered.df",
                    "spearman.test.genus.site.rclr.filtered.df",
                    "spearman.test.genus.rclr.filtered.df"))
    readr::write_rds(spearman.test.genus.rclr.filtered.lists, 
                     paste0(projectpath, "/", "usvi_spearman.rclr.sig.filtered.list-", Sys.Date(), ".rds"), 
                     compress = "gz")
    
    for(df in c("spearman.test.genus.site.time.rclr.filtered.df",
                "spearman.test.genus.site.rclr.filtered.df",
                "spearman.test.genus.rclr.filtered.df")){
      temp_df <- get0(df, mode = "any")
      readr::write_delim(temp_df, paste0(projectpath, "/", "usvi_", df, "-", Sys.Date(), ".tsv"),
                         col_names = TRUE, delim = "\t")  
    }
    # }
  }
  
}

if(!exists("padj_cutoff_list", envir = .GlobalEnv)){
  if(any(grepl("spearman.test.genus.rclr.padj_cutoff", list.files(projectpath, pattern = "usvi_.*.rds")))){
    padj_cutoff_list <- readr::read_rds(list.files(projectpath, "usvi_spearman.test.genus.rclr.padj_cutoff-.*.rds")[1])
  } else {
    
    #save the calculated FDRs as well
    
    padj_cutoff <- spearman.test.genus.rclr.full.df %>%
      dplyr::select(p_value) %>%
      tibble::deframe(.) %>% na.omit(.) %>%
      unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
      quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
      setNames(., c("q_01", "q_025", "q_05", "q_10")) #get the possible p-adj cutoffs for different q-values
    temp_df <- spearman.test.genus.site.rclr.full.df %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            droplevels %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    temp_df2 <- spearman.test.genus.site.time.rclr.full.df %>%
      split(., f = .$grouping) %>%
      map(., ~.x %>%
            dplyr::select(p_value) %>%
            tibble::deframe(.) %>% na.omit(.) %>%
            unlist %>% ashr::qval.from.lfdr(.) %>% as.matrix(.) %>%
            quantile(., probs = c(0.01, 0.025, 0.05, 0.1), na.rm = TRUE, names = FALSE,type = 7) %>%
            setNames(., c("q_01", "q_025", "q_05", "q_10"))) #get the possible p-adj cutoffs for different q-values
    
    padj_cutoff_list <- padj_cutoff %>%
      list(.) %>% 
      setNames(., c("all")) %>%
      append(., temp_df) %>%
      append(., temp_df2)
    readr::write_rds(padj_cutoff_list, paste0(projectpath, "/", "usvi_spearman.test.genus.rclr.padj_cutoff-", Sys.Date(), ".rds"))
    
    padj_cutoff <- padj_cutoff_list %>%
      dplyr::bind_rows(., .id = "model")
    
    readr::write_delim(padj_cutoff, paste0(projectpath, "/", "usvi_spearman.test.genus.rclr.padj_cutoff-", Sys.Date(), ".tsv"),
                       delim = "\t", col_names = TRUE, num_threads = nthreads)
    
  }
  
}


if(!exists("spearman.test.genus.site.time.rclr.df", envir = .GlobalEnv)){
  
  # temp_idx <- c("asv_id","simpleName", "p_value","estimate","padj_bh","grouping","padj_bh_05", 
  #               "padj_01","padj_025","padj_05","padj_10","sig")
  temp_drop <- c("filtered_estimate")
  
  spearman.test.genus.rclr.df <- spearman.test.genus.rclr.filtered.df %>%
    dplyr::select(-c(all_of(temp_drop))) %>%
    # dplyr::select(all_of(temp_idx)) %>%
    droplevels
  spearman.test.genus.site.rclr.df <- spearman.test.genus.site.rclr.filtered.df %>%
    dplyr::select(-c(all_of(temp_drop))) %>%
    # dplyr::select(all_of(temp_idx), "site") %>%
    droplevels
  spearman.test.genus.site.time.rclr.df <- spearman.test.genus.site.time.rclr.filtered.df %>%
    dplyr::select(-c(all_of(temp_drop))) %>%
    # dplyr::select(all_of(temp_idx), "site", "sampling_time") %>%
    droplevels
  
  
}


#now combine them

spearman.rclr.all.rho.p.df <- spearman.test.genus.rclr.full.df %>%
  dplyr::filter(p_value <= 0.1) %>%
  dplyr::distinct(asv_id, simpleName) %>%
  droplevels %>%
  dplyr::inner_join(., spearman.test.genus.site.time.rclr.full.df %>%
                      dplyr::filter(p_value <= 0.1) %>%
                      dplyr::distinct(asv_id, simpleName) %>%
                      droplevels, 
                    by = join_by(asv_id, simpleName)) %>%
  droplevels %>%
  dplyr::inner_join(., spearman.test.genus.site.rclr.full.df %>%
                      dplyr::filter(p_value <= 0.1) %>%
                      dplyr::distinct(asv_id, simpleName) %>%
                      droplevels, 
                    by = join_by(asv_id, simpleName)) %>%
  droplevels %>%
  dplyr::left_join(., spearman.test.genus.rclr.full.df %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels, 
                   by = join_by(asv_id, simpleName)) %>%
  dplyr::bind_rows(., spearman.test.genus.site.time.rclr.full.df %>%
                     dplyr::filter(p_value <= 0.1) %>%
                     dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels) %>%
  dplyr::bind_rows(., spearman.test.genus.site.rclr.full.df %>%
                     dplyr::filter(p_value <= 0.1) %>%
                     dplyr::distinct(asv_id, simpleName, grouping, .keep_all = TRUE) %>%
                     tidyr::pivot_longer(., cols = c(p_value, padj_bh, estimate),
                                         names_to = "metric",
                                         values_to = "value") %>%
                     droplevels) %>%
  droplevels

spearman.rclr.all.rho.p.df <- spearman.rclr.all.rho.p.df %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(asv_id, simpleName, grouping) %>%
  dplyr::filter(grepl("all", grouping), .by = c(asv_id, simpleName)) %>%
  dplyr::distinct(asv_id, simpleName) %>%
  dplyr::left_join(., spearman.rclr.all.rho.p.df,
                   by = join_by(asv_id, simpleName)) %>%
  tidyr::pivot_wider(., id_cols = c(asv_id, simpleName, grouping),
                     names_from = "metric",
                     values_from = "value") %>%
  dplyr::mutate(grouping = fct_relevel(grouping, "all", after = Inf)) %>%
  droplevels

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

spearman.test.genus.rclr.bh.dist.df <- bind_rows(spearman.test.genus.rclr.df, spearman.test.genus.site.rclr.df) %>%
  bind_rows(., spearman.test.genus.site.time.rclr.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  tidyr::drop_na(padj_10) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  # tidyr::drop_na(padj_bh_05) %>%
  # dplyr::ungroup(.) %>%
  # dplyr::arrange(desc(padj_bh_05)) %>%
  dplyr::arrange(desc(padj_bh_05), .by_group = TRUE) %>%
  dplyr::select(grouping, site, sampling_time, p_value, padj_10, padj_05, padj_025, padj_01, padj_bh_05) %>%
  dplyr::distinct(grouping, site, sampling_time, .keep_all = TRUE) %>%
  droplevels

#plot the distribution of Spearman rho's when using all samples, compared to site-specific
spearman.test.genus.rclr.rho.dist.df <- bind_rows(spearman.test.genus.rclr.filtered.df, spearman.test.genus.site.rclr.filtered.df) %>%
  bind_rows(., spearman.test.genus.site.time.rclr.filtered.df) %>%
  dplyr::mutate(site = dplyr::case_when(is.na(site) ~ "all", .default = site),
                sampling_time = dplyr::case_when(is.na(sampling_time) ~ "all", .default = sampling_time),
                grouping = dplyr::case_when(is.na(grouping) ~ "all", .default = grouping)) %>%
  dplyr::mutate(site = factor(site, levels = c(names(site_lookup), "all"))) %>%
  dplyr::mutate(sampling_time = factor(sampling_time, levels = c(names(sampling_time_lookup), "all"))) %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig, site, sampling_time), ~factor(.x))) %>%
  dplyr::distinct(asv_id, simpleName, grouping, site, sampling_time, sig, .keep_all = TRUE) %>%
  dplyr::select(asv_id, simpleName, grouping, site, sampling_time, estimate, sig) %>%
  # dplyr::mutate(sig = factor(sig, levels = names(padj_cutoff_labels))) %>%
  dplyr::group_by(asv_id, simpleName,  grouping,site, sampling_time, sig) %>%
  # dplyr::distinct(estimate, .keep_all = TRUE) %>%
  droplevels

spearman.test.genus.rclr.rho.dist.df %>%
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

spearman.test.genus.rclr.rho.dist.summary.df <- spearman.test.genus.rclr.rho.dist.df %>%
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
spearman.test.genus.rclr.rho.dist.summary.df %>%
  dplyr::filter(grepl("all", sampling_time))
# # A tibble: 16 × 5
# grouping    site         sampling_time sig      num_obs
# <fct>       <fct>        <fct>         <fct>      <int>
#   1 LB_seagrass Seagrass     all           sig_q05      205
# 2 LB_seagrass Seagrass     all           sig_q01        0
# 3 LB_seagrass Seagrass     all           maybe        522
# 4 LB_seagrass Seagrass     all           sig_q025      97
# 5 Yawzi       Yawzi Reef   all           sig_q05      348
# 6 Yawzi       Yawzi Reef   all           sig_q01      133
# 7 Yawzi       Yawzi Reef   all           maybe        808
# 8 Yawzi       Yawzi Reef   all           sig_q025     227
# 9 Tektite     Tektite Reef all           sig_q05      240
# 10 Tektite     Tektite Reef all           sig_q01      103
# 11 Tektite     Tektite Reef all           maybe        246
# 12 Tektite     Tektite Reef all           sig_q025     218
# 13 all         all          all           sig_q05      489
# 14 all         all          all           sig_q01      129
# 15 all         all          all           maybe        831
# 16 all         all          all           sig_q025     294


#now for the rCLR transformed + Spearman correlations

#note that the adjusted p-value for FDR q = 5% is now smaller than the BH-adjusted p-value from the Spearman test.
#so filtering for (!is.na(padj_bh_05) & !is.na(sig_q01)) is much more stringent

#union(maybe, sig_q01, sig_q025, sig_05)
# [1] 3147
corr_fdr01_genus_rCLR_site_sig <- spearman.test.genus.site.rclr.filtered.df %>%
  dplyr::mutate(across(c(asv_id, simpleName, grouping, sig), ~factor(.x))) %>%
  # filter(sig %in% c("sig_q01","vsig")) %>% # [1] 236 rows
  # filter(sig %in% c("sig_q01","sig_q05")) %>% # [1] 1029 rows
  # filter(sig %in% c("sig_q01","sig_q025")) %>% # [1] 778 rows
  dplyr::filter(abs(estimate) >= 0.5) %>%  filter(sig %in% c("sig_q01","sig_q025", "sig_q05")) %>% # [1] 1571
    # filter(sig %in% c("sig_q01")) %>% # [1] 236 rows
  # filter(sig %in% c("sig_q025")) %>% # [1] 542 rows
  droplevels
nrow(corr_fdr01_genus_rCLR_site_sig)

corr_fdr01_genus_rCLR_site_sig_unique_df <- corr_fdr01_genus_rCLR_site_sig %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  split(., f = .$corrlink) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::summarise(num_sites = length(unique(.[["grouping"]]))) %>%
        droplevels
  ) %>%
  bind_rows(., .id = "corrlink") %>%
  dplyr::arrange(desc(num_sites))


corr_fdr01_genus_rCLR_site_freq_sig_df <- dplyr::full_join((corr_fdr01_genus_rCLR_site_sig %>%
                                                              dplyr::ungroup(.) %>%
                                                      distinct(simpleName, asv_id) %>% #this deduplicates observations of significant ASV-metab correlations across multiple sites
                                                      count(asv_id, name = "num_unique_correlations") %>%   # Count occurrences
                                                      arrange(desc(num_unique_correlations))),
                                                   (corr_fdr01_genus_rCLR_site_sig %>%
                                                      dplyr::ungroup(.) %>%
                                                      count(asv_id, name = "num_correlations") %>%   # Count occurrences
                                                      arrange(desc(num_correlations)))) %>%
  # dplyr::filter(num_correlations >= 10) %>% #was 10 for the MH-transformed ASV-metabolite 
  dplyr::filter(num_correlations >= 2) %>%
  droplevels

#so which ASVs are observed with multiple correlations to metabolites?
corr_fdr01_genus_rCLR_site_sig_freq_unique_df <- corr_fdr01_genus_rCLR_site_sig %>%
  dplyr::mutate(corrlink = paste0(asv_id, ":", simpleName)) %>%
  dplyr::mutate(corrlink = factor(corrlink)) %>%
  dplyr::filter(asv_id %in% corr_fdr01_genus_rCLR_site_freq_sig_df[["asv_id"]]) %>%
  # dplyr::left_join(., corr_fdr01_site_freq_sig_df, by = join_by(asv_id), multiple = "all", relationship = "many-to-many") %>% 
  dplyr::distinct(corrlink, asv_id, simpleName, grouping, .keep_all = TRUE) %>%
  droplevels 

#142 pASVs-metabolite correlations are strong, significant, and observed in 2 of 3 sites 
corr_fdr01_genus_rCLR_site_sig_freq_unique_df %>%
  count(corrlink,name = "num_occ") %>%
  dplyr::arrange(desc(num_occ)) %>%
  dplyr::filter(num_occ > 1)
# # A tibble: 142 × 2
# corrlink                    num_occ
# <fct>                         <int>
#   1 pASV_00002:putrescine 2           2
# 2 pASV_00007:pantothenic acid       2


usvi_SDE_genus_rCLR_venn.list <- corr_fdr01_genus_rCLR_site_sig %>%
  dplyr::mutate(pairing = paste(asv_id, simpleName, sep = ":")) %>%
  split(., f = .$grouping) %>%
  map(., ~.x %>%
        # droplevels %>%
        split(., f = .$sig) %>%
        map(., ~.x %>%
              dplyr::distinct(pairing) %>%
              # dplyr::distinct(asv_id, simpleName) %>%
              # dplyr::distinct(asv_id) %>%
              # Filter(Negate(function(x) is.null(unlist(x))), .) %>% #thanks https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
              tibble::deframe(.) %>%
              unlist) %>%
        purrr::list_flatten(.)
  ) %>%
  purrr::list_flatten(.)




# [1] "Seagrass FDR 5%"       "Seagrass FDR 2.5%"     "Seagrass FDR 1%"      
# [4] "Tektite Reef FDR 5%"   "Tektite Reef FDR 2.5%" "Tektite Reef FDR 1%"  
# [7] "Yawzi Reef FDR 5%"     "Yawzi Reef FDR 2.5%"   "Yawzi Reef FDR 1%" 
usvi_SDE_genus_rCLR_venn.list <- usvi_SDE_genus_rCLR_venn.list %>%
  setNames(., (names(usvi_SDE_genus_rCLR_venn.list) %>%
                 stringr::str_replace_all(., site_lookup) %>%
                 stringr::str_replace_all(., padj_cutoff_labels) %>%
                 stringr::str_replace_all(., "_Q0", " FDR ") %>%
                 stringr::str_c(., "%")))

summary(usvi_SDE_genus_rCLR_venn.list)
# Length Class  Mode     
# Seagrass FDR 5%       205    -none- character
# Seagrass FDR 2.5%      97    -none- character
# Seagrass FDR 1%         0    -none- character
# Tektite Reef FDR 5%   240    -none- character
# Tektite Reef FDR 2.5% 218    -none- character
# Tektite Reef FDR 1%   103    -none- character
# Yawzi Reef FDR 5%     348    -none- character
# Yawzi Reef FDR 2.5%   227    -none- character
# Yawzi Reef FDR 1%     133    -none- character
# upset_length <- names(usvi_SDE_genus_rCLR_venn.list) %>%
#   combn(., 2) %>%
#   t() %>% as.data.frame(.) %>%
#   setNames(., c("pair1", "pair2")) %>%
#   droplevels %>%
#   dplyr::rowwise(.) %>%
#   dplyr::mutate(n_intersect = length(intersect(usvi_SDE_genus_rCLR_venn.list[[{pair1}]], usvi_SDE_genus_rCLR_venn.list[[{pair2}]]))) %>%
#   dplyr::filter(n_intersect > 0) %>%
#   dplyr::bind_rows(., data.frame("pair1" = c(names(usvi_SDE_genus_rCLR_venn.list)),
#                                  "pair2" = NA)) %>%
#   dplyr::mutate(n_intersect = dplyr::case_when(is.na(n_intersect) ~ length(usvi_SDE_genus_rCLR_venn.list[[{pair1}]]),
#                                                .default = n_intersect)) %>%
#   dplyr::arrange(pair1) %>%
#   droplevels
    

upset_length <- data.frame("pair1" = c(names(usvi_SDE_genus_rCLR_venn.list)),
                           "pair2" = NA) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(n_intersect = length(usvi_SDE_genus_rCLR_venn.list[[ {pair1} ]])) %>%
  dplyr::mutate( pair1 = factor(pair1)) %>%
  dplyr::arrange(pair1) %>%
  dplyr::bind_rows(., (names(usvi_SDE_genus_rCLR_venn.list) %>%
                         combn(., 2) %>%
                         t() %>% as.data.frame(.) %>%
                         setNames(., c("pair1", "pair2")) %>%
                         droplevels %>%
                         dplyr::rowwise(.) %>%
                         dplyr::mutate(n_intersect = length(intersect(usvi_SDE_genus_rCLR_venn.list[[{pair1}]], usvi_SDE_genus_rCLR_venn.list[[{pair2}]]))) %>%
                         dplyr::filter(n_intersect > 0) %>%
                         droplevels)) %>%
  dplyr::arrange(pair1) %>%
  dplyr::filter(n_intersect > 0) %>%
  droplevels


temp_g6 <- UpSetR::upset(UpSetR::fromList(usvi_SDE_genus_rCLR_venn.list),
                         # main.bar.color = "black",
                         mainbar.y.label = "Number of strong, significant\ngenus-metabolite correlations",
                         set.metadata = NULL,
                         sets.x.label = "Set size",
                         sets = c(names(usvi_SDE_genus_rCLR_venn.list)),
                         empty.intersections = NULL,
                         order.by = "degree", decreasing = F,
                         main.bar.color = viridis::turbo(n = nrow(upset_length)),
                         scale.intersections = "identity",
                         keep.order = T, set_size.show = FALSE,
                         query.legend = "bottom")
print(temp_g6)
# temp_g6 <- ggplotify::as.ggplot(temp_g6)
# glayout <- "
# 111111
# #22222
# "
# ggplotify::as.ggplot(temp_g2[[1]]) + ggplotify::as.ggplot(temp_g2[[2]]) + patchwork::plot_layout(design = glayout)
# g1A <- ggplotify::as.ggplot(temp_g2[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(temp_g2[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = glayout)
# g2A <- ggplotify::as.ggplot(temp_g6[[1]]) + patchwork::wrap_elements(ggplotify::as.ggplot(temp_g6[[2]]), ignore_tag = TRUE) + patchwork::plot_layout(design = glayout)

g1D <- ggplotify::as.ggplot(temp_g1) + ggtitle("Metabolite correlations\n with ASV relative abundance data")
g1A <- ggplotify::as.ggplot(temp_g2) + ggtitle("Metabolite correlations\n with ASV relative abundance or rCLR-transformed data")
g1B <- ggplotify::as.ggplot(temp_g6) + ggtitle("Genus-Metabolite correlations\n with rCLR-transformed data")
g1C <- ggplotify::as.ggplot(temp_g3) + ggtitle("ASV-Metabolite correlations\n with rCLR-transformed data")

gpatch <- (g1D | g1A | g1C | g1B) + patchwork::plot_annotation(tag_level = "A")

# gpatch <- (ggplotify::as.ggplot(temp_g2[[1]]) + ggtitle("ASV-Metabolite correlations\n with MH- or rCLR-transformed data")) /
#   ((ggplotify::as.ggplot(temp_g6) + ggtitle("ASV-Metabolite correlations\n with rCLR-transformed data")) +
#   (ggplotify::as.ggplot(temp_g3) + ggtitle("Genus-Metabolite correlations\n with rCLR-transformed data")))

gpatch

ggsave(paste0(projectpath, "/", "usvi_spearman_rCLR_vs_MH_asv_vs_genus-", Sys.Date(), ".png"),
       gpatch,
       width = 20, height =6, units = "in")



## do the ASVs in rCLR-transformed correlations fit into the genera that were significantly correlated?

#142 pASVs-metabolite correlations are strong, significant, and observed in 2 of 3 sites 
temp_df <- corr_fdr01_genus_rCLR_site_sig_freq_unique_df %>%
  count(corrlink,name = "num_occ") %>%
  dplyr::arrange(desc(num_occ)) %>%
  dplyr::filter(num_occ > 1) %>%
  dplyr::select(corrlink) %>%
  tidyr::separate_wider_delim(., corrlink, names = c("p_asv_id", "simpleName"), delim = ":", cols_remove = FALSE) %>%
  # dplyr::rename("p_asv_id" = "asv_id") %>%
  dplyr::left_join(., usvi_prok_genus.taxa %>%
                     dplyr::select(asv_id, p_asv_id), 
                   # by = join_by("asv_id" == "p_asv_id"),
                   relationship = "many-to-many", multiple = "all") %>%
  
  droplevels

temp_df2 <- temp_df %>%
  dplyr::full_join(., corr_fdr01_rCLR_site_sig_freq_unique_df %>%
                     dplyr::select(asv_id, grouping, simpleName, corrlink),
                   by = join_by(asv_id, simpleName),
                   # relationship = "many-to-many", multiple = "all",
                   suffix = c(".genus", ".rCLR")) %>%
  dplyr::full_join(., corr_fdr01_MH_site_sig_freq_unique_df %>%
                     dplyr::select(asv_id, grouping, simpleName, corrlink) %>%
                     droplevels,
                   by = join_by(asv_id, grouping, simpleName),
                   # relationship = "many-to-many", multiple = "all",
                   suffix = c(".rCLR", ".MH")) %>%
  droplevels
temp_df3 <- temp_df2 %>%
  # tidyr::drop_na(corrlink) %>%
  dplyr::ungroup(.) %>%
  tidyr::drop_na(`corrlink.rCLR`) %>%
  tidyr::drop_na(`corrlink.genus`) %>%
  droplevels

temp_df3 %>%
  tidyr::drop_na(`corrlink.rCLR`) %>%
  dplyr::distinct(asv_id, p_asv_id, simpleName, grouping, `corrlink.rCLR`) %>%
  dplyr::summarise(n_sites = length(grouping), .by = "corrlink.rCLR") %>%
  dplyr::arrange(desc(n_sites))
# # A tibble: 51 × 2
# corrlink.rCLR              n_sites
# <fct>                        <int>
#   1 ASV_06661:pantothenic acid       2
# 2 ASV_06670:pantothenic acid       2
# 3 ASV_04243:pantothenic acid       2
# 4 ASV_07040:pantothenic acid       2
# 5 ASV_10272:pantothenic acid       2
# 6 ASV_06650:pantothenic acid       2
# 7 ASV_12551:pantothenic acid       2
# 8 ASV_10264:pantothenic acid       2
# 9 ASV_11900:pantothenic acid       2
# 10 ASV_03006:pantothenic acid       2
# # ℹ 41 more rows

temp_df2 %>%
  tidyr::drop_na(`corrlink.rCLR`) %>%
  dplyr::distinct(asv_id, p_asv_id, simpleName, grouping, `corrlink.rCLR`) %>%
  dplyr::summarise(n_sig_asvs = length(asv_id), .by = "p_asv_id") %>%
  tidyr::drop_na(p_asv_id) %>%
  dplyr::arrange(desc(n_sig_asvs))
# # A tibble: 41 × 2
# p_asv_id   n_sig_asvs
# <chr>           <int>
#   1 pASV_04690          6
# 2 pASV_01689          4
# 3 pASV_02009          4
# 4 pASV_02127          4
# 5 pASV_00370          3
# 6 pASV_00524          3
# 7 pASV_01660          3
# 8 pASV_03097          3
# 9 pASV_03100          3
# 10 pASV_00359          2
# # ℹ 31 more rows
# # ℹ Use `print(n = ...)` to see more rows



# Work with the ASV-level rCLR Spearman correlations ----------------------


spearman.test.rclr.rho.dist.df %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  dplyr::reframe(bin_est = quantile(abs_estimate, probs = seq(0, 1, 0.25), names = FALSE))
# # A tibble: 5 × 1
# bin_est
# <dbl>
#   1   0.517
# 2   0.622
# 3   0.682
# 4   0.789
# 5   0.964

#how does that compare to the relative abundance ASV Spearman correlation results?
spearman.test.site.filtered.df %>%
  dplyr::filter(!(sig %in% c("maybe"))) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  dplyr::reframe(bin_est = quantile(abs_estimate, probs = seq(0, 1, 0.25), names = FALSE))
# # A tibble: 5 × 1
# bin_est
# <dbl>
#   1   0.5  
# 2   0.532
# 3   0.580
# 4   0.668
# 5   0.986

#the previous threshold was abs(rho) >= 0.5, and further refinement of abs(rho) >= 0.9 retained 18 ASVs with 19 ASV-metabolite correlations
#904 unique ASVs had 910 very strong (abs(rho) >= 0.9)) correlations with 10 unique metabolites
spearman.test.rclr.rho.dist.df %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  # dplyr::filter(abs_estimate >= 0.8) %>% # [1] 22637 correlations btw 7398 unique ASVs and 30 metabolites
  # dplyr::filter(abs_estimate >= 0.85) %>% # [1] 7341 correlations btw 4507 unique ASVs and 23 metabolites
  dplyr::filter(abs_estimate >= 0.9) %>% # [1] 910 correlations btw 904 unique ASVs and 10 metabolites
  # dplyr::distinct(asv_id) %>% # [1] 904
  # dplyr::distinct(simpleName) %>% # [1] 10
  nrow(.)

spearman.test.site.filtered.df %>%
  dplyr::filter(!(sig %in% c("maybe"))) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  dplyr::filter(abs_estimate >= 0.8) %>% # [1] 241 correlations btw 201 unique ASVs and 14 metabolites
  # dplyr::filter(abs_estimate >= 0.85) %>% # [1] 160 correlations btw 142 unique ASVs and 12 metabolites
  # dplyr::filter(abs_estimate >= 0.9) %>% # [1] 38 correlations btw 37 unique ASVs and 7 metabolites
  # dplyr::distinct(asv_id) %>% # [1] 904
  dplyr::distinct(simpleName) %>% # [1] 10
  nrow(.)

#instead of filtering by strength of correlation, what about looking at taxonomy of ASVs with strong correlations in rCLR-Spearman analysis?
#a majority of correlations are with Bacteria;Pseudomonadota, Bacteria;Bacteroidota, and Bacteria;Verrucomicrobia
#filter by abs(rho) >= 0.8
spearman.test.rclr.asv.taxonomy.df <- spearman.test.rclr.rho.dist.df %>%
  dplyr::filter(!(site == "all") & !(sampling_time == "all")) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  dplyr::filter(abs_estimate >= 0.8) %>% 
  droplevels %>%
  # dplyr::summarize(obs_asv = length(estimate), .by = "asv_id") %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(all_of(keep_tax), "asv_id") %>%
                     droplevels,
                   by = join_by(asv_id)) %>%
  # otu_to_taxonomy(., usvi_prok_asvs.taxa, level = "Genus") %>%
  # dplyr::filter(obs_asv > 0) %>%
  dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy_string")) %>%
  dplyr::mutate(taxonomy_string = stringr::str_remove_all(taxonomy_string, "^ASV_(.....): ")) %>%
  droplevels

#what does the taxonomy look like in a sankey graph?
{
  # temp_df <- spearman.test.rclr.rho.dist.df %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::mutate(abs_estimate = abs(estimate)) %>%
  #   dplyr::filter(abs_estimate >= 0.8) %>% 
  #   droplevels %>%
  #   dplyr::summarize(obs_asv = length(estimate), .by = "asv_id") %>%
  #   # dplyr::left_join(., usvi_prok_filled.taxa.df %>%
  #   #                    dplyr::select(all_of(keep_tax), "asv_id") %>%
  #   #                    droplevels,
  #   #                  by = join_by(asv_id)) %>%
  #   # dplyr::filter(grepl("Woese", Order)) %>% dplyr::summarise(num_obs = sum(obs_asv)) #3198 correlations >0.8 or < -0.8 that are with Woesearchaeales ASVs
  #   otu_to_taxonomy(., usvi_prok_asvs.taxa, level = "Genus") %>%
  #   dplyr::filter(obs_asv > 0) %>%
  #   droplevels %>%
  #   dplyr::arrange(desc(obs_asv)) %>%
  #   dplyr::mutate(across(all_of(keep_tax), ~factor(.))) %>%
  #   ggsankey::make_long(., all_of(keep_tax), value = "obs_asv") %>%
  #   dplyr::mutate(label = dplyr::case_when(value > 500 ~ node,
  #                                          .default = NA)) %>%
  #   droplevels
  # 
  # 
  # (ggplot(data = temp_df,
  #         aes(x = x, next_x = next_x,
  #             node = node, next_node = next_node,
  #             label = label,
  #             color = factor(node),
  #             fill = factor(node)))
  #   # + ggsankey::geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)
  #   # + ggsankey::geom_sankey_label(size = 3, color = "white", fill = "gray40")
  #   # + ggsankey::theme_sankey(base_size = 18)
  #   + ggsankey::geom_alluvial(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)
  #   + ggsankey::geom_alluvial_label(size = 3, color = "white", fill = "gray40")
  #   + ggsankey::theme_alluvial(base_size = 18)
  #   + theme_dark()
  #   + theme(
  #     axis.title = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks = element_blank(),
  #     panel.grid = element_blank())
  # )
  
}

#now examine the number of genus-metabolite strong correlations
#too many Archaea are represented here, when before we had no correlations with Archaea
# 
# spearman.test.rclr.asv.taxonomy.df %>%
#   dplyr::group_by(taxonomy_string, simpleName) %>%
#   dplyr::summarise(num_corr = length(estimate)) %>%
#   dplyr::arrange(desc(num_corr))
# # # A tibble: 4,812 × 3
# # # Groups:   taxonomy_string [790]
# # taxonomy_string                     simpleName              num_corr
# # <chr>                               <fct>                      <int>
# #   1 Nanoarchaeia; Woesearchaeales       sn-glycerol 3-phosphate      686
# # 2 Nanoarchaeia; Woesearchaeales       arginine                     414
# # 3 Bacteria                            sn-glycerol 3-phosphate      379
# # 4 Nanoarchaeia; Woesearchaeales       homoserine betaine           277
# # 5 Proteobacteria; Gammaproteobacteria sn-glycerol 3-phosphate      265
# # 6 Proteobacteria; Alphaproteobacteria sn-glycerol 3-phosphate      221
# # 7 Bacteria                            arginine                     217
# # 8 Nanoarchaeia; Woesearchaeales       5'UMP                        210
# #  9 Nanoarchaeia; Woesearchaeales       taurine                      210
# # 10 Bacteria                            homoserine betaine           160
# # # ℹ 4,802 more rows

#what did it look like before??
spearman.test.relabund.asv.taxonomy.df <- spearman.test.site.filtered.df %>%
  dplyr::filter(!(sig %in% c("maybe"))) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  dplyr::filter(abs_estimate >= 0.8) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(all_of(keep_tax), "asv_id") %>%
                     droplevels,
                   by = join_by(asv_id)) %>%
  # otu_to_taxonomy(., usvi_prok_asvs.taxa, level = "Genus") %>%
  # dplyr::filter(obs_asv > 0) %>%
  dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy_string")) %>%
  dplyr::mutate(taxonomy_string = stringr::str_remove_all(taxonomy_string, "^ASV_(.....): ")) %>%
  droplevels  

spearman.test.relabund.asv.taxonomy.df %>%
  dplyr::group_by(taxonomy_string, simpleName) %>%
  dplyr::summarise(num_corr = length(estimate)) %>%
  dplyr::arrange(desc(num_corr))
# # A tibble: 183 × 3
# # Groups:   taxonomy_string [93]
# taxonomy_string                         simpleName         num_corr
# <chr>                                   <chr>                 <int>
#   1 Alphaproteobacteria; SAR116 clade       GABA                      6
# 2 Alphaproteobacteria; SAR11 Clade II     3'AMP                     5
#  3 Bacteroidia; Cryomorphaceae             3'AMP                     5
# 4 Bacteroidia; NS5 marine group           3'AMP                     5
#  5 Bacteria; Marinimicrobia (SAR406 clade) GABA                      4
#  6 Bdellovibrionia; OM27 clade             5'AMP                     4
# 7 Gammaproteobacteria; Vibrio             GABA                      4
# 8 Alphaproteobacteria; Rhodobacteraceae   ciliatine                 3


#with rCLR transformation,the distribution of rare ASVs is not dissimilar to the abundant ASVs.
#0.00123 : the 50% of those observations for all ASVs Archaea;Nanoarchaeota;Nanoarchaeia;Woesearchaeales 
#-0.000162 : the 50% of those observations for all ASVs Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales
#0.00421 : the 50% of those observations for all ASVs Bacteria;Pseudomonadota;Alphaproteobacteria;SAR11 clade
#these aren't that different -_-
{
#   
# temp_df <- usvi_prok_filled.taxa.df %>%
#   dplyr::select(all_of(keep_tax), "asv_id") %>%
#   # dplyr::filter(grepl("Woese", Order)) %>% 
#   droplevels %>%
#   dplyr::left_join(., usvi_asv_rclr.tbl %>%
#                      tibble::as_tibble(rownames = "sample_id") %>%
#                      # tibble::rownames_to_column(var = "sample_id") %>%
#                      tidyr::pivot_longer(., cols = !c("sample_id"),
#                                          names_to = "asv_id",
#                                          values_to = "rclr_trans"),
#                    by = join_by("asv_id")) %>%
#   droplevels
# 
# temp_df %>%
#   dplyr::filter(grepl("Woese|SAR11|Pseudomonadales", Order)) %>%
#   dplyr::reframe(rclr_options = quantile(rclr_trans, probs = seq(0, 1, 0.25), na.rm = TRUE, names = TRUE), .by = "Order") %>%
#   droplevels
# # # A tibble: 20 × 2
# # Order                    rclr_options
# # <fct>                           <dbl>
# #   1 SAR11 clade                 -3.23    
# # 2 SAR11 clade                 -0.0261  
# # 3 SAR11 clade                  0.00421 
# # 4 SAR11 clade                  0.0944  
# # 5 SAR11 clade                  1.99    
# # 6 Pseudomonadales             -6.29    
# # 7 Pseudomonadales             -0.0271  
# # 8 Pseudomonadales             -0.000162
# # 9 Pseudomonadales              0.0315  
# # 10 Pseudomonadales              5.77    
# # 11 Woesearchaeales             -0.676   
# # 12 Woesearchaeales             -0.0341  
# # 13 Woesearchaeales              0.00123 
# # 14 Woesearchaeales              0.0365  
# # 15 Woesearchaeales              0.643   
# # 16 Candidatus Woesebacteria    -0.216   
# # 17 Candidatus Woesebacteria    -0.0490  
# # 18 Candidatus Woesebacteria     0.00209 
# # 19 Candidatus Woesebacteria     0.0753  
# # 20 Candidatus Woesebacteria     0.169  
# 
# temp_df2 <- usvi_prok_filled.taxa.df %>%
#   dplyr::select(all_of(keep_tax), "asv_id") %>%
#   droplevels %>%
#   dplyr::left_join(., usvi_asv.tbl %>%
#                      dplyr::slice(which(rowSums(.) > 0)) %>%
#                      dplyr::select(rownames(usvi_metab.tbl)) %>%
#                      apply(., 2, relabund) %>%
#                      t() %>%
#                      tibble::as_tibble(rownames = "sample_id") %>%
#                      tidyr::pivot_longer(., cols = !c("sample_id"),
#                                          names_to = "asv_id",
#                                          values_to = "rclr_trans"),
#                    by = join_by("asv_id")) %>%
#   droplevels
#   
# temp_df2 %>%
#   dplyr::filter(grepl("Woese|SAR11|Pseudomonadales", Order)) %>%
#   dplyr::reframe(rclr_options = quantile(rclr_trans, probs = seq(0, 1, 0.25), na.rm = TRUE, names = TRUE), .by = "Order") %>%
#   droplevels
# # # A tibble: 20 × 2
# # Order                    rclr_options
# # <fct>                           <dbl>
# #   1 SAR11 clade                   0      
# # 2 SAR11 clade                   0      
# # 3 SAR11 clade                   0      
# # 4 SAR11 clade                   0.0519 
# # 5 SAR11 clade                   9.68   
# # 6 Pseudomonadales               0      
# # 7 Pseudomonadales               0      
# # 8 Pseudomonadales               0      
# # 9 Pseudomonadales               0      
# # 10 Pseudomonadales              11.8    
# # 11 Woesearchaeales               0      
# # 12 Woesearchaeales               0      
# # 13 Woesearchaeales               0      
# # 14 Woesearchaeales               0      
# # 15 Woesearchaeales               0.0950 
# # 16 Candidatus Woesebacteria      0      
# # 17 Candidatus Woesebacteria      0      
# # 18 Candidatus Woesebacteria      0      
# # 19 Candidatus Woesebacteria      0      
# # 20 Candidatus Woesebacteria      0.00285
}

#for the correlations with rCLR data:

#A. filter for ASVs > 1% cumulative relative abundance: 294 ASVs
# usvi_abund_asvs_idx <- usvi_asv.tbl %>%
#   dplyr::slice(which(rowSums(.) > 0)) %>%
#   dplyr::select(rownames(usvi_metab.tbl)) %>%
#   apply(., 2, relabund) %>%
#   # tibble::as_tibble(.) %>%
#   tibble::as_tibble(., rownames = "asv_id") %>%
#   tibble::column_to_rownames("asv_id") %>%
#   # dplyr::filter(sum(.) > 1) %>%
#   dplyr::slice(which(rowSums(.) >= 1)) %>%
#   as.data.frame(.) %>%
#   droplevels %>%
#   tibble::as_tibble(., rownames = "asv_id") %>%
#   dplyr::distinct(asv_id) %>%
#   tibble::deframe(.)
# # [1] 294

#B. what about for any ASV >1 % in any sample? 57 ASVs
#C. what about for any ASV >= 0.1 % in any sample? 316 ASVs
#D. 

usvi_abund_asvs_idx <- usvi_asv.tbl %>%
# usvi_asv.tbl %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  dplyr::select(rownames(usvi_metab.tbl)) %>%
  apply(., 2, relabund) %>%
  # tibble::as_tibble(.) %>%
  tibble::as_tibble(., rownames = "asv_id") %>%
  dplyr::rowwise(.) %>% dplyr::mutate(rank = across(starts_with("Metab_")) %>% purrr::reduce(sum)) %>%
  # dplyr::mutate(rank = across(starts_with("Metab_")) %>% purrr::reduce(sum), .by = "asv_id") %>%
  # tibble::column_to_rownames("asv_id") %>%
  tidyr::pivot_longer(., cols = starts_with("Metab_"),
                      names_to = "sample_id",
                      values_to = "relabund") %>%
  # dplyr::filter(relabund >= 1) %>%
  dplyr::filter(relabund >= 0.1) %>%
  dplyr::arrange(desc(rank)) %>%
  # dplyr::slice(which(rowSums(.) >= 1)) %>%
  # as.data.frame(.) %>%
  # droplevels %>%
  # tibble::as_tibble(., rownames = "asv_id") %>%
  dplyr::distinct(asv_id) %>%
  tibble::deframe(.)
  nrow(.)



spearman.test.rclr.abund.asv.taxonomy.df <- spearman.test.rclr.asv.taxonomy.df %>% 
  dplyr::filter(asv_id %in% usvi_abund_asvs_idx) %>%
  droplevels
#when using only ASVs that are ever >1% in any the 71 samples,
# 57 ASVs remain with significant and strong (abs(rho) >= 0.8 ) correlations
#when using only ASVs that are ever >= 0.1% in any of the 71 samples,
# 316 ASVs remain with 189 significant and strong (abs(rho) >= 0.8 ) correlations


spearman.test.rclr.abund.asv.taxonomy.df %>%
  dplyr::group_by(taxonomy_string, simpleName) %>%
  dplyr::summarise(num_corr = length(estimate)) %>%
  dplyr::arrange(desc(num_corr))  
# # A tibble: 150 × 3
# # Groups:   taxonomy_string [68]
# taxonomy_string                       simpleName              num_corr
# <chr>                                 <fct>                      <int>
#   1 Bacteroidia; NS9 marine group         homoserine betaine             6
# 2 Alphaproteobacteria; SAR116 clade     taurine                        4
# 3 Alphaproteobacteria; SAR116 clade     kynurenine                     4
# 4 Bacteroidia; NS5 marine group         kynurenine                     4
# 5 Cyanobacteria; Synechococcus CC9902   kynurenine                     4
# 6 Alphaproteobacteria; SAR116 clade     arginine                       3
# 7 Bacteroidia; Cyclobacteriaceae        sn-glycerol 3-phosphate        3
# 8 Bacteroidia; NS4 marine group         kynurenine                     3
# 9 Gammaproteobacteria; OM60(NOR5) clade sn-glycerol 3-phosphate        3
# 10 Thermoplasmata; Marine Group II       sn-glycerol 3-phosphate        3
# # ℹ 140 more rows

#there are still 4 ASVs with Archaeal lineages retained in this correlation analyses 
#Yawzi: ASV_00098 and ASV_00311
#Seagrass: ASV_00292 and ASV_00438

#but they aren't dominating the results anymore

#when using only ASVs that are ever >= 0.1% in any of the 71 samples,
# 310 ASVs remain with 202 significant and strong (abs(rho) >= 0.8 ) correlations
#of these 310 ASVs with 202 strong sig relations, 
#76 ASVs had 2 or more correlations with metabolites
#8 ASVs had 5 or more correlations with metabolites
spearman.test.rclr.abund.asv.taxonomy.df %>%
  count(asv_id, name = "num_correlations") %>%   # Count occurrences
  arrange(desc(num_correlations))%>% 
  distinct()%>%
  dplyr::filter(num_correlations >= 5) %>%
  nrow(.)

#in comparison, before filtering for the more abundant ASVs with strong correlations:
#1597 ASVs had 5 or more correlations with metabolites
spearman.test.rclr.asv.taxonomy.df %>%
  count(asv_id, name = "num_correlations") %>%   # Count occurrences
  arrange(desc(num_correlations))%>% 
  distinct()%>%
  dplyr::filter(num_correlations >= 5) %>%
  nrow(.)

#in comparison, when using the relative abundance of ASVs for correlation input:
#32 ASVs had 2 or more strong, significant correlations with metabolites
#and 12 metabolites had 2 or more strong, sig correlations with ASVs
spearman.test.relabund.asv.taxonomy.df %>%
  # count(asv_id, name = "num_correlations") %>%   # Count occurrences
  count(simpleName, name = "num_correlations") %>%   # Count occurrences
  arrange(desc(num_correlations))%>% 
  distinct()%>%
  dplyr::filter(num_correlations >= 2) %>%
  nrow(.)

# # A tibble: 12 × 2
# simpleName         num_correlations
# <chr>                         <int>
#   1 3'AMP                            49
#  2 GABA                             40
#  3 histidine                        26
#  4 methionine                       25
#  5 spermidine 3                     23
#  6 homoserine betaine               17
#  7 5'AMP                            16
# 8 guanosine                        16
# 9 ciliatine                        12
# 10 putrescine 2                     11
# 11 pantothenic acid                  2
# 12 tryptamine                        2

spearman.test.rclr.abund.asv.taxonomy.df %>%
  dplyr::count(taxonomy_string, simpleName, name = "num_correlations") %>%   # Count occurrences
  dplyr::arrange(desc(num_correlations))%>% 
  dplyr::distinct()%>%
  dplyr::filter(num_correlations >= 2) %>%
  # nrow(.)
  droplevels
# # A tibble: 22 × 3
# taxonomy_string                       simpleName              num_correlations
# <chr>                                 <fct>                              <int>
#   1 Bacteroidia; NS9 marine group         homoserine betaine                     6
# 2 Alphaproteobacteria; SAR116 clade     taurine                                4
# 3 Alphaproteobacteria; SAR116 clade     kynurenine                             4
# 4 Bacteroidia; NS5 marine group         kynurenine                             4
# 5 Cyanobacteria; Synechococcus CC9902   kynurenine                             4
# 6 Alphaproteobacteria; SAR116 clade     arginine                               3
# 7 Bacteroidia; Cyclobacteriaceae        sn-glycerol 3-phosphate                3
# 8 Bacteroidia; NS4 marine group         kynurenine                             3
# 9 Gammaproteobacteria; OM60(NOR5) clade sn-glycerol 3-phosphate                3
# 10 Thermoplasmata; Marine Group II       sn-glycerol 3-phosphate                3
# # ℹ 12 more rows


spearman.test.rclr.abund.asv.taxonomy.df %>%
  dplyr::count(taxonomy_string, simpleName, site, name = "num_correlations") %>%   # Count occurrences
  dplyr::arrange(desc(num_correlations))%>% 
  dplyr::distinct()%>%
  dplyr::filter(num_correlations >= 2) %>%
  # nrow(.)
  droplevels
# # A tibble: 16 × 4
# taxonomy_string                       simpleName              site    num_correlations
# <chr>                                 <fct>                   <fct>              <int>
#   1 Bacteroidia; NS9 marine group         homoserine betaine      LB_sea…                6
# 2 Alphaproteobacteria; SAR116 clade     taurine                 Yawzi                  4
# 3 Alphaproteobacteria; SAR116 clade     kynurenine              Tektite                4
# 4 Cyanobacteria; Synechococcus CC9902   kynurenine              LB_sea…                4
# 5 Alphaproteobacteria; SAR116 clade     arginine                Yawzi                  3
# 6 Bacteroidia; NS4 marine group         kynurenine              LB_sea…                3
# 7 Bacteroidia; NS5 marine group         kynurenine              LB_sea…                3
# 8 Bacteroidia; Cryomorphaceae           homoserine betaine      LB_sea…                2
# 9 Bacteroidia; Cryomorphaceae           kynurenine              LB_sea…                2
# 10 Bacteroidia; Cyclobacteriaceae        sn-glycerol 3-phosphate LB_sea…                2
# 11 Bacteroidia; Flavobacteriaceae        putrescine 2            LB_sea…                2
# 12 Bacteroidia; NS4 marine group         homoserine betaine      LB_sea…                2
# 13 Bacteroidia; NS9 marine group         kynurenine              LB_sea…                2
# 14 Gammaproteobacteria; OM60(NOR5) clade sn-glycerol 3-phosphate Yawzi                  2
# 15 Gammaproteobacteria; SAR86 clade      sn-glycerol 3-phosphate Yawzi                  2
# 16 Thermoplasmata; Marine Group II       sn-glycerol 3-phosphate Yawzi                  2

spearman.test.rclr.abund.asv.taxonomy.df %>%
  # count(asv_id, name = "num_correlations") %>%   # Count occurrences
  count(simpleName, name = "num_correlations") %>%   # Count occurrences
  arrange(desc(num_correlations))%>% 
  distinct()%>%
  #dplyr::filter(num_correlations >= 1) %>% droplevels %>% dplyr::summarise(sum(num_correlations)) #189
  dplyr::filter(num_correlations >= 2) %>% droplevels
  # nrow(.)
# # A tibble: 12 × 2
# simpleName              num_correlations
# <fct>                              <int>
#   1 sn-glycerol 3-phosphate               46
# 2 homoserine betaine                    33
# 3 kynurenine                            30
# 4 cysteate                              15
# 5 putrescine 2                          15
# 6 arginine                              13
# 7 taurine                               12
# 8 5'UMP                                  7
#  9 GABA                                   6
# 10 aspartate                              4
# 11 tryptophan                             2
# 12 3'AMP                                  2


# Keep relative abundance + MH for NMDS -----------------------------------


#export the NMDS coordinates

nmds_asv_relabund_mh.df <- readr::read_delim(list.files(projectpath, "usvi_nmds_asv-", full.names = TRUE)[1],
                                             delim = "\t", col_names = TRUE, show_col_types = FALSE) 

#here are the most abundant ASVs: usvi_abund_asvs_idx (316)
#here are the SDA ASVs: shared_sda_asvs_idx (181)
#here are the SDA and most abundant ASVs: shared_sda_asvs_abund_idx (36)

temp_df1 <- nmds_asv_relabund_mh.df %>%
  dplyr::filter(type == "sample") %>%
  dplyr::mutate(site = factor(site, levels = names(site_lookup)),
                sampling_time = factor(sampling_time, levels = names(sampling_time_lookup))) %>%
  droplevels

temp_df2 <- nmds_asv_relabund_mh.df %>%
  dplyr::filter(type == "species") %>%
  dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"),
                   by = join_by("label" == "asv_id")) %>%
  dplyr::filter(label %in% usvi_abund_asvs_idx) %>%
  dplyr::rename(asv_id = `label`) %>% 
  dplyr::select(asv_id, NMDS1, NMDS2, taxonomy) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df %>%
                     dplyr::select(Domain, Phylum, Class, Order, Family, Genus, asv_id) %>%
                     droplevels,
                   by = join_by(asv_id)) %>%
  droplevels



# temp_df1 <- transformed_usvi_asv_nmds_keep.list[["relabund_horn"]][["points"]]
# 
# temp_df2 <- transformed_usvi_asv_nmds_keep.list[["relabund_horn"]][["species"]] %>%
#   # dplyr::filter(label %in% shared_sda_asvs_idx) %>%
#   dplyr::filter(label %in% shared_sda_asvs_abund_idx) %>%
#   dplyr::rename(asv_id = `label`) %>%
#   dplyr::select(asv_id, NMDS1, NMDS2, taxonomy) %>%
#   dplyr::left_join(., usvi_prok_filled.taxa.df %>%
# dplyr::select(Domain, Phylum, Class, Order, asv_id) %>%
#                      droplevels,
#                    by = join_by(asv_id)) %>%
#   droplevels


# annotation_taxa_order_lookup <- annotation_taxa_colors_list_v2[["color"]][["Order"]] %>%
#   tibble::enframe(name = "taxonomy_string", value = "color") %>%
#   dplyr::right_join(., temp_df2 %>%
#                       dplyr::distinct(asv_id, Domain, Phylum, Class, Order) %>%
#                       dplyr::mutate(taxonomy_string = paste(Domain, Phylum, Class, Order, sep = ";")) %>%
#                       dplyr::distinct(asv_id, taxonomy_string, Order) %>%
#                       droplevels,
#                     by = join_by(taxonomy_string))%>%
#   dplyr::distinct(asv_id, Order) %>% 
#   # dplyr::distinct(asv_id, Order) %>%
#   tibble::deframe(.)
# 
# annotation_taxa_order_colors <- annotation_taxa_colors_list[["color"]][["Order"]] %>%
#   tibble::enframe(name = "taxonomy_string", value = "color") %>%
#   dplyr::right_join(., temp_df2 %>%
#                       dplyr::distinct(asv_id, Domain, Phylum, Class, Order) %>%
#                       dplyr::left_join(., tibble::enframe(usvi_genera_relabel, name = "asv_id", value = "taxonomy"),
#                                        by = join_by(asv_id)) %>%
#                       dplyr::mutate(taxonomy_string = paste(Domain, Phylum, Class, Order, sep = ";")) %>%
#                       dplyr::distinct(asv_id, taxonomy, Order) %>%
#                       droplevels,
#                     by = join_by(taxonomy_string))%>%
#   dplyr::distinct(Order, color) %>%
#   # dplyr::distinct(asv_id, Order, color) %>% dplyr::select(Order, color) %>%
#   tibble::deframe(.)
# # droplevels 

# annotation_taxa_order_lookup <- annotation_taxa_colors_list_v2[["Order"]] %>%
#   tibble::enframe(name = "Order", value = "color") %>%
# dplyr::right_join(., temp_df2 %>%
#                     dplyr::distinct(asv_id, Phylum, Order, Genus) %>%
#                     # dplyr::distinct(asv_id, Order) %>%
#                     dplyr::distinct(asv_id, Order, Genus) %>%
#                     droplevels,
#                   # by = join_by(Order))%>%
#                   by = join_by(Genus))%>%
# dplyr::distinct(asv_id, Order) %>% 
# # dplyr::distinct(asv_id, Order) %>%
# tibble::deframe(.)

annotation_taxa_topasvs_colors <- annotation_taxa_colors_list_v2[["Genus"]] %>%
  tibble::enframe(name = "Genus", value = "color") %>%
  dplyr::right_join(., temp_df2 %>%
                      dplyr::select(-c(NMDS1, NMDS2)) %>%
                      dplyr::distinct(.) %>%
                      # dplyr::distinct(asv_id, Order, Genus) %>%
                      droplevels,
                    by = join_by(Genus))%>%
  dplyr::left_join(., annotation_taxa_colors_list_v2[["Family"]] %>%
                     tibble::enframe(name = "Family", value = "color") %>%
                     droplevels, 
                   by = join_by(Family)) %>%
  dplyr::mutate(color = across(starts_with("color")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(-c(ends_with(".x"), ends_with(".y"))) %>%
  dplyr::left_join(., annotation_taxa_colors_list_v2[["Order"]] %>%
                     tibble::enframe(name = "Order", value = "color") %>%
                     droplevels, 
                   by = join_by(Order)) %>%
  dplyr::mutate(color = across(starts_with("color")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(-c(ends_with(".x"), ends_with(".y"))) %>%
  dplyr::left_join(., annotation_taxa_colors_list_v2[["Class"]] %>%
                     tibble::enframe(name = "Class", value = "color") %>%
                     droplevels, 
                   by = join_by(Class)) %>%
  dplyr::mutate(color = across(starts_with("color")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(-c(ends_with(".x"), ends_with(".y"))) %>%
  dplyr::left_join(., annotation_taxa_colors_list_v2[["Phylum"]] %>%
                     tibble::enframe(name = "Phylum", value = "color") %>%
                     droplevels, 
                   by = join_by(Phylum)) %>%
  dplyr::mutate(color = across(starts_with("color")) %>% purrr::reduce(coalesce)) %>%
  dplyr::select(-c(ends_with(".x"), ends_with(".y"))) %>%
  dplyr::mutate(num_obs = length(asv_id), .by = "Genus")%>%
  droplevels

annotation_taxa_order_colors <- annotation_taxa_topasvs_colors %>%
  split(., f = .$num_obs) %>%
  map(., ~.x %>%
        droplevels %>%
        dplyr::select(Genus, Order, color) %>%
        # dplyr::select(Order, color) %>%
        dplyr::distinct(color, .keep_all = TRUE) %>%
        dplyr::select(Order, color) %>%
        dplyr::distinct(Order, .keep_all = TRUE) %>%
               # tibble::deframe(.))
        tibble::deframe(.) %>% scales::col_saturate(., 50) %>%
        tibble::enframe(., name = "Order", value = "color")) %>%
#   split(., f = .$Genus) %>%
#   map(., ~.x %>%
#         droplevels %>%
#         dplyr::arrange(desc(num_obs)) %>%
#   #       # dplyr::arrange(num_obs) %>%
#         dplyr::select(Genus, Order, color) %>%
#   # dplyr::select(Order, color) %>%
#         dplyr::distinct(color, .keep_all = TRUE) %>%
#         dplyr::select(Genus, color) %>%
# #        # tibble::deframe(.))
#         tibble::deframe(.) %>% scales::col_shift(., 50))
  #       scales::col_shift(., 50) %>%
  #       tibble::enframe(., name = "Order", value = "color")) %>%
  bind_rows(., .id = NULL) %>%
  # dplyr::arrange(desc(num_obs)) %>%
  dplyr::select(Order, color) %>%
  dplyr::distinct(color, .keep_all = TRUE) %>%
  dplyr::distinct(Order, .keep_all = TRUE) %>%
  tibble::deframe(.)

pal.bands(annotation_taxa_order_colors)

#export the results for Brianna
#here are the most abundant ASVs: usvi_abund_asvs_idx (316)
#here are the SDA ASVs: shared_sda_asvs_idx (181)
#here are the SDA and most abundant ASVs: shared_sda_asvs_abund_idx (36)

nmds_asv_relabund_mh_topasvs.df <- temp_df2 %>%
  dplyr::left_join(., tibble::enframe(annotation_taxa_order_colors, name = "Order", value = "color"),
                   relationship = "many-to-many", multiple = "all") %>%
  dplyr::left_join(., tibble::enframe(usvi_abund_asvs_idx, name = "abund_0.01_316", value = "asv_id"),
                   relationship = "many-to-many", multiple = "all", by = join_by(asv_id)) %>%
  dplyr::left_join(., tibble::enframe(shared_sda_asvs_idx, value = "asv_id", name = "sda_asvs_181"),
  relationship = "many-to-many", multiple = "all", by = join_by(asv_id)) %>%
  dplyr::left_join(., tibble::enframe(shared_sda_asvs_abund_idx, value = "asv_id", name = "sda_abund_asvs_36"),
                   relationship = "many-to-many", multiple = "all", by = join_by(asv_id)) %>%
  droplevels

readr::write_delim(nmds_asv_relabund_mh_topasvs.df,
                   paste0(projectpath, "/", "nmds_asv_relabund_mh_topasvs.df", ".tsv"),
                   delim = "\t", col_names = TRUE, num_threads = nthreads)

temp_df3 <- temp_df2 %>%
  dplyr::filter(asv_id %in% head(usvi_abund_asvs_idx, n = 30)) %>%
# temp_df2 <- temp_df2 %>%
  dplyr::left_join(., tibble::enframe(annotation_taxa_order_colors, name = "Order", value = "color"),
                   relationship = "many-to-many", multiple = "all") %>%
  droplevels
y <- nrow(temp_df3)

temp_nmds1 <- (ggplot(data = temp_df1 %>%
                        dplyr::filter(type == "sample") %>%
                        droplevels,
                      aes(x = NMDS1, y = NMDS2))
               + theme_bw()
               + geom_point(aes(shape = sampling_time, fill = site), color = "black",
                            size = 3, stroke = 1, alpha = 1)
               + scale_discrete_manual(aesthetics = c("shape"), 
                                       values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA), drop = TRUE)
               + scale_discrete_manual(aesthetics = c("fill"), 
                                       values = site_colors, labels = site_lookup, breaks = names(site_lookup),drop = TRUE)
               + geom_segment(data = temp_df3,
                              aes(x = 0, y = 0, group = asv_id,
                                  color = Order,
                                  xend = NMDS1, yend = NMDS2),
                              linewidth = 0.5,
                              # color = "limegreen",
                              arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                              # show.legend = TRUE,
                              lineend = "butt", linejoin = "round",
                              alpha = 0.7)
               + ggrepel::geom_text_repel(data = temp_df3,
                                          aes(label = str_wrap(taxonomy, 25),
                                              # color = Order,
                                              x = NMDS1, y = NMDS2),
                                          color = "black",
                                          size = rel(2), min.segment.length = 0,
                                          # point.size = NA, segment.color = NA,
                                          arrow = arrow(type = "open", length = unit(0.02, "npc")), segment.color = "grey",
                                          # colour = "black",
                                          direction = "both", seed = 48105, max.overlaps = 15,  show.legend =FALSE,
                                          force = 2, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"),
                                          fontface = "bold")
               + scale_discrete_manual(aesthetics = c("color"),
                                       values = annotation_taxa_order_colors,
                                       labels = names(annotation_taxa_order_colors),
                                       breaks = names(annotation_taxa_order_colors),
                                       # labels = (annotation_taxa_order_lookup),
                                       # breaks = names(annotation_taxa_order_lookup),
                                       drop = TRUE)
               + scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS2")
               + scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS1")
               + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                       panel.background = element_blank(), 
                       panel.border = element_rect(fill = "NA", colour = "grey30"),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = rel(0.9)),
                       plot.subtitle = element_text(size = rel(0.7)),
                       legend.position = "right",
                       legend.key = element_blank(),
                       legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                       legend.text = element_text(size = 12, colour = "grey30"))
               + guides(shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                             override.aes = list(color = "black", label = "", size = 2)), 
                        fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                            override.aes = list(color = "black", label = "", stroke = 1, shape = 21, size = 2)),
                        
                        color = guide_legend(order = 3, ncol = 1, title = "Bacterial Order", 
                                             override.aes = list(size = 3, label = "->"),
                                             direction = "vertical"),
                        line = "none",
                        size = "none")
               + ggtitle(label = "NMDS using Morisita-Horn on relative abundance",
                         subtitle = paste0("labels: ", {y}, " most abundant ASVs"))
)
print(temp_nmds1)

temp_nmds2 <- (ggplot()
               + geom_segment(data = temp_df3,
                              aes(x = 0, y = 0, group = asv_id,
                                  color = Order,
                                  xend = NMDS1, yend = NMDS2),
                              linewidth = 0.5,
                              # color = "limegreen",
                              arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                              # show.legend = TRUE,
                              lineend = "butt", linejoin = "round",
                              alpha = 0.7)
               + ggrepel::geom_text_repel(data = temp_df3,
                                          aes(label = str_wrap(taxonomy, 25),
                                              # color = Order,
                                              x = NMDS1, y = NMDS2),
                                          color = "black",
                                          size = rel(2), min.segment.length = 0,
                                          # point.size = NA, segment.color = NA,
                                          arrow = arrow(type = "open", length = unit(0.02, "npc")), segment.color = "grey",
                                          # colour = "black",
                                          direction = "both", seed = 48105, max.overlaps = 15,  show.legend =FALSE,
                                          force = 2, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"),
                                          fontface = "bold")
               + scale_discrete_manual(aesthetics = c("color"),
                                       values = annotation_taxa_order_colors,
                                       labels = names(annotation_taxa_order_colors),
                                       breaks = names(annotation_taxa_order_colors),
                                       # labels = (annotation_taxa_order_lookup),
                                       # breaks = names(annotation_taxa_order_lookup),
                                       drop = TRUE)
               + scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS2")
               + scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS1")
               + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                       panel.background = element_blank(), 
                       panel.border = element_rect(fill = "NA", colour = "grey30"),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = rel(0.9)),
                       plot.subtitle = element_text(size = rel(0.7)),
                       legend.position = "right",
                       legend.key = element_blank(),
                       legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                       legend.text = element_text(size = 12, colour = "grey30"))
               + guides(shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                             override.aes = list(color = "black", label = "", size = 2)), 
                        fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                            override.aes = list(color = "black", label = "", stroke = 1, shape = 21, size = 2)),
                        
                        color = guide_legend(order = 3, ncol = 1, title = "Bacterial Order", 
                                             override.aes = list(size = 3, label = "->"),
                                             direction = "vertical"),
                        line = "none",
                        size = "none")
               + ggtitle(label = "NMDS using Morisita-Horn on relative abundance",
                         subtitle = paste0("labels: ", {y}, " most abundant ASVs"))
)
temp_nmds2

ggsave(paste0(projectpath, "/", "usvi_asv_nmds_biplot_top_", y, "-", Sys.Date(), ".svg"),
       temp_nmds1,
       width = 10, height = 10, units = "in")
ggsave(paste0(projectpath, "/", "usvi_asv_only_nmds_biplot_top_", y, "-", Sys.Date(), ".svg"),
       temp_nmds2,
       width = 10, height = 10, units = "in")
