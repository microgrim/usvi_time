# 07_taxa_metabolites_nmds.R

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
library(phyloseq)
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
to_import <- c("usvi_prok_asvs.df", "usvi_prok_asvs.taxa")

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


if(!exists("ps_usvi", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))){
    # ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))
    ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))
  } else {
    cli::cli_alert_warning("Please process the USVI data through Phyloseq.")
    
  }
}
if(file.exists(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"))){
  metabolomics_sample_metadata <- readr::read_delim(paste0(projectpath, "/", "metabolomics_sample_metadata", ".tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
} else {
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

#replace NA in taxonomy with last known level

usvi_prok_filled.taxa.df <- usvi_prok_asvs.taxa %>%
  dplyr::mutate(Phylum = coalesce(Phylum, Domain)) %>%
  dplyr::mutate(Class = coalesce(Class, Phylum)) %>%
  dplyr::mutate(Order = coalesce(Order, Class)) %>%
  dplyr::mutate(Family = coalesce(Family, Order)) %>%
  dplyr::mutate(Genus = coalesce(Genus, Family)) %>%
  dplyr::mutate(Species = coalesce(Species, Genus)) %>%
  dplyr::relocate(asv_id) %>%
  dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                          grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                          .default = Phylum)) %>%
  droplevels


# Prepare color palette ---------------------------------------------------

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
  
  if(!file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
    temp_list <- list(usvi_metabolomics.df, usvi_metabolomics_long.df) %>%
      setNames(., c("usvi_metabolomics.df", "usvi_metabolomics_long.df"))
    readr::write_rds(temp_list, paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
    rm(temp_list)
  }
  
}

# Make a genus-level asv table --------------------------------------------



usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
  dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus) %>%
  dplyr::distinct(Domain, Phylum, Class, Order, Family, Genus, .keep_all = TRUE) %>%
  droplevels

usvi_sw_genus.tbl <- usvi_prok_asvs.df %>%
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

sample_relabel <- metabolomics_sample_metadata %>%
  dplyr::select(sample_id, site, sampling_day, sampling_time) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  droplevels %>%
  tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time), sep = "_", remove = TRUE)  %>%
  tibble::deframe(.)


# NMDS on NPOC/TN profiles ------------------------------------------------

#Brianna determined these were potential outliers in the metabolomes:
# LB seagrass peak: 73, 117
# LB seagrass dawn: 109,110
# Yawzi peak: 67, 69
# Yawzi dawn: none
# Tektite peak: 98
# Tektite dawn: 79

potential_metab_outliers_idx <- c(73, 117, 109, 110, 67, 69, 98, 79) %>%
  paste0("CINAR_BC_", .)
sample_relabel2 <- metabolomics_sample_metadata %>%
  dplyr::select(sample_id, metab_deriv_label, site, sampling_day, sampling_time) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  droplevels %>%
  dplyr::select(sample_id, metab_deriv_label) %>%
  dplyr::mutate(metab_deriv_label = dplyr::case_when(metab_deriv_label %in% potential_metab_outliers_idx ~ paste0(metab_deriv_label, "*"),
                                                     .default = metab_deriv_label)) %>%
  # tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time), sep = "_", remove = TRUE)  %>%
  tibble::deframe(.)
# nutrients_seawater.df <- ps_usvi %>%
#   phyloseq::sample_data(.) %>%
#   tibble::as_tibble(rownames = "sample_id") %>%
#   dplyr::filter(grepl("seawater", sample_type)) %>%
#   droplevels

cinar_npoc_tn.df <- readxl::read_excel(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/CINAR_NPOCandTNdata.2022.04.06.xlsx")) %>%
  dplyr::mutate(doc_label = dplyr::case_when(grepl("CINAR", filename) ~ stringr::str_remove_all(filename, "CINAR_"),
                                             .default = filename)) %>%
  droplevels

usvi_npoc_tn.df <- ps_usvi %>%
  phyloseq::sample_data(.) %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  droplevels %>%
  dplyr::select(sample_id, sample_type, site, sampling_date, sampling_time, sampling_day, doc_label) %>%
  # dplyr::select(sample_id, sample_type, site, sampling_date, sampling_time, sampling_day, contains("metab"), doc_label) %>%
  dplyr::mutate(doc_label = factor(doc_label)) %>%
  dplyr::left_join(., metadata %>%
                     dplyr::select(sample_id, sample_order, sample_order_all),
                   by = join_by(sample_id)) %>%
  dplyr::left_join(., cinar_npoc_tn.df %>%
                     dplyr::filter(grepl("CINAR", filename)) %>%
                     dplyr::select(filename, `NPOC(uM)`, `TN(uM)`, doc_label) %>%
                     dplyr::mutate(doc_label = factor(doc_label)) %>%
                     # dplyr::mutate(doc_label = as.numeric(doc_label)) %>%
                     droplevels,
                   by = join_by(doc_label)) %>%
  dplyr::rename(NPOC = "NPOC(uM)",
                TN = "TN(uM)") %>%
  tidyr::pivot_longer(., cols = !c(contains("sampl"), site, contains("metab"), filename, doc_label),
                      names_to = "nutrient",
                      values_to = "concentration") %>%
  dplyr::left_join(., metabolomics_sample_metadata %>%
                     dplyr::select(sample_id, metab_deriv_label),
                   by = join_by(sample_id), relationship = "many-to-many", multiple = "all") %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(.$sample_id)),
                metab_deriv_label = factor(metab_deriv_label, levels = unique(.$metab_deriv_label))) %>%
  dplyr::mutate(potential_outlier = dplyr::case_when(metab_deriv_label %in% potential_metab_outliers_idx ~ "*",
                                                     .default = NA)) %>%
  droplevels

usvi_npoc_tn.tbl <- usvi_npoc_tn.df %>%
  dplyr::distinct(sample_id, doc_label, nutrient, .keep_all = TRUE) %>%
  tidyr::pivot_wider(., id_cols = c("sample_id", "doc_label"),
                     names_from = "nutrient",
                     values_from = "concentration") %>%
  tidyr::drop_na(.) %>%
  dplyr::right_join(., metabolomics_sample_metadata %>%
                      dplyr::select(sample_id, metab_deriv_label),
                    by = join_by(sample_id), multiple = "all", relationship = "many-to-many") %>%
  tidyr::drop_na(.)
outliers::outlier(usvi_npoc_tn.tbl[,"NPOC"])
outliers::outlier(usvi_npoc_tn.tbl[,"TN"])

min(usvi_npoc_tn.tbl[,"NPOC"])
min(usvi_npoc_tn.tbl[,"TN"])

# temp_df <- usvi_npoc_tn.df %>%
usvi_npoc_tn.df <- usvi_npoc_tn.df %>%
  dplyr::mutate(potential_outlier_y = dplyr::case_when((nutrient == "NPOC" & !is.na(potential_outlier)) ~ signif(min(usvi_npoc_tn.tbl[,"NPOC"]), digits = 0),
                                                     (nutrient == "TN" & !is.na(potential_outlier)) ~ signif(min(usvi_npoc_tn.tbl[,"TN"]), digits = 0),
                                                     .default = NA)) %>%
  droplevels

g_nuts1 <- print(ggplot(data = usvi_npoc_tn.df %>%
               dplyr::filter(!grepl("73|43", metab_deriv_label)) %>%
               droplevels)
      + geom_boxplot(aes(x = interaction(sampling_day, sampling_time), y = concentration, fill = nutrient, group = interaction(sampling_day, sampling_time)))
      + geom_point(aes(x = interaction(sampling_day, sampling_time),
                       y = concentration, fill = nutrient, 
                       group = interaction(sampling_day, sampling_time)), 
                   position = position_jitterdodge(jitter.width = 0.1),
                   shape = 21)
      + theme_bw()
      + facet_grid(nutrient ~ site, drop = TRUE, scales = "free", space = "fixed", 
                   labeller = labeller(site = site_lookup))
      + scale_x_discrete(labels = sample_relabel, name = "Sample time")
      + scale_y_continuous(name = "Concentration (uM)")
      + theme(strip.text.y = element_text(angle = 0),
              axis.text.x = element_text(angle = 90), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
)

g_nuts2 <- print(ggplot(data = usvi_npoc_tn.df %>%
                          dplyr::distinct(sample_id, nutrient, .keep_all = TRUE) %>%
                    # dplyr::filter(!grepl("73|43", metab_deriv_label)) %>%
                    droplevels)
           + geom_boxplot(aes(x = sample_id, y = concentration, fill = sampling_time, group = interaction(sampling_day, sampling_time)))
           + geom_point(aes(x = sample_id,
                            y = concentration, fill = sampling_time, 
                            group = interaction(sampling_day, sampling_time)), 
                        position = position_jitterdodge(jitter.width = 0.1),
                        shape = 21)
           + theme_bw()
           + geom_text(aes(label = potential_outlier, x = sample_id, 
                           y = potential_outlier_y,
                           # y = concentration*0.7,
                           group = interaction(sampling_day, sampling_time)),
                                       stat = "identity", position = position_jitter(width = 0.1, height = 0), size =5,
                                       colour = "black", fontface = "bold")
           + scale_fill_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
           + facet_grid(nutrient ~ site, drop = TRUE, scales = "free", space = "fixed", 
                        labeller = labeller(site = site_lookup))
           + scale_x_discrete(labels = sample_relabel2, name = "Sample")
           + scale_y_continuous(name = "Concentration (uM)")
           + theme(strip.text.y = element_text(angle = 0),
                   axis.text.x = element_text(angle = 90), 
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank())
)

if(!any(grepl("nutrients", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_nutrients-", Sys.Date(), ".png"),
         g_nuts1, 
         width = 10, height = 6, units = "in")
}

# NMDS on metabolomics profiles -------------------------------------------

#here are metabolites where we don't have LODs reported:
usvi_sus_metabolites_idx <- usvi_metabolomics_long.df %>%
  dplyr::arrange(LOD) %>%
  dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
  dplyr::arrange(metabolites) %>%
  dplyr::filter(is.na(LOD)) %>%
  droplevels

#calculate Bray curtis, dispersions from site*sampling_time mean, and plot NMDS 
#first on untransformed metabolite concentrations,
#then on pseudolog (log2(x + 1)) transformed

{
  nmds_metab_df <- usvi_metabolomics.df %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
    dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
    as.matrix(.) %>%
    vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                           droplevels),
                     by = c("label" = "metab_deriv_label")) %>%
    #add labels:
    dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
    # dplyr::mutate(text_label = dplyr::case_when((grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) | label %in% potential_metab_outliers_idx)~ label,
                                                # dplyr::mutate(text_label = dplyr::case_when(grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) ~ label,
                                                .default = NA)) %>%
    droplevels
  nmds_metab_log2_df <- usvi_metabolomics.df %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
    dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
    apply(., 2, function(x) log2(x + 1)) %>%
    as.matrix(.) %>%
    vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                           droplevels),
                     by = c("label" = "metab_deriv_label")) %>%
    #add labels:
    dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
    # dplyr::mutate(text_label = dplyr::case_when(grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) ~ label,
                                                .default = NA)) %>%
    droplevels
  


  ylim <- bind_rows(nmds_metab_df, nmds_metab_log2_df) %>%
  # ylim <- nmds_metab_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS2) %>%
    simplify
  xlim <- bind_rows(nmds_metab_df, nmds_metab_log2_df) %>%
  # xlim <- nmds_metab_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS1) %>%
    simplify
  nmds_lims <- list(ylim = c(round(min(ylim), digits = 4), round(max(ylim), digits = 4)) * 1.1,
                    xlim = c(round(min(xlim), digits = 4), round(max(xlim), digits = 4)) * 1.1)
  # g3 <- (ggplot(data = nmds_metab_df,
  #              aes(x = NMDS1, y = NMDS2))
  # + theme(text = element_text(size=14))
  # + geom_point(data = (nmds_metab_df %>%
  #                                 dplyr::filter(type == "sample") %>%
  #                                 droplevels),
  #                       aes(fill = site, shape = addNA(sampling_time)), color = "black",
  #                       size = 5, stroke = 2, alpha = 1)
  # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
  #                             direction = "both", segment.color = NA, seed = 123,
  #                             force = 1, force_pull = 1,
  #                             point.padding = unit(0.01, "npc"),
  #                             box.padding = unit(0.01, "npc"),
  #                             colour = "black", fontface = "bold")
  #   # + geom_text(aes(label = text_label, x = NMDS1, vjust = "outward", hjust = "center", y = NMDS2),
  #   #             position = "dodge",
  #   #             check_overlap = TRUE, colour = "grey10", fontface = "bold")
  # + scale_shape_manual(values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
  # + scale_fill_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
  # + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                  panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                  panel.grid = element_blank(),
  #                  legend.position = "bottom",
  #                  legend.key = element_blank(),
  #                  legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                  legend.text = element_text(size = 12, colour = "grey30"))
  # + guides(color = "none",
  #                   fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
  #                                       override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                   shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                        override.aes = list(color = "black", stroke = 1, size = 2)))
  # + coord_cartesian(ylim = nmds_lims[["ylim"]],
  #                            xlim = nmds_lims[["xlim"]],
  #                            expand = TRUE)
  # )
  # g3 + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                     labeller = labeller(site = site_lookup))  + coord_cartesian(ylim = nmds_lims[["ylim"]],
  #                                                                                 xlim = nmds_lims[["xlim"]],
  #                                                                                 expand = TRUE)
  # print(g3)
  
  g3_ell <- (ggplot(data = nmds_metab_df,
                aes(x = NMDS1, y = NMDS2))
         + theme(text = element_text(size=14))
         + geom_point(data = (nmds_metab_df %>%
                                dplyr::filter(type == "sample") %>%
                                droplevels),
                      aes(shape = site, fill = addNA(sampling_time)), color = "black",
                      size = 5, stroke = 2, alpha = 1)
         + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
         # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
                                     direction = "both", segment.color = NA, seed = 123,
                                     force = 1, force_pull = 1,
                                     point.padding = unit(0.01, "npc"),
                                     box.padding = unit(0.01, "npc"),
                                     colour = "black", fontface = "bold")
         + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                        type = "norm", 
                        level = 0.95,
                        geom = "polygon", alpha = 0.05, show.legend = FALSE)
         + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
         + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
         + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
         + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                 panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                 panel.grid = element_blank(),
                 legend.position = "bottom",
                 legend.key = element_blank(),
                 legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                 legend.text = element_text(size = 12, colour = "grey30"))
         + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
                      labeller = labeller(site = site_lookup))
         + guides(color = "none",
                  fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                                      override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                  shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                       override.aes = list(color = "black", stroke = 1, size = 2)))
         + coord_cartesian(ylim = nmds_lims[["ylim"]],
                           xlim = nmds_lims[["xlim"]],
                           expand = TRUE)
  )
  
  g3_log2_ell <- (ggplot(data = nmds_metab_log2_df,
                    aes(x = NMDS1, y = NMDS2))
             + theme(text = element_text(size=14))
             + geom_point(data = (nmds_metab_log2_df %>%
                                    dplyr::filter(type == "sample") %>%
                                    droplevels),
                          aes(shape = site, fill = addNA(sampling_time)), color = "black",
                          size = 5, stroke = 2, alpha = 1)
             + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
                                             # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
                                         direction = "both", segment.color = NA, seed = 123,
                                         force = 1, force_pull = 1,
                                         point.padding = unit(0.01, "npc"),
                                         box.padding = unit(0.01, "npc"),
                                         colour = "black", fontface = "bold")
             + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                            type = "norm", 
                            level = 0.95,
                            geom = "polygon", alpha = 0.05, show.legend = FALSE)
             + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
             + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
             + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
             + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                     panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                     panel.grid = element_blank(),
                     legend.position = "bottom",
                     legend.key = element_blank(),
                     legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                     legend.text = element_text(size = 12, colour = "grey30"))
             + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
                          labeller = labeller(site = site_lookup))
             + guides(color = "none",
                      fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                                          override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                      shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                           override.aes = list(color = "black", stroke = 1, size = 2)))
             + coord_cartesian(ylim = nmds_lims[["ylim"]],
                               xlim = nmds_lims[["xlim"]],
                               expand = TRUE)
  )
  print(g3_log2_ell)
  
  # if(!any(grepl("nmds", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
  #   ggsave(paste0(projectpath, "/", "usvi_metabolomics_nmds-", Sys.Date(), ".png"),
  #          g3,
  #          width = 10, height = 8, units = "in")
  # }
}
  

#plot dispersion from centroid
usvi_metab_log2.tbl <- usvi_metabolomics.df %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  as.matrix(.) 
usvi_metab.tbl <- usvi_metabolomics.df %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  as.matrix(.) 
usvi_metab.meta <- metabolomics_sample_metadata %>%
  dplyr::filter(metab_deriv_label %in% rownames(usvi_metab.tbl)) %>%
  dplyr::select(metab_deriv_label, sample_id, sampling_time, sampling_day, site) %>%
  dplyr::mutate(grouping = interaction(site, sampling_time)) %>%
  # dplyr::select(!c(contains("label"), contains("dna_"))) %>%
  tibble::column_to_rownames(., var = "metab_deriv_label") %>%
  droplevels

dist_usvi_metab.d <- vegan::vegdist(usvi_metab.tbl, method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_metab.df <- vegan::betadisper(dist_usvi_metab.d, type = "median",
                                       usvi_metab.meta$grouping) %>%
                                        # interaction(usvi_metab.meta$site, usvi_metab.meta$sampling_time)) %>%
  purrr::pluck("distances") %>%
  tibble::enframe(value = "dissimilarity", name = "metab_deriv_label") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                         droplevels),
                   by = c("metab_deriv_label" = "metab_deriv_label")) %>%
  dplyr::mutate(potential_outlier = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", metab_deriv_label) | metab_deriv_label %in% potential_metab_outliers_idx)~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),  
  # dplyr::mutate(potential_outlier = dplyr::case_when((grepl("73|43|44|42|70|74|69|67|71|101|108|105", metab_deriv_label) | metab_deriv_label %in% potential_metab_outliers_idx)~ metab_deriv_label,
  # dplyr::mutate(potential_outlier = dplyr::case_when(metab_deriv_label %in% potential_metab_outliers_idx ~ metab_deriv_label,
                                                     .default = NA)) %>%
  dplyr::mutate(not_outlier = dplyr::case_when(is.na(potential_outlier) ~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),  
                                               .default = NA)) %>%
  droplevels
dist_usvi_metab_log2.d <- vegan::vegdist(usvi_metab_log2.tbl, method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_metab_log2.df <- vegan::betadisper(dist_usvi_metab_log2.d, 
                                             type = "median", #this is default
                                             # type = "centroid",
                                             usvi_metab.meta$grouping) %>%
                                             # interaction(usvi_metab.meta$site, usvi_metab.meta$sampling_time)) %>%
  purrr::pluck("distances") %>%
  tibble::enframe(value = "dissimilarity", name = "metab_deriv_label") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                         droplevels),
                   by = c("metab_deriv_label" = "metab_deriv_label")) %>%
  dplyr::mutate(potential_outlier = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", metab_deriv_label) | metab_deriv_label %in% potential_metab_outliers_idx)~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),  
                                                     .default = NA)) %>%
  dplyr::mutate(not_outlier = dplyr::case_when(is.na(potential_outlier) ~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),  
                                                     .default = NA)) %>%
  droplevels
g3_disp <- print(
  ggplot(data = dist_usvi_metab.df)
  + theme_bw()
  + geom_boxplot(aes(x = site, y = dissimilarity, 
                     group = interaction(site, sampling_time)), 
                 color = "black",
                 position = position_dodge2(padding = 0.2, preserve = "single"),
                 show.legend = FALSE, outliers = FALSE)
  + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
               position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
               # position = position_jitter(width = 0.75, seed = 48105, height = 0),
               alpha = 1.0, size = 3)
  # + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "Bray-Curtis dissimilarity")
  + scale_x_discrete(labels = site_lookup, name = "Site")
  + ggrepel::geom_text_repel(aes(label = not_outlier, x = site, 
                  y = dissimilarity, fill = sampling_time,
                  group = interaction(site, sampling_time)), size =4, seed = 48105, 
              stat = "identity",  position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
              colour = "grey", fontface = "bold")
  + ggrepel::geom_text_repel(aes(label = potential_outlier, x = site,
                                 y = dissimilarity, fill = sampling_time,
                                 group = interaction(site, sampling_time)), size =4,
                             stat = "identity",  position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
                             # min.segment.length = 0.1, segment.colour = "black", arrow = arrow(length = unit(0.02, "npc")), direction = "y", box.padding = 1,
                             seed = 48105, 
                             colour = "black", fontface = "bold")
  + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
          panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
          legend.text = element_text(size = 12, colour = "grey30"))
  + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
               labeller = labeller(site = site_lookup))  
  + guides(color = "none",
           fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                               override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
           shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                override.aes = list(color = "black", stroke = 1, size = 2)))
)


g3_log2_disp <- print(
  ggplot(data = dist_usvi_metab_log2.df)
  + theme_bw()
  + geom_boxplot(aes(x = site, y = dissimilarity, 
                     group = interaction(site, sampling_time)), 
                 color = "black",
                 position = position_dodge2(padding = 0.2, preserve = "single"),
                 show.legend = FALSE, outliers = FALSE)
  + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
               position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
               # position = position_jitter(width = 0.75, seed = 48105, height = 0),
               alpha = 1.0, size = 3)
  # + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Bray-Curtis dissimilarity")
  + scale_x_discrete(labels = site_lookup, name = "Site")
  + ggrepel::geom_text_repel(aes(label = not_outlier, x = site, 
                                 y = dissimilarity, fill = sampling_time,
                                 group = interaction(site, sampling_time)), size =4, seed = 48105,
                             # stat = "identity",  
                             position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
                             colour = "grey", fontface = "bold")
  + ggrepel::geom_text_repel(aes(label = potential_outlier, x = site,
                                 y = dissimilarity, fill = sampling_time,
                                 group = interaction(site, sampling_time)), size =4, seed = 48105,
                             # stat = "identity",  
                             position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
                             # min.segment.length = 0.1, segment.colour = "black", arrow = arrow(length = unit(0.02, "npc")), seed = 48105, direction = "y", box.padding = 1,
                             colour = "black", fontface = "bold")
  + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
          panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
          legend.text = element_text(size = 12, colour = "grey30"))
  + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
               labeller = labeller(site = site_lookup))  
  + guides(color = "none",
           fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                               override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
           shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                override.aes = list(color = "black", stroke = 1, size = 2)))
)


# in NMDS, CINAR_BC_81A clusters with Tektite and Yawzi samples
# and CINAR_BC_81B clusters with LB_seagrass samples


# in NMDS, CINAR_BC_73 (a peak photo LB_seagrass sample) also clusters with Tektite and Yawzi samples
# and CINAR_BC_43 (a peak photo Tektite sample) clusters closer to LB_seagrass samples than with other Tektite peak photo samples

# but in the long dataset, CINAR_BC_73 is associated with LB_seagrass Metab_280
# and CINAR_BC_43 is associated with Tektite Metab_236
#does CINAR_BC_43 look like the BC replicates for Metab_227 (CINAR_BC_44) and Metab_224 (CINAR_BC_42)?
# does CINAR_BC_73 look like BC replicates for Metab_286 (CINAR_BC_70) and Metab_289 (CINAR_BC_74)?


nmds_metab_df %>%
  dplyr::arrange(sampling_day) %>%
  # dplyr::filter(grepl("Yawzi", site))
  # dplyr::filter(grepl("LB", site))
  dplyr::filter(grepl("Tektite", site))
  # dplyr::filter(grepl("42|43|44|70|73|74", label))

#LB_seagrass: CINAR_BC_70 and CINAR_BC_74 have NMDS1 < 0 and NMDS2 < 0, whereas CINAR_BC_73 NMDS1 > 0
#Tektite: CINAR_BC_43 has NMDS2 < 0, whereas CINAR_BC_42 and CINAR_BC_44 have NMDS2 > 0

# CINAR_BC_73 is supposed to be LB peak_photo Metab_280 
# Metab_280, Metab_286, and Metab_289 are in the same quadrant together

# CINAR_BC_43 is supposed to be Tektite peak_photo Metab_236 but looks like 
# Metab_227 and Metab_236 are at least in the same quadrant

#another sample, CINAR_BC_69 (Yawzi peak photo corresponding to Metab_274), is suspiciously clustered away from its other replicates (CINAR_BC_67, Metab_271) (CINAR_BC_71, Metab_269)
#another sample, CINAR_BC_105 (Tektite dawn, Metab_321) clusters away from other dawn samples (Metab_321, Metab_315, CINAR_BC_101, CINAR_BC_108)

# #what if we arranged by NMDS1 and NMDS2 and see where BC_73 and BC_43 align?
# temp_df <- nmds_metab_df %>% 
#   dplyr::arrange(NMDS2, NMDS1)



#run NMDS on just sites' metabolomes
{
  
temp_df1 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                  dplyr::filter(grepl("LB", site)) %>%
                  droplevels %>%
                  dplyr::select(metab_deriv_label) %>%
                    tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score) 
temp_df2 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                          dplyr::filter(grepl("Yawzi", site)) %>%
                                          droplevels %>%
                                          dplyr::select(metab_deriv_label) %>%
                                          tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score)

temp_df3 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                          dplyr::filter(grepl("Tektite", site)) %>%
                                          droplevels %>%
                                          dplyr::select(metab_deriv_label) %>%
                                          tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  apply(., 2, function(x) log2(x + 1)) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score)

temp_nmds_metab_log2_df <- bind_rows(temp_df1, temp_df2, temp_df3) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                         droplevels),
                   by = c("label" = "metab_deriv_label")) %>%
  #add labels:
  dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
                                              # dplyr::mutate(text_label = dplyr::case_when(grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) ~ label,
                                              .default = NA)) %>%
  droplevels
}

#now do the untransformed:
{
temp_df1 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                          dplyr::filter(grepl("LB", site)) %>%
                                          droplevels %>%
                                          dplyr::select(metab_deriv_label) %>%
                                          tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score) 
temp_df2 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                          dplyr::filter(grepl("Yawzi", site)) %>%
                                          droplevels %>%
                                          dplyr::select(metab_deriv_label) %>%
                                          tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score)

temp_df3 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                          dplyr::filter(grepl("Tektite", site)) %>%
                                          droplevels %>%
                                          dplyr::select(metab_deriv_label) %>%
                                          tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score)

temp_nmds_metab_df <- bind_rows(temp_df1, temp_df2, temp_df3) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                         droplevels),
                   by = c("label" = "metab_deriv_label")) %>%
  #add labels:
  dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
                                              # dplyr::mutate(text_label = dplyr::case_when(grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) ~ label,
                                              .default = NA)) %>%
  droplevels
}

g3_log2_lb_ell <- (ggplot(data = temp_nmds_metab_log2_df,
                       aes(x = NMDS1, y = NMDS2))
                + theme(text = element_text(size=14))
                + geom_point(data = (temp_nmds_metab_log2_df %>%
                                       dplyr::filter(type == "sample") %>%
                                       droplevels),
                             aes(shape = site, fill = addNA(sampling_time)), color = "black",
                             size = 5, stroke = 2, alpha = 1)
                + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
                                           # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
                                           direction = "both", segment.color = NA, seed = 123,
                                           force = 1, force_pull = 1,
                                           point.padding = unit(0.01, "npc"),
                                           box.padding = unit(0.01, "npc"),
                                           colour = "black", fontface = "bold")
                + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                               type = "t", 
                               level = 0.95,
                               geom = "polygon", alpha = 0.05, show.legend = FALSE)
                + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
                + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
                + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
                + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                        panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                        panel.grid = element_blank(),
                        legend.position = "bottom",
                        legend.key = element_blank(),
                        legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                        legend.text = element_text(size = 12, colour = "grey30"))
                + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
                             labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
                + guides(color = "none",
                         fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                                             override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                         shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                              override.aes = list(color = "black", stroke = 1, size = 2)))
                # + coord_cartesian(ylim = nmds_lims[["ylim"]],
                #                   xlim = nmds_lims[["xlim"]],
                #                   expand = TRUE)
)
print(g3_log2_lb_ell)

g3_lb_ell <- (ggplot(data = temp_nmds_metab_df,
                          aes(x = NMDS1, y = NMDS2))
                   + theme(text = element_text(size=14))
                   + geom_point(data = (temp_nmds_metab_df %>%
                                          dplyr::filter(type == "sample") %>%
                                          droplevels),
                                aes(shape = site, fill = addNA(sampling_time)), color = "black",
                                size = 5, stroke = 2, alpha = 1)
                   + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
                                              # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
                                              direction = "both", segment.color = NA, seed = 123,
                                              force = 1, force_pull = 1,
                                              point.padding = unit(0.01, "npc"),
                                              box.padding = unit(0.01, "npc"),
                                              colour = "black", fontface = "bold")
                   + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                                  type = "t", 
                                  level = 0.95,
                                  geom = "polygon", alpha = 0.05, show.legend = FALSE)
                   + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
                   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
                   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
                   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                           panel.grid = element_blank(),
                           legend.position = "bottom",
                           legend.key = element_blank(),
                           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                           legend.text = element_text(size = 12, colour = "grey30"))
                   + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
                                labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
                   + guides(color = "none",
                            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                                 override.aes = list(color = "black", stroke = 1, size = 2)))
                   # + coord_cartesian(ylim = nmds_lims[["ylim"]],
                   #                   xlim = nmds_lims[["xlim"]],
                   #                   expand = TRUE)
)
print(g3_lb_ell)


# g3 <- (g3_disp + ggtitle("Community dissimilarity distance from site average") + guides(fill = "none", shape = "none")) / (g3_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed metabolite concentrations")
g3 <- (g3_disp + ggtitle("Community dissimilarity distance from site average") + guides(fill = "none", shape = "none")) / (g3_lb_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed metabolite concentrations")
# g3_log2 <- (g3_log2_disp + ggtitle("Community dissimilarity distance from site average" ) + guides(fill = "none", shape = "none")) / (g3_log2_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed metabolite concentrations")
g3_log2 <- (g3_log2_disp + ggtitle("Community dissimilarity distance from site average" ) + guides(fill = "none", shape = "none")) / (g3_log2_lb_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed metabolite concentrations")

print(g3)
print(g3_log2)

if(!any(grepl("bc_nmds", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_metabolomics_bc_nmds-", Sys.Date(), ".png"),
         g3,
         width = 14, height = 12, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_metabolomics_bc_nmds_log2-", Sys.Date(), ".png"),
         g3_log2,
         width = 14, height = 12, units = "in")
}


#which metabolites are driving the differences in each site?
temp_df1 <- usvi_metabolomics.df %>%
  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                          dplyr::filter(grepl("LB", site)) %>%
                                          droplevels %>%
                                          dplyr::select(metab_deriv_label) %>%
                                          tibble::deframe(.))) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
  as.matrix(.) %>%
  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
  ggplot2::fortify(.) %>%
  as.data.frame %>%
  dplyr::mutate(type = dplyr::case_when(score == "sites" ~ "sample",
                                        .default = "metabolite")) %>%
  dplyr::select(-score) %>%
  dplyr::mutate(site = "LB_seagrass") %>%
  droplevels %>%
  bind_rows(., (usvi_metabolomics.df %>%
                  dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                                          dplyr::filter(grepl("Tektite", site)) %>%
                                                          droplevels %>%
                                                          dplyr::select(metab_deriv_label) %>%
                                                          tibble::deframe(.))) %>%
                  tibble::column_to_rownames(var = "metab_deriv_label") %>%
                  dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
                  as.matrix(.) %>%
                  vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
                  ggplot2::fortify(.) %>%
                  as.data.frame %>%
                  dplyr::mutate(type = dplyr::case_when(score == "sites" ~ "sample",
                                                        .default = "metabolite")) %>%
                  dplyr::select(-score) %>%
                  dplyr::mutate(site = "Tektite") %>%
                  droplevels)) %>%
  bind_rows(., (usvi_metabolomics.df %>%
              dplyr::filter(metab_deriv_label %in% (metabolomics_sample_metadata %>%
                                                      dplyr::filter(grepl("Yawzi", site)) %>%
                                                      droplevels %>%
                                                      dplyr::select(metab_deriv_label) %>%
                                                      tibble::deframe(.))) %>%
              tibble::column_to_rownames(var = "metab_deriv_label") %>%
              dplyr::select(!usvi_sus_metabolites_idx[["metabolites"]]) %>%
              as.matrix(.) %>%
              vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
              ggplot2::fortify(.) %>%
              as.data.frame %>%
              dplyr::mutate(type = dplyr::case_when(score == "sites" ~ "sample",
                                                    .default = "metabolite")) %>%
              dplyr::select(-score) %>%
              dplyr::mutate(site = "Yawzi") %>%
              droplevels)) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sampling_time, sampling_day, metab_deriv_label) %>%
                         droplevels),
                   by = c("label" = "metab_deriv_label")) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(wt = dplyr::case_when(type == "metabolite" ~ round(10^(max(abs(NMDS1), abs(NMDS2))), digits = 0)),
                .default = NA) %>%
  #add labels:
  dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
                                              # dplyr::mutate(text_label = dplyr::case_when(grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) ~ label,
                                              .default = NA)) %>%
  droplevels
temp_g <- (ggplot(data = temp_df1,
                     aes(x = NMDS1, y = NMDS2))
              + theme(text = element_text(size=14))
           + geom_segment(data = (temp_df1 %>%
                                    dplyr::filter(type == "metabolite") %>%
                                    droplevels),
                          aes(x = 0, y = 0,
                              xend = NMDS1,
                              yend = NMDS2),
                          arrow = grid::arrow(angle = 20, length = unit(0.05, "npc"),
                                              ends = "last", type = "closed"),
                          linewidth = 0.5, alpha = 0.5, colour = "grey30")
           + ggrepel::geom_text_repel(data = (temp_df1 %>%
                                                dplyr::filter(type == "metabolite") %>%
                                                droplevels),
                                      aes(label = label, x = NMDS1, y = NMDS2, size = wt),
                                      # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
                                      direction = "both", segment.color = NA, seed = 123,
                                      force = 1, force_pull = 1, 
                                      point.padding = unit(0.01, "npc"),
                                      box.padding = unit(0.01, "npc"),
                                      colour = "grey", fontface = "bold")
           + geom_point(data = (temp_df1 %>%
                                  dplyr::filter(type == "sample") %>%
                                  droplevels),
                        aes(shape = site, fill = addNA(sampling_time)), color = "black",
                        size = 5, stroke = 2, alpha = 1)
           
           + ggrepel::geom_text_repel(data = (temp_df1 %>%
                                                dplyr::filter(type == "sample") %>%
                                                droplevels),
                                      aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
                                      # + ggrepel::geom_label_repel(aes(label = text_label, x = NMDS1, y = NMDS2),
                                      direction = "both", segment.color = NA, seed = 123,
                                      force = 1, force_pull = 1,
                                      point.padding = unit(0.01, "npc"),
                                      box.padding = unit(0.01, "npc"),
                                      colour = "black", fontface = "bold")
           + stat_ellipse(data = (temp_df1 %>%
                                    dplyr::filter(type == "sample") %>%
                                    droplevels),
                          aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                          type = "t", 
                          level = 0.95,
                          geom = "polygon", alpha = 0.05, show.legend = FALSE)
              + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
              + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
              + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
              + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                      panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                      panel.grid = element_blank(),
                      legend.position = "bottom",
                      legend.key = element_blank(),
                      legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                      legend.text = element_text(size = 12, colour = "grey30"))
              + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed",
                           labeller = labeller(site = site_lookup))
              + guides(color = "none",
                       size = "none",
                       fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                                           override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                       shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
                                            override.aes = list(color = "black", stroke = 1, size = 2)))
              # + coord_cartesian(ylim = nmds_lims[["ylim"]],
              #                   xlim = nmds_lims[["xlim"]],
              #                   expand = TRUE)
)
print(temp_g)

if(!any(grepl("nmds_arrows", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_metabolomics_nmds_arrows-", Sys.Date(), ".png"),
         temp_g,
         width = 16, height = 10, units = "in")
}
# NMDS on asvs ------------------------------------------------------------



#see if in the DNA based NMDS whether Metab_280 and Metab_236 look like their replicates
#LB_seagrass peak photo day 3: Metab_280, Metab_286, Metab_289
#Tektite peak photo day day 2: Metab_236, Metab_227, Metab_224

{
usvi_asv.tbl <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
  as.data.frame %>%
  apply(., 2, relabund) %>% 
  as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = -c("asv_id"),
                      names_to = "sample",
                      values_to = "abundance") %>%
  dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
                                  log2(abundance+1), #log transform abundance (with +1 pseudocount)
                                  0)) %>%
  tidyr::pivot_wider(., id_cols = "sample",
                     # values_from = "abundance",
                     values_from = "logabund",
                     names_from = "asv_id") %>%
  tibble::column_to_rownames(var = "sample")

meta.seawater <- ps_usvi %>%
  phyloseq::sample_data(.) %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  dplyr::filter(sample_id %in% rownames(usvi_asv.tbl)) %>%
  dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
  # dplyr::select(sample_id, sampling_time, sampling_day, site, contains("fcm")) %>%
  dplyr::select(!c(contains("label"), contains("dna_"))) %>%
  tibble::column_to_rownames(., var = "sample_id") %>%
  droplevels

nmds.asv <- vegan::metaMDS(usvi_asv.tbl, distance = "bray", 
                           autotransform = TRUE)

env.sw <- vegan::envfit(nmds.asv, meta.seawater, permutations = 999, 
                        na.rm = TRUE)

nmds_asv_df <- (ggplot2::fortify(env.sw)) %>%
  bind_rows(., (ggplot2::fortify(nmds.asv) %>%
                  dplyr::filter(score == "sites") %>%
                  dplyr::mutate(type = "sample") %>%
                  dplyr::select(-score))) %>%
  bind_rows(., (ggplot2::fortify(nmds.asv) %>%
                  dplyr::filter(score == "species") %>%
                  dplyr::mutate(type = "ASV") %>%
                  dplyr::select(-score))) %>%
  dplyr::left_join(., (metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, replicate) %>%
                         droplevels),
                   by = c("label" = "sample_id")) %>%
  dplyr::mutate(text_label = dplyr::case_when(type == "sample" & grepl("280|286|289|236|227|224|274|271|269|321|315|321", label) ~ label,
                                              .default = NA)) %>%
  droplevels

}

nmds_asv_df %>%
  dplyr::arrange(sampling_day) %>%
  dplyr::filter(grepl("peak_photo", sampling_time)) %>%
  # dplyr::filter(grepl("280|286|289|236|227|224", label))
  # dplyr::filter(grepl("LB", site))
  dplyr::filter(grepl("Tektite", site))
# dplyr::filter(grepl("Yawzi", site))

#if anything, Metab_227 and Metab_236 are at least in the same quadrant, and Metab_224 is not proximal to them
# and Metab_280, Metab_286, and Metab_289 are in the same quadrant together
# also, LB_seagrass peak_photo samples Metab_273 (day 4), and the three from day 5 (Metab_326, Metab_329, and Metab_332) are clustering away from other Lb samples and closer to Yawzi and Tektite samples

# in the DNA extraction log, Metab_273 and Metab_270 after extraction were originally accidentally mixed
# Metab_270 is from Tektite day 4 peak photo, and it looks like the other replicates from that same site and time (Metab_261, Metab_264)



{
  ylim <- nmds_asv_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS2) %>%
    simplify
  xlim <- nmds_asv_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS1) %>%
    simplify
  nmds_lims <- list(ylim = c(round(min(ylim), digits = 1), round(max(ylim), digits = 1)) * 1.1,
                    xlim = c(round(min(xlim), digits = 1), round(max(xlim), digits = 1)) * 1.1)
  g1 <- (ggplot(data = nmds_asv_df,
                aes(x = NMDS1, y = NMDS2))
         + theme(text = element_text(size=14))
         + geom_point(data = (nmds_asv_df %>%
                                dplyr::filter(type == "sample") %>%
                                droplevels),
                      aes(fill = site, shape = sampling_time), color = "black",
                      size = 5, stroke = 2, alpha = 1)
         + geom_text(aes(label = text_label, x = NMDS1, vjust = "outward", hjust = "center", y = NMDS2),
                     position = "dodge",
                     check_overlap = TRUE, colour = "grey10", fontface = "bold")
         + scale_shape_manual(values = c(22, 21), labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
         + scale_fill_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
         + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                 panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                 panel.grid = element_blank(),
                 legend.position = "bottom",
                 legend.key = element_blank(),
                 legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                 legend.text = element_text(size = 12, colour = "grey30"))
         + guides(color = "none",
                  fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                      override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                  shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                       override.aes = list(color = "black", stroke = 1, size = 2)))
         + coord_cartesian(ylim = nmds_lims[["ylim"]],
                           xlim = nmds_lims[["xlim"]],
                           expand = TRUE)
  )
  print(g1)
  
  if(!any(grepl("nmds", list.files(projectpath, pattern = "usvi_asv.*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_asv_nmds-", Sys.Date(), ".png"),
           g1,
           width = 10, height = 8, units = "in")
  }
  }


