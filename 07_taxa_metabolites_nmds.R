# 07_taxa_metabolites_nmds.R

#evaluate metabolome and microbiomes via NMDS

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


# g1 <- P_nmds(dataset = df %>%
#                                dplyr::filter(grepl("LB_seagrass", site)),
#                              x = taxonomy,
#                              y = relabund, 
#                              z = taxonomy,
#                              facet_form = "sampling_time ~ metric")

P_nmds <- function(dataset, facet_form) {
  
  #x: NMDS1
  #y: NMDS2
  #z: fill
  #w: shape
  
  if(exists("facet_form", inherits = TRUE)){
    if(grepl("\\.", facet_form)){
      facet_form <- deparse(substitute(facet_form))
    }
  } else {
    facet_form <- deparse(substitute(z)) %>%
      paste0(".~", .)
  }
  
  g <- ggplot(data = dataset, aes(y = NMDS2,
                                  x = NMDS1))
  g <- g + theme(text = element_text(size=14))
  g <- g + geom_point(aes(shape = site, fill = addNA(sampling_time), group = interaction(site, sampling_time)), color = "black", size = 5, stroke = 2, alpha = 1)
  g <- g + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                        type = "norm", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
  g <- g + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  g <- g + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  g <- g + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  g <- g + scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS2")
  g <- g + scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS1")
  
  g <- g + theme_bw() 
  g <- g + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                 panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                 panel.grid = element_blank(),
                 legend.position = "bottom",
                 legend.key = element_blank(),
                 legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                 legend.text = element_text(size = 12, colour = "grey30"))
  # g <- g + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #              labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  g <- g + guides(color = "none",
           fill = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                               override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
           shape = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                override.aes = list(color = "black", stroke = 1, size = 2)))
  g <- g + facet_grid(drop = TRUE, scales = "free", space = "fixed", 
                      labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup),
                      eval(facet_form))
  g
}


P_betadisp <- function(dataset, 
                       # facet_form, 
                       metric) {
  #x: site
  #y: dissimilarity
  #z: fill, sampling_time
  #w: shape, site
  
  # if(exists("facet_form", inherits = TRUE)){
  #   # if(grepl("\\.", facet_form)){
  #     facet_form <- deparse(substitute(facet_form))
  #   # }
  # } else {
  #   facet_form <- ". ~ site"
  # }
  if(exists("metric", inherits = TRUE)){
    metric <- paste0(metric, " dissimilarity")
  } else {
    metric <- "Dissimilarity"
  }
  
  g <- (ggplot(data = dataset)
  + theme_bw()
  + geom_boxplot(aes(x = site, y = dissimilarity, 
                     group = interaction(site, sampling_time)), 
                 color = "black",
                 position = position_dodge2(padding = 0.2, preserve = "single"),
                 show.legend = FALSE, outliers = FALSE)
  + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
               position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
               alpha = 1.0, size = 3)
  + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = paste0(metric))
  + scale_x_discrete(labels = site_lookup, name = "Site")
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
           fill = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                               override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
           shape = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                override.aes = list(color = "black", stroke = 1, size = 2)))
  )
g
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


if(!exists("ps_usvi", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))){
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


# Prepare color palette ---------------------------------------------------

site_lookup <- data.frame(site = c("LB_seagrass", "Tektite", "Yawzi", "control_extraction", "control_pcr", "control_seq"),
                          label = c("Lameshur Bay seagrass", "Tektite Reef", "Yawzi Reef",
                                    "Control (DNA Extraction)", "Control (PCR)", "Control (Sequencing)")) %>%
  tibble::deframe(.)
site_colors <- pals::kelly(22)[6:(5+length(site_lookup))] %>%
  setNames(., names(site_lookup))
sampling_time_lookup <- data.frame(sampling_time = c("dawn", "peak_photo"),
                                   label = c("Dawn", "Afternoon")) %>%
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


# Plot metabolite profiles ------------------------------------------------

#here are metabolites where we don't have LODs reported:
usvi_sus_metabolites_idx <- usvi_metabolomics_long.df %>%
  dplyr::arrange(LOD) %>%
  dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
  dplyr::arrange(metabolites) %>%
  dplyr::filter(is.na(LOD)) %>%
  droplevels

#spoiler: CINAR_BC_73 is suspicious, so drop it.

#also, there are some entries that are dummy values: they are reported as 1/10 the LOD
# the measurement should be removed

usvi_metab_summary.df <- usvi_metabolomics.df %>%
  dplyr::filter(!grepl("CINAR_BC_73", metab_deriv_label)) %>%
  droplevels %>%
  tidyr::pivot_longer(., cols = !c("metab_deriv_label"),
                      names_to = "metabolite",
                      values_to = "concentration") %>%
  dplyr::filter(!(metabolite %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
  dplyr::left_join(., usvi_metabolomics_long.df %>%
                     dplyr::distinct(metabolites, LODflag, adaptedDervLabel, .keep_all = FALSE),
                   by = join_by("metabolite" == "metabolites", "metab_deriv_label" == "adaptedDervLabel")) %>%
  dplyr::filter(LODflag == "TRUE") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                         droplevels),
                   by = c("metab_deriv_label")) %>%
  dplyr::select(-LODflag) %>%
  dplyr::mutate(across(c(sample_id, metabolite, sampling_time, sampling_day, site, metab_deriv_label), ~factor(.x))) %>%
  droplevels


#not recommended to drop outlier values, but here's a way to determine them:
{
  # #these values of metabolite concentrations are potentially outliers, based on all measured concentrations across sites and times:
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

usvi_metab_filtered.df <- usvi_metab_summary.df

usvi_metab_filtered.df <- usvi_metab_filtered.df %>%
  dplyr::left_join(., (usvi_metab_filtered.df %>%
                         dplyr::group_by(metabolite) %>%
                         dplyr::summarise(med_conc = median(concentration, na.rm = TRUE),
                                          conc_5 = quantile(concentration, probs = 0.05, na.rm = TRUE, type = 7),
                                          conc_95 = quantile(concentration, probs = 0.95, na.rm = TRUE, type = 7)) %>%
                         dplyr::arrange(desc(conc_95)) %>%
                         dplyr::mutate(grouping = c(rep(1:5, each = 10, length.out = 49))) %>%
                         dplyr::select(metabolite, grouping)),
                   by = join_by(metabolite), relationship = "many-to-many", multiple = "all")

metab_colors <- c("#d699cf", "#43bb4e", "#a44ec9", "#74bf3e", "#6161d3", "#adb92f", "#bf3b9e", "#47c581", "#d63b83", "#4f8a29", "#d27ddf", "#9eb553", "#9384e2", "#dcae3d", "#4e62a7", "#e48c26", "#6a94db", "#e95a36", "#3abec8", "#de334d", "#59c9ab", "#be3220", "#4aa5d4", "#be5d27", "#3d9c7e", "#e86cbc", "#82bf71", "#854a98", "#96841e", "#a776b6", "#5f6913", "#cf4160", "#4e965b", "#b6507d", "#2f7038", "#e2828c", "#277257", "#e67c63", "#56642b", "#8e4f73", "#c0af5b", "#9a4555", "#7d8c42", "#a64837", "#a6b174", "#e39155", "#827233", "#d3a06f", "#af7821") %>%
  setNames(., c(usvi_metab_filtered.df %>%
                  dplyr::distinct(metabolite, grouping) %>%
                  dplyr::arrange(grouping) %>%
                  dplyr::select(metabolite) %>%
                  tibble::deframe(.)))
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
    
    temp_g <- print(
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
                                                                                                                                subtitle = "Removed 19 outlier values and CINAR_BC_73 prior to plotting", 
                                                                                                                                tag_levels = "A")
  
  g2_metab_filtered
  if(!any(grepl("metab_filtered", list.files(projectpath, pattern = "usvi_.*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_metab_filtered-", Sys.Date(), ".png"),
           g2_metab_filtered, 
           width = 16, height = 12, units = "in")
  }
  }

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


# NMDS on metabolomics profiles -------------------------------------------

potential_metab_outliers_idx <- c(73) %>%
  paste0("CINAR_BC_", .)

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
                                                .default = NA)) %>%
    droplevels
  


  # ylim <- bind_rows(nmds_metab_df, nmds_metab_log2_df) %>%
  #   dplyr::filter(!grepl("ASV", type)) %>%
  #   dplyr::select(NMDS2) %>%
  #   simplify
  # xlim <- bind_rows(nmds_metab_df, nmds_metab_log2_df) %>%
  #   dplyr::filter(!grepl("ASV", type)) %>%
  #   dplyr::select(NMDS1) %>%
  #   simplify
  # nmds_lims <- list(ylim = c(round(min(ylim), digits = 4), round(max(ylim), digits = 4)) * 1.1,
  #                   xlim = c(round(min(xlim), digits = 4), round(max(xlim), digits = 4)) * 1.1)
  
  
  # g3_all_ell <- (ggplot(data = nmds_metab_df,
  #               aes(x = NMDS1, y = NMDS2))
  #        + theme(text = element_text(size=14))
  #        + geom_point(data = (nmds_metab_df %>%
  #                               dplyr::filter(type == "sample") %>%
  #                               droplevels),
  #                     aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                     size = 5, stroke = 2, alpha = 1)
  #        + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
  #                                    direction = "both", segment.color = NA, seed = 48105, force = 1, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "black", fontface = "bold")
  #        + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                       type = "norm", 
  #                       level = 0.95,
  #                       geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #        + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #        + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #        + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #        + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                panel.grid = element_blank(),
  #                legend.position = "bottom",
  #                legend.key = element_blank(),
  #                legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                legend.text = element_text(size = 12, colour = "grey30"))
  #        + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                     labeller = labeller(site = site_lookup))
  #        + guides(color = "none",
  #                 fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                     override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                 shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                      override.aes = list(color = "black", stroke = 1, size = 2)))
  #        + coord_cartesian(ylim = nmds_lims[["ylim"]],
  #                          xlim = nmds_lims[["xlim"]],
  #                          expand = TRUE)
  # )
  # 
  # g3_all_log2_ell <- (ggplot(data = nmds_metab_log2_df,
  #                   aes(x = NMDS1, y = NMDS2))
  #            + theme(text = element_text(size=14))
  #            + geom_point(data = (nmds_metab_log2_df %>%
  #                                   dplyr::filter(type == "sample") %>%
  #                                   droplevels),
  #                         aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                         size = 5, stroke = 2, alpha = 1)
  #            + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
  #                                        direction = "both", segment.color = NA, seed = 48105,
  #                                        force = 1, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "black", fontface = "bold")
  #            + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                           type = "norm", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #            + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #            + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #            + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #            + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                    panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                    panel.grid = element_blank(),
  #                    legend.position = "bottom",
  #                    legend.key = element_blank(),
  #                    legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                    legend.text = element_text(size = 12, colour = "grey30"))
  #            + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                         labeller = labeller(site = site_lookup))
  #            + guides(color = "none",
  #                     fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                         override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                     shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                          override.aes = list(color = "black", stroke = 1, size = 2)))
  #            + coord_cartesian(ylim = nmds_lims[["ylim"]],
  #                              xlim = nmds_lims[["xlim"]],
  #                              expand = TRUE)
  # )
  }

g3_all_log2_ell <- P_nmds(dataset = nmds_metab_log2_df, facet_form = "sampling_time ~ site")
g3_all_ell <- P_nmds(dataset = nmds_metab_df, facet_form = "sampling_time ~ site")


  # if(!any(grepl("nmds", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
  #   ggsave(paste0(projectpath, "/", "usvi_metabolomics_nmds-", Sys.Date(), ".png"),
  #          g3,
  #          width = 10, height = 8, units = "in")
  # }
  

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
  tibble::column_to_rownames(., var = "metab_deriv_label") %>%
  droplevels

dist_usvi_metab.d <- vegan::vegdist(usvi_metab.tbl, method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_metab.df <- vegan::betadisper(dist_usvi_metab.d, type = "median",
                                       usvi_metab.meta$grouping) %>%
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
dist_usvi_metab_log2.d <- vegan::vegdist(usvi_metab_log2.tbl, method = "bray", binary = FALSE, na.rm = TRUE)
dist_usvi_metab_log2.df <- vegan::betadisper(dist_usvi_metab_log2.d, 
                                             type = "median", #this is default
                                             usvi_metab.meta$grouping) %>%
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

g3_all_disp <- P_betadisp(dataset = dist_usvi_metab.df, 
           metric = "Bray-Curtis") 
# g3_disp + ggrepel::geom_text_repel(data = dist_usvi_metab.df, aes(label = not_outlier, x = site, 
#                                                                   y = dissimilarity, fill = sampling_time,
#                                                                   group = interaction(site, sampling_time)), size =4, seed = 48105, 
#                                                               stat = "identity",  position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
#                                                               colour = "grey", fontface = "bold") + 
#   ggrepel::geom_text_repel(data = dist_usvi_metab.df, aes(label = potential_outlier, x = site,
#                                y = dissimilarity, fill = sampling_time,
#                                group = interaction(site, sampling_time)), size =4,
#                            stat = "identity",  position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
#                            seed = 48105, 
#                            colour = "black", fontface = "bold")

g3_all_log2_disp <- P_betadisp(dataset = dist_usvi_metab_log2.df, 
                           metric = "Bray-Curtis")  

print(g3_all_disp)
print(g3_all_log2_disp)

{
  # g3_all_disp <- print(
  #   ggplot(data = dist_usvi_metab.df)
  #   + theme_bw()
  #   + geom_boxplot(aes(x = site, y = dissimilarity, 
  #                      group = interaction(site, sampling_time)), 
  #                  color = "black",
  #                  position = position_dodge2(padding = 0.2, preserve = "single"),
  #                  show.legend = FALSE, outliers = FALSE)
  #   + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
  #                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
  #                alpha = 1.0, size = 3)
  #   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #   + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "Bray-Curtis dissimilarity")
  #   + scale_x_discrete(labels = site_lookup, name = "Site")
  #   + ggrepel::geom_text_repel(aes(label = not_outlier, x = site, 
  #                   y = dissimilarity, fill = sampling_time,
  #                   group = interaction(site, sampling_time)), size =4, seed = 48105, 
  #               stat = "identity",  position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
  #               colour = "grey", fontface = "bold")
  #   + ggrepel::geom_text_repel(aes(label = potential_outlier, x = site,
  #                                  y = dissimilarity, fill = sampling_time,
  #                                  group = interaction(site, sampling_time)), size =4,
  #                              stat = "identity",  position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
  #                              seed = 48105, 
  #                              colour = "black", fontface = "bold")
  #   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #           panel.grid = element_blank(),
  #           legend.position = "right",
  #           legend.key = element_blank(),
  #           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           legend.text = element_text(size = 12, colour = "grey30"))
  #   + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                labeller = labeller(site = site_lookup))  
  #   + guides(color = "none",
  #            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                 override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
  # 
  # g3_all_log2_disp <- print(
  #   ggplot(data = dist_usvi_metab_log2.df)
  #   + theme_bw()
  #   + geom_boxplot(aes(x = site, y = dissimilarity, 
  #                      group = interaction(site, sampling_time)), 
  #                  color = "black",
  #                  position = position_dodge2(padding = 0.2, preserve = "single"),
  #                  show.legend = FALSE, outliers = FALSE)
  #   + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
  #                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
  #                alpha = 1.0, size = 3)
  #   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #   + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Bray-Curtis dissimilarity")
  #   + scale_x_discrete(labels = site_lookup, name = "Site")
  #   + ggrepel::geom_text_repel(aes(label = not_outlier, x = site, 
  #                                  y = dissimilarity, fill = sampling_time,
  #                                  group = interaction(site, sampling_time)), size =4, seed = 48105,
  #                              position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
  #                              colour = "grey", fontface = "bold")
  #   + ggrepel::geom_text_repel(aes(label = potential_outlier, x = site,
  #                                  y = dissimilarity, fill = sampling_time,
  #                                  group = interaction(site, sampling_time)), size =4, seed = 48105,
  #                              position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2), hjust = "outward", angle = 0, 
  #                              colour = "black", fontface = "bold")
  #   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #           panel.grid = element_blank(),
  #           legend.position = "right",
  #           legend.key = element_blank(),
  #           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           legend.text = element_text(size = 12, colour = "grey30"))
  #   + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                labeller = labeller(site = site_lookup))  
  #   + guides(color = "none",
  #            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                 override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
}


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

nmds_metab_sites_log2.df <- bind_rows(temp_df1, temp_df2, temp_df3) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                         droplevels),
                   by = c("label" = "metab_deriv_label")) %>%
  #add labels:
  dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
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

nmds_metab_sites.df <- bind_rows(temp_df1, temp_df2, temp_df3) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                         droplevels),
                   by = c("label" = "metab_deriv_label")) %>%
  #add labels:
  dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
                                              .default = NA)) %>%
  droplevels
}


g3_sites_log2_ell <- P_nmds(dataset = nmds_metab_sites_log2.df, facet_form = "sampling_time ~ site")
g3_sites_ell <- P_nmds(dataset = nmds_metab_sites.df, facet_form = "sampling_time ~ site")

# g3_sites_log2_ell <- (ggplot(data = nmds_metab_sites_log2.df,
#                        aes(x = NMDS1, y = NMDS2))
#                 + theme(text = element_text(size=14))
#                 + geom_point(data = (nmds_metab_sites_log2.df %>%
#                                        dplyr::filter(type == "sample") %>%
#                                        droplevels),
#                              aes(shape = site, fill = addNA(sampling_time)), color = "black",
#                              size = 5, stroke = 2, alpha = 1)
#                 + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
#                                            direction = "both", segment.color = NA, seed = 48105, force = 1, force_pull = 1,
#                                            point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "black", fontface = "bold")
#                 + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), type = "t", 
#                                level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
#                 + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#                 + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
#                 + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#                 + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
#                         panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#                         panel.grid = element_blank(),
#                         legend.position = "bottom",
#                         legend.key = element_blank(),
#                         legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#                         legend.text = element_text(size = 12, colour = "grey30"))
#                 + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
#                              labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
#                 + guides(color = "none",
#                          fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
#                                              override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
#                          shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
#                                               override.aes = list(color = "black", stroke = 1, size = 2)))
# )
# g3_sites_ell <- (ggplot(data = nmds_metab_sites.df,
#                           aes(x = NMDS1, y = NMDS2))
#                    + theme(text = element_text(size=14))
#                    + geom_point(data = (nmds_metab_sites.df %>%
#                                           dplyr::filter(type == "sample") %>%
#                                           droplevels),
#                                 aes(shape = site, fill = addNA(sampling_time)), color = "black",
#                                 size = 5, stroke = 2, alpha = 1)
#                    + ggrepel::geom_text_repel(aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
#                                               direction = "both", segment.color = NA, seed = 48105, force = 1, force_pull = 1,
#                                               point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "black", fontface = "bold")
#                    + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
#                                   type = "t", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
#                    + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#                    + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
#                    + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#                    + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
#                            panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#                            panel.grid = element_blank(),
#                            legend.position = "bottom",
#                            legend.key = element_blank(),
#                            legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#                            legend.text = element_text(size = 12, colour = "grey30"))
#                    + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
#                                 labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
#                    + guides(color = "none",
#                             fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
#                                                 override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
#                             shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
#                                                  override.aes = list(color = "black", stroke = 1, size = 2)))
# )
# print(g3_sites_ell)

# g3_log2_nmds <- (g3_all_log2_ell + ggtitle("All samples used in distance matrix") + guides(fill = "none", shape = "none")) / (g3_sites_log2_ell + ggtitle("Site and sampling time specific matrices")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Log2(x + 1) transformed metabolite concentrations")
# g3_nmds <- (g3_all_ell + ggtitle("All samples used in distance matrix") + guides(fill = "none", shape = "none")) / (g3_sites_ell + ggtitle("Site and sampling time specific matrices")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed metabolite concentrations")
# 
# print(g3_log2_nmds)
# print(g3_nmds)
# 
# #the relationships between samples and the outliers, are the same within each site and sampling time, regarldess of whether we use an all-sample distance matrix or distance matrices site-specific
# #use log2(x+1) transformed relative abundances, and all samples in a distance matrix, as input to NMDS


g3_all_nmds <- (g3_all_disp + ggtitle("Community dissimilarity distance from site average") + guides(fill = "none", shape = "none")) / (g3_all_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed metabolite concentrations")
g3_all_log2_nmds <- (g3_all_log2_disp + ggtitle("Community dissimilarity distance from site average" ) + guides(fill = "none", shape = "none")) / (g3_all_log2_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed metabolite concentrations")

print(g3_all_nmds)
print(g3_all_log2_nmds)

if(!any(grepl("bc_nmds", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_metabolomics_bc_nmds-", Sys.Date(), ".png"),
         g3_all_nmds,
         width = 14, height = 12, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_metabolomics_bc_nmds_log2-", Sys.Date(), ".png"),
         g3_all_log2_nmds,
         width = 14, height = 12, units = "in")
}


#site_specific:

#plot dispersion from centroid
#do it site-specific...?
usvi_metabolomics.tbl <- usvi_metabolomics.df %>%
  tibble::column_to_rownames(var = "metab_deriv_label")

meta.seawater.list <- metabolomics_sample_metadata %>%
  dplyr::filter(metab_deriv_label %in% rownames(usvi_metabolomics.tbl)) %>%
  dplyr::select(metab_deriv_label, sample_id, sampling_time, sampling_day, site) %>%
  # dplyr::select(!c(contains("label"), contains("dna_"))) %>%
  dplyr::mutate(grouping = interaction(site, sampling_time)) %>%
  tibble::column_to_rownames(., var = "metab_deriv_label") %>%
  droplevels 

meta.seawater.list <- meta.seawater.list[rownames(usvi_metabolomics.tbl),] %>%
  split(., f = .$site) 

dist_usvi_sites_metab.d <- metabolomics_sample_metadata %>%
  split(., f = .$site) %>%
  map(., ~.x %>% 
        droplevels %>%
        dplyr::select(metab_deriv_label) %>%
        dplyr::distinct(.) %>%
        dplyr::inner_join(., usvi_metabolomics.tbl %>%
                            tibble::rownames_to_column(var = "metab_deriv_label"),
                          by = join_by(metab_deriv_label)) %>%
        tibble::column_to_rownames(var = "metab_deriv_label") %>%
        as.matrix(.) %>% vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE))

dist_usvi_sites_metab.df <- map2(dist_usvi_sites_metab.d, meta.seawater.list,
                               ~vegan::betadisper(., type = "median",
                                                  add = TRUE,
                                                  .y$grouping) %>%
                                 purrr::pluck("distances") %>%
                                 tibble::enframe(value = "dissimilarity", name = "metab_deriv_label") %>%
                                 dplyr::left_join(., (metabolomics_sample_metadata %>%
                                                        dplyr::filter(grepl("seawater", sample_type)) %>%
                                                        dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                                                        dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
                                                        droplevels),
                                                  by = c("metab_deriv_label" = "metab_deriv_label")) %>%
                                 droplevels) %>%
  bind_rows(.)


dist_usvi_sites_metab_log2.d <- metabolomics_sample_metadata %>%
  split(., f = .$site) %>%
  map(., ~.x %>% 
        droplevels %>%
        dplyr::select(metab_deriv_label) %>%
        dplyr::distinct(.) %>%
        dplyr::inner_join(., (usvi_metabolomics.tbl %>%
                            apply(., 2, function(x) log2(x + 1)) %>%
                              tibble::as_tibble(rownames = "metab_deriv_label")),
                          by = join_by(metab_deriv_label)) %>%
        tibble::column_to_rownames(var = "metab_deriv_label") %>%
        as.matrix(.) %>% vegan::vegdist(., method = "bray", binary = FALSE, na.rm = TRUE))

dist_usvi_sites_metab_log2.df <- map2(dist_usvi_sites_metab_log2.d, meta.seawater.list,
                                    ~vegan::betadisper(., type = "median",
                                                       add = TRUE,
                                                       .y$grouping) %>%
                                      purrr::pluck("distances") %>%
                                      tibble::enframe(value = "dissimilarity", name = "metab_deriv_label") %>%
                                      dplyr::left_join(., (metabolomics_sample_metadata %>%
                                                             dplyr::filter(grepl("seawater", sample_type)) %>%
                                                             dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                                                             dplyr::distinct(metab_deriv_label, .keep_all = TRUE) %>%
                                                             droplevels),
                                                       by = c("metab_deriv_label" = "metab_deriv_label")) %>%
                                      droplevels) %>%
  bind_rows(.)


g3_sites_disp <- P_betadisp(dataset = dist_usvi_sites_metab.df, metric = "Bray-Curtis")
g3_sites_log2_disp <- P_betadisp(dataset = dist_usvi_sites_metab_log2.df, metric = "Bray-Curtis")


# Which metabolites drive the NMDS structure? -----------------------------

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
                                      direction = "both", segment.color = NA, seed = 48105, force = 1, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "grey", fontface = "bold")
           + geom_point(data = (temp_df1 %>%
                                  dplyr::filter(type == "sample") %>%
                                  droplevels),
                        aes(shape = site, fill = addNA(sampling_time)), color = "black",
                        size = 5, stroke = 2, alpha = 1)
           + ggrepel::geom_text_repel(data = (temp_df1 %>%
                                                dplyr::filter(type == "sample") %>%
                                                droplevels),
                                      aes(label = gsub("CINAR_BC_", "", label), x = NMDS1, y = NMDS2),
                                      direction = "both", segment.color = NA, seed = 48105, force = 1, force_pull = 1,
                                      point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "black", fontface = "bold")
           + stat_ellipse(data = (temp_df1 %>%
                                    dplyr::filter(type == "sample") %>%
                                    droplevels),
                          aes(color = sampling_time,  group = interaction(site, sampling_time)), 
                          type = "t", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
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
)
print(temp_g)

if(!any(grepl("nmds_arrows", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_metabolomics_nmds_arrows-", Sys.Date(), ".png"),
         temp_g,
         width = 16, height = 10, units = "in")
}


# NMDS on asvs ------------------------------------------------------------

#first calculate NMDS on ASVs using all sites' samples
#then do site-specific NMDS's

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
  dplyr::filter(!grepl("Metab_280", sample)) %>% #drop the sample corresponding to CINAR_BC_73: Metab_280
  dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
                                  log2(abundance+1), #log transform abundance (with +1 pseudocount)
                                  0)) 

usvi_asv_log2.tbl <- usvi_asv.tbl %>%
  tidyr::pivot_wider(., id_cols = "sample",
                     # values_from = "abundance",
                     values_from = "logabund",
                     names_from = "asv_id") %>%
  tibble::column_to_rownames(var = "sample")

usvi_asv.tbl <- usvi_asv.tbl %>%
  tidyr::pivot_wider(., id_cols = "sample",
                     values_from = "abundance",
                     # values_from = "logabund",
                     names_from = "asv_id") %>%
  tibble::column_to_rownames(var = "sample")

# meta.seawater <- ps_usvi %>%
#   phyloseq::sample_data(.) %>%
#   tibble::as_tibble(rownames = "sample_id") %>%
#   dplyr::filter(sample_id %in% rownames(usvi_asv.tbl)) %>%
#   dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
#   dplyr::select(!c(contains("label"), contains("dna_"))) %>%
#   tibble::column_to_rownames(., var = "sample_id") %>%
#   droplevels

#exploratory NMDS
#see if in the DNA based NMDS whether Metab_280 and Metab_236 look like their replicates
#LB_seagrass peak photo day 3: Metab_280, Metab_286, Metab_289
#Tektite peak photo day day 2: Metab_236, Metab_227, Metab_224

nmds_asv.df <- vegan::metaMDS(usvi_asv.tbl, 
                              # distance = "bray",
                              trymax = 100,
                              distance = "horn",
                              autotransform = TRUE)

# env.sw <- vegan::envfit(nmds.asv, meta.seawater, permutations = 999, 
#                         na.rm = TRUE)

nmds_asv.df <- (ggplot2::fortify(nmds_asv.df) %>%
                  dplyr::filter(score == "sites") %>%
                  dplyr::mutate(type = "sample") %>%
                  dplyr::select(-score)) %>%
  # bind_rows(., (ggplot2::fortify(env.sw))) %>%
  bind_rows(., (ggplot2::fortify(nmds_asv.df) %>%
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


nmds_asv_log2.df <- vegan::metaMDS(usvi_asv_log2.tbl, 
                                   # distance = "bray",
                                   trymax = 100,
                                   distance = "horn",
                                   autotransform = TRUE)

nmds_asv_log2.df <- (ggplot2::fortify(nmds_asv_log2.df) %>%
                       dplyr::filter(score == "sites") %>%
                       dplyr::mutate(type = "sample") %>%
                       dplyr::select(-score)) %>%
  # bind_rows(., (ggplot2::fortify(env.sw))) %>%
  bind_rows(., (ggplot2::fortify(nmds_asv_log2.df) %>%
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


# 
# 
# ylim <- bind_rows(nmds_asv.df, nmds_asv_log2.df) %>%
#   dplyr::filter(!grepl("ASV", type)) %>%
#   dplyr::select(NMDS2) %>%
#   simplify
# xlim <- bind_rows(nmds_asv.df, nmds_asv_log2.df) %>%
#   dplyr::filter(!grepl("ASV", type)) %>%
#   dplyr::select(NMDS1) %>%
#   simplify
# nmds_lims <- list(ylim = c(round(min(ylim), digits = 4), round(max(ylim), digits = 4)) * 1.1,
#                   xlim = c(round(min(xlim), digits = 4), round(max(xlim), digits = 4)) * 1.1)

g4_all_ell <- P_nmds(dataset = nmds_asv.df %>%
         dplyr::filter(type == "sample") %>%
         droplevels,
       facet_form = "sampling_time ~ site")

g4_all_log2_ell <- P_nmds(dataset = nmds_asv_log2.df %>%
                            dplyr::filter(type == "sample") %>%
                            droplevels,
                          facet_form = "sampling_time ~ site")


print(g4_all_log2_ell)
{
  # g4_all_ell <- (ggplot(data = nmds_asv.df %>%
  #                         dplyr::filter(type == "sample") %>%
  #                         droplevels,
  #                   aes(x = NMDS1, y = NMDS2))
  #            + theme(text = element_text(size=14))
  #            + geom_point(data = (nmds_asv.df %>%
  #                                   dplyr::filter(type == "sample") %>%
  #                                   droplevels),
  #                         aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                         size = 5, stroke = 2, alpha = 1)
  #            + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                           type = "norm", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #            + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #            + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #            + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #            + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "NMDS2")
  #            + scale_x_continuous(expand = expansion(mult = c(0.1,0.1)), name = "NMDS1")
  #            + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                    panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                    panel.grid = element_blank(),
  #                    legend.position = "bottom",
  #                    legend.key = element_blank(),
  #                    legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                    legend.text = element_text(size = 12, colour = "grey30"))
  #            + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed",
  #                         labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  #            + guides(color = "none",
  #                     fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                         override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                     shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                          override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
  # g4_all_log2_ell <- (ggplot(data = nmds_asv_log2.df %>%
  #                              dplyr::filter(type == "sample") %>%
  #                              droplevels,
  #                        aes(x = NMDS1, y = NMDS2))
  #                 + theme(text = element_text(size=14))
  #                 + geom_point(data = (nmds_asv_log2.df %>%
  #                                        dplyr::filter(type == "sample") %>%
  #                                        droplevels),
  #                              aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                              size = 5, stroke = 2, alpha = 1)
  #                 + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                                type = "norm", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #                 + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                 + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #                 + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                 + scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS2")
  #                 + scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), name = "NMDS1")
  #                 + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                         panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                         panel.grid = element_blank(),
  #                         legend.position = "bottom",
  #                         legend.key = element_blank(),
  #                         legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                         legend.text = element_text(size = 12, colour = "grey30"))
  #                 + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                              labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  #                 + guides(color = "none",
  #                          fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                              override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                          shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                               override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
}




#calculate dispersions of sites' dissimilarity metrics from the site-specific mean

dist_usvi_asv.d <- vegan::vegdist(usvi_asv.tbl, method = "horn", binary = FALSE, na.rm = TRUE)
dist_usvi_asv.df <- vegan::betadisper(dist_usvi_asv.d, type = "median",
                                      add = TRUE,
                                      meta.seawater$grouping) %>%
  purrr::pluck("distances") %>%
  tibble::enframe(value = "dissimilarity", name = "sample_id") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                         dplyr::distinct(sample_id, .keep_all = TRUE) %>%
                         droplevels),
                   by = c("sample_id" = "sample_id")) %>%
  dplyr::mutate(potential_outlier = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", metab_deriv_label) | metab_deriv_label %in% potential_metab_outliers_idx)~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),
                                                     .default = NA)) %>%
  dplyr::mutate(not_outlier = dplyr::case_when(is.na(potential_outlier) ~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),
                                               .default = NA)) %>%
  droplevels
dist_usvi_asv_log2.d <- vegan::vegdist(usvi_asv_log2.tbl, method = "horn", binary = FALSE, na.rm = TRUE)
dist_usvi_asv_log2.df <- vegan::betadisper(dist_usvi_asv_log2.d, 
                                           type = "median", #this is default
                                           add = TRUE,
                                           meta.seawater$grouping) %>%
  purrr::pluck("distances") %>%
  tibble::enframe(value = "dissimilarity", name = "sample_id") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                         dplyr::distinct(sample_id, .keep_all = TRUE) %>%
                         droplevels),
                   by = c("sample_id" = "sample_id")) %>%
  dplyr::mutate(potential_outlier = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", metab_deriv_label) | metab_deriv_label %in% potential_metab_outliers_idx)~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),
                                                     .default = NA)) %>%
  dplyr::mutate(not_outlier = dplyr::case_when(is.na(potential_outlier) ~ stringr::str_remove_all(metab_deriv_label, "CINAR_BC_"),
                                               .default = NA)) %>%
  droplevels

g4_all_disp <- P_betadisp(dataset = dist_usvi_asv.df, metric = "Morisita-Horn")
g4_all_log2_disp <- P_betadisp(dataset = dist_usvi_asv_log2.df, metric = "Morisita-Horn")
  
# g4_all_disp <- print(
#   ggplot(data = dist_usvi_asv.df)
#   + theme_bw()
#   + geom_boxplot(aes(x = site, y = dissimilarity, 
#                      group = interaction(site, sampling_time)), 
#                  color = "black",
#                  position = position_dodge2(padding = 0.2, preserve = "single"),
#                  show.legend = FALSE, outliers = FALSE)
#   + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
#                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
#                alpha = 1.0, size = 3)
#   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
#   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "Morisita-Horn dissimilarity")
#   + scale_x_discrete(labels = site_lookup, name = "Site")
#   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#           panel.grid = element_blank(),
#           legend.position = "right",
#           legend.key = element_blank(),
#           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           legend.text = element_text(size = 12, colour = "grey30"))
#   + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
#                labeller = labeller(site = site_lookup))  
#   + guides(color = "none",
#            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
#                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
#            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
#                                 override.aes = list(color = "black", stroke = 1, size = 2)))
# )
# 
# g4_all_log2_disp <- print(
#   ggplot(data = dist_usvi_asv_log2.df)
#   + theme_bw()
#   + geom_boxplot(aes(y = dissimilarity, 
#                      x = site,
#                      # x = interaction(site, sampling_time),
#                      group = interaction(site, sampling_time)), 
#                  color = "black",
#                  position = position_dodge2(padding = 0.2, preserve = "single"),
#                  show.legend = FALSE, outliers = FALSE)
#   + geom_point(aes(y = dissimilarity, 
#                    x = site,
#                    # x = interaction(site, sampling_time),
#                    fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
#                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
#                alpha = 1.0, size = 3)
#   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
#   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "Morisita-Horn dissimilarity")
#   + scale_x_discrete(labels = site_lookup, name = "Site")
#   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#           panel.grid = element_blank(),
#           legend.position = "right",
#           legend.key = element_blank(),
#           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           legend.text = element_text(size = 12, colour = "grey30"))
#   + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed",
#                labeller = labeller(site = site_lookup))
#   + guides(color = "none",
#            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
#                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
#            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
#                                 override.aes = list(color = "black", stroke = 1, size = 2)))
# )



# Site-specific ASV NMDS --------------------------------------------------

#compare the dispersions using the distances calculated specific to each site, to the distances calculated when all sites are considered



#note that for the asv table, either relative abundance or log2(+1) transformed
#that bray-curtis distance doesn't allow for convergence of a solution
#but that using Horn distance works within 50 tries

#site specific log-transformed:
{
  temp_df1 <- usvi_asv_log2.tbl %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(sample_id %in% (metabolomics_sample_metadata %>%
                                    dplyr::filter(grepl("LB", site)) %>%
                                    droplevels %>%
                                    dplyr::select(sample_id) %>%
                                    tibble::deframe(.))) %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., trymax = 50,
                   distance = "horn", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score) 
  temp_df2 <-  usvi_asv_log2.tbl %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(sample_id %in% (metabolomics_sample_metadata %>%
                                    dplyr::filter(grepl("Yawzi", site)) %>%
                                    droplevels %>%
                                    dplyr::select(sample_id) %>%
                                    tibble::deframe(.))) %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., trymax = 50,
                   distance = "horn", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score)
  
  temp_df3 <-  usvi_asv_log2.tbl %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(sample_id %in% (metabolomics_sample_metadata %>%
                                    dplyr::filter(grepl("Tektite", site)) %>%
                                    droplevels %>%
                                    dplyr::select(sample_id) %>%
                                    tibble::deframe(.))) %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., trymax = 50,
                   distance = "horn", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score)
  
  nmds_asv_sites_log2.df <- bind_rows(temp_df1, temp_df2, temp_df3) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                           dplyr::distinct(sample_id, .keep_all = TRUE) %>%
                           droplevels),
                     by = c("label" = "sample_id")) %>%
    #add labels:
    # dplyr::mutate(text_label = dplyr::case_when((grepl("_73|_43|_44|_42|_70|_74|_69|_67|_71|_101|_108|_105", label) | label %in% potential_metab_outliers_idx)~ label,
    #                                             .default = NA)) %>%
    droplevels
}

#site specific untransformed:
{
  temp_df1 <- usvi_asv.tbl %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(sample_id %in% (metabolomics_sample_metadata %>%
                                    dplyr::filter(grepl("LB", site)) %>%
                                    droplevels %>%
                                    dplyr::select(sample_id) %>%
                                    tibble::deframe(.))) %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., trymax = 50,
                   distance = "horn", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score) 
  temp_df2 <-  usvi_asv.tbl %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(sample_id %in% (metabolomics_sample_metadata %>%
                                    dplyr::filter(grepl("Yawzi", site)) %>%
                                    droplevels %>%
                                    dplyr::select(sample_id) %>%
                                    tibble::deframe(.))) %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., trymax = 50,
                   distance = "horn", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score)
  
  temp_df3 <-  usvi_asv.tbl %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(sample_id %in% (metabolomics_sample_metadata %>%
                                    dplyr::filter(grepl("Tektite", site)) %>%
                                    droplevels %>%
                                    dplyr::select(sample_id) %>%
                                    tibble::deframe(.))) %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., trymax = 50,
                   distance = "horn", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    as.data.frame %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score)
  
  nmds_asv_sites.df <- bind_rows(temp_df1, temp_df2, temp_df3) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                           dplyr::distinct(sample_id, .keep_all = TRUE) %>%
                           droplevels),
                     by = c("label" = "sample_id")) %>%
    droplevels
}

#plot dispersion from centroid
#do it site-specific...?

meta.seawater <- ps_usvi %>%
  phyloseq::sample_data(.) %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  dplyr::filter(sample_id %in% rownames(usvi_asv.tbl)) %>%
  dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
  dplyr::select(!c(contains("label"), contains("dna_"))) %>%
  dplyr::mutate(grouping = interaction(site, sampling_time)) %>%
  tibble::column_to_rownames(., var = "sample_id") %>%
  droplevels

meta.seawater.list <- meta.seawater[rownames(usvi_asv.tbl),] %>%
  split(., f = .$site) 

dist_usvi_sites_asv.d <- metabolomics_sample_metadata %>%
  split(., f = .$site) %>%
  map(., ~.x %>% 
        droplevels %>%
        dplyr::select(sample_id) %>%
        dplyr::distinct(.) %>%
        dplyr::inner_join(., usvi_asv.tbl %>%
                           tibble::rownames_to_column(var = "sample_id"),
                         by = join_by(sample_id)) %>%
          tibble::column_to_rownames(var = "sample_id") %>%
        as.matrix(.) %>% vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE))
  
dist_usvi_sites_asv.df <- map2(dist_usvi_sites_asv.d, meta.seawater.list,
                               ~vegan::betadisper(., type = "median",
                                                  add = TRUE,
                                                  .y$grouping) %>%
        purrr::pluck("distances") %>%
        tibble::enframe(value = "dissimilarity", name = "sample_id") %>%
        dplyr::left_join(., (metabolomics_sample_metadata %>%
                               dplyr::filter(grepl("seawater", sample_type)) %>%
                               dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                               dplyr::distinct(sample_id, .keep_all = TRUE) %>%
                               droplevels),
                         by = c("sample_id" = "sample_id")) %>%
        droplevels) %>%
  bind_rows(.)


dist_usvi_sites_asv_log2.d <- metabolomics_sample_metadata %>%
  split(., f = .$site) %>%
  map(., ~.x %>% 
        droplevels %>%
        dplyr::select(sample_id) %>%
        dplyr::distinct(.) %>%
        dplyr::inner_join(., usvi_asv_log2.tbl[rownames(usvi_asv.tbl),] %>%
                            tibble::rownames_to_column(var = "sample_id"),
                          by = join_by(sample_id)) %>%
        tibble::column_to_rownames(var = "sample_id") %>%
        as.matrix(.) %>% vegan::vegdist(., method = "horn", binary = FALSE, na.rm = TRUE))

dist_usvi_sites_asv_log2.df <- map2(dist_usvi_sites_asv_log2.d, meta.seawater.list,
                               ~vegan::betadisper(., type = "median",
                                                  add = TRUE,
                                                  .y$grouping) %>%
                                 purrr::pluck("distances") %>%
                                 tibble::enframe(value = "dissimilarity", name = "sample_id") %>%
                                 dplyr::left_join(., (metabolomics_sample_metadata %>%
                                                        dplyr::filter(grepl("seawater", sample_type)) %>%
                                                        dplyr::select(sample_id, metab_deriv_label, sample_type, sampling_date, sampling_time, sampling_day, site) %>%
                                                        dplyr::distinct(sample_id, .keep_all = TRUE) %>%
                                                        droplevels),
                                                  by = c("sample_id" = "sample_id")) %>%
                                 droplevels) %>%
  bind_rows(.)


g4_sites_disp <- P_betadisp(dataset = dist_usvi_sites_asv.df, metric = "Morisita-Horn")
g4_sites_log2_disp <- P_betadisp(dataset = dist_usvi_sites_asv_log2.df, metric = "Morisita-Horn")

# g4_sites_disp <- print(
#   ggplot(data = dist_usvi_sites_asv.df)
#   + theme_bw()
#   + geom_boxplot(aes(x = site, y = dissimilarity, 
#                      group = interaction(site, sampling_time)), 
#                  color = "black",
#                  position = position_dodge2(padding = 0.2, preserve = "single"),
#                  show.legend = FALSE, outliers = FALSE)
#   + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
#                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
#                alpha = 1.0, size = 3)
#   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
#   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "Morisita-Horn dissimilarity")
#   + scale_x_discrete(labels = site_lookup, name = "Site")
#   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#           panel.grid = element_blank(),
#           legend.position = "right",
#           legend.key = element_blank(),
#           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           legend.text = element_text(size = 12, colour = "grey30"))
#   + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
#                labeller = labeller(site = site_lookup))  
#   + guides(color = "none",
#            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
#                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
#            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
#                                 override.aes = list(color = "black", stroke = 1, size = 2)))
# )
# 
# g4_sites_log2_disp <- print(
#   ggplot(data = dist_usvi_sites_asv_log2.df)
#   + theme_bw()
#   + geom_boxplot(aes(x = site, y = dissimilarity, 
#                      group = interaction(site, sampling_time)), 
#                  color = "black",
#                  position = position_dodge2(padding = 0.2, preserve = "single"),
#                  show.legend = FALSE, outliers = FALSE)
#   + geom_point(aes(x = site, y = dissimilarity, fill = sampling_time, group = interaction(site, sampling_time), shape = site), 
#                position = position_jitterdodge(dodge.width = 0.75, seed = 48105, jitter.width = 0.2),
#                alpha = 1.0, size = 3)
#   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
#   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
#   + scale_y_continuous(expand = expansion(mult = c(0.1,0.1)), name = "Morisita-Horn dissimilarity")
#   + scale_x_discrete(labels = site_lookup, name = "Site")
#   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
#           panel.grid = element_blank(),
#           legend.position = "right",
#           legend.key = element_blank(),
#           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
#           legend.text = element_text(size = 12, colour = "grey30"))
#   + facet_grid(. ~ site, drop = TRUE, scales = "free", space = "fixed", 
#                labeller = labeller(site = site_lookup))  
#   + guides(color = "none",
#            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
#                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
#            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
#                                 override.aes = list(color = "black", stroke = 1, size = 2)))
# )




g4_disp <- (g4_all_disp + ggtitle("Community dissimilarity distance from site average, using all samples") + guides(fill = "none", shape = "none")) / (g4_sites_disp + ggtitle("Community dissimilarity distance from site average, site-specific distance matrix")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed relative abundances")
g4_log2_disp <- (g4_all_log2_disp + ggtitle("Community dissimilarity distance from site average, using all samples" ) + guides(fill = "none", shape = "none")) / (g4_sites_log2_disp + ggtitle("Community dissimilarity distance from site average, site-specific distance matrix")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed relative abundances")

print(g4_disp)
print(g4_log2_disp)

#while the range of MH distances is different when calculating the distance matrix using all samples, 
# compared to doing 3 site-specific matrices
#the relationship between samples is the same
#this is true for untrnsformed, and log2(x+1) transformed relative abundances



#what do the NMDS look like?

# ylim <- bind_rows(nmds_asv_sites.df, nmds_asv_sites_log2.df) %>%
#   dplyr::filter(!grepl("ASV", type)) %>%
#   dplyr::select(NMDS2) %>%
#   simplify
# xlim <- bind_rows(nmds_asv_sites.df, nmds_asv_sites_log2.df) %>%
#   dplyr::filter(!grepl("ASV", type)) %>%
#   dplyr::select(NMDS1) %>%
#   simplify
# nmds_lims <- list(ylim = c(round(min(ylim), digits = 4), round(max(ylim), digits = 4)) * 1.1,
#                   xlim = c(round(min(xlim), digits = 4), round(max(xlim), digits = 4)) * 1.1)


g4_sites_ell <- P_nmds(dataset = nmds_asv_sites.df, 
                       facet_form = "sampling_time ~ site")
g4_sites_log2_ell <- P_nmds(dataset = nmds_asv_sites_log2.df, 
                            facet_form = "sampling_time ~ site")

{
  # g4_sites_ell <- (ggplot(data = nmds_asv_sites.df,
  #                         aes(x = NMDS1, y = NMDS2))
  #                  + theme(text = element_text(size=14))
  #                  + geom_point(data = (nmds_asv_sites.df %>%
  #                                         dplyr::filter(type == "sample") %>%
  #                                         droplevels),
  #                               aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                               size = 5, stroke = 2, alpha = 1)
  #                  + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                                 type = "norm", 
  #                                 level = 0.95,
  #                                 geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #                  + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                  + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #                  + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                  + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                          panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                          panel.grid = element_blank(),
  #                          legend.position = "bottom",
  #                          legend.key = element_blank(),
  #                          legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                          legend.text = element_text(size = 12, colour = "grey30"))
  #                  + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed",
  #                               labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  #                  + guides(color = "none",
  #                           fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                               override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                           shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                                override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
  # 
  # g4_sites_log2_ell <- (ggplot(data = nmds_asv_sites_log2.df,
  #                              aes(x = NMDS1, y = NMDS2))
  #                       + theme(text = element_text(size=14))
  #                       + geom_point(data = (nmds_asv_sites_log2.df %>%
  #                                              dplyr::filter(type == "sample") %>%
  #                                              droplevels),
  #                                    aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                                    size = 5, stroke = 2, alpha = 1)
  #                       + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                                      type = "norm", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #                       + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                       + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #                       + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #                       + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                               panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #                               panel.grid = element_blank(),
  #                               legend.position = "bottom",
  #                               legend.key = element_blank(),
  #                               legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #                               legend.text = element_text(size = 12, colour = "grey30"))
  #                       + facet_grid(sampling_time ~ site, drop = TRUE, scales = "free", space = "fixed", 
  #                                    labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  #                       + guides(color = "none",
  #                                fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                                    override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #                                shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                                     override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
}


g4_all_nmds <- (g4_all_ell + ggtitle("NMDS using all samples") + guides(fill = "none", shape = "none")) / (g4_sites_ell + ggtitle("NMDS per site-specific distance matrix") + guides(fill = "none", shape = "none")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed relative abundances")
g4_all_log2_nmds <- (g4_all_log2_ell + ggtitle("NMDS using all samples" ) + guides(fill = "none", shape = "none")) / (g4_sites_log2_ell + ggtitle("NMDS per site-specific distance matrix") + guides(fill = "none", shape = "none")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed relative abundances")

print(g4_all_nmds)
print(g4_all_log2_nmds)


#confirming outliers are the same in both the site-specific distance matrices NMDS, and all-sample NMDS:
#calculate the squared product of NMDS1 and NMDS2 for each sample (cross_NMDS)
#then use outliers to identify outlier cross_NMDS for each site and sampling time 
#then filter for those samples where the ratio of cross_NMDS/outlier_cross_NMDS ~ 1 +/- 0.05
#compare which outliers were identified in the all-sample NMDS, and where they are inthe site-specific NMDS

#summary of coordinates for outliers identified:
#LB seagrass: Metab_219 and Metab_289 were identified as outliers in both metrics, but are both within the 95% ellipse in both the site-specific and all-sample NMDS
#Yawzi: Metab_228 (true outlier in both) and Metab_267 (still within the 95% ellipse in both)
#Tektite: in both metrics, Metab_241 and Metab_224 are true outliers and outside the 95% ellipse in both

{
  # temp_df1 <- nmds_asv_log2.df %>%
  #   dplyr::filter(type == "sample") %>%
  #   dplyr::select(label, NMDS1, NMDS2, site, sampling_time) %>%
  #   dplyr::distinct(label, .keep_all = TRUE) %>%
  #   dplyr::mutate(metric = "all") %>%
  #   dplyr::mutate(cross_NMDS = (NMDS1 * NMDS2)^2) %>%
  #   dplyr::group_by(site, sampling_time, metric) %>%
  #   dplyr::reframe(across(contains("cross_NMDS"),
  #                         list(outlier_upper = function(x) (outliers::outlier(x, opposite = FALSE))))) %>%
  #   dplyr::right_join(nmds_asv_log2.df %>%
  #                       dplyr::filter(type == "sample") %>%
  #                       dplyr::select(label, NMDS1, NMDS2, site, sampling_time) %>%
  #                       dplyr::distinct(label, .keep_all = TRUE) %>%
  #                       dplyr::mutate(metric = "all") %>%
  #                       dplyr::mutate(cross_NMDS = (NMDS1 * NMDS2)^2) %>%
  #                       droplevels, by = join_by(site, sampling_time, metric)) %>%
  #   dplyr::rowwise(.) %>%
  #   dplyr::mutate(ratio_NMDS_upper = cross_NMDS / cross_NMDS_outlier_upper) %>%
  #   dplyr::mutate(potential_outlier = dplyr::case_when((ratio_NMDS_upper == 1) | (ratio_NMDS_upper < 1.05 & ratio_NMDS_upper > 0.95) ~ label, .default = NA)) %>%
  #   droplevels
  # 
  # temp_df2 <- nmds_asv_sites_log2.df %>%
  #   dplyr::filter(type == "sample") %>%
  #   dplyr::select(label, NMDS1, NMDS2, site, sampling_time) %>%
  #   dplyr::distinct(label, .keep_all = TRUE) %>%
  #   dplyr::mutate(metric = "sites") %>%
  #   dplyr::mutate(cross_NMDS = (NMDS1 * NMDS2)^2) %>%
  #   dplyr::group_by(site, sampling_time, metric) %>%
  #   dplyr::reframe(across(contains("cross_NMDS"),
  #                         list(outlier_upper = function(x) (outliers::outlier(x, opposite = FALSE))))) %>%
  #   dplyr::right_join(nmds_asv_sites_log2.df %>%
  #                       dplyr::filter(type == "sample") %>%
  #                       dplyr::select(label, NMDS1, NMDS2, site, sampling_time) %>%
  #                       dplyr::distinct(label, .keep_all = TRUE) %>%
  #                       dplyr::mutate(metric = "sites") %>%
  #                       dplyr::mutate(cross_NMDS = (NMDS1 * NMDS2)^2) %>%
  #                       droplevels, by = join_by(site, sampling_time, metric)) %>%
  #   dplyr::rowwise(.) %>%
  #   dplyr::mutate(ratio_NMDS_upper = cross_NMDS / cross_NMDS_outlier_upper) %>%
  #   dplyr::mutate(potential_outlier = dplyr::case_when((ratio_NMDS_upper == 1) | (ratio_NMDS_upper < 1.05 & ratio_NMDS_upper > 0.95) ~ label, .default = NA)) %>%
  #   droplevels
  # 
  # temp_df3 <- temp_df2 %>%
  #   dplyr::select(-potential_outlier) %>%
  #   dplyr::left_join(., temp_df1 %>%
  #                      dplyr::ungroup(.) %>%
  #                      dplyr::select(site, sampling_time, label, potential_outlier) %>%
  #                      droplevels,
  #                    by = join_by(site, sampling_time, label), relationship = "many-to-many", multiple = "all") %>%
  #   bind_rows(., temp_df1) %>%
  #   droplevels
  # 
  # (ggplot(data = temp_df3 %>%
  #           dplyr::filter(grepl("LB", site)) %>%
  #           droplevels,
  #         aes(x = NMDS1, y = NMDS2))
  #   + theme(text = element_text(size=14))
  #   + geom_point(aes(shape = site, fill = addNA(sampling_time)), color = "black",
  #                size = 5, stroke = 2, alpha = 1)
  #   + ggrepel::geom_text_repel(aes(label = potential_outlier, x = NMDS1, y = NMDS2),
  #                              direction = "both", segment.color = NA, seed = 48105, force = 1, force_pull = 1, point.padding = unit(0.01, "npc"), box.padding = unit(0.01, "npc"), colour = "black", fontface = "bold")
  #   + stat_ellipse(aes(color = sampling_time,  group = interaction(site, sampling_time)), 
  #                  type = "norm", level = 0.95, geom = "polygon", alpha = 0.05, show.legend = FALSE)
  #   + scale_color_manual(name = "sampling time", values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #   + scale_shape_manual(values = c(22, 21, 23), labels = c(site_lookup, "NA"), breaks = c(names(site_lookup), NA))
  #   + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  #   + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
  #           panel.grid = element_blank(),
  #           legend.position = "bottom",
  #           legend.key = element_blank(),
  #           legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
  #           legend.text = element_text(size = 12, colour = "grey30"))
  #   + facet_grid(sampling_time ~ metric, drop = TRUE, scales = "free", space = "fixed",
  #                labeller = labeller(site = site_lookup, sampling_time = sampling_time_lookup))
  #   + guides(color = "none",
  #            fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
  #                                override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
  #            shape = guide_legend(order = 1, ncol = 1, title = "Site", direction = "vertical",
  #                                 override.aes = list(color = "black", stroke = 1, size = 2)))
  # )
}


#for clarity of methods in the manuscript, use the all-sample NMDS to make the figures



# g4 <- (g4_all_disp + ggtitle("Community dissimilarity distance from site average") + guides(fill = "none", shape = "none")) / (g4_all_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed relative abundances")
# g4_log2 <- (g4_all_log2_disp + ggtitle("Community dissimilarity distance from site average" ) + guides(fill = "none", shape = "none")) / (g4_all_log2_ell + ggtitle("NMDS by site and sampling time")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed relative abundances")

g_nmds <- (g3_all_ell + ggtitle("Metabolite concentrations") + guides(fill = "none", shape = "none")) / (g4_all_ell + ggtitle("Microbial relative abundances") + guides(fill = "none", shape = "none")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Untransformed", tag_levels = "A")
g_nmds_log2 <- (g3_all_log2_ell + ggtitle("Metabolite concentrations") + guides(fill = "none", shape = "none")) / (g4_all_log2_ell + ggtitle("Microbial relative abundances")+ guides(fill = "none", shape = "none")) + patchwork::plot_layout(guides = "collect") +  patchwork::plot_annotation(title = "Log2(x + 1) transformed", tag_levels = "A")


#more work than wanted:
# P_format_multipanel <- function(...){
#   is_valid_plot <- function(x) is.ggplot(x) || is.grob(x) #borrowed from patchwork
#   if (is_valid_plot(..1)) {
#     plots <- list(...)
#   } else if (is.list(..1)) {
#     plots <- ..1
#   } else {
#     cli_abort('Can only wrap {.cls ggplot} and/or {.cls grob} objects or a list of them')
#   }
#   lapply(plots, function(x) (x + theme(strip.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))) %>%
#     purrr::reduce(., `+`) + 
#     patchwork::plot_layout(guides = "collect", tag_level = "new") + patchwork::plot_annotation(tag_levels = "1") & guides(shape = "none")
# }
# 
# 
# g3_all_log2_ell_panels <- P_format_multipanel(P_nmds(dataset = nmds_metab_log2_df %>% dplyr::filter(grepl("LB", site)), facet_form = "sampling_time ~ site"), P_nmds(dataset = nmds_metab_log2_df %>% dplyr::filter(grepl("Tektite", site)), facet_form = "sampling_time ~ site"), P_nmds(dataset = nmds_metab_log2_df %>% dplyr::filter(grepl("Yawzi", site)), facet_form = "sampling_time ~ site")) +
#   patchwork::plot_layout(tag_level = "keep") +
#   patchwork::plot_annotation(tag_levels = list(c("A", "B", "C"))) 
# g4_all_log2_ell_panels <- P_format_multipanel(P_nmds(dataset = nmds_asv_log2.df %>% dplyr::filter(type == "sample") %>% droplevels %>% dplyr::filter(grepl("LB", site)), facet_form = "sampling_time ~ site"), P_nmds(dataset = nmds_asv_log2.df %>% dplyr::filter(type == "sample") %>% droplevels %>% dplyr::filter(grepl("Tektite", site)), facet_form = "sampling_time ~ site"), P_nmds(dataset = nmds_asv_log2.df %>% dplyr::filter(type == "sample") %>% droplevels %>% dplyr::filter(grepl("Yawzi", site)), facet_form = "sampling_time ~ site")) +
#   patchwork::plot_layout(tag_level = "keep") +
#   patchwork::plot_annotation(tag_levels = list(c("D", "E", "F")))
# 
# g_nmds_log2 <- (g3_all_log2_ell_panels + guides(fill = "none", shape = "none")) / (g4_all_log2_ell_panels + guides(fill = "none", shape = "none")) +  patchwork::plot_annotation(title = "Log2(x + 1) transformed")

g3_all_log2_disp / g3_sites_log2_disp

print(g_nmds)
print(g_nmds_log2)

g_betadisp <- (g3_all_disp + ggtitle("Metabolite concentrations") + guides(fill = "none", shape = "none")) / (g4_all_disp + ggtitle("Microbial relative abundances ") + guides(shape = "none")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Untransformed", tag_levels = "A")
g_betadisp_log2 <- (g3_all_log2_disp + ggtitle("Metabolite concentrations") + guides(fill = "none", shape = "none")) / (g4_all_log2_disp + ggtitle("Microbial relative abundances ") + guides(shape = "none")) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Log2(x + 1) transformed", tag_levels = "A")

print(g_betadisp)
print(g_betadisp_log2)



if(!any(grepl("nmds", list.files(projectpath, pattern = "usvi_metab_microb.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_metab_microb_nmds-", Sys.Date(), ".png"),
         g_nmds,
         width = 14, height = 12, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_metab_microb_log2_nmds-", Sys.Date(), ".png"),
         g_nmds_log2,
         width = 14, height = 12, units = "in")
}



if(!any(grepl("betadisp", list.files(projectpath, pattern = "usvi_metab_microb.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_metab_microb_betadisp-", Sys.Date(), ".png"),
         g_betadisp,
         width = 8, height = 12, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_metab_microb_log2_betadisp-", Sys.Date(), ".png"),
         g_betadisp_log2,
         width = 8, height = 12, units = "in")
}

