# 05_taxa_metabolites.R

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



{
 
  nmds_metab_df <- usvi_metabolomics.df %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
    as.matrix(.) %>%
    vegan::metaMDS(., distance = "bray", autotransform = TRUE) %>%
    ggplot2::fortify(.) %>%
    dplyr::filter(score == "sites") %>%
    dplyr::mutate(type = "sample") %>%
    dplyr::select(-score) %>%
    dplyr::left_join(., (metabolomics_sample_metadata %>%
                           dplyr::filter(grepl("seawater", sample_type)) %>%
                           dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, metab_deriv_label) %>%
                           droplevels),
                     by = c("label" = "metab_deriv_label")) %>%
    #add labels:
    dplyr::mutate(text_label = dplyr::case_when(grepl("73|43|44|42|70|74|69|67|71|101|108|105", label) ~ label,
                                                .default = NA)) %>%
    droplevels
  
  
  # nmds_metab_df %>%
  #   # dplyr::filter(grepl("LB", site))
  #   dplyr::filter(grepl("Tektite", site))

  ylim <- nmds_metab_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS2) %>%
    simplify
  xlim <- nmds_metab_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS1) %>%
    simplify
  nmds_lims <- list(ylim = c(round(min(ylim), digits = 4), round(max(ylim), digits = 4)) * 1.1,
                    xlim = c(round(min(xlim), digits = 4), round(max(xlim), digits = 4)) * 1.1)
  g3 <- (ggplot(data = nmds_metab_df,
               aes(x = NMDS1, y = NMDS2))
  + theme(text = element_text(size=14))
  + geom_point(data = (nmds_metab_df %>%
                                  dplyr::filter(type == "sample") %>%
                                  droplevels),
                        aes(fill = site, shape = addNA(sampling_time)), color = "black",
                        size = 5, stroke = 2, alpha = 1)
    + geom_text(aes(label = text_label, x = NMDS1, vjust = "outward", hjust = "center", y = NMDS2),
                position = "dodge",
                check_overlap = TRUE, colour = "grey10", fontface = "bold")
  + scale_shape_manual(values = c(22, 21, 23), labels = c(sampling_time_lookup, "NA"), breaks = c(names(sampling_time_lookup), NA))
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
  print(g3)
  if(!any(grepl("nmds", list.files(projectpath, pattern = "usvi_metabolomics.*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_metabolomics_nmds-", Sys.Date(), ".png"),
           g3,
           width = 10, height = 8, units = "in")
  }
}
  
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




# if(!any(grepl("seawater", list.files(projectpath, pattern = "usvi_nmds.*.png")))){
#   ggsave(paste0(projectpath, "/", "usvi_nmds_seawater-", Sys.Date(), ".png"),
#          g1,
#          width = 10, height = 10, units = "in")
# }


# Correlations between DNA profile and metabolomics -----------------------

drop <- c("CINAR_BC_73", "CINAR_BC_43", "CINAR_BC_105", "CINAR_BC_69")
usvi_metab_df <- usvi_metabolomics.df %>%
  dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  tibble::column_to_rownames(var = "metab_deriv_label") %>%
  as.matrix(.) %>%
    vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
                   autotransform = TRUE) %>%
  # as.matrix(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(var = "sample_id") %>% droplevels
  as.matrix(.)

usvi_asv_df <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, metab_deriv_label) %>%
                         droplevels),
                   by = join_by(sample_id)) %>%
    dplyr::select(-sample_id) %>%
    tibble::column_to_rownames(var = "metab_deriv_label") %>%
  t() %>%
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
  dplyr::filter(sample %in% colnames(usvi_metab_df)) %>%
  tibble::column_to_rownames(var = "sample") %>%
  tidyr::drop_na(.) %>%
  vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
                 autotransform = TRUE) %>%
  # as.matrix(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(var = "sample_id") %>% droplevels
  as.matrix(.)

usvi_mantel_res <- vegan::mantel(usvi_asv_df, usvi_metab_df, permutations = 999, parallel = nthreads,
                                 # method = "pearson") 
                                 method = "spearman")

# 
# 
# #correlate the bray-curtis NMDS of seawater samples with the FCM measurements?
# #FCM sample from Metab_280 had weird Prochlorococcus, Picoeukaryotes, and unpigmented cells
# usvi_sw_fcm_df <- ps_usvi %>%
#   phyloseq::sample_data(.) %>%
#   tibble::as_tibble(rownames = "sample_id") %>%
#   dplyr::filter(sample_id %in% rownames(usvi_asv.tbl)) %>%
#   dplyr::filter(!grepl("Metab_280", sample_id)) %>%
#   # dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
#   #   dplyr::select(sample_id, sampling_time, sampling_day, site, contains("fcm")) %>%
#   #   dplyr::select(!c(contains("label"), contains("dna_"))) %>%
#   tibble::column_to_rownames(., var = "sample_id") %>%
#   #   droplevels
#   # usvi_sw_fcm_df <- meta.seawater[!grepl("Metab_280", rownames(meta.seawater)),] %>%
#   # usvi_sw_fcm_df <- meta.seawater %>%
#   dplyr::select(starts_with("fcm_")) %>%
#   tidyr::drop_na(.) %>%
#   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
#                  autotransform = TRUE) %>%
#   as.matrix(.) %>%
#   as.data.frame(.) %>%
#   tibble::rownames_to_column(var = "sample_id") %>%
#   droplevels
# 
# usvi_sw_dist_df <- vegan::vegdist(usvi_asv.tbl, distance = "horn", binary = FALSE, upper = TRUE,
#                                   autotransform = TRUE) %>%
#   as.matrix(.) %>%
#   as.data.frame(.) %>%
#   dplyr::select(usvi_sw_fcm_df[["sample_id"]]) %>%
#   tibble::rownames_to_column(var = "sample_id") %>%
#   dplyr::filter(sample_id %in% usvi_sw_fcm_df[["sample_id"]]) %>%
#   droplevels %>%
#   tibble::column_to_rownames(var = "sample_id")
# usvi_sw_fcm_df <- usvi_sw_fcm_df %>%
#   tibble::column_to_rownames(var = "sample_id")

usvi_mantel_res <- vegan::mantel(usvi_sw_fcm_df, usvi_sw_dist_df, permutations = 999, parallel = nthreads,
                                 # method = "pearson") 
                                 method = "spearman")

#so the statistic ranges from 0.33 to 0.46 but the "signif = 0.001


#make a NMDS using the FCM measurements and see if we have anything significant.
nmds.fcm <- ps_usvi %>%
  phyloseq::sample_data(.) %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  dplyr::filter(sample_id %in% rownames(usvi_asv.tbl)) %>%
  dplyr::filter(!grepl("Metab_280", sample_id)) %>%
  tibble::column_to_rownames(., var = "sample_id") %>%
  dplyr::select(starts_with("fcm_")) %>%
  tidyr::drop_na(.) %>%
  vegan::decostand(., method = "normalize") %>%
  # dplyr::mutate(across(everything(), scale)) %>%
  vegan::metaMDS(., distance = "horn", autotransform = FALSE)


nmds_fcm_df <- ggplot2::fortify(nmds.fcm) %>%
  dplyr::filter(score == "sites") %>%
  dplyr::mutate(type = "sample") %>%
  dplyr::select(-score) %>%
  dplyr::left_join(., (metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, replicate) %>%
                         droplevels),
                   by = c("label" = "sample_id")) %>%
  droplevels

{
  ylim <- nmds_fcm_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS2) %>%
    simplify
  xlim <- nmds_fcm_df %>%
    dplyr::filter(!grepl("ASV", type)) %>%
    dplyr::select(NMDS1) %>%
    simplify
  nmds_lims <- list(ylim = c(round(min(ylim), digits = 4), round(max(ylim), digits = 4)) * 1.1,
                    xlim = c(round(min(xlim), digits = 4), round(max(xlim), digits = 4)) * 1.1)
  g3 <- ggplot(data = nmds_fcm_df,
               aes(x = NMDS1, y = NMDS2))
  g3 <- g3 + theme(text = element_text(size=14))
  g3 <- g3 + geom_point(data = (nmds_fcm_df %>%
                                  dplyr::filter(type == "sample") %>%
                                  droplevels),
                        aes(fill = site, shape = sampling_time), color = "black",
                        size = 5, stroke = 2, alpha = 1)
  g3 <- g3 + scale_shape_manual(values = c(22, 21), labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  g3 <- g3 + scale_fill_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
  g3 <- g3 + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                   panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                   panel.grid = element_blank(),
                   legend.position = "bottom",
                   legend.key = element_blank(),
                   legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                   legend.text = element_text(size = 12, colour = "grey30"))
  g3 <- g3 + guides(color = "none",
                    fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                        override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                    shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                         override.aes = list(color = "black", stroke = 1, size = 2)))
  g3 <- g3 + coord_cartesian(ylim = nmds_lims[["ylim"]],
                             xlim = nmds_lims[["xlim"]],
                             expand = TRUE)
  g3
}
