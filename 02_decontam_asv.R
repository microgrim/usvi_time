# 02_asv_analysis.R

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

library(tidyverse)
library(phyloseq)
library(decontam)
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


P_summary_taxonomy_bar <- function(dataset, x, y, z, facet_form) {
  
  
  if(exists("facet_form", inherits = TRUE)){
    # facet_rows <- stringr::str_split_i(facet_form, "~", 1)
    # facet_cols <- stringr::str_split_i(facet_form, "~", 2)
    if(grepl("\\.", facet_form)){
      facet_form <- deparse(substitute(facet_form))
    }
  } else {
    facet_form <- deparse(substitute(z)) %>%
      paste0(".~", .)
  }
  title_x <- deparse(substitute(x))
  
  g <- ggplot(data = dataset, aes(y = {{y}},
                                  x = {{x}},
                                  fill = {{z}}))
  g <- g + geom_col(width = 0.90, show.legend = TRUE, position = position_dodge2(padding = 0.1, preserve = "total"))
  g <- g + geom_col(color = "black", width = 0.90, show.legend = FALSE, position = position_dodge2(padding = 0.1, preserve = "total"))
  g <- g + theme_bw() 
  g <- g + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  g <- g + scale_fill_viridis(discrete = TRUE, option = "turbo")
  g <- g + guides(fill = guide_legend(direction = "vertical", ncol = 1))
  g <- g + facet_grid(scales = "free", space = "free",
                      # labeller = labeller(coral_species = coral_lookup),
                      # rows = vars({{facet_rows}}), cols = vars({{facet_cols}}))
                      eval(facet_form))
  g <- g + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                 strip.text.x = element_text(angle = 0),
                 strip.text.y = element_text(angle = -90))
  # g <- g + facet_wrap(scales = "free", ncol = 3, 
  #                     # labeller = labeller(coral_species = coral_lookup),
  #                     eval(facet_form))
  g <- g + labs(x = eval(title_x),
                y = "Relative abundance (%)")
  g <- g + coord_cartesian(ylim = c(0,100), expand = FALSE)
  g <- g + guides(fill = guide_legend(order = 2, ncol = 1, title = "Phylum", direction = "vertical",
                                      override.aes = list(color = "black", stroke = 1)))
  g
}


# Read in metadata --------------------------------------------------------

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
                                show_col_types = FALSE,
                                num_threads = nthreads)
} else {
  #note that the following samples were sequenced in a separate run:
  #Metab_180, Metab_182, Metab_199, Metab_219, Metab_224, Metab_231, Metab_233, Metab_235
  # c("Metab_180", "Metab_182", "Metab_199", "Metab_219", "Metab_224", "Metab_231", "Metab_233", "Metab_235")
  
sample_metadata <- readr::read_delim(paste0(projectpath, "/", "usvi_metadata.txt"),
                                     delim = "\t",
                                     quote = "",
                                     col_names = TRUE,
                                     show_col_types = FALSE,
                                     num_threads = nthreads) %>%
# sample_metadata <- metadata %>%
  dplyr::mutate(batch = dplyr::case_when((sample_ID %in% c("Metab_180", "Metab_182", "Metab_199", "Metab_219", "Metab_224", "Metab_231", "Metab_233", "Metab_235")) ~ "1",
                                         .default = "2")) %>%
  dplyr::mutate(dna_conc = dplyr::case_when(!grepl("seawater", sample_type) & is.na(dna_conc) ~ 0.001,
                                            .default = dna_conc)) %>%
  dplyr::mutate(is.neg = dplyr::case_when((sample_type == "seawater") ~ FALSE,
                                          .default = TRUE))
#maybe randomize potential dna concentrations for the 8 samples extracted during the Spike-in project that did not get DNA concentration quantified
if(any(is.na(sample_metadata$dna_conc))){
  seawater_dna_conc <- sample_metadata %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::select(dna_conc) %>%
    droplevels %>%
    unlist
  
  seawater_dna_prob <- quantile(seawater_dna_conc, probs = seq(0.25, 1, 0.25), na.rm = TRUE, names = TRUE)
  
  set.seed(48105)
  temp_conc <- data.frame(sample_id = c("Metab_180", "Metab_182", "Metab_199", "Metab_219", "Metab_224", "Metab_231", "Metab_233", "Metab_235"),
                          dna_conc = F_random_sampler(8, x = seawater_dna_conc, seawater_dna_prob)) %>%
    dplyr::mutate(dna_conc = as.numeric(dna_conc))
  sample_metadata <- sample_metadata %>%
    dplyr::left_join(., temp_conc, by = join_by("sample_ID" == "sample_id")) %>%
    dplyr::mutate(dna_conc = across(starts_with("dna_conc")) %>% purrr::reduce(coalesce)) %>%
    dplyr::select(-ends_with(c(".x", ".y"))) %>%

    droplevels
}

metadata <- sample_metadata %>%
  dplyr::select(contains("sampl"), site) %>%
  dplyr::mutate(site = dplyr::case_when(grepl("seawater", sample_type) ~ site,
                                        .default = NA)) %>%
  dplyr::mutate(sample_type = dplyr::case_when(grepl("seawater", sample_type) ~ "seawater",
                                               # (site %in% temporal_lookup[["site"]]) ~ site,
                                               grepl("pcr", sample_type) ~ "control_pcr",
                                               grepl("extraction", sample_type) ~ "control_extraction",
                                               .default = "control_seq")) %>%
  dplyr::mutate(sample_type = factor(sample_type, levels = temporal_lookup[["sample_type"]]),
                sampling_time = factor(sampling_time, levels = temporal_lookup[["sampling_time"]]),
                sampling_date = factor(sampling_date, levels = temporal_lookup[["sampling_date"]]),
                sampling_day = factor(sampling_day, levels = temporal_lookup[["sampling_day"]])) %>%
  dplyr::arrange(sample_type, site, sampling_time) %>%
  dplyr::mutate(sample_id = factor(sample_ID, levels = unique(.[["sample_ID"]]))) %>%
  droplevels

metadata <- metadata %>%
  dplyr::left_join(., bind_rows(
    (metadata %>%
       dplyr::filter(grepl("seawater", sample_type)) %>%
       droplevels %>%
       dplyr::group_by(sampling_time, sampling_day, site) %>%
       dplyr::mutate(replicate = rep(c(1:3))) %>%
       dplyr::ungroup(.) %>%
       dplyr::select(sample_ID, replicate) %>%
       droplevels),
    (sample_metadata %>%
       dplyr::filter(grepl("extraction", sample_type)) %>%
       droplevels %>%
       dplyr::arrange(sample_type) %>%
       dplyr::select(sample_ID) %>%
       dplyr::left_join(., metadata %>%
                          dplyr::select(sample_ID, sample_type)) %>%
       dplyr::group_by(sample_type) %>%
       dplyr::mutate(replicate = seq(1, length(.[["sample_type"]]))) %>%
       dplyr::ungroup(.) %>%
       dplyr::select(sample_ID, replicate) %>%
       droplevels),
    (sample_metadata %>%
       dplyr::filter(grepl("pcr", sample_type)) %>%
       droplevels %>%
       dplyr::arrange(sample_type) %>%
       dplyr::select(sample_ID) %>%
       dplyr::left_join(., metadata %>%
                          dplyr::select(sample_ID, sample_type)) %>%
       dplyr::group_by(sample_type) %>%
       dplyr::mutate(replicate = seq(1, length(.[["sample_type"]]))) %>%
       dplyr::ungroup(.) %>%
       dplyr::select(sample_ID, replicate) %>%
       droplevels)),
    by = join_by(sample_ID)) 

metadata <- metadata %>%
# temp_df2 <- metadata %>%
  dplyr::left_join(., bind_rows((metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         droplevels %>%
                           dplyr::group_by(sampling_time, sampling_day, site) %>%
                           dplyr::ungroup(.) %>%
                           dplyr::mutate(sample_order = dplyr::case_when((sampling_day == "Day1") ~ replicate * 1,
                                                                         (sampling_day == "Day2") ~ replicate + 3,
                                                                         (sampling_day == "Day3") ~ replicate + 6,
                                                                         (sampling_day == "Day4") ~ replicate + 9,
                                                                         (sampling_day == "Day5") ~ replicate + 12,
                                                                         .default = NA)) %>%
                           dplyr::mutate(sample_order_all = dplyr::case_when((sampling_day == "Day1" & sampling_time == "dawn") ~ replicate * 1,
                                                                             (sampling_day == "Day2" & sampling_time == "dawn") ~ replicate + 6,
                                                                             (sampling_day == "Day3" & sampling_time == "dawn") ~ replicate + 12,
                                                                             (sampling_day == "Day4" & sampling_time == "dawn") ~ replicate + 18,
                                                                             (sampling_day == "Day5" & sampling_time == "dawn") ~ replicate + 24,
                                                                             (sampling_day == "Day1" & sampling_time == "peak_photo") ~ replicate + 3,
                                                                             (sampling_day == "Day2" & sampling_time == "peak_photo") ~ replicate + 9,
                                                                             (sampling_day == "Day3" & sampling_time == "peak_photo") ~ replicate + 15,
                                                                             (sampling_day == "Day4" & sampling_time == "peak_photo") ~ replicate + 21,
                                                                             (sampling_day == "Day5" & sampling_time == "peak_photo") ~ replicate + 27,
                                                                             .default = NA)) %>%
                           dplyr::select(sample_ID, contains("sample_order")) %>%
                         # dplyr::group_by(site) %>%
                         # dplyr::mutate(sample_order = seq(1, length(interaction(replicate, sampling_time, sampling_day)))) %>%
                         # dplyr::ungroup(.) %>%
                         # dplyr::select(sample_ID, sample_order) %>%
                         droplevels),
                         (metadata %>%
                            dplyr::filter(!grepl("seawater", sample_type)) %>%
                            droplevels %>%
                            dplyr::arrange(sample_type, replicate) %>%
                            dplyr::group_by(sample_type) %>%
                            dplyr::mutate(sample_order_all = seq(1, length(interaction(replicate, sample_type)))) %>%
                            dplyr::mutate(sample_order = sample_order_all) %>%
                            dplyr::ungroup(.) %>%
                            dplyr::select(sample_ID, contains("sample_order")) %>%
                            droplevels)),
                   by = join_by(sample_ID))

  readr::write_delim(metadata, paste0(projectpath, "/", "usvi_metadata_tidy.txt"),
                     delim = "\t", col_names = TRUE, num_threads = nthreads)
}

# Import color scheme -----------------------------------------------------


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
  droplevels



# Make a phyloseq object --------------------------------------------------

if(!exists("ps_usvi", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))){
    ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))
  } else {
    ps_usvi <- phyloseq::phyloseq(phyloseq::otu_table((usvi_prok_asvs.df %>%
                                                         tidyr::pivot_wider(., id_cols = c("asv_id"),
                                                                            names_from = "sample_ID",
                                                                            values_from = "counts") %>%
                                                         tibble::column_to_rownames(var = "asv_id")),
                                                      taxa_are_rows=TRUE),
                                  phyloseq::sample_data(sample_metadata %>%
                                                          tibble::column_to_rownames(var = "sample_ID")),
                                  phyloseq::tax_table(usvi_prok_asvs.taxa %>%
                                                        dplyr::select(-sequence) %>%
                                                        droplevels %>%
                                                        tibble::column_to_rownames(var = "asv_id") %>%
                                                        as.matrix))
    # if(!file.exists(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))){
    readr::write_rds(ps_usvi, paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"), compress = "gz")
  }
}


#

# Decontamination ---------------------------------------------------------


set_thresholds <- c(0.1, 0.5) %>%
  setNames(., c("freq", "prev"))


##
# starting with ASV level is not broad enough
#use species-level abundances, and amplicon/DNA concentrastions


ps_usvi_samples_species <- phyloseq::phyloseq(phyloseq::otu_table((usvi_prok_asvs.df %>%
                                                                     tidyr::pivot_wider(., id_cols = c("asv_id"),
                                                                                        names_from = "sample_ID",
                                                                                        values_from = "counts") %>%
                                                                     tibble::column_to_rownames(var = "asv_id")),
                                                                  taxa_are_rows=TRUE),
                                              phyloseq::sample_data(sample_metadata %>%
                                                                      tibble::column_to_rownames(var = "sample_ID")),
                                              phyloseq::tax_table(usvi_prok_filled.taxa.df %>%
                                                                    dplyr::select(-sequence) %>%
                                                                    droplevels %>%
                                                                    tibble::column_to_rownames(var = "asv_id") %>%
                                                                    as.matrix)) %>%
  phyloseq::subset_samples(., grepl("seawater|pcr|extraction", sample_type)) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% # remove ASVs not present in any samples
  phyloseq::tax_glom(., taxrank = "Species")


decon_usvi_species_combo.df <- decontam::isContaminant(ps_usvi_samples_species, 
                                                       # method = "minimum",
                                                       method = "combined",
                                                       batch = sample_data(ps_usvi_samples_species)$batch,
                                                       threshold = set_thresholds[["freq"]],
                                                       conc = "dna_conc",
                                                       # conc = "product_conc",
                                                       neg = "is.neg")
# table(decon_usvi_species_combo.df$contaminant)
# FALSE  TRUE 
# 1170    33  #threshold = 0.1, method = "minimum", using dna_conc
# 1194     9  #threshold = 0.1, method = "combined", using dna_conc

# rownames(decon_usvi_species_combo.df)[(which(decon_usvi_species_combo.df$contaminant))]
# "ASV_00033" "ASV_00081" "ASV_00090" "ASV_00267" "ASV_00434" "ASV_01660" "ASV_03083" "ASV_05175" "ASV_06080"
contam_species_idx <- decon_usvi_species_combo.df[["contaminant"]] %>%
  setNames(., rownames(decon_usvi_species_combo.df))


potential_contam_species.taxa <- decon_usvi_species_combo.df %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  dplyr::mutate(contam_freq = dplyr::case_when(p.freq < set_thresholds[["freq"]] ~ "freq", 
                                               .default = NA)) %>%
  dplyr::mutate(contam_prev = dplyr::case_when(p.prev < set_thresholds[["prev"]] ~ "prev",
                                               .default = NA)) %>%
  dplyr::mutate(potential_contaminant = dplyr::case_when(!is.na(contam_freq) | !is.na(contam_prev) ~ TRUE,
                                                         .default = FALSE)) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df, by = join_by(asv_id)) %>%
  droplevels

{
  ps.pa <- transform_sample_counts(ps_usvi_samples_species, function(abund) 1*(abund>0))
  ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_type != "seawater", ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type == "seawater", ps.pa)
  # Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=decon_usvi_species_combo.df$contaminant)
  # ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  #   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
  df.pa.potential <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      potential_contaminant=potential_contam_species.taxa$potential_contaminant,
                      contaminant=potential_contam_species.taxa$contaminant)
  # ggplot(data=df.pa.potential, aes(x=pa.neg, y=pa.pos, color=contaminant, fill = potential_contaminant)) + geom_point(shape = 21, size = 2) +
  #   scale_fill_manual(values = c("forestgreen", "salmon")) + scale_color_manual(values = c("cyan", "magenta")) +
  #   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
  }


#not contaminant: any species that were absent in negative controls
potential_contam_species.taxa <- decon_usvi_species_combo.df %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  dplyr::mutate(contam_freq = dplyr::case_when(p.freq < set_thresholds[["freq"]] ~ "freq", 
                                               .default = NA)) %>%
  dplyr::mutate(contam_prev = dplyr::case_when(p.prev < set_thresholds[["prev"]] ~ "prev",
                                               .default = NA)) %>%
  droplevels %>%
  dplyr::left_join(., data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg)) %>%
                     tibble::rownames_to_column(var = "asv_id"), by = join_by(asv_id)) %>%
  dplyr::mutate(potential_contaminant = dplyr::case_when(pa.neg == 0 & is.na(contam_freq) & is.na(contam_prev) ~ 0,
                                                         # pa.neg == 0 & (is.na(contam_freq) | is.na(contam_prev)) ~ 1,
                                                         pa.neg > 0 & contaminant ~ 2,
                                                         pa.neg > 0 & (!is.na(contam_freq) | !is.na(contam_prev)) ~ 1,
                                                         pa.neg == 0 & !contaminant ~ 0,
                                                         .default = 0)) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df, by = join_by(asv_id)) %>%
  dplyr::relocate(asv_id, pa.pos, pa.neg, potential_contaminant, contaminant, Domain:Species) %>%
  dplyr::filter(potential_contaminant > 0 | contaminant | !is.na(contam_freq) | !is.na(contam_prev)) %>%
  droplevels
  


# Plot contaminants -------------------------------------------------------



#plot the relative abundances of these potential contaminant ASVs alongside legitimate taxa:
#" a noncontaminant sequence feature for which frequency is expected to be independent of 
# the input DNA concentration.
# a contaminant sequence feature, for which frequency is expected to be inversely proportional
# to input DNA concentration, as contaminating DNA will make up a larger fraction of the total DNA in samples with very little total DNA." 


#add three ASVs we know are legitimate taxa:
biol_asv_idx <- c("ASV_00001", "ASV_00002", "ASV_00007")

species_contam_both <- potential_contam_species.taxa %>%
  dplyr::filter(!is.na(contam_freq) & !is.na(contam_prev)) %>%
  dplyr::select(asv_id) %>%
  droplevels %>%
  unlist %>%
  as.character

species_contam_only_prev <- potential_contam_species.taxa %>%
  dplyr::filter(contaminant) %>%
  dplyr::filter(!is.na(contam_prev)) %>%
  dplyr::filter(!(asv_id %in% species_contam_both)) %>%
  tidyr::drop_na(p) %>%
  dplyr::select(asv_id) %>%
  droplevels %>%
  unlist %>%
  as.character

species_contam_only_freq <- potential_contam_species.taxa %>%
  dplyr::filter(contaminant) %>%
  dplyr::filter(!is.na(contam_freq)) %>%
  dplyr::filter(!(asv_id %in% species_contam_both)) %>%
  tidyr::drop_na(p) %>%
  dplyr::select(asv_id) %>%
  droplevels %>%
  unlist %>%
  as.character

species_contam_neither <- setdiff(potential_contam_species.taxa[["asv_id"]], c(species_contam_only_freq, species_contam_only_prev,species_contam_both))

{
  temp_df <- ps_usvi_samples_species %>% 
    phyloseq::transform_sample_counts(., function(OTU) OTU/sum(OTU)) %>%
    phyloseq::otu_table(.) %>%
    as.data.frame %>%
    tibble::as_tibble(rownames = "asv_id") %>%
    dplyr::filter(asv_id %in% c(potential_contam_species.taxa[["asv_id"]], biol_asv_idx)) %>%
    tibble::column_to_rownames(var = "asv_id") %>%
    t() %>%
    as.matrix()
  
  temp_df2 <- sample_data(ps_usvi_samples_species) %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    dplyr::select(sample_id, is.neg, contains("conc"), contains("fcm")) %>%
    dplyr::select(-c("fcm_label")) %>%
    tibble::column_to_rownames(var = "sample_id")
  
  logd <- log(seq(min(temp_df2[, "dna_conc"], na.rm = TRUE), max(temp_df2[, "dna_conc"], na.rm = TRUE), length.out=1000))
  logc <- log(seq(min(temp_df2[, "product_conc"], na.rm = TRUE), max(temp_df2[, "product_conc"], na.rm = TRUE), length.out=1000))
  
  contam_model <- NULL
  asv <- NULL
  for(asv in colnames(temp_df)){
    newdata <- data.frame(asv_id = asv, 
                          logd = logd, conc_dna = exp(logd),
                          logc = logc, conc_prod = exp(logc))
    freq <- temp_df[, asv]
    neg <- temp_df2[, "is.neg"]
    conc_dna <- temp_df2[, "dna_conc"]
    conc_prod <- temp_df2[, "product_conc"]
    
    df <- data.frame(logc=log(conc_prod), logd = log(conc_dna), logf=log(freq))
    df <- df[freq>0 & (!neg | is.na(neg)),]
    if(nrow(df) > 0 & sum(freq>0)>1) {
      lm2 <- lm(logf~offset(-1*logd), data=df)
      lm1 <- lm(logf~offset(-1*logc), data=df)
      lm0 <- lm(logf~1, data=df)
      newdata$contam.dna <- exp(predict(lm2, newdata=newdata))
      newdata$contam.pcr <- exp(predict(lm1, newdata=newdata))
      newdata$non.contam <- exp(predict(lm0, newdata=newdata))
    } else {
      newdata$contam.dna <- NA
      newdata$contam.pcr <- NA
      newdata$non.contam <- NA
    }
    contam_model[[asv]] <- newdata
  }
  contam_model <- do.call(rbind, contam_model)
}

plot_potential_contam.df <- ps_usvi_samples_species %>% 
  phyloseq::transform_sample_counts(., function(OTU) OTU/sum(OTU)) %>%
  phyloseq::otu_table(.) %>%
  as.data.frame %>%
  tibble::as_tibble(rownames = "asv_id") %>%
  dplyr::filter(asv_id %in% c(potential_contam_species.taxa[["asv_id"]], biol_asv_idx)) %>%
  tidyr::pivot_longer(., cols = !c("asv_id"),
                      names_to = "sample_id",
                      values_to = "taxon_abundance") %>%
  dplyr::left_join(., sample_data(ps_usvi_samples_species) %>%
                     tibble::as_tibble(rownames = "sample_id") %>%
                     dplyr::select(sample_id, sample_type, contains("conc"), site, contains("fcm"), batch,  is.neg) %>%
                     dplyr::mutate(modeled_dna_conc = dplyr::case_when((sample_id %in% c("Metab_180", "Metab_182", "Metab_199", "Metab_219", "Metab_224", "Metab_231", "Metab_233", "Metab_235")) ~ "modeled",
                                                                       sample_type != "seawater" ~ "modeled",
                                                                       .default = "measured")) %>%
                     droplevels, by = "sample_id") %>%
  dplyr::mutate(contaminant = dplyr::case_when(asv_id %in% biol_asv_idx ~ FALSE,
                                               .default = TRUE)) %>%
  dplyr::mutate(across(c(is.neg, contaminant), ~as.factor(.x))) %>%
  dplyr::left_join(., potential_contam_species.taxa %>%
                     dplyr::select(asv_id, contam_freq, contam_prev), by = join_by(asv_id), relationship = "many-to-many")

species_groupings <- list(biol_asv_idx, species_contam_only_prev, species_contam_only_freq, species_contam_both, species_contam_neither) %>%
  setNames(., c("biol_asv_idx", "species_contam_only_prev", "species_contam_only_freq", "species_contam_both", "species_contam_neither"))


#set up the vectors of ASVs flagegd as potential contaminants
{
# potential_contam_idx <- union(rownames(decon_usvi_prev.df)[(which(decon_usvi_prev.df$contaminant))], rownames(decon_usvi_freq.df)[(which(decon_usvi_freq.df$contaminant))])
# 
# potential_contam_usvi_asvs.taxa <- decon_usvi_prev.df %>%
#   tibble::rownames_to_column(var = "asv_id") %>%
#   dplyr::select(asv_id, prev, p.prev) %>%
#   dplyr::left_join(., decon_usvi_freq.df %>%
#                      tibble::rownames_to_column(var = "asv_id") %>%
#                      dplyr::select(asv_id, freq, p.freq), 
#                    by = join_by(asv_id)) %>%
#   dplyr::left_join(., decon_usvi_combo.df %>%
#                      tibble::rownames_to_column(var = "asv_id") %>%
#                      dplyr::select(asv_id, p, contaminant),
#                    by = join_by(asv_id)) %>%
#   dplyr::mutate(contam_freq = dplyr::case_when(p.freq < set_thresholds[["freq"]] ~ "freq", 
#                                                .default = NA)) %>%
#   dplyr::mutate(contam_prev = dplyr::case_when(p.prev < set_thresholds[["prev"]] ~ "prev",
#                                                .default = NA)) %>%
#     dplyr::filter(asv_id %in% potential_contam_idx) %>%
#   dplyr::left_join(., usvi_prok_filled.taxa.df, by = join_by(asv_id)) %>%
#   droplevels

#   asv_contam_either <- potential_contam_usvi_asvs.taxa %>%
#     dplyr::filter(asv_id %in% potential_contam_idx) %>%
#     dplyr::filter(!is.na(contam_freq) & !is.na(contam_prev)) %>%
#     # dplyr::filter(!is.na(contam_freq) | !is.na(contam_prev)) %>%
#     dplyr::filter(contaminant) %>%
#     tidyr::drop_na(p) %>%
#     dplyr::select(asv_id) %>%
#     droplevels %>%
#     unlist %>%
#     as.character
#   
#   asv_contam_only_prev <- potential_contam_usvi_asvs.taxa %>%
#     dplyr::filter(contaminant) %>%
#     dplyr::filter(!is.na(contam_prev)) %>%
#     dplyr::filter(!(asv_id %in% asv_contam_either)) %>%
#     tidyr::drop_na(p) %>%
#     dplyr::select(asv_id) %>%
#     droplevels %>%
#     unlist %>%
#     as.character
#   
# asv_contam_only_freq <- potential_contam_usvi_asvs.taxa %>%
#   dplyr::filter(contaminant) %>%
#   dplyr::filter(!is.na(contam_freq)) %>%
#   dplyr::filter(!(asv_id %in% asv_contam_either)) %>%
#   tidyr::drop_na(p) %>%
#   dplyr::select(asv_id) %>%
#   droplevels %>%
#   unlist %>%
#   as.character
# 
# # asv_contam_neither <- setdiff(potential_contam_idx, c(asv_contam_only_prev, asv_contam_only_freq, asv_contam_either))
#   
# asv_contam_neither <- potential_contam_usvi_asvs.taxa %>%
#   dplyr::filter(asv_id %in% potential_contam_idx) %>%
#   dplyr::filter(!contaminant) %>%
#   # dplyr::filter(is.na(contam_freq) & is.na(contam_prev)) %>%
#   dplyr::select(asv_id) %>%
#   droplevels %>%
#   unlist %>%
#   as.character
}


#plot the models:
{
  # temp_df <- ps_usvi_samples %>% 
  #   phyloseq::transform_sample_counts(., function(OTU) OTU/sum(OTU)) %>%
  #   phyloseq::otu_table(.) %>%
  #   as.data.frame %>%
  #   tibble::as_tibble(rownames = "asv_id") %>%
  #   dplyr::filter(asv_id %in% c(potential_contam_idx, biol_asv_idx)) %>%
  #   tibble::column_to_rownames(var = "asv_id") %>%
  #   t() %>%
  #   as.matrix()
  # 
  # temp_df2 <- sample_data(ps_usvi_samples) %>%
  #   tibble::as_tibble(rownames = "sample_id") %>%
  #   # dplyr::left_join(., temp_conc, by = join_by(sample_id)) %>%
  #   # dplyr::mutate(dna_conc = across(starts_with("dna_conc")) %>% purrr::reduce(coalesce)) %>%
  #   # dplyr::select(-ends_with(c(".x", ".y"))) %>%
  #   # dplyr::mutate(across(!c("sample_id", "is.neg"), ~dplyr::case_when((.x > 0) ~ log(.x),
  #   #                                                      .default = NA))) %>%
  #   dplyr::select(sample_id, is.neg, contains("conc"), contains("fcm")) %>%
  #   dplyr::select(-c("fcm_label")) %>%
  #   tibble::column_to_rownames(var = "sample_id")
  # 
  # logd <- log(seq(min(temp_df2[, "dna_conc"], na.rm = TRUE), max(temp_df2[, "dna_conc"], na.rm = TRUE), length.out=1000))
  # logc <- log(seq(min(temp_df2[, "product_conc"], na.rm = TRUE), max(temp_df2[, "product_conc"], na.rm = TRUE), length.out=1000))
  # 
  # contam_model <- NULL
  # asv <- NULL
  # # for(asv in c("ASV_00078", "ASV_00139")){
  # # for(asv in colnames(temp_df)[6:8]){
  # for(asv in colnames(temp_df)){
  #   newdata <- data.frame(asv_id = asv, 
  #                         logd = logd, conc_dna = exp(logd),
  #                         logc = logc, conc_prod = exp(logc))
  #   freq <- temp_df[, asv]
  #   neg <- temp_df2[, "is.neg"]
  #   conc_dna <- temp_df2[, "dna_conc"]
  #   conc_prod <- temp_df2[, "product_conc"]
  #   
  #   df <- data.frame(logc=log(conc_prod), logd = log(conc_dna), logf=log(freq))
  #   # df <- data.frame(logc=conc_prod, 
  #   #                  # logd = conc_dna, 
  #   #                  logf=freq)
  #   df <- df[freq>0 & (!neg | is.na(neg)),]
  #   # df <- df[!neg | is.na(neg),]
  #   # df <- df[freq>0,]
  #   # df <- log(df)
  #   if(nrow(df) > 0 & sum(freq>0)>1) {
  #     lm2 <- lm(logf~offset(-1*logd), data=df)
  #     lm1 <- lm(logf~offset(-1*logc), data=df)
  #     lm0 <- lm(logf~1, data=df)
  #     newdata$contam.dna <- exp(predict(lm2, newdata=newdata))
  #     newdata$contam.pcr <- exp(predict(lm1, newdata=newdata))
  #     newdata$non.contam <- exp(predict(lm0, newdata=newdata))
  #   } else {
  #     newdata$contam.dna <- NA
  #     newdata$contam.pcr <- NA
  #     newdata$non.contam <- NA
  #   }
  #   contam_model[[asv]] <- newdata
  # }
  # contam_model <- do.call(rbind, contam_model)
}

#set up dataframe with relative abundances of potential contaminants
{
# plot_potential_contam.df <- ps_usvi_samples %>% 
#   phyloseq::transform_sample_counts(., function(OTU) OTU/sum(OTU)) %>%
#   phyloseq::otu_table(.) %>%
#   as.data.frame %>%
#   tibble::as_tibble(rownames = "asv_id") %>%
#   dplyr::filter(asv_id %in% c(potential_contam_idx, biol_asv_idx)) %>%
#   tidyr::pivot_longer(., cols = !c("asv_id"),
#                       names_to = "sample_id",
#                       values_to = "taxon_abundance") %>%
#   dplyr::left_join(., sample_data(ps_usvi_samples) %>%
#                      tibble::as_tibble(rownames = "sample_id") %>%
#                      dplyr::select(sample_id, sample_type, contains("conc"), site, contains("fcm"), batch,  is.neg) %>%
#                      dplyr::mutate(modeled_dna_conc = dplyr::case_when((sample_id %in% c("Metab_180", "Metab_182", "Metab_199", "Metab_219", "Metab_224", "Metab_231", "Metab_233", "Metab_235")) ~ "modeled",
#                                                                        sample_type != "seawater" ~ "modeled",
#                                                                        .default = "measured")) %>%
#                      # dplyr::left_join(., temp_conc, by = join_by(sample_id)) %>%
#                      # dplyr::mutate(dna_conc = across(starts_with("dna_conc")) %>% purrr::reduce(coalesce)) %>%
#                      # dplyr::select(-ends_with(c(".x", ".y"))) %>%
#                      droplevels, by = "sample_id") %>%
#                      # dplyr::rowwise(.) %>%
#                      # dplyr::mutate(dna_conc = dplyr::case_when(sample_type == "seawater" & is.na(dna_conc) ~ F_random_sampler(1, seawater_dna_conc, seawater_dna_prob),
#                      #                                           sample_type != "seawater" & is.na(dna_conc) ~ 0,
#                      #                                           .default = dna_conc)), by = "sample_id") %>%
#                      # dplyr::mutate(across(contains("conc"), ~dplyr::case_when(sample_type == "seawater" & is.na(.x) ~ F_random_sampler(1, seawater_dna_conc, seawater_dna_prob),
#                      #                                                          sample_type != "seawater" & is.na(.x) ~ 0,
#                      #                                                          .default = .x))), by = "sample_id") %>%
#   dplyr::mutate(contaminant = dplyr::case_when(asv_id %in% biol_asv_idx ~ FALSE,
#                                                .default = TRUE)) %>%
#   dplyr::mutate(across(c(is.neg, contaminant), ~as.factor(.x))) %>%
#   dplyr::left_join(., potential_contam_usvi_asvs.taxa %>%
#                      dplyr::select(asv_id, contam_freq, contam_prev), by = join_by(asv_id), relationship = "many-to-many")

# asv_groupings <- list(biol_asv_idx, asv_contam_only_prev, asv_contam_only_freq, asv_contam_either, asv_contam_neither) %>%
#   setNames(., c("biol_asv_idx", "asv_contam_only_prev", "asv_contam_only_freq", "asv_contam_either", "asv_contam_neither"))
}

# for(i in names(asv_groupings)[1]){
for(i in names(species_groupings)){
  temp_idx <- get(i, inherits = TRUE)
  
  g1_temp <- (ggplot(data = plot_potential_contam.df %>%
                       dplyr::filter(asv_id %in% temp_idx))
              + geom_point(aes(x = dna_conc, y = taxon_abundance, fill = is.neg, shape = contaminant, color = modeled_dna_conc), size = 3, show.legend = TRUE)
              + geom_line(data = contam_model %>%
                            dplyr::filter(asv_id %in% temp_idx) %>%
                            droplevels, aes(x = conc_dna, y = contam.dna), color="red", linetype="solid", show.legend = FALSE)
              + geom_line(data = contam_model %>%
                            dplyr::filter(asv_id %in% temp_idx) %>%
                            droplevels, aes(x = conc_dna, y = non.contam), color="black", linetype="dashed",  show.legend = FALSE)
              + theme_bw()
              + scale_fill_manual(values = c("pink", "forestgreen"),
                                  breaks = c(TRUE, FALSE), drop = FALSE,
                                  labels = c("negative control", "seawater"))
              + scale_color_manual(values = c("black", "grey", "black"),
                                  breaks = c("measured", "modeled", NA), drop = FALSE,
                                  labels = c("measured", "modeled", NA))
              + scale_shape_manual(values = c(21, 22), drop = FALSE,
                                   breaks = c(TRUE, FALSE),
                                   labels = c("yes", "no"))
              + facet_wrap(~asv_id, scales = "fixed")
              + scale_x_continuous(transform = scales::log10_trans(),
                                   name = "DNA concentration (ng/uL)", labels = scales::number)
              + scale_y_continuous(transform = scales::log10_trans(),
                                   limits = c(NA, 1),
                                   name = "Relative abundance", labels = scales::percent)
              + theme(panel.grid.minor.x = element_blank(),
                      panel.grid.minor.y = element_blank())
              + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sample type",  direction = "vertical", override.aes = list(color = "black", shape = 23, linetype = NA)),
                       shape = guide_legend(order =2, ncol = 1, title = "Flagged as contaminant", override.aes = list(color = "black", fill = "white", linetype = NA)), 
                       color = guide_legend(order =3, ncol = 1, title = "Modeled DNA concentrations", override.aes = list(shape = 21, fill = "forestgreen", linetype = NA)))
  )
  assign(paste0("g1_dna_", i), g1_temp, envir = .GlobalEnv)
  g1_temp <- (ggplot(data = plot_potential_contam.df %>%
                       dplyr::filter(asv_id %in% temp_idx))
              + geom_point(aes(x = product_conc, y = taxon_abundance, fill = is.neg, shape = contaminant, color = modeled_dna_conc), color = "black", size = 3, show.legend = TRUE)
              + geom_line(data = contam_model %>%
                            dplyr::filter(asv_id %in% temp_idx) %>%
                            droplevels, aes(x = conc_prod, y = contam.pcr), color="red", linetype="solid", show.legend = FALSE)
              + geom_line(data = contam_model %>%
                            dplyr::filter(asv_id %in% temp_idx) %>%
                            droplevels, aes(x = conc_prod, y = non.contam), color="black", linetype="dashed",  show.legend = FALSE)
              + theme_bw()
              + scale_fill_manual(values = c("pink", "forestgreen"),
                                  breaks = c(TRUE, FALSE), drop = FALSE,
                                  labels = c("negative control", "seawater"))
              + scale_color_manual(values = c("black", "grey", "black"),
                                   breaks = c("measured", "modeled", NA), drop = FALSE,
                                   labels = c("measured", "modeled", NA))
              + scale_shape_manual(values = c(21, 22), drop = FALSE,
                                   breaks = c(TRUE, FALSE),
                                   labels = c("yes", "no"))
              + facet_wrap(~asv_id, scales = "fixed")
              + scale_x_continuous(transform = scales::log10_trans(),
                                   name = "Amplicon Concentration (ng/uL)", labels = scales::number)
              + scale_y_continuous(transform = scales::log10_trans(),
                                   limits = c(NA, 1),
                                   name = "Relative abundance", labels = scales::percent)
              + theme(panel.grid.minor.x = element_blank(),
                      panel.grid.minor.y = element_blank())
              + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sample type",  direction = "vertical", override.aes = list(color = "black", shape = 23, linetype = NA)),
                       shape = guide_legend(order =2, ncol = 1, title = "Flagged as contaminant", override.aes = list(color = "black", fill = "white", linetype = NA)), 
                       color = guide_legend(order =3, ncol = 1, title = "Modeled DNA concentrations", override.aes = list(shape = 21, fill = "forestgreen", linetype = NA)))
  )
  assign(paste0("g1_pcr_", i), g1_temp, envir = .GlobalEnv)
  rm(g1_temp)
  rm(temp_idx)
}

g_biol_asv_idx <- apropos("^g1_.*biol_asv_idx$", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `%+%`)  + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Not contaminant", tag_levels = "A")
# g_asv_contam_only_prev <- apropos("^g1_.*asv_contam_only_prev$", mode = "list") %>%
g_species_contam_only_prev <- apropos("^g1_.*species_contam_only_prev$", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `%+%`) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Flagged by prevalence", tag_levels = "A")
# g_asv_contam_either <- apropos("^g1_.*asv_contam_either$", mode = "list") %>%
g_species_contam_both <- apropos("^g1_.*species_contam_both$", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `%+%`) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Flagged in both prevalence and frequency", tag_levels = "A")
# g_asv_contam_neither <- apropos("^g1_.*asv_contam_neither$", mode = "list") %>%
g_species_contam_neither <- apropos("^g1_.*species_contam_neither$", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `%+%`) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Flagged in one but not the other", tag_levels = "A")
# g_asv_contam_only_freq <- apropos("^g1_.*asv_contam_only_freq$", mode = "list") %>%
g_species_contam_only_freq <- apropos("^g1_.*species_contam_only_freq$", mode = "list") %>%
  lapply(., get) %>%
  purrr::reduce(., `%+%`) + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(title = "Flagged by frequency", tag_levels = "A")


# ggsave(paste0(projectpath, "/", "usvi_asvs_potential_contaminants_bio-", Sys.Date(), ".png"),
ggsave(paste0(projectpath, "/", "usvi_species_potential_contaminants_bio-", Sys.Date(), ".png"),
       g_biol_asv_idx,
       width = 12, height = 6, units = "in")
# ggsave(paste0(projectpath, "/", "usvi_asvs_potential_contaminants_ei-", Sys.Date(), ".png"),
#        g_asv_contam_either,
ggsave(paste0(projectpath, "/", "usvi_species_potential_contaminants_both-", Sys.Date(), ".png"),
       g_species_contam_both,
       width = 12, height = 6, units = "in")
# ggsave(paste0(projectpath, "/", "usvi_asvs_potential_contaminants_prev-", Sys.Date(), ".png"),
#        g_asv_contam_only_prev,
ggsave(paste0(projectpath, "/", "usvi_species_potential_contaminants_prev-", Sys.Date(), ".png"),
       g_species_contam_only_prev,
       width = 16, height = 12, units = "in")
# ggsave(paste0(projectpath, "/", "usvi_asvs_potential_contaminants_nei-", Sys.Date(), ".png"),
#        g_asv_contam_neither,
ggsave(paste0(projectpath, "/", "usvi_species_potential_contaminants_nei-", Sys.Date(), ".png"),
       g_species_contam_neither,
       width = 16, height = 12, units = "in")
# ggsave(paste0(projectpath, "/", "usvi_asvs_potential_contaminants_freq-", Sys.Date(), ".png"),
#        g_asv_contam_only_freq,
ggsave(paste0(projectpath, "/", "usvi_species_potential_contaminants_freq-", Sys.Date(), ".png"),
       g_species_contam_only_freq,
       width = 16, height = 12, units = "in")



# Filter through ASVs that were contaminants but not at species le --------


ps_usvi_samples <- phyloseq::phyloseq(phyloseq::otu_table((usvi_prok_asvs.df %>%
                                                     tidyr::pivot_wider(., id_cols = c("asv_id"),
                                                                        names_from = "sample_ID",
                                                                        values_from = "counts") %>%
                                                     tibble::column_to_rownames(var = "asv_id")),
                                                  taxa_are_rows=TRUE),
                              phyloseq::sample_data(sample_metadata %>%
                                                      tibble::column_to_rownames(var = "sample_ID")),
                              phyloseq::tax_table(usvi_prok_filled.taxa.df %>%
                                                    dplyr::select(-sequence) %>%
                                                    droplevels %>%
                                                    tibble::column_to_rownames(var = "asv_id") %>%
                                                    as.matrix)) %>%
  phyloseq::subset_samples(., grepl("seawater|pcr|extraction", sample_type))
# 
# decon_usvi_freq.df <- decontam::isContaminant(ps_usvi_samples, method = "frequency",
#                                               threshold = set_thresholds[["freq"]],
#                                               batch = sample_data(ps_usvi_samples)$batch,
#                                               conc = "dna_conc")
#                                               # conc = "product_conc")
# 
#using just the concentration of amplicon products in the seawater and control samples:
# table(decon_usvi_freq.df$contaminant)
# #all togehter
# # FALSE  TRUE
# # 12928    52  #threshold = 0.1 default 
# 
# #when doing it by batches:
# #FALSE  TRUE 
# # 12917    63   #threshold = 0.1 default
# 
## using the concentration of extracted DNA (with 0.001 as dummy conc for controls, and random samplingfor the seawater extracted preivously)
# # FALSE  TRUE
# # 12898    82    #threshold = 0.1 default

# head(which(decon_usvi_freq.df$contaminant))
# 33  42  58  78  87 147 #using amplicon
#  85 134 228 264 343 359 #using DNAconcentration

# plot_frequency(ps_usvi_samples, taxa_names(ps_usvi_samples)[c(1,33)], conc="product_conc") +
#   xlab("Amplicon DNA Concentration")
# #the solid red line shows the modeled frequency of a contaminant if its proportional to the DNA conc
# #the dashed black line shows the modeled frequency of a noncontaminant with frequency independent of DNA conc
# 
# plot_frequency(ps_usvi_samples, taxa_names(ps_usvi_samples)[c(1,85,134)], conc="dna_conc") +
#   xlab("extracted DNA Concentration")


# #now look at the qualification of control samples
# decon_usvi_prev.df <- decontam::isContaminant(ps_usvi_samples, method = "prevalence",
#                                               batch = sample_data(ps_usvi_samples)$batch,
#                                               threshold = set_thresholds[["prev"]],
#                                               neg = "is.neg")
#using batches yields the same # of results as not specifying batches
# table(decon_usvi_prev.df$contaminant)
# # FALSE  TRUE 
# # 12969    11  #when thershold = 0.1
# # 12961    19 #when treshold = 0.5

decon_usvi_combo.df <- decontam::isContaminant(ps_usvi_samples,
                                               method = "minimum",
                                               # method = "combined",
                                              batch = sample_data(ps_usvi_samples)$batch,
                                              threshold = set_thresholds[["freq"]],
                                              conc = "dna_conc",
                                              # conc = "product_conc",
                                              neg = "is.neg")
table(decon_usvi_combo.df$contaminant)
# FALSE  TRUE
# 12950    30  #thresold = 0.1, method = "combined" using dna_conc
# 12896    84 # threshold = 0.1, method = "minimum" using dna_conc


{
  ps.pa <- transform_sample_counts(ps_usvi_samples, function(abund) 1*(abund>0))
  ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_type != "seawater", ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type == "seawater", ps.pa)
}

potential_contam_asvs.taxa <- decon_usvi_combo.df %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  dplyr::mutate(contam_freq = dplyr::case_when(p.freq < set_thresholds[["freq"]] ~ "freq", 
                                               .default = NA)) %>%
  dplyr::mutate(contam_prev = dplyr::case_when(p.prev < set_thresholds[["prev"]] ~ "prev",
                                               .default = NA)) %>%
  droplevels %>%
  dplyr::left_join(., data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg)) %>%
                     tibble::rownames_to_column(var = "asv_id"), by = join_by(asv_id)) %>%
  dplyr::mutate(potential_contaminant = dplyr::case_when(
                                                         pa.pos > 50 & pa.neg < 20 ~ 0,
                                                          pa.pos > 20 & pa.neg < 2 ~ 0,
                                                         pa.neg > 0 & contaminant ~ 2,
                                                         pa.neg > 0 ~ 1,
                                                         pa.neg > 0 & (!is.na(contam_freq) | !is.na(contam_prev)) ~ 1,
                                                         # pa.neg == 0 & !contaminant ~ 0,
                                                         # pa.neg == 0 & (is.na(contam_freq) | is.na(contam_prev)) ~ 1,
                                                         pa.neg == 0 & is.na(contam_freq) & is.na(contam_prev) ~ 0,
                                                         .default = 0)) %>%
  dplyr::left_join(., usvi_prok_filled.taxa.df, by = join_by(asv_id)) %>%
  dplyr::relocate(asv_id, pa.pos, pa.neg, potential_contaminant, contaminant, Domain:Species) %>%
  dplyr::filter((potential_contaminant > 0)) %>%
  # dplyr::filter(potential_contaminant > 0 | contaminant | !is.na(contam_freq) | !is.na(contam_prev)) %>%
  droplevels

#drop any ASVs that, after quality filtering for species-level abundances, were not considered contaminants
#also drop any ASVs matching taxonomy to the biological asvs control group

#these are not contaminants:
temp_df <- potential_contam_species.taxa %>%
  dplyr::filter(potential_contaminant == 0) %>%
  dplyr::select(asv_id, c(Domain:Species)) %>%
  dplyr::distinct(., .keep_all = FALSE) %>%
  dplyr::mutate(asv_id = factor(asv_id)) %>%
  dplyr::semi_join((usvi_prok_filled.taxa.df %>%
                      dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
                      droplevels), .,
                   by = join_by(Domain, Phylum, Class, Order, Family, Genus, Species)) %>%
  dplyr::bind_rows(., data.frame(asv_id = biol_asv_idx) %>%
                     dplyr::left_join((usvi_prok_filled.taxa.df %>%
                                         dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
                                         droplevels), by = join_by(asv_id)) %>%
                     dplyr::mutate(asv_id = factor(asv_id)) %>%
                     dplyr::semi_join((usvi_prok_filled.taxa.df %>%
                                         dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
                                         droplevels), .,
                                      by = join_by(Domain, Phylum, Class, Order, Family, Genus, Species))) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::mutate(asv_id = factor(asv_id)) %>%
  dplyr::arrange(asv_id)

#these are potential contaminant ASVs:
potential_contam_asvs.df <- potential_contam_species.taxa %>%
  dplyr::filter(potential_contaminant > 0) %>%
  dplyr::select(asv_id, c(Domain:Species)) %>%
  dplyr::distinct(., .keep_all = FALSE) %>%
  dplyr::mutate(asv_id = factor(asv_id)) %>%
  dplyr::semi_join((usvi_prok_filled.taxa.df %>%
                      dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
                      droplevels), .,
                   by = join_by(Domain, Phylum, Class, Order, Family, Genus, Species)) %>%
  dplyr::bind_rows(., (potential_contam_asvs.taxa %>%
  dplyr::select(asv_id, c(Domain:Species)) %>%
  dplyr::distinct(., .keep_all = FALSE) %>%
  dplyr::mutate(asv_id = factor(asv_id)) %>%
  dplyr::semi_join((usvi_prok_filled.taxa.df %>%
                      dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
                      droplevels), 
                   by = join_by(Domain, Phylum, Class, Order, Family, Genus, Species)))) %>%
  dplyr::anti_join(., temp_df, by = join_by(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species)) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  # dplyr::select(asv_id) %>% simplify %>% as.character
  droplevels
  


# Prune the contaminants --------------------------------------------------
labeller_is.neg <- c(`TRUE` = "negative control", 
                     `FALSE` = "seawater")
labeller_contaminant <- c(`TRUE` = "potential contaminant", 
                          `FALSE` = "not flagged as contaminant")

#these ASVs were flagged through DADA2 as potential contaminants.
# species_contam_both: ASV_00267 and all ASVs of the same taxonomy (through freq and prev) (1)
# species_contam_only_prev: all others in prevalence that were present in negative controls (10 total)
# species_contam_neither: all others that were present in as many negative control samples as in seawater, or flagged through prevalence because not present at all in seawater (10)
# 


# some of these other potential species may just have preferential abundance in certain sites
#27 species were flagged as potential contaminants through frequency but were not observed in negative controls
# potential_contam_freq_summary.df <- plot_potential_contam.df %>%
#   dplyr::filter(sample_type == "seawater") %>%
#   dplyr::filter(asv_id %in% species_contam_only_freq) %>%
#   dplyr::filter(taxon_abundance > 0) %>%
#   droplevels %>%
#   dplyr::group_by(asv_id, site) %>%
#   dplyr::summarise(mean_obs = mean(taxon_abundance, na.rm = TRUE),
#                    max_obs = max(taxon_abundance, na.rm = TRUE)) %>%
#   dplyr::left_join(., plot_potential_contam.df %>%
#                      dplyr::filter(sample_type == "seawater") %>%
#                      dplyr::filter(asv_id %in% species_contam_only_freq) %>%
#                      dplyr::filter(taxon_abundance > 0) %>%
#                      droplevels %>%
#                      dplyr::group_by(asv_id, site) %>%
#                      dplyr::summarise(total_obs = length(taxon_abundance)),
#                    by = join_by(asv_id, site)) %>%
#   dplyr::left_join(., potential_contam_species.taxa %>%
#                      dplyr::select(asv_id, c(Domain:Species)) %>%
#                      dplyr::distinct(., .keep_all = FALSE) %>%
#                      droplevels, by = join_by(asv_id)) %>%
#   dplyr::mutate(site = factor(site, levels = names(site_lookup))) %>%
#   droplevels
# 
# print(
#   ggplot(data = potential_contam_freq_summary.df, 
#          aes(x = asv_id, y = total_obs, fill = site, group = site))
#   + geom_col(position = "dodge")
#   + geom_col(position = "dodge", color = "black", show.legend = FALSE)
#   + theme_bw()
#   + theme(panel.grid.minor.x = element_blank(),
#           axis.text.x = element_text(angle = 90),
#           panel.grid.minor.y = element_blank())
# )
# 
# print(
#   ggplot(data = potential_contam_freq_summary.df, 
#          aes(x = asv_id, y = max_obs*100, fill = site, group = site))
#   + geom_col(position = "dodge")
#   + geom_col(position = "dodge", color = "black", show.legend = FALSE)
#   + theme_bw()
#   + theme(panel.grid.minor.x = element_blank(),
#           axis.text.x = element_text(angle = 90),
#           panel.grid.minor.y = element_blank())
# )

#many of these 27 taxa were present in mainly LB_seagrass samples, with a few occurrences in the other sites
#drop any taxa and their associated ASVs that were present in fewer than 5 samples across all sites
#drop any taxa with a max relative abundance < 0.05%
drop_species_asvs.df <- plot_potential_contam.df %>%
  dplyr::filter(sample_type == "seawater") %>%
  dplyr::filter(asv_id %in% species_contam_only_freq) %>%
  dplyr::filter(taxon_abundance > 0) %>%
  droplevels %>%
  dplyr::group_by(asv_id) %>%
  dplyr::summarise(max_obs = max(taxon_abundance, na.rm = TRUE)) %>%
  dplyr::left_join(., plot_potential_contam.df %>%
                     dplyr::filter(sample_type == "seawater") %>%
                     dplyr::filter(asv_id %in% species_contam_only_freq) %>%
                     dplyr::filter(taxon_abundance > 0) %>%
                     droplevels %>%
                     dplyr::group_by(asv_id) %>%
                     dplyr::summarise(total_obs = length(taxon_abundance)),
                   by = join_by(asv_id)) %>%
  dplyr::arrange(desc(total_obs), desc(max_obs)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(keep = dplyr::case_when((total_obs < 5) ~ 0,
                                        (max_obs < 5e-4 & total_obs < 5) ~ 0,
                                        .default = 1)) %>%
  dplyr::filter(keep < 1) %>%
  dplyr::select(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::bind_rows(., data.frame(asv_id = c(species_contam_both, species_contam_only_prev, species_contam_neither)) %>%
                     dplyr::mutate(asv_id = factor(asv_id))) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::mutate(asv_id = factor(asv_id)) %>%
  dplyr::left_join(., potential_contam_species.taxa %>%
                     dplyr::select(asv_id, c(Domain:Species)) %>%
                     dplyr::distinct(., .keep_all = FALSE) %>%
                     dplyr::mutate(asv_id = factor(asv_id)) %>%
                     droplevels, by = join_by(asv_id)) %>%
  dplyr::semi_join((usvi_prok_filled.taxa.df %>%
                      dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
                      droplevels), .,
                   by = join_by(Domain, Phylum, Class, Order, Family, Genus, Species)) %>%
  dplyr::bind_rows(., potential_contam_asvs.df) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  droplevels




#remove the contaminant asvs
drop_species_asvs_idx <- drop_species_asvs.df %>%
  dplyr::mutate(keep = FALSE) %>%
  dplyr::bind_rows(., data.frame(asv_id = usvi_prok_filled.taxa.df[["asv_id"]])) %>%
  dplyr::mutate(asv_id = factor(asv_id, levels = unique(usvi_prok_filled.taxa.df[["asv_id"]]))) %>%
  dplyr::arrange(asv_id) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::mutate(keep = dplyr::case_when(is.na(keep) ~ TRUE,
                                               .default = keep)) %>%
  droplevels %>%
  dplyr::select(asv_id, keep) %>%
  tibble::deframe(.)

ps_usvi_decontam <- phyloseq::phyloseq(phyloseq::otu_table((usvi_prok_asvs.df %>%
                                                     tidyr::pivot_wider(., id_cols = c("asv_id"),
                                                                        names_from = "sample_ID",
                                                                        values_from = "counts") %>%
                                                     tibble::column_to_rownames(var = "asv_id")),
                                                  taxa_are_rows=TRUE),
                              phyloseq::sample_data(sample_metadata %>%
                                                      tibble::column_to_rownames(var = "sample_ID")),
                              phyloseq::tax_table(usvi_prok_filled.taxa.df %>%
                                                    dplyr::select(-sequence) %>%
                                                    droplevels %>%
                                                    tibble::column_to_rownames(var = "asv_id") %>%
                                                    as.matrix)) %>%
  phyloseq::prune_taxa(drop_species_asvs_idx, .) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove ASVs not present in any samples

readr::write_rds(ps_usvi_decontam, paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"), compress = "gz")

# Summarize sequencing effort ---------------------------------------------

ps_usvi <- ps_usvi_decontam 
usvi_seq_summary.df <- ps_usvi %>%
  phyloseq::otu_table(.) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = -c("asv_id"),
                      names_to = "sample_id",
                      values_to = "abundance") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(total_seqs = sum(abundance, na.rm = TRUE)) %>%
  dplyr::left_join(., (metadata %>%
                         dplyr::select(sample_id, sample_type, sampling_time, sampling_day, site, contains("sample_order"))),
                   by = join_by(sample_id)) %>%
  dplyr::mutate(across(contains("sampl"), ~factor(.x))) %>%
  dplyr::mutate(sample_type = factor(sample_type, levels = temporal_lookup[["sample_type"]]),
                sampling_time = factor(sampling_time, levels = temporal_lookup[["sampling_time"]]),
                # sampling_date = factor(sampling_date, levels = temporal_lookup[["sampling_date"]]),
                sampling_day = factor(sampling_day, levels = temporal_lookup[["sampling_day"]])) %>%
  droplevels

temp_df <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = -c("asv_id"),
                      names_to = "sample_id",
                      values_to = "abundance") %>%
  dplyr::mutate(abundance = dplyr::case_when(abundance > 0 ~ abundance,
                                             .default = NA)) %>%
  tidyr::drop_na(abundance) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(total_ASVs = length(asv_id))

usvi_seq_summary.df %>%
  dplyr::group_by(sample_type, site, sampling_time) %>%
  dplyr::summarise(mean = mean(total_seqs),
                   max = max(total_seqs))

g1 <- print(
  ggplot(data = usvi_seq_summary.df %>%
           dplyr::filter(grepl("seawater", sample_type)) %>%
           droplevels,
         aes(x = sample_order_all, y = total_seqs, fill = sampling_time, group = site))
  + theme_bw()
  + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, alpha = 0.7,
             position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
  + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE, alpha = 0.7,
             position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
  + facet_grid(rows = NULL,
               cols = vars(site),
               scales = "free_x", drop = TRUE, space = "free",
               labeller = labeller(site = site_lookup),
               shrink = FALSE)
  + scale_fill_manual(values = sampling_time_colors, breaks = names(sampling_time_lookup), labels = sampling_time_lookup)
  + theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 90))
  + scale_y_continuous(name = "Sequences", expand = expansion(mult = c(0,0.1)),
                       labels = scales::label_number(digits =3))
  + labs(x = "Sample")
  + guides(fill = guide_legend(order = 1, ncol = 1, title = "Sampling time",  direction = "vertical", 
                               override.aes = list(color = "black", shape = 22, size = 3)),
           color = "none")
)
g2 <- print(
  ggplot(data = usvi_seq_summary.df %>%
           dplyr::filter(!grepl("seawater", sample_type)) %>%
           droplevels,
         aes(x = sample_id, y = total_seqs, fill = sample_type, group = site))
  + theme_bw()
  + geom_bar(stat = "identity", width = 0.90, show.legend = TRUE, alpha = 0.7,
             position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
  + geom_bar(color = "black", stat = "identity", width = 0.90, show.legend = FALSE, alpha = 0.7,
             position = position_dodge2(padding = 0.2, preserve = "total", reverse = TRUE))
  + facet_wrap(.~sample_type, scales = "free", drop = TRUE, shrink = FALSE,
               labeller = labeller(sample_type = site_lookup))
  # + facet_grid(rows = NULL,
  #              cols = vars(sample_type),
  #              scales = "free", drop = TRUE, space = "free",
  #              shrink = FALSE)
  # # + scale_fill_manual(values = sampling_time_colors, breaks = names(sampling_time_lookup), labels = sampling_time_lookup)
  + theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 90))
  + scale_y_continuous(name = "Sequences", expand = expansion(mult = c(0,0.1)),
                       labels = scales::label_number(digits =3))
  + labs(x = "Sample")
  + guides(fill = "none",
           color = "none")
)
gpatch_layout <- "
  AAA
  AAA
  AAA
  BBB
"

g_seq_summary <- g1 + g2 + patchwork::plot_layout(guides = "collect", design = gpatch_layout) + patchwork::plot_annotation(title = "Number of sequences after QC",
                                                                                                   tag_levels = "A")

if(!any(grepl("seqcount", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_summary_seqcount-", Sys.Date(), ".png"),
         g_seq_summary,
         width = 10, height = 10, units = "in")
}

# Broadly summarize at phylum, class, order levels ------------------------

usvi_prok_filled.taxa.df <- usvi_prok_filled.taxa.df %>%
  dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                          grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                          .default = Phylum)) %>%
  droplevels

#make relative abundance tables at specific taxonomic levels
# keep <- c("Phylum", "Genus")
keep <- c("Phylum", "Class", "Order", "Family", "Genus")
# keep <- c("Phylum", "Class", "Order")
# keep <- c("Phylum")

for(taxlevel in keep){
  namevar <- stringr::str_to_lower(taxlevel)
  temp_df <- phyloseq::transform_sample_counts(ps_usvi, function(x) x/sum(x)*100) %>%
    phyloseq::filter_taxa(., function(x) sum(x) > 0, TRUE) %>%
    phyloseq::otu_table(.) %>%
    as.data.frame %>%
    dplyr::mutate(TotAbund = rowSums(.)) %>%
    dplyr::arrange(desc(TotAbund)) %>%
    dplyr::select(-TotAbund) %>%
    tibble::rownames_to_column(var = "asv_id") %>%
    otu_to_taxonomy(dataset = .,
                    level = {{taxlevel}},
                    taxonomy = usvi_prok_filled.taxa.df) %>%
    tidyr::pivot_longer(.,
                        cols = !c(1:all_of({{taxlevel}})),
                        names_to = "sample_id", 
                        values_to = "relabund") %>%
    droplevels
  # temp_df <- phyloseq::otu_table(ps_usvi) %>%
  #   as.data.frame %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame %>%
  #   dplyr::mutate(TotAbund = rowSums(.)) %>%
  #   dplyr::arrange(desc(TotAbund)) %>%
  #   dplyr::select(-TotAbund) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   otu_to_taxonomy(dataset = .,
  #                   level = {{taxlevel}},
  #                   taxonomy = usvi_prok_filled.taxa.df) %>%
  #   tidyr::pivot_longer(.,
  #                       cols = !c(1:all_of({{taxlevel}})),
  #                       names_to = "sample_id",
  #                       values_to = "relabund") %>%
  #   droplevels
  
  # keep_tax <- colnames(usvi_prok_filled.taxa.df)[c(grep("Domain", colnames(usvi_prok_filled.taxa.df)):grep({{taxlevel}}, colnames(usvi_prok_filled.taxa.df)))]
  # 
  # temp_df2 <- temp_df %>%
  #   tidyr::pivot_wider(., id_cols = keep_tax,
  #                      names_from = "sample_id",
  #                      values_fill = 0,
  #                      values_from = "relabund") %>%
  #   tidyr::unite("taxonomy", c(all_of(keep_tax)), sep = ";", remove = TRUE) %>%
  #   dplyr::relocate(taxonomy) %>%
  #   tibble::column_to_rownames(var = "taxonomy") %>%
  #   dplyr::slice(which(rowSums(.) > 5)) %>%
  #   t()
  #   # droplevels
  # 
  # temp_df2 <- temp_df2 %>%
  #   bind_cols("sample_id" = rownames(temp_df2),
  #             ., 
  #             "Other taxa" = (100 - rowSums(temp_df2))) %>%
  #   tidyr::pivot_longer(., cols = !c("sample_id"),
  #                       names_to = "taxonomy",
  #                       values_to = "relabund") %>%
  #   dplyr::arrange(desc(relabund)) %>%
  #   dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  #   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  #   tidyr::separate_wider_delim(taxonomy, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
  #   droplevels
  
  assign(paste0("usvi_", namevar, ".df"), temp_df, envir = .GlobalEnv)
  # assign(paste0("usvi_", namevar, "_over5.df"), temp_df2, envir = .GlobalEnv)
  rm(temp_df)
}

#now work backwards to keep the class, order, family related to the top5 genera
keep_tax <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# usvi_genus_over5.df <- usvi_genus.df %>%
#   dplyr::semi_join(., usvi_seq_summary.df %>% #use only biological samples to filter taxonomic levels
#                   dplyr::filter(sample_type == "seawater") %>%
#                     dplyr::select(sample_id) %>%
#                     droplevels, by = join_by(sample_id)) %>%
#   tidyr::pivot_wider(., id_cols = all_of(keep_tax),
#                      names_from = "sample_id",
#                      values_fill = 0,
#                      values_from = "relabund") %>%
#   tidyr::unite("taxonomy", c(all_of(keep_tax)), sep = ";", remove = TRUE) %>%
#   dplyr::relocate(taxonomy) %>%
#   tibble::column_to_rownames(var = "taxonomy") %>%
#   dplyr::slice(which(rowSums(.) >= 5)) %>%
#   t()
# 
# usvi_genus_over5.df <- usvi_genus_over5.df %>%
#   # tibble::as_tibble(rownames = "sample_id") %>% #in case you don't want the "Other taxa" to confuse
#   bind_cols("sample_id" = rownames(usvi_genus_over5.df),
#             .,
#             "Other taxa" = (100 - rowSums(usvi_genus_over5.df))) %>%
#   tidyr::pivot_longer(., cols = !c("sample_id"),
#                       names_to = "taxonomy",
#                       values_to = "relabund") %>%
#   dplyr::arrange(desc(relabund)) %>%
#   dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
#   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
#   tidyr::separate_wider_delim(taxonomy, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
#   droplevels
# length(unique(usvi_genus_over5.df[["Genus"]]))
# length(unique(usvi_genus_over5.df[["Family"]]))
# 
# #the 72 unique genera in the "over 5" filter, represents 55 distinct families
# usvi_genus_over5.df %>% dplyr::filter(!grepl("Other taxa", Domain)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(sample_relabund = sum(relabund, na.rm = TRUE)) %>% dplyr::filter(grepl("Metab", sample_id)) %>% droplevels %>% dplyr::reframe(across(contains("relabund"), list(min = min, max = max)))
# #minimum of 91.3%, maximum of 99.2% of sequences in these unique genera


#what if we retained only the top 10% of genera?
usvi_genus_over10.df <- usvi_genus.df %>%
  dplyr::semi_join(., usvi_seq_summary.df %>% #use only biological samples to filter taxonomic levels
                     dplyr::filter(sample_type == "seawater") %>%
                     dplyr::select(sample_id) %>%
                     droplevels, by = join_by(sample_id)) %>%
  tidyr::pivot_wider(., id_cols = all_of(keep_tax),
                     names_from = "sample_id",
                     values_fill = 0,
                     values_from = "relabund") %>%
  tidyr::unite("taxonomy", c(all_of(keep_tax)), sep = ";", remove = TRUE) %>%
  dplyr::relocate(taxonomy) %>%
  tibble::column_to_rownames(var = "taxonomy") %>%
  dplyr::slice(which(rowSums(.) >= 10)) %>%
  t()

usvi_genus_over10.df <- usvi_genus_over10.df %>%
  # tibble::as_tibble(rownames = "sample_id") %>% #in case you don't want the "Other taxa" to confuse
  bind_cols("sample_id" = rownames(usvi_genus_over10.df),
            .,
            "Other taxa" = (100 - rowSums(usvi_genus_over10.df))) %>%
  tidyr::pivot_longer(., cols = !c("sample_id"),
                      names_to = "taxonomy",
                      values_to = "relabund") %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  tidyr::separate_wider_delim(taxonomy, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
  droplevels
length(unique(usvi_genus_over10.df[["Genus"]]))
length(unique(usvi_genus_over10.df[["Family"]]))
#43 families constitute genera with 10% or more realtive abundance
usvi_genus_over10.df %>% dplyr::filter(!grepl("Other taxa", Domain)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(sample_relabund = sum(relabund, na.rm = TRUE)) %>% dplyr::filter(grepl("Metab", sample_id)) %>% droplevels %>% dplyr::reframe(across(contains("relabund"), list(min = min, max = max)))

#use the top genera 10% or more of all sequences to further filter.
usvi_filtered_genus.df <- usvi_genus_over10.df

keep_tax2 <- rev(rev(keep_tax)[-1])
usvi_family_filtered.df <- usvi_filtered_genus.df %>%
  dplyr::filter(!grepl("Other taxa", Domain)) %>%
  dplyr::select(all_of(keep_tax2)) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::inner_join(., usvi_family.df,
                   by = join_by(!!!keep_tax2),
                   relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = keep_tax2,
                     names_from = "sample_id",
                     values_fill = 0,
                     values_from = "relabund") %>%
  tidyr::unite("taxonomy", c(all_of(keep_tax2)), sep = ";", remove = TRUE) %>%
  dplyr::relocate(taxonomy) %>%
  tibble::column_to_rownames(var = "taxonomy") %>%
  dplyr::slice(which(rowSums(.) >= 5)) %>%
  # droplevels
  t()
usvi_family_filtered.df <- usvi_family_filtered.df %>%
  # tibble::as_tibble(rownames = "sample_id") %>% #in case you don't want the "Other taxa" to confuse
  bind_cols("sample_id" = rownames(usvi_family_filtered.df),
            .,
            "Other taxa" = (100 - rowSums(usvi_family_filtered.df))) %>%
  tidyr::pivot_longer(., cols = !c("sample_id"),
                      names_to = "taxonomy",
                      values_to = "relabund") %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  tidyr::separate_wider_delim(taxonomy, names = keep_tax2, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
  droplevels
length(unique(usvi_family_filtered.df[["Family"]]))
length(unique(usvi_family_filtered.df[["Class"]]))

usvi_family_filtered.df %>% dplyr::filter(!grepl("Other taxa", Domain)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(sample_relabund = sum(relabund, na.rm = TRUE)) %>% dplyr::filter(grepl("Metab", sample_id)) %>% droplevels %>% dplyr::reframe(across(contains("relabund"), list(min = min, max = max)))

#filtere for 10% or more genera:
#75.4% to 99.4% of sequences in genera from 43 families across 12 classes

keep_tax3 <- rev(rev(keep_tax2)[-1])
usvi_order_filtered.df <- usvi_family_filtered.df %>%
  dplyr::filter(!grepl("Other taxa", Domain)) %>%
  dplyr::select(all_of(keep_tax3)) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::inner_join(., usvi_order.df,
                    by = join_by(!!!keep_tax3),
                    relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = all_of(keep_tax3),
                     names_from = "sample_id",
                     values_fill = 0,
                     values_from = "relabund") %>%
  tidyr::unite("taxonomy", c(all_of(keep_tax3)), sep = ";", remove = TRUE) %>%
  dplyr::relocate(taxonomy) %>%
  tibble::column_to_rownames(var = "taxonomy") %>%
  dplyr::slice(which(rowSums(.) >= 5)) %>%
  # droplevels
  t()
usvi_order_filtered.df <- usvi_order_filtered.df %>%
  # tibble::as_tibble(rownames = "sample_id") %>% #in case you don't want the "Other taxa" to confuse
  bind_cols("sample_id" = rownames(usvi_order_filtered.df),
            .,
            "Other taxa" = (100 - rowSums(usvi_order_filtered.df))) %>%
  tidyr::pivot_longer(., cols = !c("sample_id"),
                      names_to = "taxonomy",
                      values_to = "relabund") %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  tidyr::separate_wider_delim(taxonomy, names = keep_tax3, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
  droplevels
length(unique(usvi_order_filtered.df[["Class"]]))
length(unique(usvi_order_filtered.df[["Order"]]))

usvi_order_filtered.df %>% dplyr::filter(!grepl("Other taxa", Domain)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(sample_relabund = sum(relabund, na.rm = TRUE)) %>% dplyr::filter(grepl("Metab", sample_id)) %>% droplevels %>% dplyr::reframe(across(contains("relabund"), list(min = min, max = max)))

#10% genera filter
#28 orders in 12 classes
#90.8% to 99.4% sequences in environmental samples belong to these orders


keep_tax4 <- rev(rev(keep_tax3)[-1])
usvi_class_filtered.df <- usvi_order_filtered.df %>%
  dplyr::filter(!grepl("Other taxa", Domain)) %>%
  dplyr::select(all_of(keep_tax4)) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::inner_join(., usvi_class.df,
                    by = join_by(!!!keep_tax4),
                    relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = all_of(keep_tax4),
                     names_from = "sample_id",
                     values_fill = 0,
                     values_from = "relabund") %>%
  tidyr::unite("taxonomy", c(all_of(keep_tax4)), sep = ";", remove = TRUE) %>%
  dplyr::relocate(taxonomy) %>%
  tibble::column_to_rownames(var = "taxonomy") %>%
  dplyr::slice(which(rowSums(.) >= 5)) %>%
  # droplevels
  t()
usvi_class_filtered.df <- usvi_class_filtered.df %>%
  # tibble::as_tibble(rownames = "sample_id") %>% #in case you don't want the "Other taxa" to confuse
  bind_cols("sample_id" = rownames(usvi_class_filtered.df),
            .,
            "Other taxa" = (100 - rowSums(usvi_class_filtered.df))) %>%
  tidyr::pivot_longer(., cols = !c("sample_id"),
                      names_to = "taxonomy",
                      values_to = "relabund") %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  tidyr::separate_wider_delim(taxonomy, names = keep_tax4, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
  droplevels
length(unique(usvi_class_filtered.df[["Class"]]))
length(unique(usvi_class_filtered.df[["Phylum"]]))

usvi_class_filtered.df %>% dplyr::filter(!grepl("Other taxa", Domain)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(sample_relabund = sum(relabund, na.rm = TRUE)) %>% dplyr::filter(grepl("Metab", sample_id)) %>% droplevels %>% dplyr::reframe(across(contains("relabund"), list(min = min, max = max)))

#10% genera filter:
#12 classes in 10 phyla 
#90.8% to 100.7% of sequences are in these classes


keep_tax5 <- rev(rev(keep_tax4)[-1])
usvi_phylum_filtered.df <- usvi_class_filtered.df %>%
  dplyr::filter(!grepl("Other taxa", Domain)) %>%
  dplyr::select(all_of(keep_tax5)) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  droplevels %>%
  dplyr::inner_join(., usvi_phylum.df,
                    by = join_by(!!!keep_tax5),
                    relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = all_of(keep_tax5),
                     names_from = "sample_id",
                     values_fill = 0,
                     values_from = "relabund") %>%
  tidyr::unite("taxonomy", c(all_of(keep_tax5)), sep = ";", remove = TRUE) %>%
  dplyr::relocate(taxonomy) %>%
  tibble::column_to_rownames(var = "taxonomy") %>%
  dplyr::slice(which(rowSums(.) >= 5)) %>%
  # droplevels
  t()
usvi_phylum_filtered.df <- usvi_phylum_filtered.df %>%
  # tibble::as_tibble(rownames = "sample_id") %>% #in case you don't want the "Other taxa" to confuse
  bind_cols("sample_id" = rownames(usvi_phylum_filtered.df),
            .,
            "Other taxa" = (100 - rowSums(usvi_phylum_filtered.df))) %>%
  tidyr::pivot_longer(., cols = !c("sample_id"),
                      names_to = "taxonomy",
                      values_to = "relabund") %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  tidyr::separate_wider_delim(taxonomy, names = keep_tax5, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
  droplevels
length(unique(usvi_phylum_filtered.df[["Domain"]]))
length(unique(usvi_class_filtered.df[["Phylum"]]))

usvi_phylum_filtered.df %>% dplyr::filter(!grepl("Other taxa", Domain)) %>% dplyr::group_by(sample_id) %>% dplyr::summarise(sample_relabund = sum(relabund, na.rm = TRUE)) %>% dplyr::filter(grepl("Metab", sample_id)) %>% droplevels %>% dplyr::reframe(across(contains("relabund"), list(min = min, max = max)))
#10 phyla in 2 domains
#93.9% to 100% of sequences in environmental samples are represented in these 13 phyla




# Make color codes --------------------------------------------------------
if(!exists("annotation_taxa_colors_list", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "annotation_taxa_colors_list.rds"))){
    annotation_taxa_colors_list <- readr::read_rds(paste0(projectpath, "/", "annotation_taxa_colors_list.rds"))
  } else {
    #make colors list for plotting taxonomy of bins
    annotation_archaea_list <- usvi_filtered_genus.df %>%
      dplyr::select(Domain, Phylum, Class, Order, Family, Genus, taxonomy) %>%
      dplyr::filter(grepl("Archaea", Domain)) %>%
      dplyr::distinct(., .keep_all = TRUE) %>%
      droplevels %>%
      tidyr::pivot_longer(., cols = !c(taxonomy),
                          names_to = "taxonomic_level",
                          values_to = "value") %>%
      dplyr::mutate(taxonomic_level = factor(taxonomic_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))) %>%
      split(., f = .$taxonomic_level) %>%
      map(., ~.x %>%
            dplyr::distinct(value, .keep_all = TRUE) %>%
            tidyr::drop_na(value) %>%
            rowwise(.) %>%
            dplyr::mutate(new_t = stringr::str_split_i(taxonomy, value, 1) %>%
                            paste0(., value)) %>%
            dplyr::select(new_t, value) %>%
            droplevels %>%
            tibble::deframe(.))
    
    annotation_archaea_colors <- annotation_archaea_list  %>%
      purrr::flatten(.) 
    annotation_archaea_colors <- pals::ocean.tempo(n = length(annotation_archaea_colors)) %>%
      setNames(., names(annotation_archaea_colors)) %>%
      tibble::enframe(., name = "taxonomy_string", value = "value") %>%
      dplyr::mutate(taxonomy_string = factor(taxonomy_string))
    #what does it look like?
    {
      temp_list1 <- annotation_archaea_list %>%
        map(., ~.x %>%
              tibble::enframe(., name = "taxonomy_string", value = "classification") %>%
              tidyr::drop_na()) %>%
        map(., ~.x %>%
              dplyr::left_join(., annotation_archaea_colors, by = join_by(taxonomy_string))) %>%
        bind_rows(., .id = "taxonomic_level")
      }
    
    #now bacteria:
    annotation_bacteria_list <- usvi_filtered_genus.df %>%
      # dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
      #                                         grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
      #                                         .default = Phylum)) %>%
      dplyr::select(Domain, Phylum, Class, Order, Family, Genus, taxonomy) %>%
      dplyr::filter(!grepl("Archaea", Domain)) %>%
      dplyr::distinct(., .keep_all = TRUE) %>%
      # split(., f = .$Class) %>%
      split(., f = .$Phylum) %>%
      map(., ~.x %>%
            tidyr::pivot_longer(., cols = !c(taxonomy),
                                names_to = "taxonomic_level",
                                values_to = "value") %>%
            dplyr::mutate(taxonomic_level = factor(taxonomic_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))) %>%
            dplyr::filter(!grepl("Domain", taxonomic_level)) %>%
            split(., f = .$taxonomic_level) %>%
            map(., ~.x %>%
                  dplyr::distinct(value, .keep_all = TRUE) %>%
                  tidyr::drop_na(value) %>%
                  rowwise(.) %>%
                  dplyr::mutate(new_t = stringr::str_split_i(taxonomy, value, 1) %>%
                                  paste0(., value)) %>%
                  dplyr::select(new_t, value) %>%
                  droplevels %>%
                  tibble::deframe(.)))
    
    annotation_bacteria_colors_list <- NULL
    for(i in names(annotation_bacteria_list)){
      # { i <- names(annotation_bacteria_list)[11]
      annotation_temp_colors <- annotation_bacteria_list[[i]] %>%
        Filter(Negate(function(x) is.null(unlist(x))), .) %>%
        purrr::flatten(.)
      
      annotation_bacteria_colors_list <- c(annotation_bacteria_colors_list, list(annotation_temp_colors)) %>%
        setNames(., c(names(annotation_bacteria_colors_list), i))
      
      rm(annotation_temp_colors)
    }
    
    bacteria_options <- c(brewer.blues, brewer.gnbu,  brewer.oranges, brewer.greens, brewer.purples, 
                         
                          brewer.purd, 
                          brewer.bugn, brewer.bupu, brewer.rdpu, brewer.reds, brewer.ylgn,
                          brewer.orrd, brewer.pubu, brewer.pubugn,
                          brewer.ylgnbu, brewer.ylorbr, brewer.ylorrd)
    
    annotation_bact_phyla_solo <- NULL
    annotation_bact_phyla <- NULL
    for(i in c(1:length(annotation_bacteria_colors_list))){
      namevar <- names(annotation_bacteria_colors_list)[i]
      temp_list <- annotation_bacteria_colors_list[namevar] %>%
        purrr::flatten(.)
      if(length(temp_list) > 2){ #is there is more than one class in this phylum in the MAGs
        # if(length(temp_list) > 4){ #is there is more specific taxonomic levels for this phylum
        annotation_bact_phyla <- c(annotation_bact_phyla, namevar)
      } else {
        annotation_bact_phyla_solo <- c(annotation_bact_phyla_solo, namevar)
      }
    }
    
    annotation_bacteria_colors <- NULL
    for(i in c(1:length(annotation_bact_phyla))){
      namevar <- annotation_bact_phyla[i]
      temp_list <- annotation_bacteria_colors_list[namevar] %>%
        purrr::flatten(.)
      
      temp_palette <- bacteria_options[[i]]
      
      if(length(temp_list) > 2){
        annotation_temp_colors <- temp_palette(n = length(temp_list)+2)[1:length(temp_list)] %>%
          setNames(., rev(names(temp_list)))
        
        annotation_bacteria_colors <- c(annotation_bacteria_colors, annotation_temp_colors) %>%
          setNames(., c(names(annotation_bacteria_colors), names(temp_list)))  
      } 
    }
    
    annotation_bacteria_solo_colors <- NULL
    if(!is.null(annotation_bact_phyla_solo)){
      temp_list <- annotation_bact_phyla_solo %>%
        map(., ~annotation_bacteria_colors_list %>%
              purrr::pluck(.x)) %>%
        purrr::flatten(.)
      temp_palette <- pals::brewer.greys(n = length(temp_list)) %>%
        # temp_palette <- kovesi.isoluminant_cgo_80_c38(n = length(temp_list)+2)[1:length(temp_list)] %>%
        setNames(., rev(names(temp_list)))
      annotation_bacteria_solo_colors <- temp_palette
    }
    rm(temp_palette)
    rm(annotation_bacteria_colors_list)
    rm(bacteria_options)
    
    
    temp_df <- annotation_bacteria_colors %>%
      tibble::enframe(., name = "taxonomy_string", value = "value") %>%
      dplyr::mutate(taxonomy_string = factor(taxonomy_string))
    if(!is.null(annotation_bacteria_solo_colors)){
      temp_df2 <- annotation_bacteria_solo_colors %>%
        tibble::enframe(., name = "taxonomy_string", value = "value") %>%
        dplyr::mutate(taxonomy_string = factor(taxonomy_string))  
      temp_df <- bind_rows(temp_df, temp_df2)
    }
    
    temp_list2 <- annotation_bacteria_list %>%
      map(., ~.x %>%
            map(., ~.x %>%
                  tibble::enframe(., name = "taxonomy_string", value = "classification") %>%
                  tidyr::drop_na()) %>%
            map(., ~.x %>%
                  dplyr::left_join(., temp_df, by = join_by(taxonomy_string)) %>%
                  # dplyr::left_join(., temp_df2, by = join_by(taxonomy_string)) %>%
                  dplyr::mutate(value = across(starts_with("value")) %>% purrr::reduce(coalesce)) %>%
                  dplyr::select(-ends_with(c(".x", ".y"))))) %>%
      map(., ~.x %>%
            bind_rows(., .id = "taxonomic_level")) %>%
      bind_rows(.) %>%
      bind_rows(., data.frame(taxonomic_level = "Domain", taxonomy_string = "Bacteria", classification = "Bacteria", value = "#000000"))
    
    temp_list <- bind_rows(temp_list1, temp_list2) %>%
      split(., f = .$taxonomic_level) %>%
      map(., ~.x %>%
            dplyr::distinct(classification, .keep_all = TRUE) %>%
            tidyr::drop_na(classification) %>%
            rowwise(.) %>%
            dplyr::mutate(new_t = stringr::str_split_i(taxonomy_string, ";", -2)) %>%
            dplyr::mutate(new_t = classification) %>%
            # dplyr::mutate(new_t = dplyr::case_when(!is.na(new_t) ~ paste0(new_t, ";", classification),
            #                                        .default = classification)) %>%
            # dplyr::mutate(new_t = stringr::str_split_i(taxonomy_string, classification, 1) %>%
            # paste0(., ";", classification)) %>%
            dplyr::select(new_t, value) %>%
            droplevels %>%
            tibble::deframe(.)) %>%
      map(., ~c(.x, "#FFFFFF") %>%
            setNames(., c(names(.x), "NA")))
    temp_list[["Domain"]] <- temp_list[["Domain"]][1:2]
    annotation_taxa_colors_list <- temp_list
    
    readr::write_rds(annotation_taxa_colors_list, 
                     paste0(projectpath, "/", "annotation_taxa_colors_list.rds"))
    
    rm(list = (apropos("^temp_list*", mode = "list")))
    rm(list = (apropos("^annotation_.*_colors$", mode = "any")))
    rm(list = (apropos("^annotation_bact_phyla*", mode = "any")))
  }
}

# broad Taxonomy summary --------------------------------------------------



ASV_env_phylum_summary.tbl <- usvi_phylum_filtered.df %>%
# ASV_env_phylum_summary.tbl <- usvi_phylum_over5.df %>%
  dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
  dplyr::group_by(taxonomy, sample_type, site, sampling_time) %>%
  dplyr::summarise(across(contains("relabund"), list(mean = mean, min = min, max = max))) %>%
  tidyr::pivot_longer(., cols = contains("relabund"),
                      names_to = "metric",
                      values_to = "relabund") %>%
  dplyr::mutate(metric = gsub("relabund_", "", metric)) %>%
  droplevels

phylum_colors <- dplyr::left_join(annotation_taxa_colors_list[["Genus"]] %>%
                                    tibble::enframe(., name = "Genus", value = "color"), 
                                  usvi_filtered_genus.df %>%
                                    # dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                    #                                         grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                    #                                         .default = Phylum)) %>%
                                    dplyr::select(Genus, Domain, Phylum) %>%
                                    dplyr::distinct(Genus, Domain, Phylum, .keep_all = FALSE) %>%
                                    dplyr::mutate(Phylum = dplyr::case_when(is.na(Phylum) ~ "NA",
                                                                            .default = Phylum)),
                                  by = join_by(Genus)) %>%
  tidyr::unite("taxonomy", c(Domain, Phylum), sep = ";", remove = TRUE, na.rm = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  dplyr::select(taxonomy, color) %>%
  tibble::deframe(.)



df <- ASV_env_phylum_summary.tbl %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  droplevels
df2 <- ASV_env_phylum_summary.tbl %>%
  dplyr::filter(!grepl("seawater", sample_type)) %>%
  droplevels

{
  g1 <- P_summary_taxonomy_bar(dataset = df %>%
                                 dplyr::filter(grepl("LB_seagrass", site)),
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "sampling_time ~ metric")
  
  g2 <- P_summary_taxonomy_bar(dataset = df %>%
                                 dplyr::filter(grepl("Yawzi", site)),
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "sampling_time~metric")
  
  g3 <- P_summary_taxonomy_bar(dataset = df %>%
                                 dplyr::filter(grepl("Tektite", site)),
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "sampling_time~metric")
  
  g4 <- P_summary_taxonomy_bar(dataset = df2,
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "sample_type~metric")
  
  
  gpatch <- (g1 + ggtitle("Lameshur Bay seagrass") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")) / (g2 + ggtitle("Yawzi reef") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text.x = element_blank(), legend.position = "none")) / (g3 + ggtitle("Tektite reef") + theme(axis.text.x = element_blank(), strip.text.x = element_blank(), legend.direction = "horizontal")) / (g4 + ggtitle("Controls") + theme(axis.text.x = element_blank(), strip.text.x = element_blank(), legend.direction = "horizontal")) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right")
  gpatch <- gpatch & scale_fill_manual(values = phylum_colors, drop = TRUE, guide = "legend")
  gpatch <- gpatch + patchwork::plot_annotation(title = "Summary of relative abundances of phyla by site", 
                                                tag_level = "A")
  ggsave(paste0(projectpath, "/", "usvi_phylum_abundance_stats-", Sys.Date(), ".png"),
         gpatch,
         width = 10, height = 16, units = "in")
  gpatch_log10 <- gpatch & scale_y_continuous(trans = "pseudo_log", expand = expansion(mult = c(0,0.1)))
  ggsave(paste0(projectpath, "/", "usvi_phylum_abundance_stats_log10-", Sys.Date(), ".png"),
         gpatch_log10,
         width = 10, height = 16, units = "in")
}

# Zoom in on specific classes and orders ----------------------------------





#based on the relative abundances of phyla, focus on these:
#Proteobacteria
#Bacteroidota
#Cyanobacteria
#SAR406
#Verrucomicrobia
#Thermoplasmota
#Actinobacteria

keep <- c("Proteobacteria", "Bacteroidota", "Cyanobacteria", "Actinobacteriota", 
          "Alphaproteobacteria", "Gammaproteobacteria",
          "SAR406", "Verrucomicrobiota", "Thermoplasmatota")

taxonomy_colors <- dplyr::left_join(annotation_taxa_colors_list[["Genus"]] %>%
                                      tibble::enframe(., name = "Genus", value = "color"), 
                                    usvi_filtered_genus.df %>%
                                      dplyr::select(Genus, Domain, Phylum, Class) %>%
                                      dplyr::distinct(Genus, Domain, Phylum, Class, .keep_all = FALSE) %>%
                                      dplyr::mutate(Class = dplyr::case_when(is.na(Class) ~ "NA",
                                                                             .default = Class)),
                                    by = join_by(Genus)) %>%
  tidyr::unite("taxonomy", c(Domain, Phylum, Class), sep = ";", remove = TRUE, na.rm = TRUE) %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  dplyr::select(taxonomy, color) %>%
  tibble::deframe(.)

temp_df <- usvi_class_filtered.df %>%
  dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_class_filtered.df[["Phylum"]]), value = TRUE)) %>%
  droplevels %>%
  dplyr::arrange(desc(relabund), Class, taxonomy) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_colors))) %>%
  droplevels %>%
  dplyr::arrange(sampling_day, sampling_time, replicate) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
  dplyr::mutate(sample_order = factor(sample_order, levels = unique(.[["sample_order"]]))) %>%
  # dplyr::group_by(taxonomy, sample_type, site, sampling_day, sampling_time,replicate, sample_id, sample_order) %>%
  dplyr::distinct(sample_id, taxonomy, .keep_all = TRUE) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
  droplevels

temp_df2 <- usvi_order_filtered.df %>%
  dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_order_filtered.df[["Phylum"]]), value = TRUE)) %>%
  dplyr::arrange(desc(relabund), Class, taxonomy) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_colors))) %>%
  # dplyr::group_by(taxonomy, sample_type, site, sample_order, sampling_day, sampling_time, replicate, sample_id) %>%
  dplyr::arrange(sampling_day, sampling_time, replicate) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
  dplyr::mutate(sample_order = factor(sample_order, levels = unique(.[["sample_order"]]))) %>%
  dplyr::distinct(sample_id, taxonomy, .keep_all = TRUE) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
  droplevels

temp_df %>%
  dplyr::group_by(site, taxonomy) %>%
  # dplyr::summarise(mean_abund = sum(relabund, na.rm = TRUE)) %>%
  dplyr::summarise(mean_abund = mean(relabund, na.rm = TRUE)) %>%
  dplyr::slice_max(mean_abund, n = 10) %>%
  dplyr::select(taxonomy) %>%
  dplyr::ungroup(.) %>%
  dplyr::distinct(taxonomy, .keep_all = FALSE) %>%
  simplify %>%
  unlist %>%
  as.character


{
  print(ggplot(data = temp_df %>%
                # dplyr::filter(grepl("LB", site)) %>%
                droplevels, aes(y = relabund,
                                  x = sample_order_all,
                                group = taxonomy,
                                  fill = taxonomy))
  + geom_bar(width = 0.90, show.legend = TRUE,
                    # position = "dodge",
                    position = "stack",
                    stat = "identity")
  + geom_bar(color = "black", width = 0.90, show.legend = FALSE,
             # position = "dodge",
             position = "stack",
             stat = "identity")
  # + geom_col(width = 0.90, show.legend = TRUE, position = position_stack(reverse = FALSE))
  # + geom_col(color = "black", width = 0.90, show.legend = FALSE, position = position_stack(reverse = FALSE))
  + theme_bw() 
  + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  + scale_fill_manual(values = taxonomy_colors)
  + facet_grid(scales = "free", space = "free",
                      labeller = labeller(site = site_lookup),
                      # rows = vars({{facet_rows}}), cols = vars({{facet_cols}}))
                      . ~ site)
  + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                 strip.text.x = element_text(angle = 0),
                 strip.text.y = element_text(angle = -90))
  + labs(x = "sample",
                y = "Relative abundance (%)")
  + coord_cartesian(ylim = c(0,100), expand = FALSE)
  + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
                                      override.aes = list(color = "black", stroke = 1)))
  )
  }


