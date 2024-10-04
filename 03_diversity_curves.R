# 03_diversity_curves.R

#make rarefaction curves, etc. for ASVs in the USVI temporal series


# Resource allocation time ------------------------------------------------


if(file.exists(paste0(getwd(), "/", "00_resource_allocation.R"))){
  cat("Preparing resource allocations.")
  source(paste0(getwd(), "/", "00_resource_allocation.R"), local = FALSE)
} else {
  cat("Preparing resource allocations.")
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(data.table)
  library(BiocManager)
  library(BiocParallel)
  library(cli)
  library(furrr)
  library(progressr)
  
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
  
  
  # nthreads <- data.table::getDTthreads()
  # cluster <- multidplyr::new_cluster(n = nthreads)
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
  
}

# Load packages -----------------------------------------------------------

library(tidyverse)
library(phyloseq)

# if(!require("ggvegan", quietly = TRUE)){
#   install.packages("remotes")
#   remotes::install_github("gavinsimpson/ggvegan")
# }
library(ggvegan)
library(ggtext)
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)



# Custom functions --------------------------------------------------------

#Turn data into relative abundances (since these data are compositional = quantitative descriptions of the parts of some whole, conveying relative information.)
relabund <- function(sample) { #Write a relative abundance function
  x = sample/sum(sample)
  x = x*100
  return(x)
}

scaleFUN2 <- function(x) sprintf("%.2f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)
scaleFUN0 <- function(x) sprintf("%.0f", x)


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
  if(file.exists(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))){
    ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))
  } else {
    ps_usvi <- phyloseq::phyloseq(phyloseq::otu_table((usvi_prok_asvs.df %>%
                                                         tidyr::pivot_wider(., id_cols = c("asv_id"),
                                                                            names_from = "sample_ID",
                                                                            values_from = "counts") %>%
                                                         tibble::column_to_rownames(var = "asv_id")),
                                                      taxa_are_rows=TRUE),
                                  phyloseq::sample_data(metadata %>%
                                                          tibble::column_to_rownames(var = "sample_ID")),
                                  phyloseq::tax_table(usvi_prok_asvs.taxa %>%
                                                        dplyr::select(-sequence) %>%
                                                        droplevels %>%
                                                        tibble::column_to_rownames(var = "asv_id") %>%
                                                        as.matrix))
    
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


# Prepare alpha diversity estimates ---------------------------------------

usvi_seqdepth <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
# usvi_seqdepth <- phyloseq::otu_table(ps_usvi) %>%
  as.data.frame %>%
  t() %>%
  rowSums(.)
usvi_seqdepth_stats <- plyr::each(min, max)(usvi_seqdepth)

usvi_rarefied_stats.df <- phyloseq::otu_table(ps_usvi) %>%
  as.data.frame %>%
  t() %>%
  vegan::rarecurve(., usvi_seqdepth_stats[["min"]], tidy = TRUE) %>%
  setNames(., c("sample_id", "seqdepth", "num_ASVs")) %>%
  left_join(., metadata %>%
              dplyr::select(sample_id, sample_type, sampling_date, sampling_time, sampling_day, site, replicate) %>%
              droplevels) %>%
  dplyr::mutate(across(c(sample_id, sample_type, sampling_time, sampling_day, site, replicate), factor))

usvi_rarefied_list <- list((usvi_rarefied_stats.df %>%
                      dplyr::filter(grepl("seagrass", site)) %>%
                      droplevels),
                   (usvi_rarefied_stats.df %>%
                      dplyr::filter(grepl("Tektite", site)) %>%
                      droplevels),
                   (usvi_rarefied_stats.df %>%
                      dplyr::filter(grepl("Yawzi", site)) %>%
                      droplevels),
                   (usvi_rarefied_stats.df %>%
                      dplyr::filter(grepl("control_extraction", sample_type)) %>%
                      droplevels),
                   (usvi_rarefied_stats.df %>%
                      dplyr::filter(grepl("control_pcr", sample_type)) %>%
                      droplevels),
                   (usvi_rarefied_stats.df %>%
                      dplyr::filter(grepl("control_seq", sample_type)) %>%
                      droplevels)) %>%
  setNames(., c("LB_seagrass", "Tektite", "Yawzi", "control_extraction", "control_pcr", "control_seq"))


# Plot rarefaction curves -------------------------------------------------

# temp_list <- usvi_rarefied_list
temp_list1 <- c("LB_seagrass", "Tektite", "Yawzi") %>%
  purrr::map(., ~purrr::pluck(usvi_rarefied_list, .x)) %>%
  setNames(., c("LB_seagrass", "Tektite", "Yawzi"))
temp_list2 <- c("control_extraction", "control_pcr", "control_seq") %>%
  purrr::map(., ~purrr::pluck(usvi_rarefied_list, .x)) %>%
  setNames(., c("control_extraction", "control_pcr", "control_seq"))

xlim <- max(usvi_rarefied_stats.df %>%
              dplyr::select(seqdepth) %>%
              deframe(.) %>%
              max(.))
xlim <- round(xlim/100, digits = -1)
xlim <- xlim*100

#facetting by sampling time, for each site:
{
  # for(i in 1:length(temp_list1)){
  #   gtitle <- eval(names(temp_list1)[i])
  #   
  #   ylim <- temp_list1[[i]] %>%
  #     dplyr::select(num_ASVs) %>%
  #     deframe(.) %>%
  #     max(.)
  #   ylim <- ceiling(ylim/100)*100
  #   
  #   g <- print(
  #     ggplot(data = temp_list1[[i]],
  #            aes(x = seqdepth, y = num_ASVs, group = sample_id))
  #     + theme_bw()
  #     + geom_smooth(aes(color = sampling_day),
  #                   method = lm, formula = y ~ log2(x+1), span = 1, 
  #                   se = FALSE,
  #                   show.legend = FALSE)
  #     + geom_point(aes(fill = sampling_day, shape = site), size = 2)
  #     + scale_shape_manual(values = c(21, 22, 23),
  #                          labels = site_lookup,
  #                          breaks = names(site_lookup))
  #     + scale_fill_manual(values = sampling_day_colors, labels = sampling_day_lookup, breaks = names(sampling_day_lookup))
  #     + scale_color_manual(values = sampling_day_colors)
  #     + scale_x_continuous(name = "Number of sequences", labels = scaleFUN0)
  #     + scale_y_continuous(name = "Number of ASVs in rarefaction", labels = scaleFUN0)
  #     + guides(color = "none",
  #              fill = guide_legend(order = 2, ncol = 1, title = "Sampling day", direction = "vertical",
  #                                  override.aes = list(color = "black", stroke = 1, shape = 21)),
  #              shape = "none")
  #     + facet_grid(cols = vars(sampling_time),
  #                  scales = "free_y",
  #                  space = "fixed",
  #                  labeller = labeller(sampling_time = sampling_time_lookup))
  #     + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
  #             panel.grid.minor = element_blank())
  #     + coord_cartesian(xlim = c(0, xlim), 
  #                       ylim = c(0, ylim),
  #                       expand = FALSE)
  #     + ggtitle(gtitle)
  #     + theme(plot.title = ggtext::element_markdown())
  #   )
  #   assign(paste0("g_curve_", gtitle), g, envir = .GlobalEnv, inherits = TRUE)
  #   rm(g)
  #   rm(gtitle)
  # }
}

#faceting by sampling day for each site:
{
  for(i in 1:length(temp_list1)){
    namevar <- eval(names(temp_list1)[i])
    gtitle <- site_lookup[grepl(gsub(" ", "|", namevar), names(site_lookup))]
    
    ylim <- temp_list1[[i]] %>%
      dplyr::select(num_ASVs) %>%
      deframe(.) %>%
      max(.)
    ylim <- ceiling(ylim/100)*100
    
    g <- print(
      ggplot(data = temp_list1[[i]],
             aes(x = seqdepth, y = num_ASVs, group = sample_id))
      + theme_bw()
      + geom_smooth(aes(color = sampling_time),
                    method = lm, formula = y ~ log2(x+1), span = 1, 
                    se = FALSE,
                    show.legend = FALSE)
      + geom_point(aes(fill = sampling_time, shape = site), size = 2)
      + scale_shape_manual(values = c(21, 22, 23),
                           labels = site_lookup,
                           breaks = names(site_lookup))
      + scale_fill_manual(values = sampling_time_colors, labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
      + scale_color_manual(values = sampling_time_colors)
      + scale_x_continuous(name = "Number of sequences", labels = scaleFUN0)
      + scale_y_continuous(name = "Number of ASVs in rarefaction", labels = scaleFUN0)
      + guides(color = "none",
               fill = guide_legend(order = 2, ncol = 1, title = "Sampling time", direction = "vertical",
                                   override.aes = list(color = "black", stroke = 1, shape = 21)),
               shape = "none")
      + facet_grid(cols = vars(sampling_day),
                   scales = "free_y",
                   space = "fixed",
                   labeller = labeller(sampling_day = sampling_day_lookup))
      + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
              panel.grid.minor = element_blank())
      + coord_cartesian(xlim = c(0, xlim), 
                        ylim = c(0, ylim),
                        expand = FALSE)
      + ggtitle(gtitle)
      + theme(plot.title = ggtext::element_markdown())
    )
    assign(paste0("g_curve_", namevar), g, envir = .GlobalEnv, inherits = TRUE)
    rm(g)
    rm(gtitle)
  }
}

# gpatch <- apropos("^g_curve_.*$", mode = "list") %>%
#   lapply(., get) %>%
#   setNames(., c(names(temp_list1), names(temp_list2))) %>%
#   purrr::reduce(., `%+%`) + patchwork::plot_layout(guides = "collect")

gpatch <- (g_curve_LB_seagrass + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())) / (g_curve_Tektite + theme(axis.text.x = element_blank(), strip.text.x = element_blank(), axis.title.x = element_blank())) / (g_curve_Yawzi + theme(axis.title.y = element_blank(), strip.text.x = element_blank())) + patchwork::plot_layout(guides = "collect")
gpatch <- gpatch + patchwork::plot_annotation(title = "Richness curves of USVI temporal samples",
                                              tag_level = "A")

if(!any(grepl("seawater", list.files(projectpath, pattern = "usvi_rarecurve.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_rarecurve_seawater-", Sys.Date(), ".png"),
         gpatch,
         width = 16, height = 10, units = "in")
}

#looking at the controls now:
{
  xlim2 <- max(usvi_rarefied_stats.df %>%
                dplyr::filter(grepl("extraction", sample_type)) %>%
                droplevels %>%
                dplyr::select(seqdepth) %>%
                deframe(.) %>%
                max(.))
  xlim2 <- round(xlim2/100, digits = -1)
  xlim2 <- xlim2*100

  for(i in 1:length(temp_list2)){
    namevar <- eval(names(temp_list2)[i])
    gtitle <- site_lookup[grepl(gsub(" ", "|", namevar), names(site_lookup))]
    
    ylim <- temp_list2[[i]] %>%
      dplyr::select(num_ASVs) %>%
      deframe(.) %>%
      max(.)
    ylim <- ceiling(ylim/100)*100
    
    g <- print(
      ggplot(data = temp_list2[[i]],
             aes(x = seqdepth, y = num_ASVs, group = sample_id))
      + theme_bw()
      + geom_smooth(aes(color = sample_id),
                    method = lm, formula = y ~ log2(x+1), span = 1, 
                    se = FALSE,
                    show.legend = FALSE)
      + geom_point(aes(fill = sample_id), shape = 21, color = "black", size = 2)
      + scale_x_continuous(name = "Number of sequences", labels = scaleFUN0)
      + scale_y_continuous(name = "Number of ASVs in rarefaction", labels = scaleFUN0)
      + guides(color = "none",
               fill = guide_legend(order = 2, ncol = 1, title = "Sample ID", direction = "vertical",
                                   override.aes = list(color = "black", stroke = 1, shape = 21)),
               shape = "none")
      + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
              panel.grid.minor = element_blank())
      + coord_cartesian(xlim = c(0, xlim2), 
                        ylim = c(0, ylim),
                        expand = FALSE)
      + ggtitle(gtitle)
      + theme(plot.title = ggtext::element_markdown())
    )
    assign(paste0("g_curve_", namevar), g, envir = .GlobalEnv, inherits = TRUE)
    rm(g)
    rm(gtitle)
  }
}


gpatch2 <- apropos("^g_curve_control.*$", mode = "list") %>%
  lapply(., get) %>%
  setNames(., names(temp_list2)) %>%
  # purrr::reduce(., `%+%`) + patchwork::plot_layout(guides = "collect")
  purrr::reduce(., `/`) + patchwork::plot_layout(guides = "collect")
gpatch2 <- (g_curve_control_extraction ) / (g_curve_control_pcr + (g_curve_control_seq + coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim), expand = FALSE) + theme(axis.title.y = element_blank()))) + patchwork::plot_layout(guides = "collect")
gpatch2 <- gpatch2 + patchwork::plot_annotation(title = "Richness curves of USVI control samples",
                                              tag_level = "A")
gpatch2
if(!any(grepl("controls", list.files(projectpath, pattern = "usvi_rarecurve.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_rarecurve_controls-", Sys.Date(), ".png"),
         gpatch2,
         width = 10, height = 10, units = "in")
}


# NMDS plot ---------------------------------------------------------------

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
                   by = c("label" = "sample_id"))


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
  g1 <- ggplot(data = nmds_asv_df,
               aes(x = NMDS1, y = NMDS2))
  g1 <- g1 + theme(text = element_text(size=14))
  g1 <- g1 + geom_point(data = (nmds_asv_df %>%
                                  dplyr::filter(type == "sample") %>%
                                  droplevels),
                        aes(fill = site, shape = sampling_time), color = "black",
                        size = 5, stroke = 2, alpha = 1)
  g1 <- g1 + scale_shape_manual(values = c(22, 21), labels = sampling_time_lookup, breaks = names(sampling_time_lookup))
  g1 <- g1 + scale_fill_manual(values = site_colors, labels = site_lookup, breaks = names(site_lookup))
  g1 <- g1 + theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
                   panel.background = element_blank(), panel.border = element_rect(fill = "NA", colour = "grey30"),
                   panel.grid = element_blank(),
                   legend.position = "bottom",
                   legend.key = element_blank(),
                   legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
                   legend.text = element_text(size = 12, colour = "grey30"))
  g1 <- g1 + guides(color = "none",
                    fill = guide_legend(order = 2, ncol = 1, title = "Site", direction = "vertical",
                                        override.aes = list(color = "black", stroke = 1, shape = 21, size = 2)),
                    shape = guide_legend(order = 1, ncol = 1, title = "Sampling time", direction = "vertical",
                                         override.aes = list(color = "black", stroke = 1, size = 2)))
  g1 <- g1 + coord_cartesian(ylim = nmds_lims[["ylim"]],
                             xlim = nmds_lims[["xlim"]],
                             expand = TRUE)
  g1
  }

# g2 <- g1 + geom_segment(aes(x = 0, y = 0,
#                             xend = NMDS1,
#                             yend = NMDS2),
#                         data = (nmds_asv_df %>%
#                                   dplyr::filter(type == "Vector") %>%
#                                   droplevels),
#                         arrow = grid::arrow(angle = 20, length = unit(0.05, "npc"),
#                                             ends = "last", type = "closed"),
#                         linewidth = 0.5, alpha = 0.5, colour = "grey30") +
#   geom_text(data = (nmds_asv_df %>%
#                                   dplyr::filter(type == "Vector") %>%
#                                   dplyr::mutate(across(c(NMDS1, NMDS2), ~.x * 1.1)) %>%
#                                   droplevels),
#                         aes(label = label, x = NMDS1, vjust = "outward", hjust = "center",
#                             y = NMDS2),
#                         check_overlap = TRUE,
#                         colour = "grey30", fontface = "bold")
# g2

if(!any(grepl("seawater", list.files(projectpath, pattern = "usvi_nmds.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_nmds_seawater-", Sys.Date(), ".png"),
         g1,
         width = 10, height = 10, units = "in")
}

# g3 <- g1 + geom_segment(aes(x = 0, y = 0,
#                             xend = NMDS1,
#                             yend = NMDS2),
#                         data = (nmds_sw_df %>%
#                                   dplyr::filter(type == "ASV") %>%
#                                   dplyr::slice_max(n = 10, order_by = abs(NMDS1)) %>%
#                                   droplevels),
#                         arrow = grid::arrow(angle = 20, length = unit(0.05, "npc"),
#                                             ends = "last", type = "closed"),
#                         linewidth = 0.5, alpha = 0.5, colour = "grey30") +
#   geom_text(data = (nmds_sw_df %>%
#                       dplyr::filter(type == "ASV") %>%
#                       # dplyr::slice_max(abs(NMDS1), n = 10) %>%
#                       dplyr::slice_max(n = 10, order_by = abs(NMDS1)) %>%
#                       dplyr::mutate(across(c(NMDS1, NMDS2), ~.x * 1.1)) %>%
#                       droplevels),
#             aes(label = label, x = NMDS1, vjust = "outward", hjust = "center",
#                 y = NMDS2),
#             check_overlap = TRUE,
#             colour = "grey30", fontface = "bold")
# g3


#correlate the bray-curtis NMDS of seawater samples with the FCM measurements?
#FCM sample from Metab_280 had weird Prochlorococcus, Picoeukaryotes, and unpigmented cells
usvi_sw_fcm_df <- ps_usvi %>%
  phyloseq::sample_data(.) %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  dplyr::filter(sample_id %in% rownames(usvi_asv.tbl)) %>%
  dplyr::filter(!grepl("Metab_280", sample_id)) %>%
  # dplyr::select(sample_id, sampling_time, sampling_day, site) %>%
#   dplyr::select(sample_id, sampling_time, sampling_day, site, contains("fcm")) %>%
#   dplyr::select(!c(contains("label"), contains("dna_"))) %>%
  tibble::column_to_rownames(., var = "sample_id") %>%
#   droplevels
# usvi_sw_fcm_df <- meta.seawater[!grepl("Metab_280", rownames(meta.seawater)),] %>%
# usvi_sw_fcm_df <- meta.seawater %>%
  dplyr::select(starts_with("fcm_")) %>%
  tidyr::drop_na(.) %>%
  vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
                                 autotransform = TRUE) %>%
  as.matrix(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  droplevels

usvi_sw_dist_df <- vegan::vegdist(usvi_asv.tbl, distance = "horn", binary = FALSE, upper = TRUE,
                                  autotransform = TRUE) %>%
  as.matrix(.) %>%
  as.data.frame(.) %>%
  dplyr::select(usvi_sw_fcm_df[["sample_id"]]) %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::filter(sample_id %in% usvi_sw_fcm_df[["sample_id"]]) %>%
  droplevels %>%
  tibble::column_to_rownames(var = "sample_id")
usvi_sw_fcm_df <- usvi_sw_fcm_df %>%
  tibble::column_to_rownames(var = "sample_id")

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

