# 03_diversity_curves.R

#make rarefaction curves, etc. for ASVs in the USVI temporal series


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


# Prepare alpha diversity estimates ---------------------------------------


#note that the following samples were sequenced in a separate run:
#Metab_180, Metab_182, Metab_199, Metab_219, Metab_224, Metab_231, Metab_233, Metab_235

usvi_seqdepth <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
# usvi_seqdepth <- phyloseq::otu_table(ps_usvi) %>%
  as.data.frame %>%
  t() %>%
  rowSums(.)
#two LB samples and two Tektite samples had relatively lower seq counts
#Metab_199, Metab_237
#Metab_231, Metab_224

usvi_seqdepth[c("Metab_180", "Metab_182", "Metab_199", "Metab_219", "Metab_224", "Metab_231", "Metab_233", "Metab_235")]

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



