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

#here are metabolites where we don't have LODs reported:
if(!exists("usvi_sus_metabolites_idx", envir = .GlobalEnv)){
  usvi_sus_metabolites_idx <- data.frame(metabolites = c("2'deoxyguanosine", "HMP", "adenosine", "inosine", "pyridoxine", "4-aminobenzoic acid"))
}


# Plot heatmap of metab concentrations ------------------------------------

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
    ggsave(paste0(projectpath, "/", "usvi_metab_g3_relconc_", i, "_", namevar1 , "-", Sys.Date(), ".png"),
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
  # {i <- 1
  namevar1 <- names(usvi_sig_metab_median_list[i])
  namevar2 <- namevar1 %>%
    gsub("LB_seagrass", "Lameshur Bay", .) %>%
    gsub("peak_photo", "afternoon", .)
  temp_df1 <- usvi_sig_metab_median_list[[namevar1]] %>%
    droplevels
  fig_height <- round(usvi_sig_metab_figh[[namevar1]]/2 + 2)
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
    
    temp_df2 <- temp_df1 %>%
      dplyr::mutate(median_norm = tidyr::replace_na(median_norm, 0)) %>%
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
      + ggtitle(paste0("Normalized concentration of significant metabolites in ", namevar2))
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
    assign(paste0("g4_metab_", i, "_", namevar1), temp_g3, envir = .GlobalEnv, inherits = TRUE)
    assign(paste0("g4_metab_", i, "_", namevar1, "_med_relabund"), temp_g3a, envir = .GlobalEnv, inherits = TRUE)
    assign(paste0("g4_metab_", i, "_", namevar1, "_norm_relabund"), temp_g3b, envir = .GlobalEnv, inherits = TRUE)
    
    # # # if(!any(grepl(paste0(namevar1, "relconc", Sys.Date(), collapse = "&"), list.files(projectpath, pattern = "usvi_metab_g4_sig_relconc_*.png")))){
    ggsave(paste0(projectpath, "/", "usvi_metab_g4_sig_relconc_", i, "_", namevar1 , "-", Sys.Date(), ".png"),
           temp_g3,
           scale = 1,
           width = 20, height = fig_height, units = "in")
    # # # }
    rm(list = apropos("temp_df.*", mode = "list"))
    rm(list = apropos("namevar.*", mode = "list"))
    rm(list = apropos("temp_breaks.*", mode = "list"))
    rm(list = apropos("temp_labels.*", mode = "list"))
    rm(list = apropos("temp_pal.*", mode = "list"))
    rm(list = apropos("temp_g3.*", mode = "list"))
  }
    # rm(fig_height)
}
