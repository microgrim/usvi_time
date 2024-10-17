# 04_taxonomy_summary.R

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


# Import color scheme -----------------------------------------------------


temporal_lookup <- list(sampling_time = c("dawn", "peak_photo", NA),
                        sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5", NA),
                        sampling_date = c(20210122, 20210123, 20210124, 20210125, 20210126, NA),
                        site = c("Tektite", "Yawzi", "LB_seagrass", NA),
                        sample_type = c("seawater", "control_extraction", "control_pcr", "control_seq"))


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
  dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                          grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                          .default = Phylum)) %>%
  
  droplevels



if(!exists("ps_usvi", envir = .GlobalEnv)){
  if(file.exists(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))){
    # ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))
    ps_usvi <- readr::read_rds(paste0(projectpath, "/", "usvi_prok_decontam_phyloseq", ".rds"))
  } else {
    cli::cli_alert_warning("Please process the USVI data through Phyloseq.")
    
  }
}

# Summarize sequencing effort ---------------------------------------------

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
  # dplyr::group_by(sample_type, site, sampling_time) %>%
  dplyr::group_by(sample_type, site) %>%
  dplyr::summarise(mean = mean(total_seqs),
                   max = max(total_seqs))

#note that the following samples were sequenced in a separate run:
#Metab_180, Metab_182, Metab_199, Metab_219, Metab_224, Metab_231, Metab_233, Metab_235
#two LB samples and two Tektite samples had relatively lower seq counts
#Metab_199, Metab_237
#Metab_231, Metab_224

if(!any(grepl("seqcount", list.files(projectpath, pattern = "usvi_.*.png")))){
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
}

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
keep <- c("Phylum", "Class", "Order", "Family", "Genus")

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
 
  assign(paste0("usvi_", namevar, ".df"), temp_df, envir = .GlobalEnv)
  rm(temp_df)
}

#now work backwards to keep the class, order, family related to the top5 genera
keep_tax <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

#look at the genera that constitute 10% or more across all biological samples
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

usvi_taxonomy_filtered_list <- list(usvi_phylum_filtered.df, usvi_class_filtered.df, usvi_order_filtered.df, usvi_family_filtered.df, usvi_filtered_genus.df) %>%
  setNames(., c("usvi_phylum_filtered.df", "usvi_class_filtered.df", "usvi_order_filtered.df", "usvi_family_filtered.df", "usvi_filtered_genus.df"))
if(!file.exists(paste0(projectpath, "/", "usvi_taxonomy_filtered_list", ".rds"))){
  readr::write_rds(usvi_taxonomy_filtered_list, 
                   paste0(projectpath, "/", "usvi_taxonomy_filtered_list", ".rds"))
}


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


if(!any(grepl("phylum", list.files(projectpath, pattern = "usvi_.*_abundance_stats.*.png")))){
  
  df <- ASV_env_phylum_summary.tbl %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    droplevels
  df2 <- ASV_env_phylum_summary.tbl %>%
    dplyr::filter(!grepl("seawater", sample_type)) %>%
    droplevels
  
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


taxonomy_colors <- dplyr::left_join(annotation_taxa_colors_list[["Genus"]] %>%
                                      tibble::enframe(., name = "Genus", value = "color"), 
                                    usvi_filtered_genus.df %>%
                                      dplyr::select(Genus, Domain, Phylum, Class, Order, Family) %>%
                                      dplyr::distinct(Genus, Domain, Phylum, Class, Order, Family, .keep_all = FALSE) %>%
                                      dplyr::mutate(across(c(Class, Order, Genus), ~ dplyr::case_when(is.na(.x) ~ "NA",
                                                                             .default = .x))),
                                      # dplyr::mutate(Class = dplyr::case_when(is.na(Class) ~ "NA",
                                      #                                        .default = Class)),
                                    by = join_by(Genus)) %>%
  droplevels %>%
  tidyr::unite("taxonomy_class", c(Domain, Phylum, Class), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_order", c(Domain, Phylum, Class, Order), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_family", c(Domain, Phylum, Class, Order, Family), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_genus", c(Domain, Phylum, Class, Order, Family, Genus), sep = ";", remove = TRUE, na.rm = TRUE) %>%
  # dplyr::distinct(color, .keep_all = TRUE) %>%
  dplyr::select(color, contains("taxonomy")) %>%
  tidyr::pivot_longer(., cols = !c(color),
                      names_to = "taxonomic_level",
                      values_to = "taxonomy") %>%
  dplyr::mutate(taxonomy = dplyr::case_when(grepl("NA;NA", taxonomy) ~ "NA",
                                            .default = taxonomy)) %>%
  droplevels %>%
  dplyr::ungroup(.) %>%
  tidyr::drop_na(.) %>%
  dplyr::arrange(taxonomy) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  dplyr::select(taxonomy, color) %>%
  tibble::deframe(.)

# taxonomy_colors_outline <- rep(c("black", "grey50", "grey90"), length(taxonomy_colors)/3) %>%
#   setNames(., names(taxonomy_colors))

taxonomy_colors_lookup <- data.frame(taxonomy = names(taxonomy_colors),
                                     first = stringr::str_split_i(names(taxonomy_colors), ";", 2),
                                     second = stringr::str_split_i(names(taxonomy_colors), ";", -1)) %>%
  tidyr::unite(label, c("first", "second"), sep = ";") %>%
  dplyr::select(label, taxonomy) %>%
  tibble::deframe(.)

keep <- c("Cyanobacteria", "Alphaproteobacteria", "Bacteroidota", 
          # "Proteobacteria",
          "Gammaproteobacteria", "Actinobacteriota", 
          "SAR406", "Verrucomicrobiota", "Thermoplasmatota")
drop <- c("Rhodothermia")
{
  temp_df1 <- usvi_class_filtered.df %>%
    dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_class_filtered.df[["Phylum"]]), value = TRUE)) %>%
    dplyr::filter(!(Class %in% grep(paste0(drop, collapse = "|"), unique(usvi_class_filtered.df[["Class"]]), value = TRUE))) %>%
    dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    droplevels %>%
    dplyr::arrange(desc(relabund), Class, taxonomy) %>%
    # dplyr::mutate(Phylum = factor(Phylum, levels = unique(.[["Phylum"]]))) %>%
    dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_colors))) %>%
    droplevels %>%
    dplyr::arrange(sampling_day, sampling_time, replicate) %>%
    dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
    dplyr::mutate(sample_order = factor(sample_order, levels = unique(.[["sample_order"]]))) %>%
    # dplyr::group_by(taxonomy, sample_type, site, sampling_day, sampling_time,replicate, sample_id, sample_order) %>%
    dplyr::distinct(sample_id, taxonomy, .keep_all = TRUE) %>%
    # dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
    droplevels
  
  temp_df2 <- usvi_order_filtered.df %>%
    dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_order_filtered.df[["Phylum"]]), value = TRUE)) %>%
    dplyr::filter(!(Class %in% grep(paste0(drop, collapse = "|"), unique(usvi_order_filtered.df[["Class"]]), value = TRUE))) %>%
    dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    # dplyr::mutate(Phylum = factor(Phylum, levels = unique(.[["Phylum"]]))) %>%
    dplyr::arrange(desc(relabund), taxonomy, Order) %>%
    dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_colors))) %>%
    droplevels %>%
    dplyr::arrange(sampling_day, sampling_time, replicate) %>%
    dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
    dplyr::mutate(sample_order = factor(sample_order, levels = unique(.[["sample_order"]]))) %>%
    dplyr::distinct(sample_id, taxonomy, .keep_all = TRUE) %>%
    # dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
    droplevels
  
  temp_df3 <- usvi_family_filtered.df %>%
    dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_family_filtered.df[["Phylum"]]), value = TRUE)) %>%
    dplyr::filter(!(Class %in% grep(paste0(drop, collapse = "|"), unique(usvi_family_filtered.df[["Class"]]), value = TRUE))) %>%
    dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::arrange(desc(relabund), Family, taxonomy) %>%
    # dplyr::mutate(Phylum = factor(Phylum, levels = unique(.[["Phylum"]]))) %>%
    dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_colors))) %>%
    droplevels %>%
    dplyr::arrange(sampling_day, sampling_time, replicate) %>%
    dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
    dplyr::mutate(sample_order = factor(sample_order, levels = unique(.[["sample_order"]]))) %>%
    dplyr::distinct(sample_id, taxonomy, .keep_all = TRUE) %>%
    # dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
    droplevels
  
  temp_df4 <- usvi_filtered_genus.df %>%
    dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_filtered_genus.df[["Phylum"]]), value = TRUE)) %>%
    dplyr::filter(!(Class %in% grep(paste0(drop, collapse = "|"), unique(usvi_filtered_genus.df[["Class"]]), value = TRUE))) %>%
    dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
    dplyr::filter(grepl("seawater", sample_type)) %>%
    dplyr::arrange(desc(relabund), Genus, taxonomy) %>%
    # dplyr::mutate(Phylum = factor(Phylum, levels = unique(.[["Phylum"]]))) %>%
    dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_colors))) %>%
    droplevels %>%
    dplyr::arrange(sampling_day, sampling_time, replicate) %>%
    dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
    dplyr::mutate(sample_order = factor(sample_order, levels = unique(.[["sample_order"]]))) %>%
    dplyr::distinct(sample_id, taxonomy, .keep_all = TRUE) %>%
    # dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
    droplevels
  
  usvi_relabund_taxonomy_list <- list(temp_df1, temp_df2, temp_df3, temp_df4) %>%
    setNames(., c("Class", "Order", "Family", "Genus"))
}

# temp_df3 %>%
#   dplyr::group_by(sample_id) %>%
#   dplyr::summarise(totalabund = sum(relabund, na.rm = TRUE)) %>%
#   dplyr::arrange(desc(totalabund))

{
  temp_df <- temp_df4 %>%
    dplyr::filter(relabund > 1) %>%
    droplevels
  print(ggplot(data = temp_df %>%
                 # dplyr::filter(grepl("LB", site)) %>%
                 droplevels, aes(y = relabund,
                                 x = sample_order_all,
                                 group = taxonomy,
                                 # color = taxonomy,
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
        # + scale_y_continuous(expand = expansion(mult = c(0,0.1)), transform = "pseudo_log")
        + scale_fill_manual(values = taxonomy_colors, 
                            breaks = names(taxonomy_colors),
                            labels = names(taxonomy_colors_lookup),
                            drop = TRUE)
  # + scale_color_manual(values = taxonomy_colors_outline, 
  #                     breaks = names(taxonomy_colors_outline),
  #                     labels = names(taxonomy_colors_lookup),
  #                     drop = TRUE)
        + facet_grid(scales = "free", space = "free",
                     labeller = labeller(site = site_lookup),
                     . ~ site)
                     # Class ~ site)
        + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                strip.text.x = element_text(angle = 0),
                strip.text.y = element_text(angle = 0))
        + labs(x = "sample",
               y = "Relative abundance (%)")
        + coord_cartesian(ylim = c(0,100), expand = FALSE)
        + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
                                     override.aes = list(stroke = 1, color = "black")),
                 color = "none")
  )
}

for(taxlevel in names(usvi_relabund_taxonomy_list)){
  namevar <- taxlevel
  temp_df <- usvi_relabund_taxonomy_list[[taxlevel]] %>%
    dplyr::filter(relabund > 5) %>%
    droplevels

  g2 <- (ggplot(data = temp_df %>%
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
        + theme_bw()
        + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
        # + scale_y_continuous(expand = expansion(mult = c(0,0.1)), transform = "pseudo_log")
        + scale_fill_manual(values = taxonomy_colors, 
                            breaks = names(taxonomy_colors),
                            labels = names(taxonomy_colors_lookup),
                            drop = TRUE)
        + facet_grid(scales = "free", space = "fixed",
                     labeller = labeller(site = site_lookup),
                     # . ~ site)
                     Class ~ site)
        + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                strip.text.x = element_text(angle = 0),
                strip.text.y = element_text(angle = 0))
        + labs(x = "sample",
               y = "Relative abundance (%)")
        # + coord_cartesian(ylim = c(0,100), expand = FALSE)
        + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
                                     override.aes = list(color = "black", stroke = 1)))
  )
  assign(paste0("g_barchart_", namevar), g2, envir = .GlobalEnv)
  ggsave(paste0(projectpath, "/", "usvi_barchart_", namevar, "-", Sys.Date(), ".png"),
         g2,
         width = 14, height = 8, units = "in")
}



