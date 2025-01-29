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
                                   label = c("Dawn", "Afternoon")) %>%
  tibble::deframe(.)
sampling_time_colors <- pals::ocean.haline(n = length(sampling_time_lookup)) %>%
  setNames(., names(sampling_time_lookup))
sampling_day_lookup <- data.frame(sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5"),
                                  label = c("20210122", "20210123", "20210124", "20210125", "20210126")) %>%
  tibble::deframe(.)
sampling_day_colors <- pals::ocean.thermal(n = length(sampling_day_lookup)) %>%
  setNames(., names(sampling_day_lookup))

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

  
usvi_prok_asvs.df <- usvi_prok_asvs.df %>%
  dplyr::filter(asv_id %in% usvi_prok_asvs.taxa[["asv_id"]]) %>%
#   droplevels
# usvi_prok_asvs.df %>%
  dplyr::filter(grepl("Metab_", sample_ID)) %>%
  dplyr::filter(!grepl("Metab_K", sample_ID)) %>%
  droplevels

usvi_prok_asvs.df %>%
  dplyr::filter(grepl("Metab_", sample_ID)) %>%
  dplyr::filter(!grepl("Metab_K", sample_ID)) %>%
  dplyr::filter(counts > 0) %>%
  dplyr::summarise(num_seqs = sum(counts), num_samples = length(unique(.$sample_ID)), num_ASVs = length(unique(.$asv_id)))
# num_seqs num_samples num_ASVs
# <dbl>       <int>    <int>
#   1  7854297          90    12350

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
  phyloseq::subset_samples(., sample_type == "seawater") %>%
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

usvi_seq_summary.df %>%
  dplyr::group_by(site, sampling_time) %>%
  dplyr::summarise(across(contains("total_seqs"), list(mean = mean, min = min, max = max, median = median))) %>%
  droplevels
# sampling_time total_seqs_mean total_seqs_min total_seqs_max total_seqs_median
#   1 dawn                   85480.          13942         118596             85640
# 2 peak_photo             89229.          18535         114556             93335


#how many ASVs remain in here after filtering for non-zero ASVs in seawater samples only?

#across the dataset:
ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  dplyr::summarise(total_ASVs = length(asv_id))
# total_ASVs
# 1      12350

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
range(temp_df[["total_ASVs"]])
# [1] 169 827
usvi_seq_summary.df %>%
  # dplyr::group_by(sample_type, site, sampling_time) %>%
  dplyr::group_by(sample_type, site) %>%
  dplyr::summarise(mean = mean(total_seqs),
                   max = max(total_seqs))
usvi_seq_summary.df %>%
  dplyr::summarise(num_seqs = sum(total_seqs))
# num_seqs
# <dbl>
#   1  7854297

#are sequencing depths significantly different between sites sand sampling times?
temp_list <- usvi_seq_summary.df %>%
  dplyr::select(sample_id, site, sampling_time, total_seqs) %>%
  split(., f = interaction(.$site, .$sampling_time)) %>%
  map(., ~.x %>%
        dplyr::select(total_seqs) %>%
        tibble::deframe(.))

kruskal.test(temp_list)
# Kruskal-Wallis rank sum test
# 
# data:  temp_list
# Kruskal-Wallis chi-squared = 4.8814, df = 5, p-value = 0.4305


usvi_seq_summary.tsv <- usvi_seq_summary.df %>%
  dplyr::arrange(site, sampling_day, sampling_time, sample_order) %>%
  dplyr::select(sample_id, site, sampling_day, sampling_time, sample_order, total_seqs) %>%
  dplyr::mutate(sampling_time = recode(sampling_time, !!!sampling_time_lookup)) %>%
  # dplyr::mutate(sampling_time = dplyr::case_when(grepl("peak", sampling_time) ~ "afternoon", .default = sampling_time)) %>%
  dplyr::mutate(site = recode(site, !!!site_lookup)) %>%
  droplevels

if(!file.exists(paste0(projectpath, "/", "usvi_seq_summary.tsv", ".gz"))){
  
  readr::write_delim(usvi_seq_summary.tsv, paste0(projectpath, "/", "usvi_seq_summary.tsv", ".gz"),
                     delim = "\t", col_names = TRUE)
  
}


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


#make relative abundance tables at specific taxonomic levels
keep <- c("Phylum", "Class", "Order", "Family", "Genus")

for(taxlevel in keep){
  namevar <- stringr::str_to_lower(taxlevel)
  temp_df <- phyloseq::transform_sample_counts(ps_usvi, function(x) x/sum(x)*100) %>%
    phyloseq::subset_samples(., sample_type == "seawater") %>%
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
    #don't use the filled-in taxonomy tble for this!
    usvi_prok_asvs.taxa.mod.df <- usvi_prok_asvs.taxa %>%
      dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                              grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                              .default = Phylum)) %>%
      dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
      droplevels
    
    usvi_genus_taxonomy.df <- phyloseq::transform_sample_counts(ps_usvi, function(x) x/sum(x)*100) %>%
        phyloseq::subset_samples(., sample_type == "seawater") %>%
        phyloseq::filter_taxa(., function(x) sum(x) > 0, TRUE) %>%
        phyloseq::otu_table(.) %>%
        as.data.frame %>%
        dplyr::mutate(TotAbund = rowSums(.)) %>%
        dplyr::arrange(desc(TotAbund)) %>%
        dplyr::select(-TotAbund) %>%
        tibble::rownames_to_column(var = "asv_id") %>%
        otu_to_taxonomy(dataset = .,
                        level = "Genus",
                        taxonomy = usvi_prok_asvs.taxa.mod.df) %>%
        tidyr::pivot_longer(.,
                            cols = !c(1:all_of("Genus")),
                            names_to = "sample_id",
                            values_to = "relabund") %>%
        tidyr::unite("taxonomy", c(all_of(keep_tax)), sep = ";", remove = FALSE) %>%
      dplyr::relocate(taxonomy) %>%
      droplevels
    
    annotation_archaea_list <- usvi_genus_taxonomy.df %>%
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
    
      temp_list1 <- annotation_archaea_list %>%
        map(., ~.x %>%
              tibble::enframe(., name = "taxonomy_string", value = "classification") %>%
              tidyr::drop_na()) %>%
        map(., ~.x %>%
              dplyr::left_join(., annotation_archaea_colors, by = join_by(taxonomy_string))) %>%
        bind_rows(., .id = "taxonomic_level")
    
    #now bacteria:
    annotation_bacteria_list <- usvi_genus_taxonomy.df %>%
      dplyr::select(Domain, Phylum, Class, Order, Family, Genus, taxonomy) %>%
      dplyr::filter(!grepl("Archaea", Domain)) %>%
      dplyr::distinct(., .keep_all = TRUE) %>%
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
    
    bacteria_options <- c(brewer.blues, brewer.gnbu,  brewer.oranges, brewer.purples,  ocean.solar, 
                          ocean.haline, ocean.ice, ocean.oxy, brewer.greens, ocean.deep, ocean.dense, 
                          ocean.algae, ocean.matter, ocean.turbid, ocean.speed, ocean.amp, 
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
      # if(length(temp_list) > 2){ #is there is more than one class in this phylum in the MAGs
        if(length(temp_list) > 4){ #is there is more specific taxonomic levels for this phylum
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
    # rm(bacteria_options)
    
    
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

# Summarise taxonomic diversity -------------------------------------------


#how many different taxonomic lineages?
temp_df <- phyloseq::tax_table(ps_usvi) %>%
  as.data.frame %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  droplevels
temp_df %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  dplyr::filter(if_any(c("Domain":"Species"), ~!is.na(.x))) %>%
  dplyr::filter(if_any(c("Domain":"Species"), ~!grepl("Other", .x, ignore.case = TRUE))) %>%
  apply(., 2, function(x) length(unique(x)))
# Domain  Phylum   Class   Order  Family   Genus Species 
# 2      63     145     338     564     994    1109 

#how many non-propagated?
temp_df <- usvi_prok_asvs.taxa %>% 
  dplyr::select(asv_id, Domain:Species) %>%
  dplyr::filter(if_any(c("Domain":"Species"), ~!is.na(.x))) %>%
  droplevels
temp_df %>%
  dplyr::filter(if_any(c("Domain":"Species"), ~!grepl("Other", .x, ignore.case = TRUE))) %>%
  apply(., 2, function(x) length(unique(x)))
# Domain  Phylum   Class   Order  Family   Genus Species 
# 2      63     114     257     347     593     202 


temp_df2 <- usvi_taxonomy_filtered_list[["usvi_filtered_genus.df"]] %>%
  dplyr::group_by(Domain, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(total_abund = sum(relabund)) %>%
  dplyr::distinct(.) %>%
  dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(across(c(Domain, Phylum, Class, Order, Family, Genus), ~dplyr::case_when(!is.na(.x) ~ .x, .default = "Other taxa"))) %>%
  droplevels

temp_df2 %>%
  dplyr::filter(!grepl("Other", Domain)) %>%
  dplyr::summarise(total_abund = sum(total_abund))
temp_df2 %>%
  dplyr::filter(grepl("Other", Domain)) %>%
  dplyr::summarise(total_abund = sum(total_abund))


temp_df2 %>%
  dplyr::filter(if_any(c("Domain":"Genus"), ~!grepl("Other", .x, ignore.case = TRUE))) %>%
    dplyr::select(Domain:Genus) %>%
  apply(., 2, function(x) length(unique(x)))
# Domain Phylum  Class  Order Family  Genus 
# 2     11     12     28     43     56 


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

grep("Candidatus Nitro", names(annotation_taxa_colors_list[["Genus"]]), value = TRUE)

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

pal.bands(phylum_colors)

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

annotation_taxa_colors <- annotation_taxa_colors_list %>%
  bind_rows(., .id = "spec") %>%
  tidyr::pivot_longer(., cols = !c("spec"),
                      names_to = "taxonomy",
                      values_to = "color") %>%
  droplevels %>%
  dplyr::mutate(spec = factor(spec, levels = (c("Domain", "Phylum", "Class", "Order", "Family", "Genus")))) %>%
  dplyr::arrange(spec) %>%
  dplyr::distinct(taxonomy, color, .keep_all = TRUE) %>%
  droplevels
# annotation_taxa_colors %>%
#   dplyr::filter(grepl("Phylum", spec)) %>%
#   dplyr::select(taxonomy, color) %>%
#   droplevels %>%
#   dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
#   tibble::deframe(.) %>%
#   pal.bands
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
usvi_filtered_genus.df %>%
  dplyr::filter(!grepl("Other", taxonomy)) %>%
  dplyr::filter(relabund > 0) %>%
  droplevels %>%
  # dplyr::group_by(sample_id) %>% dplyr::summarise(total_relabund = sum(relabund)) %>% dplyr::arrange(total_relabund) #min: 89.3%
  # dplyr::summarise(num_genera = length(taxonomy)/90) #56genera in here
  dplyr::distinct(taxonomy, .keep_all = FALSE)

# usvi_relabund_taxonomy_list[["Genus"]] %>%
#   dplyr::filter(!grepl("Other", taxonomy)) %>%
#   # dplyr::filter(relabund > 0) %>%
#   droplevels %>%
#   # dplyr::group_by(sample_id) %>% dplyr::summarise(total_relabund = sum(relabund)) %>% dplyr::arrange(total_relabund) #min: 89.3%
#   dplyr::summarise(num_genera = length(taxonomy)/90) #52 here
#   # dplyr::distinct(taxonomy, .keep_all = FALSE) #52 in here

taxonomy_colors.df <- dplyr::inner_join(annotation_taxa_colors %>%
                                      tidyr::drop_na(color) %>%
                                      droplevels, 
                                    usvi_filtered_genus.df %>%
                                      dplyr::select(Genus, Domain, Phylum, Class, Order, Family) %>%
                                      dplyr::distinct(Genus, Domain, Phylum, Class, Order, Family, .keep_all = FALSE) %>%
                                      dplyr::mutate(across(c(Class, Order, Genus), ~ dplyr::case_when(is.na(.x) ~ "NA",
                                                                                                      .default = .x))),
                                    # dplyr::mutate(Class = dplyr::case_when(is.na(Class) ~ "NA",
                                    #                                        .default = Class)),
                                    by = join_by("taxonomy" == "Genus")) %>%
  droplevels %>%
  dplyr::mutate(Genus = taxonomy) %>% dplyr::relocate(Genus, .after = "Family") %>% dplyr::rename(value = "taxonomy") %>%
  tidyr::unite("taxonomy_phylum", c(Domain, Phylum), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_class", c(Domain, Phylum, Class), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_order", c(Domain, Phylum, Class, Order), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_family", c(Domain, Phylum, Class, Order, Family), sep = ";", remove = FALSE, na.rm = TRUE) %>%
  tidyr::unite("taxonomy_genus", c(Domain, Phylum, Class, Order, Family, Genus), sep = ";", remove = TRUE, na.rm = TRUE) %>%
  dplyr::select(color, contains("taxonomy")) %>%
  tidyr::pivot_longer(., cols = !c(color),
                      names_to = "taxonomic_level",
                      values_to = "taxonomy") %>%
  dplyr::mutate(taxonomy = dplyr::case_when(grepl("NA;NA", taxonomy) ~ "Other taxa",
                                            grepl("NA", taxonomy) ~ "Other taxa",
                                            .default = taxonomy)) %>%
  droplevels %>%
  dplyr::ungroup(.) %>%
  tidyr::drop_na(.) %>%
  dplyr::arrange(taxonomy, taxonomic_level) %>%
  dplyr::distinct(taxonomy, .keep_all = TRUE) %>%
  droplevels


taxonomy_colors.df %>%
  dplyr::select(taxonomy, color) %>%
  tibble::deframe(.) %>%
  pal.bands()




#apply a slgiht color filter to each phylum/class member?
# pal.bands(muted(temp_colors, l = 10, c = 70)) #adjust l and c; lower values of l make darker versions, higher values of c make it more intense
{
  # taxonomy_colors <- bind_rows(usvi_relabund_taxonomy_list[["Genus"]]) %>%
  #   dplyr::select(sample_id, c(Domain:taxonomy), relabund) %>%
  #   tidyr::pivot_wider(., id_cols = c(Domain:taxonomy),
  #                      names_from = sample_id,
  #                      values_from = relabund) %>%
  #   dplyr::mutate(total_abund = across(contains("Metab_")) %>% purrr::reduce(`+`)) %>%
  #   dplyr::relocate(total_abund, .after = "taxonomy") %>%
  #   dplyr::arrange(taxonomy) %>%
  #   dplyr::select(c(Domain:taxonomy), total_abund) %>%
  #   dplyr::distinct(., .keep_all = TRUE) %>%
  #   dplyr::left_join(., taxonomy_colors.df, by = join_by(taxonomy)) %>% droplevels %>%
  #   dplyr::select(taxonomy, color) %>%
  #   tibble::deframe(.)
  # 
  # 
  # # temp_colors <- taxonomy_colors[grepl("Alphaproteobacteria", names(taxonomy_colors))]
  # temp_colors <- taxonomy_colors[grepl("Marini", names(taxonomy_colors))]
  # 
  # 
  # y <- seq(10, 100, 10)
  # temp_pal <- temp_colors[1]
  # for(i in y){
  #   temp_pal <- c(temp_pal, muted(temp_colors[1], l = 30, c = i))
  # }
  # pal.bands(temp_pal)
  # temp_colors <- taxonomy_colors[grepl("Alphaproteobacteria", names(taxonomy_colors))]
  # y <- seq(10, 10*length(temp_colors), 10)
  # # j <- seq(10*length(temp_colors), 10, -10)
  # j <- rep_len(c(80, 50, 30), length(temp_colors))
  # temp_pal <- NULL
  # for(i in seq_len(length(y))){
  #   temp_pal <- c(temp_pal, muted(temp_colors[1], l = j[i], c = y[i]))
  # }
  # pal.bands(temp_pal)
}




# pal.bands(F_desaturate(phylum_colors, 80))


# taxonomy_rank_class_idx <- usvi_filtered_genus.df %>% #there are 12 classes here
taxonomy_rank_class_idx <- usvi_relabund_taxonomy_list[["Genus"]] %>% #there are 8 classes here that are of itnerest
  dplyr::group_by(Class) %>%
  # dplyr::summarise(mean_abund = mean(relabund)) %>%
  dplyr::summarise(mean_abund = max(relabund)) %>%
  dplyr::arrange(desc(mean_abund)) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(Class) %>%
  # dplyr::select(Class, taxonomy) %>%
  droplevels %>% tibble::deframe(.)
  # droplevels

# taxonomy_rank_genus_idx <- usvi_filtered_genus.df %>% #56 genera
taxonomy_rank_genus_idx <- usvi_relabund_taxonomy_list[["Genus"]] %>% #52
  dplyr::ungroup(.) %>%
  dplyr::distinct(Class, taxonomy) %>%
  dplyr::mutate(Class = factor(Class, levels = taxonomy_rank_class_idx)) %>%
  dplyr::arrange(Class) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
  dplyr::select(taxonomy) %>%
  droplevels %>% tibble::deframe(.)

taxonomy_filt_genus_colors.df <- c(phylum_colors, "#449C81") %>%
  setNames(., c(names(phylum_colors), "Archaea;Thermoplasmatota"))
names(taxonomy_filt_genus_colors.df)[1] <- "Archaea;Crenarchaeota"
taxonomy_filt_genus_colors.df <- taxonomy_filt_genus_colors.df %>%
  tibble::enframe(name = "first_taxonomy", value = "color")

taxonomy_filt_genus_colors <- bind_rows(usvi_relabund_taxonomy_list[["Genus"]]) %>%
  dplyr::select(sample_id, c(Domain:taxonomy), relabund) %>%
  tidyr::pivot_wider(., id_cols = c(Domain:taxonomy),
                     names_from = sample_id,
                     values_from = relabund) %>%
  dplyr::mutate(total_abund = across(contains("Metab_")) %>% purrr::reduce(`+`)) %>%
  dplyr::relocate(total_abund, .after = "taxonomy") %>%
  dplyr::arrange(desc(total_abund), taxonomy) %>%
  dplyr::select(c(Domain:taxonomy), total_abund) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::mutate(first_taxonomy = paste0(Domain, ";", Phylum)) %>%
  dplyr::left_join(., taxonomy_filt_genus_colors.df, by = join_by(first_taxonomy), relationship = "many-to-many", multiple = "all") %>% droplevels %>%
  split(., f = .$Phylum) %>%
  map(., ~.x %>%
        dplyr::select(taxonomy, color) %>%
        tibble::deframe(.)) %>%
  map(., ~F_desaturate(.x, 80)) %>%
  purrr::reduce(c)

taxonomy_filt_genus_colors <- taxonomy_filt_genus_colors[match(taxonomy_rank_genus_idx, names(taxonomy_filt_genus_colors))]

taxonomy_filt_genus_colors_lookup <- data.frame(taxonomy = names(taxonomy_filt_genus_colors),
                                                first = stringr::str_split_i(names(taxonomy_filt_genus_colors), ";", 2),
                                                second = stringr::str_split_i(names(taxonomy_filt_genus_colors), ";", -1)) %>%
  dplyr::mutate(second = dplyr::case_when(grepl("SAR11 clade", taxonomy) ~ paste0("SAR11 ", second),
                                          .default = second)) %>% #replace all of the "Alphaproteobacteria;Clade " names that are supposed to be SAR11 
  tidyr::unite(label, c("first", "second"), sep = ";") %>%
  dplyr::select(label, taxonomy) %>%
  tibble::deframe(.)

temp_df1 <- usvi_prok_filled.taxa.df %>% 
  dplyr::select(asv_id, Domain:Species) %>%
  droplevels %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(if_any(c("Family":"Genus"), ~.x %in% stringr::str_split_i(taxonomy_rank_genus_idx, ";", -1))) %>%
  # dplyr::filter(Genus %in% stringr::str_split_i(taxonomy_rank_genus_idx, ";", -1)) %>%
  tidyr::unite(col = "taxonomy", c("Domain":"Genus"), remove = FALSE, sep = ";") %>%
  dplyr::filter(taxonomy %in% taxonomy_rank_genus_idx) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = taxonomy_rank_genus_idx)) %>%
  tidyr::drop_na(taxonomy) %>%
  dplyr::arrange(taxonomy) %>%
  droplevels
# 
# temp_df3 <- temp_df2 %>%
#   dplyr::anti_join(temp_df, by = join_by(asv_id)) %>%
#   droplevels

usvi_filtered_taxonomy_rank_genus_relabund.df <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::otu_table(.) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = -c("asv_id"),
                      names_to = "sample_id",
                      values_to = "abundance") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(rel_abund = relabund(abundance)) %>%
  dplyr::mutate(abundance = dplyr::case_when(abundance > 0 ~ abundance,
                                             .default = NA)) %>%
  tidyr::drop_na(abundance)  %>%
  # dplyr::filter(asv_id %in% temp_df1[["asv_id"]]) %>%
  dplyr::right_join(temp_df1, by = join_by(asv_id), relationship = "many-to-many", multiple = "all") %>%
  droplevels

# length(unique(usvi_filtered_taxonomy_rank_genus_relabund.df$Genus))
## 52 unique genera

#these are the 56 ASVs belonging to the 52 abundant genera, that were 1% or more in at least one sample
temp_idx <- (usvi_filtered_taxonomy_rank_genus_relabund.df %>%
               dplyr::ungroup(.) %>%
               dplyr::filter(rel_abund >= 1) %>%
               dplyr::distinct(asv_id, taxonomy, .keep_all = FALSE) %>%
               dplyr::select(taxonomy, asv_id) %>%
               tibble::deframe(.))

usvi_filtered_taxonomy_rank_genus_relabund.df %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(total_relabund = sum(rel_abund)) %>%
  dplyr::arrange(total_relabund) #min: 88.8%





# taxonomy_colors <- taxonomy_colors.df %>%
#   dplyr::select(taxonomy, color) %>%
#   tibble::deframe(.)



temp_df <- usvi_relabund_taxonomy_list[["Genus"]] %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = names(taxonomy_filt_genus_colors))) %>%
  droplevels
temp_df %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(total_relabund = sum(relabund)) %>%
  dplyr::arrange(total_relabund) #min: 88.8%

g2 <- print(ggplot(data = temp_df %>%
                     droplevels, aes(y = relabund, x = sample_order_all, group = taxonomy, fill = taxonomy))
            + geom_bar(width = 0.90, show.legend = TRUE, position = "stack", stat = "identity")
            + geom_bar(color = "black", width = 0.90, show.legend = FALSE, position = "stack", stat = "identity")
            + theme_bw()
            + scale_x_continuous(name = "Sample", n.breaks = 6)
            + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Relative abundance (%)")
            + scale_fill_manual(values = taxonomy_filt_genus_colors, 
                                # breaks = names(taxonomy_filt_genus_colors),
                                breaks = taxonomy_filt_genus_colors_lookup,
                                labels = names(taxonomy_filt_genus_colors_lookup),
                                drop = TRUE)
            + facet_grid(scales = "free", space = "free",
                         labeller = labeller(site = site_lookup),
                         . ~ site)
            # Class ~ site)
            + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                    strip.text.x = element_text(angle = 0),
                    strip.text.y = element_text(angle = 0))
            + coord_cartesian(ylim = c(0,100), expand = FALSE)
            + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
                                         override.aes = list(stroke = 1, color = "black"), position = "right"),
                     color = "none")
)
# length(unique(names(temp_idx)))
# #30

#this is the relative abundace of the 30 genera that were 1% or more in any sample
temp_df2 <- temp_df %>%
  # temp_df2 <- usvi_relabund_taxonomy_list[["Genus"]] %>%
  # dplyr::filter(relabund >= 1) %>%
  # dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
  dplyr::filter(taxonomy %in% names(temp_idx)) %>%
  droplevels
# 
# #minimum relative abundance of the 30 genera in a sample: 82.4%
# temp_df2 %>%
#   dplyr::group_by(sample_id) %>%
#   dplyr::summarise(total_abund = sum(relabund)) %>%
#   dplyr::arrange(total_abund)

# #do you want to instead look at the ASVs that were 1% or more in any sample, belonging to these top 30 genera?
# temp_df2 <- usvi_filtered_taxonomy_rank_genus_relabund.df %>%
#   dplyr::filter(asv_id %in% temp_idx) %>%
#   droplevels %>%
#   dplyr::group_by(sample_id, taxonomy) %>%
#   dplyr::summarise(relabund = sum(rel_abund, na.rm = TRUE)) %>%
#   dplyr::left_join(., usvi_relabund_taxonomy_list[["Genus"]] %>%
#                      dplyr::select(starts_with("sampl"), site, replicate) %>%
#                      dplyr::distinct(.) %>%
#                      droplevels,
#                    by = join_by(sample_id), relationship = "many-to-many", multiple = "all") %>%
#   droplevels
# temp_df2 %>%
#   dplyr::group_by(sample_id) %>%
#   dplyr::summarise(total_abund = sum(relabund)) %>%
#   dplyr::arrange(total_abund) #minimum relative abundance of the top 30 genera' ASVs that were at least 1% in any sample: 68.9%

g2_oneperc <- print(ggplot(data = temp_df2 %>%
                             droplevels, aes(y = relabund, x = sample_order_all, group = taxonomy, fill = taxonomy))
                    + geom_bar(width = 0.90, show.legend = TRUE, position = "stack", stat = "identity")
                    + geom_bar(color = "black", width = 0.90, show.legend = FALSE, position = "stack", stat = "identity")
                    + theme_bw()
                    + scale_x_continuous(name = "Sample", n.breaks = 6)
                    + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Relative abundance (%)")
                    + scale_fill_manual(values = taxonomy_filt_genus_colors, 
                                        # breaks = names(taxonomy_filt_genus_colors),
                                        breaks = taxonomy_filt_genus_colors_lookup,
                                        labels = names(taxonomy_filt_genus_colors_lookup),
                                        drop = TRUE)
                    + facet_grid(scales = "free", space = "free",
                                 labeller = labeller(site = site_lookup),
                                 . ~ site)
                    # Class ~ site)
                    + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                            strip.text.x = element_text(angle = 0),
                            strip.text.y = element_text(angle = 0))
                    + coord_cartesian(ylim = c(0,100), expand = FALSE)
                    + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
                                                 override.aes = list(stroke = 1, color = "black"), position = "right"),
                             color = "none")
)

if(!any(grepl("usvi_barchart_", list.files(projectpath, pattern = "usvi_.*.png")))){
  ggsave(paste0(projectpath, "/", "usvi_barchart_filtered_genus", "-", Sys.Date(), ".png"),
         g2_oneperc,
         width = 16, height = 10, units = "in")
  ggsave(paste0(projectpath, "/", "usvi_barchart_genus", "-", Sys.Date(), ".png"),
         g2,
         width = 16, height = 10, units = "in")
}



# #which taxa were in filtered genus df, but not 1% or more abundant?
# setdiff((temp_df %>% dplyr::distinct(taxonomy) %>%  unlist %>% as.character), (temp_df2 %>% dplyr::distinct(taxonomy) %>%  unlist %>% as.character))
# #18 of them


{
  # temp_df2 <- usvi_taxonomy_filtered_list[["usvi_filtered_genus.df"]] %>%
  #   dplyr::group_by(Domain, Phylum, Class, Order, Family, Genus) %>%
  #   dplyr::summarise(total_abund = sum(relabund)) %>%
  #   dplyr::distinct(.) %>%
  #   dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::mutate(across(c(Domain, Phylum, Class, Order, Family, Genus), ~dplyr::case_when(!is.na(.x) ~ .x, .default = "Other taxa"))) %>%
  #   droplevels
  # 
  # taxonomy_colors <- rev(c("Class", "Order", "Family", "Genus")) %>%
  #   map(., ~usvi_relabund_taxonomy_list[[.x]] %>%
  #         bind_rows(., .id = "level") %>%
  #         dplyr::select(c(Domain:taxonomy)) %>%
  #         dplyr::distinct(., .keep_all = TRUE) %>%
  #         # droplevels)
  #         droplevels) %>%
  #   purrr::reduce(full_join) %>% droplevels %>%
  #   dplyr::left_join(., temp_df2 %>%
  #                      dplyr::select(Genus, total_abund), by = join_by("Genus" == "Genus")) %>%
  #   dplyr::left_join(., temp_df2 %>%
  #                      dplyr::select(Family, total_abund), by = join_by("Genus" == "Family")) %>%
  #   dplyr::left_join(., temp_df2 %>%
  #                      dplyr::select(Order, total_abund), by = join_by("Genus" == "Order")) %>%
  # dplyr::left_join(., temp_df2 %>%
  #                    dplyr::select(Class, total_abund), by = join_by("Genus" == "Class")) %>%
  #     dplyr::mutate(total_abund = across(contains("total_abund")) %>% purrr::reduce(coalesce)) %>%
  #   dplyr::select(-c(ends_with(".x"), ends_with(".y"))) %>%
  #   dplyr::distinct(., .keep_all = TRUE) %>%
  #     dplyr::left_join(., taxonomy_colors.df %>%
  #                        dplyr::distinct(color, taxonomy), by = join_by(taxonomy)) %>% droplevels %>%
  #   dplyr::arrange(taxonomy, desc(total_abund)) %>%
  #   # droplevels
  #     dplyr::select(taxonomy, color) %>%
  #     tibble::deframe(.)
  # 
  # 
  # taxonomy_colors_lookup <- data.frame(taxonomy = names(taxonomy_colors),
  #                                      first = stringr::str_split_i(names(taxonomy_colors), ";", 2),
  #                                      second = stringr::str_split_i(names(taxonomy_colors), ";", -1)) %>%
  #   tidyr::unite(label, c("first", "second"), sep = ";") %>%
  #   dplyr::select(label, taxonomy) %>%
  #   tibble::deframe(.)


# keep <- c("Cyanobacteria", "Alphaproteobacteria", "Bacteroidota", 
#           # "Proteobacteria",
#           "Gammaproteobacteria", "Actinobacteriota", 
#           "SAR406", "Verrucomicrobiota", "Thermoplasmatota")
# drop <- c("Rhodothermia")
# 
# 
# for(taxlevel in names(usvi_relabund_taxonomy_list)){
#   namevar <- taxlevel
#   temp_df <- usvi_relabund_taxonomy_list[[taxlevel]] %>%
#     dplyr::filter(relabund > 5) %>%
#     droplevels
# 
#   g2 <- (ggplot(data = temp_df %>%
#                  # dplyr::filter(grepl("LB", site)) %>%
#                  droplevels, aes(y = relabund,
#                                  x = sample_order_all,
#                                  group = taxonomy,
#                                  fill = taxonomy))
#         + geom_bar(width = 0.90, show.legend = TRUE,
#                    # position = "dodge",
#                    position = "stack",
#                    stat = "identity")
#         + geom_bar(color = "black", width = 0.90, show.legend = FALSE,
#                    # position = "dodge",
#                    position = "stack",
#                    stat = "identity")
#         + theme_bw()
#         + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#         # + scale_y_continuous(expand = expansion(mult = c(0,0.1)), transform = "pseudo_log")
#         + scale_fill_manual(values = taxonomy_colors, 
#                             breaks = names(taxonomy_colors),
#                             labels = names(taxonomy_colors_lookup),
#                             drop = TRUE)
#         + facet_grid(scales = "free", space = "fixed",
#                      labeller = labeller(site = site_lookup),
#                      # . ~ site)
#                      Class ~ site)
#         + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
#                 strip.text.x = element_text(angle = 0),
#                 strip.text.y = element_text(angle = 0))
#         + labs(x = "sample",
#                y = "Relative abundance (%)")
#         # + coord_cartesian(ylim = c(0,100), expand = FALSE)
#         + guides(fill = guide_legend(order = 2, ncol = 1, title = "Taxonomy", direction = "vertical",
#                                      override.aes = list(color = "black", stroke = 1)))
#   )
#   assign(paste0("g_barchart_", namevar), g2, envir = .GlobalEnv)
#   # ggsave(paste0(projectpath, "/", "usvi_barchart_", namevar, "-", Sys.Date(), ".png"),
#   #        g2,
#   #        width = 14, height = 8, units = "in")
# }

}


# How many ASVs in these 52 dominant genera? -----------------------------

temp_df2 <- usvi_prok_asvs.taxa %>%
  dplyr::select(asv_id, Domain:Species) %>%
  # dplyr::filter(if_any(c("Domain":"Species"), ~!is.na(.x))) %>%
  dplyr::filter(if_any(c("Domain":"Genus"), ~!is.na(.x))) %>%
  droplevels %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(if_any(c("Domain":"Genus"), ~.x %in% stringr::str_split_i(taxonomy_rank_genus_idx, ";", -1))) %>%
  # dplyr::filter(if_any(c("Family":"Genus"), ~.x %in% stringr::str_split_i(taxonomy_rank_genus_idx, ";", -1))) %>%
  droplevels



# temp_idx <- usvi_filtered_taxonomy_rank_genus_relabund.df %>%
#   dplyr::group_by(taxonomy, asv_id) %>%
#   dplyr::summarise(weight = mean(rel_abund, na.rm = TRUE)) %>%
#   dplyr::mutate(weight = nafill(weight, fill = 0)) %>%
#   dplyr::group_by(taxonomy) %>%
#   dplyr::mutate(weight = relabund(weight)/100) %>%
#   # dplyr::mutate(weight = scale(weight, center = TRUE, scale = FALSE)) %>%
#   # dplyr::summarise(sum_weight = sum(weight, na.rm = TRUE)) %>%
#   # droplevels
#   dplyr::ungroup(.) %>% dplyr::select(asv_id, weight) %>%  tibble::deframe(.)
# 
# temp_df <- usvi_filtered_taxonomy_rank_genus_relabund.df %>%
#   split(., f = .$sample_id) %>%
#   map(., ~.x %>%
#         dplyr::left_join(., tibble::enframe(temp_idx, name = "asv_id", value = "weight"), by = join_by(asv_id)) %>%
#         dplyr::mutate(w_mean = (rel_abund * weight)) %>%
#         # dplyr::group_by(taxonomy) %>%
#         # dplyr::summarise(w_mean = mean(w_mean, na.rm= TRUE)) %>%
#         droplevels) %>%
#   bind_rows(.) %>%
#   dplyr::group_by(taxonomy) %>%
#   dplyr::summarise(w_mean = mean(w_mean, na.rm= TRUE)) %>%
#   droplevels

#what is the average relative abundance of the ASVs assigned to these genera-level taxonomies
#if we look at all ASVs (~2058)
#compared to only those ASVs 1% or more in a sample
temp_df <- usvi_filtered_taxonomy_rank_genus_relabund.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::group_by(taxonomy) %>%
                     dplyr::summarise(`Mean relative abundance of all ASVs` = mean(rel_abund, na.rm= TRUE)) %>%
                     droplevels %>%
  dplyr::left_join(., usvi_filtered_taxonomy_rank_genus_relabund.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::filter(asv_id %in% temp_idx) %>%
                     # dplyr::mutate(keep = dplyr::case_when(rel_abund >= 1 ~ "Abundant",
                     #                                       .default = NA)) %>%
                     tidyr::drop_na(.) %>%
                     dplyr::group_by(taxonomy) %>%
                     dplyr::summarise(`Mean relative abundance of abundant ASVs` = mean(rel_abund, na.rm= TRUE)) %>%
                     droplevels, 
                   by= join_by(taxonomy)) %>%
  tidyr::drop_na(taxonomy) %>%
  dplyr::mutate(across(contains("relative abundance"), nafill, fill = 0)) %>%
  droplevels

#forexample, the average relative abundance of the 4 abundant Synechococcus ASVs was 5.69%
#in contrast, the average relative abundance of all 26 Synechoccocus ASVs was 2.95%
usvi_filtered_taxonomy_rank_genus_relabund.df %>%
  tidyr::drop_na(.) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(total_relabund = sum(rel_abund,na.rm = TRUE)) %>%
  dplyr::left_join(., usvi_filtered_taxonomy_rank_genus_relabund.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::filter(taxonomy %in% names(temp_idx)) %>%
                     droplevels %>%
                     dplyr::group_by(sample_id) %>%
                     dplyr::summarise(most_relabund = sum(rel_abund,na.rm = TRUE)),
                   by = join_by(sample_id)) %>%
  dplyr::left_join(., usvi_filtered_taxonomy_rank_genus_relabund.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::filter(asv_id %in% temp_idx) %>%
                     droplevels %>%
                     # dplyr::filter(rel_abund >= 1) %>%
                     # tidyr::drop_na(.) %>%
                     dplyr::group_by(sample_id) %>%
                     dplyr::summarise(asv_one_relabund = sum(rel_abund,na.rm = TRUE)),
                   by = join_by(sample_id)) %>%
  # dplyr::arrange(most_relabund)
  dplyr::arrange(asv_one_relabund)
  # dplyr::arrange(total_relabund)

usvi_filtered_taxonomy_rank_genus_summary.tbl <- usvi_filtered_taxonomy_rank_genus_relabund.df %>%
  dplyr::ungroup(.) %>%
  # dplyr::filter(rel_abund >= 1) %>%
  dplyr::distinct(asv_id, .keep_all = TRUE) %>%
  dplyr::group_by(taxonomy, Genus) %>%
  dplyr::summarise(`Total ASVs` = length(asv_id)) %>%
  tidyr::drop_na(.) %>%
  dplyr::left_join(., usvi_filtered_taxonomy_rank_genus_relabund.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::filter(asv_id %in% temp_idx) %>%
                     # dplyr::filter(taxonomy %in% names(temp_idx)) %>%
                     dplyr::distinct(asv_id, .keep_all = TRUE) %>%
                     dplyr::group_by(taxonomy, Genus) %>%
                     dplyr::summarise(`Abundant ASVs` = length(asv_id)) %>%
                     tidyr::drop_na(.),
                   by = join_by(taxonomy, Genus)) %>%
  dplyr::left_join(., (usvi_filtered_taxonomy_rank_genus_relabund.df %>% #what is the mean relative abundance of the most abundant ASV in each taxon?
      dplyr::group_by(asv_id, taxonomy) %>%
      dplyr::summarise(`Mean relative abundance of top ASV` = mean(rel_abund)) %>%
      dplyr::arrange(desc(`Mean relative abundance of top ASV`)) %>% 
      dplyr::ungroup(.) %>%
        dplyr::select(taxonomy, `Mean relative abundance of top ASV`) %>%
      dplyr::distinct(taxonomy, .keep_all = TRUE)),
      by = join_by(taxonomy)) %>%
  dplyr::mutate(`Abundant ASVs` = nafill(`Abundant ASVs`, fill = 0)) %>%
  dplyr::left_join(., temp_df, by = join_by(taxonomy)) %>%
  dplyr::mutate(Class = stringr::str_extract(taxonomy, paste0(taxonomy_rank_class_idx, collapse = "|"))) %>%
  dplyr::relocate(Class, .before = Genus) %>%
  dplyr::ungroup(.) %>%
  # dplyr::arrange(desc(`Mean relative abundance of abundant ASVs`)) %>%
  dplyr::arrange(desc(`Mean relative abundance of top ASV`)) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(across(contains("relative abundance"), scales::percent, scale = 1, accuracy = 0.01)) %>%
  droplevels

readr::write_delim(usvi_filtered_taxonomy_rank_genus_summary.tbl, paste0(projectpath, "/", "usvi_filtered_genus_summary", ".tsv"),
                   delim = "\t", col_names = TRUE)

usvi_filtered_taxonomy_rank_genus_summary.df <- usvi_filtered_taxonomy_rank_genus_summary.tbl %>%
  dplyr::select(-contains("relative abundance")) %>%
  tidyr::pivot_longer(., cols = contains("ASVs"),
                      names_to = "parameter",
                      values_to = "count") %>%
  dplyr::mutate(count = nafill(count, fill = 0)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = taxonomy_rank_genus_idx)) %>%
  dplyr::arrange(taxonomy) %>%
  dplyr::mutate(Genus = factor(Genus))

print(
  ggplot(data = usvi_filtered_taxonomy_rank_genus_summary.df,
         aes(x = Genus, y = count, fill = taxonomy, group_by = taxonomy))
  + geom_bar(width = 0.90, show.legend = TRUE, position = "stack", stat = "identity")
  + geom_bar(color = "black", width = 0.90, show.legend = FALSE, position = "stack", stat = "identity")
  + theme_bw()
  + scale_y_continuous(expand = expansion(mult = c(0,0.1)), name = "Number of ASVs", transform = "pseudo_log")
  + scale_fill_manual(values = taxonomy_filt_genus_colors, 
                      # breaks = names(taxonomy_filt_genus_colors),
                      breaks = taxonomy_filt_genus_colors_lookup,
                      labels = names(taxonomy_filt_genus_colors_lookup),
                      drop = TRUE)
  + facet_grid(scales = "free", space = "fixed",
               # labeller = labeller(site = site_lookup),
               parameter ~ .)
  )


# Make fasta file of sequences to make phylogenetic tree ------------------


#in R:

usvi_seqs_idx <- ps_usvi %>%
  phyloseq::subset_samples(., sample_type == "seawater") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>% # remove ASVs not present in any samples
  phyloseq::otu_table(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  dplyr::select(asv_id) %>%
  dplyr::right_join(., usvi_prok_asvs.taxa, by = join_by(asv_id)) %>%
# usvi_seqs_idx <- usvi_prok_asvs.taxa %>%
#   dplyr::filter(asv_id %in% usvi_seqs_idx ) %>%
  droplevels %>%
  dplyr::select(asv_id, sequence) %>%
  dplyr::mutate(abundance = 1)
  # dplyr::select(asv_id, sequence) %>%
  # simplify
  # tibble::deframe(.)
library(dada2)
dada2::uniquesToFasta(usvi_seqs_idx, paste0(projectpath, "/", "usvi_prok_asvs.fna"), ids = usvi_seqs_idx[["asv_id"]], mode = "w")

