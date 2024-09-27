# 02_asv_analysis.R

# Load packages -----------------------------------------------------------

# if(!require("ggvegan", quietly = TRUE)){
#   install.packages("remotes")
#   remotes::install_github("gavinsimpson/ggvegan")
# }
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(tidyverse)
library(data.table)
library(BiocManager)
library(BiocParallel)
library(dada2)
library(phyloseq)
# library(ggvegan)
library(viridis)
library(DESeq2)
library(cli)
library(furrr)
library(progressr)
library(patchwork)

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


# Custom functions --------------------------------------------------------

#Turn data into relative abundances (since these data are compositional = quantitative descriptions of the parts of some whole, conveying relative information.)
relabund <- function(sample) { #Write a relative abundance function
  x = sample/sum(sample)
  x = x*100
  return(x)
}


F_trans_abund <- function(counts, mdata, suffix){
  #this function will take the ASV matrix and your subsetted metadata,
  #and do 3 things:
  # filter out the samples retained in the metadata
  # drop the ASVs with 0 occurrences in the subset of samples
  # log2(x + 1) transform the relative abundances of ASVs
  
  if(!missing(suffix)){
    suffix <- rlang::parse_expr(suffix)
  } else {
    suffix <- "filtered"
  }
  namevar <- deparse(substitute(counts))
  metadata <- get0("mdata", inherits = TRUE)
  lookup <- c(sample = "sample_id", sample = "sample_name", sample = "samplename")
  metadata <- metadata %>%
    dplyr::rename(., any_of(lookup))
  
  if(any(colnames(counts) %in% levels(metadata$sample))){
    temp_counts <- counts %>%
      dplyr::select(as.character(metadata[["sample"]])) %>%
      tibble::rownames_to_column(var = "ASV")
  } else {
    temp_counts <- counts %>%
      dplyr::filter(sample %in% levels(metadata$sample)) %>%
      droplevels %>%
      data.table::transpose(., make.names = "sample", keep.names = "ASV")
  }
  temp_counts <- temp_counts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(c_across(!matches("ASV")))) %>%
    dplyr::slice(which(total != 0)) %>% #drop any ASVs that are nonoccurring in the filtered dataset
    dplyr::select(-total) %>%
    tibble::column_to_rownames(., var = "ASV") %>%
    apply(., 2, relabund) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., var = "ASV") %>%
    tidyr::pivot_longer(., cols = -c("ASV"),
                        names_to = "sample",
                        values_to = "abundance") %>%
    dplyr::mutate(logabund = ifelse(!(is.na(abundance) | (abundance < 0)),
                                    log2(abundance+1), #log transform abundance (with +1 pseudocount)
                                    0)) %>%
    tidyr::pivot_wider(., id_cols = "sample",
                       names_from = "ASV",
                       values_from = "logabund")
  assign(paste0(namevar, "_", suffix), temp_counts, envir = .GlobalEnv)
}

otu_to_taxonomy <- function(dataset, taxonomy, level){
  #dataset: col1 is asv_id, cols 2-end are the samples
  #taxonomy: col1 is asv_id, cols 2-end are taxonomic levels
  #level: text string of the specific taxonomic level you want to examine
  
  tax_levels <- c("Domain", "Kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_levels <- grep(paste0(tax_levels, collapse = "|"), colnames(taxonomy), ignore.case = TRUE, value = TRUE)
  tax_levels <- tax_levels[1:grep( {level}, tax_levels, ignore.case = TRUE)]
  
  sum_dataset <- taxonomy %>%
    dplyr::select(1, all_of(tax_levels)) %>%
    left_join(., dataset, by = "asv_id") %>%
    tidyr::pivot_longer(., cols = !c("asv_id", all_of(tax_levels)),
                        names_to = "sample_id",
                        values_to = "abundance") %>%
    dplyr::summarise(relabund = sum(abundance, na.rm = TRUE),
                     .by = c(sample_id, all_of(tax_levels))) %>%
    droplevels %>%
    tidyr::pivot_wider(.,
                       names_from = "sample_id",
                       values_from = "relabund", 
                       values_fill = 0)
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



sample_metadata <- readr::read_delim(paste0(projectpath, "/", "usvi_metadata.txt"),
                                     delim = "\t",
                                     quote = "",
                                     col_names = TRUE,
                                     show_col_types = FALSE,
                                     num_threads = nthreads)
temporal_lookup <- list(sampling_time = c("dawn", "peak_photo", NA),
                        sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5", NA),
                        sampling_date = c(20210122, 20210123, 20210124, 20210125, 20210126, NA),
                        site = c("Tektite", "Yawzi", "LB_seagrass", NA),
                        sample_type = c("seawater", "control_extraction", "control_pcr", "control_seq"))

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
                         dplyr::arrange(site, sampling_time, sampling_day) %>%
                         dplyr::group_by(site) %>%
                         dplyr::mutate(sample_order = seq(1, length(interaction(replicate, sampling_time, sampling_day)))) %>%
                         dplyr::ungroup(.) %>%
                         dplyr::select(sample_ID, sample_order) %>%
                         droplevels),
                         (metadata %>%
                            dplyr::filter(!grepl("seawater", sample_type)) %>%
                            droplevels %>%
                            dplyr::arrange(sample_type, replicate) %>%
                            dplyr::group_by(sample_type) %>%
                            dplyr::mutate(sample_order = seq(1, length(interaction(replicate, sample_type)))) %>%
                            dplyr::ungroup(.) %>%
                            dplyr::select(sample_ID, sample_order) %>%
                            droplevels)),
                   by = join_by(sample_ID))
  

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



# Broadly summarize at phylum, class, order levels ------------------------

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

#make relative abundance tables at specific taxonomic levels
keep <- c("Phylum", "Class", "Order")
# keep <- c("Phylum")

for(taxlevel in keep){
  namevar <- stringr::str_to_lower(taxlevel)
  temp_df <- phyloseq::otu_table(ps_usvi) %>%
    as.data.frame %>%
    dplyr::slice(which(rowSums(.) > 0)) %>%
    apply(., 2, relabund) %>%
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
  
  keep_tax <- colnames(usvi_prok_filled.taxa.df)[c(grep("Domain", colnames(usvi_prok_filled.taxa.df)):grep({{taxlevel}}, colnames(usvi_prok_filled.taxa.df)))]
  
  temp_df2 <- temp_df %>%
    tidyr::pivot_wider(., id_cols = keep_tax,
                       names_from = "sample_id",
                       values_from = "relabund") %>%
    tidyr::unite("taxonomy", c(all_of(keep_tax)), sep = ";", remove = TRUE) %>%
    dplyr::relocate(taxonomy) %>%
    tibble::column_to_rownames(var = "taxonomy") %>%
    dplyr::slice(which(rowSums(.) > 5)) %>%
    t()
    # droplevels
  
  temp_df2 <- temp_df2 %>%
    bind_cols("sample_id" = rownames(temp_df2),
              ., 
              "Other taxa" = (100 - rowSums(temp_df2))) %>%
    tidyr::pivot_longer(., cols = !c("sample_id"),
                        names_to = "taxonomy",
                        values_to = "relabund") %>%
    dplyr::arrange(desc(relabund)) %>%
    dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
    dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
    tidyr::separate_wider_delim(taxonomy, names = keep_tax, delim = ";", too_few = "align_start", too_many = "drop", cols_remove = FALSE) %>%
    droplevels
  
  assign(paste0("usvi_", namevar, ".df"), temp_df, envir = .GlobalEnv)
  assign(paste0("usvi_", namevar, "_over5.df"), temp_df2, envir = .GlobalEnv)
  rm(temp_df)
}


#non looped version:
{
  # usvi_order.tbl <- phyloseq::otu_table(ps_usvi) %>%
  #   as.data.frame %>%
  #   dplyr::slice(which(rowSums(.) > 0)) %>%
  #   apply(., 2, relabund) %>%
  #   as.data.frame %>%
  #   dplyr::mutate(TotAbund = rowSums(.)) %>%
  #   dplyr::arrange(desc(TotAbund)) %>%
  #   dplyr::select(-TotAbund) %>%
  #   tibble::rownames_to_column(var = "asv_id") %>%
  #   otu_to_taxonomy(dataset = .,
  #                   level = "Order",
  #                   taxonomy = usvi_prok_filled.taxa.df) %>%
  #   tidyr::pivot_longer(.,
  #                       cols = !c(1:all_of("Order")),
  #                       names_to = "sample_id",
  #                       values_to = "relabund") %>%
  #   droplevels
  # 
  # temp_df <- usvi_order.tbl %>%
  #   tidyr::pivot_wider(., id_cols = c(Domain, Phylum, Class, Order),
  #                      names_from = "sample_id",
  #                      values_from = "relabund") %>%
  #   dplyr::mutate(taxonomy = paste(Domain, Phylum, Class, Order, sep = ";")) %>%
  #   dplyr::select(-c(Domain, Phylum, Class, Order)) %>%
  #   tibble::column_to_rownames(var = "taxonomy") %>%
  #   # droplevels #354 distinct orders:class:phyla
  #   dplyr::slice(which(rowSums(.) > 5)) %>%
  #   t() #49 orders:class:phyla that are over 5% relative abundance over entire dataset of 109 samples
  # temp_df <- temp_df %>%
  #   bind_cols("sample_id" = rownames(temp_df),
  #             ., 
  #             "Other taxa" = (100 - rowSums(temp_df))) %>%
  #   tidyr::pivot_longer(., cols = !c("sample_id"),
  #                       names_to = "taxonomy",
  #                       values_to = "relabund") %>%
  #   dplyr::arrange(desc(relabund)) %>%
  #   dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  #   dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  #   # dplyr::left_join(., sample_metadata %>%
  #   #                    dplyr::select(sample_source, sample_id, sample_type, reef_origin, reef_outplant, coral_species, Colony_origin_ID, replant, replicate) %>%
  #   #                    droplevels) %>%
  #   droplevels
  # 
  # #sanity check:
  # temp_df %>%
  #   dplyr::group_by(sample_id) %>%
  #   dplyr::summarise(totalabund = sum(relabund, na.rm = TRUE))
  # 
}

length(unique(usvi_phylum_over5.df[["taxonomy"]]))
#17 phyla have cumulative relative abundance over all samples 5% or more
length(unique(usvi_class_over5.df[["taxonomy"]]))
length(unique(usvi_class_over5.df[["Phylum"]]))
#26 classes
#representing 17 phyla

length(unique(usvi_order_over5.df[["taxonomy"]]))
length(unique(usvi_order_over5.df[["Class"]]))
#49 orders
#representing 22 classes

# broad Taxonomy summary --------------------------------------------------



ASV_env_phylum_summary.tbl <- usvi_phylum_over5.df %>%
  dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
  dplyr::group_by(taxonomy, sample_type, site, sampling_time) %>%
  dplyr::summarise(across(contains("relabund"), list(mean = mean, min = min, max = max))) %>%
  tidyr::pivot_longer(., cols = contains("relabund"),
                      names_to = "metric",
                      values_to = "relabund") %>%
  dplyr::mutate(metric = gsub("relabund_", "", metric))

total_phylum <- length(unique(usvi_phylum_over5.df[["taxonomy"]]))
phylum_colors <- viridisLite::turbo(n = total_phylum) %>%
  # setNames(., as.character(unique(ASV_env_phylum.df$taxonomy))[c(seq(1, total_phylum, by = 2), seq(2, total_phylum, by = 2))])
  setNames(., as.character(unique(usvi_phylum_over5.df$taxonomy)))

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
  gpatch <- gpatch & scale_y_continuous(trans = "pseudo_log", expand = expansion(mult = c(0,0.1)))
  # gpatch <- gpatch & scale_fill_manual(values = phylum_colors, drop = TRUE, guide = "legend")
  gpatch <- gpatch + patchwork::plot_annotation(title = "Summary of relative abundances of phyla by site", 
                                                tag_level = "A")
  ggsave(paste0(projectpath, "/", "usvi_phylum_abundance_stats-", Sys.Date(), ".png"),
         gpatch,
         width = 10, height = 16, units = "in")
  # gpatch & scale_fill_manual(values = phylum_colors, drop = TRUE, guide = "legend")
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
          "SAR406", "Verrucomicrobiota", "Thermoplasmatota")

temp_df <- usvi_class_over5.df %>%
  dplyr::left_join(., metadata, by = join_by(sample_id)) %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  # dplyr::filter(Phylum %in% keep) %>%
  dplyr::filter(Phylum %in% grep(paste0(keep, collapse = "|"), unique(usvi_class_over5.df[["Phylum"]]), value = TRUE)) %>%
  dplyr::arrange(desc(relabund), Class, taxonomy) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.[["taxonomy"]]))) %>%
  dplyr::arrange(sampling_day, sampling_time, replicate) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(.[["sample_id"]]))) %>%
  # dplyr::group_by(taxonomy, sample_type, site, sampling_time) %>%
  # dplyr::summarise(across(contains("relabund"), list(mean = mean, min = min, max = max))) %>%
  dplyr::group_by(taxonomy, sample_type, site, sampling_day, sampling_time,replicate, sample_id, sample_order) %>%
  dplyr::summarise(across(contains("relabund"), list(mean = mean))) %>%
  tidyr::pivot_longer(., cols = contains("relabund"),
                      names_to = "metric",
                      values_to = "relabund") %>%
  dplyr::mutate(metric = gsub("relabund_", "", metric))

{
  g <- ggplot(data = temp_df, aes(y = relabund,
                                  x = sample_order,
                                  # x = interaction(replicate, sampling_day),
                                  fill = taxonomy))
  g <- g + geom_col(width = 0.90, show.legend = TRUE, position = "stack")
  g <- g + geom_col(color = "black", width = 0.90, show.legend = FALSE, position = "stack")
  g <- g + theme_bw() 
  g <- g + scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  g <- g + scale_fill_viridis(discrete = TRUE, option = "turbo")
  g <- g + guides(fill = guide_legend(direction = "vertical", ncol = 1))
  g <- g + facet_grid(scales = "free", space = "free",
                      # labeller = labeller(coral_species = coral_lookup),
                      # rows = vars({{facet_rows}}), cols = vars({{facet_cols}}))
                      . ~ site)
  g <- g + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
                 strip.text.x = element_text(angle = 0),
                 strip.text.y = element_text(angle = -90))
  g <- g + labs(x = "sample",
                y = "Relative abundance (%)")
  g <- g + coord_cartesian(ylim = c(0,100), expand = FALSE)
  g <- g + guides(fill = guide_legend(order = 2, ncol = 1, title = "Class", direction = "vertical",
                                      override.aes = list(color = "black", stroke = 1)))
  g
  }

#26 classes
#representing 17 phyla

length(unique(usvi_order_over5.df[["taxonomy"]]))
length(unique(usvi_order_over5.df[["Class"]]))
