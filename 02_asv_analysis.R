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
  transplant_lookup <- c(Initial = "Initial",
                         CRock_CRock = "Channel Rock to Channel Rock",
                         CRock_CBay = "Channel Rock to Cane Bay",
                         CBay_CBay = "Cane Bay to Cane Bay",
                         CBay_CRock = "Cane Bay to Channel Rock")
  coral_lookup <- c(OANN = "Orbicella annularis",
                    OFAV = "Orbicella faveolata",
                    PSTR = "Pseudodiploria strigosa",
                    ACER = "Acorpora cervicornis",
                    PPOR = "Porites porites",
                    seawater = "Seawater")
  colony_lookup <- coral_lookup %>%
    as_tibble_col(column_name = "scientific_coral") %>%
    dplyr::mutate(coral_species = names(coral_lookup))
  colony_lookup <- sample_metadata %>%
    dplyr::select(Colony_origin_ID, coral_species) %>%
    tidyr::drop_na(Colony_origin_ID) %>%
    droplevels %>%
    distinct(., .keep_all = FALSE) %>%
    dplyr::left_join(., colony_lookup) %>%
    dplyr::mutate(across(everything(), as.character))
  colony_lookup <- colony_lookup[["scientific_coral"]] %>%
    setNames(., colony_lookup[["Colony_origin_ID"]])
  
  if(exists("facet_form", inherits = TRUE)){
    facet_form <- deparse(substitute(facet_form))
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
  g <- g + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))
  g <- g + scale_fill_viridis(discrete = TRUE, option = "turbo")
  g <- g + guides(fill = guide_legend(direction = "vertical", ncol = 1))
  g <- g + facet_wrap(scales = "free", ncol = 3,
                      eval(facet_form),
                      labeller = labeller(coral_species = coral_lookup))
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

# Replace NA with last known taxonomy -------------------------------------

usvi_prok_filled.taxa.df <- usvi_prok_asvs.taxa %>%
  dplyr::mutate(Phylum = coalesce(Phylum, Kingdom)) %>%
  dplyr::mutate(Class = coalesce(Class, Phylum)) %>%
  dplyr::mutate(Order = coalesce(Order, Class)) %>%
  dplyr::mutate(Family = coalesce(Family, Order)) %>%
  dplyr::mutate(Genus = coalesce(Genus, Family)) %>%
  dplyr::mutate(Species = coalesce(Species, Genus)) %>%
  dplyr::relocate(asv_id) %>%
  droplevels


# Make a phyloseq object --------------------------------------------------

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
if(!file.exists(paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"))){
  readr::write_rds(ps_usvi, paste0(projectpath, "/", "usvi_prok_phyloseq", ".rds"), compress = "gz")
}


# Broadly summarize at phylum, class, order levels ------------------------

usvi_phylum.tbl <- phyloseq::otu_table(ps_usvi) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  apply(., 2, relabund) %>%
  as.data.frame %>%
  dplyr::mutate(TotAbund = rowSums(.)) %>%
  dplyr::arrange(desc(TotAbund)) %>%
  dplyr::select(-TotAbund) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  otu_to_taxonomy(dataset = .,
                  taxonomy = usvi_prok_filled.taxa.df,
                  level = "Phylum") %>%
  tidyr::pivot_longer(.,
                      cols = !c(1:Phylum),
                      names_to = "sample_id",
                      values_to = "relabund") %>%
  droplevels

usvi_class.tbl <- phyloseq::otu_table(ps_usvi) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  apply(., 2, relabund) %>%
  as.data.frame %>%
  dplyr::mutate(TotAbund = rowSums(.)) %>%
  dplyr::arrange(desc(TotAbund)) %>%
  dplyr::select(-TotAbund) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  otu_to_taxonomy(dataset = .,
                  taxonomy = usvi_prok_filled.taxa.df,
                  level = "Class") %>%
  tidyr::pivot_longer(.,
                      cols = !c(1:Class),
                      names_to = "sample_id",
                      values_to = "relabund") %>%
  droplevels

usvi_order.tbl <- phyloseq::otu_table(ps_usvi) %>%
  as.data.frame %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  apply(., 2, relabund) %>%
  as.data.frame %>%
  dplyr::mutate(TotAbund = rowSums(.)) %>%
  dplyr::arrange(desc(TotAbund)) %>%
  dplyr::select(-TotAbund) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  otu_to_taxonomy(dataset = .,
                  taxonomy = usvi_prok_filled.taxa.df,
                  level = "Order") %>%
  tidyr::pivot_longer(.,
                      cols = !c(1:Order),
                      names_to = "sample_id",
                      values_to = "relabund") %>%
  droplevels


#drop all < 5% order
temp_df <- usvi_order.tbl %>%
  tidyr::pivot_wider(., id_cols = c(Kingdom, Phylum, Class, Order),
                     names_from = "sample_id",
                     values_from = "relabund") %>%
  dplyr::mutate(taxonomy = paste(Kingdom, Phylum, Class, Order, sep = ";")) %>%
  dplyr::select(-c(Kingdom, Phylum, Class, Order)) %>%
  tibble::column_to_rownames(var = "taxonomy") %>%
  # droplevels #354 distinct orders:class:phyla
  dplyr::slice(which(rowSums(.) > 5)) %>%
  t() #49 orders:class:phyla that are over 5% relative abundance over entire dataset of 109 samples
temp_df <- temp_df %>%
  bind_cols("sample_id" = rownames(temp_df),
          ., 
          "Other taxa" = (100 - rowSums(temp_df))) %>%
  tidyr::pivot_longer(., cols = !c("sample_id"),
                      names_to = "taxonomy",
                      values_to = "relabund") %>%
  dplyr::arrange(desc(relabund)) %>%
  dplyr::mutate(across(c(sample_id, taxonomy), factor)) %>%
  dplyr::mutate(taxonomy = factor(taxonomy, levels = unique(.$taxonomy))) %>%
  # dplyr::left_join(., sample_metadata %>%
  #                    dplyr::select(sample_source, sample_id, sample_type, reef_origin, reef_outplant, coral_species, Colony_origin_ID, replant, replicate) %>%
  #                    droplevels) %>%
  droplevels

temp_df %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(totalabund = sum(relabund, na.rm = TRUE))






# broad Taxonomy summary --------------------------------------------------

ASV_env_phylum_summary.tbl <- ASV_env_phylum.df %>%
  dplyr::group_by(taxonomy, sample_type, coral_species) %>%
  dplyr::summarise(across(contains("relabund"), list(mean = mean, min = min, max = max))) %>%
  tidyr::pivot_longer(., cols = contains("relabund"),
                      names_to = "metric",
                      values_to = "relabund") %>%
  dplyr::mutate(metric = gsub("relabund_", "", metric))
total_phylum <- length(unique(ASV_env_phylum.df[["taxonomy"]]))
phylum_colors <- viridisLite::turbo(n = total_phylum) %>%
  # setNames(., as.character(unique(ASV_env_phylum.df$taxonomy))[c(seq(1, total_phylum, by = 2), seq(2, total_phylum, by = 2))])
  setNames(., as.character(unique(ASV_env_phylum.df$taxonomy)))

df <- ASV_env_phylum_summary.tbl %>%
  dplyr::filter(grepl("coral", sample_type)) %>%
  droplevels

{
  g1 <- P_summary_taxonomy_bar(dataset = df %>%
                                 dplyr::filter(grepl("OANN", coral_species)),
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "metric~.")
  
  g2 <- P_summary_taxonomy_bar(dataset = df %>%
                                 dplyr::filter(grepl("OFAV", coral_species)),
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "metric~.")
  
  g3 <- P_summary_taxonomy_bar(dataset = df %>%
                                 dplyr::filter(grepl("PPOR", coral_species)),
                               x = taxonomy,
                               y = relabund, 
                               z = taxonomy,
                               facet_form = "metric~.")
  
  gpatch <- (g1 + ggtitle("Orbicella annularis") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")) / (g2 + ggtitle("Orbicella faveolata") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")) / (g3 + ggtitle("Porites porites") + theme(axis.text.x = element_blank(), legend.direction = "horizontal")) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right")
  gpatch <- gpatch & scale_y_continuous(trans = "pseudo_log", expand = expansion(mult = c(0,0.1)))
  # gpatch <- gpatch & scale_fill_manual(values = phylum_colors, drop = TRUE, guide = "legend")
  gpatch <- gpatch + patchwork::plot_annotation(title = "Summary of relative abundances of phyla by coral species", 
                                                tag_level = "A")
  ggsave(paste0(getwd(), "/", Sys.Date(), "-coral_phylum_abundance_stats.png"),
         gpatch,
         width = 10, height = 8, units = "in")
  # gpatch & scale_fill_manual(values = phylum_colors, drop = TRUE, guide = "legend")
}
