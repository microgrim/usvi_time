## 00_custom_functions.R

#custom functions for your R analyses


# Generate a random number using provided observation and quantiles -------


F_random_sampler <- function(n, x, prob){
  #n is how many observations you want
  #x is your vector (or range) of potential numbers
  #prob is your quantile desired the range of observations
  
  stopifnot(is.numeric(x), is.numeric(prob))
  
  #if you provided a vector from quantile:
  mod_melt <- NULL
  if(length(prob) > 1){
    for(i in seq_len(length(prob))){
      a1 <- min(prob[[i]], x, na.rm = TRUE) #make this your minimum
      b1 <- min(prob[[i]], max(x, na.rm = TRUE), na.rm = TRUE) #Make this the maximum below the quantile
      # b1 <- prob[[i]] #Make this the maximum below the quantile
      c1 <- seq(from = a1, to = b1, length.out = 100)
      mod_melt[[i]] <- c1
    }
    mod_melt <- do.call(rbind, mod_melt)
  } else {
    #if you provided one quantile value < 1, and want just a random assortment of values around that:
    p1 <- min(prob*0.9, prob)
    p2 <- min(prob*1.1, 1)
    p1 <- quantile(x, probs = p1, na.rm = TRUE)
    p2 <- quantile(x, probs = p2, na.rm = TRUE)
    c1 <- seq(from = p1, to = p2, length.out = 100)
    mod_melt <- c1
  }
  sample(mod_melt, n, replace = TRUE)
}


# Relabund ----------------------------------------------------------------


relabund <- function(sample) { #Write a relative abundance function
  x = sample/sum(sample)
  x = x*100
  return(x)
}


# Make a data table from OTUs and common taxonomies -----------------------

otu_to_taxonomy <- function(dataset, taxonomy, level){
  #dataset: col1 is asv_id, cols 2-end are the samples
  #taxonomy: col1 is asv_id, cols 2-end are taxonomic levels
  #level: text string of the specific taxonomic level you want to examine
  
  tax_levels <- c("Domain", "Kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_levels <- grep(paste0(tax_levels, collapse = "|"), colnames(taxonomy), ignore.case = TRUE, value = TRUE)
  tax_levels <- tax_levels[1:grep( {level}, tax_levels, ignore.case = TRUE)]
  
  sum_dataset <- taxonomy %>%
    dplyr::select(1, all_of(tax_levels)) %>%
    dplyr::left_join(., dataset, by = "asv_id") %>%
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



# Transform a counts matrix -----------------------------------------------

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


# Gimme that legend -------------------------------------------------------

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Plot taxonomy bar charts ------------------------------------------------


color_resample_list <- c("#e6d1eb",
                          "#77bd86",
                          "#d3adf3",
                          "#b6fcc0",
                          "#ef9ccd",
                          "#86f6d7",
                          "#ffaad8",
                          "#7ce3bd",
                          "#f79695",
                          "#7efeec",
                          "#ffa3a8",
                          "#5aedef",
                          "#e19b76",
                          "#1abee8",
                          "#f7ba82",
                          "#4ed3ff",
                          "#d4a266",
                          "#7eb5fa",
                          "#fff1a3",
                          "#a1a9f2",
                          "#e2f1a4",
                          "#bfb6ff",
                          "#b1c97f",
                          "#c29ddd",
                          "#acdc98",
                          "#d6beff",
                          "#a4c780",
                          "#5eb3ee",
                          "#fad08c",
                          "#4cb6e7",
                          "#e4be7a",
                          "#03bddb",
                          "#ffb29c",
                          "#3adbeb",
                          "#e3998a",
                          "#34d4da",
                          "#ffacbd",
                          "#2ac4be",
                          "#e095b8",
                          "#97ebb8",
                          "#ffbfeb",
                          "#70c697",
                          "#ffb8d5",
                          "#8ad39d",
                          "#c99dcb",
                          "#82b97c",
                          "#f6d5ff",
                          "#b1b267",
                          "#a4a8de",
                          "#d2cb7f",
                          "#8bd1ff",
                          "#c4ac65",
                          "#72f4ff",
                          "#e09997",
                          "#96ffed",
                          "#de98ab",
                          "#94fffb",
                          "#d19bbf",
                          "#d2ffc9",
                          "#bba3c9",
                          "#f1ffc0",
                          "#d1cdff",
                          "#dfcc81",
                          "#97acd3",
                          "#ffe5ae",
                          "#5fb8cb",
                          "#ffd2ac",
                          "#8de5ff",
                          "#cea289",
                          "#94f5ff",
                          "#ffbbc1",
                          "#6ebaa2",
                          "#ffcbe9",
                          "#8cb686",
                          "#ffd6f4",
                          "#a0b37d",
                          "#dddaff",
                          "#f9ffcb",
                          "#a6a9d0",
                          "#e7ffd5",
                          "#c9b4cc",
                          "#cdffde",
                          "#c6a681",
                          "#b8f0ff",
                          "#bda984",
                          "#c9fff5",
                          "#ffcfc4",
                          "#6bb9af",
                          "#ffdbb7",
                          "#88b4ab",
                          "#fcdbee",
                          "#9cb38b",
                          "#ffe6ce",
                          "#8bb59b",
                          "#fff1d2",
                          "#b0ad85",
                          "#dbf7e1",
                          "#bbb095",
                          "#dee2c5",
                          "#d6c6ac")


# Global variables --------------------------------------------------------


set.alpha <- 0.05
set.seed(48105)

model_dispersion_lookup <- data.frame(v1 = c("manual_dispersion", 
                                             "manual_dawn", "manual_photo", "auto_dawn", "auto_photo",
                                             "auto_dispersion"),
                                      label = c("All LB seagrass, manual dispersion",
                                                "Dawn LB seagrass, manual dispersion",
                                                "Peak photo LB seagrass, manual dispersion",
                                                "Dawn LB seagrass, auto dispersion",
                                                "Peak photo LB seagrass, auto dispersion",
                                                "All LB seagrass, auto dispersion")) %>%
  tibble::deframe(.)


group_labels_lookup <- c("Tektite_d" = "Tektite dawn",
                         "Tektite_p" = "Tektite peak photo",
                         "Yawzi_d" = "Yawzi dawn",
                         "Yawzi_p" = "Yawzi peak photo",
                         "LB_seagrass_d" = "LB seagrass dawn",
                         "LB_seagrass_p" = "LB seagrass peak photo",
                         "Tektite" = "Tektite",
                         "Yawzi" = "Yawzi",
                         "LB_seagrass" = "LB seagrass")

group_labels_colors <- viridisLite::turbo(length(group_labels_lookup)) %>%
  setNames(., names(group_labels_lookup))

site_lookup <- data.frame(site = c("LB_seagrass", "Tektite", "Yawzi", "control_extraction", "control_pcr", "control_seq"),
                          label = c("Lameshur Bay seagrass", "Tektite Reef", "Yawzi Reef",
                                    "Control (DNA Extraction)", "Control (PCR)", "Control (Sequencing)")) %>%
  dplyr::mutate(label = stringr::str_wrap(label, width = 16)) %>%
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

# usvi_genera_relabel <- usvi_prok_filled.taxa.df %>%
#   dplyr::mutate(across(everything(), ~stringr::str_replace_all(.x, " clade", ""))) %>%
#   dplyr::rowwise(.) %>%
#   dplyr::mutate(first = dplyr::case_when((Order == Genus) ~ Class,
#                                          .default = Order)) %>%
#   dplyr::mutate(first = dplyr::case_when((Class == Order) ~ NA,
#                                          (Class != Phylum) ~ Phylum,
#                                          grepl(paste0(c("Synechococcales"), collapse = "|"), first) ~ NA,
#                                          .default = first)) %>%
#   dplyr::mutate(second = dplyr::case_when((Order == Genus) ~ Order,
#                                           .default = Genus)) %>%
#   dplyr::mutate(taxa_label = dplyr::case_when(!is.na(first) ~ paste0(first, "; ", second),
#                                               .default = second)) %>%
#   dplyr::select(asv_id, taxa_label) %>%
#   tibble::deframe(.)
# 
global_labeller <- labeller(
  # model2 = model_dispersion_lookup,
  model = model_dispersion_lookup,
  Combo = group_labels_lookup,
  site = site_lookup,
  sampling_day = sampling_day_lookup,
  sampling_time = sampling_time_lookup,
  # asv_id = usvi_genera_relabel,
  .multi_line = TRUE,
  .default = label_both
)
