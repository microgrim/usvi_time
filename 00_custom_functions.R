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


# Plot taxonomy bar charts ------------------------------------------------


