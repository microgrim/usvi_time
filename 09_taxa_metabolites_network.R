## 09_taxa_metabolites_network.R



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
# library(phyloseq)
if (!require("NetCoMi", quietly = TRUE)){
  install.packages("devtools")
  devtools::install_github("zdk123/SpiecEasi")
  devtools::install_github("GraceYoon/SPRING")
  devtools::install_github("stefpeschel/NetCoMi", 
                           repos = c("https://cloud.r-project.org/",
                                     BiocManager::repositories()))
  NetCoMi::installNetCoMiPacks()
}
library(NetCoMi)
library(igraph)
# library(ggvegan)
# library(ggtext)
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

set.seed(48105)

# find existing processed files -------------------------------------------

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

to_import <- c("usvi_prok_asvs.df", "usvi_prok_asvs.taxa", 
               "metabolomics_sample_metadata")

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
      if(file.exists(paste0(projectpath, "/", namevar, ".tsv"))){
        cli::cli_alert_info("Importing this dataset: {namevar}")
        temp_df <- readr::read_delim(paste0(projectpath, "/", namevar, ".tsv"),
                                     col_names = TRUE, show_col_types = FALSE, delim = "\t", num_threads = nthreads)
        assign(paste0(namevar), temp_df, envir = .GlobalEnv)
        rm(temp_df)
        to_import <- to_import[grep(namevar, to_import, value = FALSE, invert = TRUE)]
      } else {
        cli::cli_alert_warning("Please prepare this dataset: {namevar}")  
      }
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
  dplyr::relocate(asv_id) %>%
  dplyr::mutate(Phylum = dplyr::case_when(grepl("Gammaproteobacteria", Class) ~ "Gammaproteobacteria",
                                          grepl("Alphaproteobacteria", Class) ~ "Alphaproteobacteria",
                                          .default = Phylum)) %>%
  droplevels

if(file.exists(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))){
  temp_list <- readr::read_rds(paste0(projectpath, "/", "usvi_metabolomics_dfs_list", ".rds"))
  list2env(temp_list, envir = .GlobalEnv)
  rm(temp_list)
} else {
  cli::cli_alert_warning("Please tidy the metabolomics datasets.")
}

sample_relabel <- metabolomics_sample_metadata %>%
  dplyr::select(sample_id, site, sampling_day, sampling_time) %>%
  dplyr::distinct(., .keep_all = TRUE) %>%
  dplyr::arrange(site, sampling_time, sampling_day) %>%
  droplevels %>%
  tidyr::unite("relabeled_sample", c(site, sampling_day, sampling_time), sep = "_", remove = TRUE)  %>%
  tibble::deframe(.)


# Tidy datasets for analyses ----------------------------------------------

#try first at the genus level

usvi_sw_genus.taxa.df <- usvi_prok_filled.taxa.df %>%
  dplyr::select(asv_id, Domain, Phylum, Class, Order, Family, Genus) %>%
  dplyr::distinct(Domain, Phylum, Class, Order, Family, Genus, .keep_all = TRUE) %>%
  droplevels

usvi_sw_genus.tbl <- usvi_prok_asvs.df %>%
  dplyr::right_join(., metadata %>%
                      dplyr::filter(grepl("seawater", sample_type)) %>%
                      dplyr::distinct(sample_ID) %>%
                      droplevels,
                    by = join_by("sample_ID"), relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "asv_id",
                     names_from = "sample_ID",
                     values_from = "counts",
                     values_fill = 0) %>%
  otu_to_taxonomy(., usvi_prok_filled.taxa.df, level = "Genus") %>%
  dplyr::left_join(., usvi_sw_genus.taxa.df,
                   by = join_by(Domain, Phylum, Class, Order, Family, Genus)) %>%
  dplyr::relocate(asv_id) %>%
  dplyr::select(-c(Domain, Phylum, Class, Order, Family, Genus)) %>%
  # tibble::column_to_rownames(var = "asv_id") %>%
  droplevels


#here are metabolites where we don't have LODs reported:
usvi_sus_metabolites_idx <- usvi_metabolomics_long.df %>%
  dplyr::arrange(LOD) %>%
  dplyr::distinct(metabolites, LOD, LOQ, .keep_all = TRUE) %>%
  dplyr::arrange(metabolites) %>%
  dplyr::filter(is.na(LOD)) %>%
  droplevels

drop <- c("CINAR_BC_73")
usvi_metab.tbl <- usvi_metabolomics.df %>%
  dplyr::filter(!grepl(paste0(drop, collapse = "|"), metab_deriv_label)) %>%
  dplyr::left_join(., (metabolomics_sample_metadata %>%
                         dplyr::filter(grepl("seawater", sample_type)) %>%
                         dplyr::select(sample_id, metab_deriv_label) %>%
                         droplevels),
                   by = join_by(metab_deriv_label), multiple = "all", relationship = "many-to-many") %>%
  dplyr::relocate(sample_id) %>%
  tidyr::pivot_longer(., cols = !c(sample_id, metab_deriv_label),
                      names_to = "simpleName",
                      values_to = "conc") %>%
  dplyr::filter(!(simpleName %in% usvi_sus_metabolites_idx[["metabolites"]])) %>%
  droplevels %>%
  dplyr::mutate(across(c(sample_id, metab_deriv_label, simpleName), ~factor(.x))) %>%
  dplyr::ungroup(.) %>%
  dplyr::group_by(sample_id, simpleName) %>%
  # dplyr::summarise(num = length(conc)) %>%
  dplyr::summarise(mean_conc = mean(conc, na.rm = TRUE), 
                   num = length(conc),
                   # .by = c(sample_id, simpleName),
                   .groups = "keep",
                   sd = sd(conc, na.rm = TRUE)) %>%
  dplyr::rename(conc = "mean_conc") %>%
  dplyr::mutate(log_conc = ifelse(!(is.na(conc) | (conc < 0)),
                                  log2(conc+1), #log transform abundance (with +1 pseudocount)
                                  0)) %>%
  # dplyr::select(sample_id, simpleName, log_conc) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "sample_id",
                     # values_from = "log_conc",
                     values_from = "conc",
                     names_from = "simpleName") %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  droplevels

usvi_top100_genus.tbl <- usvi_sw_genus.tbl %>%
  tibble::column_to_rownames(var = "asv_id") %>%
  apply(., 2, relabund) %>%
  as.data.frame(.) %>%
  dplyr::slice(which(rowSums(.) > 0)) %>%
  dplyr::mutate(TotAbund = rowSums(.)) %>%
  dplyr::arrange(desc(TotAbund)) %>%
  dplyr::slice_head(., n = 100) %>%
  dplyr::select(-TotAbund) %>%
  tibble::rownames_to_column(var = "asv_id") %>%
  tidyr::pivot_longer(., cols = -c("asv_id"),
                      names_to = "sample",
                      values_to = "abundance")  %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = "sample",
                     values_from = "abundance",
                     names_from = "asv_id") %>%
  tibble::column_to_rownames(var = "sample") %>%
  tidyr::drop_na(.) %>%
  droplevels


usvi_top100_genus_with_metab.tbl <- usvi_top100_genus.tbl[rownames(usvi_metab.tbl),]

usvi_metab.groups <- metabolomics_sample_metadata %>%
  dplyr::distinct(sample_id, sampling_time) %>%
  dplyr::filter(sample_id %in% rownames(usvi_top100_genus_with_metab.tbl)) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = rownames(usvi_top100_genus_with_metab.tbl))) %>%
  droplevels %>%
  dplyr::mutate(sampling_time = dplyr::case_when(sampling_time == "dawn" ~ 1,
                                                 sampling_time == "peak_photo" ~ 2,
                                                 .default = NA)) %>%
  dplyr::select(sampling_time) %>%
  tibble::deframe(.) %>%
  as.integer(.)


# spiec-easi and other methods for asv-asv analyses -----------------------



#make a correlation matrix for ASVs

## usvi_top100_genus.cor <- cor(usvi_top100_genus_with_metab.tbl, method = "spearman")
#try using spieceasi or Spring for microbial correlations

# #this is sample to sample:
# if(!file.exists(paste0(projectpath, "/", "se.usvi_top100_genus_t.gl", ".rds"))){
#   se.usvi_top100_genus_t.gl <- spiec.easi(t(as.matrix(usvi_top100_genus_with_metab.tbl)), method = "glasso", 
#                                         pulsar.select = TRUE, 
#                                         pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, ncores = nthreads))
#   readr::write_rds(se.usvi_top100_genus_t.gl, paste0(projectpath, "/", "se.usvi_top100_genus_t.gl", ".rds"), compress = "gz")
# } else {
#   se.usvi_top100_genus_t.gl <- readr::read_rds(paste0(projectpath, "/", "se.usvi_top100_genus_t.gl", ".rds"))
# }

#this is genera to genera:


if(!file.exists( paste0(projectpath, "/", "se.usvi_top100_genus.mb", ".rds"))){
  se.usvi_top100_genus.mb <- spiec.easi(as.matrix(usvi_top100_genus_with_metab.tbl), 
                                        method = "mb", 
                                        pulsar.select = TRUE, 
                                        pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 50, seed = 48105, ncores = nthreads))
  
  readr::write_rds(se.usvi_top100_genus.mb, paste0(projectpath, "/", "se.usvi_top100_genus.mb", ".rds"), compress = "gz")
  
} else {
  se.usvi_top100_genus.mb <- readr::read_rds(paste0(projectpath, "/", "se.usvi_top100_genus.mb", ".rds"))
}

#it looks like "glasso" method with spiec-easi isn't as good as "mb"

if(!file.exists( paste0(projectpath, "/", "se.usvi_top100_genus.gl", ".rds"))){
  se.usvi_top100_genus.gl <- spiec.easi(as.matrix(usvi_top100_genus_with_metab.tbl), method = "glasso", 
                                        pulsar.select = TRUE, 
                                        lambda.min.ratio = 0.01, 
                                        nlambda = 20, 
                                        pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, ncores = nthreads))
  readr::write_rds(se.usvi_top100_genus.gl, paste0(projectpath, "/", "se.usvi_top100_genus.gl", ".rds"), compress = "gz")
} else {
  se.usvi_top100_genus.gl <- readr::read_rds(paste0(projectpath, "/", "se.usvi_top100_genus.gl", ".rds"))
}


#can we do slr method? "sparse + low-rank decomposition"
if(!file.exists( paste0(projectpath, "/", "se.usvi_top100_genus.slr", ".rds"))){
  se.usvi_top100_genus.slr <- spiec.easi(as.matrix(usvi_top100_genus_with_metab.tbl), method = "slr", 
                                         pulsar.select = TRUE, 
                                         lambda.min.ratio = 0.01, 
                                         nlambda = 100,
                                         pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 50, seed = 48105, ncores = nthreads))
  readr::write_rds(se.usvi_top100_genus.slr, paste0(projectpath, "/", "se.usvi_top100_genus.slr", ".rds"), compress = "gz")
} else {
  se.usvi_top100_genus.slr <- readr::read_rds(paste0(projectpath, "/", "se.usvi_top100_genus.slr", ".rds"))
}


#examine ROC via huge and via stars:
{
  # #make graph topology matrix:
  # # usvi_genus_graph <- SpiecEasi::make_graph('cluster', 
  # #                                           D = ncol(usvi_top100_genus_with_metab.tbl), 
  # #                                           # e = nrow(usvi_top100_genus_with_metab.tbl))
  # #                                           e = ncol(usvi_top100_genus_with_metab.tbl))
  # usvi_genus_graph <- igraph::sample_k_regular(ncol(usvi_top100_genus_with_metab.tbl),
  #                                              ncol(usvi_top100_genus_with_metab.tbl),
  #                                              # nrow(usvi_top100_genus_with_metab.tbl),
  #                                              directed = FALSE, multiple = TRUE)
  # 
  # # pulsar(data = X, fun = match.fun(estFun), fargs = args, criterion = "stars", 
  # #        thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, 
  # #        ncores = c(system = 10))
  # stars.roc(getOptMerge(se.usvi_top100_genus.gl), usvi_genus_graph, verbose=FALSE)
  # # True Postive Rate: from 0 to 0.1587251 
  # # False Positive Rate: from 0 to 0.159461 
  # # Area under Curve: 0.01323045 
  # # Maximum F1 Score: 0.2408236 
  # stars.pr(getOptMerge(se.usvi_top100_genus.gl), usvi_genus_graph, verbose=FALSE)
  # # True Postive Rate: from 0 to 0.1587251 
  # # False Positive Rate: from 0 to 0.159461 
  # # Area under Curve: 0.2980994 
  # # Maximum F1 Score: 0.2408236 
  # huge::huge.roc(se.usvi_top100_genus.gl$est$path, usvi_genus_graph, verbose = FALSE)
  # # True Postive Rate: from 0 to 0.5258757 
  # # False Positive Rate: from 0 to 0.5421112 
  # # Area under Curve: 0.1455523 
  # # Maximum F1 Score: 0.5085871 
  # 
  # # usvi_genus_graph <- SpiecEasi::make_graph('cluster', D = ncol(usvi_top100_genus_with_metab.tbl), e = ncol(usvi_top100_genus_with_metab.tbl))
  # # usvi_genus_graph <- SpiecEasi::make_graph('scale_free',  D = ncol(usvi_top100_genus_with_metab.tbl), e = ncol(usvi_top100_genus_with_metab.tbl))
  # #                                           # e = nrow(usvi_top100_genus_with_metab.tbl))
  # # pulsar(data = X, fun = match.fun(estFun), fargs = args, criterion = "stars", 
  # #        thresh = 0.05, subsample.ratio = 0.8, rep.num = 50, seed = 48105, 
  # #        ncores = c(system = 10))
  # stars.roc(getOptMerge(se.usvi_top100_genus.mb), usvi_genus_graph, verbose=FALSE, plot = TRUE)
  # # True Postive Rate: from 0 to 0.1110406 
  # # False Positive Rate: from 0 to 0.121802 
  # # Area under Curve: 0.006192195 
  # # Maximum F1 Score: 0.1801375 
  # stars.pr(getOptMerge(se.usvi_top100_genus.mb), usvi_genus_graph, verbose=FALSE, plot = TRUE)
  # # True Postive Rate: from 0 to 0.1110406 
  # # False Positive Rate: from 0 to 0.121802 
  # # Area under Curve: 0.2521146 
  # # Maximum F1 Score: 0.1801375 
  # huge::huge.roc(se.usvi_top100_genus.mb$est$path, usvi_genus_graph, verbose = FALSE)
  # # True Postive Rate: from 0 to 0.7737944 
  # # False Positive Rate: from 0 to 0.7775306 
  # # Area under Curve: 0.2968418 
  # # Maximum F1 Score: 0.6065824 
  # 
  #
  # usvi_genus_graph <- igraph::sample_k_regular(ncol(usvi_top100_genus_with_metab.tbl),
  #                                              ncol(usvi_top100_genus_with_metab.tbl),
  #                                              # nrow(usvi_top100_genus_with_metab.tbl),
  #                                              directed = FALSE, multiple = TRUE)
  # # pulsar(data = X, fun = match.fun(estFun), fargs = args, criterion = "stars", 
  # #        thresh = 0.05, subsample.ratio = 0.8, rep.num = 50, seed = 48105, 
  # #        ncores = c(system = 10))
  # stars.roc(getOptMerge(se.usvi_top100_genus.slr), usvi_genus_graph, verbose=FALSE)
  # # True Postive Rate: from 0 to 0.1189766 
  # # False Positive Rate: from 0 to 0.1263814 
  # # Area under Curve: 0.00729227 
  # # Maximum F1 Score: 0.1910721 
  # stars.pr(getOptMerge(se.usvi_top100_genus.slr), usvi_genus_graph, verbose=FALSE)
  # # True Postive Rate: from 0 to 0.1189766 
  # # False Positive Rate: from 0 to 0.1263814 
  # # Area under Curve: 0.2463619 
  # # Maximum F1 Score: 0.1910721 
  # huge::huge.roc(se.usvi_top100_genus.slr$est$path, usvi_genus_graph, verbose = FALSE)
  # # True Postive Rate: from 0 to 0.1148956 
  # # False Positive Rate: from 0 to 0.1241145 
  # # Area under Curve: 0.006992314 
  # # Maximum F1 Score: 0.1854636 
  # # 
  # 
  
}


# #the StARS summary statistic is:
# getStability(se.usvi_top100_genus.gl)
# # 0.04768485
# getStability(se.usvi_top100_genus.mb)
# # 0.04688646
# # getStability(usvi_top100_genus.se.gl2)
# # # 0.0994901
# 
# #the stability threshold can be better attained by adjusting nlambda
# #how many connections?
# # sum(getRefit(usvi_top100_genus.se.gl2))/2
# # #883
# sum(getRefit(se.usvi_top100_genus.gl))/2
# #475
# sum(getRefit(se.usvi_top100_genus.mb))/2
# #193
# 

# getStability(se.usvi_top100_genus.slr)
# # 0.04817584
# #how many connections?
# sum(getRefit(se.usvi_top100_genus.slr))/2
# #291
# 


#network via sparcc:
usvi_top100_genus.sparcc <- SpiecEasi::sparcc(as.matrix(usvi_top100_genus_with_metab.tbl))
usvi_sparcc.graph <- abs(usvi_top100_genus.sparcc$Cor) >= 0.3
diag(usvi_sparcc.graph) <- 0
usvi_sparcc.graph <- Matrix::Matrix(usvi_sparcc.graph, sparse=TRUE) #dsCMatrix

#testing out plots of networks:
{
  # #plot the networks made via sparcc, spiec-easi gl, and spiec-easi mb
  # ig.gl <- adj2igraph(getRefit(se.usvi_top100_genus.gl))
  # ig.sparcc <- adj2igraph(usvi_sparcc.graph)
  # 
  # # # temp_ig.gl <- getRefit(se.usvi_top100_genus.gl) #equivalent to:
  # # temp_ig.gl <- Matrix::drop0(se.usvi_top100_genus.gl$refit$stars) #makes a  dsCMatrix
  # # temp_ig.gl <- igraph::graph_from_adjacency_matrix(temp_ig.gl, mode = "undirected", weighted = TRUE, diag=FALSE)
  # 
  # #the default for adj2igraph() invokes igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected", weighted = TRUE, diag=diag)
  # #but for the mb model, it errors out.
  # 
  # ## ig.mb <- adj2igraph(getRefit(se.usvi_top100_genus.mb)) #no!
  # # temp_ig.mb <- getOptBeta(se.usvi_top100_genus.mb) #makes a dgCMatrix
  # temp_ig.mb <- Matrix::drop0(se.usvi_top100_genus.mb$refit$stars) #makes a dgCMatrix
  # # class(as(temp_ig.mb, "symmetricMatrix")) #no!
  # # as(temp_ig.mb, "dsCMatrix") #also no!
  # temp_ig.mb <- Matrix::forceSymmetric(temp_ig.mb, uplo = "U") #makes it a dsCMatrix!
  # ig.mb <- adj2igraph(temp_ig.mb)
  # 
  # 
  # 
  # 
  # ## set size of vertex proportional to clr-mean
  # vsize <- rowMeans(clr(as.matrix(usvi_top100_genus_with_metab.tbl), 1))+6
  # # am.coord <- igraph::layout_with_fr(ig.gl)
  # am.coord <- igraph::layout_with_fr(ig.mb)
  # # 
  # # par(mfrow=c(1,3))
  # # plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  # # plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
  # # plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
  # # 
  # # dev.off()
  # # 
  
  # #now add taxon names
  # ig2.gl <- adj2igraph(getRefit(se.usvi_top100_genus.gl),
  #                      rmEmptyNodes = TRUE, 
  #                      vertex.attr = list(name = se.usvi_top100_genus.gl[["est"]][["data"]] %>% colnames(.)))
  # ig2.sparcc <- adj2igraph(usvi_sparcc.graph, 
  #                          rmEmptyNodes = TRUE,
  #                          vertex.attr = list(name = colnames(usvi_top100_genus_with_metab.tbl)))
  # ig2.mb <- adj2igraph(temp_ig.mb, 
  #                      rmEmptyNodes = TRUE, 
  #                      vertex.attr = list(name = se.usvi_top100_genus.mb[["est"]][["data"]] %>% colnames(.)))
  # ig2.slr <- adj2igraph(getRefit(se.usvi_top100_genus.slr),
  #                       rmEmptyNodes = TRUE,
  #                       vertex.attr = list(name = se.usvi_top100_genus.slr[["est"]][["data"]] %>% colnames(.)))
  # phyloseq::plot_network(ig2.slr)
  # phyloseq::plot_network(ig2.gl)
  # phyloseq::plot_network(ig2.sparcc)
  # phyloseq::plot_network(ig2.mb)
  # 
  # 
  # 
  # 
  # #edge dissimilarities
  # #compare edge dissimilarities between the two types of networks:
  # otu1 <- se.usvi_top100_genus.gl[["est"]][["data"]] %>% colnames(.)
  # otu2 <- colnames(usvi_top100_genus_with_metab.tbl)
  # 
  # #for reference, 1 = most similar; 0 = most dissimilar
  # edge.diss(temp_ig.mb, temp_ig.mb, 'jaccard', otu1, otu1)
  # # 1
  # edge.diss(getRefit(se.usvi_top100_genus.gl), getRefit(se.usvi_top100_genus.gl), 'jaccard', otu1, otu1)
  # # 1
  # edge.diss(getRefit(se.usvi_top100_genus.gl), temp_ig.mb, 'jaccard', otu1, otu1)
  # # 0.2393321
  # edge.diss(getRefit(se.usvi_top100_genus.gl), getRefit(se.usvi_top100_genus.slr), 'jaccard', otu1, otu1)
  # # 0.1857585
  # edge.diss(getRefit(se.usvi_top100_genus.slr), usvi_sparcc.graph, 'jaccard', otu1, otu2)
  # # 0.1567398
  # edge.diss(temp_ig.mb, getRefit(se.usvi_top100_genus.slr), 'jaccard', otu1, otu1)
  # # 0.1415094
  # edge.diss(getRefit(se.usvi_top100_genus.gl), usvi_sparcc.graph, 'jaccard', otu1, otu2)
  # # 0.1355236
  # edge.diss(temp_ig.mb, usvi_sparcc.graph, 'jaccard', otu1, otu2)
  # # 0.0627451
  # 
  # 
}



# Make a species-species matrix, and species-metabolite matrix ------------

# usvi_samples_with_metab.mat <- usvi_top100_genus_with_metab.tbl %>%
#   tidyr::drop_na(.) %>%
#   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
#                  autotransform = TRUE) %>%
#   as.matrix(.)
# 
# usvi_metab.mat <- usvi_metab.tbl %>% 
#   as.matrix(.) %>%
#   vegan::vegdist(., binary = FALSE, upper = TRUE,
#                  distance = "horn", autotransform = TRUE) %>%
#   as.matrix(.)
# 
# #make a dissimilarity matrix for the ASVs
# usvi_top100_genus_with_metab.mat <- usvi_top100_genus_with_metab.tbl %>%
#   tibble::rownames_to_column("asv_id") %>%
#   data.table::transpose(., keep.names = "asv_id", make.names = TRUE) %>% droplevels %>%
#   tibble::column_to_rownames(var = "asv_id") %>%
#   tidyr::drop_na(.) %>%
#   vegan::vegdist(., distance = "horn", binary = FALSE, upper = TRUE,
#                  autotransform = TRUE) %>%
#   as.matrix(.)
# 
# 
# 
# #make a correlation matrix between genera and metabolites
# 
# usvi_top100_metabolites.cor <- matrix(nrow = ncol(usvi_top100_genus_with_metab.tbl), ncol = ncol(usvi_metab.tbl))
# colnames(usvi_top100_metabolites.cor) <- colnames(usvi_metab.tbl)
# rownames(usvi_top100_metabolites.cor) <- colnames(usvi_top100_genus_with_metab.tbl)
# 
# usvi_top100_metabolites.cor.rho <- usvi_top100_metabolites.cor
# 
# y <- length(colnames(usvi_top100_metabolites.cor))
# for(j in seq_len(y)){
#   vector_metab <- usvi_metab.tbl[, j]
#   # for(i in seq_len(2)){
#   for(i in seq_len(nrow(usvi_top100_metabolites.cor))){
#     usvi_top100_metabolites.cor[i, j] <- cor.test(vector_metab, usvi_top100_genus_with_metab.tbl[,i], method = "spearman", exact = FALSE) %>%
#       purrr::pluck(., "p.value")
#     usvi_top100_metabolites.cor.rho[i, j] <- cor.test(vector_metab, usvi_top100_genus_with_metab.tbl[,i], method = "spearman", exact = FALSE) %>%
#       purrr::pluck(., "estimate")
#   }
# }
# {
#   q_value <- 0.1
#   padj_cutoff <- usvi_top100_metabolites.cor %>%
#     apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
#     apply(., 2, function(x) ashr::qval.from.lfdr(x)) %>%
#     unlist %>%
#     as.matrix(.) %>%
#     quantile(., probs = q_value, na.rm = TRUE, names = FALSE,type = 7)
#   
# }
# 
# usvi_top100_metabolites.cor.corrected <- usvi_top100_metabolites.cor %>%
#   apply(., 2, function(x) p.adjust(x, method = "BH")) %>% #multiple testing corrections
#   # apply(., 2, function(x) ifelse(x < 0.05, x, NA)) #drop the p.values > 0.05 or did not compute
#   apply(., 2, function(x) ifelse(x < padj_cutoff, x, NA)) #drop the p.values > the adjusted p-value or did not compute



#try spiecieasi with the metabolites data
usvi_metab_lc <- usvi_metab.tbl %>%
  dplyr::mutate(across(everything(), ~ifelse(!(is.na(.x) | (.x < 0)),
                                  log2(.x+1), #log transform abundance (with +1 pseudocount)
                                  0))) %>%
  droplevels %>%
  apply(., 2, function(x) scale(x, center = FALSE, scale = FALSE), simplify = TRUE)
rownames(usvi_metab_lc) <- rownames(usvi_metab.tbl)
if(!file.exists(paste0(projectpath, "/","se.usvi_metab_genus.gl", ".rds"))){
  se.usvi_metab_genus.gl <- spiec.easi(list(as.matrix(usvi_top100_genus_with_metab.tbl), usvi_metab_lc), 
                                       method = "glasso", 
                                       pulsar.select = TRUE, 
                                       lambda.min.ratio = 0.01, 
                                       nlambda = 100, 
                                       pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, ncores = nthreads))
  readr::write_rds(se.usvi_metab_genus.gl, paste0(projectpath, "/","se.usvi_metab_genus.gl", ".rds"), compress = "gz")
} else {
  se.usvi_metab_genus.gl <- readr::read_rds(paste0(projectpath, "/","se.usvi_metab_genus.gl", ".rds"))
}

# getStability(se.usvi_metab_genus.gl)
# # 0.04865092
# sum(getRefit(se.usvi_metab_genus.gl))/2
# #1049
# # 
#            
# ig.metab.gl <- adj2igraph(getRefit(se.usvi_metab_genus.gl))
# # plot(ig.metab.gl, vertex.color = dtype+1)
# # phyloseq::plot_network(ig.metab.gl)
# 
# ig2.metab.gl <- adj2igraph(getRefit(se.usvi_metab_genus.gl),
#                       rmEmptyNodes = TRUE,
#                       vertex.attr = list(name = se.usvi_metab_genus.gl[["est"]][["data"]] %>% colnames(.)))
# dtype <- c(rep(1,ncol(usvi_top100_genus_with_metab.tbl)), rep(2,ncol(usvi_metab_lc)))
# plot(ig2.metab.gl, vertex.color = dtype+1)
# phyloseq::plot_network(ig2.metab.gl)

if(!file.exists(paste0(projectpath, "/","se.usvi_metab_genus.mb", ".rds"))){
  se.usvi_metab_genus.mb <- spiec.easi(list(as.matrix(usvi_top100_genus_with_metab.tbl), usvi_metab_lc), 
                                       method = "mb", 
                                       pulsar.select = TRUE, 
                                       lambda.min.ratio = 0.01, 
                                       nlambda = 20, 
                                       pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, ncores = nthreads))
  readr::write_rds(se.usvi_metab_genus.mb, paste0(projectpath, "/","se.usvi_metab_genus.mb", ".rds"), compress = "gz")
} else {
  se.usvi_metab_genus.mb <- readr::read_rds(paste0(projectpath, "/","se.usvi_metab_genus.mb", ".rds"))
}
# getStability(se.usvi_metab_genus.mb)
# # 0.047807
# sum(getRefit(se.usvi_metab_genus.mb))/2
# #453

se.metab.mb <- getOptBeta(se.usvi_metab_genus.mb) #dgCmatrix
#which one to use?

se.metab.mb_cor <- se.metab.mb %>%
  as.matrix(.) %>% symBeta(., mode = "ave") #dsyMatrix


{
  temp_m <- as.matrix(se.metab.mb)
  temp_u <- Matrix::triu(temp_m) %>% as.matrix(.)
  temp_l <- Matrix::tril(temp_m) %>% as.matrix(.) %>% t(.)
  
  #ave method:
  temp_beta_ave <- (temp_m+t(temp_m))/2
  # all(temp_beta_ave == as.matrix(se.metab.mb_cor))
  # TRUE
  colnames(temp_beta_ave) <- rownames(temp_beta_ave) <- colnames(se.usvi_metab_genus.gl[["est"]][["data"]])
  
  #maxabs method:
  {
    maxt <- pmax(abs(temp_u), abs(temp_l))
    uptind <- Matrix::which(maxt == abs(temp_u))
    lotind <- Matrix::which(maxt == abs(temp_l))
    
    if (length(uptind != 0)) maxt[uptind] <- maxt[uptind]*sign(temp_u[uptind])
    if (length(lotind != 0)) maxt[lotind] <- maxt[lotind]*sign(temp_l[lotind])
    temp_beta_ma <-  maxt + t(maxt)
    
  }
  colnames(temp_beta_ma) <- rownames(temp_beta_ma) <- colnames(se.usvi_metab_genus.gl[["est"]][["data"]])
  # all(temp_beta_ma == as.matrix(se.metab.mb_cor))
  
  
  
  #upper method:
  temp_beta_u <- temp_u + t(temp_u)
  # all(temp_beta_u == as.matrix(se.metab.mb_cor))
  colnames(temp_beta_u) <- rownames(temp_beta_u) <- colnames(se.usvi_metab_genus.gl[["est"]][["data"]])
  
  
  #lower method:
  temp_beta_l <- temp_l + t(temp_l)
  # all(temp_beta_l == as.matrix(se.metab.mb_cor))
  colnames(temp_beta_l) <- rownames(temp_beta_l) <- colnames(se.usvi_metab_genus.gl[["est"]][["data"]])
  # 
  # 
  # temp_beta <- temp_u + t(temp_l)
  # temp_beta <- as.matrix(temp_beta)
  # temp_beta <- Matrix::forceSymmetric(temp_beta)
  
  }
# NetCoMi::plotHeat(temp_beta_ave, type = "upper")
# corrplot::corrplot(temp_beta_ave, type = "upper", method = "color", add = FALSE)


# Prepare for network construction ----------------------------------------
nodeCols <- c(rep("lightblue", length(colnames(as.matrix(usvi_top100_genus_with_metab.tbl)))), rep("orange", length(colnames(usvi_metab_lc))))
names(nodeCols) <- c(colnames(as.matrix(usvi_top100_genus_with_metab.tbl)), (colnames(usvi_metab_lc)))

se.metab.mb_cor <- temp_beta_ave
# se.metab.mb_cor <- temp_beta_u
# se.metab.mb_cor <- temp_beta_l
# se.metab.mb_cor <- temp_beta_ma

quantile(se.metab.mb_cor[se.metab.mb_cor != 0], probs = seq(0, 1,0.25))
quantile(abs(se.metab.mb_cor[se.metab.mb_cor != 0]), probs = seq(0, 1,0.25))
quantile(abs(se.metab.mb_cor[se.metab.mb_cor != 0]), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)

# temp_cutoff <- quantile(abs(se.metab.mb_cor), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)
# temp_mat[(abs(temp_mat) < 0.3)] <- 0
se.metab.mb_cor[(abs(se.metab.mb_cor) < temp_cutoff)] <- 0
diag(se.metab.mb_cor) <- 0


net_se.metab.mb <- NetCoMi::netConstruct(data = se.metab.mb_cor, 
                                         dataType = "condDependence",
                                         sparsMethod = "none",
                                         cores = nthreads,
                                         verbose = 1,
                                         seed = 48105)
netcom_se.metab.mb <- netAnalyze(net_se.metab.mb, 
                                 # centrLCC = FALSE, sPathAlgo = "automatic",
                                 # hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
                                 # normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
                                 clustMethod = "hierarchical")
summary(netcom_se.metab.mb)

netcom_se.metab.mb %>%
  plot(.,
       layout = "spring",
       repulsion = 1.2,
       shortenLabels = "none",
       labelScale = TRUE,
       rmSingles = TRUE,
       nodeSize = "eigenvector",
       nodeSizeSpread = 2,
       # nodeColor = "cluster",
       nodeColor = "colorVec", colorVec = nodeCols,
       hubBorderCol = "gray60",
       cexNodes = 1.8,
       cexLabels = 2,
       cexHubLabels = 2.2,
       title1 = "Network for Genera/Metabolite data", 
       showTitle = TRUE,
       cexTitle = 2.3)

#trythe glasso networks:
se.metab.gl <- getOptCov(se.usvi_metab_genus.gl)
se.metab.gl_cov <- se.metab.gl %>%
  as.matrix(.) %>%
  cov2cor(.)
rownames(se.metab.gl_cov) <- colnames(se.usvi_metab_genus.gl[["est"]][["data"]])
colnames(se.metab.gl_cov) <- colnames(se.usvi_metab_genus.gl[["est"]][["data"]])

net_se.metab.gl <- NetCoMi::netConstruct(data = se.metab.gl_cov, 
                                         dataType = "condDependence",
                                         sparsMethod = "none",
                                         cores = nthreads,
                                         verbose = 1,
                                         seed = 48105)
netcom_se.metab.gl <- netAnalyze(net_se.metab.gl, 
                                 # centrLCC = FALSE, sPathAlgo = "automatic",
                                 # hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
                                 # normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
                                 clustMethod = "hierarchical")
summary(netcom_se.metab.gl)


#what if we removed all values that fell below threshold?
temp_mat <- se.metab.gl_cov
temp_cutoff <- quantile(abs(temp_mat), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)
# temp_mat[(abs(temp_mat) < 0.3)] <- 0
temp_mat[(abs(temp_mat) < temp_cutoff)] <- 0
diag(temp_mat) <- 0

net_se.metab.gl_filter <- NetCoMi::netConstruct(data = temp_mat, 
                                         # data = se.metab.gl_cov, 
                                         dataType = "condDependence",
                                         sparsMethod = "none",
                                         cores = nthreads,
                                         verbose = 1,
                                         seed = 48105)
netcom_se.metab.gl_filter <- netAnalyze(net_se.metab.gl_filter, 
                                 # centrLCC = FALSE, sPathAlgo = "automatic",
                                 # hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
                                 # normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
                                 # clustPar = list(method = "average", k = 3),
                                 hubPar = "eigenvector",
                                 clustMethod = "hierarchical")
                                 # clustMethod = "cluster_fast_greedy")

plot(netcom_se.metab.gl_filter,
     layout = "spring",
     # repulsion = 1.2,
     # nodeColor = "cluster",
     nodeColor = "colorVec", colorVec = nodeCols,
     nodeSize = "fix",
     nodeSizeSpread = 2,
     hubTransp = 40,
     highlightHubs = TRUE,
     edgeTranspLow = 60,
     rmSingles = TRUE,
     # edgeFilter = "threshold", edgeFilterPar = temp_cutoff,
     shortenLabels = "simple",
     labelLength = 6,
     mar = c(1, 3, 3, 5))

# netcom_se.metab.gl_filter %>%
netcom_se.metab.gl %>%
  plot(.,
       layout = "spring",
       repulsion = 1.2,
       shortenLabels = "simple",
       labelLength = 6,
       # labelScale = TRUE,
       rmSingles = TRUE,
       nodeSize = "fix",
       nodeSizeSpread = 2,
       # nodeColor = "cluster",
       nodeColor = "colorVec", colorVec = nodeCols,
       hubBorderCol = "gray60",
       hubTransp = 40,
       highlightHubs = TRUE,
       cexNodes = 1.8, cexLabels = 2, cexHubLabels = 2.2,
       edgeFilter = "threshold", edgeFilterPar = temp_cutoff,
       # edgeInvisFilter = "threshold", edgeInvisPar = temp_cutoff,
       title1 = "Network for Genera/Metabolite data", 
       showTitle = TRUE,
       cexTitle = 2.3)


net1 <- NetCoMi::netConstruct(data = usvi_top100_genus.cor,
                              data2 = usvi_top100_metabolites.cor,
                              dataType = "correlation",
                                  measure = "spring",
                              measurePar = list(nlambda=10,
                                                rep.num=10,
                                                Rmethod = "approx"),
                              filtTax = "none",
                                  # filtTax = "highestFreq", filtTaxPar = list(highestFreq = 50),
                                  zeroMethod = "none",
                                  normMethod = "none",
                                  sparsMethod = "none",
                              dissFunc = "signed",
                              cores = 1,
                                  # cores = nthreads,
                                  verbose = 1,
                                  seed = 48105)


# net2 <- NetCoMi::netConstruct(data = usvi_samples_with_metab.mat,
#                               # data2 = usvi_metab.mat,
#                               # group = usvi_metab.groups,
#                               # measure = "spearman", zeroMethod = "none", normMethod = "none", sparsMethod = "none",
#                                   dataType = "dissimilarity", 
#                               measure = "aitchison", zeroMethod = "multRepl", sparsMethod = "knn", kNeighbor = 3,
#                               filtTax = "none",
#                                   # filtTax = "highestFreq", filtTaxPar = list(highestFreq = 50),
#                                   
#                               dissFunc = "signed",
#                               cores = 1,
#                                   # cores = nthreads,
#                                   verbose = 1,
#                                   seed = 48105)

net2 <- NetCoMi::netConstruct(data = usvi_top100_genus_with_metab.tbl, dataType = "counts",
  # data = usvi_top100_genus_with_metab.mat, dataType = "dissimilarity", 
                              # data2 = usvi_metab.mat,
                              group = usvi_metab.groups,
  measure = "spring", measurePar = list(nlambda = 10, rep.num = 10, Rmethod = "approx"), 
                              # measure = "spearman",
  zeroMethod = "none", normMethod = "none", sparsMethod = "none",
                              # measure = "aitchison", zeroMethod = "multRepl", sparsMethod = "knn", kNeighbor = 3,
                              filtTax = "none",
                              # filtTax = "highestFreq", filtTaxPar = list(highestFreq = 50),
                              
                              dissFunc = "signed",
                              cores = 1,
                              # cores = nthreads,
                              verbose = 1,
                              seed = 48105)

props2 <- netAnalyze(net2, 
                     centrLCC = FALSE, sPathAlgo = "automatic",
                     hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
                     normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
                     clustMethod = "cluster_fast_greedy")
summary(props2)
plot(props2)
