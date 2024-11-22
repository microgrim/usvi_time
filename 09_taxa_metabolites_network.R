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

# #run PCA on these data
## not enough info to suggest pcaMethods works for the transformed data out of spiec-easi
{
  # se.mb.data <- se.usvi_metab_genus.mb$est$data #a 71x149. (num of samples x (100 genera + 49 metabolites))
  # temp_pca <- pcaMethods::pca(se.mb.data, method = "svd", nPcs = 2, cv = "q2", scale = "none", center = TRUE)
  # summary(temp_pca)
  # # plot(temp_pca)
  # pcaMethods::slplot(temp_pca)
  # 
  # se.mb.data_t <- t(se.usvi_metab_genus.mb$est$data) #a 149x71. ((100 genera + 49 metabolites) x num of samples )
  # temp_pca2 <- pcaMethods::pca(se.mb.data_t, method="nlpca", nPcs = 2, maxSteps = 100, cv = "q2")
  # summary(temp_pca2)
  # pcaMethods::slplot(temp_pca2)
  # pcaMethods::scores(temp_pca2)
  # # pcaMethods::loadings(temp_pca2)
  # # fittedData <- fitted(temp_pca)
  # plot(fitted(temp_pca2))
  
  
}

#can we do robust PCA on these data?
#no.

{
 
# X <- se.usvi_metab_genus.mb$est$data #a 71x149. (num of samples x (100 genera + 49 metabolites))
# L <- se.usvi_metab_genus.mb$est$beta[[getOptInd(se.usvi_metab_genus.mb)]] # a 149x149 matrix (100 genera + 49 metabolites)
# se.usvi_metab_genus.mb.pca <- robustPCA(X, L)
# 
# age <- as.numeric(as.character(sample_data(amgut2.filt.phy)$AGE))
# bmi <- as.numeric(as.character(sample_data(amgut2.filt.phy)$BMI))
# depth <- colSums(otu_table(amgut2.filt.phy))
# 
# cor(age, se.usvi_metab_genus.mb.pca$scores, use='pairwise')
# cor(bmi, se.usvi_metab_genus.mb.pca$scores, use='pairwise')

}

se.metab.mb <- getOptBeta(se.usvi_metab_genus.mb) #dgCmatrix
#which one to use?

# se.metab.mb_cor <- se.metab.mb %>%
#   as.matrix(.) %>% symBeta(., mode = "ave") #dsyMatrix

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
  
  }
# NetCoMi::plotHeat(temp_beta_ave, type = "upper")
# corrplot::corrplot(temp_beta_ave, type = "upper", method = "color", add = FALSE)


# Use NetCoMi -------------------------------------------------------------

nodeCols <- c(rep("lightblue", length(colnames(as.matrix(usvi_top100_genus_with_metab.tbl)))), rep("orange", length(colnames(usvi_metab_lc))))
names(nodeCols) <- c(colnames(as.matrix(usvi_top100_genus_with_metab.tbl)), (colnames(usvi_metab_lc)))

se.metab.mb_cor <- temp_beta_ave
# #with average beta, in the node graph dataframe 130 genera and metabolites are retained, in one cluster only and weight ranging from 0.21 to 0.57

# se.metab.mb_cor <- temp_beta_u
# se.metab.mb_cor <- temp_beta_l
# se.metab.mb_cor <- temp_beta_ma
# #with maxabs, in the node graph dataframe 130 genera and metabolites are retained, in 3 clusters and weight ranged from 0.14-0.60

se.mb_cutoff <- quantile(abs(se.metab.mb_cor[se.metab.mb_cor != 0]), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)

# se.mb_cutoff <- quantile(abs(se.metab.mb_cor), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)
# temp_mat[(abs(temp_mat) < 0.3)] <- 0
# se.metab.mb_cor[(abs(se.metab.mb_cor) < se.mb_cutoff)] <- 0
# diag(se.metab.mb_cor) <- 0


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
# summary(netcom_se.metab.mb)
# 
# p1 <- netcom_se.metab.mb %>%
#   plot(.,
#        layout = "spring",
#        repulsion = 1.2,
#        shortenLabels = "none",
#        labelScale = TRUE,
#        rmSingles = TRUE,
#        nodeSize = "eigenvector",
#        nodeSizeSpread = 2,
#        # nodeColor = "cluster",
#        nodeColor = "colorVec", colorVec = nodeCols,
#        # edgeFilter = "threshold", edgeFilterPar = se.mb_cutoff,
#        hubBorderCol = "gray60",
#        cexNodes = 1.8,
#        cexLabels = 2,
#        cexHubLabels = 2.2,
#        title1 = "Network for Genera/Metabolite data",
#        showTitle = TRUE,
#        cexTitle = 2.3)

#extract the edges and weights for plotting:
{
  net_mb.edges <- dplyr::select(net_se.metab.mb$edgelist1, v1, v2)
  # net_mb.edges$Source <- as.numeric(factor(net_mb.edges$v1))
  # net_mb.edges$Target <- as.numeric(factor(net_mb.edges$v2))
  # net_mb.edges$Type <- "Undirected"
  # net_mb.edges$Weight <- net_se.metab.mb$edgelist1$adja
  # Add Source and Target variables (as IDs)
  net_mb.nodes <- data.frame(Label = union(net_mb.edges$v1, net_mb.edges$v2),
                             Id = seq_len(length(union(net_mb.edges$v1, net_mb.edges$v2))))
  net_mb.edges <- net_mb.edges %>%
    dplyr::inner_join(net_mb.nodes, by = join_by("v1" == "Label")) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(net_mb.nodes, by = join_by("v2" == "Label")) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = net_se.metab.mb$edgelist1$adja)
  
  
  # net_mb.nodes <- unique(net_mb.edges[,c('v1','Source')])
  # colnames(net_mb.nodes) <- c("Label", "Id")
  # 
  # Add category with clusters (can be used as node colors in Gephi)
  net_mb.nodes$Category <- netcom_se.metab.mb$clustering$clust1[net_mb.nodes$Label]
  net_mb.edges <- dplyr::select(net_mb.edges, Source, Target, Type, Weight)
  
  readr::write_delim(net_mb.edges, paste0(projectpath, "/", "net_se.mb.edges", ".csv"),
                     delim = ",", col_names = TRUE)
  readr::write_delim(net_mb.nodes, paste0(projectpath, "/", "net_se.mb.nodes", ".csv"),
                     delim = ",", col_names = TRUE)
  
}




#try the glasso networks:
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

{
  net_gl.edges <- dplyr::select(net_se.metab.gl$edgelist1, v1, v2)
  
  # Add Source and Target variables (as IDs)
  # net_gl.edges$Source <- as.numeric(factor(net_gl.edges$v1))
  # net_gl.edges$Target <- as.numeric(factor(net_gl.edges$v2))
  # net_gl.edges$Type <- "Undirected"
  # net_gl.edges$Weight <- net_se.metab.gl$edgelist1$adja
  # 
  # net_gl.nodes <- unique(net_gl.edges[,c('v1','Source')])
  # colnames(net_gl.nodes) <- c("Label", "Id")
  # 
  net_gl.nodes <- data.frame(Label = union(net_gl.edges$v1, net_gl.edges$v2),
                             Id = seq_len(length(union(net_gl.edges$v1, net_gl.edges$v2))))
  net_gl.edges <- net_gl.edges %>%
    dplyr::inner_join(net_gl.nodes, by = join_by("v1" == "Label")) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(net_gl.nodes, by = join_by("v2" == "Label")) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = net_se.metab.gl$edgelist1$adja)
  
  # Add category with clusters (can be used as node colors in Gephi)
  net_gl.nodes$Category <- netcom_se.metab.gl$clustering$clust1[net_gl.nodes$Label]
  net_gl.edges <- dplyr::select(net_gl.edges, Source, Target, Type, Weight)
  
  readr::write_delim(net_gl.edges, paste0(projectpath, "/", "net_se.gl.edges", ".csv"),
                     delim = ",", col_names = TRUE)
  readr::write_delim(net_gl.nodes, paste0(projectpath, "/", "net_se.gl.nodes", ".csv"),
                     delim = ",", col_names = TRUE)
  
}


# #what if we removed all values that fell below threshold?
# temp_mat <- se.metab.gl_cov
# se.gl_cutoff <- quantile(abs(temp_mat), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)
# # temp_mat[(abs(temp_mat) < 0.3)] <- 0
# temp_mat[(abs(temp_mat) < se.gl_cutoff)] <- 0
# diag(temp_mat) <- 0
# 
# net_se.metab.gl_filter <- NetCoMi::netConstruct(data = temp_mat, 
#                                          # data = se.metab.gl_cov, 
#                                          dataType = "condDependence",
#                                          sparsMethod = "none",
#                                          cores = nthreads,
#                                          verbose = 1,
#                                          seed = 48105)
# netcom_se.metab.gl_filter <- netAnalyze(net_se.metab.gl_filter, 
#                                  # centrLCC = FALSE, sPathAlgo = "automatic",
#                                  # hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
#                                  # normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
#                                  # clustPar = list(method = "average", k = 3),
#                                  hubPar = "eigenvector",
#                                  clustMethod = "hierarchical")
#                                  # clustMethod = "cluster_fast_greedy")
# 
# plot(netcom_se.metab.gl_filter,
#      layout = "spring",
#      # repulsion = 1.2,
#      # nodeColor = "cluster",
#      nodeColor = "colorVec", colorVec = nodeCols,
#      nodeSize = "fix",
#      nodeSizeSpread = 2,
#      hubTransp = 40,
#      highlightHubs = TRUE,
#      edgeTranspLow = 60,
#      rmSingles = TRUE,
#      # edgeFilter = "threshold", edgeFilterPar = se.gl_cutoff,
#      shortenLabels = "simple",
#      labelLength = 6,
#      mar = c(1, 3, 3, 5))
# 
# # netcom_se.metab.gl_filter %>%
# netcom_se.metab.gl %>%
#   plot(.,
#        layout = "spring",
#        repulsion = 1.2,
#        shortenLabels = "simple",
#        labelLength = 6,
#        # labelScale = TRUE,
#        rmSingles = TRUE,
#        nodeSize = "fix",
#        nodeSizeSpread = 2,
#        # nodeColor = "cluster",
#        nodeColor = "colorVec", colorVec = nodeCols,
#        hubBorderCol = "gray60",
#        hubTransp = 40,
#        highlightHubs = TRUE,
#        cexNodes = 1.8, cexLabels = 2, cexHubLabels = 2.2,
#        edgeFilter = "threshold", edgeFilterPar = temp_cutoff,
#        # edgeInvisFilter = "threshold", edgeInvisPar = temp_cutoff,
#        title1 = "Network for Genera/Metabolite data", 
#        showTitle = TRUE,
#        cexTitle = 2.3)

# Try using two groups in network analyses --------------------------------

#first, LB seagrass vs Reefs
sample_idx_seagrass <- metabolomics_sample_metadata %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::filter(grepl("seagrass", site)) %>%
  dplyr::distinct(sample_id) %>%
  unlist %>%
  as.character(.)
sample_idx_reefs <- metabolomics_sample_metadata %>%
  dplyr::filter(grepl("seawater", sample_type)) %>%
  dplyr::filter(!grepl("seagrass", site)) %>%
  dplyr::distinct(sample_id) %>%
  unlist %>%
  as.character(.)

metab_seagrass_lc <- usvi_metab_lc[rownames(usvi_metab_lc) %in% sample_idx_seagrass,]
top100_genus_seagrass.tbl <- usvi_top100_genus_with_metab.tbl[rownames(usvi_top100_genus_with_metab.tbl) %in% sample_idx_seagrass,]

metab_reefs_lc <- usvi_metab_lc[rownames(usvi_metab_lc) %in% sample_idx_reefs,]
top100_genus_reefs.tbl <- usvi_top100_genus_with_metab.tbl[rownames(usvi_top100_genus_with_metab.tbl) %in% sample_idx_reefs,]


#mb method
#use nlamba = 50 for seagrass, and nlambda = 200 for reefs

if(!file.exists(paste0(projectpath, "/","se.usvi_seagrass.mb", ".rds"))){
  se.usvi_seagrass.mb <- spiec.easi(list(as.matrix(top100_genus_seagrass.tbl), metab_seagrass_lc), 
                                       method = "mb", 
                                       pulsar.select = TRUE, 
                                       lambda.min.ratio = 0.01, 
                                       nlambda = 50, 
                                       pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, ncores = nthreads))
  readr::write_rds(se.usvi_seagrass.mb, paste0(projectpath, "/","se.usvi_seagrass.mb", ".rds"), compress = "gz")
} else {
  se.usvi_seagrass.mb <- readr::read_rds(paste0(projectpath, "/","se.usvi_seagrass.mb", ".rds"))
}




if(!file.exists(paste0(projectpath, "/","se.usvi_reefs.mb", ".rds"))){
  se.usvi_reefs.mb <- spiec.easi(list(as.matrix(top100_genus_reefs.tbl), metab_reefs_lc), 
                                    method = "mb", 
                                    pulsar.select = TRUE, 
                                    lambda.min.ratio = 0.01, 
                                    nlambda = 200, 
                                    pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, rep.num = 20, seed = 48105, ncores = nthreads))
  readr::write_rds(se.usvi_reefs.mb, paste0(projectpath, "/","se.usvi_reefs.mb", ".rds"), compress = "gz")
} else {
  se.usvi_reefs.mb <- readr::read_rds(paste0(projectpath, "/","se.usvi_reefs.mb", ".rds"))
}

#assess quality of the network:
{
  
  # #the StARS summary statistic is:
  # getStability(se.usvi_seagrass.mb)
  # # 0.04665427 with nlambda = 20
  # #0.04827861 with nlambda = 50
  # 
  # sum(getRefit(se.usvi_seagrass.mb))/2
  # #287 with nlambda = 20
  # #299 with nlambda = 50
  # 
  # 
  # getStability(se.usvi_reefs.mb)
  # # 0.04231906 with nlambda = 20
  # # 0.04793397 with nlambda = 100
  # # 0.04878016 with nlambda = 200
  # sum(getRefit(se.usvi_reefs.mb))/2
  # # 358 with nlambda = 20
  # # 397 with nlambda = 100
  # # 400 with nlambda = 200
}


#pull. out the beta matrices for each site type, and compile in netconstruct
# se.reefs.mb <- getOptBeta(se.usvi_reefs.mb)
{
  temp_m <- getOptBeta(se.usvi_reefs.mb) %>% as.matrix(.)
  temp_beta_ave <- (temp_m+t(temp_m))/2
  # all(temp_beta_ave == as.matrix(symBeta(temp_m, mode = "ave")))
  # # TRUE
  colnames(temp_beta_ave) <- rownames(temp_beta_ave) <- colnames(se.usvi_reefs.mb[["est"]][["data"]])
  se.reefs.mb <- temp_beta_ave
}
# se.seagrass.mb <- getOptBeta(se.usvi_seagrass.mb) 
{
  temp_m <- getOptBeta(se.usvi_seagrass.mb) %>% as.matrix(.)
  temp_beta_ave <- (temp_m+t(temp_m))/2
  # all(temp_beta_ave == as.matrix(symBeta(temp_m, mode = "ave")))
  # # TRUE
    colnames(temp_beta_ave) <- rownames(temp_beta_ave) <- colnames(se.usvi_seagrass.mb[["est"]][["data"]])
  se.seagrass.mb <- temp_beta_ave
  
}

#trying to input covariance matrices to netconstruct--not possible from spiec-easi mb method
{
  # #so we can provide covariance matrices to netconstruct
  # #how to convert from the mb correlation matrix to covariance?
  # getOptInd(se.usvi_seagrass.mb)
  # # temp_array <- unlist(se.usvi_seagrass.mb[["est"]][["beta"]]) %>% map(., as.matrix)
  # temp_array <- unlist(se.usvi_seagrass.mb[["est"]][["beta"]]) %>% map(., as.matrix)  %>% array(., c(149, 149, 50))
  # lapply(temp_array, 1, sd)
  # apply(temp_array, c(1, 2), sd)
  # apply(array(se.usvi_seagrass.mb[["est"]][["beta"]], c(149, 149, 50)), c(1,2), sd)
  # 
  # se.seagrass.mb.sds <- lapply(se.usvi_seagrass.mb[["est"]][["beta"]][getOptInd(se.usvi_seagrass.mb):50], sd)
  # se.usvi_seagrass.mb[["est"]][["beta"]]
  # # corpcor::cor2pcor(se.seagrass.mb)
  # 
  # cor2cov(se.seagrass.mb)
}
#for plotting the networks, can use a threshold (edgeFilter) to not show any connections with abs(asso) < threshold
#what should that threshold be?
# se.mb_cutoff <- quantile(c(abs(se.reefs.mb[se.reefs.mb != 0]), abs(se.seagrass.mb[se.seagrass.mb != 0])), probs = 0.75, na.rm = TRUE, names = FALSE, type = 7)

# temp_idx <- (c((se.reefs.mb[se.reefs.mb != 0]), (se.seagrass.mb[se.seagrass.mb != 0])))
# range(temp_idx)
# #in this pair of matrices from symBeta, range of values is [-0.3408765,0.6629936]
# range(abs(temp_idx))
# #[1.698641e-14, 6.629936e-01]
# quantile(temp_idx, probs = seq(0, 1, 0.25), names = FALSE, type = 7)
# # -3.408765e-01  1.808638e-05  4.049252e-02  1.278752e-01  6.629936e-01
# quantile(abs(temp_idx), probs = seq(0, 1, 0.25), names = FALSE, type = 7)
# # 1.698641e-14 2.194431e-02 6.178854e-02 1.422803e-01 6.629936e-01
# 
# # hist(temp_idx, breaks = c(-0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8))

#the value at the 75th percentile, in both the abs(asso) and true(asso) is approximately 0.128-0.142
#at the 50th percentile, it's 0.004 to 0.006
#lets try 0.005


se.mb_cutoff <- 0.005

net_se.sites.mb <- NetCoMi::netConstruct(data = se.reefs.mb,
                                         data2 = se.seagrass.mb,
                                         dataType = "condDependence",
                                         sparsMethod = "none", #no filtering yet, just all the data
                                         cores = nthreads,
                                         verbose = 2,
                                         seed = 48105)

netprops_se.sites.mb <- NetCoMi::netAnalyze(net_se.sites.mb, 
                                 clustMethod = "hierarchical")

# summary(netprops_se.sites.mb)




#plotting:
if(!exists("p1", envir = .GlobalEnv)){
  p1 <- (
    plot(netprops_se.sites.mb,
         layout = "spring",
         repulsion = 1.2,
         shortenLabels = "none",
         labelScale = TRUE,
         rmSingles = TRUE,
         nodeSize = "eigenvector",
         nodeSizeSpread = 2,
         # nodeColor = "cluster",
         nodeColor = "colorVec", colorVec = nodeCols,
         # edgeFilter = "threshold", edgeFilterPar = se.mb_cutoff,
         hubBorderCol = "gray60",
         cexNodes = 1.8,
         cexLabels = 2,
         cexHubLabels = 2.2,
         title1 = "Network for Genera/Metabolite data",
         showTitle = TRUE,
         groupNames = c("Reef sites", "Seagrass"),
         cexTitle = 2.3)
  )
  p1.cutoff <- netprops_se.sites.mb %>%
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
         edgeFilter = "threshold", edgeFilterPar = se.mb_cutoff,
         hubBorderCol = "gray60",
         cexNodes = 1.8,
         cexLabels = 2,
         cexHubLabels = 2.2,
         title1 = "Network for Genera/Metabolite data",
         showTitle = TRUE,
         groupNames = c("Reef sites", "Seagrass"),
         cexTitle = 2.3)
  
  # ggsave(paste0(projectpath, "/", "net_sites.mb.adjgraph-", Sys.Date(), ".png"),
  #        p1,
  #        width = 10, height = 10, units = "in")
}

#in netCompare, can't do permutation tests with the two association matrices from spiec-easi (condDependence method)

netcom_se.sites.mb <- NetCoMi::netCompare(netprops_se.sites.mb, 
                                 adjust = "adaptBH",
                                 cores = nthreads,
                                 testRand = TRUE,
                                 permTest = FALSE,
                                 # permTest = TRUE, nPerm = 100, storeAssoPerm = TRUE, fileStoreAssoPerm = "netcom_se.sites.mb.assoPerm",
                                 seed = 48105)

# summary(netcom_se.sites.mb, groupNames = c("Reef sites", "Seagrass"))


#save the association matrices
if(!file.exists(paste0(projectpath, "/", "usvi_se.asso_mat.reefs.mb", ".tsv"))){
  temp_mat1 <- se.reefs.mb[names(nodeCols), names(nodeCols)] %>%
    tibble::as_tibble(rownames = "var1")
  temp_mat2 <- se.seagrass.mb[names(nodeCols), names(nodeCols)] %>%
    tibble::as_tibble(rownames = "var1")
  readr::write_delim(temp_mat1, paste0(projectpath, "/", "usvi_se.asso_mat.reefs.mb", ".tsv"), col_names = TRUE,
                     delim = "\t")
  readr::write_delim(temp_mat2, paste0(projectpath, "/", "usvi_se.asso_mat.seagrass.mb", ".tsv"), col_names = TRUE,
                     delim = "\t")
}

# Work with the association matrices --------------------------------------


#visualize the asso matrices as heatmaps
se.reefs.mb.dend <- net_se.sites.mb[["assoMat1"]] %>%
  dist(t(.), method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

se.seagrass.mb.dend <- net_se.sites.mb[["assoMat2"]] %>%
  dist(t(.), method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

se.reefs.mb.df <- se.reefs.mb %>%
  tibble::as_tibble(rownames = "v1") %>%
  tidyr::pivot_longer(., cols = !c("v1"),
                      names_to = "v2",
                      values_to = "asso") %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
  dplyr::mutate(asso = dplyr::case_when(asso != 0 ~ asso,
                                        .default = NA)) %>%
  dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                        .default = "taxon")) %>%
  droplevels
# se.reefs.mb.df %>%
#   dplyr::filter(!is.na(asso)) %>%
#   dplyr::group_by(type, v1) %>%
#   dplyr::summarise(num_asso = length(is.na(asso))) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::arrange(desc(num_asso)) %>%
#   dplyr::reframe(distro = quantile(num_asso, probs = seq(0,1, 0.25)), .by = "type")
se.reefs.mb.ranked.df <- se.reefs.mb.df %>%
  dplyr::filter(!is.na(asso)) %>%
  dplyr::group_by(type, v1) %>%
  dplyr::summarise(num_asso = length(is.na(asso))) %>%
  dplyr::left_join(., se.reefs.mb.df %>%
                     dplyr::filter(!is.na(asso)) %>%
                     dplyr::group_by(type, v1) %>%
                     dplyr::summarise(num_asso = length(is.na(asso))) %>%
                     dplyr::ungroup(.) %>%
                     dplyr::arrange(desc(num_asso)) %>%
                     dplyr::reframe(distro = quantile(num_asso, probs = seq(0,1, 0.25)), .by = "type") %>%
                     dplyr::slice_tail(prop = 0.8, by = "type") %>%
                     dplyr::mutate(rank = rep(paste0("group", seq_len(4)), 2)) %>%
                     droplevels, by = join_by("type"), relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  dplyr::filter(num_asso <= distro) %>%
  dplyr::arrange(type, v1, desc(num_asso)) %>%
  dplyr::distinct(v1, .keep_all = TRUE) %>%
  dplyr::select(type, v1, rank) %>%
  dplyr::mutate(v2 = v1) %>%
  droplevels
  
temp_df <- se.reefs.mb.df %>%
  dplyr::left_join(., se.reefs.mb.ranked.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::select(type, v1, rank) %>%
                     dplyr::rename(rank1 = "rank"), 
                  by = join_by(type, v1), relationship = "many-to-many", multiple = "all") %>%
  dplyr::left_join(., se.reefs.mb.ranked.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::select(v2, rank) %>%
                     dplyr::rename(rank2 = "rank"),
                   by = join_by(v2), relationship = "many-to-many", multiple = "all") %>%
  dplyr::mutate(asso = dplyr::case_when(v1 == v2 ~ 1,
                                        .default = asso)) %>%
  # dplyr::filter(!is.na(asso)) %>%
  droplevels

g_hm1 <- print(
  ggplot(data = temp_df %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(x = v1, y = v2, fill = asso), 
              show.legend = TRUE)
  + scale_fill_gradientn(colors = colorRampPalette(pals::ocean.haline(n = 3))(100),
                         aesthetics = "fill", expand = expansion(1.1,1.1),
                         na.value = NA)
                         # limits = c(-1,1),
                         # na.value = "black")
  # + facet_wrap(.~rank2,
  # + facet_grid(rank2~rank1, space = "free",
  #              scales = "free", 
  #              shrink = TRUE, drop = TRUE)
  + theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())
  + labs(x = "var1", y = "var2")
  + guides(color = "none",
           fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = "Association value",
                                 direction = "vertical"))
)
g_hm1 + facet_grid(rank2~rank1, space = "free",
                   scales = "free", 
                   shrink = TRUE, drop = TRUE)

#what if we use pheatmap?
#with just the reefs' samples association matrix:
{
  # temp_mat <- se.reefs.mb[labels(se.reefs.mb.dend), labels(se.reefs.mb.dend)] %>%
  #   as.data.frame(.) %>%
  #   dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
  #                                                        .default = .x))) %>%
  #   droplevels
}

#with both:
{
  temp_mat1 <- se.reefs.mb[names(nodeCols), names(nodeCols)] %>% #cols are from reefs
    Matrix::triu(.) %>%
    as.matrix(.)
  
  temp_mat2 <- se.seagrass.mb[names(nodeCols), names(nodeCols)] %>% #rows are from seagrass
    Matrix::tril(.) %>%
    as.matrix(.)
  
  temp_mat <- temp_mat1 + temp_mat2
  # temp_mat <- temp_mat[labels(se.reefs.mb.dend), labels(se.seagrass.mb.dend)] %>%
  temp_mat <- temp_mat %>%
      as.data.frame(.) %>%
      dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
                                                           .default = .x))) %>%
      droplevels %>%
    as.matrix(.)
  diag(temp_mat) <- 1

}
temp_df2 <- temp_mat1 %>%
  as.data.frame(.) %>%
  dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
                                                       .default = .x))) %>%
  tibble::as_tibble(rownames = "v1") %>%
  tidyr::pivot_longer(., cols = !c("v1"),
                      names_to = "v2",
                      values_to = "asso") %>%
  dplyr::mutate(site = "reef") %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
  dplyr::filter(asso != 0) %>%
  droplevels %>%
  bind_rows(., (temp_mat2 %>%
                  as.data.frame(.) %>%
                  dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
                                                                       .default = .x))) %>%
                  tibble::as_tibble(rownames = "v1") %>%
                  tidyr::pivot_longer(., cols = !c("v1"),
                                      names_to = "v2",
                                      values_to = "asso") %>%
                  dplyr::mutate(site = "seagrass") %>%
                  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
                  dplyr::filter(asso != 0) %>%
                  droplevels)) %>%
  dplyr::arrange(site, v1, v2)


se.sites.mb.ranked.df <- temp_df2 %>%
  dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                        .default = "taxon")) %>%
  # dplyr::filter(!is.na(asso)) %>%
  dplyr::group_by(type, v1, site) %>%
  dplyr::summarise(num_asso = length(is.na(asso))) %>%
  dplyr::left_join(., temp_df2 %>%
                     dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                                           .default = "taxon")) %>%
                     # dplyr::filter(!is.na(asso)) %>%
                     dplyr::group_by(type, v1, site) %>%
                     dplyr::summarise(num_asso = length(is.na(asso))) %>%
                     dplyr::ungroup(.) %>%
                     dplyr::arrange(desc(num_asso)) %>%
                     dplyr::reframe(distro = quantile(num_asso, probs = seq(0,1, 0.25)), .by = c("type", "site")) %>%
                     dplyr::slice_tail(prop = 0.8, by = c("site", "type")) %>%
                     dplyr::mutate(rank = rep(paste0("group", seq_len(4)), 4)) %>%
                     droplevels, by = join_by("type"), relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  dplyr::filter(num_asso <= distro) %>%
  dplyr::arrange(type, v1, desc(num_asso)) %>%
  dplyr::distinct(v1, .keep_all = TRUE) %>%
  dplyr::select(type, v1, rank) %>%
  dplyr::mutate(v2 = v1) %>%
  droplevels

se.sites.mb.ranked.df <- temp_df2 %>%
  dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                        .default = "taxon")) %>%
  dplyr::left_join(., se.sites.mb.ranked.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::select(type, v1, rank), 
                   by = join_by(type, v1), relationship = "many-to-many", multiple = "all") %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols))))

se.sites.mb.ranked.full.df <- se.sites.mb.ranked.df %>%
  tidyr::complete(v2, site, nesting(v1, rank, type)) %>%
  # dplyr::full_join(., se.sites.mb.ranked.df) %>%
  # tidyr::complete(v1, v2, site) %>%
  # dplyr::left_join(., se.sites.mb.ranked.df %>%
  #                    dplyr::ungroup(.) %>%
  #                    # dplyr::select(v1, type, rank) %>%
  #                    tidyr::expand(v1, v2, site) %>%
  #                    droplevels, relationship = "many-to-many", multiple = "all") %>%
  dplyr::select(site, v1, v2, asso, type, rank) %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
  dplyr::arrange(site, v1, rev(v2), rank) %>%
  droplevels

if(!file.exists(paste0(projectpath, "/", "usvi_se.asso_mat.sites.df", ".tsv"))){
  usvi_se.asso_mat.sites.df <- se.sites.mb.ranked.full.df %>%
    dplyr::select(site, v1, v2, asso, type)
  readr::write_delim(usvi_se.asso_mat.sites.df, paste0(projectpath, "/", "usvi_se.asso_mat.sites.df", ".tsv"), col_names = TRUE,
                     delim = "\t")
}

g_hm2 <- print(
  ggplot(data = se.sites.mb.ranked.df %>%
           droplevels)
  # ggplot(data = se.sites.mb.ranked.df %>%
  # ggplot(data = se.sites.mb.ranked.full.df %>%
           # dplyr::filter(grepl("group1", rank)))
  + theme_bw() 
  + geom_tile(aes(x = v1, y = v2, group = site, fill = asso), na.rm = FALSE,
              show.legend = TRUE)
  + scale_x_discrete(name = "var1", drop = FALSE)
  + scale_y_discrete(name = "var2", drop = FALSE)
  + scale_fill_gradientn(colors = colorRampPalette(pals::ocean.haline(n = 3))(100),
                         aesthetics = "fill", expand = expansion(1.1,1.1),
                         na.value = NA)
  # + facet_wrap(rank~site, ncol = 2, scales = "free",
  + theme(axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90),
          strip.text.y = element_text(angle = 90),
          strip.text.x = element_text(angle = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  # + labs(x = "var1", y = "var2")
  + guides(color = "none",
           fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = "Association value",
                                 direction = "vertical"))
)
gpatch <- g_hm2 + theme(axis.text.y = element_blank()) + facet_grid(.~site, space = "free", scales = "fixed",
                   shrink = TRUE, drop = FALSE)
gpatch <- g_hm2  + gpatch + patchwork::plot_layout(guides = "collect")
gpatch <- gpatch + patchwork::plot_annotation(title = "Association matrices between 49 metabolites and 100 top genera",
                                              tag_levels = "A",
                                              subtitle = "In just reef samples compared to just seagrass samples")
ggsave(paste0(projectpath, "/", "usvi_se_sites.mb-", Sys.Date(), ".png"),
       gpatch,
       width = 22, height = 8, units = "in")

#what if I made v2 from reefs, the v1?
#and keep v1 from seagrass

temp_df3 <- temp_mat1 %>%
  as.data.frame(.) %>%
  # dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
  #                                                      .default = .x))) %>%
  tibble::as_tibble(rownames = "v1") %>%
  tidyr::pivot_longer(., cols = !c("v1"),
                      names_to = "v2",
                      values_to = "asso") %>%
  dplyr::mutate(site = "reef") %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
  # dplyr::filter(asso != 0) %>%
  droplevels %>%
  bind_rows(., (temp_mat2 %>%
                  as.data.frame(.) %>%
                  # dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
                  #                                                      .default = .x))) %>%
                  tibble::as_tibble(rownames = "v2") %>%
                  tidyr::pivot_longer(., cols = !c("v2"),
                                      names_to = "v1",
                                      values_to = "asso") %>%
                  dplyr::mutate(site = "seagrass") %>%
                  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
                  # dplyr::filter(asso != 0) %>%
                  droplevels)) %>%
  # dplyr::filter(v1 != v2) %>%
  # dplyr::filter(asso != 0, .by = c("v1", "v2", "site")) %>%
  dplyr::arrange(site, v1, v2)
temp_df3 %>%
  dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                        .default = "taxon")) %>%
  # dplyr::filter(!is.na(asso)) %>%
  dplyr::filter(asso != 0) %>%
  dplyr::group_by(type, v1, site) %>%
  dplyr::summarise(num_asso = length(is.na(asso))) %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(desc(num_asso)) %>%
  dplyr::reframe(distro = quantile(num_asso, probs = seq(0,1, 0.25)), .by = c("type", "site")) %>%
  dplyr::slice_tail(prop = 0.8, by = c("site", "type")) %>%
  dplyr::mutate(rank = rep(paste0("group", seq_len(4)), 4))

temp_ranked.df <- temp_df3 %>%
  dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                        .default = "taxon")) %>%
  # dplyr::filter(!is.na(asso)) %>%
  dplyr::filter(asso != 0) %>%
  dplyr::group_by(type, v1, site) %>%
  dplyr::summarise(num_asso = length(is.na(asso))) %>%
  dplyr::left_join(., temp_df3 %>%
                     dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                                           .default = "taxon")) %>%
                     dplyr::filter(asso != 0) %>%
                     dplyr::group_by(type, v1, site) %>%
                     dplyr::summarise(num_asso = length(is.na(asso))) %>%
                     dplyr::ungroup(.) %>%
                     dplyr::arrange(desc(num_asso)) %>%
                     dplyr::reframe(distro = quantile(num_asso, probs = seq(0,1, 0.25)), .by = c("type", "site")) %>%
                     dplyr::slice_tail(prop = 0.8, by = c("site", "type")) %>%
                     dplyr::mutate(rank = rep(paste0("group", seq_len(4)), 4)) %>%
                     droplevels, by = join_by("type"), relationship = "many-to-many", multiple = "all") %>%
  droplevels %>%
  dplyr::filter(num_asso <= distro) %>%
  dplyr::arrange(type, v1, desc(num_asso)) %>%
  dplyr::distinct(v1, .keep_all = TRUE) %>%
  dplyr::select(type, v1, rank) %>%
  # dplyr::mutate(v2 = v1) %>%
  droplevels

temp_ranked.df <- temp_df3 %>%
  dplyr::mutate(type = dplyr::case_when(v1 %in% colnames(usvi_metab_lc) ~ "metabolite",
                                        .default = "taxon")) %>%
  dplyr::left_join(., temp_ranked.df %>%
                     dplyr::ungroup(.) %>%
                     dplyr::select(type, v1, rank), 
                   by = join_by(type, v1), relationship = "many-to-many", multiple = "all") %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols)))) %>%
  droplevels

temp_group4_mat <- temp_ranked.df %>%
  dplyr::filter(grepl("group4", rank)) %>%
  droplevels %>%
  tidyr::pivot_wider(., id_cols = c("v1", "site", "type", "rank"),
                     names_from = "v2",
                     values_from = "asso",
                     values_fill = NA) %>%
  dplyr::mutate(across(everything(), ~dplyr::case_when(.x == 0 ~ NA,
                                                       .default = .x))) %>%
  tidyr::pivot_longer(., cols = !c("v1", "site", "type", "rank"),
                      names_to = "v2",
                      values_to = "asso") %>%
  dplyr::mutate(across(c("v1", "v2"), ~factor(.x, levels = names(nodeCols))))

g_hm2 <- print(
  ggplot(data = temp_group4_mat %>%
  # ggplot(data = temp_ranked.df %>%
  # ggplot(data = se.sites.mb.ranked.df %>%
           # dplyr::filter(grepl("group4", rank)) %>%
           droplevels)
  + theme_bw() 
  + geom_tile(aes(x = v1, y = v2, group = site, fill = asso), 
              show.legend = TRUE)
  + scale_fill_gradientn(colors = colorRampPalette(pals::ocean.haline(n = 3))(100),
                         aesthetics = "fill", expand = expansion(1.1,1.1),
                         na.value = NA)
  # + facet_wrap(.~site, ncol = 2, scales = "free",
  # # + facet_grid(rows = NULL, cols = vars(site), scales = "free", space = "free",
  #              shrink = TRUE, drop = FALSE)
  + theme(axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90),
          strip.text.y = element_text(angle = 90),
          strip.text.x = element_text(angle = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  + labs(x = "var1", y = "var2")
  + guides(color = "none",
           fill = guide_colorbar(order = 1, ncol = 1, 
                                 title = "Association value",
                                 direction = "vertical"))
)

color_hm <- rev(colorRampPalette(pals::cubehelix(n = 30)[5:25])(30))
g_hm2_ph <- pheatmap::pheatmap(temp_mat,
                               color = colorRampPalette(rev(pals::brewer.spectral(n = 7)))(100),
               # color = color_hm,
               border_color = "black",
               scale = "none",
               na_col = "white",
               cluster_cols = F,
               cluster_rows = F,
               # annotation_col = dend_colorbar_col,
               # annotation_row = dend_colorbar_row,
               # annotation_colors = color_names,
               # cutree_rows = 3,
               # gaps_col = c(6, 12), #this was for just the 13H samples
               # gaps_col = c(22, 44), #for all fractions
               cutree_cols = 3,
               angle_col = 90,
               # main = paste0("Expression of MAG bins, grouped by vent and phylogeny"),
               fontsize = 8,
               # breaks = temp_breaks,
               # legend_breaks = temp_legend,
               # legend_labels = temp_legendLabels,
               legend = FALSE,
               # legend = TRUE,
               annotation_legend = FALSE,
               silent = TRUE,
               show_rownames = T
) %>%
  ggplotify::as.ggplot(.)
g_hm1_ph




# More network stuff ------------------------------------------------------



#do differential network comparison
#first, since we have asso matrices as input, can use only Fisher z-test
#but on first pass, keep all pvalues
net_se.sites.mb_diff1 <- diffnet(net_se.sites.mb, 
        fisherTrans = TRUE,
        diffMethod = "fisherTest", n1 = 149, n2 = 149, #can only use Fisher z-test for asso matrices as input to netConstruct
        alpha = 0.05, 
        adjust = "adaptBH", trueNullMethod = "lfdr",
        # adjust = "lfdr", lfdrThresh = 0.5,
        verbose = TRUE, cores = nthreads, seed = 48105)
# #next use those p-vals as input to diffnet again
# net_se.sites.mb_diff2 <- diffnet(net_se.sites.mb, 
#         fisherTrans = FALSE,
#         diffMethod = "fisherTest", n1 = 149, n2 = 149, #can only use Fisher z-test for asso matrices as input to netConstruct
#         # diffMethod = "permute", nPerm = 100, permPvalsMethod = "pseudo",
#         adjust = "none", pvalsVec = net_se.sites.mb_diff1$pvalsVec, alpha = 0.05,
#         verbose = TRUE, cores = nthreads, seed = 48105)


plot(net_se.sites.mb_diff1,
     layout = "spring", repulsion = 0.8, adjusted = FALSE,
     shortenLabels = "none", legendArgs = list(bty = "n"),
     rmSingles = FALSE,
     nodeColor = nodeCols,
     # legend = FALSE,
     legend = TRUE, legendPos = "topright", legendTitle = "Sites", legendGroupnames = c("Reef sites", "Seagrass"),
     labelScale = TRUE)
{
  
  net_sites.mb.edges <- dplyr::select(net_se.sites.mb$edgelist1, v1, v2)
  
  net_sites.mb.nodes <- data.frame(Label = union(net_sites.mb.edges$v1, net_sites.mb.edges$v2),
                                   Id = seq_len(length(union(net_sites.mb.edges$v1, net_sites.mb.edges$v2))))
  net_sites.mb.edges1 <- net_sites.mb.edges %>%
    dplyr::inner_join(net_sites.mb.nodes, by = join_by("v1" == "Label")) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(net_sites.mb.nodes, by = join_by("v2" == "Label")) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = net_se.sites.mb$edgelist1$adja)%>%
    dplyr::mutate(cluster = "reefs")
  # dplyr::rename_with(cols = everything(), ~paste0("reefs_", .x))
  
  net_sites.mb.edges2 <- dplyr::select(net_se.sites.mb$edgelist2, v1, v2) %>%
    dplyr::inner_join(net_sites.mb.nodes, by = join_by("v1" == "Label")) %>%
    dplyr::rename(Source = "Id") %>%
    dplyr::inner_join(net_sites.mb.nodes, by = join_by("v2" == "Label")) %>%
    dplyr::rename(Target = "Id") %>%
    dplyr::mutate(Type = "Undirected") %>%
    dplyr::mutate(Weight = net_se.sites.mb$edgelist2$adja) %>%
    dplyr::mutate(cluster = "seagrass")
  # dplyr::rename_with(cols = everything(), ~paste0("seagrass_", .x))
  
  # Add category with clusters (can be used as node colors in Gephi)
  net_sites.mb.nodes$reefs_hub <- netcom_se.sites.mb[["properties"]][["clust1"]][net_sites.mb.nodes$Label]
  net_sites.mb.nodes$seagrass_hub <- netcom_se.sites.mb[["properties"]][["clust2"]][net_sites.mb.nodes$Label]
  
  net_sites.mb.edges <- bind_rows(net_sites.mb.edges1, net_sites.mb.edges2) %>%
    dplyr::select(contains("cluster"), contains("Source"), contains("Target"), contains("Type"), contains("Weight"))
  
  readr::write_delim(net_sites.mb.edges, paste0(projectpath, "/", "net_sites.mb.edges", ".csv"),
                     delim = ",", col_names = TRUE)
  readr::write_delim(net_sites.mb.nodes, paste0(projectpath, "/", "net_sites.mb.nodes", ".csv"),
                     delim = ",", col_names = TRUE)
  
}

# More network manipulation -----------------------------------------------


# #what if we pre-filter in the network construction step?
# net_se.sites.mb.tt <- NetCoMi::netConstruct(data = se.reefs.mb,
#                                             data2 = se.seagrass.mb,
#                                             dataType = "condDependence",
#                                             sparsMethod = "t-test", alpha = 0.05, adjust = "adaptBH", sampleSize = c(149, 149),
#                                             # sparsMethod = "none",
#                                             cores = nthreads,
#                                             verbose = 2,
#                                             seed = 48105)
# #this limits the edgelist to 22 from se.reefs.mb, and 23 from seagrass
# 
# netprops_se.sites.mb.tt <- netAnalyze(net_se.sites.mb.tt, 
#                                       # centrLCC = FALSE, sPathAlgo = "automatic",
#                                       # hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
#                                       # normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
#                                       clustMethod = "hierarchical")
# 
# 
# p2 <- netprops_se.sites.mb.tt %>%
#   plot(.,
#        layout = "spring",
#        repulsion = 1.2,
#        shortenLabels = "none",
#        labelScale = TRUE,
#        rmSingles = TRUE,
#        nodeSize = "eigenvector",
#        nodeSizeSpread = 2,
#        # nodeColor = "cluster",
#        nodeColor = "colorVec", colorVec = nodeCols,
#        # edgeFilter = "threshold", edgeFilterPar = se.mb_cutoff,
#        hubBorderCol = "gray60",
#        cexNodes = 1.8,
#        cexLabels = 2,
#        cexHubLabels = 2.2,
#        title1 = "Network for Genera/Metabolite data",
#        showTitle = TRUE,
#        groupNames = c("Reef sites", "Seagrass"),
#        cexTitle = 2.3)
# 
# 
# #if we use threshold method in the netConstruct step:
# net_se.sites.mb.thresh <- NetCoMi::netConstruct(data = se.reefs.mb,
#                                                 data2 = se.seagrass.mb,
#                                                 dataType = "condDependence",
#                                                 sparsMethod = "threshold", thresh = se.mb_cutoff,
#                                                 # sparsMethod = "none",
#                                                 cores = nthreads,
#                                                 verbose = 2,
#                                                 seed = 48105)
# 
# # #this limits the edgelist to 379 from se.reefs.mb, and 274 from seagrass
# # all(abs(net_se.sites.mb.thresh[["edgelist1"]][["asso"]]) > 0.005)
# # all(abs(net_se.sites.mb.thresh[["edgelist1"]][["asso"]]) > se.mb_cutoff)
# 
# netprops_se.sites.mb.thresh <- netAnalyze(net_se.sites.mb.thresh, 
#                                           # centrLCC = FALSE, sPathAlgo = "automatic",
#                                           # hubPar = c("degree", "eigenvector"), gcmHeat = FALSE, hubQuant = 0.9, sPathNorm = FALSE, lnormFit = FALSE,
#                                           # normDeg = FALSE, normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
#                                           clustMethod = "hierarchical")
# 
# p3 <- netprops_se.sites.mb.thresh %>%
#   plot(.,
#        layout = "spring",
#        repulsion = 1.2,
#        shortenLabels = "none",
#        labelScale = TRUE,
#        rmSingles = TRUE,
#        nodeSize = "eigenvector",
#        nodeSizeSpread = 2,
#        # nodeColor = "cluster",
#        nodeColor = "colorVec", colorVec = nodeCols,
#        # edgeFilter = "threshold", edgeFilterPar = 0.1,
#        hubBorderCol = "gray60",
#        cexNodes = 1.8,
#        cexLabels = 2,
#        cexHubLabels = 2.2,
#        title1 = "Network for Genera/Metabolite data",
#        showTitle = TRUE,
#        groupNames = c("Reef sites", "Seagrass"),
#        cexTitle = 2.3)
# 
# # length(which(abs(net_se.sites.mb.thresh[["edgelist1"]][["asso"]]) > 0.1))
# # #145 associations had abs(asso) > 0.1 
# 
# 
# net_se.sites.mb_diff2 <- diffnet(net_se.sites.mb.thresh, 
#                                  fisherTrans = TRUE,
#                                  diffMethod = "fisherTest", n1 = 149, n2 = 149, #can only use Fisher z-test for asso matrices as input to netConstruct
#                                  alpha = 0.05, 
#                                  adjust = "adaptBH", trueNullMethod = "lfdr",
#                                  # adjust = "lfdr", lfdrThresh = 0.5,
#                                  verbose = TRUE, cores = nthreads, seed = 48105)
# 
# plot(net_se.sites.mb_diff2,
#      layout = "spring", repulsion = 0.8, adjusted = FALSE,
#      shortenLabels = "none", legendArgs = list(bty = "n"),
#      rmSingles = FALSE, nodeColor = nodeCols,
#      # legend = FALSE,
#      legend = TRUE, legendPos = "topright", legendTitle = "Sites", legendGroupnames = c("Reef sites", "Seagrass"),
#      labelScale = TRUE)



