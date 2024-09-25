## 00_data_ingest.R

# import the amplicon fastqs for processing through DADA2

# Load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
library(BiocParallel)
library(dada2)
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

# read in files -----------------------------------------------------------


sample_metadata <- readr::read_delim(paste0(projectpath, "/", "usvi_metadata.txt"),
                                     delim = "\t",
                                     quote = "",
                                     col_names = TRUE,
                                     show_col_types = FALSE,
                                     num_threads = nthreads)

# Set up data import manifest ---------------------------------------------

# #raw files are here:
# /proj/omics/apprill/MiSeq/LibraryW2
# #some are also here:
# /proj/omics/apprill/MiSeq/LibraryK2/

usvi_temporal_k2 <- data.frame(indir = "/proj/omics/apprill/MiSeq/LibraryK2/",
                               read1 = sort(grep("plus", invert = TRUE, list.files(path = "/proj/omics/apprill/MiSeq/LibraryK2/", pattern = "R1_001", recursive = FALSE, include.dirs = TRUE), value = TRUE)),
                               read2 = sort(grep("plus",invert = TRUE, list.files(path = "/proj/omics/apprill/MiSeq/LibraryK2/", pattern = "R2_001", recursive = FALSE, include.dirs = TRUE), value = TRUE))) %>%
  dplyr::mutate(sample_name = gsub("(^[[:alnum:]]*)(_.*)", "\\1", read1)) %>%
  dplyr::mutate(read1 = paste0(.$indir, read1),
                read2 = paste0(.$indir, read2)) %>%
  dplyr::filter(sample_name %in% sample_metadata[["original_ID"]]) %>%
  droplevels
#rename the controls:

control_rename_lookup <- data.frame(sample_name = c("MockCommunity_Even_HM_782D", "neg_29cycles", "neg_30cy_PCR2",
                                                    "neg_30cy_PCR3", "neg_32cycles", "neg_34cycles", 
                                                    "Bowhead_3", "C-102"),
                                    relabel = c("mock_community", "pcr_control_01_29_cyc", "pcr_control_02_30_cyc", 
                                                "pcr_control_03_30_cyc", "pcr_control_04_32_cyc", "pcr_control_05_34_cyc",
                                                "bowhead_3", "control_102_jb")) %>%
  droplevels
# tibble::deframe(.)

usvi_temporal_w2 <- data.frame(indir = "/proj/omics/apprill/MiSeq/LibraryW2/",
                               read1 = sort(grep("R1_001", list.files(path = "/proj/omics/apprill/MiSeq/LibraryW2/", pattern = "R1_001", recursive = FALSE, include.dirs = TRUE), value = TRUE)),
                               read2 = sort(grep("R2_001", list.files(path = "/proj/omics/apprill/MiSeq/LibraryW2/", pattern = "R2_001", recursive = FALSE, include.dirs = TRUE), value = TRUE))) %>%
  dplyr::mutate(sample_name = gsub("((^[[:alnum:]_]*){1,})((_[A-Z]{6,8})(-)([A-Z]{6,8}.*))", "\\1", read1)) %>%
  dplyr::mutate(read1 = paste0(.$indir, read1),
                read2 = paste0(.$indir, read2)) %>%
  dplyr::left_join(., control_rename_lookup, by = join_by("sample_name")) %>%
  dplyr::mutate(sample_ID = dplyr::case_when(is.na(relabel) ~ sample_name,
                                             !is.na(relabel) ~ relabel,
                                             .default = sample_name)) %>%
  dplyr::select(-relabel) %>%
  droplevels

usvi_temporal_import <- bind_rows(usvi_temporal_w2, usvi_temporal_k2) %>%
  dplyr::mutate(sample_name = dplyr::case_when(is.na(sample_name) ~ sample_ID,
                                               .default =sample_name))

rm(usvi_temporal_w2)
rm(usvi_temporal_k2)

# Start DADA2 filtering ---------------------------------------------------


#look at quality profiles since that seems to be standard protocol
# 250bp PE reads

# dada2::plotQualityProfile(usvi_temporal_import$read1)
# #read1 quality generally stays above Q30 throughout the length of the read
# #in some samples, quality has more variability and slightly drops at ~cycle 60-70: Metab_K9, neg_34cycles, Metab_K10, either Metab_K5 or Metab_232


# dada2::plotQualityProfile(usvi_temporal_import$read2)
#plot just those for read2:
temp_quality_plot_df <- usvi_temporal_import %>%
  dplyr::slice_sample(n = 10) %>%
  bind_rows(., usvi_temporal_import %>%
              dplyr::filter(sample_name %in% c("neg_34cycles", "Metab_K9", "Metab_K5", "Metab_232", "Metab_K10"))) %>%
  dplyr::distinct(sample_name, .keep_all = TRUE) %>%
  droplevels

dada2::plotQualityProfile(temp_quality_plot_df$read1)
#from library W2, actual samples read1 drop in quality below Q30 as early as cycle 80
#negative controls read2 drop in quality below Q30 at cycle 80


dada2::plotQualityProfile(temp_quality_plot_df$read2)
#from library K2, sample 235 read2 has a drop in quality around 150
#from library W2, actual samples read2 drop in quality below Q30 as late as cycle 210 but as early as cycle 45
#negative controls read2 drop in quality below Q30 at cycle 45

#this is the 515-806 primer combo for the hypervariable v4 region of the 16S rRNA gene,so product would be ~250-300bp
#if we use both r1 and r2, r2 should be trimmed after 160 as per tutorial recs: https://benjjneb.github.io/dada2/tutorial.html
#trim read1 to 210bp (standard is at most 240bp)


# filter and trim
#output filtered to scratch for now: ~/scratch/projects/apprill/usvi/
usvi_qc_list <- usvi_temporal_import %>%
  dplyr::mutate(outdir = paste0("~/scratch/projects/apprill/usvi/", "filtered/"),
                filt1 = paste0(outdir, sample_name, "_R1_filt.fastq.gz"),
                filt2 = paste0(outdir, sample_name, "_R2_filt.fastq.gz"))

usvi_list <- usvi_qc_list %>%
  # head(n = 2) %>% #for testing out
  dplyr::select(sample_ID, read1, read2, filt1, filt2) %>%
  tidyr::pivot_longer(., cols = -c("sample_ID"),
                      names_to = "readtype",
                      values_to = "filepath") %>%
  split(., f = .$readtype) %>%
  map(., ~.x %>%
        dplyr::select(-c(readtype)) %>%
        droplevels %>%
        tidyr::unnest(., cols = filepath) %>%
        droplevels) %>%
  map(., ~.x[["filepath"]] %>%
        setNames(., .x[["sample_ID"]]))


usvi_qc_summary <- dada2::filterAndTrim(fwd = usvi_list$read1, 
                                        filt = usvi_list$filt1, 
                                        rev = usvi_list$read2,
                                        filt.rev = usvi_list$filt2,
                                        truncLen = c(210,160), 
                                        maxN = 0, 
                                        matchIDs = TRUE,
                                        multithread = nthreads,
                                        maxEE = c(2,2), 
                                        truncQ = 2,
                                        rm.phix = TRUE, 
                                        compress = TRUE)

##when first running filterAndTrim in a new directory, it will automatically make the directory for output:
## Creating output directory: /vortexfs1/home/sharon.grim/scratch/projects/apprill/usvi/filtered

usvi_qc_summary_df <- usvi_qc_summary %>%
  tibble::as_tibble(., rownames = "read1") %>%
  dplyr::mutate(sample_name = gsub("((^[[:alnum:]_]*){1,})((_[A-Z]{6,8})(-)([A-Z]{6,8}.*))", "\\1", read1)) %>%
  dplyr::select(-read1) %>%
  dplyr::left_join(., control_rename_lookup, by = join_by("sample_name")) %>%
  dplyr::left_join(., sample_metadata %>%
                     dplyr::select(sample_ID, original_ID),
                   by = join_by(sample_name == original_ID)) %>%
  dplyr::mutate(sample_ID = dplyr::case_when(is.na(sample_ID) & !is.na(relabel) ~ relabel,
                                             is.na(sample_ID) & is.na(relabel) ~ sample_name,
                                             .default = sample_ID)) %>%
  dplyr::select(sample_ID, sample_name, reads.in, reads.out) %>%
  droplevels

#error learning

usvi_qc_err_list <- c("filt1", "filt2") %>%
  map(., ~usvi_list[[.x]] %>%
        # simplify(.)) %>% #for testing out
        dada2::learnErrors(., multithread = nthreads)) %>%
  setNames(., c("err1", "err2"))

##with just two samples:
##56099610 total bases in 267141 reads from 2 samples will be used for learning the error rates.

#with all samples:
# 100609110 total bases in 479091 reads from 4 samples will be used for learning the error rates.
# 107766560 total bases in 673541 reads from 6 samples will be used for learning the error rates.

# dada2::plotErrors(usvi_qc_err_list$err1, nominalQ = TRUE)
# dada2::plotErrors(usvi_qc_err_list$err2, nominalQ = TRUE)
# #error rates of base substitution, decrease as consensus quality score increases

#sample inference

# as.list(c("filt1", "filt2"))
usvi_dada_list <- map2(c("filt1", "filt2"), c("err1", "err2"),
                       ~dada2::dada(derep = usvi_list[[.x]],
                                    err = usvi_qc_err_list[[.y]],
                                    multithread = nthreads)) %>%
  setNames(., c("dadas1", "dadas2"))

# Sample 1 - 154100 reads in 39822 unique sequences.
# Sample 2 - 113041 reads in 33496 unique sequences.
# Sample 3 - 106173 reads in 32302 unique sequences.
# Sample 4 - 105777 reads in 34737 unique sequences.
# Sample 5 - 107241 reads in 33617 unique sequences.
# Sample 6 - 87209 reads in 26467 unique sequences.
# Sample 7 - 126802 reads in 33962 unique sequences.
# Sample 8 - 96097 reads in 28993 unique sequences.
# Sample 9 - 69196 reads in 21766 unique sequences.
# Sample 10 - 90607 reads in 27630 unique sequences.
# Sample 11 - 119505 reads in 36684 unique sequences.
# Sample 12 - 118145 reads in 34818 unique sequences.
# Sample 13 - 88785 reads in 27110 unique sequences.
# Sample 14 - 106856 reads in 33002 unique sequences.
#...
# Sample 107 - 16491 reads in 8663 unique sequences.
# Sample 108 - 98483 reads in 36802 unique sequences.
# Sample 109 - 129348 reads in 42942 unique sequences.


usvi_dada_list <- usvi_dada_list %>%
  map(., ~.x %>%
        setNames(., names(usvi_list[["filt1"]])))

#inpsect the results:
# usvi_dada_list[["dadas2"]][1]
## generally for each sample, there are fewer reverse reads due to trimming and inference steps


#intermediate saving:

{
  temp_list <- list(usvi_list, usvi_qc_err_list,usvi_qc_list, usvi_dada_list, usvi_qc_summary, usvi_temporal_import) %>%
    setNames(., c("usvi_list", "usvi_qc_err_list","usvi_qc_list", "usvi_dada_list", "usvi_qc_summary", "usvi_temporal_import"))
  readr::write_rds(temp_list, paste0(projectpath, "/", "usvi_dada_res", ".rds"))
  rm(temp_list)
  }


#merge pairs of reads
usvi_dada_merged <- dada2::mergePairs(dadaF = usvi_dada_list[["dadas1"]],
                                      derepF = usvi_list[["filt1"]],
                                      dadaR = usvi_dada_list[["dadas2"]],
                                      derepR = usvi_list[["filt2"]],
                                      minOverlap = 12,
                                      verbose = TRUE)
#the verbose output indicates that over 95% of read pairs merged
#intermediate file saving:
if(!file.exists(paste0(projectpath, "/", "usvi_dada_merged", ".rds"))){
  readr::write_rds(usvi_dada_merged, paste0(projectpath, "/", "usvi_dada_merged", ".rds"))  
}

#construct sequence table of ASVs

usvi_asvs <- usvi_dada_merged %>%
  dada2::makeSequenceTable(., orderBy = "abundance")

# 
# temp_df <- tibble::enframe(rowSums(usvi_asvs), value = "reads.merged", name = "sample_name") %>%
#   dplyr::mutate(sample_ID = dplyr::case_when(grepl("_", sample_name) ~ sample_name,
#                                              # grepl("bowhead", sample_name) ~ stringr::str_to_title(sample_name),
#                                              .default = paste0("Metab_", sample_name))) %>%
#   dplyr::select(sample_ID, reads.merged)




#if you want to check the length distribution of the ASVs:
table(nchar(dada2::getSequences(usvi_asvs)))
#a majority of ASVs are 253bp in length.
#a small handful are <247bp, and >258bp

# retain only the seqs within that range.
usvi_asvs <- usvi_asvs[,nchar(colnames(usvi_asvs)) %in% 247:258]
table(nchar(dada2::getSequences(usvi_asvs)))

usvi_asvs_key <- usvi_asvs %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "sequence") %>%
  tibble::rowid_to_column(var = "rowid") %>%
  dplyr::mutate(rowid = sprintf("%05.f", rowid)) %>%
  dplyr::mutate(asv_id = paste0("ASV_", rowid)) %>%
  dplyr::select(asv_id, sequence) %>%
  droplevels
#chimera checking


usvi_nochim_asvs <- usvi_asvs %>%
  dada2::removeBimeraDenovo(., method = "consensus", multithread = nthreads, verbose = TRUE)
# Identified 318 bimeras out of 14300 input sequences.

sum(usvi_nochim_asvs)/sum(usvi_asvs)
#99% of sequences retained
# 
# temp_df <- tibble::enframe(rowSums(usvi_nochim_asvs), value = "reads.nochim", name = "sample_name") %>%
#     dplyr::mutate(sample_ID = dplyr::case_when(grepl("_", sample_name) ~ sample_name,
#                                                # grepl("bowhead", sample_name) ~ stringr::str_to_title(sample_name),
#                                                .default = paste0("Metab_", sample_name))) %>%
#     dplyr::select(sample_ID, reads.nochim)
# temp_df <- usvi_qc_summary_df %>%
usvi_qc_summary_df <- usvi_qc_summary_df %>%
  dplyr::left_join(., tibble::enframe(rowSums(usvi_asvs), value = "reads.merged", name = "sample_name") %>%
                     dplyr::mutate(sample_ID = dplyr::case_when(grepl("_", sample_name) ~ sample_name,
                                                                # grepl("bowhead", sample_name) ~ stringr::str_to_title(sample_name),
                                                                .default = paste0("Metab_", sample_name))) %>%
                     dplyr::select(sample_ID, reads.merged),
                   by = join_by(sample_ID)) %>%
  dplyr::left_join(., tibble::enframe(rowSums(usvi_nochim_asvs), value = "reads.nochim", name = "sample_name") %>%
                     dplyr::mutate(sample_ID = dplyr::case_when(grepl("_", sample_name) ~ sample_name,
                                                                # grepl("bowhead", sample_name) ~ stringr::str_to_title(sample_name),
                                                                .default = paste0("Metab_", sample_name))) %>%
                     dplyr::select(sample_ID, reads.nochim),
                   by = join_by(sample_ID))

#13982 ASVs
usvi_nochim_asvs_key <- usvi_nochim_asvs %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "sequence") %>%
  dplyr::left_join(., usvi_asvs_key, by = join_by(sequence)) %>% #since we already assigned ASV identifiers before chimera checking...
  # tibble::rowid_to_column(var = "rowid") %>%
  # dplyr::mutate(rowid = sprintf("%05.f", rowid)) %>%
  # dplyr::mutate(asv_id = paste0("ASV_", rowid)) %>%
  dplyr::select(asv_id, sequence) %>%
  droplevels

usvi_dada_nochim <- list(usvi_asvs, usvi_asvs_key, usvi_nochim_asvs, usvi_nochim_asvs_key) %>%
  setNames(., c("usvi_asvs", "usvi_asvs_key", "usvi_nochim_asvs", "usvi_nochim_asvs_key"))
if(!file.exists(paste0(projectpath, "/", "usvi_dada_nochim", ".rds"))){
  readr::write_rds(usvi_dada_nochim, file = paste0(projectpath, "/", "usvi_dada_nochim", ".rds"))
}



if(!file.exists(paste0(projectpath, "/", "usvi_qc_summary_df", ".tsv"))){
  readr::write_delim(usvi_qc_summary_df, paste0(projectpath, "/", "usvi_qc_summary_df", ".tsv"),
                     delim  = "\t", col_names = TRUE, num_threads = nthreads)
}