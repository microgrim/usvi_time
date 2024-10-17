# 00_resource_allocation.R

#prepare resources for R on your computer

# Load packages -----------------------------------------------------------

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("furrr", quietly = TRUE)){
  install.packages("furrr")
  install.packages("future")
}
library(data.table)
library(BiocManager)
library(BiocParallel)
library(cli)
library(furrr)
library(progressr)

# Plan for resource allocation --------------------------------------------
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
  register(bpparam_multi)
} else {
  bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                workers = nthreads,
                                                stop.on.error = TRUE,
                                                RNGseed = 48105,
                                                progressbar = TRUE)
  register(bpparam_multi)
}


f_projectpath <- function() {
  ANSWER <- readline("Type in your full project path. ")
  if (substr(ANSWER, 1, 1) == "/"){
    if(file.exists(ANSWER)){
      cli::cli_alert_info("Setting path to {ANSWER} \n") 
      assign("projectpath", ANSWER, inherits = TRUE, envir = .GlobalEnv)  
    } else
      stop("Your provided path does not seem to exist.")
  } else {
    stop("Please provide a path.\n") 
  }
}

# if(interactive()) try(f_projectpath())
# try(f_projectpath())


