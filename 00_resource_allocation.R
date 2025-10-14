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
if (!require("svglite", quietly = TRUE)){
  install.packages("svglite")
}
library(BiocManager)
library(BiocParallel)
library(data.table)
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
    if(ANSWER == "getwd()"){
      assign("projectpath", getwd(), inherits = TRUE, envir = .GlobalEnv)  
    } else {
      stop("Please provide a path.\n")   
    }
  }
}

# if(interactive()) try(f_projectpath())
# try(f_projectpath())

# User input to proceed with code chunk -----------------------------------


f_run_chunk <- function() {
  #run 'try(f_run_chunk())' to set variable 'execchunk' to TRUE or FALSE
  #after running this function, use 'if(execchunk){}' to run a code chunk within the curly braces
  ANSWER <- readline("Do you want to run this chunk? ")
  if (substr(ANSWER, 1, 1) != "y"){
    assign("execchunk", FALSE, envir = .GlobalEnv, inherits = TRUE)
    stop("Skipping the code chunk.")
  } else {
    #do the code chunk
    cli::cli_alert_info("Running the code chunk... \n") 
    assign("execchunk", TRUE, envir = .GlobalEnv, inherits = TRUE)
  }
}