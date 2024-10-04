# 00_process_env.R

# Load packages -----------------------------------------------------------

if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
  devtools::install_github("dempsey-CMAR/sensorstrings")
}
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

library(sensorstrings)

library(tidyverse)
library(data.table)
library(BiocManager)
library(BiocParallel)
# library(phyloseq)
library(viridis)
library(cli)
library(furrr)
library(progressr)
library(patchwork)
library(viridisLite)
library(pals)

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

# Prepare lookup vectors --------------------------------------------------

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
                                   label = c("Dawn", "Peak photosynthesis")) %>%
  tibble::deframe(.)
sampling_time_colors <- pals::ocean.haline(n = length(sampling_time_lookup)) %>%
  setNames(., names(sampling_time_lookup))
sampling_day_lookup <- data.frame(sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5"),
                                  label = c("20210122", "20210123", "20210124", "20210125", "20210126")) %>%
  tibble::deframe(.)
sampling_day_colors <- pals::ocean.thermal(n = length(sampling_day_lookup)) %>%
  setNames(., names(sampling_day_lookup))


# Import the tide gauge record from Lameshur Bay --------------------------

noaa_lb_tide_df <- readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "TideElev_LameshurBay_Jan2021_hourly_m.txt"),
                             delim = "\t", quote = "", skip = 14, comment = "#", show_col_types = FALSE, col_names = FALSE, col_select = c(1:4)) %>%
  setNames(., c("ymd", "day", "time", "height")) %>%
  dplyr::mutate(fulldate = paste0(ymd, " ", time) %>%
                  lubridate::ymd_hms(., tz = "America/Virgin")) %>%
  dplyr::mutate(site = "LB_seagrass") %>%
  droplevels

# Import HOBOdata ---------------------------------------------------------

#time inthe HOBO output looks like: 01/21/21 04:32:30 PM

hobo_df <- readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Tektite_Temporal_Study_SJ_Jan21.csv"),
                             delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
  setNames(., c("order", "date", "temp", "lumens")) %>%
  dplyr::mutate(separadate = lubridate::parse_date_time(date, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
  # dplyr::mutate(separadate = lubridate::parse_date_time(date, c("%m/%d/%y %I:%M:%S Op", "mdy"))) %>%
  # dplyr::mutate(separadate = lubridate::ymd_hms(separadate, tz = "UTC")) %>%
  # tidyr::separate_wider_delim(., separadate, names = c("day", "time"), delim = " ", too_few = "align_start", too_many = "merge", cols_remove = FALSE) %>%
  droplevels %>%
  dplyr::mutate(site = "Tektite") %>%
  bind_rows(., (readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Yawzi_Temporal_Study_SJ_Jan21.csv"),
                                  delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
                  setNames(., c("order", "date", "temp", "lumens")) %>%
                  dplyr::mutate(separadate = lubridate::parse_date_time(date, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
                  # tidyr::separate_wider_delim(., separadate, names = c("day", "time"), delim = " ", too_few = "align_start", too_many = "merge", cols_remove = FALSE) %>%
                  droplevels %>%
                  dplyr::mutate(site = "Yawzi")),
            (readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Seagrass_Temporal_Study_SJ_Jan21.csv"),
                              delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
              setNames(., c("order", "date", "temp", "lumens")) %>%
               dplyr::mutate(separadate = lubridate::parse_date_time(date, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
              # tidyr::separate_wider_delim(., separadate, names = c("day", "time"), delim = " ", too_few = "align_start", too_many = "merge", cols_remove = FALSE) %>%
              droplevels %>%
              dplyr::mutate(site = "LB_seagrass"))) %>%
  dplyr::mutate(separadate = lubridate::with_tz(separadate, tz = "America/Virgin")) %>%
  tidyr::separate_wider_delim(., separadate, names = c("day", "time"), delim = " ", too_few = "align_start", too_many = "merge", cols_remove = FALSE)



# Plot tide data from LB --------------------------------------------------
noaa_lb_tide_filtered_df <- noaa_lb_tide_df %>%
  dplyr::filter(grepl(paste0(c("2021-01-22", "2021-01-23", "2021-01-24", "2021-01-25", "2021-01-26"), collapse = "|"), ymd)) %>%
  droplevels

g <- print(
  ggplot(data = noaa_lb_tide_filtered_df,
         aes(x = fulldate, y = height, fill = site))
  + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

# Plot HOBO temp data -----------------------------------------------------


#ALl of the hobo time seems to be in GMT, given that lumens begins to increase around 11:00:00 each day which is ~7:00am in USVI
#so offset the times by 4hours:

hobo_filtered_df <- hobo_df %>%
  # dplyr::mutate(separadate = lubridate::with_tz(separadate, tz = "America/Virgin")) %>%
  # dplyr::mutate(separadate = lubridate::with_tz(separadate, tz = "GMT")) %>%
  # dplyr::mutate(separadate = lubridate::ymd_hms(separadate, tz = "UTC")) %>%
  # dplyr::mutate(separadate = lubridate::force_tz(separadate, tzones = "UTC", tzone_out = "America/New York")) %>%
  dplyr::filter(day %in% c("2021-01-22", "2021-01-23", "2021-01-24", "2021-01-25", "2021-01-26")) %>%
  dplyr::arrange(order) %>%
  dplyr::group_by(site) %>%
  dplyr::slice_head(prop = 0.9) %>%
  droplevels 

start_time <- hobo_filtered_df %>%
  dplyr::ungroup(.) %>%
  dplyr::select(separadate) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::distinct(., .keep_all = FALSE) %>%
  dplyr::mutate(separadate = lubridate::parse_date_time(separadate, c("ymd HMS"))) %>%
  droplevels

g <- print(
  ggplot(data = hobo_filtered_df %>%
           dplyr::filter(grepl("LB", site)),
         aes(x = separadate, y = temp, fill = site))
  + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)


# Calculate changes in temp and lumens ------------------------------------

temp_list <- hobo_filtered_df %>%
  split(., f = .$site) %>%
  purrr::map(., ~.x %>%
               # dplyr::filter(site == "LB_seagrass") %>%
               dplyr::ungroup(.) %>%
               dplyr::select(order, separadate, lumens, temp) %>%
               droplevels)
temp_list2 <- temp_list %>%
  purrr::map(., ~bind_rows(.x %>%
                      dplyr::slice(., -(1:1)),
                    (tibble::tribble(
                      ~order, ~separadate, ~lumens, ~temp,
                      max(hobo_df[["order"]]), NA, NA, NA))) %>%
        dplyr::rename_with(., ~paste0(.x, ".next"), everything()))

hobo_filtered_delta_df <- names(temp_list2) %>%
  map(., ~temp_list[[.x]] %>%
        bind_cols(., temp_list2[[.x]])) %>%
  setNames(., names(temp_list2)) %>%
  imap(., ~.x %>%
         dplyr::mutate(site = .y)) %>%
  purrr::list_flatten(., name_spec = "{inner}") %>%
  map(., ~.x %>%
        dplyr::rowwise(.) %>%
  dplyr::mutate(interval = lubridate::time_length(interval(separadate, separadate.next), "minute"),
                del_temp = (temp.next - temp),
                del_lumens = (lumens.next - lumens)) %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(time_lapsed = purrr::accumulate(interval, `+`)) %>%
  droplevels) %>%
  bind_rows(., .id = NULL)

g2 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           droplevels,
         aes(x = time_lapsed, y = temp - 26.8, color = site))
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_color_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

#null any temp entries where abs(del_lumens) > 5

# temp_df2 <- hobo_filtered_df %>%
#   dplyr::filter(site == "LB_seagrass") %>%
#   dplyr::ungroup(.) %>%
#   dplyr::select(order, separadate, lumens, temp) %>%
#   droplevels %>%
#   bind_cols(., bind_rows(hobo_filtered_df %>%
#                   dplyr::filter(site == "LB_seagrass") %>%
#                   dplyr::ungroup(.) %>%
#                   dplyr::select(order, separadate, lumens, temp) %>%
#                   droplevels %>%
#                   dplyr::slice(., -(1:1)),
#                   (tibble::tribble(
#                     ~order, ~separadate, ~lumens, ~temp,
#                     max(hobo_df[["order"]]), NA, NA, NA))) %>%
#       dplyr::rename_with(., ~paste0(.x, ".next"), everything())) %>%
#   dplyr::mutate(site = "LB_seagrass") %>%
#   dplyr::rowwise(.) %>%
#   dplyr::mutate(interval = lubridate::time_length(interval(separadate, separadate.next), "minute"),
#                 del_temp = (temp.next - temp),
#                 del_lumens = (lumens.next - lumens)) %>%
#   dplyr::group_by(site) %>%
#   dplyr::mutate(time_lapsed = purrr::accumulate(interval, `+`)) %>%
#   droplevels

g3 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           droplevels)
  + geom_point(aes(x = time_lapsed, y = del_temp, fill = site), shape = 22)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

g3 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           droplevels)
  + geom_point(aes(x = time_lapsed, y = del_lumens, fill = site), shape = 23)
  # + geom_point(aes(x = time_lapsed, y = del_lumens/5, fill = site), shape = 23)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

#now pair with the NOAA tide table predictions
g4 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels, aes(x = separadate, y = temp, color = site))
# + geom_point(shape = 21)
+ geom_line(show.legend = FALSE)
+ facet_wrap(.~site, labeller = labeller(site = site_lookup))
+ scale_color_manual(values = site_colors)
+ theme(axis.text.x = element_text(angle = 90))
+ geom_line(data = noaa_lb_tide_filtered_df, aes(x = fulldate, y = height + 27), color = "cyan")
+ scale_y_continuous(name = "Temperature ˚C", 
                     sec.axis = dup_axis(~. - 27, 
                                         name = "Predicted tide height (m)"))
+ scale_x_datetime("Date", date_labels = "%y%m%d %I%M")
)

g5 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels, aes(x = separadate, y = lumens, color = site))
  # + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_color_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = fulldate, y = height *4000), color = "cyan")
  + scale_y_continuous(name = "Lumens", 
                       sec.axis = dup_axis(~. /4000, 
                                           name = "Predicted tide height (m)"))
  + scale_x_datetime("Date", date_labels = "%y%m%d %I%M")
)

temp_df3 <- hobo_filtered_delta_df %>%
  dplyr::mutate(del_temp = dplyr::case_when(abs(del_lumens) < 100 ~ del_temp,
                                            .default = NA)) %>%
  dplyr::ungroup(.) %>%
  # dplyr::mutate(del_temp = dplyr::case_when(abs(del_temp) < 0.1 ~ NA,
  #                                           .default = del_temp)) %>%
  dplyr::mutate(coswave = cos((time_lapsed+200)*pi/750)*0.25) %>%
  dplyr::mutate(sinewave = sin((time_lapsed+1000)*pi/750)*0.25) %>%
  droplevels

# g3 <- print(
#   ggplot(data = temp_df3,
#          aes(x = time_lapsed, y = del_temp, fill = site))
#          # aes(x = separadate.next, y = del_temp, fill = site))
#   + geom_point(shape = 21)
#   # + geom_line(show.legend = FALSE)
#   # + geom_smooth(formula = y~sin(x*pi/3600),
#   #               method = "lm")
#   + facet_wrap(.~site, labeller = labeller(site = site_lookup))
#   + scale_fill_manual(values = site_colors)
#   + theme(axis.text.x = element_text(angle = 90))
#   + geom_line(data = temp_df3, aes(x = time_lapsed, y = sinewave), shape = 22)
# )
g3 <- print(
  ggplot(data = temp_df3 %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels,
         aes(x = time_lapsed, y = temp, fill = site))
  + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + geom_point(data = temp_df3, aes(x = time_lapsed, y = coswave+26.8), shape = 23)
  + theme(axis.text.x = element_text(angle = 90))
)

g3 <- print(
  ggplot(data = temp_df3 %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels)
  + geom_point(shape = 21, aes(x = separadate, y = temp, fill = site))
  + geom_line(show.legend = FALSE, aes(x = separadate, y = temp))
  + geom_point(data = temp_df3, aes(x = separadate.next, y = coswave+26.8), shape = 23)
  + theme(axis.text.x = element_text(angle = 90))
  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = fulldate, y = height +26.8 ), color = "cyan")
  + scale_y_continuous(name = "Temperature ˚C", 
                       sec.axis = dup_axis(~. - 26.8, 
                                           name = "Predicted tide height (m)"))
  # + scale_x_datetime("Date", date_labels = "%y%m%d %I%M")
)




# Play with sine and cosine -----------------------------------------------


temp_df4 <- hobo_filtered_delta_df %>%
  dplyr::mutate(del_temp = dplyr::case_when(abs(del_lumens) < 100 ~ del_temp,
                                            .default = NA)) %>%
  dplyr::ungroup(.) %>%
  # dplyr::mutate(del_temp = dplyr::case_when(abs(del_temp) < 0.1 ~ NA,
  #                                           .default = del_temp)) %>%
  dplyr::mutate(coswave = cos((time_lapsed+200)*pi/725)*0.25) %>%
  # dplyr::mutate(sinewave = sin((time_lapsed+1000)*pi/725)*1) %>%
  dplyr::mutate(sinewave = sin((time_lapsed+1000)*pi/725)*0.75) %>%
  dplyr::mutate(coswave = coswave + 26.8,
                sinewave = sinewave + 27.2) %>%
  # dplyr::mutate(combinewave = (sinewave/ coswave + 26.8)) %>%
  droplevels

print(
  ggplot(data = temp_df4 %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels)
  + geom_line(shape = 21, aes(x = separadate, y = temp), color = "black")
  + geom_line(aes(x = separadate, y = coswave), color = "purple") #this is the little peak
  + geom_line(aes(x = separadate, y = sinewave), color = "orange") #this is  the big peak 
  # + geom_line(aes(x = separadate, y = combinewave), color = "blue") #this is  the big peak
  + theme(axis.text.x = element_text(angle = 90))
  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = fulldate, y = height +26.8 ), color = "cyan")
  + scale_y_continuous(name = "Temperature ˚C", 
                       sec.axis = dup_axis(~. - 26.8, 
                                           name = "Predicted tide height (m)"))
  # + scale_x_datetime("Date", date_labels = "%y%m%d %I%M")
)
