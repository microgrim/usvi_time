# 00_process_env.R

# Resource allocation time ------------------------------------------------

if(file.exists(paste0(getwd(), "/", "00_resource_allocation.R"))){
  cat("Preparing resource allocations.")
  source(paste0(getwd(), "/", "00_resource_allocation.R"), local = FALSE,
         echo = TRUE, verbose = getOption("verbose"), prompt.echo = getOption("prompt"))
  try(f_projectpath())
} else {
  cat("Preparing resource allocations.")
  
  #load packages
  # if (!require("devtools", quietly = TRUE)){
  #   install.packages("devtools")
  #   devtools::install_github("dempsey-CMAR/sensorstrings")
  # }
  # library(sensorstrings)
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
library(viridis)
library(patchwork)
library(viridisLite)
library(pals)

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
  dplyr::mutate(date_ast = paste0(ymd, " ", time) %>%
                  lubridate::ymd_hms(., tz = "America/Virgin")) %>%
  dplyr::mutate(site = "LB_seagrass") %>%
  droplevels

# Import HOBOdata ---------------------------------------------------------

#time inthe HOBO output looks like: 01/21/21 04:32:30 PM
#ALl of the hobo time seems to be in GMT, given that lumens begins to increase around 11:00:00 each day which is ~7:00am in USVI
#so offset the times by 4hours:

#it looks like at MIS, we used HOBOs deployed at multiple depths to 23m below surface to measure light
#Hobos record in Lux or Lumens, which is not easily translated to PAR
#to get the full light field, we also used multiple light profilers: Hyperspectral, LiCor, C-OPs
#the CONVERSION of the reported Hobo light in Lux, to uE was as follows:

#at 23m depth, only blue-green PAR was available as determined by the Blackbird hyperspectral profiler. 
#in July 2016, Blackbird measured ~75uE of available light almost entirely in Blue-Green wavelengths
#HOBOs measured light at 23m depth on those same days as ~996 lux
#so the light-profile/depth specific conversion was:
#Lux/54/0.245142572 = available PAR at 23m depth
#in this case, 996/54/.245142572 = 75 uE

#in MIS the attenuation coefficients averaged -0.132840405 (range: -0.161106824, -0.101399784)


hobo_df <- readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Tektite_Temporal_Study_SJ_Jan21_0-lux.csv"),
                             delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
  setNames(., c("order", "date_ast", "temp", "lux")) %>%
  # dplyr::mutate(date_ast = lubridate::parse_date_time(date_gmt, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "AST")) %>%
  droplevels %>%
  dplyr::mutate(site = "Tektite") %>%
  bind_rows(., (readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Yawzi_Temporal_Study_SJ_Jan21_0-lux.csv"),
                                  delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
                  setNames(., c("order", "date_ast", "temp", "lux")) %>%
                  # dplyr::mutate(date_ast = lubridate::parse_date_time(date_gmt, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "AST")) %>%
                  droplevels %>%
                  dplyr::mutate(site = "Yawzi")),
            (readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Seagrass_Temporal_Study_SJ_Jan21_0-lux.csv"),
                               delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
               setNames(., c("order", "date_ast", "temp", "lux")) %>%
               # dplyr::mutate(date_ast = lubridate::parse_date_time(date_gmt, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "AST")) %>%
               droplevels %>%
               dplyr::mutate(site = "LB_seagrass"))) %>%
  dplyr::mutate(date_ast = lubridate::parse_date_time(date_ast, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
  # dplyr::mutate(date_ast = lubridate::with_tz(date_ast, tz = "America/Virgin")) %>%
  tidyr::separate_wider_delim(., date_ast, names = c("day", "time"), delim = " ", too_few = "align_start", too_many = "merge", cols_remove = FALSE)

hobo_df2 <- readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Tektite_Temporal_Study_SJ_Jan21.csv"),
                             delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
  setNames(., c("order", "date_gmt", "temp", "lumens")) %>%
  dplyr::mutate(date_ast = lubridate::parse_date_time(date_gmt, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
  droplevels %>%
  dplyr::mutate(site = "Tektite") %>%
  bind_rows(., (readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Yawzi_Temporal_Study_SJ_Jan21.csv"),
                                  delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
                  setNames(., c("order", "date_gmt", "temp", "lumens")) %>%
                  dplyr::mutate(date_ast = lubridate::parse_date_time(date_gmt, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
                  droplevels %>%
                  dplyr::mutate(site = "Yawzi")),
            (readr::read_delim(paste0("/Users/sharongrim/projects/apprill/usvi_temporal/", "Seagrass_Temporal_Study_SJ_Jan21.csv"),
                              delim = ",", quote = "", skip = 2, comment = "#", show_col_types = FALSE, col_names = FALSE) %>%
              setNames(., c("order", "date_gmt", "temp", "lumens")) %>%
               dplyr::mutate(date_ast = lubridate::parse_date_time(date_gmt, c("%m/%d/%y %I:%M:%S Op", "mdy"), tz = "GMT")) %>%
              droplevels %>%
              dplyr::mutate(site = "LB_seagrass"))) %>%
  dplyr::mutate(date_ast = lubridate::with_tz(date_ast, tz = "America/Virgin")) %>%
  tidyr::separate_wider_delim(., date_ast, names = c("day", "time"), delim = " ", too_few = "align_start", too_many = "merge", cols_remove = FALSE)

hobo_df <- hobo_df %>%
  dplyr::left_join(., hobo_df2 %>% dplyr::select(order, site, lumens, temp))


#to look at intervals of time:
# interval <- lubridate::as.interval(duration(hours = 4), ymd_hms("2021-01-21 11:00:00", tz = "America/Virgin")) 
# interval_vector_peak <- list(interval,
#                              lubridate::int_shift(interval, days(1)),
#                              lubridate::int_shift(interval, days(2)),
#                              lubridate::int_shift(interval, days(3)),
#                              lubridate::int_shift(interval, days(4)),
#                              lubridate::int_shift(interval, days(5)))
# 
# interval <- lubridate::as.interval(duration(hours = 4), ymd_hms("2021-01-21 04:00:00", tz = "America/Virgin")) 
# interval_vector_dawn <- list(interval,
#                              lubridate::int_shift(interval, days(1)),
#                              lubridate::int_shift(interval, days(2)),
#                              lubridate::int_shift(interval, days(3)),
#                              lubridate::int_shift(interval, days(4)),
#                              lubridate::int_shift(interval, days(5)))

if(!file.exists(paste0(projectpath, "/", "usvi_hobo_light_temp", ".tsv"))){
  readr::write_delim(hobo_df, paste0(projectpath, "/", "usvi_hobo_light_temp", ".tsv", ".gz"), delim = "\t", col_names = TRUE)
}

# Plot tide data from LB --------------------------------------------------
noaa_lb_tide_filtered_df <- noaa_lb_tide_df %>%
  dplyr::filter(grepl(paste0(c("2021-01-22", "2021-01-23", "2021-01-24", "2021-01-25", "2021-01-26"), collapse = "|"), ymd)) %>%
  droplevels

g1 <- (
  ggplot(data = noaa_lb_tide_filtered_df,
         aes(x = date_ast, y = height, fill = site))
  + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)
if(interactive()){
  print(g1)
}

# Plot HOBO temp data -----------------------------------------------------

hobo_filtered_df <- hobo_df %>%
  dplyr::filter(day %in% c("2021-01-22", "2021-01-23", "2021-01-24", "2021-01-25", "2021-01-26")) %>%
  dplyr::arrange(order) %>%
  dplyr::group_by(site) %>%
  dplyr::slice_head(prop = 0.9) %>%
  droplevels 

# start_time <- hobo_filtered_df %>%
#   dplyr::ungroup(.) %>%
#   dplyr::select(date_ast) %>%
#   dplyr::slice_head(n = 1) %>%
#   dplyr::distinct(., .keep_all = FALSE) %>%
#   dplyr::mutate(date_ast = lubridate::parse_date_time(date_ast, c("ymd HMS"))) %>%
#   droplevels

g2 <- (
  ggplot(data = hobo_filtered_df %>%
           dplyr::filter(grepl("LB", site)),
         aes(x = date_ast, y = temp, fill = site))
  + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)
if(interactive()){
  print(g2)
}

# Calculate changes in temp and lumens ------------------------------------

temp_list <- hobo_filtered_df %>%
  split(., f = .$site) %>%
  purrr::map(., ~.x %>%
               # dplyr::filter(site == "LB_seagrass") %>%
               dplyr::ungroup(.) %>%
               dplyr::select(order, date_ast, lumens, temp) %>%
               droplevels)
temp_list2 <- temp_list %>%
  purrr::map(., ~bind_rows(.x %>%
                      dplyr::slice(., -(1:1)),
                    (tibble::tribble(
                      ~order, ~date_ast, ~lumens, ~temp,
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
  dplyr::mutate(interval = lubridate::time_length(interval(date_ast, date_ast.next), "minute"),
                del_temp = (temp.next - temp),
                del_lumens = (lumens.next - lumens)) %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(time_lapsed = purrr::accumulate(interval, `+`)) %>%
  droplevels) %>%
  bind_rows(., .id = NULL)

g3a <- (
  ggplot(data = hobo_filtered_delta_df %>%
           droplevels,
         aes(x = time_lapsed, y = temp - 26.8, color = site))
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_color_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)
g3b <- (
  ggplot(data = hobo_filtered_delta_df %>%
           droplevels)
  + geom_point(aes(x = time_lapsed, y = del_temp, fill = site), shape = 22)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

g3c <- (
  ggplot(data = hobo_filtered_delta_df %>%
           droplevels)
  + geom_point(aes(x = time_lapsed, y = del_lumens, fill = site), shape = 23)
  # + geom_point(aes(x = time_lapsed, y = del_lumens/5, fill = site), shape = 23)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_fill_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
)

if(interactive()){
  g3 <- g3a / g3b / g3c + patchwork::plot_layout(guides = "collect")
  print(g3)
}
#null any temp entries where abs(del_lumens) > 5

# temp_df2 <- hobo_filtered_df %>%
#   dplyr::filter(site == "LB_seagrass") %>%
#   dplyr::ungroup(.) %>%
#   dplyr::select(order, date_ast, lumens, temp) %>%
#   droplevels %>%
#   bind_cols(., bind_rows(hobo_filtered_df %>%
#                   dplyr::filter(site == "LB_seagrass") %>%
#                   dplyr::ungroup(.) %>%
#                   dplyr::select(order, date_ast, lumens, temp) %>%
#                   droplevels %>%
#                   dplyr::slice(., -(1:1)),
#                   (tibble::tribble(
#                     ~order, ~date_ast, ~lumens, ~temp,
#                     max(hobo_df[["order"]]), NA, NA, NA))) %>%
#       dplyr::rename_with(., ~paste0(.x, ".next"), everything())) %>%
#   dplyr::mutate(site = "LB_seagrass") %>%
#   dplyr::rowwise(.) %>%
#   dplyr::mutate(interval = lubridate::time_length(interval(date_ast, date_ast.next), "minute"),
#                 del_temp = (temp.next - temp),
#                 del_lumens = (lumens.next - lumens)) %>%
#   dplyr::group_by(site) %>%
#   dplyr::mutate(time_lapsed = purrr::accumulate(interval, `+`)) %>%
#   droplevels


# Tie together the site data with NOAA tide table -------------------------



#now pair with the NOAA tide table predictions
g4 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels, aes(x = date_ast, y = temp, color = site))
# + geom_point(shape = 21)
+ geom_line(show.legend = FALSE)
+ facet_wrap(.~site, labeller = labeller(site = site_lookup))
+ scale_color_manual(values = site_colors)
+ theme(axis.text.x = element_text(angle = 90))
+ geom_line(data = noaa_lb_tide_filtered_df, aes(x = date_ast, y = height + 27), color = "cyan")
+ scale_y_continuous(name = "Temperature ˚C", 
                     sec.axis = dup_axis(~. - 27, 
                                         name = "Predicted tide height (m)"))
+ scale_x_datetime("Date", date_labels = "%y%m%d %I%M")
)

g5 <- print(
  ggplot(data = hobo_filtered_delta_df %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels, aes(x = date_ast, y = lumens, color = site))
  # + geom_point(shape = 21)
  + geom_line(show.legend = FALSE)
  + facet_wrap(.~site, labeller = labeller(site = site_lookup))
  + scale_color_manual(values = site_colors)
  + theme(axis.text.x = element_text(angle = 90))
  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = date_ast, y = height *4000), color = "cyan")
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
#          # aes(x = date_ast.next, y = del_temp, fill = site))
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
  + geom_point(shape = 21, aes(x = date_ast, y = temp, fill = site))
  + geom_line(show.legend = FALSE, aes(x = date_ast, y = temp))
  + geom_point(data = temp_df3, aes(x = date_ast.next, y = coswave+26.8), shape = 23)
  + theme(axis.text.x = element_text(angle = 90))
  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = date_ast, y = height +26.8 ), color = "cyan")
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
  dplyr::mutate(sinewave = sin((time_lapsed+1000)*pi/725)*1.1) %>% #this is the big peak
  dplyr::mutate(coswave = cos((time_lapsed+400)*pi/725)*0.22) %>% #this is the little peak
  # dplyr::mutate(coswave = coswave + 26.9,
  #               sinewave = sinewave + 26.8) %>%
  dplyr::mutate(tide_sinewave = sin((time_lapsed+00)*pi/740)*0.2) %>% #this is the tide's big peak
  dplyr::mutate(tide_coswave = cos((time_lapsed+700)*pi/725)*0.09) %>% #this is the tide's little peak
  dplyr::mutate(tide_shift = sin((time_lapsed-200)*pi/12600)*0.2+0.92) %>%
  droplevels

print(
  ggplot(data = temp_df4 %>%
           dplyr::filter(grepl("LB", site)) %>%
           droplevels)
  + geom_line(aes(x = date_ast, y = temp), color = "black")
  + geom_line(aes(x = date_ast, y = lumens/4000+26.8), color = "yellow")
  + geom_line(aes(x = date_ast, y = sinewave + 26.8), color = "tan", linetype = 2) #this is  the big peak
  + geom_line(aes(x = date_ast, y = coswave + 26.9), color = "magenta", linetype = 2) #this is the little peak

  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = date_ast, y = height +26.8 ), color = "cyan", linewidth = 3)
  + scale_y_continuous(name = "Temperature ˚C",
                       sec.axis = dup_axis(~. - 26.8,
                                           name = "Predicted tide height (m)"))
  + geom_line(aes(x = date_ast, y = tide_sinewave*tide_shift + 26.8), color = "orange", linetype = 2) #this is  the big peak
  + geom_line(aes(x = date_ast, y = tide_coswave/tide_shift+ 26.8), color = "purple", linetype = 2) #this is the little peak
  # + geom_line(aes(x = date_ast, y = tide_shift), color = "salmon") #this is the slight shift in tide peak
)

#it looks like on these 4 days, every day between noon-one there is an inflection in the tide height
#also a local dip between the lumens max 
#also a dip in temp

#in the modeled regression lines of Lameshur Bay
#where the secondary (magenta)Temperature profile intersects the primary (orange) tide profile
#seems to precede the dip in light and temperature

#so solve for the intersections of coswave and (tide_sinewave*tide_shift + 26.8)
#cos((time_lapsed+400)*pi/725)*0.22 + 26.9
#(sin((time_lapsed+00)*pi/740)*0.2)*(sin((time_lapsed-200)*pi/12600)*0.2+0.92) + 26.8


#what do the other sites look like?
print(
  ggplot(data = temp_df4 %>%
           dplyr::filter(grepl("Tektite", site)) %>%
           droplevels)
  + geom_line(aes(x = date_ast, y = temp), color = "black")
  + geom_line(aes(x = date_ast, y = lumens/4000+26.8), color = "yellow")
  + geom_line(aes(x = date_ast, y = sinewave/1.7 + 26.8), color = "tan", linetype = 2) #this is  the big peak
  + geom_line(aes(x = date_ast, y = coswave/1.7 + 26.9), color = "magenta", linetype = 2) #this is the little peak
  
  + geom_line(data = noaa_lb_tide_filtered_df, aes(x = date_ast, y = height +26.8 ), color = "cyan", linewidth = 3)
  + scale_y_continuous(name = "Temperature ˚C",
                       sec.axis = dup_axis(~. - 26.8,
                                           name = "Predicted tide height (m)"))
  # + geom_line(aes(x = date_ast, y = tide_sinewave*tide_shift + 26.8), color = "orange", linetype = 2) #this is  the big peak
  # + geom_line(aes(x = date_ast, y = tide_coswave/tide_shift+ 26.8), color = "purple", linetype = 2) #this is the little peak
  # + geom_line(aes(x = date_ast, y = tide_shift), color = "salmon") #this is the slight shift in tide peak
)
