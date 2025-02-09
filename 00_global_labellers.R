## 00_global_labellers.R

#project-wide renamers and lookups



model_dispersion_lookup <- data.frame(v1 = c("manual_dispersion", 
                                             "manual_dawn", "manual_photo", "auto_dawn", "auto_photo",
                                             "auto_dispersion"),
                                      label = c("All LB seagrass, manual dispersion",
                                                "Dawn LB seagrass, manual dispersion",
                                                "Peak photo LB seagrass, manual dispersion",
                                                "Dawn LB seagrass, auto dispersion",
                                                "Peak photo LB seagrass, auto dispersion",
                                                "All LB seagrass, auto dispersion")) %>%
  bind_rows(., data.frame(v1 = c("manual_lameshur", "manual_yawzi", "manual_tektite", "manual_reef",
                                 "auto_lameshur", "auto_yawzi",
                                 "auto_tektite", "auto_reef"),
                          label = c("Dawn LB seagrass, manual dispersion", "Dawn Yawzi, manual dispersion", "Dawn Tektite, manual dispersion", "Tektite, manual dispersion",
                                    "Dawn LB seagrass, auto dispersion", "Dawn Yawzi, auto dispersion", "Dawn Tektite, auto dispersion", "Tektite, auto dispersion"))) %>%
  tibble::deframe(.)


group_labels_lookup <- c("Tektite_d" = "Tektite dawn",
                         "Tektite_p" = "Tektite afternoon",
                         "Yawzi_d" = "Yawzi dawn",
                         "Yawzi_p" = "Yawzi afternoon",
                         "LB_seagrass_d" = "Lameshur Bay seagrass dawn",
                         "LB_seagrass_p" = "Lameshur Bay seagrass afternoon",
                         "Tektite" = "Tektite",
                         "Yawzi" = "Yawzi",
                         "LB_seagrass" = "Lameshur Bay seagrass")

group_labels_colors <- viridisLite::turbo(length(group_labels_lookup)) %>%
  setNames(., names(group_labels_lookup))

variable_labels_lookup <- c("Tektite - LB_seagrass" = "Tektite and Lameshur Bay seagrass",
                         "Yawzi - LB_seagrass" = "Yawzi and Lameshur Bay seagrass",
                         "Yawzi - Tektite" = "Yawzi and Tektite Reefs",
                         "peak_photo - dawn" = "afternoon and dawn")

variable_labels_colors <- viridisLite::turbo(length(variable_labels_lookup)) %>%
  setNames(., names(variable_labels_lookup))

contrast_labels_lookup <- data.frame(v1 = c("LB_seagrass (peak_photo - dawn)", 
                                            "Yawzi (peak_photo - dawn)", 
                                            "Tektite (peak_photo - dawn)", 
                                            "dawn (Yawzi - LB_seagrass)",
                                            "dawn (Tektite - LB_seagrass)",
                                            "dawn (Yawzi - Tektite)",
                                            "peak_photo (Yawzi - LB_seagrass)",
                                            "peak_photo (Tektite - LB_seagrass)",
                                            "peak_photo (Yawzi - Tektite)"),
                                     label = c("Lameshur Bay (afternoon - dawn)", 
                                               "Yawzi (afternoon - dawn)", 
                                               "Tektite (afternoon - dawn)", 
                                               "dawn (Yawzi - Lameshur Bay)",
                                               "dawn (Tektite - Lameshur Bay)",
                                               "dawn (Yawzi - Tektite)",
                                               "afternoon (Yawzi - Lameshur Bay)",
                                               "afternoon (Tektite - Lameshur Bay)",
                                               "afternoon (Yawzi - Tektite)")) %>%
  tibble::deframe(.)

site_lookup <- data.frame(site = c("LB_seagrass", "Yawzi", "Tektite",  "control_extraction", "control_pcr", "control_seq"),
                          label = c("Lameshur Bay seagrass", "Yawzi Reef", "Tektite Reef", 
                                    "Control (DNA Extraction)", "Control (PCR)", "Control (Sequencing)")) %>%
  tibble::deframe(.)
site_colors <- pals::kelly(22)[6:(5+length(site_lookup))] %>%
  # site_colors <- viridisLite::cividis(n = length(site_lookup), direction = 1) %>%
  setNames(., names(site_lookup))
sampling_time_lookup <- data.frame(sampling_time = c("dawn", "peak_photo"),
                                   label = c("Dawn", "Afternoon")) %>%
  tibble::deframe(.)
sampling_time_colors <- pals::ocean.haline(n = length(sampling_time_lookup)) %>%
  setNames(., names(sampling_time_lookup))
sampling_day_lookup <- data.frame(sampling_day = c("Day1", "Day2", "Day3", "Day4", "Day5"),
                                  label = c("20210122", "20210123", "20210124", "20210125", "20210126")) %>%
  tibble::deframe(.)
sampling_day_colors <- pals::ocean.thermal(n = length(sampling_day_lookup)) %>%
  setNames(., names(sampling_day_lookup))

enriched_labels_lookup <- c(site_lookup, sampling_time_lookup, group_labels_lookup)

usvi_sus_metabolites_idx <- data.frame(metabolites = c("2'deoxyguanosine", "HMP", "adenosine", "inosine", "pyridoxine", "4-aminobenzoic acid"))