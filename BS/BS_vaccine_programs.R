
## SIMULATING X-YEAR TIME SERIES FOR EACH CLUSTER/COUNTRY IN INF A+B ##

## RUNNING WITH VARIOUS VACCINATION SCENARIOS ##

# vaccine types:
vacc_type_list <- list(
  current = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), # c(Matched u65, 65+, mismatched u65, 65+)
    mismatch = T
  ),
  improved_minimal = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = T
  ),
  improved_efficacy = list(
    imm_duration = 2,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.7, 0.4),
    mismatch = T
  ),
  improved_breadth = list(
    imm_duration = 3,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.7, 0.46),
    mismatch = F
  ),
  universal = list(
    imm_duration = 5,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.9, 0.7),
    mismatch = F
  ),
  current2 = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), # c(Matched u65, 65+, mismatched u65, 65+)
    mismatch = T
  ),
  improved_minimal2 = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = T
  ),
  improved_efficacy2 = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.7, 0.4),
    mismatch = T
  ),
  improved_breadth2 = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.7, 0.46),
    mismatch = F
  ),
  universal2 = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.9, 0.7, 0.9, 0.7),
    mismatch = F
  ),
  current3 = list(
    imm_duration = 0.5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28), # c(Matched u65, 65+, mismatched u65, 65+)
    mismatch = T
  ),
  improved_minimal3 = list(
    imm_duration = 1,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = T
  ),
  improved_efficacy3 = list(
    imm_duration = 2,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = T
  ),
  improved_breadth3 = list(
    imm_duration = 3,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = F
  ),
  universal3 = list(
    imm_duration = 5,
    coverage_pattern = 1,
    VE = c(0.7, 0.46, 0.42, 0.28),
    mismatch = F
  )
)

# default inputs for some vaccine program properties:
default_v <- list(
  weeks_vaccinating = 12,
  first_year_all = T,
  NH_vacc_date = '01-10',
  SH_vacc_date = '01-04',
  init_vaccinated = c(0, 0, 0, 0),
  relative_vacc_infec = 1
)

#coverage_vec <- c(0.5, 0.5, 0.0, 0.5)
coverage_vec <- c(0.6, 0.6, 0, 0.6)

# vaccine programs, list unfinished:
vaccine_programs <- list(
  V1 = list(pop_coverage = coverage_vec,
            vacc_type = 'current',
            weeks_vaccinating = default_v$weeks_vaccinating,
            first_year_all = default_v$first_year_all,
            NH_vacc_date = default_v$NH_vacc_date,
            SH_vacc_date = default_v$SH_vacc_date,
            init_vaccinated = default_v$init_vaccinated,
            relative_vacc_infec = default_v$relative_vacc_infec),
  V2 = list(pop_coverage = coverage_vec,
            vacc_type = 'improved_minimal',
            weeks_vaccinating = default_v$weeks_vaccinating,
            first_year_all = default_v$first_year_all,
            NH_vacc_date = default_v$NH_vacc_date,
            SH_vacc_date = default_v$SH_vacc_date,
            init_vaccinated = default_v$init_vaccinated,
            relative_vacc_infec = default_v$relative_vacc_infec),
  V3 = list(pop_coverage = coverage_vec,
            vacc_type = 'improved_efficacy',
            weeks_vaccinating = default_v$weeks_vaccinating,
            first_year_all = default_v$first_year_all,
            NH_vacc_date = default_v$NH_vacc_date,
            SH_vacc_date = default_v$SH_vacc_date,
            init_vaccinated = default_v$init_vaccinated,
            relative_vacc_infec = default_v$relative_vacc_infec),
  V4 = list(pop_coverage = coverage_vec,
            vacc_type = 'improved_breadth',
            weeks_vaccinating = default_v$weeks_vaccinating,
            first_year_all = default_v$first_year_all,
            NH_vacc_date = default_v$NH_vacc_date,
            SH_vacc_date = default_v$SH_vacc_date,
            init_vaccinated = default_v$init_vaccinated,
            relative_vacc_infec = default_v$relative_vacc_infec),
  V5 = list(pop_coverage = coverage_vec,
            vacc_type = 'universal',
            weeks_vaccinating = default_v$weeks_vaccinating,
            first_year_all = default_v$first_year_all,
            NH_vacc_date = default_v$NH_vacc_date,
            SH_vacc_date = default_v$SH_vacc_date,
            init_vaccinated = default_v$init_vaccinated,
            relative_vacc_infec = default_v$relative_vacc_infec)
  # V12 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'current2',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V22 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'improved_minimal2',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V32 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'improved_efficacy2',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V42 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'improved_breadth2',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V52 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'universal2',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V13 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'current3',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V23 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'improved_minimal3',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V33 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'improved_efficacy3',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V43 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'improved_breadth3',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec),
  # V53 = list(pop_coverage = coverage_vec,
  #           vacc_type = 'universal3',
  #           weeks_vaccinating = default_v$weeks_vaccinating,
  #           first_year_all = default_v$first_year_all,
  #           NH_vacc_date = default_v$NH_vacc_date,
  #           SH_vacc_date = default_v$SH_vacc_date,
  #           init_vaccinated = default_v$init_vaccinated,
  #           relative_vacc_infec = default_v$relative_vacc_infec)
)

















