
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
  )
)

# default inputs for some vaccine program properties:
default_v <- list(
  weeks_vaccinating = 12,
  first_year_all = T,
  NH_vacc_date = '01-10',
  SH_vacc_date = '01-04',
  init_vaccinated = c(0, 0, 0, 0)
)

cov1 <- 0.5
cov2 <- 0.2
cov3 <- 0.8

coverage_vector <- function(target, cov){
  if(! target %in% 1:5){
    stop("Target scheme must be in 1:5")
  }
  case_when(
    target == 1 ~ c(cov, 0, 0, 0),
    target == 2 ~ c(cov, cov*(6/15), 0, 0),
    target == 3 ~ c(cov, cov*(13/15), 0, 0),
    target == 4 ~ c(0, 0, 0, cov),
    target == 5 ~ c(cov, cov*(13/15), 0, cov)
  )
}

fcn_vacc_prog <- function(target_input, cov_input, relative_vacc_infec_input){
  output_list <- list(
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'current',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'improved_minimal',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'improved_efficacy',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'improved_breadth',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v),
    c(list(pop_coverage = coverage_vector(target_input, cov_input),
              vacc_type = 'universal',
              relative_vacc_infec = relative_vacc_infec_input),
      default_v)
  )
  names(output_list) <- paste0('vt_', 1:5, '_ct_', target_input)
  output_list
}

# vaccine programs, list unfinished:
vaccine_programs_base <- c(
  fcn_vacc_prog(1, cov1, 1),
  fcn_vacc_prog(2, cov1, 1),
  fcn_vacc_prog(3, cov1, 1),
  fcn_vacc_prog(4, cov1, 1),
  fcn_vacc_prog(5, cov1, 1)
)

vaccine_programs_low_cov <- c(
  fcn_vacc_prog(1, cov2, 1),
  fcn_vacc_prog(2, cov2, 1),
  fcn_vacc_prog(3, cov2, 1),
  fcn_vacc_prog(4, cov2, 1),
  fcn_vacc_prog(5, cov2, 1)
)

vaccine_programs_high_cov <- c(
  fcn_vacc_prog(1, cov3, 1),
  fcn_vacc_prog(2, cov3, 1),
  fcn_vacc_prog(3, cov3, 1),
  fcn_vacc_prog(4, cov3, 1),
  fcn_vacc_prog(5, cov3, 1)
)

vaccine_programs_rel_inf <- c(
  fcn_vacc_prog(1, cov1, 0.5),
  fcn_vacc_prog(2, cov1, 0.5),
  fcn_vacc_prog(3, cov1, 0.5),
  fcn_vacc_prog(4, cov1, 0.5),
  fcn_vacc_prog(5, cov1, 0.5)
)

vaccine_programs_merged <- list(
  vaccine_programs_base = vaccine_programs_base,
  vaccine_programs_low_cov = vaccine_programs_low_cov,
  vaccine_programs_high_cov = vaccine_programs_high_cov,
  vaccine_programs_rel_inf = vaccine_programs_rel_inf
)















