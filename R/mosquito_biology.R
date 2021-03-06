#' @title Calculate equilibrium solution for vector counts
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param foim equilibrium foim
#' @param m (optional) the total number of female adult mosquitos
#' @noRd
initial_mosquito_counts <- function(parameters, foim = 0, m = NULL) {
  if (is.null(m)) {
    m <- parameters$total_M
  }
  omega <- calculate_omega(parameters)
  n_F <- 3 
  n_L <- 2 * parameters$mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_P <- 2 * parameters$dpl * parameters$mum * m

  n_T <-10

  n_Sm <- m * parameters$mum / (foim + parameters$mum)

  incubation_survival <- exp(-parameters$mum * parameters$dem)

  n_Pm <- m * foim / (foim + parameters$mum) * (
    1. - incubation_survival
  )

  n_Im <- m * foim / (foim + parameters$mum) * incubation_survival

  c(n_F, n_L, n_P, n_T, n_Sm, n_Pm, n_Im)
}

initial_mosquito_counts_ld <- function(parameters) {

  print(paste("KF ",parameters$KF," Amax ",parameters$Amax," me ", parameters$me," beta ", parameters$beta," G0 ",parameters$G0," mum ",parameters$mum))

  KF<-parameters$KF; Amax<-parameters$Amax; G0<-parameters$G0; beta<-parameters$beta; me<-parameters$me; mum<-parameters$mum 

  F_eq<-KF/(Amax*((1/me)*
  log(0.5*beta/mum))/(5*Amax)-1)
  print(paste("F_eq " , F_eq))
  L_eq<-G0*((1/me)*log(0.5*beta/mum))/(5*Amax)
  print(paste("L_eq " , L_eq))
  P_eq<-0.3*me*G0*(1/me)*log(0.5*beta/mum)/(5*Amax)/(beta-mum)
  print(paste("P_eq " , P_eq))
  T_eq<-(1/me)*log(0.5*beta/mum)
  print(paste("T_eq ", T_eq))
  c(F_eq,L_eq,P_eq,T_eq)

}


#' @title Calculate omega value
#' @description useful value for calculating equilibrium solutions for vectors
#' taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @noRd
calculate_omega <- function(parameters) {
  sub_omega <- parameters$gamma * parameters$ml / parameters$me - (
    parameters$del / parameters$dl
  ) + (
    (parameters$gamma - 1) * parameters$ml * parameters$del
  )

  -.5 * sub_omega + sqrt(
    .25 * sub_omega**2 +
      .5 * parameters$gamma * parameters$beta * parameters$ml * parameters$del /
      (parameters$me * parameters$mum * parameters$dl * (
        1. + parameters$dpl * parameters$mup
      ))
  )
}

#' @title Calculate the vector carrying capacity
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @noRd
calculate_carrying_capacity <- function(parameters) {
  m <- parameters$total_M
  omega <- calculate_omega(parameters)

  m * 2 * parameters$dl * parameters$mum * (
    1. + parameters$dpl * parameters$mup
  ) * parameters$gamma * (omega + 1) / (
    omega / (parameters$ml * parameters$del) - (
      1. / (parameters$ml * parameters$dl)
    ) - 1.
  )
}

#' @title Calculate the mean rainfall throughout the year
#' @param parameters model parameters
#' @noRd
calculate_R_bar <- function(parameters) {
  mean(vnapply(1:365, function(t) rainfall(
		t,
		parameters$days_per_timestep,
    parameters$g0,
    parameters$g,
    parameters$h
	)))
}

#' @title Calculate equilibrium total_M from parameters
#'
#' @param parameters to work from
#' @param EIR equilibrium to use, bites per person per year
#' @noRd
equilibrium_total_M <- function(parameters, EIR) {
  if (EIR == 0) {
    return(0)
  }
  if (parameters$init_foim == 0) {
    stop('init_foim must be > 0 to calculate a non-zero equilibrium total_M')
  }
  total_daily_eir <- EIR * parameters$human_population / 365
  lifetime <- parameters$init_foim * exp(-parameters$mum * parameters$dem) / (
    parameters$init_foim + parameters$mum
  )
  total_daily_eir / sum(
    parameters$species_proportions * parameters$blood_meal_rates * parameters$Q0 * lifetime
  )
}

#' @title Calculate the yearly offset (in timesteps) for the peak mosquito
#' season
#'
#' @param parameters to work from
#' @export
peak_season_offset <- function(parameters) {
  K0 <- calculate_carrying_capacity(parameters)
  R_bar <- calculate_R_bar(parameters)
  which.max(vnapply(seq(365), function(t) {
    carrying_capacity(
      t,
      parameters$model_seasonality,
      parameters$days_per_timestep,
      parameters$g0,
      parameters$g,
      parameters$h,
      K0,
      R_bar
    )
  }))[[1]]
}

#' @title Calculate the effects of biting on mosquito individuals
#'
#' @param variables a list of variables in this simulation
#' @param human_infectivity the infectivity for each human
#' @param lambda the effective biting rate for this species on each human
#' @param mosquito_infection an event for mosquito infection
#' @param species the index of the species to calculate for
#' @param susceptible_species the indices of susceptible mosquitos of the
#' species
#' @param adult_species the indices of adult mosquitos of the
#' species
#' @param W the mean probability that a mosquito feeds and survives
#' @param Z the mean probability that a mosquito is repelled
#' @param f the feeding rate for this species of mosquito
#' @param renderer the model renderer object
#' @param timestep the current timestep
#' @param parameters the model parameters
#' @noRd
calculate_mosquito_effects <- function(
    variables,
    human_infectivity,
    lambda,
    events,
    species,
    susceptible_species,
    adult_species,
    W,
    Z,
    f,
    parameters,
    renderer,
    timestep
  ) {
  # deal with mosquito infections
  lambda <- sum(human_infectivity * lambda)
  renderer$render(paste0('FOIM_', species), lambda, timestep)
  target <- sample_bitset(susceptible_species, lambda)
  variables$mosquito_state$queue_update('Pm', target)
  events$mosquito_infection$schedule(
    target,
    log_uniform(target$size(), parameters$dem)
  )

  # deal with mosquito deaths
  p1_0 <- exp(-parameters$mum * parameters$foraging_time)
  gonotrophic_cycle <- get_gonotrophic_cycle(species, parameters)
  p2 <- exp(-parameters$mum * gonotrophic_cycle)
  p1 <- p1_0 * W / (1 - Z * p1_0)
  mu <- -f * log(p1 * p2)
  died <- sample_bitset(adult_species, mu)
  renderer$render(paste0('mu_', species), mu, timestep)

  events$mosquito_death$schedule(died, 0)
}

get_gonotrophic_cycle <- function(v, parameters) {
  f <- parameters$blood_meal_rates[[v]]
  gonotrophic_cycle <- 1 / f - parameters$foraging_time
}
