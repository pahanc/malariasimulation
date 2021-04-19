ODE_INDICES <- c(Fl = 1, L = 2, P = 3, Tl = 4)


parameterise_ode <- function(parameters) {
  lapply(
    parameters$species_proportions,
    function(p) {
      m <- p * parameters$total_M
      create_mosquito_model(
        initial_mosquito_counts_ld(parameters)[ODE_INDICES],
	parameters$beta,
        parameters$del,
        parameters$me,
        p * calculate_carrying_capacity(parameters),
        parameters$gamma,
        parameters$dl,
        parameters$ml,
        parameters$dpl,
        parameters$mup,
        m,
        parameters$model_seasonality,
        parameters$days_per_timestep,
        parameters$g0,
        parameters$g,
        parameters$h,
	parameters$history_f,
        parameters$history_m,
	parameters$G0,
        parameters$KF,
        parameters$Amax,
	parameters$mum,
	calculate_R_bar(parameters)
	)
    }
  )
}

create_ode_rendering_process <- function(renderer, odes) {
  function(timestep) {
    counts <- rep(0, length(ODE_INDICES))
    for (ode in odes) {
      row <- mosquito_model_get_states(ode)
      counts <- counts + row
    }
    for (i in seq_along(ODE_INDICES)) {
      renderer$render(
        paste0('mosquito_', names(ODE_INDICES)[[i]], '_count'),
        counts[[i]],
        timestep
      )
    }
  }
}

#' @title Step mosquito ODE
#' @description calculates total_M per species and updates the vector ode
#'
#' @param odes the models to step, one for each species
#' @param state the mosquito state variable
#' @param species the mosquito species variable
#' @param species_names the names of the mosquito species
#' @param renderer the model renderer
#' @noRd
create_ode_stepping_process <- function(odes, state, species, species_names, renderer) {
  function(timestep) {
    adult <- state$get_index_of("NonExistent")$not()
    for (s_i in seq_along(species_names)) {
      total_M <- species$get_index_of(species_names[[s_i]])$and(adult)$size()
      renderer$render(paste0('total_M_', s_i), total_M, timestep)
      mosquito_model_step(odes[[s_i]], total_M)
    }
    renderer$render('total_M', adult$size(), timestep)
  }
}
