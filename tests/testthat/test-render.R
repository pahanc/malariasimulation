test_that('that default rendering works', {
  timestep <- 0
  parameters <- get_parameters()
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S'),
    c('U', 'A', 'D', 'S')
  )
  birth <- individual::DoubleVariable$new(
    -c(2, 5, 10, 11) * 365
  )
  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    rep('no', 4)
  )

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
    is_severe,
    parameters,
    renderer
  )

  process(timestep)

  mockery::expect_args(
    renderer$render,
    1,
    'pv_730_3650',
    2/3,
    timestep
  )
})

test_that('that default rendering works when no one is in the age range', {
  timestep <- 0
  parameters <- get_parameters()
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S'),
    rep('S', 4)
  )
  birth <- individual::DoubleVariable$new(
    -c(1, 11, 21, 11) * 365
  )
  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    rep('no', 4)
  )

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
    is_severe,
    parameters,
    renderer
  )
  process(timestep)

  mockery::expect_args(
    renderer$render,
    1,
    'pv_730_3650',
    0,
    timestep
  )
})

test_that('that severe rendering works', {
  timestep <- 0
  year <- 365
  parameters <- get_parameters(list(
    severe_prevalence_rendering_min_ages = c(0, 2) * year,
    severe_prevalence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S'),
    c('U', 'D', 'D', 'S')
  )
  birth <- individual::DoubleVariable$new(
    -c(2, 5, 10, 11) * 365
  )
  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    c('no', 'yes', 'no', 'no')
  )
  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
    is_severe,
    parameters,
    renderer
  )
  process(timestep)

  mockery::expect_args(
    renderer$render,
    1,
    'pv_severe_0_1825',
    1/2,
    timestep
  )

  mockery::expect_args(
    renderer$render,
    2,
    'pv_severe_730_3650',
    1/3,
    timestep
  )
})

test_that('that incidence rendering works', {
  timestep <- 0
  year <- 365
  parameters <- get_parameters(list(
    incidence_rendering_min_ages = c(0, 2) * year,
    incidence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))

  birth <- individual::DoubleVariable$new(
    -c(2, 5, 10, 11) * 365
  )

  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    c('no', 'yes', 'no', 'no')
  )

  renderer <- mock_render(1)
  process <- create_incidence_renderer(birth, is_severe, parameters, renderer)

  process(timestep, individual::Bitset$new(4)$insert(c(1, 2, 4)))

  mockery::expect_args(
    renderer$render,
    1,
    'inc_0_1825',
    1,
    timestep
  )

  mockery::expect_args(
    renderer$render,
    2,
    'inc_730_3650',
    2/3,
    timestep
  )
})

test_that('that severe incidence rendering works', {
  timestep <- 0
  year <- 365
  parameters <- get_parameters(list(
    severe_incidence_rendering_min_ages = c(0, 2) * year,
    severe_incidence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))

  birth <- individual::DoubleVariable$new(
    -c(2, 6, 10, 11) * 365
  )

  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    c('no', 'yes', 'no', 'no')
  )

  renderer <- mock_render(1)
  process <- create_incidence_renderer(birth, is_severe, parameters, renderer)

  process(timestep, individual::Bitset$new(4)$insert(c(1, 2, 4)))

  mockery::expect_args(
    renderer$render,
    1,
    'inc_severe_0_1825',
    0,
    timestep
  )

  mockery::expect_args(
    renderer$render,
    2,
    'inc_severe_730_3650',
    1/3,
    timestep
  )
})
