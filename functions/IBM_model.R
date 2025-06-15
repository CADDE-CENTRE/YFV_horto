### Individual model
run_simulation2 <- function(seed, steps, dt, N, initial_infections, death_obs_prop, 
                            beta, initial_run, overall_run_length, transmission_type,
                            importation_rate,
                            EIP_gamma_shape,
                            EIP_gamma_rate,
                            latent_period_gamma_shape, 
                            latent_period_gamma_rate,
                            infectious_period_gamma_shape, 
                            infectious_period_gamma_rate,
                            death_observation_mixture_prob,
                            death_observation_mixture_exponential_rate,
                            death_observation_mixture_gamma_shape, 
                            death_observation_mixture_gamma_rate,
                            state) {
  
  ## Setting the seed
  set.seed(seed)
  
  ## Defining transmission type
  if (!(transmission_type %in% c("frequency_dependent", "density_dependent"))) {
    stop("specify transmission_type correctly")
  } 
  
  ## Define disease states
  health <- CategoricalVariable$new(categories = c("S", "E", "I", "D", "D_unobs", "Dobs"),
                                    initial_values = c(rep("S", N - initial_infections), 
                                                       rep("E", initial_infections)))
  
  ## Exposure process moving individuals S->E
  if (transmission_type == "density_dependent") {
    exposure_process <- function(t){
      
      ## Getting the index of all those still susceptible
      S_all <- health$get_index_of("S")  # get index of all monkeys still in susceptible state
      
      ## Infections from other PEL monkeys
      I <- health$get_index_of("I")
      I_inf <- I$size()                     # calculating the number of infectious monkeys contributing to the FOI
      foi <- beta * I_inf                   # calculate FOI experienced by susceptible monkeys
      p_inf <- 1 - exp(-foi * dt)           # converting the FOI (instantaneous rate) to a probability of being infected in the timestep
      local_infections <- S_all$copy()      # copy the susceptible set
      if (p_inf < 0) {
        print("p_inf is less than 0")
      }
      if (p_inf > 1) {
        print("p_inf is greater than 1")
      }
      local_infections$sample(rate = p_inf) # pick who gets infected locally
      
      ## Infections representing importations
      num_imports <- rpois(1, importation_rate * dt)
      if (num_imports > 0) {
        import_infections <- S_all$copy()
        import_infections$set_difference(local_infections)
        if (import_infections$size() > 0) {
          num_to_import <- min(num_imports, import_infections$size())
          if (is.na(num_to_import)) {
            num_to_import <- 0
          } else {
            health_render$render('num_to_import', 0, t)
          }
          health_render$render('num_to_import', num_to_import, t)
          import_infections$choose(num_to_import)
        } else {
          health_render$render('num_to_import', 0, t)
        }
        newly_infected <- local_infections
        newly_infected$or(import_infections)
        health$queue_update(value = "E",index = newly_infected)      # updating the health categorical variable
      } else {
        health_render$render('num_to_import', 0, t)
        health$queue_update(value = "E",index = local_infections)    # updating the health categorical variable
      }
    }
  } else if (transmission_type == "frequency_dependent") {
    exposure_process <- function(t){
      
      ## Getting the index of all those still susceptible
      S_all <- health$get_index_of("S")  # get index of all monkeys still in susceptible state
      
      ## Infections from other PEL monkeys
      I <- health$get_index_of("I")
      I_inf <- I$size()                            # calculating the number of infectious monkeys contributing to the FOI
      N <- health$get_size_of(c("S", "E", "I"))    # maybe this should be "N" instead (at which point we're assuming something semi-density dependent I think)
      if (N < 1) {
        foi <- 0                                   # calculate FOI experienced by susceptible monkeys
      } else {
        foi <- beta * I_inf / N                    # calculate FOI experienced by susceptible monkeys
      }
      p_inf <- 1 - exp(-foi * dt)                  # converting the FOI (instantaneous rate) to a probability of being infected in the timestep
      local_infections <- S_all$copy()             # copy the susceptible set
      local_infections$sample(rate = p_inf)        # pick who gets infected locally
      
      ## Infections representing importations
      num_imports <- rpois(1, importation_rate * dt)
      if (num_imports > 0) {
        import_infections <- S_all$copy()
        import_infections$set_difference(local_infections)
        if (import_infections$size() > 0) {
          num_to_import <- min(num_imports, import_infections$size())
          if (is.na(num_to_import)) {
            num_to_import <- 0
          }
          health_render$render('num_to_import', num_to_import, t)
          import_infections$choose(num_to_import)
        }
        newly_infected <- local_infections
        newly_infected$or(import_infections)
        health$queue_update(value = "E",index = newly_infected)      # updating the health categorical variable
      } else {
        num_to_import <- 0
        health_render$render('num_to_import', num_to_import, t)
        health$queue_update(value = "E",index = local_infections)    # updating the health categorical variable
      }
    }
  } else {
    stop("something wrong in transmission type specification")
  }
  
  ## Define infection event and process to schedule moves from E->I
  exposed_infectious_event <- TargetedEvent$new(population_size = N)
  exposed_infectious_event$add_listener(function(t, target) {
    health$queue_update("I", target)
  })
  exposed_progression_process <- function(t){
    E <- health$get_index_of("E")
    exposed_infectious_already_scheduled <- exposed_infectious_event$get_scheduled()
    E$and(exposed_infectious_already_scheduled$not(inplace = TRUE))
    infection_times <- rgamma(n = E$size(), 
                              shape = latent_period_gamma_shape, 
                              rate = latent_period_gamma_rate) / dt
    EIP_times <- rgamma(n = E$size(), 
                        shape = EIP_gamma_shape, 
                        rate = EIP_gamma_rate) / dt
    exposed_infection_progression_times <- infection_times + EIP_times
    exposed_infectious_event$schedule(target = E, delay = round(exposed_infection_progression_times))
    
  }
  
  ## Define death of I event and process to schedule moves from I->D
  infectious_death_event <- TargetedEvent$new(population_size = N)
  infectious_death_event$add_listener(function(t, target) {
    health$queue_update("D", target)
  })
  infectious_death_process <- function(t){
    I <- health$get_index_of("I")
    already_scheduled <- infectious_death_event$get_scheduled()
    I$and(already_scheduled$not(inplace = TRUE))
    death_times <- rgamma(n = I$size(), 
                          shape = infectious_period_gamma_shape, 
                          rate = infectious_period_gamma_rate) / dt
    infectious_death_event$schedule(target = I, delay = death_times)
  }
  
  ## Define event for individuals in D to go to either Dobs or D_unobs
  observation_event <- TargetedEvent$new(population_size = N)
  observation_event$add_listener(function(t, target) {
    health$queue_update("Dobs", target)
  })
  unobserved_event <- TargetedEvent$new(population_size = N)
  unobserved_event$add_listener(function(t, target) {
    health$queue_update("D_unobs", target)
  })
  observation_process <- function(t){
    
    ## Getting index of individuals in D state
    D <- health$get_index_of("D")
    
    ## Getting individuals already scheduled
    Dobs_already_scheduled <- observation_event$get_scheduled()
    Dunobs_already_scheduled <- unobserved_event$get_scheduled()
    
    ## Subsetting D inplace by those who aren't already scheduled
    D$and(Dobs_already_scheduled$not(inplace = TRUE))
    D$and(Dunobs_already_scheduled$not(inplace = TRUE))
    
    unscheduled_to_Dobs <- D$copy()
    unscheduled_to_Dobs <- unscheduled_to_Dobs$sample(death_obs_prop)
    unscheduled_to_D_unobs <- D$and(unscheduled_to_Dobs$not(inplace = FALSE))
    
    if (unscheduled_to_Dobs$size() > 0) {
      death_observation_mixture_dist_indicator <- rbinom(n = unscheduled_to_Dobs$size(), 
                                                         size = 1, 
                                                         prob = death_observation_mixture_prob)
      observation_times <- vector(mode = "numeric", length = length(death_observation_mixture_dist_indicator))
      for (i in 1:length(death_observation_mixture_dist_indicator)) {
        if (death_observation_mixture_dist_indicator[i] == 1) {
          observation_times[i] <- rexp(n = 1, 
                                       rate = death_observation_mixture_exponential_rate) / dt
        } else if (death_observation_mixture_dist_indicator[i] == 0) {
          observation_times[i] <- rgamma(n = 1, 
                                         shape = death_observation_mixture_gamma_shape, 
                                         rate = death_observation_mixture_gamma_rate) / dt
        } else {
          stop("something has gone wrong with death_observation_mixture_dist_indicator")
        }
      }
      observation_event$schedule(target = unscheduled_to_Dobs, delay = observation_times)
    }
    unobserved_event$schedule(target = unscheduled_to_D_unobs, delay = 1)
    
    ## Define render
    # health_render$render('D_obs_new', unscheduled_to_Dobs$size(), t)
    
  }
  
  ## Define render
  health_render <- Render$new(timesteps = steps)
  health_render_process <- categorical_count_renderer_process(
    renderer = health_render,
    variable = health,
    categories =  c("Dobs", "S")
  )
  
  ## Run simulation loop
  final_state <- simulation_loop(
    variables = list(health),
    events = list(exposed_infectious_event, observation_event, unobserved_event, infectious_death_event),
    processes = list(exposure_process, exposed_progression_process, 
                     infectious_death_process, observation_process, health_render_process),
    timesteps = steps,
    state = state)
  
  return(list(result = health_render$to_dataframe(), state = final_state))
}

