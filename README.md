# increment_cov_derivative
## Output (do_optimisation)
The function returns:
- $locs - locations *in the end* of algorithm
- $history_locs - locations at every step: side optimisation, rest locs optimisation, repeated rest locs optimisation.
- $history_det - Sigma_2 determinants after every step
- $iters - amount of iterations after first try (1. line) and after repetitions (2. line)
- $codes - codes indicating why algorithm stopped
  - 10 - reached left bound.
  - 11 - reached right bound.
  - 5 - marking time (←→ / →←)
  - 2 - same as 5, but different tracking method
  - 3.. - location has overlaped another point (not working, because R cheats: it gives you 2 points like |loc_1 - loc_2| = 1.386457e-15). Still try to fix it.

## Variables / parameters (do_optimisation)
### Obligatory to define
- xi_1: initial locations
- xi_2: increment locations
- sigma2: variance (squared standard deviation)
- phi: nuicance parameter for correlation function. Higher phi - lower influence of sigma2.


### Not obligatory to define
- l_rate: learning rate. Default value = 0.05
- bounds: restriction of coordiantes of the locations. Default value = c(0,20)
- iteration_stop: Iterations per every location. Default value = ((bounds[2]-bounds[1])/l_rate)*1.5
- repetitions: Amount of algorithm repeats. Default value = 2.
- debug_mode_global: debug-modes for every step. Obligated to have 3 elements. Default value = c(F,F,F).
- show_fails: Extra debug mode (specific message) for the same steps. Obligated to have 3 elements. Default value = c(F,F,F).
- draw_plots: draws plots at every step (shows xi-1, old xi_2 and new xi_2 locations). Default values = T (True).
