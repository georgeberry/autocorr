source("sim_code.R")
#### Run and save #############################################################

# main
print('Main node')
df_main_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_main_node
)
print('Main edge')
write_csv(
  df_main_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_runs.csv'
)
df_main_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_main_edge
)
write_csv(
  df_main_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_edge_runs.csv'
)
# main_perf
print('Main perf node')
df_node_perf = fn_par_bootstrap(
  fn_performance,
  N_SIMS,
  fn_dgp_main_node
)
write_csv(
  df_node_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_perf.csv'
)
print('Main perf edge')
df_edge_perf = fn_par_bootstrap(
  fn_performance,
  N_SIMS,
  fn_dgp_main_edge
)
write_csv(
  df_edge_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_edge_perf.csv'
)

# indep
print('Indep node')
df_indep_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_indep_node
)
write_csv(
  df_indep_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/indep_node_runs.csv'
)
print('Indep edge')
df_indep_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_indep_edge
)
write_csv(
  df_indep_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/indep_edge_runs.csv'
)

# degcorr
print('Degcorr node')
df_degcorr_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_degcorr_node
)
write_csv(
  df_degcorr_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/degcorr_node_runs.csv'
)
print('Degcorr edge')
df_degcorr_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_degcorr_edge
)
write_csv(
  df_degcorr_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/degcorr_edge_runs.csv'
)

# unobs
print('Unobs node')
df_unobs_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_unobs_node
)
write_csv(
  df_unobs_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/unobs_node_runs.csv'
)
print('Unobs edge')
df_unobs_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_unobs_edge
)
write_csv(
  df_unobs_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/unobs_edge_runs.csv'
)

# sampling
print('Sampling node')
df_sampling_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_sampling_node
)
write_csv(
  df_sampling_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_node_runs.csv'
)
print('Sampling edge')
df_sampling_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_sampling_edge
)
write_csv(
  df_sampling_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_edge_runs.csv'
)