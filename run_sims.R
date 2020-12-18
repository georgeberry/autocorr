source("sim_code.R")

N_SIMS = 100

#### NEW GRAPH #################################################################

print('New graph main node nodesamp')
df_new_graph_main_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_main_node
)
write_csv(
  df_new_graph_main_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_main_node_runs.csv'
)

print('New graph main node degsamp')
df_new_graph_main_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_main_node
)
write_csv(
  df_new_graph_main_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_main_node_runs.csv'
)

print('New graph main edge edgesamp')
df_new_graph_main_edge_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_main_edge
)
write_csv(
  df_new_graph_main_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_main_edge_runs.csv'
)

if (FALSE) {
# indep
print('New graph indep node')
df_new_graph_indep_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_indep_node
)
write_csv(
  df_new_graph_indep_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_indep_node_runs.csv'
)

print('New graph indep edge')
df_new_graph_indep_edge_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_indep_edge
)
write_csv(
  df_new_graph_indep_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_indep_edge_runs.csv'
)

# equal
print('New graph equal node')
df_new_graph_equal_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_equal_node
)
write_csv(
  df_new_graph_equal_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_equal_node_runs.csv'
)

print('New graph equal edge')
df_new_graph_equal_edge_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_new_dgp_equal_edge
)
write_csv(
  df_new_graph_equal_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_equal_edge_runs.csv'
)

# main_perf
print('Main perf node')
df_new_node_perf = fn_par_run(
  fn_performance,
  N_SIMS,
  fn_new_dgp_main_node
)
write_csv(
  df_new_node_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_main_node_perf.csv'
)
print('Main perf edge')
df_new_edge_perf = fn_par_run(
  fn_performance,
  N_SIMS,
  fn_new_dgp_main_edge
)
write_csv(
  df_new_edge_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_main_edge_perf.csv'
)

#### OLD GRAPH #################################################################

# main
print('Main node')
df_main_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_dgp_main_node
)
print('Main edge')
write_csv(
  df_main_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_runs.csv'
)
df_main_edge_runs = fn_par_run(
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
df_node_perf = fn_par_run(
  fn_performance,
  N_SIMS,
  fn_dgp_main_node
)
write_csv(
  df_node_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_perf.csv'
)
print('Main perf edge')
df_edge_perf = fn_par_run(
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
df_indep_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_dgp_indep_node
)
write_csv(
  df_indep_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/indep_node_runs.csv'
)
print('Indep edge')
df_indep_edge_runs = fn_par_run(
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
df_degcorr_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_dgp_degcorr_node
)
write_csv(
  df_degcorr_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/degcorr_node_runs.csv'
)
print('Degcorr edge')
df_degcorr_edge_runs = fn_par_run(
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
df_unobs_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_dgp_unobs_node
)
write_csv(
  df_unobs_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/unobs_node_runs.csv'
)
print('Unobs edge')
df_unobs_edge_runs = fn_par_run(
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
df_sampling_node_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_dgp_sampling_node
)
write_csv(
  df_sampling_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_node_runs.csv'
)
print('Sampling edge')
df_sampling_edge_runs = fn_par_run(
  fn_oos,
  N_SIMS,
  fn_dgp_sampling_edge
)
write_csv(
  df_sampling_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_edge_runs.csv'
)



#### Bootstrap??
#### 
#### this code creates a graph, and does N_SIMS resamples from the graph

if (FALSE) {
BOOT_REPS = 100

df_boot = foreach(
  i=1:BOOT_REPS,
  .combine=rbind
) %dopar% {
  g = fn_dgp_indep_edge()
  res = fn_naive_bootstrap_edge(g)
  res$rep = i
  return(res)
}

vars_to_summarize = vars(
  Y_true,
  Y_hat_edge
  # Y_hat_node
)

df_boot %>%
  group_by(rep) %>%
  summarize_at(
    vars_to_summarize,
    list(
      ~ quantile(., probs=c(0.05)),
      ~ quantile(., probs=c(0.95))
    )
  ) %>%
  mutate(
    in_range = Y_hat_edge_quantile..1 <= Y_true_quantile..1 &
    Y_true_quantile..1 <= Y_hat_edge_quantile..2,
  ) %>%
  summarize(
    mean(in_range)
  )



df_boot = foreach(
  i=1:500,
  .combine=rbind
) %dopar% {
  g = fn_dgp_indep_node()
  res = fn_naive_bootstrap_node(g)
  res$rep = i
  return(res)
}

vars_to_summarize = vars(
  Y_true,
  Y_hat_node
)

df_boot %>%
  group_by(rep) %>%
  summarize_at(
    vars_to_summarize,
    list(
      ~ quantile(., probs=c(0.1)),
      ~ quantile(., probs=c(0.9))
    )
  ) %>%
  mutate(
    in_range = Y_hat_node_quantile..1 <= Y_true_quantile..1 &
      Y_true_quantile..1 <= Y_hat_node_quantile..2,
  ) %>%
  summarize(
    mean(in_range)
  )
}
}