source("sim_code.R")

registerDoParallel(cores=12)
N_SIMS = 100
N_NODES = 4000
GROUND_TRUTH_LABELING_BUDGET = 500
ALPHA = 0.8  # powerlaw exponent
BETA = 0.7  # strength of homophily
EDGES_PER_NEW_NODE = 5
MAJORITY_GROUP_FRAC = 0.7

#### NEW GRAPH #################################################################

df_new_runs_no_homophily = fn_par_run(
  fn_run_sims,
  N_SIMS,

  # sim params
  N_NODES,
  EDGES_PER_NEW_NODE,
  MAJORITY_GROUP_FRAC,
  ALPHA,
  0.0
)
write_csv(
  df_new_runs_no_homophily,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_sims_no_homophily.csv'
)

df_new_runs_low_homophily = fn_par_run(
  fn_run_sims,
  N_SIMS,

  # sim params
  N_NODES,
  EDGES_PER_NEW_NODE,
  MAJORITY_GROUP_FRAC,
  ALPHA,
  0.4
)
write_csv(
  df_new_runs_low_homophily,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_sims_low_homophily.csv'
)

df_new_runs_main = fn_par_run(
  fn_run_sims,
  N_SIMS,

  # sim params
  N_NODES,
  EDGES_PER_NEW_NODE,
  MAJORITY_GROUP_FRAC,
  ALPHA,
  BETA
)
write_csv(
  df_new_runs_main,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_sims_main.csv'
)


df_new_runs_high_homophily = fn_par_run(
  fn_run_sims,
  N_SIMS,

  # sim params
  N_NODES,
  EDGES_PER_NEW_NODE,
  MAJORITY_GROUP_FRAC,
  ALPHA,
  1.2
)
write_csv(
  df_new_runs_high_homophily,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_sims_high_homophily.csv'
)

#### Performance ###############################################################

df_perf = fn_par_run(
  fn_run_perf,
  N_SIMS,

  # sim params
  N_NODES,
  EDGES_PER_NEW_NODE,
  MAJORITY_GROUP_FRAC,
  ALPHA,
  BETA
)
write_csv(
  df_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_sims_perf.csv'
)




if (FALSE) {
# indep
print('New graph indep node')
df_new_graph_indep_node_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_new_dgp_indep_node
)
write_csv(
  df_new_graph_indep_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_indep_node_runs.csv'
)

print('New graph indep edge')
df_new_graph_indep_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_new_dgp_equal_node
)
write_csv(
  df_new_graph_equal_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/new_graph_equal_node_runs.csv'
)

print('New graph equal edge')
df_new_graph_equal_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_dgp_main_node
)
print('Main edge')
write_csv(
  df_main_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_runs.csv'
)
df_main_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_dgp_indep_node
)
write_csv(
  df_indep_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/indep_node_runs.csv'
)
print('Indep edge')
df_indep_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_dgp_degcorr_node
)
write_csv(
  df_degcorr_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/degcorr_node_runs.csv'
)
print('Degcorr edge')
df_degcorr_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_dgp_unobs_node
)
write_csv(
  df_unobs_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/unobs_node_runs.csv'
)
print('Unobs edge')
df_unobs_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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
  fn_gen_graph_and_fit_models,
  N_SIMS,
  fn_dgp_sampling_node
)
write_csv(
  df_sampling_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_node_runs.csv'
)
print('Sampling edge')
df_sampling_edge_runs = fn_par_run(
  fn_gen_graph_and_fit_models,
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