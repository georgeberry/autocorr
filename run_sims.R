library(tidyverse)
library(ggraph)
library(tidygraph)
library(foreach)
library(doParallel)
library(pROC)
library(caret)

registerDoParallel(cores=12)
N_SIMS = 500
N_NODES = 4000
NODE_SAMP_FRAC = 0.2
EDGE_SAMP_FRAC = 0.05
POWERLAW_EXPONENT = 0.8
EDGES_PER_NEW_NODE = 5

#### Functions ################################################################

fn_par_bootstrap = function(fn, n_reps, ...) {
  df_results = foreach(
    i=1:n_reps,
    .combine=rbind
  ) %dopar% {
    result = fn(...)
    result$rep = i
    return(result)
  }
  return(df_results)
}

######## Main sims ############################################################

fn_dgp_main_node = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      # ZZ = rnorm(n()), 
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          max(as_tibble(neighborhood, active='nodes')$X)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)), #  - log1p(outdeg))),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      ),
      gt = rbinom(n(), 1, NODE_SAMP_FRAC)
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
      gt = as.numeric(gt_ego & gt_nbr),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    )
  return(g)
}


fn_dgp_main_edge = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      # ZZ = rnorm(n()), 
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          max(as_tibble(neighborhood, active='nodes')$X)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)), #  - log1p(outdeg))),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      )
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      # sample at edge level
      gt = rbinom(n(), 1, EDGE_SAMP_FRAC),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    ) %>%
    # for any node who is part of a gt edge, give that node gt=1
    activate(nodes) %>%
    mutate(
      # if any edge containing the focal node is 1, the node is gt=1
      gt=map_local_dbl(
        .f = function(neighborhood, node, ...) {
          df_tmp = as_tibble(neighborhood, active='edges') %>%
            filter(gt==1)
          min(nrow(df_tmp), 1) # if 0 in gt, return 0 else 1
        }
      )
    ) %>%
    # now go back and label the gt_ego/gt_nbr based on the node being included
    # in any ground-truth edge
    activate(edges) %>%
    mutate(
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
    )
    
  return(g)
}


######## Indep sims ###########################################################

fn_dgp_indep_node = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      Z = rnorm(n()),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      ),
      gt = rbinom(n(), 1, NODE_SAMP_FRAC)
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
      gt = as.numeric(gt_ego & gt_nbr),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    )
  return(g)
}


fn_dgp_indep_edge = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      Z = rnorm(n()), 
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      )
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      # sample at edge level
      gt = rbinom(n(), 1, EDGE_SAMP_FRAC),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    ) %>%
    # for any node who is part of a gt edge, give that node gt=1
    activate(nodes) %>%
    mutate(
      # if any edge containing the focal node is 1, the node is gt=1
      gt=map_local_dbl(
        .f = function(neighborhood, node, ...) {
          df_tmp = as_tibble(neighborhood, active='edges') %>%
            filter(gt==1)
          min(nrow(df_tmp), 1) # if 0 in gt, return 0 else 1
        }
      )
    ) %>%
    # now go back and label the gt_ego/gt_nbr based on the node being included
    # in any ground-truth edge
    activate(edges) %>%
    mutate(
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
    )
    
  return(g)
}


######## Degree correlation sims ##############################################


fn_dgp_degcorr_node = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - 0.03 * outdeg)),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      ),
      gt = rbinom(n(), 1, NODE_SAMP_FRAC)
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
      gt = as.numeric(gt_ego & gt_nbr),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    )
  return(g)
}


fn_dgp_degcorr_edge = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - 0.03 * outdeg)),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      )
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      # sample at edge level
      gt = rbinom(n(), 1, EDGE_SAMP_FRAC),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    ) %>%
    # for any node who is part of a gt edge, give that node gt=1
    activate(nodes) %>%
    mutate(
      # if any edge containing the focal node is 1, the node is gt=1
      gt=map_local_dbl(
        .f = function(neighborhood, node, ...) {
          df_tmp = as_tibble(neighborhood, active='edges') %>%
            filter(gt==1)
          min(nrow(df_tmp), 1) # if 0 in gt, return 0 else 1
        }
      )
    ) %>%
    # now go back and label the gt_ego/gt_nbr based on the node being included
    # in any ground-truth edge
    activate(edges) %>%
    mutate(
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
    )
    
  return(g)
}

######## Unobs sims ###########################################################

fn_dgp_corr_unobs_node = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      ZZ = rnorm(n()), 
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          max(as_tibble(neighborhood, active='nodes')$ZZ)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      ),
      gt = rbinom(n(), 1, NODE_SAMP_FRAC)
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
      gt = as.numeric(gt_ego & gt_nbr),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    )
  return(g)
}


fn_dgp_corr_unobs_edge = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      ZZ = rnorm(n()), 
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          max(as_tibble(neighborhood, active='nodes')$ZZ)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      )
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      # sample at edge level
      gt = rbinom(n(), 1, EDGE_SAMP_FRAC),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    ) %>%
    # for any node who is part of a gt edge, give that node gt=1
    activate(nodes) %>%
    mutate(
      # if any edge containing the focal node is 1, the node is gt=1
      gt=map_local_dbl(
        .f = function(neighborhood, node, ...) {
          df_tmp = as_tibble(neighborhood, active='edges') %>%
            filter(gt==1)
          min(nrow(df_tmp), 1) # if 0 in gt, return 0 else 1
        }
      )
    ) %>%
    # now go back and label the gt_ego/gt_nbr based on the node being included
    # in any ground-truth edge
    activate(edges) %>%
    mutate(
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
    )
    
  return(g)
}



######## Sampling sims ########################################################

fn_dgp_sampling_node = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      # ZZ = rnorm(n()), 
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          max(as_tibble(neighborhood, active='nodes')$X)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)), #  - log1p(outdeg))),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      ),
      # probably we see the true value
      Q_prob = 1 / (1 + exp(1.6 - 0.03 * outdeg - 0.2 * X)),
      gt = rbinom(graph_order(), 1, Q_prob)
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
      gt = as.numeric(gt_ego & gt_nbr),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    )
  return(g)
}


fn_dgp_sampling_edge = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_barabasi_albert(N_NODES, POWERLAW_EXPONENT, EDGES_PER_NEW_NODE)
  g = g %>%
    bind_edges(
      g %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(tmp=to, to=from, from=tmp) %>%
        select(to, from)
    ) %>%
    activate(nodes) %>%
    mutate(
      X = rnorm(n()),
      # ZZ = rnorm(n()), 
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          max(as_tibble(neighborhood, active='nodes')$X)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X - Z)), #  - log1p(outdeg))),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      Y = rbinom(
        n(),
        1,
        Y_prob_true
      )
    ) %>%
    activate(edges) %>%
    mutate(
      Y_ego = .N()$Y[from],
      Y_nbr = .N()$Y[to],
      Y = Y_ego * Y_nbr,
      indeg_ego = .N()$indeg[from],
      indeg_nbr = .N()$indeg[to],
      outdeg_ego = .N()$outdeg[from],
      outdeg_nbr = .N()$outdeg[to],
      indeg_inv_ego = .N()$indeg_inv[from],
      indeg_inv_nbr = .N()$indeg_inv[to],
      outdeg_inv_ego = .N()$outdeg_inv[from],
      outdeg_inv_nbr = .N()$outdeg_inv[to],
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
      # sample at edge level
      # probably we see the true value
      # 1 / (1 + exp(1.6 - 0.03 * outdeg - 0.2 * X))
      Q_prob = 1 / (1 + exp(
        7.4 -
        0.03 * outdeg_ego -
        0.03 * indeg_nbr - 
        0.1 * X_ego - 
        0.1 * X_nbr
      )),
      gt = rbinom(graph_size(), 1, Q_prob)
    ) %>%
    # for any node who is part of a gt edge, give that node gt=1
    activate(nodes) %>%
    mutate(
      # if any edge containing the focal node is 1, the node is gt=1
      gt=map_local_dbl(
        .f = function(neighborhood, node, ...) {
          df_tmp = as_tibble(neighborhood, active='edges') %>%
            filter(gt==1)
          min(nrow(df_tmp), 1) # if 0 in gt, return 0 else 1
        }
      )
    ) %>%
    # now go back and label the gt_ego/gt_nbr based on the node being included
    # in any ground-truth edge
    activate(edges) %>%
    mutate(
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
    )
    
  return(g)
}

######## Compute performance metrics ##########################################

fn_compute_performance_metrics = function(true_vals, probs) {
  auc_val = as.numeric(auc(roc(true_vals, probs)))
  # zero-one predictions
  preds = probs > 0.5
  prec_val = sum(true_vals * preds) / sum(preds)
  rec_val = sum(true_vals * preds) / sum(true_vals)
  acc_val = sum(true_vals == preds) / length(true_vals)
  return(
    data.frame(
      auc=auc_val,
      precision=prec_val,
      recall=rec_val,
      accuracy=acc_val
    )
  )
}

fn_performance = function(fn_g) {
  # fit models, evlauate performance metrics
  # node no network features (baseline)
  # node
  # edge
  # ego-then-alter
  # ego-then-alter (with prob passed to alter model)
  
  g = fn_g()
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()

  # node no network features
  node_mod = glm(
    Y ~ X, # + log1p(outdeg),
    family='binomial',
    data=df_nodes %>%
      filter(gt == 1)
  )
  node_preds = predict(
    node_mod,
    newdata=df_nodes %>%
      filter(gt == 0),
    type='response',
  )
  node_nonetwork_perf = fn_compute_performance_metrics(
    df_nodes %>%
      filter(gt == 0) %>%
      .$Y,
    node_preds
  ) %>%
    mutate(
      model='node_nonetwork'
    )
  
  # nodemod
  node_mod = glm(
    Y ~ X + outdeg_inv + outdeg,
    family='binomial',
    data=df_nodes %>%
      filter(gt == 1)
  )
  node_preds = predict(
    node_mod,
    newdata=df_nodes %>%
      filter(gt == 0),
    type='response',
  )
  node_perf = fn_compute_performance_metrics(
    df_nodes %>%
      filter(gt == 0) %>%
      .$Y,
    node_preds
  ) %>%
    mutate(
      model='node'
    )
  
  # ego then alter
  ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1)
  )
  df_edges$ego_preds = predict(
    ego_mod,
    newdata=df_edges,
    type='response',
  )
  df_edges = df_edges %>%
    group_by(from) %>%
    summarize(
      gt_ego=mean(gt_ego),
      ego_preds=mean(ego_preds),
      Y_ego=mean(Y_ego)
    )
  print(head(df_edges))
  egoalter_perf = fn_compute_performance_metrics(
    df_edges %>%
      filter(gt_ego == 0) %>%
      .$Y_ego,
    df_edges %>%
      filter(gt_ego == 0) %>%
      .$ego_preds
  ) %>%
    mutate(
      model='egoalter'
    )

  return(
    bind_rows(
      node_perf,
      node_nonetwork_perf,
      egoalter_perf
    )
  )
}


#### Run and save #############################################################

# main
print('Main node')
df_main_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_main_node,
)
print('Main edge')
write_csv(
  df_main_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_runs.csv'
)
df_main_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_main_edge,
)
write_csv(
  df_main_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_edge_runs.csv'
)
# main perf
print('Main perf node')
df_node_perf = fn_par_bootstrap(fn_performance, N_SIMS, fn_dgp_main_node)
write_csv(
  df_node_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_node_perf.csv'
)
print('Main perf edge')
df_edge_perf = fn_par_bootstrap(fn_performance, N_SIMS, fn_dgp_main_edge)
write_csv(
  df_edge_perf,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_edge_perf.csv'
)

# indep
print('Indep node')
df_indep_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_indep_node,
)
write_csv(
  df_indep_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/indep_node_runs.csv'
)
print('Indep edge')
df_indep_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_indep_edge,
)
write_csv(
  df_indep_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/main_edge_runs.csv'
)

# degcorr
print('Degcorr node')
df_degcorr_node_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_degcorr_node,
)
write_csv(
  df_degcorr_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/degcorr_node_runs.csv'
)
print('Degcorr edge')
df_degcorr_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_degcorr_edge,
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
  fn_dgp_unobs_node,
)
write_csv(
  df_unobs_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/unobs_node_runs.csv'
)
print('Unobs edge')
df_unobs_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_unobs_edge,
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
  fn_dgp_sampling_node,
)
write_csv(
  df_sampling_node_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_node_runs.csv'
)
print('Sampling edge')
df_sampling_edge_runs = fn_par_bootstrap(
  fn_oos,
  N_SIMS,
  fn_dgp_sampling_edge,
)
write_csv(
  df_sampling_edge_runs,
  '/Users/georgeberry/Dropbox/project-autocorr/data/sampling_edge_runs.csv'
)