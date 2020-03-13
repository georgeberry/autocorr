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
EDGE_SAMP_FRAC = 0.025
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

fn_oos = function(fn_g) {
  # fit 4 models
  # node
  # edge
  # ego-then-alter
  # ego-then-alter (with prob passed to alter model)
  
  g = fn_g()
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()

  # node no network features
  node_mod = glm(
    Y ~ X,
    family='binomial',
    data=df_nodes %>%
      filter(gt == 1)
  )
  node_preds = predict(node_mod, newdata=df_nodes, type='response')
  Y_hat_node_nonetwork = sum(
    df_edges$outdeg_inv_ego *
    node_preds[df_edges$from] *
    node_preds[df_edges$to]
  )
  node_mod = glm(
    Y ~ X + outdeg + outdeg_inv,
    family='binomial',
    data=df_nodes %>%
      filter(gt == 1)
  )
  node_preds = predict(node_mod, newdata=df_nodes, type='response')
  Y_hat_node = sum(
    df_edges$outdeg_inv_ego *
    node_preds[df_edges$from] *
    node_preds[df_edges$to]
  )
  
  # edge
  edge_mod = glm(
    Y ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt == 1)
  )
  Y_hat_edge = sum(
    df_edges$outdeg_inv_ego *
    predict(edge_mod, newdata=df_edges, type='response')
  )
  
  # ego then alter
  ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1)
  )
  df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
  nbr_mod = glm(
    Y_nbr ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_nbr == 1)
  )
  df_edges$Y_nbr_hat = predict(nbr_mod, newdata=df_edges, type='response')
  Y_hat_egoalter = sum(
    df_edges$outdeg_inv_ego *
    df_edges$Y_ego_hat *
    df_edges$Y_nbr_hat
  )
  
  # clear
  df_edges$Y_ego_hat = NULL
  df_edges$Y_nbr_hat = NULL
  
  # ego then alter with prob passed
  ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1)
  )
  df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
  nbr_mod = glm(
    Y_nbr ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego + Y_ego_hat,
    family='binomial',
    data=df_edges %>%
      filter(gt_nbr == 1)
  )
  df_edges$Y_nbr_hat = predict(nbr_mod, newdata=df_edges, type='response')
  Y_hat_egoalter_passed = sum(
    df_edges$outdeg_inv_ego *
    df_edges$Y_ego_hat *
    df_edges$Y_nbr_hat
  )

  # ego then alter hardmode
  ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1) %>%
      sample_frac(0.5) # this sampling makes it hard
  )
  df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
  nbr_mod = glm(
    Y_nbr ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego + Y_ego_hat,
    family='binomial',
    data=df_edges %>%
      filter(gt_nbr == 1) %>%
      sample_frac(0.5)
  )
  df_edges$Y_nbr_hat = predict(nbr_mod, newdata=df_edges, type='response')
  Y_hat_egoalter_hardmode = sum(
    df_edges$outdeg_inv_ego *
    df_edges$Y_ego_hat *
    df_edges$Y_nbr_hat
  )

  # with weights
  if (FALSE) {
    ego_wghts = glm(
      gt_ego ~ X_ego + X_nbr + outdeg_inv_ego + outdeg_inv_nbr,
      family='binomial',
      data=df_edges
    )
    nbr_wghts = glm(
      gt_nbr ~ X_ego + X_nbr + outdeg_inv_ego + outdeg_inv_nbr,
      family='binomial',
      data=df_edges
    )
    
    ego_mod = glm(
      Y_ego ~ X_ego + X_nbr +  outdeg_inv_ego,
      family='binomial',
      data=df_edges %>%
        filter(gt_ego == 1),
      weights=1/predict(ego_wghts, type='response', newdata=df_edges[which(df_edges$gt_ego == 1),])
        # 1/outdeg_ego
      #df_edges %>%
      #  filter(gt_ego == 1) %>%
      #  group_by(from) %>%
      #  summarize(
      #    Y_ego=mean(Y_ego),
      #    X_ego=mean(X_ego),
      #    X_nbr=mean(X_nbr),
      #    outdeg_inv_ego=mean(outdeg_inv_ego),
      #    outdeg_ego=mean(outdeg_ego),
      #    indeg_nbr=mean(indeg_nbr)
    )
    df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
    nbr_mod = glm(
      Y_nbr ~ X_ego + X_nbr + outdeg_inv_ego + Y_ego_hat,
      family='binomial',
      data=df_edges %>%
        filter(gt_nbr == 1),
      weights=1/predict(nbr_wghts, type='response', newdata=df_edges[which(df_edges$gt_nbr == 1),])
      #df_edges %>%
      #  filter(gt_nbr == 1) %>%
      #  group_by(to) %>%
      #  summarize(
      #    Y_nbr=mean(Y_nbr),
      #    X_ego=mean(X_ego),
      #    X_nbr=mean(X_nbr),
      #    outdeg_inv_ego=mean(outdeg_inv_ego),
      #    outdeg_ego=mean(outdeg_ego),
      #    indeg_nbr=mean(indeg_nbr),
      #    Y_ego_hat=mean(Y_ego_hat)
      #  )
    )
  }
  # True values in the network
  Y_true = sum(df_edges$outdeg_inv_ego * df_edges$Y)

  Y_no_model = sum(
    df_edges %>% filter(gt == 1) %>% .$outdeg_inv_ego * 
    df_edges %>% filter(gt == 1) %>% .$Y
  ) * (nrow(df_edges) / sum(df_edges$gt))

  return(
    data.frame(
      Y_true=Y_true,
      Y_no_model=Y_no_model,
      Y_hat_node_nonetwork=Y_hat_node_nonetwork,
      Y_hat_node=Y_hat_node,
      Y_hat_edge=Y_hat_edge,
      Y_hat_egoalter=Y_hat_egoalter,
      Y_hat_egoalter_passed=Y_hat_egoalter_passed,
      Y_hat_egoalter_hardmode=Y_hat_egoalter_hardmode
    )
  )
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
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X)),
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
      outdeg = centrality_degree(mode='out'),
      Y_prob_true = 1 / (1 + exp(- 2 * X)),
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
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          mean(as_tibble(neighborhood, active='nodes')$outdeg)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
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
      Z = map_local_dbl(
        .f = function(neighborhood, ...) {
          mean(as_tibble(neighborhood, active='nodes')$outdeg)
        }
      ),
      Z = (Z - mean(Z)) / sd(Z),
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

######## Unobs sims ###########################################################

fn_dgp_unobs_node = function() {
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


fn_dgp_unobs_edge = function() {
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
      Q_prob = 1 / (1 + exp(1.2 - 0.05 * outdeg - 0.2 * X)),
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
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
      Z_ego = .N()$Z[from],
      Z_nbr = .N()$Z[to],
      # sample at edge level
      # probably we see the true value
      # 1 / (1 + exp(1.6 - 0.03 * outdeg - 0.2 * X))
      Q_prob = 1 / (1 + exp(
        12.0 -
        0.05 * outdeg_ego +
        0.05 * indeg_nbr -
        0.2 * X_ego - 
        0.2 * X_nbr
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
  df_edges_perf = df_edges %>%
    group_by(from) %>%
    summarize(
      gt_ego=mean(gt_ego),
      ego_preds=mean(ego_preds),
      Y_ego=mean(Y_ego)
    )
  egoalter_perf = fn_compute_performance_metrics(
    df_edges_perf %>%
      filter(gt_ego == 0) %>%
      .$Y_ego,
    df_edges_perf %>%
      filter(gt_ego == 0) %>%
      .$ego_preds
  ) %>%
    mutate(
      model='egoalter'
    )

  # hard mode egoalter perf: randomly remove some ground truth cases
  ego_mod_hardmode = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1) %>%
      sample_frac(0.5)
  )
  df_edges$ego_preds_hardmode = predict(
    ego_mod_hardmode,
    newdata=df_edges,
    type='response',
  )
  df_edges_perf_hardmode = df_edges %>%
    group_by(from) %>%
    summarize(
      gt_ego=mean(gt_ego),
      ego_preds_hardmode=mean(ego_preds_hardmode),
      Y_ego=mean(Y_ego)
    )
  egoalter_hardmode_perf = fn_compute_performance_metrics(
    df_edges_perf_hardmode %>%
      filter(gt_ego == 0) %>%
      .$Y_ego,
    df_edges_perf_hardmode %>%
      filter(gt_ego == 0) %>%
      .$ego_preds_hardmode
  ) %>%
    mutate(
      model='egoalter_hardmode'
    )

  return(
    bind_rows(
      node_perf,
      node_nonetwork_perf,
      egoalter_perf,
      egoalter_hardmode_perf
    )
  )
}


######## New graph sims ########################################################


gen_undir_pref_attach_homophily_graph = function(
  n_nodes,
  edges_per_node,
  majority_group_frac,
  same_grp_coef
) {
  
  majority_N = round(n_nodes * majority_group_frac)
  majority_Y = rep(1, majority_N)
  majority_X = rnorm(majority_N, 1, 1)
  
  minority_N = round(n_nodes * (1 - majority_group_frac))
  minority_Y = rep(0, minority_N)
  minority_X = rnorm(minority_N, -1, 1)
  
  # Group
  Y = c(majority_Y, minority_Y)
  
  X = c(majority_X, minority_X)
  
  nodes = 1:n_nodes
  
  node_df = data.frame(
    Y=Y,
    X=X
  )
  
  # Degree
  D = list()
  for (node in nodes) {
    D[[node]] = 0
  }
  
  # seed nodes are chose at random and completely connected (clique)
  # edges_per_node is also the number of seeds
  seeds = sample(nodes, edges_per_node, replace=FALSE)
  
  for (seed in seeds) {
    D[[seed]] = edges_per_node
  }
  
  edge_df = expand.grid(seeds, seeds)
  colnames(edge_df) = c('from', 'to')
  
  for (node in setdiff(nodes, seeds)) {
    node_grp = Y[node]
    for (edge_to_add in 1:edges_per_node) {
      same_grp_nodes = as.integer(Y == node_grp)
      
      # Determining a connection is a function of
      # degree, group Y
      
      # we first take the standard pref attachment probs, turn into
      # logits, and then add some factor for same group
      logits = log(
        (unlist(D) + 0.00001) / sum(unlist(D) + 0.00001)
      ) + same_grp_coef * same_grp_nodes
      probs = 1 / (1 + exp(-logits))
      probs = probs / sum(probs)
      
      draw = sample(nodes, 1, replace=FALSE, prob=probs)
      while (draw == node) {
        # if you sample yourself, redraw
        draw = sample(nodes, 1, replace=FALSE, prob=probs)
      }
      
      edge_df = bind_rows(
        edge_df,
        data.frame(from=c(node, draw), to=c(draw, node))
      )
      D[[node]] = D[[node]] + 1
      D[[draw]] = D[[draw]] + 1
    }
  }
  return(
    as_tbl_graph(edge_df) %>%
      activate(nodes) %>%
      mutate(
        Y=Y[as.integer(name)],
        X=X[as.integer(name)]
      )
  )
}


fn_dgp_new_graph_node = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = gen_undir_pref_attach_homophily_graph(
    N_NODES,
    EDGES_PER_NEW_NODE,
    0.8,
    0.5
  )
  g = g %>%
    activate(nodes) %>%
    mutate(
      outdeg = centrality_degree(mode='out'),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
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


fn_dgp_new_graph_edge = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = gen_undir_pref_attach_homophily_graph(
    N_NODES,
    EDGES_PER_NEW_NODE,
    0.8,
    0.5
  )
  g = g %>%
    activate(nodes) %>%
    mutate(
      outdeg = centrality_degree(mode='out'),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
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