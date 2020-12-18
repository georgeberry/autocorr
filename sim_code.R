library(tidyverse)
library(ggraph)
library(tidygraph)
library(foreach)
library(doParallel)
library(pROC)
library(caret)

registerDoParallel(cores=12)
N_SIMS = 100
N_NODES = 4000
N_GROUND_TRUTH_NODES = 500
POWERLAW_EXPONENT = 0.8
EDGES_PER_NEW_NODE = 5
MAJORITY_GROUP_FRAC = 0.7
SAME_GRP_COEF = 0.7

# The basic idea is that there's a graph gen function, a function which does 
# the sampling, and then a function which fits the models. The basic structure
# is: fn_model_fitting(fn_sampling), and the fn_sampling knows how to 
# gen the graph.

#### Graph gen #################################################################

play_undir_powerlaw_homophily_graph = function(
    n_nodes,
    edges_per_node,
    majority_group_frac,
    powerlaw_exponent,
    same_grp_coef
) {

    majority_N = round(n_nodes * majority_group_frac)
    majority_Y = rep(1, majority_N)
    majority_Z = rnorm(majority_N, 0, 1)
    majority_X = 1 + majority_Z

    minority_N = round(n_nodes * (1 - majority_group_frac))
    minority_Y = rep(0, minority_N)
    minority_Z = rnorm(minority_N, 0, 1)
    minority_X = -1 + minority_Z

    # Group
    Y = c(majority_Y, minority_Y)

    X = c(majority_X, minority_X)

    Z = c(majority_Z, minority_Z)

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

    node_iter_order = sample(
      setdiff(nodes, seeds),
      length(setdiff(nodes, seeds)),
      replace=FALSE
    )

    for (node in node_iter_order) {
        node_grp = Y[node]
        same_grp_nodes = as.integer(Y == node_grp)
        for (edge_to_add in 1:edges_per_node) {
            # Determining a connection is a function of
            # degree, group_Y, Z

            # Using the conditional logit expression in Overgoor et al 2020
            deg_probs = exp(
              powerlaw_exponent * log(unlist(D)) +
              same_grp_coef * same_grp_nodes
            )
            probs = deg_probs / sum(deg_probs)

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
        as_tbl_graph(edge_df)  %>%
            activate(nodes) %>%
            mutate(
                Y=Y[as.integer(name)],
                X=X[as.integer(name)],
                Z=Z[as.integer(name)]
            )
    )
}


#### Samplers ##################################################################


fn_new_dgp_main_node_randsamp = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_undir_powerlaw_homophily_graph(
    N_NODES,
    EDGES_PER_NEW_NODE,
    MAJORITY_GROUP_FRAC,
    POWERLAW_EXPONENT,
    SAME_GRP_COEF
  )
  g = g %>%
    activate(nodes) %>%
    mutate(
      outdeg = centrality_degree(mode='out'),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
      gt = sample(
        c(rep(1, N_GROUND_TRUTH_NODES), rep(0, N_NODES - N_GROUND_TRUTH_NODES)),
        N_NODES,
        replace=FALSE
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
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
      gt = as.numeric(gt_ego & gt_nbr),
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
    )
  return(g)
}


fn_new_dgp_main_node_degsamp = function() {
  # Creates a graph and draws Y = f(X, Z), where X ~ N(0, 1) and
  # Z is the ego-network avg of the observed X
  g = play_undir_powerlaw_homophily_graph(
    N_NODES,
    EDGES_PER_NEW_NODE,
    MAJORITY_GROUP_FRAC,
    POWERLAW_EXPONENT,
    SAME_GRP_COEF
  )
  g = g %>%
    activate(nodes) %>%
    mutate(
      outdeg = centrality_degree(mode='out'),
      outdeg_inv = 1 / outdeg,
      indeg = centrality_degree(mode='in'),
      indeg_inv = 1 / indeg,
    )

  tmp = g %>%
    activate(nodes) %>%
    as_tibble() %>%
    select(name, outdeg) %>%
    mutate(pr = outdeg / sum(outdeg))
  gt_nodes = sample(
    tmp$name,
    N_GROUND_TRUTH_NODES,
    replace=FALSE,
    prob=tmp$pr
  )
    
  g = g %>%
    mutate(
      gt = as.integer(name %in% gt_nodes)
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


fn_new_dgp_main_edge_randsamp = function() {
  # The sampling here is tricky b/c we need to sample a random sample of edges
  # given a fixed node labeling budget
  # 
  # One way to do this would be to assign each edge a random number between 1
  # and E, and then set all edges to GT which are < the # which gives about
  # the right number of nodes

  g = play_undir_powerlaw_homophily_graph(
    N_NODES,
    EDGES_PER_NEW_NODE,
    MAJORITY_GROUP_FRAC,
    POWERLAW_EXPONENT,
    SAME_GRP_COEF
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
      X_ego = .N()$X[from],
      X_nbr = .N()$X[to],
      rand_num = sample(1:n(), n(), replace=FALSE),
    )

  tmp = g %>%
    activate(edges) %>%
    as_tibble() %>%
    select(from, to, rand_num)

  # literally just iterate through until the N of GT nodes is within 1 of
  # the desired
  rand_num_cutoff = NULL
  for (i in 1:nrow(tmp)) {
    tmp2 = tmp %>%
      filter(rand_num <= i)
    n_unique_nodes = length(unique(c(tmp2$from, tmp2$to)))
    if (
      # break first time we see this condition, should always be within 1 of
      # correct answer
      n_unique_nodes >= N_GROUND_TRUTH_NODES
    ) {
      rand_num_cutoff = i
      break
    }
  }
  gt_nodes = unique(c(tmp2$from, tmp2$to))

  g = g %>%
    activate(edges) %>%
    mutate(
      gt = ifelse(rand_num <= rand_num_cutoff, 1, 0)
    ) %>%
    # for any node who is part of a gt edge, give that node gt=1
    activate(nodes) %>%
    mutate(
      # if any edge containing the focal node is 1, the node is gt=1
      gt=as.integer(name %in% gt_nodes)
    ) %>%
    # now go back and label the gt_ego/gt_nbr based on the node being included
    # in any ground-truth edge
    # Note you can't actually use the other edges that by chance end up
    # in the ground truth set or it messes up the edge sample
    activate(edges) %>%
    mutate(
      gt_ego = .N()$gt[from],
      gt_nbr = .N()$gt[to],
    )
  return(g)
}



#### Helpers ###################################################################

# A simple helper to run stuff in parallel
fn_par_run = function(fn, n_reps, ...) {
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

#### Model fitting #############################################################

fn_gen_graph_and_fit_models = function(fn_g) {
  # fit 7 models
  # node (no network features)
  # node (ctrl for ego_deg_inv)
  # node (ctrl for ego_deg_inv and ego_deg)
  # edge (ctrl for X_ego, X_nbr, ego_deg_inv)
  # edge (ctrl for X_ego, X_nbr, ego_deg_inv, nbr_deg_inv, d_ego, d_nbr)
  # ego-alter (ctrl for X_ego, X_nbr, ego_deg_inv)
  # ego-alter (ctrl for X_ego, X_nbr, ego_deg_inv, nbr_deg_inv, d_ego, d_nbr)

  g = fn_g()

  #### True values in the network
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  Y_true = sum(df_edges$outdeg_inv_ego * df_edges$Y)
  Y_no_model = sum(
    df_edges %>% filter(gt == 1) %>% .$outdeg_inv_ego * 
    df_edges %>% filter(gt == 1) %>% .$Y
  ) * (nrow(df_edges) / sum(df_edges$gt))

  #### node no network features
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
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

  #### node with basic network features
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  node_mod = glm(
    Y ~ X + outdeg_inv,
    family='binomial',
    data=df_nodes %>%
      filter(gt == 1)
  )
  node_preds = predict(node_mod, newdata=df_nodes, type='response')
  Y_hat_node_basic = sum(
    df_edges$outdeg_inv_ego *
    node_preds[df_edges$from] *
    node_preds[df_edges$to]
  )

  #### node with full network features
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  node_mod = glm(
    Y ~ X + outdeg + outdeg_inv,
    family='binomial',
    data=df_nodes %>%
      filter(gt == 1)
  )
  node_preds = predict(node_mod, newdata=df_nodes, type='response')
  Y_hat_node_full = sum(
    df_edges$outdeg_inv_ego *
    node_preds[df_edges$from] *
    node_preds[df_edges$to]
  )
  
  #### Predict edge categories directly (basic features)
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  edge_mod = glm(
    Y ~ X_ego + X_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt == 1)
  )
  Y_hat_edge_basic = sum(
    df_edges$outdeg_inv_ego *
    predict(edge_mod, newdata=df_edges, type='response')
  )

  #### Predict edge categories directly (full features)
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  edge_mod = glm(
    Y ~ X_ego +
      X_nbr +
      outdeg_inv_ego +
      outdeg_inv_nbr +
      outdeg_ego +
      outdeg_nbr,
    family='binomial',
    data=df_edges %>%
      filter(gt == 1)
  )
  Y_hat_edge_full = sum(
    df_edges$outdeg_inv_ego *
    predict(edge_mod, newdata=df_edges, type='response')
  )

  #### ego-alter (basic features)
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1)
  )
  df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
  df_edges %>%
    left_join(
      df_edges %>%
        select(from, to, Y_nbr_hat=Y_ego_hat),
      by=c('from'='to', 'to'='from')
    ) ->
    df_edges
  Y_hat_egoalter_basic = sum(
    df_edges$outdeg_inv_ego *
    df_edges$Y_ego_hat *
    df_edges$Y_nbr_hat
  )

  #### ego-alter (full features)
  df_nodes = g %>% activate(nodes) %>% as_tibble()
  df_edges = g %>% activate(edges) %>% as_tibble()
  ego_mod = glm(
    Y_ego ~ X_ego +
      X_nbr +
      outdeg_inv_ego +
      outdeg_inv_nbr +
      outdeg_ego +
      outdeg_nbr,
    family='binomial',
    data=df_edges %>%
      filter(gt_ego == 1)
  )
  df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
  df_edges %>%
    left_join(
      df_edges %>%
        select(from, to, Y_nbr_hat=Y_ego_hat),
      by=c('from'='to', 'to'='from')
    ) ->
    df_edges
  Y_hat_egoalter_full = sum(
    df_edges$outdeg_inv_ego *
    df_edges$Y_ego_hat *
    df_edges$Y_nbr_hat
  )

  return(
    data.frame( 
      Y_true=Y_true,
      Y_no_model=Y_no_model,

      Y_hat_node_nonetwork=Y_hat_node_nonetwork,
      Y_hat_node_basic=Y_hat_node_basic,
      Y_hat_node_full=Y_hat_node_full,

      Y_hat_edge_basic=Y_hat_edge_basic,
      Y_hat_edge_full=Y_hat_edge_full,

      Y_hat_egoalter_basic=Y_hat_egoalter_basic,
      Y_hat_egoalter_full=Y_hat_egoalter_full
    )
  )
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

