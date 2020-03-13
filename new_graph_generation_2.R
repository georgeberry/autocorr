library(tidyverse)
library(tidygraph)

# X = observed covariate
# Z = latent covariate
# Y = attribute to predict

# Algorithm
# - Pick fraction of Y = 1
# - Draw Y's
# - Draw X | Y
# - For each node, add link to another node probability based on edges and group

n_nodes = 500
edges_per_node = 5
majority_group_frac = 0.8
same_grp_coef = 0.5

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
                Y=Y[name],
                X=X[name]
            )
    )
}


g = gen_undir_pref_attach_homophily_graph(
    n_nodes,
    edges_per_node,
    majority_group_frac,
    same_grp_coef
)


# sanity checks: log log plot
df = g %>%
    activate(nodes) %>%
    mutate(deg=centrality_degree()) %>%
    as_tibble()

df %>% group_by(deg) %>% summarize(n=n()) %>% ggplot() + aes(x=log(deg), y=log(n)) + geom_point()

# sanity check: homophily

df_edge = g %>%
    activate(edges) %>%
    mutate(
        Y_ego = .N()$Y[from],
        Y_nbr = .N()$Y[to],
        YY = as.integer(Y_ego == 1 & Y_nbr == 1)
    ) %>%
    as_tibble()

H = (mean(df_edge$YY) / mean(df_edge$Y_ego) - majority_group_frac) / (1 - majority_group_frac)

