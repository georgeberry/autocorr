library(tidyverse)
library(tidygraph)
library(keras)

# X = observed covariate
# Z = latent covariate
# Y = attribute to predict

# Algorithm
# - Pick fraction of Y = 1
# - Draw Y's
# - Draw X | Y
# - For each node, add link to another node probability based on edges and group

N_SIMS = 100
N_NODES = 4000
NODE_SAMP_FRAC = 0.2
EDGE_SAMP_FRAC = 0.025
POWERLAW_EXPONENT = 0.8
EDGES_PER_NEW_NODE = 5
MAJORITY_GROUP_FRAC = 0.8
SAME_GRP_COEF = 0.8
NOISE_COEF = 0.2


g = play_undir_powerlaw_homophily_graph(
    N_NODES,
    EDGES_PER_NEW_NODE,
    MAJORITY_GROUP_FRAC,
    POWERLAW_EXPONENT,
    SAME_GRP_COEF,
    NOISE_COEF
)




# sanity checks: log log plot
df = g %>%
    activate(nodes) %>%
    mutate(deg=centrality_degree()) %>%
    as_tibble()

df %>% group_by(deg) %>% summarize(n=n()) %>% ggplot() + aes(x=log(deg), y=log(n)) + geom_point()

# sanity check: homophily

df_node = g %>%
    activate(nodes) %>%
    as_tibble()

df_edge = g %>%
    activate(edges) %>%
    mutate(
        Y_ego = .N()$Y[from],
        Y_nbr = .N()$Y[to],
        YY = as.integer(Y_ego == 1 & Y_nbr == 1),
        Z_ego = .N()$Z[from],
        Z_nbr = .N()$Z[to]
    ) %>%
    as_tibble()

# frac of edges origining in Y = 1 which end up in Y = 1
(mean(df_edge$YY) / mean(df_edge$Y_ego) - MAJORITY_GROUP_FRAC) / (1 - MAJORITY_GROUP_FRAC)

cor(df_edge$Z_ego, df_edge$Z_nbr)



# for node embs

el = g %>%
  activate(edges) %>%
  as_tibble()

write_delim(
    el,
    '/Users/georgeberry/Dropbox/project-autocorr/el.edgelist',
    delim=' ',
    col_names=F
)

# do embedding #

embs = read_delim(
    '/Users/georgeberry/Dropbox/project-autocorr/tmp.embeddings',
    delim=' ',
    col_names=c(
        'idx',
        'X1', 'X2', 'X3', 'X4'
    ),
    skip=1
) %>%
    mutate(idx=as.character(idx))


g = g %>%
    activate(nodes) %>%
    mutate(
        outdeg = centrality_degree(mode='out'),
        outdeg_inv = 1 / outdeg,
        indeg = centrality_degree(mode='in'),
        indeg_inv = 1 / indeg,
        gt = rbinom(n(), 1, NODE_SAMP_FRAC),
    ) %>%
    left_join(
        embs,
        by=c('name'='idx')
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
        X1_ego = .N()$X1[from],
        X2_ego = .N()$X2[from],
        X3_ego = .N()$X3[from],
        X4_ego = .N()$X4[from],
        X1_nbr = .N()$X1[to],
        X2_nbr = .N()$X2[to],
        X3_nbr = .N()$X3[to],
        X4_nbr = .N()$X4[to],
    )

df_edges = g %>% activate(edges) %>% as_tibble()


ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + log1p(outdeg_ego) + log1p(indeg_nbr) + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
        filter(gt_ego == 1)
)
df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
nbr_mod = glm(
    Y_nbr ~ X_ego + X_nbr + log1p(outdeg_ego) + log1p(indeg_nbr) + outdeg_inv_ego,
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
print(Y_hat_egoalter)

Y_true = sum(df_edges$outdeg_inv_ego * df_edges$Y)
print(Y_true)

ego_mod = glm(
    Y_ego ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + X1_ego + X2_ego + X3_ego + X4_ego + X1_nbr + X2_nbr + X3_nbr + X4_nbr + outdeg_inv_ego,
    family='binomial',
    data=df_edges %>%
        filter(gt_ego == 1)
)
df_edges$Y_ego_hat = predict(ego_mod, newdata=df_edges, type='response')
nbr_mod = glm(
    Y_nbr ~ X_ego + X_nbr + outdeg_ego + indeg_nbr + X1_ego + X2_ego + X3_ego + X4_ego + X1_nbr + X2_nbr + X3_nbr + X4_nbr + outdeg_inv_ego,
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
print(Y_hat_egoalter)


X = df_edges %>%
    select(indeg_ego:outdeg_inv_nbr, X_ego:X4_nbr) %>%
    as.matrix()
X_train_ego = df_edges %>%
    filter(gt_ego == 1) %>%
    select(indeg_ego:outdeg_inv_nbr, X_ego:X4_nbr) %>%
    as.matrix()
X_train_nbr = df_edges %>%
    filter(gt_nbr == 1) %>%
    select(indeg_ego:outdeg_inv_nbr, X_ego:X4_nbr) %>%
    as.matrix()
y_ego = to_categorical(df_edges %>% filter(gt_ego == 1) %>% .$Y_ego)
y_nbr = to_categorical(df_edges %>% filter(gt_nbr == 1) %>% .$Y_nbr)


mdl = keras_model_sequential() 
mdl %>%
    layer_dense(units=8) %>%
    layer_dense(units=8) %>%
    layer_dense(units=2, activation='sigmoid')

mdl %>% compile(
    loss = 'binary_crossentropy',
    optimizer = optimizer_adam(lr=0.01),
    metrics = c('accuracy')
)

hist = mdl %>% fit(
    X_train_ego, y_ego, epochs=10, batch_size=32
)

y_hat_ego = predict(mdl, X)[,2]

mdl = keras_model_sequential() 
mdl %>%
    layer_dense(units=8) %>%
    layer_dense(units=8) %>%
    layer_dense(units=2, activation='sigmoid')

mdl %>% compile(
    loss = 'binary_crossentropy',
    optimizer = optimizer_adam(lr=0.01),
    metrics = c('accuracy')
)

hist = mdl %>% fit(
    X_train_nbr, y_nbr, epochs=10, batch_size=32
)

y_hat_nbr = predict(mdl, X)[,2]

sum(
    df_edges$outdeg_inv_ego *
        y_hat_ego * y_hat_nbr
)