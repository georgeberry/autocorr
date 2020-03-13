library(tidyverse)
source('/Users/georgeberry/Dropbox/project-autocorr/autocorr/sim_code.R')

g = fn_dgp_main_node()

g %>% activate(edges) %>%
  as_tibble(.) %>%
  write_csv(., '/Users/georgeberry/Dropbox/project-autocorr/autocorr/edgelist.csv')

g %>% activate(nodes) %>%
  as_tibble(.) %>%
  write_csv(., '/Users/georgeberry/Dropbox/project-autocorr/autocorr/node_attrs.csv')