#' Constructs a dataframe identifying targets
#' Use drake to check dependency structure

library(drake)
dir.create("arrayOutput")

# set.seed(4387943)
# sample(1e9, 4)
# [1] 882966155 304646912 197903977 404958566

my_plan <- plan(
  inits1 = initialize_chain(882966155, "stickBreaking"),
  inits2 = initialize_chain(304646912, "stickBreaking"),
  inits3 = initialize_chain(197903977, "stickBreaking"),
  inits4 = initialize_chain(404958566, "stickBreaking"),
  chain1_sb = sample_bnp_model(inits1),
  chain2_sb = sample_bnp_model(inits2),
  chain3_sb = sample_bnp_model(inits3),
  chain4_sb = sample_bnp_model(inits4),
  strings_in_dots = "literals"
)

check(my_plan)
# plot_graph(my_plan)
save(my_plan, file = "targets.RData")