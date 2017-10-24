#' Constructs a dataframe identifying targets
#' Use drake to check dependency structure

library(drake)
dir.create("arrayOutput")

# set.seed(4387943)
# sample(1e9, 4)
# [1] 882966155 304646912 197903977 404958566

my_plan <- plan(
  inits1_1024 = initialize_chain(882966155, "stickBreaking", 20000, 20000, 100000),
  inits2_1024 = initialize_chain(304646912, "stickBreaking", 20000, 20000, 100000),
  inits3_1024 = initialize_chain(197903977, "stickBreaking", 20000, 20000, 100000),
  inits4_1024 = initialize_chain(404958566, "stickBreaking", 20000, 20000, 100000),
  chain1_sb_1024 = sample_bnp_model(inits1_1024),
  chain2_sb_1024 = sample_bnp_model(inits2_1024),
  chain3_sb_1024 = sample_bnp_model(inits3_1024),
  chain4_sb_1024 = sample_bnp_model(inits4_1024),
  strings_in_dots = "literals"
)

check(my_plan)
# plot_graph(my_plan)
save(my_plan, file = "targets.RData")