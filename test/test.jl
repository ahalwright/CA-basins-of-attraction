# Run without loading the Basins package
#   To load the Basins package, include ../src/run.jl
include("../src/types.jl")
include("../src/rule.jl")
include("../src/basins.jl")
n=3; w=4; R = 30
(basin_counts,attractor_lengths) = find_all_basins(n,w,R)
println("basin_counts: ",basin_counts,"  attractor_lengths: ",attractor_lengths)
result_df = analyze_all_rules(n,w)
println("result_df: ",result_df)

