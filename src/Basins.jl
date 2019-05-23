module Basins
#export init, create_B, find_attractor, create_preimage_dict, preimage_list, preimage_tree, find_basin
#export create_rule, apply_rule
include("../src/types.jl")
include("../src/rule.jl")
include("../src/basins.jl")
end
using Basins
