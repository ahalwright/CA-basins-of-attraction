# Create and update the basin of attraction data structure B
#
using DataFrames
using Printf
#include("rule.jl")

export init, create_B, find_attractor, create_preimage_dict, preimage_list, preimage_tree, find_basin, find_all_basins, B, Binv_df, Binv_dict, analyze_all_rules

#=
Example
 n = 3
 w = 5
 R = 30
 B = create_B( n, w, R )
 find_attractor(0x0007,B,1)
 Binv_df, Binv_dict = create_preimage_dict(B);
 preimage_list( 0x000, Binv_df, Binv_dict )
=#

function init( n, w, R )
  global B = create_B( n, w, R )
  #find_attractor(0x0007,B,1)
  global Binv_df
  global Binv_dict
  Binv_df, Binv_dict = create_preimage_dict(B);
end

function create_B( n::Int64, w::Int64, R::Int64 )
r = create_rule(n,R)
  zero = convert(STATE_TYPE,0)
  last_state = convert(STATE_TYPE,2^w-1)
  state_list = collect(zero:last_state)
  next_state_list = map(s->apply_rule(n,w,r,s),state_list)
  B = DataFrame(state=state_list,next_state=next_state_list, basin=zeros(Int64,2^w), distance=zeros(Int64,2^w) )
end

function find_attractor(start_state::STATE_TYPE, B::DataFrame, basin_index::Int64=0 )
  state = start_state
  visited = Set{STATE_TYPE}(state)
  #next_state = apply_rule(n,w,rule,start_state)
  next_state = B[:next_state][state+1]
  while !in(next_state,visited)
    push!(visited,next_state)
    next_state = B[:next_state][next_state+1]
    #@printf("next_state: 0X%2X\n",next_state)
  end
  attractor_start_state = next_state
    @printf("attractor_start_state: 0X%2X\n",attractor_start_state)
  attractor = Array{STATE_TYPE}{1}()
  push!(attractor,attractor_start_state)
  #println("attractor: ",attractor)
  if basin_index != 0
    B[:basin][attractor_start_state+1] = basin_index
  end
  next_state = B[:next_state][next_state+1]
  while next_state != attractor_start_state
    push!(attractor,next_state)
    if basin_index != 0
      B[:basin][next_state+1] = basin_index
    end
    next_state = B[:next_state][next_state+1]
    #@printf("next_state: 0X%2X\n",next_state)
  end
  attractor
end

function create_preimage_dict(B::DataFrame)
  Binv_df = sort!(deepcopy(B),[order(:next_state),order(:state)])
  Binv_dict = Dict(Binv_df[:next_state][i]=>i for i = size(B)[1]:-1:1)
  (Binv_df,Binv_dict)
end

function preimage_list( state::STATE_TYPE, B::DataFrame, Binv_df::DataFrame, Binv_dict::Dict )
  result = Array{STATE_TYPE}{1}()
  row = get(Binv_dict,state,-1)
  if row == -1
    return result
  end
  push!(result,Binv_df[:state][row])
  row += 1
  #@printf("row: %d   Binv_df[:state][row]: %X\n",row,Binv_df[:state][row])
  while row <= size(B)[1] && Binv_df[:next_state][row] == state
    push!(result,Binv_df[:state][row])
    row += 1
  end
  result
end

"""
@doc function preimage_tree()
  Traverses the tree of preimages starting from state s
  Sets the :basin and :distance fields of B for every state in the tree.
  Returns the number of states in the tree.
"""

function preimage_tree( s::STATE_TYPE, B::DataFrame, Binv_df::DataFrame, Binv_dict::Dict )
  #count += 1
  #printf("count: ",count)
  #println("B[:state][s+1]: ",B[:state][s+1] )
  @assert B[:state][s+1] == s
  if B[:basin][s+1] == 0
    error("The basin of the start state have been set in B.")
  end
  basin = B[:basin][s+1]
  distance = B[:distance][s+1]
  #@printf("state: %X   distance: %d   ",s,distance)
  #println("basin: ",basin,"  distance: ",distance)
  preimages = preimage_list( s, B, Binv_df, Binv_dict )
  #println("preimage: ",preimages)
  if length(preimages) == 0
    return 1
  end
  if distance > 15
    return -1
  end
  count = 1
  for p in preimages
    if B[:basin][p+1] == 0
      B[:basin][p+1] = basin
      B[:distance][p+1] = distance+1
      count += preimage_tree( p, B, Binv_df, Binv_dict )
      #@printf("state: %X   distance: %d  count: %d \n",s,distance,count)
    end
  end
  #println("returned count: ",count)
  return count
end

function find_basin( start_state::STATE_TYPE, B::DataFrame, Binv_df::DataFrame, Binv_dict::Dict, basin_index::Int64 )
  attractor = find_attractor(start_state, B, basin_index)
  count = 0
  for a in attractor
    count += preimage_tree(a, B, Binv_df, Binv_dict )
  end
  #println("find_basin: len attr: ",length(attractor),"  count: ",count)
  return length(attractor),count
end

function find_all_basins( n::Int64, w::Int64, R:: Int64 )
  println("find_all_basins: n: ",n,"  w: ",w,"  R: ",R)
  B = create_B( n, w, R)
  Binv_df, Binv_dict = create_preimage_dict( B )
  zero = convert(STATE_TYPE,0)
  len = size(B)[1]
  basin_ind = 0
  basin_counts = Int64[]
  attractor_lengths = Int64[]
  for i = 1 : len
    #println("i: ",i,"  b: ",B[:basin][i],"  state: ",B[:state][i])
    if B[:basin][i] == 0
      basin_ind += 1
      (alen,count) = find_basin( B[:state][i], B, Binv_df, Binv_dict, basin_ind )
      push!(basin_counts,count)
      push!(attractor_lengths,alen)
      #println("basin_ind: ",basin_ind,"  alen: ",alen,"  count: ",count)
    end
  end
  #println("sum(basin_counts): ",sum(basin_counts)) 
  #println("attractor_lengths: ",attractor_lengths)
  (basin_counts,attractor_lengths)
end
    
function analyze_all_rules(n::Int64, w::Int64 )
  len = 2^2^n   # number of rules
  result_df = DataFrame( rule=collect(0:(len-1)), basin_sizes=Vector{Int64}[Int64[] for i = 1:len], 
      attractor_lengths=Vector{Int64}[Int64[] for i = 1:len] )
  for R = 0:len-1
    println("R: ",R)
    basin_counts, attractor_lengths = find_all_basins( n, w, R)
    println("basin_counts: ",basin_counts,"  attractor_lengths: ",attractor_lengths)
    result_df[:basin_sizes][R+1] = basin_counts
    result_df[:attractor_lengths][R+1] = attractor_lengths
  end
  result_df
end
