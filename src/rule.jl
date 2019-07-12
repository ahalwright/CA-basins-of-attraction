#  Implement and apply a regular CA rule for an elementary periodic CA where states are represented
#  as a list of bits which is represented by the binary representation of an integer.  
#  https://en.wikipedia.org/wiki/Cellular_automaton  gives examples of rule 30 and rule 110.
using Printf
#include("types.jl")

export create_rule, apply_rule

# Assumptions:  
#   binary:  states are a bit string represented as an unsigned integer
#   periodic boundary conditions
#   regular:  rule widths are 3, 5, 7 and always centered at the state component being updated
# Notation:
#   w   width (number of bits) of a state
#   n   rule width (should be 3, 5, 7, etc.)


#Rules are functions from the range 0:2^w-1 to {0,1}.
#These functions are represented by an array of {0,1}  indexed by 0:2^w-1

# Example:  
#> include("types.jl")
#> include("rule.jl")
#> include("basins.jl")
#> n = 3    # rule width 3
#> w = 8    # state width 8
#> test_state = 0x00b8    # Equivalent to convert(STATE_TYPE,184)
#> R = 30   # one of the standard rules
#> r = create_rule(n,R)
# 0x0000
# 0x0001
# 0x0001
# 0x0001
# 0x0001
# 0x0000
# 0x0000
# 0x0000
#> apply_rule(n,w,r,test_state)
# 0x00a5

# Constructs a rule of width n whose binary representation is the binary representation of R
# R is the standard rule numberin used by Wolfram and Wuensche & Lesser.
# The result is a rule represented by an array of {0,1}  indexed by 0:2^w-1
function create_rule( n::Int64, R::Int64 )
  r = zeros(STATE_TYPE,2^n)
  m = convert(STATE_TYPE,1)
  shift = 0
  for i = 1:2^n
    r[i] = (R & m) >> shift
    m <<= 1
    shift += 1
    #println("m:",m)
  end
  r
end

# Applies rule r to state s
# Returns the new state
function apply_rule( n::Int64, w::Int64, r::Array{STATE_TYPE}{1}, s::STATE_TYPE )
  if n == 0    # a non-cellular-automaton rule
    return nonCA_apply_rule(n,w,r,s)
  end
  one = convert(STATE_TYPE,1)      
  mask = one
  for j = 2:n
    mask = (mask << one) | one
  end
  es = extended_state( s, n, w )  
  ns = convert(STATE_TYPE,0)    # new state is all zeros
  #@printf("state:%3X  mask:%3X  es:%3X  ns: %0X\n",s,mask,es,ns)
  shift = 0
  for i = 1:w
    j = (es & (mask<<shift)) >>> shift  # Extract w bits of extended state, shift to low order
    #print("i:",i,"  j: ",j)
    #@printf("  r[j+1]:%3X  r[j+1]<<shift:%3X\n",r[j+1],r[j+1]<<shift)
    ns |= r[j+1] << shift  # apply a bit from rule to next state
    shift += 1
    #@printf("ns: %4X\n",ns)
  end
  ns
end

# Extends a state by replicating the first n/2 bits on the right end and the last n/2 bits on the left end
function extended_state( s::STATE_TYPE, n::Int64, w::Int64 )
  one = convert(STATE_TYPE,1)      
  mask = one
  n_2 = convert(STATE_TYPE,n >>> 1)     # n/2
  for j = 2:n_2
    mask = (mask << one) | one
  end
  center = s << n_2
  right = (center & (mask << w)) >> w
  left = (s & mask) << (w + n_2)
  es = right | center | left
  return es
end

# Example of a nonCA rule for w=2 (array format)
# r = [0x0003, 0x0002, 0x0001, 0x0000]  # flips bits
# Thus, in binary notation: r[00]=11, r[01]=10, r[10]=01, r[11]=00.
# In integer format, concatenate the components
# R = 0xE4 = 11100100 (binary) = 228 (decimal)
# Another rule is R=108 which is cyclic.

# Converts a rule from a decimal integer representation to a list representation
function nonCA_convert_rule(w::Int64,R::Int64)
  r = zeros(STATE_TYPE,2^w)
  one = convert(STATE_TYPE,1)      
  mask = one
  for j = 2:w
    mask = (mask << one) | one
  end
  @printf("mask:%3X\n",mask)
  shift = 0
  for i = 2^w:-1:1
    r[i] = mask & (R>>>shift)
    @printf("i:%d  shift:%d  R>>>shift:%3X  r[i]:%3X\n",i,shift,R>>>shift,r[i])
    shift += w
  end
  r
end
  
# Applies rule r to state s.
# Rule r is in list representation where r[s+1] is the result of applying r to state s.
# Note that Julia lists are base 1 rather than base 0, so "s+1" converts the state to base 1.
function nonCA_apply_rule( n::Int64, w::Int64, r::Array{STATE_TYPE}{1}, s::STATE_TYPE )
  if n != 0
    error("n must be zero for function non_CA_rule_apply")
  end
  @assert length(r) == 2^w
  return r[s+1]
end

function concatenate_rules( w1::Int64, w2::Int64, r1::RULE_TYPE, r2::RULE_TYPE )
  len1 = 2^w1
  len2 = 2^w2
  len = len1*len2
  println("len1: ",len1,"  len2: ",len2,"  len: ",len)
  r = zeros(STATE_TYPE,len)
  for i2 = 0:(len2-1)
    for i1 = 0:(len1-1)
      i = i2*len1 + i1
      print("i2: ",i2,"  i1: ",i1,"  i: ",i)
      r[i+1] = r2[i2+1]*len1 + r1[i1+1]
      @printf("  r2[i2+1]:%3X  r1[i1+1]:%3X  r[i+1]:%4X \n",r2[i2+1],r1[i1+1],r[i+1])
    end
  end
  r
end
