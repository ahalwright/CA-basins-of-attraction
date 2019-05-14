#  Implement and apply a regular CA rule for an elementary periodic CA where states are represented
#  as a list of bits which is represented by the binary representation of an integer.  
#  https://en.wikipedia.org/wiki/Cellular_automaton  gives examples of rule 30 and rule 110.

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
#> n = 3    # rule width 3
#> w = 8    # state width 8
#> test_state = 0x00b8    # Equivalent to convert(STATE_TYPE,184)
#> R = 30   # one of the standard rules
#> r = rule_construct(n,R)
# 0x0000
# 0x0001
# 0x0001
# 0x0001
# 0x0001
# 0x0000
# 0x0000
# 0x0000
#> rule_apply(n,w,r,test_state)
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
    #print("m: ",m)
  end
  r
end

# Applies rule r to state s
# Returns the new state
function apply_rule( n::Int64, w::Int64, r::Array{STATE_TYPE}{1}, s::STATE_TYPE )
  one = convert(STATE_TYPE,1)      
  mask = one
  for j = 2:n
    mask = (mask << one) | one
  end
  #@printf("mask: %4X\n",mask)  
  es = extended_state( s, n, w )
  ns = convert(STATE_TYPE,0)    # new state
  #@printf("mask: %4X  es: %4X  ns: %0X\n",mask,es,ns)
  shift = 0
  for i = 1:w
    j = (es & (mask<<shift)) >>> shift
    #println("i:",i,"  j: ",j,"  r[j+1]: ",r[j+1],"  r[j+1]<<shift:",r[j+1]<<shift)
    ns |= r[j+1] << shift
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
  
  

