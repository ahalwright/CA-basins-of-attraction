# States are represented as the binary representation of unsigned integers
# The number of bits in the representation must be at least w + n - 1
export STATE_TYPE, RULE_TYPE
const STATE_TYPE = UInt16
const RULE_TYPE = Array{STATE_TYPE}{1}

