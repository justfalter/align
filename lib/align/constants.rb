module Align
  module Constants
    CELL_FLAG_BITS = 8
    CELL_FLAG_MASK = (1 << CELL_FLAG_BITS) - 1

    CELL_FLAG_ALIGN  = (1 << 0)
    CELL_FLAG_SHIFT1 = (1 << 1)
    CELL_FLAG_SHIFT2 = (1 << 2)
  end
end
