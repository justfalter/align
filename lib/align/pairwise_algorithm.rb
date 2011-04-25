module Align
  # Provides a base for algorithms that align two sequences.
  class PairwiseAlgorithm
    attr_reader :seq1, :seq2, :scoring

    def initialize(seq1, seq2, scoring)
      @seq1 = seq1
      @seq2 = seq2
      @scoring = scoring
    end

    # Determines the maximum value of three variables. 3-4 times faster than
    # [v1,v2,v3].max.
    def max3(v1,v2,v3)
      (v1 >= v2) ? ((v1 >= v3)? v1 : v3) : ((v2 >= v3)? v2 : v3)
    end
  end
end
