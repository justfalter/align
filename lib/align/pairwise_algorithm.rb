module Align
  # Provides a base for algorithms that align two sequences.
  class PairwiseAlgorithm
    attr_reader :seq1, :seq2, :scoring

    def initialize(seq1, seq2, scoring)
      @seq1 = seq1
      @seq2 = seq2
      @scoring = scoring
    end

    # Max of 2
    def max2(a,b)
      a >= b ? a : b
    end

    # Determines the maximum value of three variables. 3-4 times faster than
    # [a,b,c].max.
    def max3(a,b,c)
      (a >= b) ? ((a >= c)? a : c) : ((b >= c)? b : c)
    end

    # Returns the max of 4 integers
    def max4(a,b,c,d)
      x = a >= b ? a : b
      y = c >= d ? c : d
      (x >= y) ? x : y
    end



    # Returns the sequences in aligned arrays. Gaps are filled with :skip_obj
    # @return Two arrays containing the sequences, and their elements.
    def align
      raise NotImplementedError.new("#{self.class}#align")
    end
  end
end
