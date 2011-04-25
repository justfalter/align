module Align
  # Basic Scoring interface
  class BasicScoring
    # @param [Numeric] align_match Price for alignment.
    # @param [Numeric] align_mismatch Penalty for misalignment
    # @param [Numeric] gap_penalty Gap penalty for insert/delete
    def initialize(align_match, align_mismatch, gap_penalty)
      @align_match = align_match
      @align_mismatch = align_mismatch
      @gap_penalty = gap_penalty
    end

    def score_align(a,b)
      (a == b) ? @align_match : @align_mismatch
    end

    def score_insert(a)
      @gap_penalty
    end

    def score_delete(a)
      @gap_penalty
    end
  end
end
