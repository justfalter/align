require 'align/scoring_interface'
require 'align/pairwise_algorithm'

module Align
  # Default scoring for the Needleman-Wunsch algorithm.
  class NeedlemanWunschScoring < ScoringInterface
    MATCH = 1
    NOMATCH = 0
    GAP_PENALTY = 0

    def score_align(a,b)
      (a == b) ? MATCH : NOMATCH
    end

    def score_insert(a)
      GAP_PENALTY
    end

    def score_delete(a)
      GAP_PENALTY
    end
  end

  # Align two sequences via [NeedlemanWunsch.align]
  # References:
  # [http://www.avatar.se/molbioinfo2001/dynprog/dynamic.html]
  class NeedlemanWunsch < PairwiseAlgorithm
    attr_reader :highest_score, :highest_score_loc
    attr_reader :rows, :cols

    DEFAULT_SCORING = NeedlemanWunschScoring.new

    # @param [#[], #size] seq1 The first sequence
    # @param [#[], #size] seq2 The second sequence
    # @param [Hash] opts Options
    # @option opts [NeedlemanWunschScoring] :scoring (NeedlemanWunschScoring) An instance of a scoring object.
    # @option opts [Object] :skip_obj (nil) An object to shove into the gaps of
    #  the aligned sequences
    def initialize(seq1, seq2, opts = {})
      super(seq1, seq2, opts[:scoring] || DEFAULT_SCORING)

      @highest_score = nil
      @highest_score_loc = nil

      @cols = @seq1.size + 1
      @rows = @seq2.size + 1

      @skip_obj = opts[:skip_obj] || nil

      @matrix = Array.new(@cols) do 
        Array.new(@rows)
      end

      fill()
    end

    # Returns the score matrix # @return [Array<Array<Fixnum>>] A matrix where each cell is the score.
    def to_score_matrix
      @matrix
    end

    # Returns the value at col, row
    # @param [Fixnum] col
    # @param [Fixnum] row
    # @return [Object]
    def [](col, row)
      if col >= @cols || row >= @rows
        raise ArgumentError.new("out of bounds (col: #{col} >= #{@cols} || row: #{row} >= #{@rows})")
      end
      @matrix[col][row]
    end

    # Fills the matrix with the alignment map.
    def fill
      0.upto(@cols-1) {|i| @matrix[i][0] = 0}
      @matrix[0].fill(0)

      1.upto(@cols-1) do |i|
        prv_col = @matrix[i-1]
        cur_col = @matrix[i]

        1.upto(@rows-1) do |j|
          
          seq1_obj = @seq1[i-1]
          seq2_obj = @seq2[j-1]

          # Calculate the score.
          score_align = prv_col[j-1] + @scoring.score_align(seq1_obj, seq2_obj)
          score_delete = prv_col[j] + @scoring.score_delete(seq1_obj)
          score_insert = cur_col[j-1] + @scoring.score_insert(seq2_obj)
          max = max3(score_align, score_delete, score_insert)

          # Store the highest score and where we've seen it.
          if @highest_score.nil? || max >= @highest_score
            @highest_score = max
            @highest_score_loc = [i,j]
          end

          @matrix[i][j] = max
        end
      end
    end # fill

    # Traces backward, finding the alignment.
    # @param [Integer] i The column to start tracing back from
    # @param [Integer] j The row to start tracing back from
    # @yield [i,j,step] 
    # @yieldparam i [Integer] The location in sequence one
    # @yieldparam j [Integer] The location in sequence two
    # @yieldparam step [Integer] The direction we took
    def traceback
      i = @highest_score_loc[0]
      j = @highest_score_loc[1]

      while (i > 0 && j > 0)
        score = @matrix[i][j]

        seq1_obj = @seq1[i-1]
        seq2_obj = @seq2[j-1]

        score_align = @matrix[i-1][j-1] + @scoring.score_align(seq1_obj, seq2_obj)
        score_delete = @matrix[i-1][j] + @scoring.score_delete(seq1_obj)
        score_insert = @matrix[i][j-1] + @scoring.score_insert(seq2_obj)

        flags = 0
        need_select = false

        if score == score_align
          flags = :align
          i-=1
          j-=1
        elsif score == score_delete
          flags = :delete
          i-=1
        else
          flags = :insert
          j-=1
        end

        yield(i,j,flags) 
      end # while

      while i > 0
        i-=1
        yield(i,j,:delete) 
      end

      while j > 0
        j-=1
        yield(i,j,:insert) 
      end
    end # traceback

    # Like traceback, but returns an array of the traceback instead of 
    # yielding blocks.
    def traceback_array
      trace = []
      traceback do |i,j,flags|
        trace << [i,j,flags]
      end
      trace
    end # traceback_array

    # Returns the sequences in aligned arrays. Gaps are filled with :skip_obj
    # @return Two arrays containing the sequences, and their elements.
    # @param [Integer] i The column to start tracing back from
    # @param [Integer] j The row to start tracing back from
    def align
      alignment_1 = []
      alignment_2 = []

      traceback do |i, j, flags|
        seq1_val = seq2_val = @skip_obj
        case flags
        when :align
          seq1_val = @seq1[i]
          seq2_val = @seq2[j]
        when :insert
          seq2_val = @seq2[j]
        when :delete
          seq1_val = @seq1[i]
        end
        alignment_1.unshift seq1_val
        alignment_2.unshift seq2_val
      end

      [alignment_1, alignment_2]
    end

    # Aligns two sequences together.
    # @param see (#initialize)
    def self.align(seq1, seq2, opts = {})
      self.new(seq1, seq2, opts).align
    end

  end
end
