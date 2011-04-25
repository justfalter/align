require 'align/basic_scoring'
require 'align/pairwise_algorithm'

module Align
  # Align two sequences via [NeedlemanWunsch.align]
  # References:
  # [http://www.avatar.se/molbioinfo2001/dynprog/dynamic.html]
  class NeedlemanWunsch < PairwiseAlgorithm
    attr_reader :cols, :rows, :matrix

    # Default scoring for 
    SCORING_DEFAULT = Align::BasicScoring.new(1,0,0)
    SCORING_ALT1    = Align::BasicScoring.new(1,-1,-1)

    # @param [#[], #size] seq1 The first sequence
    # @param [#[], #size] seq2 The second sequence
    # @param [Hash] opts Options
    # @option opts [NeedlemanWunschScoring] :scoring (NeedlemanWunschScoring) An instance of a scoring object.
    # @option opts [Object] :skip_obj (nil) An object to shove into the gaps of
    #  the aligned sequences
    def initialize(seq1, seq2, opts = {})
      super(seq1, seq2, opts[:scoring] || SCORING_DEFAULT)

      @highest_score = nil
      @highest_score_loc = nil

      @rows = @seq1.size + 1
      @cols = @seq2.size + 1

      @skip_obj = opts[:skip_obj] || nil

      @matrix = Array.new(@rows) do 
        Array.new(@cols)
      end

      fill()
    end

    def score
      @matrix[@rows-1][@cols-1]
    end

    # Fills the matrix with the alignment map.
    def fill
      @matrix[0][0] = 0
      # Set up the first column on each row.
      1.upto(@rows-1) {|i| @matrix[i][0] = @matrix[i-1][0] + @scoring.score_delete(@seq1[i])}
      # Set up the first row 
      1.upto(@cols-1) {|j| @matrix[0][j] = @matrix[0][j-1] + @scoring.score_insert(@seq2[j])}

      1.upto(@rows-1) do |i|
        prv_row = @matrix[i-1]
        cur_row = @matrix[i]

        1.upto(@cols-1) do |j|
          
          seq1_obj = @seq1[i-1]
          seq2_obj = @seq2[j-1]

          # Calculate the score.
          score_align = prv_row[j-1] + @scoring.score_align(seq1_obj, seq2_obj)
          score_delete = prv_row[j] + @scoring.score_delete(seq1_obj)
          score_insert = cur_row[j-1] + @scoring.score_insert(seq2_obj)
          max = max3(score_align, score_delete, score_insert)

          @matrix[i][j] = max
        end
      end
    end # fill

    # Traces backward, finding the alignment.
    # @yield [i,j,step] 
    # @yieldparam i [Integer] The location in sequence one
    # @yieldparam j [Integer] The location in sequence two
    # @yieldparam step [Integer] The direction we took
    def traceback
      i = @rows - 1
      j = @cols - 1

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
