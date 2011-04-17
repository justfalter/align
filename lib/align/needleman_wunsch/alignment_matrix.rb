require 'align/needleman_wunsch/constants'
require 'align/needleman_wunsch/alignment_matrix_helpers'

module Align
  module NeedlemanWunsch
    # Obviously, this would be loads faster if this were done in C. But, I'm
    # leaving it as is, for now, for portability reasons.
    class AlignmentMatrix
      include Align::NeedlemanWunsch::Constants
      include Align::NeedlemanWunsch::AlignmentMatrixHelpers
      #extend Align::NeedlemanWunsch::AlignmentMatrixHelpers

      attr_reader :rows, :cols, :gap_penalty, :score_proc, :select_alignment_proc

      # @param [#[], #size] seq1 The first sequence
      # @param [#[], #size] seq2 The second sequence
      # @param [Hash] opts Options
      # @option opts [Fixnum] :gap_penalty (0) The gap penalty
      # @option opts [Proc] :score_proc A proc that compares two sequence items,
      #  and returns a score. Defaults to :default_score_proc
      # @option opts [Proc] :select_alignment_proc  A proc that allows for the
      #  controlling of which way to go about aligning if there are multiple
      #  paths available. Defaults to default_select_alignment_proc
      def initialize(seq1, seq2, opts = {})
        @cols = seq1.size + 1
        @rows = seq2.size + 1
        @gap_penalty = opts[:gap_penalty] || 0
        @score_proc = opts[:score_proc] || method(:default_score_proc)
        @select_alignment_proc = opts[:select_alignment_proc] || 
          method(:default_select_alignment_proc)

        @matrix = Array.new(@cols) do 
          Array.new(@rows)
        end

        0.upto(@cols-1) {|i| @matrix[i][0] = create_cell(0)}
        0.upto(@rows-1) {|j| @matrix[0][j] = create_cell(0)}

        fill(seq1, seq2)
      end

      # Returns the score matrix
      def to_score_matrix
        @matrix.map do |col|
          col.map do |row_item|
            row_item >> CELL_FLAG_BITS
          end
        end
      end

      # Returns the score matrix
      def to_traceback_matrix
        @matrix.map do |col|
          col.map do |row_item|
            row_item & CELL_FLAG_MASK
          end
        end
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
      def fill(seq1, seq2)
        1.upto(@cols-1) do |i|
          1.upto(@rows-1) do |j|
            
            score_align = parse_score_from_cell(@matrix[i-1][j-1]) + @score_proc.call(seq1[i-1], seq2[j-1])
            score_shift2 = parse_score_from_cell(@matrix[i-1][j]) + @gap_penalty
            score_shift1 = parse_score_from_cell(@matrix[i][j-1]) + @gap_penalty
            @matrix[i][j] = create_cell(score_align, score_shift1, score_shift2)
          end
        end
      end # fill

      # @return [Fixnum] The global alignment score.
      def score
        self[@cols - 1, @rows - 1] >> CELL_FLAG_BITS
      end

      # Traces backward, finding the alignment. This method is influenced by the
      # [#select_alignment_proc].
      # @return [Array] An array of tuples, pointing to coordinates in the
      #  Matrix. Note that for use in constructing the alignment of the provided
      #  sequences, you will want to substract 1 from each value. This is because
      #  the matrix is offset by 1 to account for the alignment penalty.
      def traceback
        i = @cols - 1
        j = @rows - 1

        trace = []

        last_move = CELL_FLAG_ALIGN

        while (i > 0 && j > 0)
          trace << [i,j, last_move]
          yield(i,j, last_move) if block_given?

          last_move = @select_alignment_proc.call(@matrix[i][j])
          case last_move
          when CELL_FLAG_SHIFT2; i -= 1
          when CELL_FLAG_SHIFT1;   j -= 1
          when CELL_FLAG_ALIGN; i -= 1; j -= 1
          else
            raise "invalid return from select_alignment_proc: #{last_move.inspect}"
          end
        end # while

        while i > 0
          trace << [i,j,last_move]
          yield(i,j,last_move) if block_given?
          last_move = CELL_FLAG_SHIFT2 
          i-=1
        end

        while j > 0
          trace << [i,j,last_move]
          yield(i,j,last_move) if block_given?
          last_move = CELL_FLAG_SHIFT1
          j-=1
        end

        trace
      end # traceback

    end # AlignmentMatrix
  end
end
