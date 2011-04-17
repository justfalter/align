require 'align/needleman_wunsch/constants'

module Align
  module NeedlemanWunsch
    # Helpers for working with the AlignmentMatrix
    module AlignmentMatrixHelpers
      include Align::NeedlemanWunsch::Constants

      # Creates a cell for the alignment matrix with the purpose of tracking
      # the score and source of the score. 
      # @param [Fixnum] align The calculated diag score.
      # @param [Fixnum] shift1 The calculated 'up' score
      # @param [Fixnum] shift2 The calculated 'left' score
      # @return [Fixnum] The score, with flags set on the lower CELL_FLAG_BITS
      # bits.
      def create_cell(align=0, shift1=0, shift2=0)
        max = [align,shift1,shift2].max
        val = max << CELL_FLAG_BITS
        val |= CELL_FLAG_ALIGN  if align  == max
        val |= CELL_FLAG_SHIFT1 if shift1 == max
        val |= CELL_FLAG_SHIFT2 if shift2 == max
        val
      end

      # Parses the score from a cell value
      # @param [Fixnum] cell The value in the cell
      # @return [Fixnum] The score for that cell.
      def parse_score_from_cell(cell)
        cell >> CELL_FLAG_BITS
      end

      # The default for 'select_alignment_proc'. Its order of preference
      # is align, shift2, shift1.
      # @param [AlignmentScore] score 
      # @return [CELL_FLAG_ALIGN, CELL_FLAG_SHIFT2, CELL_FLAG_SHIFT1] The symbol representing where to proceed.
      def default_select_alignment_proc(score)
        return CELL_FLAG_ALIGN if (score & CELL_FLAG_ALIGN) != 0
        return CELL_FLAG_SHIFT2 if (score & CELL_FLAG_SHIFT2) != 0
        CELL_FLAG_SHIFT1
      end

      # The default scoring proc. 
      # @param [Object] seq1_val
      # @param [Object] seq2_val
      # @return 1 if seq1_val and seq2_val are equal, and 0 if they are not.
      def default_score_proc(seq1_val, seq2_val)
        (seq1_val == seq2_val) ? 1 : 0
      end
    end # AlignmentMatrixHelpers
  end # NeedlemanWunsch
end # Align

