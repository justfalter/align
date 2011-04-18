require 'align/alignment_matrix'

module Align
  # References:
  # [http://www.avatar.se/molbioinfo2001/dynprog/dynamic.html]
  module NeedlemanWunsch
    # Aligns two sequences together.
    # @param (see AlignmentMatrix.initialize)
    # @param [Hash] opts
    # @option opts [Object] :skip_obj (nil) The object that will be inserted
    #  when it is found that we must skip.
    def self.align(seq1, seq2, opts = {})
      skip_obj = opts[:skip_obj] || nil
      matrix = Align::AlignmentMatrix.new(seq1, seq2, opts)

      alignment_1 = []
      alignment_2 = []

      last_i = last_j = 0
      traceback = matrix.traceback.reverse
      traceback.each do |i,j,last_move|
        if last_i != i && last_j != j
          alignment_1 << seq1[i-1]
          alignment_2 << seq2[j-1]
        elsif last_i != i
          alignment_1 << seq1[i-1]
          alignment_2 << skip_obj
        elsif last_j != j
          alignment_1 << skip_obj
          alignment_2 << seq2[j-1]
        end
        last_i = i
        last_j = j
      end

      [alignment_1, alignment_2]
    end
  end
end
