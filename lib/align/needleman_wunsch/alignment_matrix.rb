module Align
  module NeedlemanWunsch
    # Stores a score in the AlignmentMatrix. This makes our code an awful
    # lot easier to use when it comes to caluculating the 
    class AlignmentScore
      attr_reader :shift2, :align, :shift1, :value

      # @param [Fixnum] shift2
      # @param [Fixnum] align
      # @param [Fixnum] shift1
      def initialize(shift2=0, align=0, shift1=0)
        @shift2 = shift2
        @align = align
        @shift1   = shift1
        @value = [@shift2, @align, @shift1].max
      end

      # @return [Boolean] True if the value came from the shift2
      def shift2?; @shift2 == @value; end

      # @return [Boolean] True if the value came from the align
      def align?; @align == @value; end

      # @return [Boolean] True if the value came from above
      def shift1?;   @shift1   == @value; end

      # Basic test for equality.
      def ==(other)
        @value == other
      end
    end

    # Obviously, this would be loads faster if this were done in C. But, I'm
    # leaving it as is, for now, for portability reasons.
    class AlignmentMatrix
      attr_reader :rows, :cols, :gap_penalty, :score_proc, :select_alignment_proc

      # @param [#[], #size] seq1 The first sequence
      # @param [#[], #size] seq2 The second sequence
      # @param [Hash] opts Options
      # @option opts [Fixnum] :gap_penalty (0) The gap penalty
      # @option opts [Proc] :score_proc A proc that compares two sequence items,
      #  and returns a score. Defaults to :default_score_proc
      # @option opts [Proc] :select_alignment_proc  A proc that allows for the
      #  controlling of which way to go about aligning if there are multiple
      #  paths available. Defaults to AlignmentMatrix.default_select_alignment_proc
      def initialize(seq1, seq2, opts = {})
        @cols = seq1.size + 1
        @rows = seq2.size + 1
        @gap_penalty = opts[:gap_penalty] || 0
        @score_proc = opts[:score_proc] || AlignmentMatrix.method(:default_score_proc)
        @select_alignment_proc = opts[:select_alignment_proc] || 
          AlignmentMatrix.method(:default_select_alignment_proc)

        @matrix = Array.new(@cols) do 
          Array.new(@rows)
        end

        0.upto(@cols-1) {|i| @matrix[i][0] = AlignmentScore.new(0)}
        0.upto(@rows-1) {|j| @matrix[0][j] = AlignmentScore.new(0)}

        fill(seq1, seq2)
      end

      # Returns the matrix data. 
      def to_a
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
      def fill(seq1, seq2)
        1.upto(@cols-1) do |i|
          1.upto(@rows-1) do |j|
            score_align = @matrix[i-1][j-1].value + @score_proc.call(seq1[i-1], seq2[j-1])
            score_shift2 = @matrix[i-1][j].value + @gap_penalty
            score_shift1 = @matrix[i][j-1].value + @gap_penalty
            @matrix[i][j] = AlignmentScore.new(score_shift2, score_align, score_shift1)
          end
        end
      end # fill

      # @return [Fixnum] The global alignment score.
      def score
        self[@cols - 1, @rows - 1]
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

        last_move = :align

        while (i > 0 && j > 0)
          trace << [i,j, last_move]
          yield(i,j, last_move) if block_given?

          last_move = @select_alignment_proc.call(@matrix[i][j])
          case last_move
          when :shift2; i -= 1
          when :shift1;   j -= 1
          when :align; i -= 1; j -= 1
          else
            raise "invalid return from select_alignment_proc: #{last_move.inspect}"
          end
        end # while

        while i > 0
          trace << [i,j,last_move]
          yield(i,j,last_move) if block_given?
          last_move = :shift2
          i-=1
        end

        while j > 0
          trace << [i,j,last_move]
          yield(i,j,last_move) if block_given?
          last_move = :shift1
          j-=1
        end

        trace
      end # traceback

      # The default for 'select_alignment_proc'. Its order of preference
      # is align, shift2, shift1.
      # @param [AlignmentScore] score 
      # @return [:align, :shift2, :shift1] The symbol representing where to proceed.
      def AlignmentMatrix.default_select_alignment_proc(score)
        return :align if score.align?
        return :shift2 if score.shift2?
        :shift1
      end

      # The default scoring proc. 
      # @param [Object] seq1_val
      # @param [Object] seq2_val
      # @return 1 if seq1_val and seq2_val are equal, and 0 if they are not.
      def AlignmentMatrix.default_score_proc(seq1_val, seq2_val)
        (seq1_val == seq2_val) ? 1 : 0
      end
    end # AlignmentMatrix
  end
end
