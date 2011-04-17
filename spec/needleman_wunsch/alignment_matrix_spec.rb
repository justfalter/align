require_relative '../spec_helper'
require 'align/needleman_wunsch/alignment_matrix'
require 'benchmark'

require 'matrix'
require 'pp'

describe Align::NeedlemanWunsch::AlignmentMatrix do
  include Align::NeedlemanWunsch
  include Align::NeedlemanWunsch::Constants
  before :each do
    @seq1 = "GAATTCAGTTA"
    @seq2 = "GGATCGA"
  end

## Commented out stuff that I use to help me assess performance.
#  it "should bench" do
#    Benchmark.bm do |x|
#      x.report("100     ") {100.times {Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2)}}
#      x.report("1000    ") {1000.times {Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2)}}
#      x.report("10000   ") {10000.times {Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2)}}
#    end
#  end
#
#  it "should bench" do
#    matrix = Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2)
#    Benchmark.bm do |x|
#      x.report("100     ") {100.times {matrix.traceback}}
#      x.report("1000    ") {1000.times {matrix.traceback}}
#      x.report("10000   ") {10000.times {matrix.traceback}}
#    end
#  end
#
#  it "should bench" do
#    Benchmark.bm do |x|
#      x.report("100     ") {100.times {Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2).traceback}}
#      x.report("1000    ") {1000.times {Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2).traceback}}
#      x.report("10000   ") {10000.times {Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2).traceback}}
#    end
#  end
  
  context "with two similar sequences" do
    before :each do
      @matrix = Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2)
    end

    # Actual matrix
    #     #0 1 2 3 4 5 6 7 8 9 A B 
    #     #- G A A T T C A G T T A
    # 0 - [0,0,0,0,0,0,0,0,0,0,0,0]
    # 1 G [0,1,1,1,1,1,1,1,1,1,1,1]
    # 2 G [0,1,1,1,1,1,1,1,2,2,2,2]
    # 3 A [0,1,2,2,2,2,2,2,2,2,2,3]
    # 4 T [0,1,2,2,3,3,3,3,3,3,3,3]
    # 5 C [0,1,2,2,3,3,4,4,4,4,4,4]
    # 6 G [0,1,2,2,3,3,4,4,5,5,5,5]
    # 7 A [0,1,2,3,3,3,4,5,5,5,5,6]

    # Actual matrix
    #     #0 1 2 3 4 5 6 7 8 9 A B 
    #     #- G A A T T C A G T T A
    # 0 - [0,0,0,0,0,0,0,0,0,0,0,0]
    # 1 G [0,X,1,1,1,1,1,1,1,1,1,1]
    # 2 G [0,1,X,1,1,1,1,1,2,2,2,2]
    # 3 A [0,1,2,X,X,2,2,2,2,2,2,3]
    # 4 T [0,1,2,2,3,X,3,3,3,3,3,3]
    # 5 C [0,1,2,2,3,3,X,X,4,4,4,4]
    # 6 G [0,1,2,2,3,3,4,4,X,X,X,5]
    # 7 A [0,1,2,3,3,3,4,5,5,5,5,X]

    subject {@matrix}
    its(:cols) { should == @seq1.size + 1}
    its(:rows) { should == @seq2.size + 1}
    its(:score) { should == 6}
    its(:gap_penalty) {should == 0}

    describe "#[]" do
      it "should raise an error when accessing an out of bound column" do
        lambda do
          @matrix[@matrix.cols+1,0]
        end.should raise_error(ArgumentError, "out of bounds (col: 13 >= 12 || row: 0 >= 8)")
      end
      it "should raise an error when accessing an out of bound row" do
        lambda do
          @matrix[0,@matrix.rows+1]
        end.should raise_error(ArgumentError, "out of bounds (col: 0 >= 12 || row: 9 >= 8)")
      end
    end

    it "should have a properly built score matrix" do
      @matrix.to_score_matrix.should == 
        [
        [0,0,0,0,0,0,0,0],
        [0,1,1,1,1,1,1,1],
        [0,1,1,2,2,2,2,2],
        [0,1,1,2,2,2,2,3],
        [0,1,1,2,3,3,3,3],
        [0,1,1,2,3,3,3,3],
        [0,1,1,2,3,4,4,4],
        [0,1,1,2,3,4,4,5],
        [0,1,2,2,3,4,5,5],
        [0,1,2,2,3,4,5,5],
        [0,1,2,2,3,4,5,5],
        [0,1,2,3,3,4,5,6],
      ]
    end

    it "should have a properly built traceback matrix" do
      @matrix.to_traceback_matrix.should == 
        [
        [7,7,7,7,7,7,7,7],
        [7,1,3,2,2,2,3,2],
        [7,4,7,1,2,2,2,3],
        [7,4,7,5,7,7,7,1],
        [7,4,7,4,1,2,2,6],
        [7,4,7,4,5,7,7,7],
        [7,4,7,4,4,1,2,2],
        [7,4,7,5,4,4,7,1],
        [7,5,1,6,4,4,1,6],
        [7,4,4,7,5,4,4,7],
        [7,4,4,7,5,4,4,7],
        [7,4,4,1,6,4,4,1],
      ]
    end
    
    it "should be possible to return the proper traceback" do
      @matrix.traceback.should == [
        [11,7, CELL_FLAG_ALIGN],
        [10,6, CELL_FLAG_ALIGN],
        [9,6, CELL_FLAG_SHIFT2],
        [8,6, CELL_FLAG_SHIFT2],
        [7,5, CELL_FLAG_ALIGN],
        [6,5, CELL_FLAG_SHIFT2],
        [5,4, CELL_FLAG_ALIGN],
        [4,3, CELL_FLAG_ALIGN],
        [3,3, CELL_FLAG_SHIFT2],
        [2,2, CELL_FLAG_ALIGN],
        [1,1, CELL_FLAG_ALIGN]
      ]
    end

    it "should be possible to construct the proper traceback via yield" do
      traceback = []
      @matrix.traceback do |i,j, last_move|
        traceback << [i,j, last_move]
      end
      traceback.should == [
        [11,7, CELL_FLAG_ALIGN],
        [10,6, CELL_FLAG_ALIGN],
        [9,6, CELL_FLAG_SHIFT2],
        [8,6, CELL_FLAG_SHIFT2],
        [7,5, CELL_FLAG_ALIGN],
        [6,5, CELL_FLAG_SHIFT2],
        [5,4, CELL_FLAG_ALIGN],
        [4,3, CELL_FLAG_ALIGN],
        [3,3, CELL_FLAG_SHIFT2],
        [2,2, CELL_FLAG_ALIGN],
        [1,1, CELL_FLAG_ALIGN]
      ]
    end
  end # "with two similar sequences"

  describe "#traceback" do
    it "should raise an exception if the select_alignment_proc returns something other than CELL_FLAG_ALIGN, CELL_FLAG_SHIFT1, or CELL_FLAG_SHIFT2" do
      select_proc = lambda do |score|
        :foo
      end

      lambda do
        matrix = Align::NeedlemanWunsch::AlignmentMatrix.new(@seq1, @seq2, 
                                          :select_alignment_proc => select_proc)
        matrix.traceback
      end.should raise_error(StandardError, "invalid return from select_alignment_proc: :foo")
    end

  end
end
