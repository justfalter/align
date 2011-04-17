require_relative '../spec_helper'
require 'align/needleman_wunsch/alignment_matrix'
require 'benchmark'

require 'matrix'
require 'pp'

describe Align::NeedlemanWunsch::AlignmentMatrix do
  include Align::NeedlemanWunsch
  before :each do
    @seq1 = "GAATTCAGTTA"
    @seq2 = "GGATCGA"
  end

  describe :default_score_proc do
    it "should return 1 if two things are equal" do
      AlignmentMatrix.default_score_proc("a", "a").should == 1
    end

    it "should return 0 if two things are not equal" do
      AlignmentMatrix.default_score_proc("a", "b").should == 0
    end
  end # :default_score_proc

  describe :default_select_slignment_proc do
    it "should return :align all are equal" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(1,1,1)).should == :align
    end

    it "should return :align if shift2/align are the max" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(1,1,0)).should == :align
    end

    it "should return :shift2 if shift2/shift1 are the max" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(1,0,1)).should == :shift2
    end

    it "should return :shift2 if shift2 is the max" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(1,0,0)).should == :shift2
    end
    
    it "should return :align if align/shift1 are the max" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(0,1,1)).should == :align
    end

    it "should return :align if align is the max" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(0,1,0)).should == :align
    end
    
    it "should return :shift1 if shift1 is the max" do
      AlignmentMatrix.default_select_alignment_proc(AlignmentScore.new(0,0,1)).should == :shift1
    end
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

    it "should have the properly built matrix" do
      @matrix.to_a.should == 
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
    
    it "should be possible to return the proper traceback" do
      @matrix.traceback.should == [
        [11,7, :align],
        [10,6, :align],
        [9,6, :shift2],
        [8,6, :shift2],
        [7,5, :align],
        [6,5, :shift2],
        [5,4, :align],
        [4,3, :align],
        [3,3, :shift2],
        [2,2, :align],
        [1,1, :align]
      ]
    end

    it "should be possible to construct the proper traceback via yield" do
      traceback = []
      @matrix.traceback do |i,j, last_move|
        traceback << [i,j, last_move]
      end
      traceback.should == [
        [11,7, :align],
        [10,6, :align],
        [9,6, :shift2],
        [8,6, :shift2],
        [7,5, :align],
        [6,5, :shift2],
        [5,4, :align],
        [4,3, :align],
        [3,3, :shift2],
        [2,2, :align],
        [1,1, :align]
      ]
    end

  end # "with two similar sequences"
end
