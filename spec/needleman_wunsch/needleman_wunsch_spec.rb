require_relative '../spec_helper'
require 'align/needleman_wunsch'
require 'benchmark'
require 'pp'

shared_examples_for "a needleman-wunsch alignment algorithm" do
  before :each do
    @seq1 = "GAATTCAGTTA"
    @seq2 = "GGATCGA"
  end

  context "for 'GATTCAGTTA' and 'GGATCGA'" do

    context "when aligning 'azzz' and 'zzz'" do
      before :each do
        @align1, @align2 = @align_proc.call('azzz', 'zzz', :skip_obj => '-')
      end
      
      it "should align 1 as 'azzz'" do
        @align1.join('').should == "azzz"
      end
      it "should align 2 as '-zzz'" do
        @align2.join('').should == "-zzz"
      end

    end
    context "by default" do
      before :each do
        @align1, @align2 = @align_proc.call(@seq1, @seq2, :skip_obj => '-')
      end

      it "should return arrays" do
        @align1.should be_a(Array)
        @align2.should be_a(Array)
      end
      
      it "should align 1 as 'GAATTCAGTTA'" do
        @align1.join('').should == "GAATTCAGTTA"
      end
      it "should align 2 as 'GGA-TC-G--A'" do
        @align2.join('').should == "GGA-TC-G--A"
      end

    end # by default
  end
end

describe Align::NeedlemanWunsch do
  before :each do
    @seq1 = "GAATTCAGTTA"
    @seq2 = "GGATCGA"
  end

## Commented out stuff that I use to help me assess performance.
#  it "should bench" do
#    Benchmark.bm do |x|
#      x.report("100     ") {100.times {Align::NeedlemanWunsch.new(@seq1, @seq2)}}
#      x.report("1000    ") {1000.times {Align::NeedlemanWunsch.new(@seq1, @seq2)}}
#      x.report("10000   ") {10000.times {Align::NeedlemanWunsch.new(@seq1, @seq2)}}
#    end
#  end
#
#  it "should bench" do
#    matrix = Align::NeedlemanWunsch.new(@seq1, @seq2)
#    Benchmark.bm do |x|
#      x.report("100     ") {100.times {matrix.traceback}}
#      x.report("1000    ") {1000.times {matrix.traceback}}
#      x.report("10000   ") {10000.times {matrix.traceback}}
#    end
#  end
#
#  it "should bench" do
#    Benchmark.bm do |x|
#      x.report("100     ") {100.times {Align::NeedlemanWunsch.new(@seq1, @seq2).traceback}}
#      x.report("1000    ") {1000.times {Align::NeedlemanWunsch.new(@seq1, @seq2).traceback}}
#      x.report("10000   ") {10000.times {Align::NeedlemanWunsch.new(@seq1, @seq2).traceback}}
#    end
#  end
  
  context "with two similar sequences" do
    before :each do
      @matrix = Align::NeedlemanWunsch.new(@seq1, @seq2)
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
    its(:rows) { should == @seq1.size + 1}
    its(:cols) { should == @seq2.size + 1}
    its(:highest_score) { should == 6}
    its(:highest_score_loc) { should == [11,7] }

    describe "#[]" do
      it "should raise an error when accessing an out of bound column" do
        lambda do
          @matrix[@matrix.rows+1,0]
        end.should raise_error(ArgumentError, "out of bounds (col: 13 >= 12 || row: 0 >= 8)")
      end
      it "should raise an error when accessing an out of bound row" do
        lambda do
          @matrix[0,@matrix.cols+1]
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

    it "should be possible to return the proper traceback" do
      @matrix.traceback_array.should == [
        [10,6, :align],
        [9,6, :delete],
        [8,6, :delete],
        [7,5, :align],
        [6,5, :delete],
        [5,4, :align],
        [4,3, :align],
        [3,3, :delete],
        [2,2, :align],
        [1,1, :align],
        [0,0, :align]
      ]
    end

    it "should be possible to construct the proper traceback via yield" do
      traceback = []
      @matrix.traceback do |i,j, last_move|
        traceback << [i,j, last_move]
      end
      traceback.should == [
        [10,6, :align],
        [9,6, :delete],
        [8,6, :delete],
        [7,5, :align],
        [6,5, :delete],
        [5,4, :align],
        [4,3, :align],
        [3,3, :delete],
        [2,2, :align],
        [1,1, :align],
        [0,0, :align]
      ]
    end
  end # "with two similar sequences"
end

describe Align::NeedlemanWunsch, :align do
  before :all do
    @align_proc = Align::NeedlemanWunsch.method(:align)
  end

  it_should_behave_like "a needleman-wunsch alignment algorithm" 
end

