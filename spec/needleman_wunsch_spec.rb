require_relative 'spec_helper'
require 'align/needleman_wunsch'

shared_examples_for "a needleman-wunsch alignment algorithm" do
  examples = [ 
    ['bcefg', 'abcdef', '-bc-efg', 'abcdef-'],
    ['GAATTCAGTTA', 'GGATCGA', 'GAATTCAGTTA', 'GGA-TC-G--A']
  ]

  examples.each do |seq1, seq2, exp1, exp2|
    context "when aligning '#{seq1}' and '#{seq2}'" do
      before :each do
        @opts[:skip_obj] = "-"
        @align1, @align2 = @klass.align(seq1, seq2, @opts)
      end
      
      it "should align 1 as '#{exp1}'" do
        @align1.join('').should == exp1
      end
      it "should align 2 as '#{exp2}'" do
        @align2.join('').should == exp2
      end
    end

    context "when aligning '#{seq2}' and '#{seq1}'" do
      before :each do
        @opts[:skip_obj] = "-"
        @align1, @align2 = @klass.align(seq2, seq1, @opts)
      end
      
      it "should align 1 as '#{exp2}'" do
        @align1.join('').should == exp2
      end
      it "should align 2 as '#{exp1}'" do
        @align2.join('').should == exp1
      end
    end

  end

  context "with two similar sequences" do
    before :each do
      @seq1 = "GAATTCAGTTA"
      @seq2 = "GGATCGA"
      @nw = @klass.new(@seq1, @seq2, @opts)
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

    subject {@nw}
    its(:rows) { should == @seq1.size + 1}
    its(:cols) { should == @seq2.size + 1}

    it "should use the score for the lower, right" do
      @nw.score.should == @score
    end

#    it "should have a properly built score matrix" do
#      @nw.matrix.should == 
#        [
#        [0,0,0,0,0,0,0,0],
#        [0,1,1,1,1,1,1,1],
#        [0,1,1,2,2,2,2,2],
#        [0,1,1,2,2,2,2,3],
#        [0,1,1,2,3,3,3,3],
#        [0,1,1,2,3,3,3,3],
#        [0,1,1,2,3,4,4,4],
#        [0,1,1,2,3,4,4,5],
#        [0,1,2,2,3,4,5,5],
#        [0,1,2,2,3,4,5,5],
#        [0,1,2,2,3,4,5,5],
#        [0,1,2,3,3,4,5,6],
#      ]
#    end

    it "should be possible to return the proper traceback" do
      @nw.traceback_array.should == [
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
      @nw.traceback do |i,j, last_move|
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

describe Align::NeedlemanWunsch do
  before :each do
    @klass = Align::NeedlemanWunsch
    @opts = {}
    @score = 6
  end
  context "by default" do
    it_should_behave_like "a needleman-wunsch alignment algorithm" 
  end

  context "using SCORING_ALT1" do
    before :each do
      @klass = Align::NeedlemanWunsch
      @opts[:scoring] = Align::NeedlemanWunsch::SCORING_ALT1
      @score = 1
    end

    it_should_behave_like "a needleman-wunsch alignment algorithm" 
  end
end

