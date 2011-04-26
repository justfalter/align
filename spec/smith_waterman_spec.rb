require_relative 'spec_helper'
require 'align/smith_waterman'
require 'benchmark'

describe Align::SmithWaterman do
  before :each do
    @klass = Align::SmithWaterman
    @opts = {}
  end

  examples = [ 
    ['CAGCCUCGCUUAG', 'AAUGCCAUUGACGG', 'GCC-UCG', 'GCCAUUG'],
    ['AAUGCCAUUGACGG', 'CAGCCUCGCUUAG', 'GCCAUUG', 'GCC-UCG'],
    ['ACACACTA', 'AGCACACA', 'CACAC', 'CACAC']
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
  end

  context "with two similar sequences" do
    before :each do
      @seq1 = "CAGCCUCGCUUAG"
      @seq2 = "AAUGCCAUUGACGG"
      @nw = @klass.new(@seq1, @seq2, @opts)
    end

    subject {@nw}
    its(:rows) { should == @seq1.size + 1}
    its(:cols) { should == @seq2.size + 1}
    its(:max_score) { should == 6 }
    its(:max_row) {should == 8 }
    its(:max_col) {should == 10}

    it "should be possible to return the proper traceback" do
      @nw.traceback_array.should == 
      [
        [7, 9, :align],
        [6, 8, :align],
        [5, 7, :align],
        [5, 6, :insert],
        [4, 5, :align],
        [3, 4, :align],
        [2, 3, :align]
      ]
    end

    it "should be possible to construct the proper traceback via yield" do
      traceback = []
      @nw.traceback do |i,j, last_move|
        traceback << [i,j, last_move]
      end
      traceback.should == 
      [
        [7, 9, :align],
        [6, 8, :align],
        [5, 7, :align],
        [5, 6, :insert],
        [4, 5, :align],
        [3, 4, :align],
        [2, 3, :align]
      ]
    end
  end # "with two similar sequences"

end
