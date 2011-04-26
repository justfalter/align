require_relative 'spec_helper'
require 'align/pairwise_algorithm'
require 'align/basic_scoring'

describe Align::PairwiseAlgorithm do
  before :each do
    @seq1 = "GAATTCAGTTA"
    @seq2 = "GGATCGA"
    @scoring = Align::BasicScoring.new(1,-1,-1)
  end
  subject {Align::PairwiseAlgorithm.new(@seq1, @seq2, @scoring)}
  its(:seq1) {should == @seq1}
  its(:seq2) {should == @seq2}
  its(:scoring) {should == @scoring}
  describe :align do
    it "should raise a NotImplementedError" do
      lambda {subject.align}.should raise_error(NotImplementedError)
    end
  end

  it "#max2 should find the max amongst two numbers" do
    [1,2].permutation.each do |a,b|
      subject.max2(a,b).should == 2
    end
  end

  it "#max3 should find the max amongst three numbers" do
    [1,2,3].permutation.each do |a,b,c|
      subject.max3(a,b,c).should == 3
    end
  end

  it "#max4 should find the max amongst four numbers" do
    [1,2,3,4].permutation.each do |a,b,c,d|
      subject.max4(a,b,c,d).should == 4
    end
  end
end
