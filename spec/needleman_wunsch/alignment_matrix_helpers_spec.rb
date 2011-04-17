require_relative '../spec_helper'
require 'align/needleman_wunsch/alignment_matrix_helpers'
require 'benchmark'


describe Align::NeedlemanWunsch::AlignmentMatrixHelpers do
  include Align::NeedlemanWunsch::Constants
  include Align::NeedlemanWunsch::AlignmentMatrixHelpers

  describe :default_score_proc do
    it "should return 1 if two things are equal" do
      default_score_proc("a", "a").should == 1
    end

    it "should return 0 if two things are not equal" do
      default_score_proc("a", "b").should == 0
    end
  end # :default_score_proc

  describe :default_select_slignment_proc do
    it "should return CELL_FLAG_ALIGN all are equal" do
      default_select_alignment_proc(create_cell(1,1,1)).should == CELL_FLAG_ALIGN
    end

    it "should return CELL_FLAG_ALIGN if shift2/align are the max" do
      default_select_alignment_proc(create_cell(1,0,1)).should == CELL_FLAG_ALIGN
    end

    it "should return CELL_FLAG_SHIFT2 if shift2/shift1 are the max" do
      default_select_alignment_proc(create_cell(0,1,1)).should == CELL_FLAG_SHIFT2
    end

    it "should return CELL_FLAG_SHIFT2 if shift2 is the max" do
      default_select_alignment_proc(create_cell(0,0,1)).should == CELL_FLAG_SHIFT2
    end
    
    it "should return CELL_FLAG_ALIGN if align/shift1 are the max" do
      default_select_alignment_proc(create_cell(1,1,0)).should == CELL_FLAG_ALIGN
    end

    it "should return CELL_FLAG_ALIGN if align is the max" do
      default_select_alignment_proc(create_cell(1,0,0)).should == CELL_FLAG_ALIGN
    end
    
    it "should return CELL_FLAG_SHIFT1 if shift1 is the max" do
      default_select_alignment_proc(create_cell(0,1,0)).should == CELL_FLAG_SHIFT1
    end
  end
end
