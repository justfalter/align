require_relative '../spec_helper'
require 'align/needleman_wunsch'
require 'pp'

describe Align::NeedlemanWunsch::AlignmentMatrix do
  before :each do
    @seq1 = "GAATTCAGTTA"
    @seq2 = "GGATCGA"
  end

  context "for 'GATTCAGTTA' and 'GGATCGA'" do

    describe :align do
      context "when aligning 'azzz' and 'zzz'" do
        before :each do
          @align1, @align2 = Align::NeedlemanWunsch.align('azzz', 'zzz', :skip_obj => '-')
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
          @align1, @align2 = Align::NeedlemanWunsch.align(@seq1, @seq2, :skip_obj => '-')
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

      context "by changing the select_alignment_proc to favor shift2, shift1, align" do
        before :each do
          select_align_proc = lambda do |score|
            if score.shift2?
              :shift2
            elsif score.shift1?
              :shift1
            else
              :align
            end
          end

          @align1, @align2 = Align::NeedlemanWunsch.align(@seq1, @seq2, :skip_obj => '-', 
                                                  :select_alignment_proc => select_align_proc)
        end # before


        it "should align 1 as 'G-AATTCAGTTA'" do
          @align1.join('').should == "G-AATTCAGTTA"
        end
        it "should align 2 as 'GGA-T-C-G--A'" do
          @align2.join('').should == "GGA-T-C-G--A"
        end
      end

      context "by changing the select_alignment_proc to favor shift2, align, shift1" do
        before :each do
          select_align_proc = lambda do |score|
            if score.shift2?
              :shift2
            elsif score.align?
              :align
            else
              :shift1
            end
          end

          @align1, @align2 = Align::NeedlemanWunsch.align(@seq1, @seq2, :skip_obj => '-', 
                                                  :select_alignment_proc => select_align_proc)
        end # before


        it "should align 1 as '-GAATTCAGTTA'" do
          @align1.join('').should == "-GAATTCAGTTA"
        end
        it "should align 2 as 'GGA-T-C-G--A'" do
          @align2.join('').should == "GGA-T-C-G--A"
        end
      end

      context "by changing the select_alignment_proc to foo" do
        before :each do
          select_align_proc = lambda do |score|
            if score.shift2?
              :shift2
            elsif score.align?
              :align
            else
              :shift1
            end
          end
          score_proc = lambda do |a,b|
            a == b ? 2 : -1
          end

          @align1, @align2 = Align::NeedlemanWunsch.align(@seq1, @seq2, :skip_obj => '-', 
                                                  :gap_penalty => -2, :score_proc => score_proc,
                                                  :select_alignment_proc => select_align_proc)
        end # before


        it "should align 1 as 'GAATTCAGTTA'" do
          @align1.join('').should == "GAATTCAGTTA"
        end
        it "should align 2 as 'GGAT-C-G--A'" do
          @align2.join('').should == "GGAT-C-G--A"
        end
      end
    end

  end
end
