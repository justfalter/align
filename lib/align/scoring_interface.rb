module Align
  # Basic Scoring interface
  class ScoringInterface
    def score_align(a,b)
      raise NotImplementedError.new("#{self.class}#score_align")
    end

    def score_insert(a)
      raise NotImplementedError.new("#{self.class}#score_align")
    end

    def score_delete(a)
      raise NotImplementedError.new("#{self.class}#score_align")
    end
  end
end
