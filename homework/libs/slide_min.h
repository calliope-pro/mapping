#ifndef LIBS_SLIDE_MIN_H_
#define LIBS_SLIDE_MIN_H_

#include <deque>

namespace libs
{
  using Index = int;
  using Value = int;

  class SlideMin
  {
  private:
    std::deque<std::pair<Index, Value>> deq_;
    const int section_length_;

  public:
    SlideMin(Value init_value, const int section_length)
        : deq_(), section_length_(section_length)
    {
      deq_.emplace_back(std::pair<Index, Value>(0, init_value));
    }

    std::pair<Index, Value> add(Value new_value)
    {
      Index new_index = deq_.back().first + 1;
      while (!deq_.empty() && deq_.back().second >= new_value)
      {
        deq_.pop_back();
      }
      if (!deq_.empty() && deq_.front().first <= new_index - section_length_)
      {
        deq_.pop_front();
      }
      deq_.emplace_back(std::pair<Index, Value>(new_index, new_value));
      return min();
    }
    Value last_value()
    {
      return deq_.back().second;
    }
    std::pair<Index, Value> min()
    {
      return deq_.front();
    }
  };
}

#endif // LIBS_SLIDE_MIN_H_
