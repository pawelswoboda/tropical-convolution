#ifndef MIN_SUM_HXX
#define MIN_SUM_HXX

#include <tuple>
#include <cassert>

namespace tropical_convolution {

      template<typename ITERATOR_1, typename ITERATOR_2>
      typename std::iterator_traits<ITERATOR_1>::value_type min_sum(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const std::size_t sum)
      {
         using VALUE_TYPE = typename std::iterator_traits<ITERATOR_1>::value_type;
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<ITERATOR_2>::value_type>::value, "iterators must have same value type");
         VALUE_TYPE val = std::numeric_limits<VALUE_TYPE>::infinity();

         const auto a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const auto b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);

         assert(sum <= a_size-1 + b_size-1);

         for(auto i=std::max(0, int(sum)-int(b_size)+1); i<std::min(sum+1, std::size_t(a_size)); ++i) {
            val = std::min(val, a_begin[i] + b_begin[sum-i]);
         }
         return val;
      }

      // find out the two indices (left,right) whose sum is minimal
      template<typename ITERATOR_1, typename ITERATOR_2>
      std::tuple<typename std::iterator_traits<ITERATOR_1>::value_type,std::size_t,std::size_t> 
      arg_min_sum(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const std::size_t sum)
      {
         using VALUE_TYPE = typename std::iterator_traits<ITERATOR_1>::value_type;
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<ITERATOR_2>::value_type>::value, "iterators must have same value type");
         VALUE_TYPE val = std::numeric_limits<VALUE_TYPE>::infinity();

         std::size_t a = std::numeric_limits<std::size_t>::max();
         std::size_t b = std::numeric_limits<std::size_t>::max();

         const auto a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const auto b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);

         assert(sum <= a_size-1 + b_size-1);


         for(auto i=std::max(0, int(sum)-int(b_size)+1); i<std::min(sum+1, std::size_t(a_size)); ++i) {
           const auto cur_val = a_begin[i] + b_begin[sum-i];
           if(cur_val <= val) {
             a = i;
             b = sum - i;
             val = cur_val;
           }
         }
         assert(a < a_size && b < b_size);
         return std::make_tuple(val,a,b);;
      }

}

#endif
