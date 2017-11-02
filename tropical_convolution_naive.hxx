#ifndef TROPICAL_CONVOLUTION_NAIVE_HXX
#define TROPICAL_CONVOLUTION_NAIVE_HXX

#include <limits>
#include <functional>
#include <algorithm>
#include <tuple>
#include <cassert>

// implement naive min convolution with O(n^2) runtime.
// to do: use SIMD

namespace tropical_convolution {

      template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR>
      void min_conv_naive(INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, OUTPUT_ITERATOR result_begin, OUTPUT_ITERATOR result_end)
      {
         using VALUE_TYPE = typename std::iterator_traits<INPUT_ITERATOR_1>::value_type;
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<INPUT_ITERATOR_2>::value_type>::value, "input iterators must have same value type");
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<OUTPUT_ITERATOR>::value_type>::value, "input and output iterators must have same value type");

         const auto a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const auto b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);
         const auto result_size = std::distance(result_begin, result_end);
         assert(result_size <= a_size + b_size - 1);
         std::fill(result_begin, result_end, std::numeric_limits<VALUE_TYPE>::infinity());

         for(auto i=0; i<std::min(result_size, a_size); ++i) {
            for(auto j=0; j<std::min(b_size, result_size-1); ++j) {
               result_begin[i+j] = std::min(result_begin[i+j], a_begin[i] + b_begin[j]); 
            }
         }
      }

      // additionally return the index coming from the first vector for the optimum convolution. 
      template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR_VAL, typename OUTPUT_ITERATOR_INDEX>
      void min_conv_naive(
          INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, 
          OUTPUT_ITERATOR_VAL result_begin, OUTPUT_ITERATOR_VAL result_end, OUTPUT_ITERATOR_INDEX result_index_a_begin)
      {
         using VALUE_TYPE = typename std::iterator_traits<INPUT_ITERATOR_1>::value_type;
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<INPUT_ITERATOR_2>::value_type>::value, "input iterators must have same value type");
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<OUTPUT_ITERATOR_VAL>::value_type>::value, "input and output iterators must have same value type");

         const auto a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const auto b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);
         const auto result_size = std::distance(result_begin, result_end);
         assert(result_size <= a_size + b_size - 1);
         std::fill(result_begin, result_end, std::numeric_limits<VALUE_TYPE>::infinity());

         for(auto i=0; i<std::min(result_size, a_size); ++i) {
            for(auto j=0; j<std::min(b_size, result_size - i); ++j) {
               const auto cur_val = a_begin[i] + b_begin[j];
               if(cur_val <= result_begin[i+j]) {
                  result_begin[i+j] = a_begin[i] + b_begin[j];
                  result_index_a_begin[i+j] = i;
               }
            }
         } 
      }

} // end namespace tropical_convolution

#endif // TROPICAL_CONVOLUTION_NAIVE_HXX
