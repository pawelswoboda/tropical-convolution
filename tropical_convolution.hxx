#ifndef TROPICAL_CONVOLUTION_HXX
#define TROPICAL_CONVOLUTION_HXX

#include <tuple>
#include <cassert>
#include "tropical_convolution_naive.hxx"
#include "tropical_convolution_bussieck.hxx"

// do zrobienia: possibly make own project out of min convolution, but then used vector class and associated memory allocator also needs to be a separate project

// perform min convolution. Choose automatically betweeen the naive version which is fast on small inputs and the method by Bussieck et al that has better asymptotic running time.

// do zrobienia: no return type, write directly into output iterators.

namespace tropical_convolution{

      // automatically choose between naive and efficient version of min convolution
      // when more elements than indicated by threshold need to be computed, use heuristic, otherwise use naive implementation.
      // a good value can be gleaned from the results of the benchmark program.
      static const std::size_t min_conv_threshold = 500; 

      template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR>
      void min_conv(INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, OUTPUT_ITERATOR result_begin, OUTPUT_ITERATOR result_end)
      {
         const auto result_size = std::distance(result_begin, result_end);
         if(result_size < min_conv_threshold) {
            return min_conv_naive(a_begin, a_end, b_begin, b_end, result_begin, result_end);
         } else {
            return min_conv_Bussieck_et_al(a_begin, a_end, b_begin, b_end, result_begin, result_end);
         }
      }

      // additionally return the index coming from the first vector for the optimum convolution. 
      template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR_VAL, typename OUTPUT_ITERATOR_INDEX>
      void min_conv(
          INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, 
          OUTPUT_ITERATOR_VAL result_begin, OUTPUT_ITERATOR_VAL result_end, OUTPUT_ITERATOR_INDEX result_index_a_begin)
      {
         const auto result_size = std::distance(result_begin, result_end);
         if(result_size < min_conv_threshold) {
            return min_conv_naive(a_begin, a_end, b_begin, b_end, result_begin, result_end, result_index_a_begin);
         } else {
            return min_conv_Bussieck_et_al(a_begin, a_end, b_begin, b_end, result_begin, result_end, result_index_a_begin);
         }
      } 

} // end namespace tropical_convolution

#endif // TROPICAL_CONVOLUTION_HXX
