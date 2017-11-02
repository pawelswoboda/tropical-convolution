#ifndef TROPICAL_CONVOLUTION_BUSSIECK_HXX
#define TROPICAL_CONVOLUTION_BUSSIECK_HXX

#include "min_sum.hxx"
#include <limits>
#include <cmath>
#include <functional>
#include <queue>
#include <algorithm>
#include <vector>
#include <cassert>

// implement the tropical convolution algorithm described in
// M. Bussieck, H. Hassler, G. J. Woeginger, and U. T. Zimmermann:
// Fast algorithms for the maximum convolution problem.
// Oper. Res. Let.  , 15:1–5, 1994

namespace tropical_convolution {

  namespace detail {

       template<typename VALUE_TYPE>
       struct index_value_pair {
         VALUE_TYPE val;
         std::size_t idx;
       };

       template<typename ITERATOR>
       std::vector<index_value_pair<typename std::iterator_traits<ITERATOR>::value_type>> sort_indices(ITERATOR begin, ITERATOR end)
       {
         using VALUE_TYPE = typename std::iterator_traits<ITERATOR>::value_type;
         std::vector<index_value_pair<VALUE_TYPE>> idx(std::distance(begin,end));
         for(auto i=0; i<idx.size(); ++i) {
           idx[i].idx = i;
           idx[i].val = begin[i];
         }
         auto compare = [&](const auto& x, const auto& y){ return x.val < y.val; };
         std::sort(idx.begin(), idx.end(), compare);
         return idx; 
       }

       template<typename INDICES_ITERATOR>
       bool is_covered(std::size_t i, std::size_t j, INDICES_ITERATOR begin, INDICES_ITERATOR end)
       {
         for(auto it=begin; it!=end; ++it) {
           if(i >= (*it).idx[0] && j >= (*it).idx[1]) {
             return true;
           }
         }
         return false;
       }

       // fill missing values explicit min sum computation
       template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR_VAL, typename OUTPUT_ITERATOR_INDEX>
       void min_sum_missing(
           INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, 
           OUTPUT_ITERATOR_VAL result_begin, OUTPUT_ITERATOR_VAL result_end, OUTPUT_ITERATOR_INDEX result_index_a_begin)
       {
         using VALUE_TYPE = typename std::iterator_traits<INPUT_ITERATOR_1>::value_type;
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<INPUT_ITERATOR_2>::value_type>::value, "input iterators must have same value type");
         static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<OUTPUT_ITERATOR_VAL>::value_type>::value, "input and output iterators must have same value type");

         const auto a_size = std::distance(a_begin, a_end);
         const auto b_size = std::distance(b_begin, b_end);
         const auto result_size = std::distance(result_begin, result_end);

         for(auto i=0; i<result_size; ++i) {
           std::size_t a_idx, b_idx;
           if(result_begin[i] == std::numeric_limits<VALUE_TYPE>::infinity()) { // not assigned yet
             std::tie(result_begin[i], a_idx, b_idx) = arg_min_sum(a_begin, a_end, b_begin, b_end, i);
             result_index_a_begin[i] = a_idx;
           } 
         }
       }

     } // end namespace detail

      // return convolution values and additionally the index coming from the first vector for the optimum convolution. 
      template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR_VAL, typename OUTPUT_ITERATOR_INDEX>
      void min_conv_Bussieck_et_al(
          INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, 
          OUTPUT_ITERATOR_VAL result_begin, OUTPUT_ITERATOR_VAL result_end, OUTPUT_ITERATOR_INDEX result_index_a_begin)
      {
        using VALUE_TYPE = typename std::iterator_traits<INPUT_ITERATOR_1>::value_type;
        static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<INPUT_ITERATOR_2>::value_type>::value, "input iterators must have same value type");
        static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<OUTPUT_ITERATOR_VAL>::value_type>::value, "input and output iterators must have same value type");
        //VALUE_TYPE val = std::numeric_limits<VALUE_TYPE>::infinity();

       const auto a_size = std::distance(a_begin, a_end);
       const auto b_size = std::distance(b_begin, b_end);
       const auto result_size = std::distance(result_begin, result_end);
       assert(result_size <= a_size + b_size - 1);

       // store sorted indices of vector a and b respectively
       auto idx_a = detail::sort_indices(a_begin, a_end);
       auto idx_b = detail::sort_indices(b_begin, b_end);

       std::fill(result_begin, result_end, std::numeric_limits<VALUE_TYPE>::infinity()); // the output. Note: infinity is not really the best marker, since the input can have infinities as well! (but then the algorithm can terminate early once it has reached such a value.

       std::size_t open = result_size;

       struct indices_value_pair {
         std::size_t i;
         std::size_t j;
         VALUE_TYPE val;
       };
       auto compare = [](const indices_value_pair& x, const indices_value_pair& y) {
         return x.val > y.val; // std::priority_queue keeps the largest element on top, hence we have to invert sorting.
       };
       std::priority_queue<indices_value_pair,std::vector<indices_value_pair>,decltype(compare) > queue(compare);
       std::vector<unsigned char> cover_a(a_size,0);
       std::vector<unsigned char> cover_b(b_size,0);

       auto add_cover = [&](const std::size_t i, const std::size_t j) {
         assert(cover_a[i] <= 1 && cover_b[j] <= 1);
         if( cover_a[i] == 0 && cover_b[j] == 0) {
           cover_a[i] = 1;
           cover_b[j] = 1;
           queue.push({i,j, idx_a[i].val + idx_b[j].val});
         } 
       };

       auto remove_cover = [&](const std::size_t i, const std::size_t j) {
         assert(cover_a[i] <= 1);
         assert(cover_b[j] <= 1);
         if( cover_a[i] > 0 ) { cover_a[i] = 0; }
         if( cover_b[j] > 0 ) { cover_b[j] = 0; } 
       };

       auto add_cover_i = [&](const std::size_t i, const std::size_t j) {
         if(cover_a[i+1] == 0 && cover_b[j] == 0) {
           for(std::size_t inc=i+1; inc<a_size; ++inc) {
             const std::size_t underlying_sum = idx_a[inc].idx + idx_b[j].idx;
             if(underlying_sum < result_size && result_begin[underlying_sum] == std::numeric_limits<VALUE_TYPE>::infinity()) {
               add_cover(inc,j);
               return;
             }
             if(cover_a[inc] > 0) { return; }
           }
         }
       };

       auto add_cover_j = [&](const std::size_t i, const std::size_t j) {
         if(cover_a[i] == 0 && cover_b[j+1] == 0) {
           for(std::size_t inc=j+1; inc<b_size; ++inc) {
             const std::size_t underlying_sum = idx_a[i].idx + idx_b[inc].idx;
             if(underlying_sum < result_size && result_begin[underlying_sum] == std::numeric_limits<VALUE_TYPE>::infinity()) { 
               add_cover(i,inc);
               return;
             }
             if(cover_b[inc] > 0) { return; }
           }
         }
       };

       add_cover(0,0); 

       while(open > 0) {
         assert(!queue.empty());

         // check whether number of queue elements is larger than number of unassigned elements. If so, explicitly compute missing values
         if(queue.size() > 0.08*open) { // 0.08 seems to be a good value based on preliminary experiments
           detail::min_sum_missing(a_begin, a_end, b_begin, b_end, result_begin, result_end, result_index_a_begin);
           break;
         }

         const auto i = queue.top().i;
         assert(i < a_size);
         const auto j = queue.top().j;
         assert(j < b_size);
         const auto val = queue.top().val;
         queue.pop();
         remove_cover(i,j);
         if(i+j >= result_size) { continue; }
         const auto underlying_i = idx_a[i].idx;
         assert(underlying_i < a_size);
         const auto underlying_j = idx_b[j].idx;
         assert(underlying_j < b_size);
         const auto underlying_sum = underlying_i + underlying_j;

         // new minimum found
         if(underlying_sum < result_size && result_begin[underlying_sum] == std::numeric_limits<VALUE_TYPE>::infinity()) {
           open--;
           result_begin[underlying_sum] = val;
           result_index_a_begin[underlying_sum] = underlying_i;
         }

         if( i+1 < a_size && j+1 < b_size ) {
           if( cover_a[i+1] == 0 && cover_b[j] == 0 && cover_a[i] == 0 && cover_b[j+1] == 0) {
             add_cover(i+1,j);
             add_cover(i,j+1);
             //add_cover_i(i,j); // doing so and leaving out add_cover(i+1,j) is proposed by Bussieck et al, but in my experiments gives slower runtime.
           } else if( cover_a[i+1] == 0 && cover_b[j] == 0 ) {
             add_cover_i(i,j);
             //add_cover(i+1,j);
           } else if( cover_a[i] == 0 && cover_b[j+1] == 0 ){
             add_cover_j(i,j);
             //add_cover(i,j+1);
           }
         } else if( i+1 < a_size ) {
           add_cover_i(i,j);
         } else if( j+1 < b_size ) {
           add_cover_j(i,j);
         }
       }
     }

      // only return values of min convolution.
     template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR>
     void min_conv_Bussieck_et_al(INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, OUTPUT_ITERATOR result_begin, OUTPUT_ITERATOR result_end)
     {
       using VALUE_TYPE = typename std::iterator_traits<INPUT_ITERATOR_1>::value_type;
       // define output iterators that will do nothing (for not recording optimal indices).
       struct null_value {
         null_value& operator=(VALUE_TYPE) { return *this; }
       };
       struct null_iterator {
         null_value operator*() { return null_value(); }
         null_value operator[](std::size_t) { return null_value(); }
       } null_it;

       min_conv_Bussieck_et_al(a_begin, a_end, b_begin, b_end, result_begin, result_end, null_it);
     } 

} // end namespace tropical_convolution

#endif // TROPICAL_CONVOLUTION_BUSSIECK_HXX
