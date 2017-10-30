#ifndef TROPICAL_CONVOLUTION_BUSSIECK_HXX
#define TROPICAL_CONVOLUTION_BUSSIECK_HXX

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

       template<class T, class C = std::vector<T>, class P = std::less<typename C::value_type> >
       struct iterable_priority_queue :
         std::priority_queue<T,C,P> {
           using std::priority_queue<T,C,P>::priority_queue;
           typename C::iterator begin() { return std::priority_queue<T, C, P>::c.begin(); }
           typename C::iterator end() { return std::priority_queue<T, C, P>::c.end(); }
         };

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
         auto compare = [](const auto x, const auto y){ return x.val < y.val; };
         std::sort(idx.begin(), idx.end(), compare);
         return idx; 
       }

       struct min_conv_index : public std::array<std::size_t,2> {};
       
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

       // note: this cover heuristic can only guarantee that if it returns that an element is covered that it is. 
       // However it may occur that an index is returned as not covered while it actually is
       struct cover : public std::vector<std::size_t> {
         cover(const std::size_t _n, const std::size_t _m) : 
           vector(_n + _m, 0),
           n(_n)
         {} 

         const std::size_t n; // number of left entries;

         void remove_cover(const std::size_t i, const std::size_t j)
         {
           //assert((*this)[i] > 0);
           //assert((*this)[n+j] > 0);
           //if( (*this)[i] > 0 )
           { (*this)[i]--; }
           //if( (*this)[n+j] > 0 )
           { (*this)[n+j]--; }
         }

         void add_cover(const std::size_t i, const std::size_t j)
         {
           //if( (*this)[i] == 0 && (*this)[n+j] == 0 ){
             (*this)[i]++; 
             (*this)[n+j]++;
           //}
         }

         bool covered(const std::size_t i, const std::size_t j) const
         {
           return !( (*this)[i] == 0 || (*this)[n+j] == 0 );
         }

         /*
         auto add_cover_I = [&](Index i,Index j){
           if( cover[i+1] == 0 && cover[n+j] == 0 ){
             INDEX inc = 1;
             while( i+inc < n ){
               INDEX idx = op(idxa_[i+inc],idxb_[j]);//idxa_[i+inc]+idxb_[j];
               idx = std::min(idx, INDEX(c_.size()-1));
               assert(idx<cp_.size());
               if( cp_[idx] == 0 ){
                 add_cover(i+inc,j);
                 break;
               }
               if( cover[i+inc] > 0 ){ break; }
               inc++;
             }
           }
         };

         auto addCoverJ = [&](Index i,Index j){
           if( cover[i] == 0 && cover[n+j+1] == 0 ){
                  Index inc = 1;
                  while( j+inc < m ){
                     Index idx = op(idxa_[i],idxb_[j+inc]);//idxa_[i]+idxb_[j+inc];
                     idx = std::min(idx, Index(c_.size()-1));
                     assert(idx<cp_.size());
                     if( cp_[idx] == 0 ){
                       add_cover(i,j+inc);
                        break;
                     }
                     if( cover[n+j+inc] > 0 ){ break; }
                     inc++;
                  }
               }
            }
       */
       };

     } // end namespace detail

      // additionally return the index coming from the first vector for the optimum convolution. 
      template<typename INPUT_ITERATOR_1, typename INPUT_ITERATOR_2, typename OUTPUT_ITERATOR_VAL, typename OUTPUT_ITERATOR_INDEX>
      void min_conv_Bussieck_et_al(
          INPUT_ITERATOR_1 a_begin, INPUT_ITERATOR_1 a_end, INPUT_ITERATOR_2 b_begin, INPUT_ITERATOR_2 b_end, 
          OUTPUT_ITERATOR_VAL result_begin, OUTPUT_ITERATOR_VAL result_end, OUTPUT_ITERATOR_INDEX result_index_a_begin)
      {
        using VALUE_TYPE = typename std::iterator_traits<INPUT_ITERATOR_1>::value_type;
        static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<INPUT_ITERATOR_2>::value_type>::value, "iterators must have same value type");
        static_assert(std::is_same<VALUE_TYPE, typename std::iterator_traits<OUTPUT_ITERATOR_VAL>::value_type>::value, "iterators must have same value type");
        VALUE_TYPE val = std::numeric_limits<VALUE_TYPE>::infinity();

       const auto a_size = std::distance(a_begin, a_end);
       const auto b_size = std::distance(b_begin, b_end);
       const auto result_size = std::distance(result_begin, result_end);
       assert(result_size <= a_size + b_size - 1);

       // store sorted indices according to vector a,b respectively
       auto idx_a = detail::sort_indices(a_begin, a_end);
       auto idx_b = detail::sort_indices(b_begin, b_end);

       std::fill(result_begin, result_end, std::numeric_limits<VALUE_TYPE>::infinity()); // the output. Note: infinity is not really the best marker, since the input can have infinities as well! (but then the algorithm can terminate early once it has reached such a value.

       std::size_t open = result_size;

       struct indices_value_pair {
         std::array<std::size_t,2> idx;
         VALUE_TYPE val;
       };
       auto compare = [](const indices_value_pair& x, const indices_value_pair& y) {
         return x.val > y.val; // prioirty_queue keeps the largest element, hence we have to invert sorting!
       };

       detail::iterable_priority_queue<indices_value_pair,std::vector<indices_value_pair>,decltype(compare) > queue(compare);
       queue.push(indices_value_pair({0, 0, idx_a[0].val+idx_b[0].val})); 

       detail::cover cover(a_size,b_size);
       cover.add_cover(0,0);


       while(open > 0) {
         assert(!queue.empty());

         const auto i = queue.top().idx[0];
         const auto j = queue.top().idx[1];
         const auto val = queue.top().val;
         queue.pop();
         cover.remove_cover(i,j);
         const auto underlying_i = idx_a[i].idx;
         const auto underlying_j = idx_b[j].idx;
         const auto underlying_sum = underlying_i + underlying_j;

         //std::cout << "val = " << val << ", sorted indices = (" << i << "," << j << "), true indices = (" << underlying_i << "," << underlying_j << ")\n";
         //for(auto it=queue.begin(); it!=queue.end(); ++it) {
         //  std::cout << "(" << (*it).idx[0] << "," << (*it).idx[1] << "); ";
         //}
         //std::cout << "\n";

         // new minimum found
         if(underlying_sum < result_size && result_begin[underlying_sum] == std::numeric_limits<VALUE_TYPE>::infinity()) {
           open--;
           result_begin[underlying_sum] = val;
           result_index_a_begin[underlying_sum] = underlying_i;
         }

         // update cover
         //if(idx_a[i[0]] + idx_b[i[1]] >= ) { continue; } // no need to update cover

         // FILL: cover elements
         //cover.add_cover(i[0] + 1, i[1]);
         //queue.push_back(indices({i[0]+1, i[1]}));
         //cover.add_cover(i[0], i[1] + 1);
         //queue.push_back(indices({i[0], i[1]+1}));

         // FILL1: only cover uncovered elements
         if(i + 1 < a_size) { 
           if(!cover.covered(i+1,j)) {
           //if(!detail::is_covered(i+1, j, queue.begin(), queue.end()))
             queue.push({i+1, j, idx_a[i+1].val + idx_b[j].val});
             cover.add_cover(i+1,j);
           }
         }

         if(j + 1 < b_size) {
           //if(!detail::is_covered(i, j+1, queue.begin(), queue.end()))
           if(!cover.covered(i,j+1)) {
             queue.push({i, j+1, idx_a[i].val + idx_b[j+1].val});
             cover.add_cover(i,j+1);
           }
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
