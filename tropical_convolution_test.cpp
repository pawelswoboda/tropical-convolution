#include "tropical_convolution.hxx"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <random>

inline void
test(const bool& pred)
{
  if (!pred)
    throw std::runtime_error("Test failed.");
}

template<typename VEC1, typename VEC2>
void test_naive_Bussieck(const VEC1 a, const VEC2 b)
{
  std::vector<double> result_naive(a.size() + b.size()-1);
  std::vector<double> result_Bussieck_et_al(a.size() + b.size()-1);
  tropical_convolution::min_conv_naive(a.begin(), a.end(), b.begin(), b.end(), result_naive.begin(), result_naive.end());
  tropical_convolution::min_conv_Bussieck_et_al(a.begin(), a.end(), b.begin(), b.end(), result_Bussieck_et_al.begin(), result_Bussieck_et_al.end());
  test(result_naive.size() == result_Bussieck_et_al.size());
  for(auto i=0; i<result_naive.size(); ++i) {
    test(result_naive[i] == result_Bussieck_et_al[i]);
  }

  std::vector<std::size_t> result_idx_naive(result_naive.size());
  std::vector<std::size_t> result_idx_Bussieck_et_al(result_Bussieck_et_al.size());
  tropical_convolution::min_conv_naive(a.begin(), a.end(), b.begin(), b.end(), result_naive.begin(), result_naive.end(), result_idx_naive.begin());
  tropical_convolution::min_conv_Bussieck_et_al(a.begin(), a.end(), b.begin(), b.end(), result_Bussieck_et_al.begin(), result_Bussieck_et_al.end(), result_idx_Bussieck_et_al.begin());
  test(result_naive.size() == result_Bussieck_et_al.size());
  for(auto i=0; i<result_naive.size(); ++i) {
    test(result_naive[i] == result_Bussieck_et_al[i]);
    test(a[result_idx_naive[i]] + b[i - result_idx_naive[i]] == result_naive[i]);
    test(a[result_idx_Bussieck_et_al[i]] + b[i - result_idx_Bussieck_et_al[i]] == result_Bussieck_et_al[i]);
  }

}

// test whether naive implementation and Bussieck et al algorithms return same results
int main()
{

  // artificial input
  {
   std::vector<double> a {0.1, 0.2, 0.05, 1};
   std::vector<double> b {0.1, 0.2, 0.05, 1};
   std::reverse(b.begin(), b.end());

   test_naive_Bussieck(a,b);
  }

  // random input
  {
    // initialize seed explicitly to make unit test reproducible
    // random numbers for size of underlying vectors
    std::mt19937 gen(1); //Standard mersenne_twister_engine seeded with 1
    std::uniform_int_distribution<> dis_int(2, 100);

    // random numbers for vector values
    std::uniform_real_distribution<> dis_real(1.0, 2.0);
    for(auto run=0; run<1000; ++run) {
      const auto a_size = dis_int(gen); 
      const auto b_size = dis_int(gen);

      std::vector<double> a(a_size);
      for(auto i=0; i<a.size(); ++i) {
        a[i] = dis_real(gen);
      }

      std::vector<double> b(b_size);
      for(auto i=0; i<b.size(); ++i) {
        b[i] = dis_real(gen);
      }

      test_naive_Bussieck(a,b);

    } 
  }
}

