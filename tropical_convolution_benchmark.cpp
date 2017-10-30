#include "tropical_convolution.hxx"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

void benchmark(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b) 
{
  const auto runs = a.size();
  const auto size = a[0].size();

  std::vector<double> result(2*size-1);

  auto begin_naive = std::chrono::high_resolution_clock::now();
  for(auto i=0; i<runs; ++i) {
    tropical_convolution::min_conv_naive(a[i].begin(), a[i].end(), b[i].begin(), b[i].end(), result.begin(), result.end());
  }
  auto end_naive = std::chrono::high_resolution_clock::now();

  auto begin_Bussieck_et_al = std::chrono::high_resolution_clock::now();
  for(auto i=0; i<runs; ++i) {
    tropical_convolution::min_conv_Bussieck_et_al(a[i].begin(), a[i].end(), b[i].begin(), b[i].end(), result.begin(), result.end());
  }
  auto end_Bussieck_et_al = std::chrono::high_resolution_clock::now();

  std::cout << "Input size = " << size 
    << "; naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_naive-begin_naive).count() 
    << "(ms); Bussieck et al: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_Bussieck_et_al-begin_Bussieck_et_al).count() << "(ms)" << std::endl;
}

// compare runtimes of tropical convolution implementations
int main()
{
  // random numbers for vector values
  std::mt19937 gen(1); //Standard mersenne_twister_engine seeded with 1
  std::uniform_real_distribution<> dis_real(1.0, 2.0);

  // compare runtimes for various vector sizes
  auto max_size = 10000;
  for(auto vec_size=10; vec_size<max_size; vec_size*=2) {

    std::vector<std::vector<double>> a(max_size/vec_size);
    for(auto& vec : a) {
      vec.resize(vec_size);
      for(auto i=0; i<vec.size(); ++i) {
        vec[i] = dis_real(gen);
      }
    }
    std::vector<std::vector<double>> b(max_size/vec_size);
    for(auto& vec : b) {
      vec.resize(vec_size);
      for(auto i=0; i<vec.size(); ++i) {
        vec[i] = dis_real(gen);
      }
    }

    benchmark(a,b);
  }
}

