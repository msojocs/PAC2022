// #include "arrayMD.h"
#include <complex>
#include <chrono>
#include <cmath>
#include <iostream>
#include <immintrin.h>

using namespace std;
using namespace chrono;

using namespace std::chrono;

#define nstart 0
#define nend 3

using DataType = double;
using DataType_VEC = __m512d;

#define ComplexType std::complex<DataType>

// ArrayMD definitions
// #define ARRAY2D Array2D<ComplexType>
// #define ARRAY1D Array1D<ComplexType>
// #define ARRAY1D_int Array1D<int>
// #define ARRAY1D_DataType Array1D<DataType>

// Function Definitions

void
noflagOCC_solver(size_t number_bands,
                 size_t ngpown,
                 size_t ncouls,
                 int* inv_igp_index,
                 int* indinv,
                 DataType* wx_array,
                 DataType* wtilde_array_a,
                 DataType* wtilde_array_b,
                 DataType* aqsmtemp_a,
                 DataType* aqsmtemp_b,
                 DataType* aqsntemp_a,
                 DataType* aqsntemp_b,
                 DataType* I_eps_array_a,
                 DataType* I_eps_array_b,
                 DataType* vcoul,
                 ComplexType* achtemp);

inline void ComplexType_print(ComplexType &src)
{
  printf("(%f,%f) \n",src.real(),src.imag());
}

// 共轭复数 如果我们将复数z表示为(real，img)，则其共轭为(real，-img)
inline ComplexType ComplexType_conj(ComplexType& src)
{
  return (std::conj(src));
}
