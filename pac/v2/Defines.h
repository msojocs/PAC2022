#include "arrayMD.h"
#include <complex>
#include <chrono>
#include <cmath>
#include <iostream>
using namespace std;
using namespace chrono;

using namespace std::chrono;

#define nstart 0
#define nend 3

using DataType = double;

#define ComplexType std::complex<DataType>

// ArrayMD definitions
#define ARRAY2D Array2D<ComplexType>
#define ARRAY1D Array1D<ComplexType>
#define ARRAY1D_int Array1D<int>
// #define ARRAY1D_DataType Array1D<DataType>

// Function Definitions

void
noflagOCC_solver(size_t number_bands,
                 size_t ngpown,
                 size_t ncouls,
                 ARRAY1D_int& inv_igp_index,
                 ARRAY1D_int& indinv,
                 DataType* wx_array,
                 ARRAY2D& wtilde_array,
                 ARRAY2D& aqsmtemp,
                 ARRAY2D& aqsntemp,
                 ARRAY2D& I_eps_array,
                 DataType* vcoul,
                 ARRAY1D& achtemp);

inline void ComplexType_print(ComplexType &src)
{
  printf("(%f,%f) \n",src.real(),src.imag());
}

// 共轭复数 如果我们将复数z表示为(real，img)，则其共轭为(real，-img)
inline ComplexType ComplexType_conj(ComplexType& src)
{
  return (std::conj(src));
}
