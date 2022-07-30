#ifndef _ARRAYMDCPU_H
#define _ARRAYMDCPU_H

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

template <typename T>
struct Array1D
{
  unsigned n1;
  unsigned size;
  T *dptr;

  /**
   * 对()进行重载
   * 取第 i1 个元素
   *
   * */
  inline T &operator()(unsigned i1)
  {
    return dptr[i1];
  }

  Array1D() = default;

  Array1D(const Array1D &p)
  {
    n1 = p.n1;
    size = 0;
    dptr = p.dptr;
  }

  Array1D(int in1)
  {
    n1 = in1;
    size = n1;
    // 创建一个大小为in1的数组
    dptr = new T[size];
  }

  ~Array1D()
  {
    if (size && dptr)
      delete[] dptr;
  }
  /**
   * 获取占用空间大小
   * 元素个数 * 每个元素大小
   * 
   * */
  // unsigned getSize() { return size * sizeof(T); } // NB: in bytes
};

template <typename T>
struct Array2D
{
  unsigned y_len, x_len, b1;
  unsigned size;
  T *dptr;

  /**
   * 
   * */
  inline T &operator()(unsigned y, unsigned x)
  {
    return dptr[x + (x_len * y)];
  }

  // 正方形
  Array2D()
  {
    y_len = x_len = 0, b1 = 0;
    size = 0;
    dptr = NULL;
  }

  Array2D(const Array2D &p)
  {
    y_len = p.y_len;
    x_len = p.x_len;
    size = 0;
    dptr = p.dptr;
  }

  Array2D(int y, int x)
  {
    y_len = y;
    x_len = x;
    size = y_len * x_len;
    dptr = new T[size];
  }

  ~Array2D()
  {
    if (size && dptr)
      delete[] dptr;
  }

  // unsigned getSize() { return size * sizeof(T); } // NB: in bytes
};

// template<typename T>
// struct Array3D
// {
//   unsigned n1, n2, n3;
//   unsigned size;
//   T* dptr;

//   inline T& operator()(unsigned i1, unsigned i2, unsigned i3)
//   {
//     return dptr[i3 + i2 * f2 + i1 * f1];
//   }

//   Array3D() = default;

//   Array3D(const Array3D& p)
//   {
//     n1 = p.n1;
//     n2 = p.n2;
//     n3 = p.n3;
//     size = 0;
//     dptr = p.dptr;
//     f2 = n3;
//     f1 = f2 * n2;
//   }

//   Array3D(unsigned in1, unsigned in2, unsigned in3)
//   {
//     n1 = in1;
//     n2 = in2;
//     n3 = in3;
//     size = n1 * n2 * n3;
//     f2 = n3;
//     f1 = f2 * n2;
//     dptr = new T[size];
//   }

//   ~Array3D()
//   {
//     if (size && dptr)
//       delete[] dptr;
//   }

//   unsigned getSize() { return size * sizeof(T); }  // NB: in bytes

// private:
//   unsigned f2, f1, b1, b2;
// };

// template<typename T>
// struct Array4D
// {
//   unsigned n1, n2, n3, n4;
//   unsigned size;
//   T* dptr;

//   inline T& operator()(unsigned i1, unsigned i2, unsigned i3, unsigned i4)
//   {
//     return dptr[i4 + i3 * f3 + i2 * f2 + i1 * f1];
//   }

//   Array4D() = default;

//   Array4D(const Array4D& p)
//   {
//     n1 = p.n1;
//     n2 = p.n2;
//     n3 = p.n3;
//     n4 = p.n4;
//     size = 0;
//     dptr = p.dptr;
//     f3 = n4;
//     f2 = f3 * n3;
//     f1 = f2 * n2;
//   }

//   Array4D(unsigned in1, unsigned in2, unsigned in3, unsigned in4)
//   {
//     n1 = in1;
//     n2 = in2;
//     n3 = in3;
//     n4 = in4;
//     size = n1 * n2 * n3 * n4;
//     f3 = n4;
//     f2 = f3 * n3;
//     f1 = f2 * n2;
//     dptr = new T[size];
//   }

//   ~Array4D()
//   {
//     if (size && dptr)
//       delete[] dptr;
//   }

// private:
//   unsigned f3, f2, f1, b1, b2, b3;
// };

#endif
