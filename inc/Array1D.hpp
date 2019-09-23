/* inc/Array1D.hpp
* 
* Copyright (C) 2019 Bihemo Kimasa
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef ARRAY1D_HPP
#define ARRAY1D_HPP

#include <cstddef>
#include <algorithm>
#include <utility>

#include "DynProg.hpp"
#include "Grid.hpp"


namespace DTDP
{
template <class T>
class Array1D : public Grid<T>
{
public:
  Array1D() = default;
  Array1D(const Array1D<T> &) = default;
  Array1D<T> &operator=(const Array1D<T> &) = default;
  Array1D(Array1D<T> &&) = default;
  Array1D<T> &operator=(Array1D<T> &&) = default;
  ~Array1D() = default;

  explicit Array1D(Index n);
  explicit Array1D(const T *elems, Index n);

  Index dim() const noexcept { return Grid<T>::GridSize(); }
  Index Size() const noexcept { return Grid<T>::GridSize(); }

  T operator()(Index i) const { return GetElem(i); }
  T &operator()(Index i) { return GetElem(i); }

  T operator[](Index i) const { return GetElem(i); }
  T &operator[](Index i) { return GetElem(i); }

  T GetElem(Index i) const { return Grid<T>::GetElem(0, i); }
  T &GetElem(Index i) { return Grid<T>::GetElem(0, i); }

  T const *Elems(Index pos = 0) const { return Grid<T>::Elems(0, pos); };
  T *Elems(Index pos = 0) { return Grid<T>::Elems(0, pos); }

  virtual std::istream &Read(std::istream &is) override;
  virtual std::ostream &Write(std::ostream &os) const override;

  Array1D<T> Slice(Index initpos, Index endpos) const;

  T Max() const { return *std::max_element(Grid<T>::begin(), Grid<T>::end()); }
  T Min() const { return *std::min_element(Grid<T>::begin(), Grid<T>::end()); }
  std::pair<T, T> MinMax() const;
  T Sum() const { return std::accumulate(Grid<T>::begin(), Grid<T>::end(), static_cast<T>(0)); }
  Array1D<T> CumSum() const;

private:
};
} // namespace DTDP

namespace DTDP
{
template <class T>
Array1D<T>::Array1D(Index n)
    : Grid<T>(std::vector<Index>{n})
{
  // std::vector<Index> el_dim = {n};
  // // Index n_elem{1};
  // // Index el_dim[n_elem] = {n};
  // Grid<T>::operator=(Grid<T>(el_dim));
}

template <class T>
Array1D<T>::Array1D(const T *elems, Index n)
    : Grid<T>(&elems[0], std::vector<Index>{n})
{
  // std::vector<Index> el_dim = {n};
  // Grid<T>::operator=(Grid<T>(elems, el_dim));
}

template <class T>
std::istream &Array1D<T>::Read(std::istream &is)
{
  Index N{0};
  is >> N;
  if (!(N == this->size()))
    *this = Array1D<T>(N);
  for (Index i = 0; i < N; i++)
    is >> this->GetElem(i);
  return is;
}

template <class T>
std::ostream &Array1D<T>::Write(std::ostream &os) const
{
  os << this->size() << std::endl;
  for (Index i = 0; i < this->size(); i++)
  {
    os << std::setprecision(DOUBLE_PRECISION) << this->GetElem(i);
    if (i < this->size() - 1)
      os << " ";
  }
  os << std::endl;
  return os;
}

template <class T>
Array1D<T> Array1D<T>::Slice(Index initpos, Index endpos) const
{
  Array1D<T> B(endpos - initpos + 1);
  for (Index i{0}; i < B.dim(); i++)
    B(i) = GetElem(initpos + i);
  return B;
}

// MATLAB like linspace function:
// r is a spacing factor, controls the width size: r>1 -> greater density near a
template <class T>
Array1D<T> LinearSpace(T a, T b, Index N, T r = static_cast<T>(1.0))
{
  Array1D<T> A(N);
  for (Index i{0}; i < N; i++)
    A(i) = a + (b - a) * std::pow(i / static_cast<T>(N - 1), r);
  return A;
}
// Chebyshev interpolation nodes

template <class T>
Array1D<T> ChebySpace(T a, T b, Index N)
{
  Array1D<T> A(N);
  A(0) = a;
  A(N - 1) = b;
  for (Index i{N - 2}, j{1}; i > 0; i--, j++)
    A(j) = (a + b) / 2 + ((b - a) / 2) * std::cos((2 * i + 1) * PI / (2 * (N - 1) + 2));
  return A;
}

template <class T>
std::pair<T, T> Array1D<T>::MinMax() const
{
  auto mmax = std::minmax_element(Grid<T>::begin(), Grid<T>::end());
  return std::make_pair(*mmax.first, *mmax.second);
}

template <class T>
Array1D<T> Array1D<T>::CumSum() const
{
  Array1D<T> A(Size());
  A(0) = GetElem(0);
  for (Index i{1}; i < Size(); i++)
    A(i) = A(i - 1) + GetElem(i);
  return A;
}
} // namespace DTDP

#endif