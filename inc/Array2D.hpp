/* inc/Array2D.hpp
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

#ifndef ARRAY2D_HPP
#define ARRAY2D_HPP

#include "Grid.hpp"
#include "Array1D.hpp"

namespace DTDP
{
template <class T>
class Array2D : public Grid<T>
{
public:
  using Index = int;
  Array2D() = default;
  Array2D(const Array2D<T> &) = default;
  Array2D<T> &operator=(const Array2D<T> &) = default;
  Array2D(Array2D<T> &&) = default;
  Array2D<T> &operator=(Array2D<T> &&) = default;
  ~Array2D() = default;

  explicit Array2D(Index rows, Index cols);
  explicit Array2D(const T *elems, Index rows, Index cols);

  Index RowSize() const noexcept { return aNumRows; }
  Index ColSize() const noexcept { return aNumCols; }
  Index Size() const noexcept { return Grid<T>::GridSize(); }
  // Index size() const noexcept { return Grid<T>::GridSize(); }

  // Index num_rows() const noexcept { return aNumRows; }
  // Index num_cols() const noexcept { return aNumCols; }

  // std::tuple<Index, Index> size() const noexcept { return std::make_tuple(aNumRows, aNumCols); }

  virtual T GetElem(Index i_row, Index j_col) const override { return Grid<T>::GetElem(0, LinearInd(i_row, j_col)); } // { return this->aElems[LinearInd(i_row, j_col)]; } // { return this->aElems[LinearInd(i_row, j_col)]; }
  virtual T &GetElem(Index i_row, Index j_col) override { return Grid<T>::GetElem(0, LinearInd(i_row, j_col)); }

  virtual T const *Elems(Index i_row, Index j_col) const override { return Grid<T>::Elems(0, LinearInd(i_row, j_col)); };
  virtual T *Elems(Index i_row, Index j_col) override {return Grid<T>::Elems(0, LinearInd(i_row, j_col)); }

  Index LinearInd(Index i_row, Index j_col) const noexcept
  {
    return (i_row * aNumCols + j_col);
  }
  Index RowIndex(Index k) const noexcept { return (k / aNumCols); }
  Index ColIndex(Index k) const noexcept { return k % aNumCols; }

  bool DimsAreEqual(Index rows, Index cols) const { return aNumRows == rows && aNumCols == cols; }

  // linear index
  Array1D<T> Row(Index i) const;
  Array1D<T> Col(Index j) const;

  // Array1D<T> row(Index i) const { return Row(i); };
  // Array1D<T> col(Index j) const { return Col(j); };
  // IO
  virtual std::istream &Read(std::istream &is) override;
  virtual std::ostream &Write(std::ostream &os) const override;

  T Max(Index beginRow, Index beginCol, Index blockLength, Index stride) const;
  T Min(Index beginRow, Index beginCol, Index blockLength, Index stride) const;
  T Sum(Index beginRow, Index beginCol, Index blockLength, Index stride) const;
  Array1D<T> CumSum(Index beginRow, Index beginCol, Index blockLength, Index stride) const;

  T Max() const { return *std::max_element(Grid<T>::Begin(), Grid<T>::End()); }
  T RowMax(Index i_row) const { return Max(i_row, 0, aNumCols, 1); }
  T ColMax(Index j_col) const { return Max(0, j_col, aNumRows, aNumCols); }
  Array1D<T> RowMax() const;
  Array1D<T> ColMax() const;

  T Min() const { return *std::min_element(Grid<T>::Begin(), Grid<T>::End()); }
  T RowMin(Index i_row) const { return Min(i_row, 0, aNumCols, 1); };
  T ColMin(Index j_col) const { return Min(0, j_col, aNumRows, aNumCols); };
  Array1D<T> RowMin() const;
  Array1D<T> ColMin() const;

  std::pair<T, T> MinMax() const;

  T Sum() const { return std::accumulate(Grid<T>::Begin(), Grid<T>::End(), static_cast<T>(0)); }
  T RowSum(Index i_row) const { return Sum(i_row, 0, aNumCols, 1); };
  T ColSum(Index j_col) const { return Sum(0, j_col, aNumRows, aNumCols); };
  Array1D<T> RowSum() const;
  Array1D<T> ColSum() const;

  Array1D<T> RowCumSum(Index i_row) const { return CumSum(i_row, 0, aNumCols, 1); };
  Array1D<T> ColCumSum(Index j_col) const { return CumSum(0, j_col, aNumRows, aNumCols); }

  Array2D<T> RowCumSum() const;
  Array2D<T> ColCumSum() const;

private:
  Index aNumRows{0}, aNumCols{0};
};
} // namespace DTDP

//
// Implementation of class Array2D
//
namespace DTDP
{
template <class T>
Array2D<T>::Array2D(Index rows, Index cols)
    : Grid<T>(std::vector<Index>{rows * cols}),
      aNumRows{rows}, aNumCols{cols}
{
  // std::vector<Index> el_dim = {rows * cols};
  // Index n_elem{1};
  // Index el_dim[n_elem] = {rows * cols};
  // Grid<T>::operator=(Grid<T>(el_dim));
}

template <class T>
Array2D<T>::Array2D(const T *elems, Index rows, Index cols)
    : Grid<T>(&elems[0], std::vector<Index>{rows * cols}),
      aNumRows{rows}, aNumCols{cols}
{
  // std::vector<Index> el_dim = {rows * cols};
  // Grid<T>::operator=(Grid<T>(elems, el_dim));
}

template <class T>
Array1D<T> Array2D<T>::Row(Index i) const
{
  Array1D<T> A(aNumCols);
  for (Index j = 0; j < aNumCols; j++)
    A(j) = GetElem(i, j);
  return A;
}

template <class T>
Array1D<T> Array2D<T>::Col(Index j) const
{
  Array1D<T> A(aNumRows);
  for (Index i = 0; i < aNumRows; i++)
    A(i) = GetElem(i, j);
  return A;
}

// IO
template <class T>
std::istream &Array2D<T>::Read(std::istream &is)
{
  Index M{0}, N{0};
  is >> M >> N;
  if (!(M == aNumRows && N == aNumCols))
    *this = Array2D<T>(M, N);
  for (Index i = 0; i < M; i++)
    for (Index j = 0; j < N; j++)
      is >> GetElem(i, j);
  return is;
}

template <class T>
std::ostream &Array2D<T>::Write(std::ostream &os) const
{
  os << aNumRows << " " << aNumCols << std::endl;
  for (Index i = 0; i < aNumRows; i++)
  {
    for (Index j = 0; j < aNumCols; j++)
    {
      os << std::setprecision(DOUBLE_PRECISION) << GetElem(i, j);
      if (j < ((aNumCols - 1)))
        os << " ";
    }
    os << std::endl;
  }
  return os;
}

template <class T>
T Array2D<T>::Max(Index beginRow, Index beginCol, Index blockLength, Index stride) const
{
  Index k{LinearInd(beginRow, beginCol)};
  T val = GetElem(RowIndex(k), ColIndex(k));
  for (Index j{0}; j < blockLength; j++, k += stride)
    val = std::max(val, GetElem(RowIndex(k), ColIndex(k)));
  return val;
}

template <class T>
T Array2D<T>::Min(Index beginRow, Index beginCol, Index blockLength, Index stride) const
{
  Index k{LinearInd(beginRow, beginCol)};
  T val = GetElem(RowIndex(k), ColIndex(k));
  for (Index j{0}; j < blockLength; j++, k += stride)
    val = std::min(val, GetElem(RowIndex(k), ColIndex(k)));
  return val;
}

template <class T>
T Array2D<T>::Sum(Index beginRow, Index beginCol, Index blockLength, Index stride) const
{
  Index k{LinearInd(beginRow, beginCol)};
  T val = GetElem(RowIndex(k), ColIndex(k));
  for (Index j{0}; j < blockLength; j++, k += stride)
    val += GetElem(RowIndex(k), ColIndex(k));
  return val;
}

template <class T>
Array1D<T> Array2D<T>::RowMax() const
{
  Array1D<T> A(RowSize());
  for (Index i{0}; i < RowSize(); i++)
    A(i) = RowMax(i);
  return A;
}

template <class T>
Array1D<T> Array2D<T>::ColMax() const
{
  Array1D<T> A(ColSize());
  for (Index j{0}; j < ColSize(); j++)
    A(j) = ColMax(j);
  return A;
}

template <class T>
Array1D<T> Array2D<T>::RowMin() const
{
  Array1D<T> A(RowSize());
  for (Index i{0}; i < RowSize(); i++)
    A(i) = RowMin(i);
  return A;
}

template <class T>
Array1D<T> Array2D<T>::ColMin() const
{
  Array1D<T> A(ColSize());
  for (Index j{0}; j < ColSize(); j++)
    A(j) = ColMin(j);
  return A;
}

template <class T>
std::pair<T, T> Array2D<T>::MinMax() const
{
  auto mmax = std::minmax_element(Grid<T>::Begin(), Grid<T>::End());
  return std::make_pair(*mmax.first, *mmax.second);
}

template <class T>
Array1D<T> Array2D<T>::RowSum() const
{
  Array1D<T> A(RowSize());
  for (Index i{0}; i < RowSize(); i++)
    A(i) = RowSum(i);
  return A;
}

template <class T>
Array1D<T> Array2D<T>::ColSum() const
{
  Array1D<T> A(ColSize());
  for (Index j{0}; j < ColSize(); j++)
    A(j) = ColSum(j);
  return A;
}

template <class T>
Array1D<T> Array2D<T>::CumSum(Index beginRow, Index beginCol, Index blockLength, Index stride) const
{
  Index k{LinearInd(beginRow, beginCol)};
  Array1D<T> A(blockLength);
  A(0) = GetElem(RowIndex(k), ColIndex(k));
  k += stride;
  for (Index j{1}; j < blockLength; j++, k += stride)
    A(j) = A(j - 1) + GetElem(RowIndex(k), ColIndex(k));
  return A;
}

template <class T>
Array2D<T> Array2D<T>::RowCumSum() const
{
  Array2D<T> A(RowSize(), ColSize());
  for (Index i{0}; i < RowSize(); i++)
  {
    A(i, 0) = GetElem(i, 0);
    for (Index j{1}; j < ColSize(); j++)
      A(i, j) = A(i, j - 1) + GetElem(i, j);
  }
  return A;
}
template <class T>
Array2D<T> Array2D<T>::ColCumSum() const
{
  Array2D<T> A(RowSize(), ColSize());
  for (Index i{0}; i < ColSize(); i++)
  {
    A(0, i) = GetElem(0, i);
    for (Index j{1}; j < RowSize(); j++)
      A(j, i) = A(j - 1, i) + GetElem(j, i);
  }
  return A;
}

} // namespace DTDP

#endif