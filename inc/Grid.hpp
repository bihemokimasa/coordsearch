/* inc/Grid.hpp
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

// Simply put, a grid is a set of arrays (of doubles or floats)
// the arrays may be unequal dimensions


#ifndef GRID_HPP
#define GRID_HPP

#include <cstddef>
#include <numeric>
#include <algorithm>
#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <iomanip>

// #include <iostream>

#include "DynProg.hpp"
#include "Exception.hpp"

namespace DTDP
{

template <class T>
class Grid
{
public:
  Grid(){};
  Grid(const Grid &);
  Grid &operator=(const Grid &);
  Grid(Grid &&) noexcept;
  Grid &operator=(Grid &&) noexcept;
  virtual ~Grid() noexcept;

  explicit Grid(const Index *el_dim, Index size);
  explicit Grid(const std::vector<Index> &el_dim);
  explicit Grid(const T *elems, const std::vector<Index> &el_dim);
  /* Constructing a matrix or vector */
  //explicit Grid(Index dim1, Index dim2 = 0 /*for vector*/);

  void SetDims(const Index *el_dim, Index sz);
  void CopyElemDims(const Index *el_dim, Index sz);
  void SetDims(const std::vector<Index> &el_dim);
  void CopyElemDims(const std::vector<Index> &el_dim);

  Index GridSize(Index i_row) const noexcept;
  Index GridSize() const noexcept;
  Index NumElems() const noexcept { return aNumElems; }
  Index size() const noexcept { return GridSize(); }

  Index *ElemsDims() const noexcept { return &aElemsDims[0]; };

  // access
  T operator()(Index i_row, Index j_col) const { return GetElem(i_row, j_col); }
  T &operator()(Index i_row, Index j_col) { return GetElem(i_row, j_col); }

  virtual T GetElem(Index i_row, Index j_col) const;
  virtual T &GetElem(Index i_row, Index j_col);
  void GetElems(T *elems, Index size) const;

  virtual T const *Elems(Index i_row = 0, Index j_col = 0) const;
  virtual T *Elems(Index i_row = 0, Index j_col = 0);

  T const *begin() const { return &aElems[0]; };
  T *begin() { return &aElems[0]; };

  T const *end() const { return begin() + aGridSize; };
  T *end() { return begin() + aGridSize; };

  // dimension based comparison
  virtual bool VerifyDims(const Grid &g) const { return VerifyDims(g.aElemsDims, g.aNumElems, g.aGridSize); }
  virtual bool VerifyDims(Grid &g) const { return VerifyDims(g.aElemsDims, g.aNumElems, g.aGridSize); }
  virtual bool VerifyDims(Grid &&g) const { return VerifyDims(std::move(g)); };
  bool VerifyDims(const Index *elemsDim, Index numElems, Index gridSize) const;
  //bool operator==(const Grid &g) const;
  // bool operator==(Grid &&) const;

  Index StartInd(Index i_row) const;
  void CheckRange(Index i_row, Index j_col = 0) const;

  // IO

  virtual std::istream &Read(std::istream &is);
  virtual std::ostream &Write(std::ostream &os) const;

  friend std::istream &operator>>(std::istream &is, Grid<T> &g) { return g.Read(is); }
  friend std::ostream &operator<<(std::ostream &os, const Grid<T> &g) { return g.Write(os); };

  // T Max(Index i_row, Index beginColIndex, Index blockLength, Index stride) const;
  // T Min(Index i_row, Index beginColIndex, Index blockLength, Index stride) const;
  // std::pair<T, T> MinMax() const;
  // T Sum(Index i_row, Index beginColIndex, Index blockLength, Index stride) const;
  // Array1D<T> CumSum(Index i_row, Index beginColIndex, Index blockLength, Index stride) const;

protected:
  // void ReInitialize() noexcept;
  void Allocate();
  void DeAllocate() noexcept;

  void Copy(const T *elems, Index size);
  void Copy(const Grid &);
  // void Swap(Grid &) noexcept;
  // void Swap(Grid &&) noexcept;

private:
  T *aElems{nullptr};         //  pointer to elements of the grid, each a vector
  Index *aElemsDims{nullptr}; // dimensions of elements of the grid
  Index *aStartInd{nullptr};  // gets us to the first element of a grid element
  Index aNumElems{0};         // number of elements in a grid
  Index aGridSize{0};         // grid size, the total of aElemsDims
};

} // namespace DTDP

//
// Implementation of class Grid
//
namespace DTDP
{

template <class T>
Grid<T>::Grid(const Index *el_dim, Index sz) : Grid()
{
  SetDims(el_dim, sz);
  Allocate();
  CopyElemDims(el_dim, sz);
}

template <class T>
Grid<T>::Grid(const std::vector<Index> &el_dim) : Grid()
{
  SetDims(el_dim);
  Allocate();
  CopyElemDims(el_dim);
}

template <class T>
Grid<T>::Grid(const T *elems, const std::vector<Index> &el_dim)
{
  SetDims(el_dim);
  Allocate();
  CopyElemDims(el_dim);
  Copy(&elems[0], std::accumulate(el_dim.begin(), el_dim.end(), 0));
}

template <class T>
Grid<T>::Grid(const Grid &g) : Grid()
{
  // std::cout << g << std::endl;
  aNumElems = g.aNumElems;
  aGridSize = g.aGridSize;
  Allocate();
  Copy(g);
  // std::copy(&g.aElemsDims[0], &g.aElemsDims[0] + g.aNumElems, &aElemsDims[0]);
  // std::copy(&g.aStartInd[0], &g.aStartInd[0] + g.aNumElems, &aStartInd[0]);
  // std::copy(&g.aElems[0], &g.aElems[0] + g.aGridSize, &aElems[0]);
}

template <class T>
Grid<T> &Grid<T>::operator=(const Grid &g)
{
  if (this->VerifyDims(g.aElemsDims, g.aNumElems, g.aGridSize))
    Copy(g);
  else
  {
    DeAllocate();
    aNumElems = g.aNumElems;
    aGridSize = g.aGridSize;
    Allocate();
    Copy(g);
    // swap(g);
  }
  return *this;
}

template <class T>
Grid<T>::Grid(Grid &&g) noexcept : Grid() { *this = std::move(g); }

template <class T>
Grid<T> &Grid<T>::operator=(Grid &&g) noexcept
{
  if (this != &g)
  {
    this->DeAllocate();
    aElems = &g.aElems[0];
    aElemsDims = &g.aElemsDims[0];
    aStartInd = &g.aStartInd[0];
    aNumElems = g.aNumElems;
    aGridSize = g.aGridSize;

    g.aElems = nullptr;
    g.aElemsDims = nullptr;
    g.aStartInd = nullptr;
    g.aNumElems = 0;
    g.aGridSize = 0;

    // Allocate();
    // Swap(g);
    // g.ReInitialize();
  }
  return *this;
}

template <class T>
Grid<T>::~Grid() noexcept { DeAllocate(); }

// template <class T>
// void Grid<T>::ReInitialize() noexcept
// {
//   aElems = nullptr;
//   aElemsDims = nullptr;
//   aStartInd = nullptr;
//   aNumElems = 0;
//   aGridSize = 0;
// }

template <class T>
void Grid<T>::SetDims(const Index *el_dim, Index sz)
{
  aNumElems = sz;
  aGridSize = std::accumulate(&el_dim[0], &el_dim[0] + sz, 0);
}

template <class T>
void Grid<T>::CopyElemDims(const Index *el_dim, Index sz)
{
  std::copy(&el_dim[0], &el_dim[0] + sz, &aElemsDims[0]);
  aStartInd[0] = 0;
  for (Index i{1}; i < sz; i++)
    aStartInd[i] = aStartInd[i - 1] + el_dim[i - 1];
}

template <class T>
void Grid<T>::SetDims(const std::vector<Index> &el_dim)
{
  aNumElems = el_dim.size();
  aGridSize = std::accumulate(el_dim.begin(), el_dim.end(), 0);
}

template <class T>
void Grid<T>::CopyElemDims(const std::vector<Index> &el_dim)
{
  std::copy(el_dim.begin(), el_dim.end(), &aElemsDims[0]);
  aStartInd[0] = 0;
  for (Index i{1}; i < static_cast<Index>(el_dim.size()); i++)
    aStartInd[i] = aStartInd[i - 1] + el_dim[i - 1];
}

template <class T>
void Grid<T>::Allocate()
{
  aElems = new T[aGridSize]();
  aElemsDims = new Index[aNumElems]();
  aStartInd = new Index[aNumElems]();
}

template <class T>
void Grid<T>::DeAllocate() noexcept
{
  delete[] aElems;
  delete[] aElemsDims;
  delete[] aStartInd;
  aElems = nullptr;
  aElemsDims = nullptr;
  aStartInd = nullptr;
  aNumElems = 0;
  aGridSize = 0;
}

template <class T>
void Grid<T>::Copy(const Grid &g)
{
  std::copy(&g.aElemsDims[0], &g.aElemsDims[0] + g.aNumElems, &aElemsDims[0]);
  std::copy(&g.aStartInd[0], &g.aStartInd[0] + g.aNumElems, &aStartInd[0]);
  std::copy(&g.aElems[0], &g.aElems[0] + g.aGridSize, &aElems[0]);
}

template <class T>
void Grid<T>::Copy(const T *elems, Index size)
{
  std::copy(&elems[0], &elems[0] + size, &aElems[0]);
}

// template <class T>
// void Grid<T>::Swap(Grid &g) noexcept
// {
//   std::swap(aElems, g.aElems);
//   std::swap(aElemsDims, g.aElemsDims);
//   std::swap(aStartInd, g.aStartInd);
//   std::swap(aNumElems, g.aNumElems);
//   std::swap(aGridSize, g.aGridSize);
// }

// template <class T>
// void Grid<T>::Swap(Grid &&g) noexcept
// {
//   std::swap(aElems, g.aElems);
//   std::swap(aElemsDims, g.aElemsDims);
//   std::swap(aStartInd, g.aStartInd);
//   std::swap(aNumElems, g.aNumElems);
//   std::swap(aGridSize, g.aGridSize);
// }

template <class T>
Index Grid<T>::GridSize(Index i_row) const noexcept
{
  CheckRange(i_row);
  return aElemsDims[i_row];
}

template <class T>
Index Grid<T>::GridSize() const noexcept { return aGridSize; }

template <class T>
Index Grid<T>::StartInd(Index i_row) const
{
  CheckRange(i_row);
  return aStartInd[i_row];
}

template <class T>
T Grid<T>::GetElem(Index i_row, Index j_col) const
{
  CheckRange(i_row, j_col);
  return aElems[aStartInd[i_row] + j_col];
}

template <class T>
T &Grid<T>::GetElem(Index i_row, Index j_col)
{
  CheckRange(i_row, j_col);
  return aElems[aStartInd[i_row] + j_col];
}

template <class T>
void Grid<T>::GetElems(T *elems, Index size) const
{
  std::copy(&aElems[0], &aElems[0] + size, &elems[0]);
}

template <class T>
T const *Grid<T>::Elems(Index i_row, Index j_col) const
{
  CheckRange(i_row, j_col);
  return &aElems[aStartInd[i_row] + j_col];
}

template <class T>
T *Grid<T>::Elems(Index i_row, Index j_col)
{
  CheckRange(i_row, j_col);
  return &aElems[aStartInd[i_row] + j_col];
}

template <class T>
bool Grid<T>::VerifyDims(const Index *elemsDim, Index numElems, Index gridSize) const
{
  if (aGridSize != gridSize || aNumElems != numElems)
    return false;
  for (Index i{0}; i < aNumElems; i++)
    if (aElemsDims[i] != elemsDim[i])
      return false;
  return true;
}

template <class T>
void Grid<T>::CheckRange(Index i_row, Index j_col) const
{
  if (i_row >= aNumElems || j_col >= aElemsDims[i_row])
    Exception::error(__FILE__, __LINE__,
                     std::string{"RangeError at ("} +
                         std::to_string(i_row) +
                         std::string{", "} +
                         std::to_string(j_col) +
                         std::string{") where size = "} +
                         std::to_string(GridSize(i_row)),
                     static_cast<int>(ERROR_NUMBER::IGNORE));
}

template <class T>
std::istream &Grid<T>::Read(std::istream &is)
{
  std::string line, token;
  std::getline(is, line); // the first line has dimensions
  // extract the dimensions from the line, store in a vector
  std::istringstream ss(line);
  std::vector<Index> dims; //
  while (std::getline(ss, token, ' '))
    dims.push_back(std::stoi(token));

  // store the dimensions in an array of Index type
  Index numElems = dims.size();                             // number of elements
  Index elemDims[numElems];                                 // the array with the number of elements
  std::copy(dims.begin(), dims.end(), &elemDims[0]);        // copy
  Index gsz = std::accumulate(dims.begin(), dims.end(), 0); // calculate total grid size

  *this = Grid(dims);

  // Index i = 0;
  for (Index i{0}; i < gsz; i++)
    is >> aElems[i];
  // while ((is >> aElems[i]) && (i < gsz))
  //   i++;
  return is;
}

template <class T>
std::ostream &Grid<T>::Write(std::ostream &os) const
{
  for (Index i{0}; i < aNumElems; i++)
  {
    os << aElemsDims[i];
    if (i < aNumElems - 1)
      os << " ";
  }
  os << std::endl;

  for (Index i{0}; i < aNumElems; i++)
  {
    for (Index j{0}; j < aElemsDims[i]; j++)
    {
      os << std::setprecision(DOUBLE_PRECISION) << GetElem(i, j);
      if (j < aElemsDims[i] - 1)
        os << " ";
    }
    os << std::endl;
  }
  return os;
}
} // namespace DTDP

#endif