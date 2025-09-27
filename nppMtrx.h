#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <type_traits>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <functional>
#include <numeric>
#include "nppMath.h"
#include "nppCmplx.h"
#include "nppVctr.h"

namespace npp
{
  template <class _T>
  class Mtrx
  {
    static_assert(
      std::is_arithmetic_v<_T> || std::is_same_v<_T, npp::Cmplx>,
      "The npp::Mtrx template class requires an arithmetic type like int, float, double, npp::Cmplx, etc."
    );

  public:
    using ValueType = _T;
    using SizeType = typename std::vector<_T>::size_type;
    using Iterator = typename std::vector<_T>::iterator;
    using ConstIterator = typename std::vector<_T>::const_iterator;

    class RowProxy
    {
    private:
      _T *m_rowData;
      SizeType m_cols;

    public:
      RowProxy(_T *, const SizeType);
      _T &operator[](const SizeType) noexcept;
      const _T &operator[](const SizeType) const noexcept;
    };

    class ConstRowProxy
    {
    private:
      const _T *m_rowData;
      SizeType m_cols;

    public:
      ConstRowProxy(const _T *, const SizeType);
      const _T &operator[](const SizeType) const noexcept;
    };

  private:
    std::vector<_T> m_data;
    SizeType m_rows, m_cols;

  public:
    Mtrx();
    explicit Mtrx(const SizeType, const SizeType);
    explicit Mtrx(const SizeType, const SizeType, const _T &);
    Mtrx(std::initializer_list<std::initializer_list<_T>>);
    Mtrx(const Mtrx &);
    Mtrx(Mtrx &&) noexcept;
    ~Mtrx() noexcept;

  private:
    SizeType index(const SizeType, const SizeType) const noexcept;

    void checkRowBounds(const SizeType) const;
    void checkColBounds(const SizeType) const;
    void checkBounds(const SizeType, const SizeType) const;

    void checkSizeCompatibility(const SizeType, const SizeType, const SizeType, const SizeType) const;
    void checkSizeCompatibility(const Mtrx &) const;
    
    void checkMultCompatibility(const SizeType, const SizeType, const SizeType, const SizeType) const;
    void checkMultCompatibility(const Mtrx &) const;

    void checkSquareMatrix(const SizeType, const SizeType) const;
    void checkSquareMatrix(const Mtrx &) const;
    
  public:
    Mtrx &operator=(const Mtrx &);
    Mtrx &operator=(Mtrx &&) noexcept;
    void fill(const _T &);

    _T &at(const SizeType, const SizeType);
    const _T &at(const SizeType, const SizeType) const;

    RowProxy operator[](const SizeType) noexcept;
    ConstRowProxy operator[](const SizeType) const noexcept;

    SizeType rows() const noexcept;
    SizeType cols() const noexcept;
    std::vector<SizeType> dims() const noexcept;
    SizeType elements() const noexcept;

    bool isEmpty() const noexcept;
    bool isSquare() const noexcept;

    void resize(const SizeType, const SizeType);
    void resize(const SizeType, const SizeType, const _T &);

    Iterator begin() noexcept;
    ConstIterator begin() const noexcept;
    ConstIterator cbegin() const noexcept;

    Iterator end() noexcept;
    ConstIterator end() const noexcept;
    ConstIterator cend() const noexcept;

    Mtrx &operator+=(const Mtrx &);
    Mtrx &operator+=(const _T &);

    Mtrx &operator-=(const Mtrx &);
    Mtrx &operator-=(const _T &);

    Mtrx &operator*=(const Mtrx &);
    Mtrx &operator*=(const _T &);

    Mtrx &operator/=(const Mtrx &);
    Mtrx &operator/=(const _T &);

    const Mtrx &operator+() const;
    Mtrx operator+(const Mtrx &) const;
    Mtrx operator+(const _T &) const;

    Mtrx operator-() const;
    Mtrx operator-(const Mtrx &) const;
    Mtrx operator-(const _T &) const;

    Mtrx operator*(const Mtrx &) const;
    Mtrx operator*(const _T &) const;

    Mtrx operator/(const Mtrx &) const;
    Mtrx operator/(const _T &) const;

    bool operator==(const Mtrx &) const;
    bool operator==(const _T &) const;

    bool operator!=(const Mtrx &) const;
    bool operator!=(const _T &) const;

    Vctr<_T> getRow(const SizeType) const;
    Vctr<_T> getCol(const SizeType) const;

    Mtrx T() const;
    Mtrx t() const;
    Mtrx H() const;
    Mtrx h() const;

    Mtrx dot(const Mtrx &) const;

    Mtrx minor(const SizeType, const SizeType) const;

    ValueType det() const;
    static ValueType det(const Mtrx &);

    Mtrx adj() const;

    Mtrx inv() const;
    Mtrx lInv(const Mtrx &) const;
    Mtrx rInv(const Mtrx &) const;

    ValueType sum() const;
    ValueType min() const;
    ValueType max() const;

    template <class UnaryFunc>
    Mtrx &apply(UnaryFunc);

    template <class UnaryFunc>
    Mtrx map(UnaryFunc) const;

    void print() const;
    void print(const std::string &) const;

    static Mtrx zeros(const SizeType, const SizeType);
    static Mtrx ones(const SizeType, const SizeType);
    static Mtrx diag(const SizeType, const ValueType &);
    static Mtrx diag(const Vctr<_T> &);
  };
}

template <class _T>
npp::Mtrx<_T>::Mtrx() : m_data(), m_rows(0), m_cols(0) {}

template <class _T>
npp::Mtrx<_T>::Mtrx(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols) : m_data(rows * cols), m_rows(rows), m_cols(cols) {}

template <class _T>
npp::Mtrx<_T>::Mtrx(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols, const _T &value) : m_data(rows * cols, value), m_rows(rows), m_cols(cols) {}

template <class _T>
npp::Mtrx<_T>::Mtrx(std::initializer_list<std::initializer_list<_T>> rows) : m_rows(rows.size()), m_cols(rows.size() > 0 ? rows.begin()->size() : 0)
{
  for (const auto &row : rows)
    if (row.size() != m_cols)
      throw std::invalid_argument("Column numbers are not consistent");

  m_data.reserve(m_rows * m_cols);

  for (const auto &row : rows)
    m_data.insert(m_data.end(), row.begin(), row.end());
}

template <class _T>
npp::Mtrx<_T>::Mtrx(const npp::Mtrx<_T> &other) = default;

template <class _T>
npp::Mtrx<_T>::Mtrx(npp::Mtrx<_T> &&other) noexcept = default;

template <class _T>
npp::Mtrx<_T>::~Mtrx() noexcept = default;

template <class _T>
typename npp::Mtrx<_T>::SizeType npp::Mtrx<_T>::index(const npp::Mtrx<_T>::SizeType row, const npp::Mtrx<_T>::SizeType col) const noexcept
{
  return col + row * cols();
}

template <class _T>
void npp::Mtrx<_T>::checkRowBounds(const npp::Mtrx<_T>::SizeType row) const
{
  if (row >= rows())
    throw std::invalid_argument("Row " + std::to_string(row) + " exceeds matrix row bounds " + std::to_string(rows()));
}

template <class _T>
void npp::Mtrx<_T>::checkColBounds(const npp::Mtrx<_T>::SizeType col) const
{
  if (col >= cols())
    throw std::invalid_argument("Column " + std::to_string(col) + " exceeds matrix column bounds " + std::to_string(cols()));
}

template <class _T>
void npp::Mtrx<_T>::checkBounds(const npp::Mtrx<_T>::SizeType row, const npp::Mtrx<_T>::SizeType col) const
{
  checkRowBounds(row);
  checkColBounds(col);
}

template <class _T>
void npp::Mtrx<_T>::checkSizeCompatibility(const npp::Mtrx<_T>::SizeType r0, const npp::Mtrx<_T>::SizeType c0, const npp::Mtrx<_T>::SizeType r1, const npp::Mtrx<_T>::SizeType c1) const
{
  if (r0 != r1 || c0 != c1)
    throw std::invalid_argument("Matrix dimension mismatch: " + std::to_string(r0) + 'x' + std::to_string(c0) + " vs " + std::to_string(r1) + 'x' + std::to_string(c1));
}

template <class _T>
void npp::Mtrx<_T>::checkSizeCompatibility(const npp::Mtrx<_T> &other) const
{
  checkSizeCompatibility(rows(), cols(), other.rows(), other.cols());
}

template <class _T>
void npp::Mtrx<_T>::checkMultCompatibility(const npp::Mtrx<_T>::SizeType r0, const npp::Mtrx<_T>::SizeType c0, const npp::Mtrx<_T>::SizeType r1, const npp::Mtrx<_T>::SizeType c1) const
{
  if (c0 != r1)
    throw std::invalid_argument("Matrix dimension mismatch for multiplication: " + std::to_string(r0) + 'x' + std::to_string(c0) + " vs " + std::to_string(r1) + 'x' + std::to_string(c1));
}

template <class _T>
void npp::Mtrx<_T>::checkMultCompatibility(const npp::Mtrx<_T> &other) const
{
  checkMultCompatibility(rows(), cols(), other.rows(), other.cols());
}

template <class _T>
void npp::Mtrx<_T>::checkSquareMatrix(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols) const
{
  if (rows != cols)
    throw std::invalid_argument("Matrix is not a square matrix: it has " + std::to_string(rows) + " rows and " + std::to_string(cols) + " columns");
}

template <class _T>
void npp::Mtrx<_T>::checkSquareMatrix(const npp::Mtrx<_T> &m) const
{
  checkSquareMatrix(m.rows(), m.cols());
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator=(const npp::Mtrx<_T> &other) = default;

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator=(npp::Mtrx<_T> &&other) noexcept = default;

template <class _T>
void npp::Mtrx<_T>::fill(const _T &value)
{
  std::fill(begin(), end(), value);
}

template <class _T>
_T &npp::Mtrx<_T>::at(const npp::Mtrx<_T>::SizeType i, const npp::Mtrx<_T>::SizeType j)
{
  return m_data.at(index(i, j));
}

template <class _T>
const _T &npp::Mtrx<_T>::at(const npp::Mtrx<_T>::SizeType i, const npp::Mtrx<_T>::SizeType j) const
{
  return m_data.at(index(i, j));
}

template <class _T>
npp::Mtrx<_T>::RowProxy::RowProxy(_T *rowData, const npp::Mtrx<_T>::SizeType cols)
: m_rowData(rowData), m_cols(cols) {}

template <class _T>
_T &npp::Mtrx<_T>::RowProxy::operator[](const npp::Mtrx<_T>::SizeType col) noexcept
{
  return m_rowData[col];
}

template <class _T>
const _T &npp::Mtrx<_T>::RowProxy::operator[](const npp::Mtrx<_T>::SizeType col) const noexcept
{
  return m_rowData[col];
}

template <class _T>
typename npp::Mtrx<_T>::RowProxy npp::Mtrx<_T>::operator[](const npp::Mtrx<_T>::SizeType row) noexcept
{
  return RowProxy(&m_data[index(row, 0)], m_cols);
}

template <class _T>
npp::Mtrx<_T>::ConstRowProxy::ConstRowProxy(const _T *rowData, const npp::Mtrx<_T>::SizeType cols)
: m_rowData(rowData), m_cols(cols) {}

template <class _T>
const _T &npp::Mtrx<_T>::ConstRowProxy::operator[](const npp::Mtrx<_T>::SizeType col) const noexcept
{
  return m_rowData[col];
}

template <class _T>
typename npp::Mtrx<_T>::ConstRowProxy npp::Mtrx<_T>::operator[](const npp::Mtrx<_T>::SizeType row) const noexcept
{
  return ConstRowProxy(&m_data[index(row, 0)], m_cols);
}

template <class _T>
typename npp::Mtrx<_T>::SizeType npp::Mtrx<_T>::rows() const noexcept
{
  return m_rows;
}

template <class _T>
typename npp::Mtrx<_T>::SizeType npp::Mtrx<_T>::cols() const noexcept
{
  return m_cols;
}

template <class _T>
std::vector<typename npp::Mtrx<_T>::SizeType> npp::Mtrx<_T>::dims() const noexcept
{
  return {m_rows, m_cols};
}

template <class _T>
typename npp::Mtrx<_T>::SizeType npp::Mtrx<_T>::elements() const noexcept
{
  return m_rows * m_cols;
}

template <class _T>
bool npp::Mtrx<_T>::isEmpty() const noexcept
{
  return m_data.empty();
}

template <class _T>
bool npp::Mtrx<_T>::isSquare() const noexcept
{
  return m_rows == m_cols;
}

template <class _T>
void npp::Mtrx<_T>::resize(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols)
{
  m_rows = rows;
  m_cols = cols;

  m_data.resize(rows * cols);
}

template <class _T>
void npp::Mtrx<_T>::resize(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols, const _T &value)
{
  m_rows = rows;
  m_cols = cols;

  m_data.resize(rows * cols, value);
}

template <class _T>
typename npp::Mtrx<_T>::Iterator npp::Mtrx<_T>::begin() noexcept
{
  return m_data.begin();
}

template <class _T>
typename npp::Mtrx<_T>::ConstIterator npp::Mtrx<_T>::begin() const noexcept
{
  return m_data.begin();
}

template <class _T>
typename npp::Mtrx<_T>::ConstIterator npp::Mtrx<_T>::cbegin() const noexcept
{
  return m_data.cbegin();
}

template <class _T>
typename npp::Mtrx<_T>::Iterator npp::Mtrx<_T>::end() noexcept
{
  return m_data.end();
}

template <class _T>
typename npp::Mtrx<_T>::ConstIterator npp::Mtrx<_T>::end() const noexcept
{
  return m_data.end();
}

template <class _T>
typename npp::Mtrx<_T>::ConstIterator npp::Mtrx<_T>::cend() const noexcept
{
  return m_data.cend();
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator+=(const npp::Mtrx<_T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::plus<_T>{});

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator+=(const _T &other)
{
  checkSizeCompatibility(other);

  for (auto &element : m_data)
    element += other;

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator-=(const npp::Mtrx<_T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::minus<_T>{});

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator-=(const _T &other)
{
  checkSizeCompatibility(other);

  for (auto &element : m_data)
    element -= other;

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator*=(const npp::Mtrx<_T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::multiplies<_T>{});

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator*=(const _T &other)
{
  checkSizeCompatibility(other);

  for (auto &element : m_data)
    element *= other;

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator/=(const npp::Mtrx<_T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::divides<_T>{});

  return *this;
}

template <class _T>
npp::Mtrx<_T> &npp::Mtrx<_T>::operator/=(const _T &other)
{
  checkSizeCompatibility(other);

  for (auto &element : m_data)
    element /= other;

  return *this;
}

template <class _T>
const npp::Mtrx<_T> &npp::Mtrx<_T>::operator+() const
{
  return *this;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator+(const npp::Mtrx<_T> &other) const
{
  npp::Mtrx<_T> result(*this);

  result += other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator+(const _T &other) const
{
  npp::Mtrx<_T> result(*this);

  result += other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator-() const
{
  npp::Mtrx<_T> result(rows(), cols());

  std::transform(begin(), end(), result.begin(), std::negate<_T>{});

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator-(const npp::Mtrx<_T> &other) const
{
  npp::Mtrx<_T> result(*this);

  result -= other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator-(const _T &other) const
{
  npp::Mtrx<_T> result(*this);

  result -= other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator*(const npp::Mtrx<_T> &other) const
{
  npp::Mtrx<_T> result(*this);

  result *= other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator*(const _T &other) const
{
  npp::Mtrx<_T> result(*this);

  result *= other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator/(const npp::Mtrx<_T> &other) const
{
  npp::Mtrx<_T> result(*this);

  result /= other;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::operator/(const _T &other) const
{
  npp::Mtrx<_T> result(*this);

  result /= other;

  return result;
}

template <class _T>
bool npp::Mtrx<_T>::operator==(const npp::Mtrx<_T> &other) const
{
  return
    rows() == other.rows()
    && cols() == other.cols()
    && std::equal(begin(), end(), other.begin());
}

template <class _T>
bool npp::Mtrx<_T>::operator==(const _T &other) const
{
  return std::find(begin(), end(), other) != end();
}

template <class _T>
bool npp::Mtrx<_T>::operator!=(const npp::Mtrx<_T> &other) const
{
  return !(*this == other);
}

template <class _T>
bool npp::Mtrx<_T>::operator!=(const _T &other) const
{
  return !(*this == other);
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::T() const
{
  npp::Mtrx<_T> result(cols(), rows());

  for (SizeType i = 0; i < rows(); i++)
    for (SizeType j = 0; j < cols(); j++)
      result[j][i] = m_data[index(i, j)];

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::t() const
{
  return this->_T();
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::H() const
{
  return this->_T();
}

template <>
npp::Mtrx<npp::Cmplx> npp::Mtrx<npp::Cmplx>::H() const
{
  npp::Mtrx<npp::Cmplx> result(this->T()); // AGGIUNGERE conj()

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::h() const
{
  return this->H();
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::dot(const npp::Mtrx<_T> &other) const
{
  checkMultCompatibility(other);

  npp::Mtrx<_T> result(rows(), other.cols());

  for (SizeType i = 0; i < rows(); i++)
    for (SizeType j = 0; j < other.cols(); j++)
    {
      result[i][j] = 0;
      for (SizeType k = 0; k < cols(); k++)
        result[i][j] += m_data[index(i, k)] * other[k][j];
    }

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::minor(const npp::Mtrx<_T>::SizeType i, const npp::Mtrx<_T>::SizeType j) const
{
  npp::Mtrx<_T> result(rows() - 1, cols() - 1);

  for (npp::Mtrx<_T>::SizeType r = 0; r < i; r++)
  {
    for (npp::Mtrx<_T>::SizeType c = 0; c < j; c++)
    {
      result[r, c] = m_data[index(r, c)];
    }
    for (npp::Mtrx<_T>::SizeType c = j; c < cols() - 1; c++)
    {
      result[r, c] = m_data[index(r, c + 1)];
    }
  }
  for (npp::Mtrx<_T>::SizeType r = i; r < rows() - 1; r++)
  {
    for (npp::Mtrx<_T>::SizeType c = 0; c < j; c++)
    {
      result[r, c] = m_data[index(r + 1, c)];
    }
    for (npp::Mtrx<_T>::SizeType c = j; c < cols() - 1; c++)
    {
      result[r, c] = m_data[index(r + 1, c + 1)];
    }
  }

  return result;
}

template <class _T>
typename npp::Mtrx<_T>::ValueType npp::Mtrx<_T>::det() const
{
  checkSquareMatrix(*this);

  if (rows() == 1)
  {
    return m_data[index(0, 0)];
  }
  else if (rows() == 2)
  {
    return m_data[index(0, 0)] * m_data[index(1, 1)] - m_data[index(0, 1)] * m_data[index(1, 0)];
  }
  else
  {
    npp::Mtrx<_T>::ValueType result = 0;

    for (npp::Mtrx<_T>::SizeType k = 0; k < cols(); k++)
      result += std::pow(-1.0, k) * m_data[index(0, k)] * minor(0, k).det();
  
    return result;
  }

  return -1;
}

template <class _T>
typename npp::Mtrx<_T>::ValueType npp::Mtrx<_T>::det(const npp::Mtrx<_T> &m)
{
  return m.det();
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::adj() const
{
  npp::Mtrx<_T> result(rows(), cols());

  for (npp::Mtrx<_T>::SizeType i = 0; i < rows(); i++)
  {
    for (npp::Mtrx<_T>::SizeType j = 0; j < cols(); j++)
    {
      result[i][j] = std::pow(-1.0, i + j) * minor(i, j).det();
    }
  }
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::inv() const
{
  return adj() / det();
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::lInv(const npp::Mtrx<_T> &other) const
{
  return inv().dot(other);
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::rInv(const npp::Mtrx<_T> &other) const
{
  return other.dot(inv());
}

template <class _T>
typename npp::Mtrx<_T>::ValueType npp::Mtrx<_T>::sum() const
{
  return std::accumulate(begin(), end(), _T{});
}

template <class _T>
typename npp::Mtrx<_T>::ValueType npp::Mtrx<_T>::min() const
{
  if (isEmpty())
    throw std::runtime_error("min() called on empty matrix");

  return *std::min_element(begin(), end());
}

template <class _T>
typename npp::Mtrx<_T>::ValueType npp::Mtrx<_T>::max() const
{
  if (isEmpty())
    throw std::runtime_error("max() called on empty matrix");

  return *std::max_element(begin(), end());
}

template <class _T>
template <class UnaryFunc>
npp::Mtrx<_T> &npp::Mtrx<_T>::apply(UnaryFunc func)
{
  std::transform(begin(), end(), begin(), func);

  return *this;
}

template <class _T>
template <class UnaryFunc>
npp::Mtrx<_T> npp::Mtrx<_T>::map(UnaryFunc func) const
{
  npp::Mtrx<_T> result(*this);

  result.apply(func);

  return result;
}

template <class _T>
void npp::Mtrx<_T>::print() const
{
  for (npp::Mtrx<_T>::SizeType i = 0; i < rows(); i++)
  {
    std::cout << '\n';
    for (npp::Mtrx<_T>::SizeType j = 0; j < cols(); j++)
      std::cout << m_data[index(i, j)] << '\t';
  }
  std::cout << '\n';
}

template <class _T>
void npp::Mtrx<_T>::print(const std::string &name) const
{
  std::cout << '\n' << name << " =";

  print();
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::zeros(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols)
{
  npp::Mtrx<_T> result(rows, cols);

  for (auto &element : result)
    element = 0.0;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::ones(const npp::Mtrx<_T>::SizeType rows, const npp::Mtrx<_T>::SizeType cols)
{
  npp::Mtrx<_T> result(rows, cols);

  for (auto &element : result)
    element = 1.0;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::diag(const npp::Mtrx<_T>::SizeType dim, const npp::Mtrx<_T>::ValueType &value)
{
  npp::Mtrx<_T> result(dim, dim);

  for (npp::Mtrx<_T>::SizeType k = 0; k < dim; k++)
    result[k][k] = value;

  return result;
}

template <class _T>
npp::Mtrx<_T> npp::Mtrx<_T>::diag(const npp::Vctr<_T> &values)
{
  npp::Mtrx<_T> result(values.size(), result.size());

  for (npp::Mtrx<_T>::SizeType k = 0; k < values.size(); k++)
    result[k][k] = values[k];

  return result;
}