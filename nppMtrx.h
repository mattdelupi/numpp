#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
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
      operator _T &();
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
      operator const _T &() const;
      const _T &operator[](const SizeType) const noexcept;
    };

  protected:
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

  protected:
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

    Mtrx &operator^=(const Mtrx &);
    Mtrx &operator^=(const _T &);

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

    Mtrx operator^(const Mtrx &) const;
    Mtrx operator^(const _T &) const;

    bool operator==(const Mtrx &) const;
    bool operator==(const _T &) const;

    bool operator!=(const Mtrx &) const;
    bool operator!=(const _T &) const;

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

    Mtrx vectorize() const;

    static Mtrx zeros(const SizeType, const SizeType);
    static Mtrx ones(const SizeType, const SizeType);
    static Mtrx diag(const SizeType, const ValueType &);
    static Mtrx diag(const Mtrx &);
  };
}

template <class _T>
npp::Mtrx<_T> operator+(const _T &, const npp::Mtrx<_T> &);

template <class _T>
npp::Mtrx<_T> operator-(const _T &, const npp::Mtrx<_T> &);

template <class _T>
npp::Mtrx<_T> operator*(const _T &, const npp::Mtrx<_T> &);

template <class _T>
npp::Mtrx<_T> operator/(const _T &, const npp::Mtrx<_T> &);

template <class _T>
npp::Mtrx<_T> operator^(const _T &, const npp::Mtrx<_T> &);

#include "nppMtrx.tpp"