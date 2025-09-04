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

namespace npp
{
  template <class T>
  class Mtrx
  {
    static_assert(
      std::is_arithmetic_v<T> || std::is_same_v<T, npp::Cmplx> || std::is_same_v<T, npp::Vctr>,
      "The npp::Mtrx template class requires an arithmetic type like int, float, double, npp::Cmplx, etc."
    );

  public:
    using ValueType = T;
    using SizeType = typename std::vector<T>::size_type;
    using Iterator = typename std::vector<T>::iterator;
    using ConstIterator = typename std::vector<T>::const_iterator;

    class RowProxy
    {
    private:
      T *m_rowData;
      SizeType m_cols;

    public:
      RowProxy(const T *, const SizeType);
      T &operator[](const SizeType) noexcept;
      const T &operator[](const SizeType) const noexcept;
    };

    class ConstRowProxy
    {
    private:
      const T *m_rowData;
      SizeType m_cols;

    public:
      ConstRowProxy(const T *, const SizeType);
      const T &operator[](const SizeType) const noexcept;
    };

  private:
    std::vector<T> m_data;
    SizeType m_rows, m_cols;

  public:
    Mtrx();
    explicit Mtrx(const SizeType, const SizeType);
    explicit Mtrx(const SizeType, const SizeType, const T &);
    Mtrx(const Vctr<T> &);
    Mtrx(std::initializer_list<std::initializer_list<T>>);
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
    
  public:
    Mtrx &operator=(const Mtrx &);
    Mtrx &operator=(Mtrx &&) noexcept;
    void fill(const T &);

    T &at(const SizeType, const SizeType);
    const T &at(const SizeType, const SizeType) const;

    RowProxy operator[](const SizeType) noexcept;
    ConstRowProxy operator[](const SizeType) const noexcept;

    SizeType rows() const noexcept;
    SizeType cols() const noexcept;
    std::vector<SizeType> dims() const noexcept;
    SizeType elements() const noexcept;

    bool isEmpty() const noexcept;
    bool isSquare() const noexcept;

    void resize(const SizeType, const SizeType);
    void resize(const SizeType, const SizeType, const T &);

    Iterator begin() noexcept;
    ConstIterator begin() const noexcept;
    ConstIterator cbegin() const noexcept;

    Iterator end() noexcept;
    ConstIterator end() const noexcept;
    ConstIterator cend() const noexcept;

    Mtrx &operator+=(const Mtrx &);
    Mtrx &operator+=(const T &);

    Mtrx &operator-=(const Mtrx &);
    Mtrx &operator-=(const T &);

    Mtrx &operator*=(const Mtrx &);
    Mtrx &operator*=(const T &);

    Mtrx &operator/=(const Mtrx &);
    Mtrx &operator/=(const T &);

    const Mtrx &operator+() const;
    Mtrx operator+(const Mtrx &) const;
    Mtrx operator+(const T &) const;

    Mtrx operator-() const;
    Mtrx operator-(const Mtrx &) const;
    Mtrx operator-(const T &) const;

    Mtrx operator*(const Mtrx &) const;
    Mtrx operator*(const T &) const;

    Mtrx operator/(const Mtrx &) const;
    Mtrx operator/(const T &) const;

    bool operator==(const Mtrx &) const;
    bool operator==(const Mtrx &) const;

    bool operator!=(const Mtrx &) const;
    bool operator!=(const T &);

    Vctr<T> getRow(const SizeType) const;
    Vctr<T> getCol(const SizeType) const;

    Mtrx T() const;
    Mtrx t() const;
    Mtrx H() const;
    Mtrx h() const;

    Mtrx dot(const Mtrx &) const;

    T det() const;

    Mtrx inv() const;
    Mtrx lInv(const Mtrx &) const;
    Mtrx rInv(const Mtrx &) const;

    T sum() const;
    T min() const;
    T max() const;

    template <class UnaryFunc>
    Mtrx &apply(UnaryFunc);

    template <class UnaryFunc>
    Mtrx map(UnaryFunc) const;

    static Mtrx zeros(const SizeType, const SizeType);
    static Mtrx ones(const SizeType, const SizeType);
    static Mtrx diag(const T &) const;
    static Mtrx diag(const Vctr<T> &) const;
  };
}