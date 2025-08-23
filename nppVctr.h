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

namespace npp
{
  template <class T>
  class Vctr
  {
    static_assert(std::is_arithmetic_v<T>, "The npp::Vctr template class requires an arithmetic type like int, float, double, npp::Cmplx, etc.");

  public:
    using ValueType = T;
    using SizeType = typename std::vector<T>::size_type;
    using Iterator = typename std::vector<T>::iterator;
    using ConstIterator = typename std::vector<T>::const_iterator;

  private:
    std::vector<T> m_data;
    bool m_isColumn;

  public:
    Vctr();
    explicit Vctr(const SizeType);
    explicit Vctr(const SizeType, const T &);
    Vctr(std::initializer_list<T>);
    Vctr(const Vctr &);
    Vctr(Vctr &&) noexcept;
    ~Vctr();

    template <class Iter>
    Vctr(Iter, Iter);

  private:
    void checkBounds(const SizeType) const;
    void checkSizeCompatibility(const Vctr &) const;

  public:
    Vctr &operator=(const Vctr &);
    Vctr &operator=(Vctr &&) noexcept;
    void fill(const T &);

    T &at(const SizeType);
    const T &at(const SizeType) const;

    T &operator[](const SizeType) noexcept;
    const T &operator[](const SizeType) const noexcept;

    T &front();
    const T &front() const;

    T &back();
    const T &back() const;

    SizeType size() const noexcept;
    bool isEmpty() const noexcept;
    void reserve(const SizeType);
    SizeType capacity() const noexcept;

    void pushBack(const T &);
    void pushBack(T &&);

    void popBack();
    void clear() noexcept;

    void resize(const SizeType);
    void resize(const SizeType, const T &);

    Iterator begin() noexcept;
    ConstIterator begin() const noexcept;
    ConstIterator cbegin() const noexcept;

    Iterator end() noexcept;
    ConstIterator end() const noexcept;
    ConstIterator cend() const noexcept;

    Vctr &operator+=(const Vctr &);
    Vctr &operator+=(const T &);

    Vctr &operator-=(const Vctr &);
    Vctr &operator-=(const T &);

    Vctr &operator*=(const Vctr &);
    Vctr &operator*=(const T &);

    Vctr &operator/=(const Vctr &);
    Vctr &operator/=(const T &);

    const Vctr &operator+() const;
    Vctr operator+(const Vctr &) const;
    Vctr operator+(const T &) const;

    template <class U>
    friend Vctr<U> operator+(const U &, const Vctr<U> &);

    Vctr operator-() const;
    Vctr operator-(const Vctr &) const;
    Vctr operator-(const T &) const;

    template <class U>
    friend Vctr<U> operator-(const U &, const Vctr<U> &);

    Vctr operator*(const Vctr &) const;
    Vctr operator*(const T &) const;

    template <class U>
    friend Vctr<U> operator*(const U &, const Vctr<U> &);

    Vctr operator/(const Vctr &) const;
    Vctr operator/(const T &) const;

    template <class U>
    friend Vctr<U> operator/(const U &, const Vctr<U> &);

    bool operator==(const Vctr &) const;
    bool operator==(const T &) const;

    template <class U>
    friend bool operator==(const U &, const Vctr<U> &);

    bool operator!=(const Vctr &) const;
    bool operator!=(const T &) const;

    template <class U>
    friend bool operator!=(const U &, const Vctr<U> &);

    T dot(const Vctr &) const;

    T squaredNorm() const;
    T norm() const;

    Vctr &normalize();
    Vctr normalized() const;

    T sum() const;
    T min() const;
    T max() const;

    template <class UnaryFunc>
    Vctr &apply(UnaryFunc);

    template <class UnaryFunc>
    Vctr map(UnaryFunc) const;

    bool isColumn() const;
  };
}