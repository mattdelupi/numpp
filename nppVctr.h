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
  class Vctr
  {
    static_assert(
      std::is_arithmetic_v<T> || std::is_same_v<T, npp::Cmplx>,
      "The npp::Vctr template class requires an arithmetic type like int, float, double, npp::Cmplx, etc."
    );

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
    ~Vctr() noexcept;

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

    Vctr operator-() const;
    Vctr operator-(const Vctr &) const;
    Vctr operator-(const T &) const;

    Vctr operator*(const Vctr &) const;
    Vctr operator*(const T &) const;

    Vctr operator/(const Vctr &) const;
    Vctr operator/(const T &) const;

    bool operator==(const Vctr &) const;
    bool operator==(const T &) const;

    bool operator!=(const Vctr &) const;
    bool operator!=(const T &) const;

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

    static Vctr range(const T &, const T &);
    static Vctr range(const T &, const T &, const T &);

    static Vctr linSpace(const T &, const T &, const SizeType &);

    bool isColumn() const;
  };
}

template <class T, class U>
npp::Vctr<T> operator+(const U &, const npp::Vctr<T> &);

template <class T, class U>
npp::Vctr<T> operator-(const U &, const npp::Vctr<T> &);

template <class T, class U>
npp::Vctr<T> operator*(const U &, const npp::Vctr<T> &);

template <class T, class U>
npp::Vctr<T> operator/(const U &, const npp::Vctr<T> &);

template <class T, class U>
bool operator==(const U &, const npp::Vctr<T> &);

template <class T, class U>
bool operator!=(const U &, const npp::Vctr<T> &);

template <class T>
npp::Vctr<T>::Vctr() : m_data(), m_isColumn(true) {}

template <class T>
npp::Vctr<T>::Vctr(const npp::Vctr<T>::SizeType s) : m_data(s), m_isColumn(true) {}

template <class T>
npp::Vctr<T>::Vctr(const npp::Vctr<T>::SizeType s, const T &value) : m_data(s, value), m_isColumn(true) {}

template <class T>
npp::Vctr<T>::Vctr(std::initializer_list<T> l) : m_data(l), m_isColumn(true) {}

template <class T>
npp::Vctr<T>::Vctr(const npp::Vctr<T> &other) = default;

template <class T>
npp::Vctr<T>::Vctr(npp::Vctr<T> &&other) noexcept = default;

template <class T>
template <class Iter>
npp::Vctr<T>::Vctr(Iter first, Iter last) : m_data(first, last), m_isColumn(true) {}

template <class T>
npp::Vctr<T>::~Vctr() noexcept = default;

template <class T>
void npp::Vctr<T>::checkBounds(const npp::Vctr<T>::SizeType i) const
{
  if (i >= m_data.size())
    throw (std::out_of_range("Index " + std::to_string(i) + " out of range for vector of size " + std::to_string(size())));
}

template <class T>
void npp::Vctr<T>::checkSizeCompatibility(const npp::Vctr<T> &other) const
{
  if (size() != other.size())
    throw (std::invalid_argument("Vector size mismatch: " + std::to_string(size()) + " vs " + std::to_string(other.size())));
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator=(const npp::Vctr<T> &other) = default;

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator=(npp::Vctr<T> &&other) noexcept = default;

template <class T>
void npp::Vctr<T>::fill(const T &value)
{
  std::fill(begin(), end(), value);
}

template <class T>
T &npp::Vctr<T>::at(const npp::Vctr<T>::SizeType i)
{
  return m_data.at(i);
}

template <class T>
const T &npp::Vctr<T>::at(const npp::Vctr<T>::SizeType i) const
{
  return m_data.at(i);
}

template <class T>
T &npp::Vctr<T>::operator[](const npp::Vctr<T>::SizeType i) noexcept
{
  return m_data[i];
}

template <class T>
const T &npp::Vctr<T>::operator[](const npp::Vctr<T>::SizeType i) const noexcept
{
  return m_data[i];
}

template <class T>
T &npp::Vctr<T>::front()
{
  return m_data.front();
}

template <class T>
const T &npp::Vctr<T>::front() const
{
  return m_data.front();
}

template <class T>
T &npp::Vctr<T>::back()
{
  return m_data.back();
}

template <class T>
const T &npp::Vctr<T>::back() const
{
  return m_data.back();
}

template <class T>
typename npp::Vctr<T>::SizeType npp::Vctr<T>::size() const noexcept
{
  return m_data.size();
}

template <class T>
bool npp::Vctr<T>::isEmpty() const noexcept
{
  return m_data.empty();
}

template <class T>
void npp::Vctr<T>::reserve(const npp::Vctr<T>::SizeType s)
{
  m_data.reserve(s);
}

template <class T>
typename npp::Vctr<T>::SizeType npp::Vctr<T>::capacity() const noexcept
{
  return m_data.capacity();
}

template <class T>
void npp::Vctr<T>::pushBack(const T &value)
{
  m_data.push_back(value);
}

template <class T>
void npp::Vctr<T>::pushBack(T &&value)
{
  m_data.push_back(std::move(value));
}

template <class T>
void npp::Vctr<T>::popBack()
{
  m_data.pop_back();
}

template <class T>
void npp::Vctr<T>::clear() noexcept
{
  m_data.clear();
}

template <class T>
void npp::Vctr<T>::resize(const npp::Vctr<T>::SizeType s)
{
  m_data.resize(s);
}

template <class T>
void npp::Vctr<T>::resize(const npp::Vctr<T>::SizeType s, const T &value)
{
  m_data.resize(s, value);
}

template <class T>
typename npp::Vctr<T>::Iterator npp::Vctr<T>::begin() noexcept
{
  return m_data.begin();
}

template <class T>
typename npp::Vctr<T>::ConstIterator npp::Vctr<T>::begin() const noexcept
{
  return m_data.begin();
}

template <class T>
typename npp::Vctr<T>::ConstIterator npp::Vctr<T>::cbegin() const noexcept
{
  return m_data.cbegin();
}

template <class T>
typename npp::Vctr<T>::Iterator npp::Vctr<T>::end() noexcept
{
  return m_data.end();
}

template <class T>
typename npp::Vctr<T>::ConstIterator npp::Vctr<T>::end() const noexcept
{
  return m_data.end();
}

template <class T>
typename npp::Vctr<T>::ConstIterator npp::Vctr<T>::cend() const noexcept
{
  return m_data.cend();
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator+=(const npp::Vctr<T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::plus<T>{});

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator+=(const T &value)
{
  for (auto &element : m_data)
    element += value;

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator-=(const npp::Vctr<T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::minus<T>{});

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator-=(const T &value)
{
  for (auto &element : m_data)
    element -= value;

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator*=(const npp::Vctr<T> &other)
{
  checkSizeCompatibility(other);

  std::transform(begin(), end(), other.begin(), other.end(), std::multiplies<T>{});

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator*=(const T &value)
{
  for (auto &element : m_data)
    element *= value;

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator/=(const npp::Vctr<T> &other)
{
  checkSizeCompatibility(other);
  npp::checkDivision(other);

  std::transform(begin(), end(), other.begin(), begin(), std::divides<T>{});

  return *this;
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::operator/=(const T &value)
{
  npp::checkDivision(value);

  for (auto &element : m_data)
    element /= value;

  return *this;
}

template <class T>
const npp::Vctr<T> &npp::Vctr<T>::operator+() const
{
  return *this;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::operator+(const npp::Vctr<T> &other) const
{
  npp::Vctr<T> result(*this);

  result += other;

  return result;
}

template <class T>npp::Vctr<T> npp::Vctr<T>::operator+(const T &other) const
{
  npp::Vctr<T> result(*this);

  result += other;

  return result;
}

template <class T, class U>
npp::Vctr<T> operator+(const U &other, const npp::Vctr<T> &v)
{
  npp::Vctr<T> result(v);

  result += other;

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::operator-() const
{
  npp::Vctr<T> result(size());

  std::transform(begin(), end(), result.begin(), std::negate<T>{});

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::operator-(const npp::Vctr<T> &other) const
{
  npp::Vctr<T> result(*this);

  result -= other;

  return result;
}

template <class T>npp::Vctr<T> npp::Vctr<T>::operator-(const T &other) const
{
  npp::Vctr<T> result(*this);

  result -= other;

  return result;
}

template <class T, class U>
npp::Vctr<T> operator-(const U &other, const npp::Vctr<T> &v)
{
  npp::Vctr<T> result(v);

  for (auto &element : result)
    element = other - element;

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::operator*(const npp::Vctr<T> &other) const
{
  npp::Vctr<T> result(*this);

  result *= other;

  return result;
}

template <class T>npp::Vctr<T> npp::Vctr<T>::operator*(const T &other) const
{
  npp::Vctr<T> result(*this);

  result *= other;

  return result;
}

template <class T, class U>
npp::Vctr<T> operator*(const U &other, const npp::Vctr<T> &v)
{
  npp::Vctr<T> result(v);

  result *= other;

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::operator/(const npp::Vctr<T> &other) const
{
  npp::Vctr<T> result(*this);

  result /= other;

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::operator/(const T &other) const
{
  npp::Vctr<T> result(*this);

  result /= other;

  return result;
}

template <class T, class U>
npp::Vctr<T> operator/(const U &other, const npp::Vctr<T> &v)
{
  npp::checkDivision(v);

  npp::Vctr<T> result(v);

  for (auto &element : result)
    element = other / element;

  return result;
}

template <class T>
bool npp::Vctr<T>::operator==(const npp::Vctr<T> &other) const
{
  return
    size() == other.size()
    && std::equal(begin(), end(), other.begin());
}

template <class T>
bool npp::Vctr<T>::operator==(const T &other) const
{
  return std::find(begin(), end(), other) != end();
}

template <class T, class U>
bool operator==(const U &other, const npp::Vctr<T> &v)
{
  return v == other;
}

template <class T>
bool npp::Vctr<T>::operator!=(const npp::Vctr<T> &other) const
{
  return !(*this == other);
}

template <class T>
bool npp::Vctr<T>::operator!=(const T &other) const
{
  return !(*this == other);
}

template <class T, class U>
bool operator!=(const U &other, const npp::Vctr<T> &v)
{
  return v != other;
}

template <class T>
T npp::Vctr<T>::dot(const npp::Vctr<T> &other) const
{
  checkSizeCompatibility(other);

  return std::inner_product(begin(), end(), other.begin(), T{});
}

template <class T>
T npp::Vctr<T>::squaredNorm() const
{
  return std::inner_product(begin(), end(), begin(), T{});
}

template <class T>
T npp::Vctr<T>::norm() const
{
  return std::sqrt(squaredNorm());
}

template <class T>
npp::Vctr<T> &npp::Vctr<T>::normalize()
{
  T mag = norm();

  if (mag == T{})
    throw (std::runtime_error("Cannot normalize a zero vector"));
  

  *this /= mag;
  return *this;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::normalized() const
{
  npp::Vctr<T> result(*this);

  return result.normalize();
}

template <class T>
T npp::Vctr<T>::sum() const
{
  return std::accumulate(begin(), end(), T{});
}

template <class T>
T npp::Vctr<T>::min() const
{
  if (isEmpty())
    throw (std::runtime_error("min() called on empty vector"));

  return *std::max_element(begin(), end());
}

template <class T>
T npp::Vctr<T>::max() const
{
  if (isEmpty())
    throw(std::runtime_error("max() called on empty vector"));

  return *std::min_element(begin(), end());
}

template <class T>
template <class UnaryFunc>
npp::Vctr<T> &npp::Vctr<T>::apply(UnaryFunc func)
{
  std::transform(begin(), end(), begin(), func);

  return *this;
}

template <class T>
template <class UnaryFunc>
npp::Vctr<T> npp::Vctr<T>::map(UnaryFunc func) const
{
  npp::Vctr<T> result(*this);

  result.apply(func);

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::range(const T &first, const T &last)
{
  T unit = last >= first ? 1.0 : -1.0;

  npp::Vctr<T>::SizeType N = std::floor(std::abs(last - first)) + 1;

  npp::Vctr<T> result(N);

  for (npp::Vctr<T>::SizeType i = 0; i < N; i++)
    result[i] = first + unit * i;

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::range(const T &first, const T &last, const T &step)
{
  T unit = last >= first ? 1.0 : -1.0;

  npp::Vctr<T>::SizeType N = std::floor(std::abs(last - first) / std::abs(step)) + 1;

  npp::Vctr<T> result(N);

  for (npp::Vctr<T>::SizeType i = 0; i < N; i++)
    result[i] = first + unit * std::abs(step) * i;

  return result;
}

template <class T>
npp::Vctr<T> npp::Vctr<T>::linSpace(const T &first, const T &last, const npp::Vctr<T>::SizeType &length)
{
  npp::Vctr<T> result(length);

  for (npp::Vctr<T>::SizeType i = 0; i < length; i++)
    result[i] = first + (last - first) * i / (length - 1);

  return result;
}

template <class T>
bool npp::Vctr<T>::isColumn() const
{
  return m_isColumn;
}