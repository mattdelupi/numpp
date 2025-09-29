#pragma once

#include "nppMtrx.h"

namespace npp
{
  template <class _T>
  class Vctr : public Mtrx<_T>
  {
  public:
    using ValueType = _T;
    using SizeType = typename std::vector<_T>::size_type;
    using Iterator = typename std::vector<_T>::iterator;
    using ConstIterator = typename std::vector<_T>::const_iterator;

  public:
    Vctr();
    explicit Vctr(const SizeType);
    explicit Vctr(const SizeType, const ValueType &);
    Vctr(std::initializer_list<_T>);
    Vctr(const Vctr &);
    Vctr(Vctr &&) noexcept;
    Vctr(const Mtrx<_T> &);
    Vctr(Mtrx<_T> &&) noexcept;
    ~Vctr() noexcept;

  public:
    _T &operator[](const SizeType);
    const _T &operator[](const SizeType) const;

    static Vctr range(const _T &, const _T &);
    static Vctr range(const _T &, const _T &, const _T &);
    static Vctr linSpace(const _T &, const _T &, const SizeType &);
  };

  using VctrD = Vctr<double>;
  using VctrF = Vctr<float>;
  using VctrC = Vctr<Cmplx>;
  using VctrI = Vctr<int>;
  using VctrL = Vctr<long>;
  using VctrLL = Vctr<long long>;
  using VctrU = Vctr<unsigned>;
  using VctrB = Vctr<bool>;
  using VctrS = Vctr<std::string>;
  using VctrCh = Vctr<char>;
}

#include "nppVctr.tpp"