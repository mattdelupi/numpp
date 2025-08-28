#pragma once

#include <iostream>
#include <cmath>
#include <stdexcept>
#include "nppMath.h"

namespace npp
{
  class Cmplx
  {
  private:
    double m_real;
    double m_imag;

  public:
    Cmplx();
    Cmplx(const double);
    explicit Cmplx(const double, const double);
    Cmplx(const Cmplx &);

    static Cmplx cartesian(const double, const double);
    static Cmplx polar(const double, const double);

    static Cmplx i();
    static Cmplx zero();
    static Cmplx one();

  public:
    double &re();
    const double &re() const;
    double &im();
    const double &im() const;

    double mag() const;
    double arg() const;

    Cmplx &operator=(const Cmplx &);
    Cmplx &operator=(const double);

    const Cmplx &operator+() const;
    Cmplx operator+(const Cmplx &) const;
    Cmplx operator+(const double) const;

    Cmplx operator-() const;
    Cmplx operator-(const Cmplx &) const;
    Cmplx operator-(const double) const;

    Cmplx operator*(const Cmplx &) const;
    Cmplx operator*(const double) const;

    Cmplx operator/(const Cmplx &) const;
    Cmplx operator/(const double) const;

    bool operator==(const Cmplx &) const;
    bool operator==(const double) const;

    bool operator!=(const Cmplx &) const;
    bool operator!=(const double) const;

    Cmplx conj() const;

    Cmplx pow(const double) const;

    Cmplx exp() const;
    static Cmplx exp(const Cmplx &);

    void print() const;
    static void print(const Cmplx &);
    void printPolar() const;
    static void printPolar(const Cmplx &);
  };
}

npp::Cmplx operator+(const double, const npp::Cmplx &);
npp::Cmplx operator-(const double, const npp::Cmplx &);
npp::Cmplx operator*(const double, const npp::Cmplx &);
npp::Cmplx operator/(const double, const npp::Cmplx &);
bool operator==(const double, const npp::Cmplx &);
bool operator!=(const double, const npp::Cmplx &);