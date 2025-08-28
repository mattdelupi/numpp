#include "nppCmplx.h"

npp::Cmplx::Cmplx()
: m_real(0.0), m_imag(0.0) {}

npp::Cmplx::Cmplx(const double other)
: m_real(other), m_imag(0.0) {}

npp::Cmplx::Cmplx(const double re, const double im)
: m_real(re), m_imag(im) {}

npp::Cmplx::Cmplx(const npp::Cmplx &other)
: m_real(other.re()), m_imag(other.im()) {}

npp::Cmplx npp::Cmplx::cartesian(const double re, const double im)
{
  return Cmplx(re, im);
}

npp::Cmplx npp::Cmplx::polar(const double r, const double t)
{
  if (r < 0)
    std::invalid_argument("complex number magnitude must be greater than zero");

  return Cmplx(r * std::cos(t), r * std::sin(t));
}

npp::Cmplx npp::Cmplx::i()
{
  return Cmplx(0.0, 1.0);
}

npp::Cmplx npp::Cmplx::zero()
{
  return Cmplx(0.0, 0.0);
}

npp::Cmplx npp::Cmplx::one()
{
  return Cmplx(1.0, 0.0);
}

double &npp::Cmplx::re()
{
  return m_real;
}

const double &npp::Cmplx::re() const
{
  return m_real;
}

double &npp::Cmplx::im()
{
  return m_imag;
}

const double &npp::Cmplx::im() const
{
  return m_imag;
}

double npp::Cmplx::mag() const
{
  return std::sqrt(m_real * m_real + m_imag * m_imag);
}

double npp::Cmplx::arg() const
{
  return std::atan2(m_imag, m_real);
}

npp::Cmplx &npp::Cmplx::operator=(const npp::Cmplx &other)
{
  if (&other == this)
    return *this;

  m_real = other.re();
  m_imag = other.im();
  return *this;
}

npp::Cmplx &npp::Cmplx::operator=(const double other)
{
  m_real = other;
  m_imag = 0;
  return *this;
}

const npp::Cmplx &npp::Cmplx::operator+() const
{
  return *this;
}

npp::Cmplx npp::Cmplx::operator+(const npp::Cmplx &other) const
{
  return Cmplx(m_real + other.re(), m_imag + other.im());
}

npp::Cmplx npp::Cmplx::operator+(const double other) const
{
  return Cmplx(m_real + other, m_imag);
}

npp::Cmplx operator+(const double other, const npp::Cmplx &z)
{
  return z + other;
}

npp::Cmplx npp::Cmplx::operator-() const
{
  return Cmplx(-m_real, -m_imag);
}

npp::Cmplx npp::Cmplx::operator-(const npp::Cmplx &other) const
{
  return Cmplx(m_real - other.re(), m_imag - other.im());
}

npp::Cmplx npp::Cmplx::operator-(const double other) const
{
  return Cmplx(m_real - other, m_imag);
}

npp::Cmplx operator-(const double other, const npp::Cmplx &z)
{
  return -(z - other);
}

npp::Cmplx npp::Cmplx::operator*(const npp::Cmplx &other) const
{
  return Cmplx
  (
    m_real * other.re() - m_imag * other.im(),
    m_real * other.im() + m_imag * other.re()
  );
}

npp::Cmplx npp::Cmplx::operator*(const double other) const
{
  return Cmplx(m_real * other, m_imag * other);
}

npp::Cmplx operator*(const double other, const npp::Cmplx &z)
{
  return z * other;
}

npp::Cmplx npp::Cmplx::operator/(const npp::Cmplx &other) const
{
  npp::checkDivision(other);

  return Cmplx
  (
    (m_real * other.re() + m_imag * other.im()) / (other.re() * other.re() + other.im() * other.im()),
    (m_imag * other.re() - m_real * other.im()) / (other.re() * other.re() + other.im() * other.im())
  );
}

npp::Cmplx npp::Cmplx::operator/(const double other) const
{
  npp::checkDivision(other);

  return Cmplx(m_real / other, m_imag / other);
}

npp::Cmplx operator/(const double other, const npp::Cmplx &z)
{
  npp::checkDivision(other);

  return npp::Cmplx
  (
    (z.re() * other) / (z.re() * z.re() + z.im() * z.im()),
    (-z.im() * other) / (z.re() * z.re() + z.im() * z.im())
  );
}

npp::Cmplx npp::Cmplx::operator+=(const npp::Cmplx &other)
{
  *this = *this + other;

  return *this;
}

npp::Cmplx npp::Cmplx::operator+=(const double &other)
{
  *this = *this + other;

  return *this;
}

npp::Cmplx operator+=(const double &other, const npp::Cmplx &z)
{
  npp::Cmplx result = other + z;

  return result;
}

npp::Cmplx npp::Cmplx::operator-=(const npp::Cmplx &other)
{
  *this = *this - other;

  return *this;
}

npp::Cmplx npp::Cmplx::operator-=(const double &other)
{
  *this = *this - other;

  return *this;
}

npp::Cmplx operator-=(const double &other, const npp::Cmplx &z)
{
  npp::Cmplx result = other - z;

  return result;
}

npp::Cmplx npp::Cmplx::operator*=(const npp::Cmplx &other)
{
  *this = *this * other;

  return *this;
}

npp::Cmplx npp::Cmplx::operator*=(const double &other)
{
  *this = *this * other;

  return *this;
}

npp::Cmplx operator*=(const double &other, const npp::Cmplx &z)
{
  npp::Cmplx result = other * z;

  return result;
}

npp::Cmplx npp::Cmplx::operator/=(const npp::Cmplx &other)
{
  npp::checkDivision(other);

  *this = *this / other;

  return *this;
}

npp::Cmplx npp::Cmplx::operator/=(const double &other)
{
  npp::checkDivision(other);

  *this = *this / other;

  return *this;
}

npp::Cmplx operator/=(const double &other, const npp::Cmplx &z)
{
  npp::checkDivision(z);

  npp::Cmplx result = other / z;

  return result;
}

bool npp::Cmplx::operator==(const npp::Cmplx &other) const
{
  if (m_real != other.re() || m_imag != other.im())
    return false;
  else
    return true;
}

bool npp::Cmplx::operator==(const double other) const
{
  if (m_real != other || m_imag != 0.0)
    return false;
  else
    return true;
}

bool operator==(const double other, const npp::Cmplx &z)
{
  return z == other;
}

bool npp::Cmplx::operator!=(const npp::Cmplx &other) const
{
  return !(*this == other);
}

bool npp::Cmplx::operator!=(const double other) const
{
  return !(*this == other);
}

bool operator!=(const double other, const npp::Cmplx &z)
{
  return z != other;
}

npp::Cmplx npp::Cmplx::conj() const
{
  return Cmplx(m_real, -m_imag);
}

npp::Cmplx npp::Cmplx::pow(const double e) const
{
  return Cmplx::polar(std::pow(mag(), e), e * arg());
}

npp::Cmplx npp::Cmplx::exp() const
{
  return Cmplx
  (
    std::exp(m_real) * std::cos(m_imag),
    std::exp(m_real) * std::sin(m_imag)
  );
}

npp::Cmplx npp::Cmplx::exp(const npp::Cmplx &z)
{
  return z.exp();
}

void npp::Cmplx::print() const
{
  std::cout << m_real;

  if (m_imag >= 0.0)
    std::cout << " + i" << std::abs(m_imag) << '\n';
  else
    std::cout << " - i" << -m_imag << '\n';
}

void npp::Cmplx::print(const npp::Cmplx &z)
{
  z.print();
}

void npp::Cmplx::printPolar() const
{
  std::cout << mag() << '<' << arg() << '\n';
}

void npp::Cmplx::printPolar(const Cmplx &z)
{
  z.printPolar();
}