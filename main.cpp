#include <iostream>
#include <cmath>
#include "nppCmplx.h"

int main()
{
  npp::Cmplx a;
  npp::Cmplx b(-1.0, 2.0);
  npp::Cmplx c(2.5);
  npp::Cmplx d(b);

  auto e = npp::Cmplx::cartesian(3.0, 4.0);
  auto f = npp::Cmplx::polar(5.0, std::atan2(4.0, 3.0));

  auto g = npp::Cmplx::i();
  auto h = npp::Cmplx::zero();
  auto i = npp::Cmplx::one();

  npp::Cmplx j = e;
  npp::Cmplx k; k = -2.3;

  npp::Cmplx l = e.conj();

  npp::Cmplx m = f.pow(1.5);

  npp::Cmplx n = m.exp();
  auto o = npp::Cmplx::exp(m);

  // Now prints
  a.print(); a.printPolar();
  b.print(); b.printPolar();
  c.print(); c.printPolar();
  d.print(); d.printPolar();

  e.print(); e.printPolar();
  f.print(); f.printPolar();

  g.print(); g.printPolar();
  h.print(); h.printPolar();
  i.print(); i.printPolar();

  j.print(); j.printPolar();
  k.print(); k.printPolar();

  l.print(); l.printPolar();

  m.print(); m.printPolar();

  n.print(); n.printPolar();
  o.print(); o.printPolar();

  return 0;
}