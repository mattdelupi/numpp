#pragma once

#include <stdexcept>

namespace npp
{
  template <class T>
  void checkDivision(const T &value)
  {
    if (value == 0.0)
    throw (std::runtime_error("division by zero occurred"));
  }
}