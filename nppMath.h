#pragma once

#include <stdexcept>

namespace npp
{
  template <class T>
  void checkDivision(const T &value)
  {
    if (value == T{})
    throw (std::runtime_error("division by zero occurred"));
  }
}