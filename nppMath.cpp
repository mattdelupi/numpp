#include "nppMath.h"

template <class T>
void npp::checkDivision(const T &value)
{
  if (value == T{})
    throw (std::runtime_error("division by zero occurred"));
}