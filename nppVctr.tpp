template <class _T>
npp::Vctr<_T>::Vctr()
: npp::Mtrx<_T>() {}

template <class _T>
npp::Vctr<_T>::Vctr(const npp::Vctr<_T>::SizeType s)
: npp::Mtrx<_T>(s, 1) {}

template <class _T>
npp::Vctr<_T>::Vctr(const npp::Vctr<_T>::SizeType s, const npp::Vctr<_T>::ValueType &x)
: npp::Mtrx<_T>(s, 1, x) {}

template <class _T>
npp::Vctr<_T>::Vctr(std::initializer_list<_T> elements)
: npp::Mtrx<_T>(elements.size(), 1)
{
  for (npp::Vctr<_T>::SizeType i = 0; i < elements.size(); i++)
    *(this->begin() + i) = *(elements.begin() + i);
}

template <class _T>
npp::Vctr<_T>::Vctr(const npp::Vctr<_T> &other) = default;

template <class _T>
npp::Vctr<_T>::Vctr(npp::Vctr<_T> &&other) noexcept = default;

template <class _T>
npp::Vctr<_T>::Vctr(const npp::Mtrx<_T> &other)
: npp::Mtrx<_T>(other) {}

template <class _T>
npp::Vctr<_T>::Vctr(npp::Mtrx<_T> &&other) noexcept
: npp::Mtrx<_T>(other) {}

template <class _T>
npp::Vctr<_T>::~Vctr() noexcept = default;

template <class _T>
_T &npp::Vctr<_T>::operator[](const npp::Vctr<_T>::SizeType i)
{
  return *(this->begin() + i);
}

template <class _T>
const _T &npp::Vctr<_T>::operator[](const npp::Vctr<_T>::SizeType i) const
{
  return *(this->begin() + i);
}

template <class _T>
npp::Vctr<_T> npp::Vctr<_T>::range(const _T &first, const _T &last)
{
  _T unit = last >= first ? 1.0 : -1.0;

  npp::Vctr<_T>::SizeType N = std::floor(std::abs(last - first)) + 1;

  npp::Vctr<_T> result(N);

  for (npp::Vctr<_T>::SizeType i = 0; i < N; i++)
    result[i] = first + unit * i;

  return result;
}

template <class _T>
npp::Vctr<_T> npp::Vctr<_T>::range(const _T &first, const _T &last, const _T &step)
{
  _T unit = last >= first ? 1.0 : -1.0;

  npp::Vctr<_T>::SizeType N = std::floor(std::abs(last - first) / std::abs(step)) + 1;

  npp::Vctr<_T> result(N);

  for (npp::Vctr<_T>::SizeType i = 0; i < N; i++)
    result[i] = first + unit * std::abs(step) * i;

  return result;
}

template <class _T>
npp::Vctr<_T> npp::Vctr<_T>::linSpace(const _T &first, const _T &last, const npp::Vctr<_T>::SizeType &length)
{
  npp::Vctr<_T> result(length);

  for (npp::Vctr<_T>::SizeType i = 0; i < length; i++)
    result[i] = first + (last - first) * i / (length - 1);

  return result;
}