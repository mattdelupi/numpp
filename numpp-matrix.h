#pragma once

#include <iostream>

namespace dlp {
	template <class _T>
	class matrix {
	private:
		std::size_t rows, cols;
		_T* ref;
	public:
		matrix();
		matrix(const std::size_t&, const std::size_t&);
		matrix(const std::size_t&, const std::size_t&, const _T&);
		matrix(const std::size_t&, const std::size_t&, const _T*);
		matrix(const matrix<_T>&);
		~matrix();
		void print() const;
		std::size_t nrows() const;
		std::size_t ncols() const;
		const _T& at(const std::size_t&, const std::size_t&) const;
		_T& at(const std::size_t&, const std::size_t&);
		/*matrix<_T> dot(const matrix<_T>&) const;
		matrix<_T>& dot(const matrix<_T>&, matrix<_T>&);
		matrix<_T> add(const matrix<_T>&) const;
		matrix<_T>& add(const matrix<_T>&, matrix<_T>&);
		matrix<_T> add(const _T&) const;
		matrix<_T>& add(const _T&, matrix<_T>&);
		matrix<_T> negate() const;
		matrix<_T>& negate(matrix<_T>&);
		matrix<_T> multiply(const matrix<_T>&) const;
		matrix<_T>& multiply(const matrix<_T>&, matrix<_T>&);
		matrix<_T> multiply(const _T&) const;
		matrix<_T>& multiply(const _T&, matrix<_T>&);
		void reset();
		void reset(const std::size_t&, const std::size_t&);
		void reset(const std::size_t&, const std::size_t&, const _T&);
		void reset(const matrix<_T>&);*/
	};

	template <class _T>
	matrix<_T>::matrix() : rows(0), cols(0), ref(nullptr) {}

	template <class _T>
	matrix<_T>::matrix(const std::size_t& m, const std::size_t& n)
		: rows(m > 0 ? m : 0), cols(n > 0 ? n : 0), ref(rows * cols > 0 ? new _T[rows * cols] : nullptr) {}

	template <class _T>
	matrix<_T>::matrix(const std::size_t& m, const std::size_t& n, const _T& value)
		: rows(m > 0 ? m : 0), cols(n > 0 ? n : 0), ref(rows * cols > 0 ? new _T[rows * cols] : nullptr) {
		for (std::size_t i = 0; i < rows; i++)
			for (std::size_t j = 0; j < cols; j++)
				ref[cols * i + j] = value;
	}

	template <class _T>
	matrix<_T>::matrix(const std::size_t& m, const std::size_t& n, const _T* array)
		: rows(m > 0 ? m : 0), cols(n > 0 ? n : 0), ref(rows * cols > 0 ? new _T[rows * cols] : nullptr) {
		for (std::size_t i = 0; i < rows; i++)
			for (std::size_t j = 0; j < cols; j++)
				ref[cols * i + j] = *array++;
	}

	template <class _T>
	matrix<_T>::matrix(const matrix<_T>& rhs) : rows(rhs.nrows()), cols(rhs.ncols()), ref(rows * cols > 0 ? new _T[rows * cols] : nullptr) {
		for (std::size_t i = 0; i < rows; i++)
			for (std::size_t j = 0; j < cols; j++)
				ref[cols * i + j] = rhs.at(i, j);
	}

	template <class _T>
	matrix<_T>::~matrix() {
		if (ref != nullptr)
			delete[] ref;
	}

	template <class _T>
	void matrix<_T>::print() const {
		for (std::size_t i = 0; i < rows; i++) {
			for (std::size_t j = 0; j < cols; j++)
				std::cout << ref[cols * i + j] << " ";
			std::cout << "\n";
		}
	}

	template <class _T>
	std::size_t matrix<_T>::nrows() const {
		return rows;
	}

	template <class _T>
	std::size_t matrix<_T>::ncols() const {
		return cols;
	}

	template <class _T>
	const _T& matrix<_T>::at(const std::size_t& i, const std::size_t& j) const {
		return ref[cols * i + j];
	}

	template <class _T>
	_T& matrix<_T>::at(const std::size_t& i, const std::size_t& j) {
		return ref[cols * i + j];
	}
}
