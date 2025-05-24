#pragma once

#include <iostream>
#include <vector>

namespace dlp {
	template <class _T>
	class vector {
	private:
		std::size_t length;
		_T* ref;
	public:
		vector();
		vector(const std::size_t&);
		vector(const std::size_t&, const _T&);
		vector(const std::size_t&, const _T*);
		vector(const vector<_T>&);
		~vector();
		void print() const;
		std::size_t size() const;
		const _T& operator[](const std::size_t&) const;
		_T& operator[](const std::size_t&);
		const _T& at(const std::size_t&) const;
		_T& at(const std::size_t&);
		_T dot(const vector<_T>&) const;
		_T& dot(const vector<_T>&, _T&) const;
		vector<_T> add(const vector<_T>&) const;
		vector<_T>& add(const vector<_T>&, vector<_T>&);
		vector<_T> add(const _T&) const;
		vector<_T>& add(const _T&, vector<_T>&);
		vector<_T> negate() const;
		vector<_T>& negate(vector<_T>&);
		vector<_T> multiply(const vector<_T>&) const;
		vector<_T>& multiply(const vector<_T>&, vector<_T>&);
		vector<_T> multiply(const _T&) const;
		vector<_T>& multiply(const _T&, vector<_T>&);
		void reset();
		void reset(const std::size_t&);
		void reset(const std::size_t&, const _T&);
		void reset(const vector<_T>&);
		void append(const _T&);
	};

	template <class _T>
	vector<_T>::vector() : length(0), ref(nullptr) {}

	template <class _T>
	vector<_T>::vector(const std::size_t& size) : length(size > 0 ? size : 0), ref(size > 0 ? new _T[size] : nullptr) {}

	template <class _T>
	vector<_T>::vector(const std::size_t& size, const _T& value) : length(size > 0 ? size : 0), ref(size > 0 ? new _T[size] : nullptr) {
		for (std::size_t i = 0; i < length; i++)
			ref[i] = value;
	}

	template <class _T>
	vector<_T>::vector(const std::size_t& size, const _T* array) : length(size > 0 ? size : 0), ref(size > 0 ? new _T[size] : nullptr) {
		for (std::size_t i = 0; i < length; i++)
			ref[i] = array[i];
	}

	template <class _T>
	vector<_T>::vector(const vector<_T>& rhs) : length(rhs.size()), ref(length > 0 ? new _T[length] : nullptr) {
		for (std::size_t i = 0; i < length; i++)
			ref[i] = rhs[i];
	}

	template <class _T>
	vector<_T>::~vector() {
		if (ref != nullptr)
			delete[] ref;
	}

	template <class _T>
	void vector<_T>::print() const {
		for (std::size_t i = 0; i < length; i++)
			std::cout << ref[i] << " ";
		std::cout << "\n";
	}

	template <class _T>
	std::size_t vector<_T>::size() const {
		return length;
	}

	template <class _T>
	const _T& vector<_T>::operator[](const std::size_t& i) const {
		return ref[i];
	}

	template <class _T>
	_T& vector<_T>::operator[](const std::size_t& i) {
		return ref[i];
	}

	template <class _T>
	const _T& vector<_T>::at(const std::size_t& i) const {
		return ref[i];
	}

	template <class _T>
	_T& vector<_T>::at(const std::size_t& i) {
		return ref[i];
	}

	template <class _T>
	_T vector<_T>::dot(const vector<_T>& other) const {
		if (length == other.size()) {
			_T temp = 0;
			for (std::size_t i = 0; i < length; i++)
				temp += ref[i] * other[i];
			return temp;
		}
		else {
			throw ("Sizes must match in order to get the dot product.");
		}
	}

	template <class _T>
	_T& vector<_T>::dot(const vector<_T>& other, _T& out_var) const {
		if (length == other.size()) {
			out_var = 0;
			for (std::size_t i = 0; i < length; i++)
				out_var += ref[i] * other[i];
			return out_var;
		}
		else {
			throw ("Sizes must match in order to get the dot product.");
		}
	}
	
	template <class _T>
	vector<_T> vector<_T>::add(const vector<_T>& other) const {
		if (length == other.size()) {
			vector<_T> temp(length);
			for (std::size_t i = 0; i < length; i++)
				temp[i] = ref[i] + other[i];
			return temp;
		}
		else {
			throw ("Sizes must match in order to add the vectors.");
		}
	}

	template <class _T>
	vector<_T>& vector<_T>::add(const vector<_T>& other, vector<_T>& out_var) {
		if (length == other.size()) {
			for (std::size_t i = 0; i < length; i++)
				out_var[i] = ref[i] + other[i];
			return out_var;
		}
	}

	template <class _T>
	vector<_T> vector<_T>::add(const _T& other) const {
		vector<_T> temp(length);
		for (std::size_t i = 0; i < length; i++)
			temp[i] = ref[i] + other;
		return temp;
	}

	template <class _T>
	vector<_T>& vector<_T>::add(const _T& other, vector<_T>& out_var) {
		for (std::size_t i = 0; i < length; i++)
			out_var[i] = ref[i] + other;
		return out_var;
	}

	template <class _T>
	vector<_T> vector<_T>::negate() const {
		vector<_T> temp(length);
		for (std::size_t i = 0; i < length; i++)
			temp[i] = -ref[i];
		return temp;
	}

	template <class _T>
	vector<_T>& vector<_T>::negate(vector<_T>& out_var) {
		for (std::size_t i = 0; i < length; i++)
			out_var[i] = -ref[i];
		return out_var;
	}

	template <class _T>
	vector<_T> vector<_T>::multiply(const vector<_T>& other) const {
		if (length == other.size()) {
			vector<_T> temp(length);
			for (std::size_t i = 0; i < length; i++)
				temp[i] = ref[i] * other[i];
			return temp;
		}
		else {
			throw ("Sizes must match in order to multiply the vectors.");
		}
	}

	template <class _T>
	vector<_T>& vector<_T>::multiply(const vector<_T>& other, vector<_T>& out_var) {
		if (length == other.size()) {
			for (std::size_t i = 0; i < length; i++)
				out_var[i] = ref[i] * other[i];
			return out_var;
		}
	}

	template <class _T>
	vector<_T> vector<_T>::multiply(const _T& other) const {
		vector<_T> temp(length);
		for (std::size_t i = 0; i < length; i++)
			temp[i] = ref[i] * other;
		return temp;
	}

	template <class _T>
	vector<_T>& vector<_T>::multiply(const _T& other, vector<_T>& out_var) {
		for (std::size_t i = 0; i < length; i++)
			out_var[i] = ref[i] * other;
		return out_var;
	}

	template <class _T>
	void vector<_T>::reset() {
		if (ref != nullptr)
			delete[] ref;
		ref = new _T[length];
	}

	template <class _T>
	void vector<_T>::reset(const std::size_t& size) {
		if (ref != nullptr)
			delete[] ref;
		length = size;
		ref = new _T[length];
	}

	template <class _T>
	void vector<_T>::reset(const std::size_t& size, const _T& value) {
		this->reset(size);
		for (std::size_t i = 0; i < length; i++)
			ref[i] = value;
	}

	template <class _T>
	void vector<_T>::reset(const vector<_T>& rhs) {
		this->reset(rhs.size());
		for (std::size_t i = 0; i < length; i++)
			ref[i] = rhs[i];
	}

	template <class _T>
	void vector<_T>::append(const _T& value) {
		std::vector<_T> temp;
		for (std::size_t i = 0; i < length; i++)
			temp.push_back(ref[i]);
		temp.push_back(value);
		this->reset(++length);
		for (std::size_t i = 0; i < length; i++)
			ref[i] = temp[i];
	}
}
