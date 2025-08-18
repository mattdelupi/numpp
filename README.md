# numpp

library structure:

- complex numbers = Cmplx

- static vector = StackVector
  - constructor
    - content using std::array
    - content using a scalar
    - length
    - copy constructor
  - print
  - length
  - operator[]
    - without reference
    - with reference
  - operator=
    - with another vector
    - with a scalar
  - operator+
    - with another vector
    - with a scalar
  - operator-
    - with another vector
    - with a scalar
    - unary
  - operator*
    - with another vector
    - with a scalar
  - operator/
    - with another vector
    - with a scalar
  - dot
  - cross
  - t, T

- dynamic vector = HeapVector

- static matrix = StackMatrix

- dynamic matrix = HeapMatrix

- static tensor = StackTensor

- dynamic tensor = HeapTensor