//Artur Lukianov
#include <iostream>
#include <cmath>
#include <iomanip>


// Gnuplot
#define GNUPLOT_NAME "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist"


//#define ADJUST_ZERO(x) (fabs((x)+0.0) < 0.005f ? 0 : (x))
#define ADJUST_ZERO(x) (x)

const size_t MAXN = 100; // Max matrix size

// Dimensional exception defines exception from dimension mismatch
class DimensionalException: public std::exception {
public:
  // Override the "what" method to return error message
  char * what() {
    return (char *)"Error: the dimensional problem occurred\r\n";
  }
};

// Matrix representation
class Matrix {
public:
  double content[MAXN][MAXN]; // 2D array to store the matrix values
  size_t n,
    m; // Dimensions of matrix. n - number of rows, m - number of columns

  // Default constructor
  Matrix(): n(0),
            m(0) {}

  // Constructor with parameters to initialize the dimensions
  Matrix(size_t n, size_t m): n(n),
                              m(m) {}

  // Assign (copy) operator override
  Matrix & operator = (const Matrix & other) {
    n = other.n;
    m = other.m;

    // Copy the content of the other matrix
    for (size_t i = 0; i < other.n; i++) {
      for (size_t j = 0; j < other.m; j++) {
        content[i][j] = other.content[i][j];
      }
    }

    return * this;
  }

  // Matrix transposition method
  Matrix transpose() {
    Matrix transposed;
    transposed.n = m;
    transposed.m = n;

    // Transpose the content of the matrix
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < m; j++) {
        transposed.content[j][i] = content[i][j];
      }
    }

    return transposed;
  }
};

// Square matrix implementation as child class of Matrix
// Square matrix is basically the same matrix, but with equal dimensions
class SquareMatrix : public Matrix {
public:
  // Default constructor
  SquareMatrix() : Matrix() {}

  // Constructor with dimenstion parameter
  SquareMatrix(size_t n) : Matrix(n, n) {}


  // Assign (copy) operator override
  SquareMatrix & operator = (const SquareMatrix & other) {
    n = other.n;
    m = other.m;

    // Copy the content of the other matrix
    for (size_t i = 0; i < other.n; i++) {
      for (size_t j = 0; j < other.m; j++) {
        content[i][j] = other.content[i][j];
      }
    }

    return * this;
  }


  // Assign (copy) operator override
  SquareMatrix & operator = (const Matrix & other) {
    n = other.n;
    m = other.m;

    // Copy the content of the other matrix
    for (size_t i = 0; i < other.n; i++) {
      for (size_t j = 0; j < other.m; j++) {
        content[i][j] = other.content[i][j];
      }
    }

    return * this;
  }
};


// Implemenation of identity matrix (1 on main diagonal)
class IdentityMatrix: public SquareMatrix {
public:
  // Default constructor for IdentityMatrix
  IdentityMatrix() : SquareMatrix() {}

  // Constructor for IdentityMatrix with specified size
  IdentityMatrix(size_t n) : SquareMatrix(n) {
    for(size_t i = 0; i < n; i++)
      for(size_t j = 0; j < m; j++)
        content[i][j] = 0;

    // Set diagonal elements to 1
    for (size_t i = 0; i < n; i++) {
      content[i][i] = 1;
    }
  }
};

// Implementation of permutation matrix
// By multiplying it by matrix of same size, the specified rows will be swapped
class PermutationMatrix: public SquareMatrix {
public:
  // Default constructor for PermutationMatrix
  PermutationMatrix() : SquareMatrix() {}

  // Constructor for PermutationMatrix with specified size and row indices
  PermutationMatrix(size_t n, size_t r1, size_t r2) : SquareMatrix(n) {
    for(size_t i = 0; i < n; i++)
      for(size_t j = 0; j < m; j++)
        content[i][j] = 0;

    // Set diagonal elements to 1
    for (size_t i = 0; i < n; i++) {
      content[i][i] = 1;
    }

    // Convert row indices to indexation from 0
    r1--; r2--;

    // Swap row r1 and row r2 in permutation matrix
    content[r1][r1] = 0;
    content[r1][r2] = 1;
    content[r2][r2] = 0;
    content[r2][r1] = 1;
  }
};

// Implementation of elimination matrix
// By multiplying it by specified matrix, element specified in constructor will be eliminated
class EliminationMatrix : public SquareMatrix {
public:
  // Default constructor for EliminationMatrix
  EliminationMatrix() : SquareMatrix() {}

  // Constructor for EliminationMatrix with SquareMatrix, row index r, and column index c
  EliminationMatrix(Matrix m, size_t r, size_t c) : SquareMatrix(m.n) {
    for(size_t i = 0; i < m.n; i++)
      for(size_t j = 0; j < m.m; j++)
        content[i][j] = 0;

    // We cannot just be inherited from identity matrix, as it would be a logic error
    // Set diagonal elements to 1
    for (size_t i = 0; i < m.n; i++)
      content[i][i] = 1;

    // Convert row and column indices to indexation from 0
    r--; c--;

    // Find the row with the first nonzero element in column c, except for row r
    for (int i = r - 1; i >= 0; i--) {
      if (m.content[i][c] != 0) {
        // Set the element at position (r,c) to -m[r,c]/m[i,c]
        content[r][c] = -m.content[r][c] / m.content[i][c];
        break;
      }
    }
  }
};


// Implementation of backward elimination matrix
// By multiplying it by specified matrix, element specified in constructor will be eliminated
class BackwardEliminationMatrix : public SquareMatrix {
public:
  // Default constructor for backward elimination matrix
  BackwardEliminationMatrix() : SquareMatrix() {}

  // Constructor for backward elimination matrix with SquareMatrix, row index r, and column index c
  BackwardEliminationMatrix(Matrix m, size_t r, size_t c) : SquareMatrix(m.n) {
    for(size_t i = 0; i < m.n; i++)
      for(size_t j = 0; j < m.m; j++)
        content[i][j] = 0;

    // We cannot just be inherited from identity matrix, as it would be a logic error
    // Set diagonal elements to 1
    for (size_t i = 0; i < m.n; i++)
      content[i][i] = 1;

    // Convert row and column indices to indexation from 0
    r--; c--;

    // Find the row with the first nonzero element in column c, except for row r
    for (int i = r + 1; i < m.n; i++) {
      if (m.content[i][c] != 0) {
        // Set the element at position (r,c) to -m[r,c]/m[i,c]
        content[r][c] = -m.content[r][c] / m.content[i][c];
        break;
      }
    }
  }
};


// Matrix containing two matrices - for gaussian elemnation
class AugmentedMatrix {
public:
  Matrix right; // Right part of augmented matrix
  Matrix left; // Left part of augmented matrix

  AugmentedMatrix() {}
  AugmentedMatrix(Matrix _left, Matrix _right) : right(_right), left(_left) {}
};

// Matrix in column form (1 column, N rows)
class ColumnVector : public Matrix {
public:
  // Default constructor
  ColumnVector() : Matrix() {}
  // Full contructor with number of rows
  ColumnVector(size_t n) : Matrix(n, 1) {}
};

// Plus operator override
Matrix operator + (const Matrix & a,
                   const Matrix & b) {
  Matrix c(a.n, a.m);

  // Check if matrices have the same dimensions
  if (!(a.n == b.n && a.m == b.m)) throw DimensionalException();

  // Add the content of matrices a and b
  for (size_t i = 0; i < a.n; i++) {
    for (size_t j = 0; j < a.m; j++) {
      c.content[i][j] = a.content[i][j] + b.content[i][j];
    }
  }

  return c;
}

// Minus operator override
Matrix operator - (const Matrix & a,
                   const Matrix & b) {
  Matrix c(a.n, a.m);

  // Check if matrices have the same dimensions
  if (a.n != b.n || a.m != b.m) throw DimensionalException();

  // Subtract the content of matrix b from a
  for (size_t i = 0; i < a.n; i++) {
    for (size_t j = 0; j < a.m; j++) {
      c.content[i][j] = a.content[i][j] - b.content[i][j];
    }
  }

  return c;
}

// Matrix multiplication
Matrix operator * (const Matrix & a,
                   const Matrix & b) {
  Matrix c(a.n, b.m);

  // Check if the number of columns in matrix a is equal to the number of rows in matrix b
  if (a.m != b.n) throw DimensionalException();

  // Perform matrix multiplication
  for (size_t i = 0; i < a.n; i++) {
    for (size_t j = 0; j < b.m; j++) {
      c.content[i][j] = 0;
      for (size_t k = 0; k < a.m; k++) {
        c.content[i][j] += a.content[i][k] * b.content[k][j];
      }
      c.content[i][j] = ADJUST_ZERO(c.content[i][j]);
    }
  }

  return c;
}



AugmentedMatrix operator * (const SquareMatrix & a,
                            const AugmentedMatrix & b) {
  return AugmentedMatrix(a * b.left, a * b.right);
}




// IO stream overloading
// Overloads the ">>" operator for the Matrix class to allow input from a std::istream object
std::istream & operator >> (std::istream & in, Matrix & m) {
  // Reads the number of rows and columns of the Matrix object from the std::istream object
  in >> m.n >> m.m;
  // Loop through each row and column of the Matrix object and read in the value from the std::istream object
  for (size_t i = 0; i < m.n; i++) {
  for (size_t j = 0; j < m.m; j++) {
  in >> m.content[i][j];
}
}
  // Returns the std::istream object to allow for chaining of input operations
  return in;
}

// Overloads the "<<" operator for the Matrix class to allow output to a std::ostream object
std::ostream & operator << (std::ostream & out, Matrix & m) {
  // Loop through each row and column of the Matrix object and output the value to the std::ostream object
  for (size_t i = 0; i < m.n; i++) {
    for (size_t j = 0; j < m.m; j++) {
      out << ADJUST_ZERO(m.content[i][j]);
      // Adds a space between each value in a row, except for the last one
      if (j != m.m - 1) std::cout << " ";
    }
    // Starts a new line for each row
    out << "\r\n";
  }
  // Returns the std::ostream object to allow for chaining of output operations
  return out;
}

// Overloads the ">>" operator for the SquareMatrix class to allow input from a std::istream object
std::istream & operator >> (std::istream & in, SquareMatrix & m) {
  // Reads the number of rows and columns of the SquareMatrix object from the std::istream object
  // As dismensions are equal, we can just read one of them
  in >> m.n;
  m.m = m.n;
  // Loop through each row and column of the SquareMatrix object and read in the value from the std::istream object
  for (size_t i = 0; i < m.n; i++) {
  for (size_t j = 0; j < m.m; j++) {
  in >> m.content[i][j];
}
}
  // Returns the std::istream object to allow for chaining of input operations
  return in;
}

// Overloads the "<<" operator for the SquareMatrix class to allow output to a std::ostream object
std::ostream & operator << (std::ostream & out, SquareMatrix & m) {
  // Loop through each row and column of the SquareMatrix object and output the value to the std::ostream object
  for (size_t i = 0; i < m.n; i++) {
    for (size_t j = 0; j < m.m; j++) {
      out << ADJUST_ZERO(m.content[i][j]);
      // Adds a space between each value in a row, except for the last one
      if (j != m.m - 1) std::cout << " ";
    }
    // Starts a new line for each row
    out << "\r\n";
  }
  // Returns the std::ostream object to allow for chaining of output operations
  return out;
}


// Overloads the "<<" operator for the AugmentedMatrix class to allow output to a std::ostream object
std::ostream & operator << (std::ostream & out, AugmentedMatrix & m) {
  // Loop through each row and column of the SquareMatrix object and output the value to the std::ostream object
  for (size_t i = 0; i < m.left.n; i++) {
    // Left part of adjastency matrix
    for (size_t j = 0; j < m.left.m; j++) {
      out << ADJUST_ZERO(m.left.content[i][j]);
      // Adds a space between each value in a row
      std::cout << " ";
    }
    // Right part of adjastency matrix
    for (size_t j = 0; j < m.right.m; j++) {
      out << ADJUST_ZERO(m.right.content[i][j]);
      // Adds a space between each value in a row
      std::cout << " ";
    }
    // Starts a new line for each row
    out << "\r\n";
  }
  // Returns the std::ostream object to allow for chaining of output operations
  return out;
}

// Overloads the ">>" operator for the ColumnVector class to allow input from a std::istream object
std::istream & operator >> (std::istream & in, ColumnVector & m) {
  // Reads the number of rows and columns of the ColumnVector object from the std::istream object
  // As dismensions are equal, we can just read one of them
  in >> m.n;
  m.m = 1;
  // Loop through each row and column of the ColumnVector object and read in the value from the std::istream object
  for (size_t i = 0; i < m.n; i++) {
    for (size_t j = 0; j < m.m; j++) {
      in >> m.content[i][j];
    }
  }
  // Returns the std::istream object to allow for chaining of input operations
  return in;
}


Matrix findInverse(Matrix x) {
  IdentityMatrix ident(x.m);

  AugmentedMatrix a(x, ident);
  unsigned int step = 1;

  // Perform Gaussian elimination
  for (size_t i = 0; i < a.left.n - 1; i++) {
    // Permutate rows to pivot maximal element
    double mx = fabs(a.left.content[i][i]);
    size_t mxi = i;
    for (size_t j = i + 1; j < a.left.n; j++) {
      if (fabs(a.left.content[j][i]) > mx) {
        mx = fabs(a.left.content[j][i]);
        mxi = j;
      }
    }

    // If we need to do permutation to pivot, then perform it
    if (mxi != i) {
      PermutationMatrix pm(a.left.n, i+1, mxi+1);
      a = pm * a;
    }

    // If permutation did not succeded to get non-zero, skip (because everywhere is already 0 on all rows)
    if (a.left.content[i][i] == 0) continue;

    // Eliminate all rows bellow current
    for (size_t j = i + 1; j < a.left.n; j++) {
      if (a.left.content[j][i] == 0) continue;
      EliminationMatrix e(a.left, j+1, i+1);
      a = e * a;
    }
  }

  // Perform Gaussian elimination
  for (int i = a.left.n-1; i >= 0; i--) {
    // Eliminate all rows above current
    for (int j = i - 1; j >= 0; j--) {
      if (a.left.content[j][i] == 0) continue;
      BackwardEliminationMatrix e(a.left, j+1, i+1);
      a = e * a;
    }
  }

  // Perform diagonal normalization
  for (int i = 0; i < a.left.n; i++) {
    if (a.left.content[i][i] == 0) continue;
    for (int j = 0; j < a.left.m; j++) {
      a.right.content[i][j] /= a.left.content[i][i];
    }
    a.left.content[i][i] = 1;
  }

  return a.right;
}


int main() {
    // Interact with gnuplot by file descriptor
    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    // Input data required for LSA
    int m, n;

    std::cout << std::setprecision(4) << std::fixed;

    std::cin >> m;

    Matrix A(m, 2);
    ColumnVector b(m);

    for (int i = 0; i < m; i++) {
        A.content[i][0] = 1;
        std::cin >> A.content[i][1] >> b.content[i][0];
    }

    // Initialize matrix with polynomial values
    std::cin >> n;
    A.m = n + 1;
    for (int i = 0; i < m; i++) {
        for (int j = 2; j < n + 1; j++) {
            A.content[i][j] = A.content[i][j - 1] * A.content[i][1];
        }
    }

    // Output needed values and calculate (AtA)^-1*AtB=x
    std::cout << "A:\n";
    std::cout << A;

    Matrix A2 = A.transpose() * A;
    std::cout << "A_T*A:\n";
    std::cout << A2;

    Matrix A3 = findInverse(A2);
    std::cout << "(A_T*A)^-1:\n";
    std::cout << A3;

    Matrix A4 = A.transpose() * b;
    std::cout << "A_T*b:\n";
    std::cout << A4;

    Matrix x = A3 * A4;
    std::cout << "x~:\n";
    std::cout << x;

    std::string polynomial = "";
    for (int i = 0; i < n + 1; i++) {
        polynomial += std::to_string(x.content[i][0]) + "*x**" + std::to_string(i);
        if (i != n) polynomial += " + ";
    }

    fprintf(pipe, "plot [-5: 25] [-5: 25] %s, '-' using 1:2 with points\n", polynomial.c_str());
    for (int i = 0; i < m; i++) {
        fprintf(pipe, "%f\t%f\n", A.content[i][1], b.content[i][0]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);

    _pclose(pipe);
}
