A high-performance library for computing differential operators on scalar and vector fields, implemented in Fortran with Python bindings.

## Overview

This library provides efficient implementations of common differential operators used in computational fluid dynamics, physics simulations, and other scientific computing applications. The core functionality is implemented in Fortran for performance, with Python bindings for ease of use.

## Features

- **Fast and accurate finite difference methods**:
  - First and second-order derivatives
  - Central, forward, and backward differencing schemes
  - Support for 1D, 2D, and 3D domains

- **Vector field operators**:
  - Gradient
  - Divergence
  - Curl
  - Laplacian and vector Laplacian
  - Jacobian matrices
  - Hessian matrices

- **Parallel computing**:
  - MPI support for distributed computing
  - Domain decomposition

- **Python Interface**:
  - Access all operators from Python
  - NumPy integration
  - Easy-to-use API

## Building the Library

### Prerequisites

- Fortran compiler (gfortran or Intel Fortran)
- MPI implementation (OpenMPI or MPICH)
- CMake (3.22 or newer)
- Python 3.x with NumPy

### Build Instructions

```bash
# Create a build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build
cmake --build .

# Install
cmake --install .
```

## Usage Examples

### Using the Python API

```python
import numpy as np
from operators_lib import OperatorsLib

# Create 3D scalar field
nx, ny, nz = 100, 100, 100
dx = dy = dz = 0.01
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
z = np.linspace(0, 1, nz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Create scalar field f(x,y,z) = x^2 + y^2 + z^2
scalar_field = X**2 + Y**2 + Z**2

# Calculate gradient
gradient = OperatorsLib.calc_gradient(scalar_field, dx, dy, dz)

# Calculate Laplacian
laplacian = OperatorsLib.calc_laplacian(scalar_field, dx, dy, dz)

# Create vector field
vector_field = np.zeros((3, nx, ny, nz))
vector_field[0] = X**2 + Y**2
vector_field[1] = Y**2 + Z**2
vector_field[2] = Z**2 + X**2

# Calculate divergence
divergence = OperatorsLib.calc_divergence(vector_field, dx, dy, dz)

# Calculate curl
curl = OperatorsLib.calc_curl(vector_field, dx, dy, dz)
```

### Using Fortran Directly

```fortran
program example
    use serial_vector_field_operators
    use iso_fortran_env, only: real64
    implicit none
    
    integer, parameter :: nx = 100, ny = 100, nz = 100
    real(real64), parameter :: dx = 0.01, dy = 0.01, dz = 0.01
    real(real64) :: f(nx, ny, nz), df(nx, ny, nz)
    real(real64) :: x, y, z
    integer :: i, j, k
    
    ! Initialize scalar field
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                x = real(i, real64) * dx
                y = real(j, real64) * dy
                z = real(k, real64) * dz
                f(i,j,k) = x**2 + y**2 + z**2
            end do
        end do
    end do
    
    ! Calculate x-derivative
    df = calculate_derivative(f, dx, dy, dz, 1, 1)
    
    ! Calculate Laplacian
    df = laplacian(f, dx, dy, dz)
end program example
```

## Testing

Run the test suite to verify correct operation:

```bash
cd build
ctest
# or directly
./bin/test_operators
```

## License

This library is provided under the MIT License.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
