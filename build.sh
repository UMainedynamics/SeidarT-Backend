#!/bin/bash

# Detect OS
UNAME_S=$(uname -s)

# Set paths based on the operating system
if [[ "$UNAME_S" == "Linux" ]]; then
    INCLUDE_PATH="/usr/local/include"
    LIB_PATH="/usr/local/lib"
    BIN_PATH="src/executables/linux"
elif [[ "$UNAME_S" == "Darwin" ]]; then
    JSON_FORTRAN_PREFIX="/opt/homebrew/Cellar/json-fortran/9.0.2"
    INCLUDE_PATH=$JSON_FORTRAN_PREFIX/include
    LIB_PATH=$JSON_FORTRAN_PREFIX/lib
    BIN_PATH=$JSON_FORTRAN_PREFIX/bin
else
    echo "Unsupported OS: $UNAME_S"
    exit 1
fi

# Compiler and flags
FC=$(which gfortran)
echo "Using GFortran from $FC"

FFLAGS="-I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -g -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -O0 -fcheck=all"

# OpenMP support (optional)
if [[ "$1" == "--openmp" ]]; then
    FFLAGS="$FFLAGS -fopenmp"
fi

# Source files
SRC_DIR="src/fortran"
SOURCES="$SRC_DIR/constants.f08 $SRC_DIR/seidart_types.f08 $SRC_DIR/seidartio.f08 $SRC_DIR/cpmlfdtd.f08 $SRC_DIR/main.f08"
EXECUTABLE="seidartfdtd"

# Create executable directory if it doesn't exist
mkdir -p "$BIN_PATH"

# Compile the executable
echo "Compiling the SeidarT CPML FDTD executable..."
$FC $FFLAGS -o $EXECUTABLE $SOURCES > compile_output.txt 2>&1

if [[ $? -ne 0 ]]; then
    echo "Compilation failed. Check compile_output.txt for details."
    cat compile_output.txt
    exit 1
fi

# Move the executable to the correct folder
echo "Moving the executable to $BIN_PATH..."
mv $EXECUTABLE $BIN_PATH/

echo "Build and installation complete. Executable is located in $BIN_PATH."
