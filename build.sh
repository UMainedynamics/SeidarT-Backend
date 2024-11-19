#!/bin/bash

# -----------------------------------------------------------------------------
user_defined_builddir=false
openmp_compile=false
clean=false
# Loop through the arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -b|--builddir)
      BIN_PATH="$2"
      user_defined_builddir=true
      shift 2
      ;;
    -o|--openmp)
      openmp_compile=true
      shift
      ;;
    -c|--clean)
      clean=true
      shift
      ;;
    -*|--*)
      echo "Unknown option: $1"
      echo "Usage: cmd [-b install_directory] [-c|--clean]"
      exit 1
      ;;
    *)
      echo "Unexpected argument: $1"
      exit 1
      ;;
  esac
done
# -----------------------------------------------------------------------------

# Detect OS
UNAME_S=$(uname -s)
arch=$(uname -m | tr '[:upper:]' '[:lower:]')
os=$(uname -s | tr '[:upper:]' '[:lower:]')

# INCLUDE_PATH="/usr/local/jsonfortran-gnu-9.0.2/lib"
# LIB_PATH="/usr/local/jsonfortran-gnu-9.0.2/lib"

# Set paths based on the operating system
# if [[ "$UNAME_S" == "Linux" ]]; then
#     INCLUDE_PATH="$CONDA_PREFIX/include"
#     LIB_PATH="$CONDA_PREFIX/lib"
#     BIN_PATH="$CONDA_PREFIX/linux"
# elif [[ "$UNAME_S" == "Darwin" ]]; then
#     JSON_FORTRAN_PREFIX="/opt/homebrew/Cellar/json-fortran/9.0.2"
#     INCLUDE_PATH=$JSON_FORTRAN_PREFIX/include
#     LIB_PATH=$JSON_FORTRAN_PREFIX/lib
#     BIN_PATH=$JSON_FORTRAN_PREFIX/bin
# else
#     echo "Unsupported OS: $UNAME_S"
#     exit 1
# fi
# BASE_PATH=`conda info --base`
# INCLUDE_PATH=`echo $BASE_PATH/envs/seidart/include`
# LIB_PATH=`echo $BASE_PATH/envs/seidart/lib`
INCLUDE_PATH="$CONDA_PREFIX/include"
LIB_PATH=`echo $CONDA_PREFIX/lib` 

echo "Include path is set to $INCLUDE_PATH"
echo "Library path is set to $LIB_PATH"
export LD_LIBRARY_PATH=/usr/local/jsonfortran-gnu-9.0.2/lib:$LD_LIBRARY_PATH

if [[ "$user_defined_builddir" = false ]]; then
    # BIN_PATH=`echo $BASE_PATH/envs/seidart/bin`
    BIN_PATH=`echo $CONDA_PREFIX/bin`
fi

echo "Install path is set to $BIN_PATH"

# -----------------------------------------------------------------------------
# Compiler and flags
FC=$(which gfortran)
echo "Using GFortran from $FC"

FFLAGS="-I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -g -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -O0 -fcheck=all -fPIC"

# OpenMP support (optional)
if [[ "$openmp_compile" = true ]]; then
    FFLAGS="$FFLAGS -fopenmp"
fi

# Source files
SRC_DIR="src/fortran"
SOURCES="$SRC_DIR/constants.f08 $SRC_DIR/seidart_types.f08 $SRC_DIR/seidartio.f08 $SRC_DIR/cpmlfdtd.f08 $SRC_DIR/main.f08"
EXECUTABLE="seidartfdtd-$os-$arch"

# Create executable directory if it doesn't exist
mkdir -p "$BIN_PATH"

# Compile the executable
echo "Compiling the SeidarT CPML FDTD executable..."
# $FC -v $FFLAGS -o $EXECUTABLE $SOURCES > compile_output.txt 2>&1
/usr/bin/gfortran -v -I/usr/local/jsonfortran-gnu-9.0.2/lib \
-L/usr/local/jsonfortran-gnu-9.0.2/lib -g -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow \
-O0 -fcheck=all -o seidartfdtd-linux-x86_64 \
src/fortran/constants.f08 src/fortran/seidart_types.f08 src/fortran/seidartio.f08 \
src/fortran/cpmlfdtd.f08 src/fortran/main.f08 \
-ljsonfortran -lgfortran -lm -shared-libgcc

if [[ $? -ne 0 ]]; then
    echo "Compilation failed. Check compile_output.txt for details."
    cat compile_output.txt
    exit 1
fi

# -----------------------------------------------------------------------------
# Move the executable to the correct folder
echo "Moving the executable to $BIN_PATH..."
mv $EXECUTABLE $BIN_PATH/

echo "Build and installation complete. Executable is located in $BIN_PATH."

if [[ "$clean" = true ]]; then
    echo "Cleaning up intermediate files..."
	rm -f $EXECUTABLE *.o *.mod
	echo "Cleanup complete."
fi