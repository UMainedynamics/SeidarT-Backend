#!/bin/bash

# -----------------------------------------------------------------------------
user_defined_builddir=false
openmp_compile=false
clean=false
debug=false
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
    -d|--debug)
      debug=true
      shift
      ;;
    -*|--*)
      echo "Unknown option: $1"
      echo "Usage: cmd [-b install_directory] [-c|--clean] [-d|--debug]"
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

INCLUDE_PATH="$CONDA_PREFIX/include"
LIB_PATH=`echo $CONDA_PREFIX/lib` 

echo "Include path is set to $INCLUDE_PATH"
echo "Library path is set to $LIB_PATH"
# export LD_LIBRARY_PATH=/usr/local/jsonfortran-gnu-9.0.2/lib:$LD_LIBRARY_PATH

echo "Listed files and directories in the include path"
echo "Listed files and directories in the library path"

if [[ "$user_defined_builddir" = false ]]; then
    # BIN_PATH=`echo $BASE_PATH/envs/seidart/bin`
    BIN_PATH=`echo $CONDA_PREFIX/bin`
fi

echo "Install path is set to $BIN_PATH"

# Check if the OS is macOS
if [[ "$(uname)" == "Darwin" ]]; then
    # Ensure that CONDA_PREFIX is set (it should be if your environment is activated)
    if [ -n "$CONDA_PREFIX" ]; then
        export DYLD_LIBRARY_PATH="$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH"
        echo "DYLD_LIBRARY_PATH set to: $DYLD_LIBRARY_PATH"
    else
        echo "Warning: CONDA_PREFIX is not set. Are you sure the conda environment is activated?"
    fi
fi

# -----------------------------------------------------------------------------
# Compiler and flags
FC=$(which gfortran)
echo "Using GFortran from $FC"

# Debugging options are the first FFLAGS. Uncomment when necessary
# FFLAGS=" -I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -march=native -funroll-loops -g -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -fPIC -O2 -ffast-math -ffp-contract=fast"
# These are the production flags
# if [[ "$debug" = true ]]; then 
#   FFLAGS="-I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -std=f2008 -g -Og \
#                                           -fbacktrace -fcheck=all \
#                                           -ffpe-trap=invalid,zero,overflow \
#                                           -finit-real=snan -finit-integer=0 \
#                                           -Wall"
# else
FFLAGS="-I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -O3 -march=native -funroll-loops -ffast-math -ffp-contract=fast -fomit-frame-pointer -fPIC -Wall"
# fi 

# OpenMP support (optional)
if [[ "$openmp_compile" = true ]]; then
    FFLAGS="$FFLAGS -fopenmp"
fi

# Source files
SRC_DIR="src/fortran"
SOURCES="$SRC_DIR/constants.f08 \
         $SRC_DIR/seidart_types.f08 \
         $SRC_DIR/seidartio.f08 \
         $SRC_DIR/density_averaging.f08 \
         $SRC_DIR/cpmlfdtd.f08 \
         $SRC_DIR/main.f08"
EXECUTABLE="seidartfdtd-$os-$arch"

# Create executable directory if it doesn't exist
mkdir -p "$BIN_PATH"

# Compile the executable
echo "Compiling the SeidarT CPML FDTD executable..."
$FC -v $FFLAGS -o $EXECUTABLE $SOURCES > compile_output.txt 2>&1

if [[ $? -ne 0 ]]; then
    echo "Compilation failed. Check compile_output.txt for details."
    cat compile_output.txt
    exit 1
fi

# If debug mode, exit after successful build
if [[ "$debug" = true ]]; then
    echo "Debug build complete; skipping installation."
    exit 0
fi

# -----------------------------------------------------------------------------
# Move the executable to the correct folder
echo "Moving the executable to $BIN_PATH..."
mv $EXECUTABLE $BIN_PATH/seidartfdtd

echo "Build and installation complete. Executable is located in $BIN_PATH."

if [[ "$clean" = true ]]; then
    echo "Cleaning up intermediate files..."
	rm -f $EXECUTABLE *.o *.mod
	echo "Cleanup complete."
fi