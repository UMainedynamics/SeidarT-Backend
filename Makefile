# 'make USE_OPENMP=1' will compile with the openmp directives. Default 'make' 
# will compile without OpenMP. 
#
UNAME_S := $(shell uname -s)


CONDA_PREFIX := $(shell conda info --base)/envs/seidart
INCLUDE_PATH := $(CONDA_PREFIX)/include
LIB_PATH := $(CONDA_PREFIX)/lib
BIN_PATH := $(CONDA_PREFIX)/bin

# Compiler and flags
FC := $(shell which gfortran)
$(info Using GFortran from $(FC))

# ifeq ($(UNAME_S), Darwin)
#     UNSET_CONDA_ENV_VARS := $(shell unset CONDA_PREFIX && unset LD_LIBRARY_PATH && unset DYLD_LIBRARY_PATH)
# 	JSON_FORTRAN_PREFIX := $(shell brew --prefix json-fortran)
# 	LIB_PATH := $(JSON_FORTRAN_PREFIX)/lib
# 	INCLUDE_PATH := $(JSON_FORTRAN_PREFIX)/include/jsonfortran
	
# endif

ifeq ($(UNAME_S), Linux)
    CONDA_PREFIX := $(shell conda info --base)/envs/seidart
    INCLUDE_PATH := $(CONDA_PREFIX)/include
    LIB_PATH := $(CONDA_PREFIX)/lib
    BIN_PATH := $(CONDA_PREFIX)/bin
else ifeq ($(UNAME_S), Darwin)
    JSON_FORTRAN_PREFIX := $(shell brew --prefix json-fortran) 
	#/opt/homebrew/Cellar/json-fortran/9.0.2/ #$(shell brew --prefix json-fortran)
    INCLUDE_PATH := $(JSON_FORTRAN_PREFIX)/include
    LIB_PATH := $(JSON_FORTRAN_PREFIX)/lib
    BIN_PATH := $(JSON_FORTRAN_PREFIX)/bin
endif

FFLAGS := -I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -g -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -O0 -fcheck=all

# Source files
SRC_DIR := src/fortran
SOURCES := $(SRC_DIR)/constants.f08 $(SRC_DIR)/seidart_types.f08 $(SRC_DIR)/seidartio.f08 $(SRC_DIR)/cpmlfdtd.f08 $(SRC_DIR)/main.f08
EXECUTABLE := seidartfdtd 

# OpenMP support (conditionally add -fopenmp flag)
USE_OPENMP ?= 0
ifeq ($(USE_OPENMP), 1)
    FFLAGS += -fopenmp
endif

# Default target
all: $(EXECUTABLE) install clean
	# export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/json-fortran/9.0.2/lib:$$DYLD_LIBRARY_PATH
	@echo "To uninstall run 'make uninstall'"

# Rule to build the executable
$(EXECUTABLE): $(SOURCES)
	@echo "Compiling the SeidarT CPML FDTD executable..."
	# export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/json-fortran/9.0.2/lib:$$DYLD_LIBRARY_PATH
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(SOURCES) -lgfortran #> compile_output.txt 2>&1
	@echo "Compilation finished."

# Install the executable into the correct Miniconda folder
install: $(EXECUTABLE)
	@echo "Installing the executable to $(BIN_PATH)..."
	# export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/json-fortran/9.0.2/lib:$$DYLD_LIBRARY_PATH
	mv $(EXECUTABLE) $(BIN_PATH)
	@echo "SeidarT CPML FDTD executable is in the $(BIN_PATH)"

	
# Uninstall
uninstall:
	@echo "Removing $(BIN_PATH)$(EXECUTABLE)"
	rm -f $(BIN_PATH)/$(EXECUTABLE)
	@echo "Uninstall successful."

# Reinstall target depends on uninstall, install, and clean
reinstall: uninstall install

# Clean up
clean:
	@echo "Cleaning up intermediate files..."
	rm -f $(EXECUTABLE) *.o *.mod
	@echo "Cleanup complete."

# Uninstall and clean
clean_all: uninstall clean