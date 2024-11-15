# 'make USE_OPENMP=1' will compile with the openmp directives. Default 'make' 
# will compile without OpenMP. 
#

# The executable will by default be placed in the BIN_PATH. If a different path 
# is passed at the command line, we want to use that one

BIN_PATH ?= 
UNAME_S := $(shell uname -s)

# Compiler and flags
FC := gfortran 
# FCFLAGS = -O3
# LDFLAGS = -ljsonfortran

ifeq ($(UNAME_S), Linux)
    CONDA_PREFIX := $(shell conda info --base)/envs/seidart
    INCLUDE_PATH := $(CONDA_PREFIX)/include
    LIB_PATH := $(CONDA_PREFIX)/lib
    BIN_PATH ?= $(CONDA_PREFIX)/bin
else ifeq ($(UNAME_S), Darwin)
    JSON_FORTRAN_PREFIX := /opt/homebrew/Cellar/json-fortran/9.0.2
    INCLUDE_PATH := $(JSON_FORTRAN_PREFIX)/include
    LIB_PATH := $(JSON_FORTRAN_PREFIX)/lib
    BIN_PATH ?= $(JSON_FORTRAN_PREFIX)/bin
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
	@echo "To uninstall run 'make uninstall'"

# Rule to build the executable
$(EXECUTABLE): $(SOURCES)
	@echo "Compiling the SeidarT CPML FDTD executable..."
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(SOURCES) -lgfortran
	@echo "Compilation finished."

# Install the executable into the correct Miniconda folder
install: $(EXECUTABLE)
	@echo "Installing the executable to $(BIN_PATH)..."
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
