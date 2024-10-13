CONDA_PREFIX := $(shell conda info --base)/envs/seidart
INCLUDE_PATH := $(CONDA_PREFIX)/include
LIB_PATH := $(CONDA_PREFIX)/lib
BIN_PATH := $(CONDA_PREFIX)/bin

# Compiler and flags
FC := gfortran
FFLAGS := -I${INCLUDE_PATH} -L${LIB_PATH} -ljsonfortran -g -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -O0 -fcheck=all

# Source files
SRC_DIR := src/fortran
SOURCES := $(SRC_DIR)/constants.f08 $(SRC_DIR)/seidart_types.f08 $(SRC_DIR)/seidartio.f08 $(SRC_DIR)/cpmlfdtd.f08 $(SRC_DIR)/main.f08
EXECUTABLE := seidartfdtd 

# Default target
all: $(EXECUTABLE) install clean
	@echo "To uninstall run 'make uninstall'"

# Rule to build the executable
$(EXECUTABLE): $(SOURCES)
	@echo "Compiling the SeidarT CPML FDTD executable..."
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(SOURCES) -lgfortran > compile_output.txt 2>&1
	@echo "Compilation finished."

# Install the executable into the correct Miniconda folder
install: $(EXECUTABLE)
	@echo "Installing the executable to $(BIN_PATH)..."
	cp $(EXECUTABLE) $(BIN_PATH)
	@echo "SeidarT CPML FDTD executable is in the $(BIN_PATH)"
	
# Uninstall
uninstall:
	@echo "Removing $(BIN_PATH)$(EXECUTABLE)"
	rm -f $(BIN_PATH)/$(EXECUTABLE)
	@echo "Uninstall successful."
	

# Clean up
clean:
	@echo "Cleaning up intermediate files..."
	rm -f $(EXECUTABLE) *.o *.mod
	@echo "Cleanup complete."

# Uninstall and clean
clean_all: uninstall clean