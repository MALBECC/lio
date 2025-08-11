# PyLIO - Python Wrapper for LIO
# Makefile for building Python extensions

.PHONY: all clean test install wrapper libs

# Default target
all: libs wrapper

# Build LIO libraries first
libs:
	@echo "Building LIO libraries..."
	$(MAKE) -C .. g2g
	$(MAKE) -C .. liblio
	@echo "✓ LIO libraries built"

# Create Python wrapper
wrapper: libs
	@echo "Setting up Python wrapper..."
	python3 -c "import numpy; print('✓ NumPy available')"
	python3 -c "import ctypes; print('✓ ctypes available')"
	chmod +x *.py
	@echo "✓ Python wrapper ready"

# Test the wrapper
test: wrapper
	@echo "Testing PyLIO wrapper..."
	python3 simple_pylio.py

# Install Python dependencies
install-deps:
	@echo "Installing Python dependencies..."
	pip3 install numpy

# Clean build artifacts
clean:
	rm -f *.pyc *.pyo
	rm -rf __pycache__
	rm -f *.inp *.out *.bas
	rm -f lio_fortran.*.so
	rm -f f2py_lio_interface.f90

# Advanced: build f2py extension (requires proper Fortran interface)
f2py-build:
	@echo "Building f2py extension..."
	python3 f2py_wrapper.py

# Create example files
examples:
	@echo "Creating example calculations..."
	mkdir -p examples
	python3 -c "
import numpy as np
from simple_pylio import LIOInputGenerator

# Water molecule
gen = LIOInputGenerator()
atoms = ['O', 'H', 'H']
coords = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])
gen.create_input_file('examples/water.inp', atoms, coords, basis='sto-3g', functional='pbe')

# Methane molecule  
atoms = ['C', 'H', 'H', 'H', 'H']
coords = np.array([[0.0, 0.0, 0.0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629], [-0.629, 0.629, -0.629], [0.629, -0.629, -0.629]])
gen.create_input_file('examples/methane.inp', atoms, coords, basis='sto-3g', functional='b3lyp')

print('✓ Example files created in examples/')
"

# Show usage information
help:
	@echo "PyLIO Makefile"
	@echo "=============="
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build libraries and wrapper (default)"
	@echo "  libs         - Build LIO libraries only"
	@echo "  wrapper      - Setup Python wrapper"
	@echo "  test         - Test the wrapper"
	@echo "  install-deps - Install Python dependencies"
	@echo "  examples     - Create example input files"
	@echo "  f2py-build   - Build f2py extension (experimental)"
	@echo "  clean        - Clean build artifacts"
	@echo "  help         - Show this help"
	@echo ""
	@echo "Usage:"
	@echo "  make all       # Build everything"
	@echo "  make test      # Test the wrapper"
	@echo "  make examples  # Create examples"
