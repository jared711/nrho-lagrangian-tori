# Makefile for nrho-lagrangian-tori
# Computes invariant tori around NRHOs using parameterization method

# Directories
SRC_DIR = src/cpp
HEADER_DIR = src/cpp/headers
BUILD_DIR = build
BIN_DIR = bin

# Compiler flags
OPT = -g -Wall
#OPT = -O3 -Wall
# Use -iquote for local headers (searched for #include "...")
# but not -I (which affects #include <...>)
CFLAGS = $(OPT) -fopenmp -ffast-math -fdiagnostics-color=always -iquote$(HEADER_DIR)
CXXFLAGS = $(OPT) -fopenmp -ffast-math -fdiagnostics-color=always -I$(HEADER_DIR)

# Target executable
TARGET = $(BIN_DIR)/param

# Source files
CXX_SOURCES = $(SRC_DIR)/param.cc
C_SOURCES = $(SRC_DIR)/seccp.c \
            $(SRC_DIR)/fluxvp.c \
            $(SRC_DIR)/rk78vp.c \
            $(SRC_DIR)/rtbphp.c \
            $(SRC_DIR)/campvp.c \
            $(SRC_DIR)/scread.c \
            $(SRC_DIR)/vbprintf.c

# Object files (in build directory)
OBJECTS = $(BUILD_DIR)/seccp.o \
          $(BUILD_DIR)/fluxvp.o \
          $(BUILD_DIR)/rk78vp.o \
          $(BUILD_DIR)/rtbphp.o \
          $(BUILD_DIR)/campvp.o \
          $(BUILD_DIR)/scread.o \
          $(BUILD_DIR)/vbprintf.o

# Header files (for dependency tracking)
HEADERS = $(HEADER_DIR)/complex.h \
          $(HEADER_DIR)/grid.h \
          $(HEADER_DIR)/matrix.h \
          $(HEADER_DIR)/utils.h \
          $(HEADER_DIR)/seccp.h \
          $(HEADER_DIR)/fluxvp.h \
          $(HEADER_DIR)/rk78vp.h \
          $(HEADER_DIR)/rtbphp.h \
          $(HEADER_DIR)/campvp.h \
          $(HEADER_DIR)/scread.h \
          $(HEADER_DIR)/vbprintf.h

.PHONY: all clean realclean help

all: $(TARGET)

# Link executable
$(TARGET): $(CXX_SOURCES) $(OBJECTS) | $(BIN_DIR)
	g++ -o $@ $(CXXFLAGS) $(CXX_SOURCES) $(OBJECTS) -lm
	@echo ""
	@echo "Build complete: $(TARGET)"
	@echo "Run with: ./$(TARGET) data/initial_conditions/approxQPO.csv"

# Compile C source files to object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	gcc -c $(CFLAGS) $< -o $@

# Create directories if they don't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean build artifacts
clean:
	rm -f $(BUILD_DIR)/*.o
	@echo "Build artifacts cleaned"

# Clean everything including executable
realclean: clean
	rm -f $(TARGET)
	@echo "Executable removed"

# Help target
help:
	@echo "Makefile for nrho-lagrangian-tori"
	@echo ""
	@echo "Targets:"
	@echo "  all (default) - Build the param executable"
	@echo "  clean         - Remove object files from build/"
	@echo "  realclean     - Remove object files and executable"
	@echo "  help          - Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make          # Build param"
	@echo "  make clean    # Clean build artifacts"
	@echo "  ./bin/param data/initial_conditions/approxQPO.csv"
	@echo ""
	@echo "Directory structure:"
	@echo "  src/cpp/        - C/C++ source files"
	@echo "  src/cpp/headers/ - Header files"
	@echo "  build/          - Object files (.o)"
	@echo "  bin/            - Executables"
	@echo "  data/           - Input data files"
	@echo "  output/         - Generated torus output"
	@echo "  results/        - Test results and logs"
