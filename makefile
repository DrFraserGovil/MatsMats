## This is a custom makefile built with the help of an LLM because I don't know how makefiles work and have been using the same one for 10 years.
## Fingers crossed.

# Project Name
PROJECT = matsmats

# Compiler
CC = g++

# Compiler and Linker Options
CXXFLAGS = -std=c++20 -pthread -O3 -Wall
LDFLAGS = -lpthread

# Dependency Flags
DEPFLAGS = -MMD -MP

# Source and Build Directories
SRC_DIR = src
BUILD_DIR = .build

### DON'T EDIT BELOW HERE

# Find all source files for the main application
# This assumes all .cpp files in SRC_DIR and its subdirectories are part of the main application.
SOURCE_FILES := $(shell find $(SRC_DIR) -name "*.cpp")

# Generate object files for the main application
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCE_FILES))

# Dependency Files for Main
MAIN_DEPS := $(OBJECTS:.o=.d)

# Default Target: Build the main project
.PHONY: all
all: $(PROJECT)

# Build Main Project executable
$(PROJECT): $(OBJECTS)
	@echo "Linking $(PROJECT)..."
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile Main Source Files
# This rule creates directories as needed for object files, maintaining the source directory structure in .build/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo "Compiling $< for $(PROJECT)..."
	$(CC) $(CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

# Include Dependencies
-include $(MAIN_DEPS)

# Run Main Project
.PHONY: run
run: $(PROJECT)
	./$(PROJECT)

# Clean Targets
.PHONY: clean
clean:
	@echo "Cleaning up build..."
	rm -rf $(PROJECT) $(BUILD_DIR)

.PHONY: depclean
depclean:
	@echo "Removing dependency files..."
	rm -f $(MAIN_DEPS)

.PHONY: clean-all # Kept for convenience, effectively same as 'clean depclean' now
clean-all: clean depclean