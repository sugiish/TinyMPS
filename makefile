# Build all source files.
#
# Targets in this file are:
# all     Compiles and links all source files
# clean   Removes all intermediate files
#
# Copyright (c) 2017 Shota SUGIHARA
# Distributed under the MIT License.

TARGET := mps.out
SOURCE_DIR := src
OBJECT_DIR := obj
INCLUDE_DIR := include

CC := g++
CXX := g++
DEBUGFLAGS := -g -O0
CXXFLAGS := $(DEBUGFLAGS) -std=c++11 -Wall -Wextra -MP -MMD
CPPFLAGS := -I $(INCLUDE_DIR)

MKDIR := mkdir -p
MV := mv -f
RM := rm -f
SED := sed
TEST := test

# Creates an object directory if it does not exist.
to_create_dir := $(OBJECT_DIR)
create-object-directory :=                  \
	$(shell for f in $(to_create_dir);      \
		do                                  \
			$(TEST) -d $$f | $(MKDIR) $$f;  \
		done)

sources := $(wildcard $(SOURCE_DIR)/*.cpp)
objects := $(addprefix $(OBJECT_DIR)/, $(notdir $(sources:.cpp=.o)))
dependencies := $(objects:.o=.d)

.PHONY: all
all: $(TARGET)

$(TARGET): $(objects)
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

.PHONY: clean
clean:
	$(RM) $(TARGET) $(objects) $(dependencies)

ifneq "$(MAKECMDGOALS)" "clean"
	-include $(dependencies)
endif