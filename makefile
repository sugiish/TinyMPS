CC := g++
CXX := g++
TARGET := mps.out
CXXFLAGS := -g -O0 -Wall -MP -MMD

#COMPILE.C = $(COMPILE.cc)
#COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
#OUTPUT_OPTION = -o $@

SOURCE_DIR := src
OBJECT_DIR := obj
INCLUDE_DIR := include

#VPATH = $(SOURCE_DIR) $(INCLUDE_DIR)
#vpath %.h $(INCLUDE_DIR)
#vpath %.cpp $(SOURCE_DIR)
#vpath %.o $(OBJECT_DIR)

CPPFLAGS := -I $(INCLUDE_DIR)

MKDIR := mkdir -p
MV := mv -f
RM := rm -f
SED := sed
TEST := test

sources := $(wildcard $(SOURCE_DIR)/*.cpp)
objects := $(addprefix $(OBJECT_DIR)/, $(notdir $(sources:.cpp=.o)))
dependencies := $(objects:.cpp=.d)

create-object-directory :=					\
	$(shell for f in $(OBJECT_DIR);			\
		do									\
			$(TEST) -d $$f | $(MKDIR) $$f;	\
		done)


#mps.out: main.o particles.o timer.o
#	g++ -Wall $(CPPFLAGS) $^ -o $@

#main.o: main.cpp
#particles.o: particles.cpp
#particles.o: particles.h
#timer.o: timer.cpp
#grid.o: grid.cpp

#$(TARGET): $(objects)
#	$(LINK.o) $^ $(LOADLIBES) $(LDLIBS) -o $@


#%.out: %.o
#	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

#all:$(TARGET)


#$(TARGET): $(objects)
#	$(LINK.o) $^ $(LOADLIBES) $(LDLIBS) -o $@

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

.PHONY: clean
clean:
	$(RM) $(TARGET) $(objects) $(dependencies)

-include $(dependencies)