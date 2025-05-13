# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++23

# Target executable
TARGET = 3bp

# Source files
SRCS = 3bp.cpp fileutils.cpp stringutils.cpp gnuplot.cpp

# Header files
HDRS = fileutils.h stringutils.h gnuplot.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile source files into object files
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean