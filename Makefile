CC = g++
CXXFLAGS = -Wall -g
TARGET = program
SRCS = classes.cpp
OBJS = $(SRCS:.c=.o)

all: $(TARGET)
	@echo "Running test:"
	./$(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
