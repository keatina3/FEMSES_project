CC = g++

CFLAGS = -g -Wall -W -O3

LDFLAGS =

OBJECTS = *.o
FILES = *.cc

TARGET = prog

all: $(FILES)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

#.PHONY: test mem_test clean
#
#test: all
#	./$(TARGET) -m 95000 -t 350 -s 17 1> results.out 2> error.log
#
#mem_test: all
#	valgrind ./$(TARGET) -m 900 -t 350 -s 17

clean:
	$(RM) $(OBJECTS) $(TARGET) results.out error.log
