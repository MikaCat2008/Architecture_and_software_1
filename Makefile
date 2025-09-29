CC = gcc
CFLAGS = -g
OBJS = tests.o mat_int.o mat_float.o mat_double.o
TARGET = program

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
