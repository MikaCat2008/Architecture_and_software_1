CC = gcc
OBJS = tests.o mat_int.o mat_float.o mat_double.o
CFLAGS = -g
SRC_DIR = src
BIN_DIR = bin
TARGET = $(BIN_DIR)/program

all: $(TARGET)

$(TARGET): $(addprefix $(SRC_DIR)/,$(OBJS)) | $(BIN_DIR)
	$(CC) $(addprefix $(SRC_DIR)/,$(OBJS)) -o $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
