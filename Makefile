CC = g++
OBJS = main.o
CFLAGS = -g -fexec-charset=cp1251 -finput-charset=cp1251
SRC_DIR = src
BIN_DIR = bin
TARGET = $(BIN_DIR)/main

all: $(TARGET)

$(TARGET): $(addprefix $(SRC_DIR)/,$(OBJS)) | $(BIN_DIR)
	$(CC) $(addprefix $(SRC_DIR)/,$(OBJS)) -o $(TARGET) & del /s "src\*.o"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
