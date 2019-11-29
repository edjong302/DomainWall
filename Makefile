CC=g++
TARGET_DIR ?= ./bin
BUILD_DIR ?= ./build
SRC_DIR ?= ./src
INC_DIR := ./includes
OUTPUT_DIR := ./output

INC_FLAGS := $(addprefix -I,$(INC_DIR))
OPTFLAG = -O3 -ffast-math
OMPFLAG = -fopenmp

all: $(BUILD_DIR)/main.o $(BUILD_DIR)/domain.o $(BUILD_DIR)/potential.o $(BUILD_DIR)/initial_conditions.o
	$(CC) $(OPTFLAG) $(OMPFLAG) $(INC_FLAGS) -o domain $(BUILD_DIR)/main.o $(BUILD_DIR)/domain.o $(BUILD_DIR)/potential.o $(BUILD_DIR)/initial_conditions.o

$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp
	$(CC) $(OPTFLAG) $(OMPFLAG) $(INC_FLAGS) -c $(SRC_DIR)/main.cpp -o $(BUILD_DIR)/main.o

$(BUILD_DIR)/domain.o: $(SRC_DIR)/domain.cpp
	$(CC) $(OPTFLAG) $(OMPFLAG) $(INC_FLAGS) -c $(SRC_DIR)/domain.cpp -o $(BUILD_DIR)/domain.o

$(BUILD_DIR)/potential.o: $(SRC_DIR)/potential.cpp
	$(CC) $(OPTFLAG) $(OMPFLAG) $(INC_FLAGS) -c $(SRC_DIR)/potential.cpp -o $(BUILD_DIR)/potential.o

$(BUILD_DIR)/initial_conditions.o: $(SRC_DIR)/initial_conditions.cpp
	$(CC) $(OPTFLAG) $(OMPFLAG) $(INC_FLAGS) -c $(SRC_DIR)/initial_conditions.cpp -o $(BUILD_DIR)/initial_conditions.o

clean:
	rm $(BUILD_DIR)/*.o
	rm domain

data-clean:
	rm $(OUTPUT_DIR)/*.txt

gif:
	python domain.py