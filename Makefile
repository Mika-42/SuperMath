CXX = g++
FLAGS = -Wall -Wextra
CPP_VERSION = -std=c++23
CPP_FILE = src/*.cpp
OUTPUT_DIR = bin
OUTPUT_NAME = out
DEBUG = -O1 -g3 -time

all:
	${CXX} ${CPP_VERSION} ${CPP_FILE} -o ${OUTPUT_DIR}/${OUTPUT_NAME} ${DEBUG} ${FLAGS}
	${OUTPUT_DIR}/${OUTPUT_NAME}
	