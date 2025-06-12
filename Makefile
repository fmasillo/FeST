# Set variables
COMPILER = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -O3 -g -std=c++20 -march=native -ffast-math # -fsanitize=address

all: enhanced_splay_tree 

enhanced_splay_tree: FeST.o
	${COMPILER} ${CXXFLAGS} FeST.o -o FeST



clean:
	rm -f FeST *.o

