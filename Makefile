# Compile the file called splay_tree.c

# Set variables

COMPILER = g++
CXXFLAGS = -Wall -O3 -g -std=c++20 -march=native -ffast-math # -fsanitize=address

all: splay_tree  enhanced_splay_tree enhanced_splay_tree_with_reversed

splay_tree: splay_tree.o
	${COMPILER} ${CXXFLAGS} splay_tree.o -o splay_tree


enhanced_splay_tree: enhanced_splay_tree.o
	${COMPILER} ${CXXFLAGS} enhanced_splay_tree.o -o enhanced_splay_tree


enhanced_splay_tree_with_reversed: enhanced_splay_tree_with_reversed.o
	${COMPILER} ${CXXFLAGS} enhanced_splay_tree_with_reversed.cpp -o enhanced_splay_tree_with_reversed


clean:
	rm -f splay_tree enhanced_splay_tree enhanced_splay_tree_with_reversed *.o

