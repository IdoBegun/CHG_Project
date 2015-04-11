#!/bin/tcsh

rm -f main.o LDTree.o exe
g++ -g -Wall -c LDTree.cpp
g++ -g -Wall -c main.cpp
g++ -g -Wall -o exe main.o LDTree.o

