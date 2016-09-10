test: test.o
	g++ -o test test.o -lgmp -lgmpxx -g

test.o: test.cpp polyring.h modring.h
	g++ -c test.cpp -std=c++11 -g

## Remove all the compilation and debugging files
clean: 
	rm -f core test main.o file-processing.o threads.o compute.o *~
