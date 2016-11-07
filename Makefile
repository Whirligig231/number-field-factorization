test: test.o alg.o modring.o numbers.o
	g++ -o test test.o alg.o modring.o numbers.o -lgmp -lgmpxx -g

test.o: test.cpp polyring.h modring.h polymodring.h alg.h
	g++ -c test.cpp -std=c++11 -g -isystem /usr/include/eigen3/
	
alg.o: alg.cpp alg.h polyring.h modring.h typedefs.h numbers.h
	g++ -c alg.cpp -std=c++11 -g -isystem /usr/include/eigen3/
	
modring.o: modring.cpp modring.h numbers.h
	g++ -c modring.cpp -std=c++11 -g -isystem /usr/include/eigen3/
	
numbers.o: numbers.cpp numbers.h
	g++ -c numbers.cpp -std=c++11 -g -isystem /usr/include/eigen3/

## Remove all the compilation and debugging files
clean:
	rm -f core test *.o *~
