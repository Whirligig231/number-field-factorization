test: test.o alg.o modring.o
	g++ -o test test.o alg.o modring.o -lgmp -lgmpxx -g

test.o: test.cpp polyring.h modring.h alg.h
	g++ -c test.cpp -std=c++11 -g
	
alg.o: alg.cpp alg.h polyring.h modring.h typedefs.h
	g++ -c alg.cpp -std=c++11 -g
	
modring.o: modring.cpp modring.h
	g++ -c modring.cpp -std=c++11 -g

## Remove all the compilation and debugging files
clean:
	rm -f core test *.o *~
