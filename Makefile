test: test.o alg.o modring.o numbers.o complex.o numberfield.o polyio.o
	g++ -o test test.o alg.o modring.o numbers.o complex.o numberfield.o polyio.o -lgmp -lgmpxx -g -pg

test.o: test.cpp polyring.h modring.h polymodring.h alg.h numberfield.h polyio.h
	g++ -c test.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg
	
alg.o: alg.cpp alg.h polyring.h modring.h polymodring.h typedefs.h numbers.h numberfield.h
	g++ -c alg.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg
	
modring.o: modring.cpp modring.h numbers.h
	g++ -c modring.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg
	
numbers.o: numbers.cpp numbers.h
	g++ -c numbers.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg
	
complex.o: complex.cpp complex.h numbers.h
	g++ -c complex.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg
	
numberfield.o: numberfield.cpp numberfield.h polyring.h modring.h typedefs.h numbers.h
	g++ -c numberfield.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg
	
polyio.o: polyio.cpp polyring.h modring.h polymodring.h alg.h numberfield.h polyio.h
	g++ -c polyio.cpp -std=c++11 -g -isystem /usr/include/eigen3/ -pg

## Remove all the compilation and debugging files
clean:
	rm -f core test *.o *~
