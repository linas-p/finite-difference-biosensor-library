main.o: explicit.o utils.o
	g++ -pedantic main.cpp calculator.o utils.o -o main

explicit.o: utils.o calculator.cpp calculator.h
	g++ -pedantic -c -fPIC calculator.cpp -o calculator.o

utils.o: utils.cpp utils.h
	g++ -pedantic -c -fPIC utils.cpp -o utils.o

clean:
	rm calculator.o utils.o
