main.o: explicit.o utils.o
	g++ main.cpp calculator.o utils.o -o main

explicit.o: utils.o calculator.cpp calculator.h
	g++ -c -fPIC calculator.cpp -o calculator.o

utils.o: utils.cpp utils.h
	g++ -c -fPIC utils.cpp -o utils.o

clean:
	rm calculator.o calculator.so utils.o
