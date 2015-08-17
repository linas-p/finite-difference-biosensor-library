explicit.o: utils.o explicit_calculator.cpp explicit_calculator.h 
	g++ -c -fPIC explicit_calculator.cpp -o explicit_calculator.o

utils.o: utils.cpp utils.h
	g++ -c -fPIC utils.cpp -o utils.o

clean:
	rm explicit_calculator.o explicit_calculator.so utils.o
