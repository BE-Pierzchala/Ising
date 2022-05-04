main: main.cpp myFunctions.h
	clang++ -std=c++11 main.cpp -o base
	clang++ -std=c++11 main_h.cpp -o base_h

