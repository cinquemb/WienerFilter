all:
	clang++ -std=gnu++11 wiener.cpp -o wiener

debug:
	clang++ -g -std=gnu++11 wiener.cpp -o wiener