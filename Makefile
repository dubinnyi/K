
all:
	g++ -g ncs.cpp test_ncs.cpp -o test_ncs

test:
	./test_ncs
