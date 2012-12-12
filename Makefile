all:main.cpp evaluate.cpp
	g++ -g main.cpp evaluate.cpp -o main
clean:
	rm -rf *.o
