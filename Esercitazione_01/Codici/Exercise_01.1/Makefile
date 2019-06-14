CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o TLC.o
	$(CC) random.o main.o TLC.o -o main.exe
main.o : main.cpp TLC.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

TLC.o : TLC.cpp TLC.h random.h
	$(CC) -c TLC.cpp -o TLC.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out 
