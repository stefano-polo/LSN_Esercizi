CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o Buffon.o TLC.o
	$(CC) random.o main.o Buffon.o TLC.o -o main.exe
main.o : main.cpp Buffon.h random.h TLC.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

Buffon.o : Buffon.cpp Buffon.h random.h
	$(CC) -c Buffon.cpp -o Buffon.o $(CFLAGS)
TLC.o : TLC.cpp TLC.h 
	$(CC) -c TLC.cpp -o TLC.o $(CFLAGS)


clean :
	rm *.o main.exe seed.out 
