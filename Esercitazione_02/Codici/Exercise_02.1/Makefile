CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o Funzioni.o Integrazione.o TLC.o
	$(CC) random.o main.o Funzioni.o Integrazione.o TLC.o -o main.exe
main.o : main.cpp random.h Funzioni.h Integrazione.h TLC.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h 
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Funzioni.o : Funzioni.cpp Funzioni.h 
	$(CC) -c Funzioni.cpp -o Funzioni.o $(CFLAGS)
Integrazione.o : Integrazione.cpp Integrazione.h random.h Funzioni.h
	$(CC) -c Integrazione.cpp -o Integrazione.o $(CFLAGS)
TLC.o : TLC.cpp TLC.h 
	$(CC) -c TLC.cpp -o TLC.o $(CFLAGS)


clean :
	rm *.o main.exe seed.out 
