CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o TLC.o Funzioni.o Metropolis.o
	$(CC) main.o random.o TLC.o Funzioni.o Metropolis.o -o main.exe
main.o : main.cpp Funzioni.h random.h Metropolis.h TLC.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h 
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Metropolis.o : Metropolis.cpp random.h Funzioni.h
	$(CC) -c Metropolis.cpp -o Metropolis.o $(CFLAGS)
TLC.o : TLC.cpp TLC.h 
	$(CC) -c TLC.cpp -o TLC.o $(CFLAGS)
Funzioni.o : Funzioni.cpp Funzioni.h 
	$(CC) -c Funzioni.cpp -o Funzioni.o $(CFLAGS)
clean:
	\rm main.exe *.o seed.out
esegui:
	./main.exe
