CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o emcee.o TLC.o 
	$(CC) random.o main.o emcee.o TLC.o  -o main.exe
main.o : main.cpp random.h emcee.h TLC.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h 
	$(CC) -c random.cpp -o random.o $(CFLAGS)
emcee.o : emcee.cpp emcee.h random.h
	$(CC) -c emcee.cpp -o emcee.o $(CFLAGS)
TLC.o : TLC.cpp TLC.h 
	$(CC) -c TLC.cpp -o TLC.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

esegui:
	./main.exe
