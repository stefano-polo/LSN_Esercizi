CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o Salesman.o
	$(CC) random.o main.o Salesman.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Salesman.o : Salesman.cpp Salesman.h
	$(CC) -c Salesman.cpp -o Salesman.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.out


esegui:
	make && ./main.exe 
