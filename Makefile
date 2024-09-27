CC = gcc
CFLAGS = -lm -Wall -Wextra -Wno-unused-variable -Wno-error 

main : main.c functions.o
	$(CC) $(CFLAGS) -o main main.c functions.o -lm


functions.o: functions.c functions.h
	$(CC) $(CFLAGS) -c functions.c -lm
 
clean:
	rm -f functions.o main
