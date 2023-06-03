all:
	gcc -o computeBondLength computeBondLength.c -lm -Wall
	gcc -o msd2 msd2.c -lm -Wall
	gcc -o msd msd.c -lm -Wall
