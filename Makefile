all:
	gcc -o computeBondLength computeBondLength.c -lm -Wall
	./computeBondLength bond_info.dump
