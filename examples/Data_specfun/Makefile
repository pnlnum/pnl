LDFLAGS=-L/opt/local/lib  -lgsl
CPPFLAGS=-I/opt/local/include
CFLAGS=-Wall
CC=gcc

all:expint_gsl hyperg_gsl

data:
	./expint_gsl

expint_gsl: expint_gsl.o
	$(CC) $(LDFLAGS) -o $@ $<


hyperg_gsl: hyperg_gsl.o
	$(CC) $(LDFLAGS) -o $@ $<
