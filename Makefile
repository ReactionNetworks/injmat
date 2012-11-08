#####################################################
#
#  Makefile for the whole system:
#  Note that certain features are for readability
#  and not for compactness of code
#  
######################################################
CC = gcc -Wall
CXX = g++ -Wall
CFLAGS = -g -O2
CXXFLAGS= -g -O2
VERSION = 1.0
MODELNAME = CRNLIB$(VERSION)

HEADERS =	src/crn.h

MAINLIB = 	src/lib/libcrn.a

MAINLINK =	-lcrn -lcln -lginac

OBJDIR = 	src/obj

AR      =       ar crs

OBJECTS =	src/obj/crn.o


project: bin/test;

clean: 
	rm bin/* src/lib/* src/obj/*

src/obj/%.o: src/%.c $(HEADERS)
	$(CXX) -c -o $@ $<

src/lib/libcrn.a: $(HEADERS) $(OBJECTS)
	$(AR) src/lib/libcrn.a $(OBJECTS)



##################

bin/test: src/test.c $(MAINLIB)
	$(CXX) -o bin/test src/test.c -lm -Lsrc/lib $(MAINLINK)


