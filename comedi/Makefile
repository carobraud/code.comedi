noinst_PROGRAMS = cmd mmap sample20k6ch

LIBS = -lcomedi -lm
LIBCImg = -I.. -I../NetCDF/include -lnetcdf_c++ -lnetcdf -Wall -W -lm -lcomedi -ansi -pedantic -O0 -check -stacktrace -fno-tree-pre -Dcimg_use_vt100 -I/usr/X11R6/include -Dcimg_use_xshm -Dcimg_use_xrandr -L/usr/X11R6/lib -lpthread -lX11 -lXext -lXrandr
EXTRA = common.c
CC = gcc
CPP = g++

all: $(noinst_PROGRAMS)

cmd: cmd.cc
	$(CPP) $(EXTRA) $@.cc $(LIBS) -o $@
mmap: mmap.cpp common.c examples.h
	$(CPP) $(EXTRA) $@.cpp $(LIBCImg) -o $@
sample20k6ch: sample20k6ch.cc common.c examples.h
	$(CPP) $(EXTRA) $@.cc $(LIBCImg) -o $@

clean:
	rm -f *.o
	rm -rf .libs
	@list='$(noinst_PROGRAMS)'; for p in $$list; do \
	  rm -f $$p ; \
	done


