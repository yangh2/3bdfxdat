CC=gcc
CFLAGS=-lm -g

SOURCES=3bdfxdat.c 3bdfxdat_sum.c 3bdf2s3.c 3bdf_view_theta.c 3bdfMsup.c 3bdf_print_sup.c pdf2supp.c hist2supp.c 3bdf2s3_eps.c

TARGETS=$(SOURCES:.c=)
HEADERS=hist.h 3bdf.h integrate.h supp.h random_16807.h
HSOURCES=$(HEADERS:.h=.c)

all : $(TARGETS)

.PHONY : all

$(TARGETS) : $(SOURCES) $(HEADERS) $(HSOURCES)
	$(CC)  -o $@ $(@).c  $(HSOURCES) $(CFLAGS)

clean:
	rm *.o
cleanall:
	rm  $(TARGETS)
