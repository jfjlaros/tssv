SWIG := swig

CC := gcc
CFLAGS := -O3 -fPIC -I/usr/include/python2.7
DFLAGS := -Wall -ggdb

LD  := ld
LFLAGS := -shared

INPUT_FILES := sg_align.c
PRECIOUS := .py
TEMPORARY := _wrap.c .o _wrap.o

INPUT := $(basename $(INPUT_FILES))
TEMPORARY_FILES := $(foreach I, $(TEMPORARY), $(addsuffix $I, $(INPUT)))
PRECIOUS_FILES := $(foreach I, $(PRECIOUS), $(addsuffix $I, $(INPUT)))
OUTPUT_FILES := _$(INPUT).so

.SUFFIXES:
.PHONY: all clean distclean
.PRECIOUS: $(PRECIOUS_FILES)

all: $(OUTPUT_FILES)

clean:
	rm -f $(TEMPORARY_FILES)

distclean: clean
	rm -f $(PRECIOUS_FILES) $(OUTPUT_FILES)

%_wrap.c %.py: %.i
	$(SWIG) -python $<
%.py: %.i
	@

%.o: %.c %_wrap.c
	$(CC) $(CFLAGS) -c $^
%_wrap.o: %.c %_wrap.c
	@

_%.so: %.o %_wrap.o
	$(LD) $(LFLAGS) $^ -o $@

sg_align.o: sg_align.c sg_align.h
	$(CC) $(DFLAGS) -c $^

debug: debug.c sg_align.o
	$(CC) $(DFLAGS) -o $@ $^
