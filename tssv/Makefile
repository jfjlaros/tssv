EXEC := d
MAIN := debug
TESTS :=
LIBS := sg_align
FIXTURES :=


CC := g++
INCLUDE :=
CC_ARGS := -O3
V_ARGS := -U __SSE2__


OBJS := $(addsuffix .o, $(TESTS) $(FIXTURES) $(LIBS))
V_OBJS := $(addsuffix _v.o, $(TESTS) $(FIXTURES) $(LIBS))

.PHONY: clean distclean


all: $(EXEC) $(EXEC)_v

$(EXEC): $(MAIN).c $(OBJS)
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ $^

$(EXEC)_v: $(MAIN).c $(V_OBJS)
	$(CC) $(V_ARGS) $(CC_ARGS) $(INCLUDE) -o $@ $^

%.o: %.cc
%.o: %.c
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

%_v.o: %.cc
%_v.o: %.c
	$(CC) $(V_ARGS) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

clean:
	rm -f $(OBJS) $(V_OBJS)

distclean: clean
	rm -f $(EXEC) v$(EXEC)
