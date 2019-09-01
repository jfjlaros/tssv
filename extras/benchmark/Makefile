EXEC := benchmark
MAIN := benchmark
TESTS :=
LIBS := ../../tssv/sgAlign ../../tssv/sgAlignSSE
FIXTURES :=


CC := gcc
INCLUDE :=
CC_ARGS := -O3


OBJS := $(addsuffix .o, $(TESTS) $(FIXTURES) $(LIBS))

.PHONY: clean distclean


all: $(EXEC)

$(EXEC): $(MAIN).c $(OBJS)
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ $^

%.o: %.c
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

clean:
	rm -f $(OBJS) $(V_OBJS)

distclean: clean
	rm -f $(EXEC)
