EXEC := run_tests
MAIN := test_lib
TESTS := test_sgAlign test_sgAlignSSE
LIBS := ../tssv/sgAlign ../tssv/sgAlignSSE
FIXTURES :=


CC := g++
INCLUDE :=
CC_ARGS := -O3 -std=c++11


OBJS := $(addsuffix .o, $(TESTS) $(FIXTURES) $(LIBS))

.PHONY: all clean distclean


all: $(EXEC)

$(EXEC): $(MAIN).cc $(OBJS)
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ $^

%.o: %.cc
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

%.o: %.c
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

check: all
	./$(EXEC)

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(EXEC)
