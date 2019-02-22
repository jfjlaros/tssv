EXEC := run_tests
MAIN := test_sgAlign
TESTS :=
LIBS := ../tssv/sgAlign
FIXTURES :=


CC := g++
INCLUDE :=
CC_ARGS := -O3
V_ARGS := -U __SSE2__


OBJS := $(addsuffix .o, $(TESTS) $(FIXTURES) $(LIBS))
V_OBJS := $(addsuffix _v.o, $(TESTS) $(FIXTURES) $(LIBS))

.PHONY: all clean distclean


all: $(EXEC) $(EXEC)_v

$(EXEC): $(MAIN).cc $(OBJS)
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ $^

$(EXEC)_v: $(MAIN).cc $(V_OBJS)
	$(CC) $(V_ARGS) $(CC_ARGS) $(INCLUDE) -o $@ $^

%.o: %.cc
%.o: %.c
	$(CC) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

%_v.o: %.cc
%_v.o: %.c
	$(CC) $(V_ARGS) $(CC_ARGS) $(INCLUDE) -o $@ -c $<

check: all
	./$(EXEC); ./$(EXEC)_v

clean:
	rm -f $(OBJS) $(V_OBJS)

distclean: clean
	rm -f $(EXEC) $(EXEC)_v
