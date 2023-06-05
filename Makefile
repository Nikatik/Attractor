CC := cc
CXX := c++
SHELL := $(shell whereis bash | awk '{print $$2}')

COMPLIER := 0
CSOURCES := $(shell find . -name '*.c')
CPPSOURCES := $(shell find . -name '*.cpp')

ifeq ($(CC),gcc)
    COMPLIER := 1
else ifeq ($(CC), clang)
    COMPLIER := 2
else
endif

ifeq ($(CXX),g++)
    COMPLIER := 1
else ifeq ($(CXX), clang++)
    COMPLIER := 2
else
endif

FLAGS := -Wall -Wextra -O -Wfloat-equal -MMD -Wundef -Wshadow -Wcast-align -Wconversion -Wunreachable-code -ftrapv -Wformat=2 -Wreturn-type -Wdouble-promotion -Wstrict-overflow=5 -Wconversion
LDLIBS := -lm -lstdc++ -lquadmath 

all: FLAGS += -O3
debug: FLAGS += -g -Og
test: FLAGS += -Werror -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=null -fstack-protector-all -Og -g
test: LDFLAGS := -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=null -fstack-protector-all
suck: FLAGS := -Wall -Wextra -Wfloat-equal -O -MMD -Og -g 

ifeq ($(COMPLIER),1)
    FLAGS += -Wshift-overflow=2 -Wformat-overflow=2
    test: LDFLAGS += -fsanitize=bounds-strict
    test: FLAGS += -fsanitize=bounds-strict
else ifeq ($(COMPLIER),2)
    FLAGS += -Wshift-overflow -Wpedantic -pedantic
else
endif

CFLAGS := -std=c2x
CXXFLAGS := -std=c++2a

all: main
debug: main
test: main
suck: main

main: CFLAGS += $(FLAGS)
main: CXXFLAGS += $(FLAGS)
main: $(CSOURCES:%.c=%.o) $(CPPSOURCES:%.cpp=%.o) 

DEPS := $(shell find -name "*.d")
-include $(DEPS)

clean:
	@rm -vf main *.o *.d trajectory*
