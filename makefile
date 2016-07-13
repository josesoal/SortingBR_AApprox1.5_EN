CC = g++
CFLAGS = -g -Wall # flags when *compiling*
LFLAGS = -g -Wall # flags when *linking*
LIBS = -lm # math library
INCLUDE = -I boost_1_48_0 
SOURCES = breakpointGraph.c cycleDecomposition.c reversalGraph.c eliminate.c main.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = approx

all: $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@ $(LIBS)

%.o:%.c
	$(CC) $(INCLUDE) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
