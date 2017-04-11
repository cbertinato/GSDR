CC=g++
CFLAGS=-c -Wall
LDFLAGS=-lfftw3 -lm
SOURCES=acquire.cpp generateCA.cpp vector_math.cpp
LIBDIR=/usr/local/lib
INCLUDEDIR=/usr/local/include
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=acquire

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) -I$(INCLUDEDIR) -L$(LIBDIR) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) -I$(INCLUDEDIR) -L$(LIBDIR)  $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)
