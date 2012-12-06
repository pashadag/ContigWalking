COPT =  -Wall -Wno-sign-compare -fPIC -fexceptions 
#COPT =  -Wall -Wno-sign-compare -O3 -fPIC -fexceptions 

# redefine inbuilt macros and define new ones
CC    = g++
CFLAGS = $(COPT) 


#.SUFFIXES:       # remove inbuilt definitions
#.SUFFIXES: .cpp .o # define default compilation procedure
#.cpp.o:            # .o files depend on .c files
#	$(CC) $(CFLAGS) $*.c # start with a tab character!!

# Default target
all: walks
#	@echo "compilation done"

debug: COPT += -g
debug: all

walks: walks.o
	$(CC) walks.o -o walks
walks.o : walks.cpp 
	$(CC) -c $(COPT) walks.cpp -o walks.o

# Target deleting unwanted files
clean:
	rm -f *.o *~ walks core mppcore
