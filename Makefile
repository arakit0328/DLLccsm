CC = c++
CFLAGS = -Wall -g
FLAGS = -Wall -g
LIBS = -lm
OBJS = SCPv.o skcp.o


mmas_ml: $(OBJS)
	$(CC) $(FLAGS) -o skcp $(OBJS) $(LIBS)
.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	/bin/rm -rf *.o *~ skcp $(OBJS) $(TARGET)
