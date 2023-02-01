CC = c++
CFLAGS = -Wall # -g
FLAGS = -Wall -O2
LIBS = -lm
OBJS = SCPv.o skcp_main.o


mmas_ml: $(OBJS)
	$(CC) $(FLAGS) -o skcp_main $(OBJS) $(LIBS)
.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	/bin/rm -rf *.o *~ skcp_main $(OBJS) $(TARGET)
