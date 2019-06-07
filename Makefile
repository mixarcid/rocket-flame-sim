CPP = gcc -Wall 
SRCS = main.c solver3d.c perlin.c

all:
	$(CPP) $(SRCS) -o fluid -lm

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
