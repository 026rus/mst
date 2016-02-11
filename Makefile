CC=mpic++ -Wall
CLEAN=rm -f
# Name of the program
PM=mst1
# Object files
OF=*.o
# Source files
SF= main.cpp


$(PM):  $(OF)
	$(CC) -o $@ $^

$(OF): $(SF)
	$(CC) -c $^

clean:
	$(CLEAN) *.o
	$(CLEAN) $(PM)
