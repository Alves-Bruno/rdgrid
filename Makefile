all:
	g++ rdgrid.cpp -g -DNORCPP -I./ -o rdgrid
clean:
	rm rdgrid out.*
