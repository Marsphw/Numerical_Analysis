all: run

run: compile execute

compile:
	g++ ppForm_linear.cpp -o ppForm_linear -Wall
	g++ BSpline_linear.cpp -o BSpline_linear -Wall

execute:
	./ppForm_linear < ../data/ppForm_linear.in | python draw.py
	./BSpline_linear < ../data/BSpline_linear.in | python draw.py
		

clean:
	rm -f *.exe
