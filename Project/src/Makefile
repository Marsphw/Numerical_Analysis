all: run

run: compile execute

compile:
	g++ ppForm_linear.cpp -o ppForm_linear
	g++ BSpline_linear.cpp -o BSpline_linear
	g++ -o ppForm_cubic ppForm_cubic.cpp -llapacke -llapack -lblas
	@echo "All programs have been complied."

execute:
	./ppForm_linear < ../data/ppForm_linear.in | python draw.py
	./BSpline_linear < ../data/BSpline_linear.in | python draw.py
	./ppForm_cubic | python draw_TaskA.py
		

clean:
	rm -f *.exe
