all: run

run: compile execute

compile:
	g++ -o ppForm_linear.exe ppForm_linear.cpp
	g++ -o BSpline_linear.exe BSpline_linear.cpp
	g++ -o Task_A.exe Task_A.cpp -llapacke -llapack -lblas
	g++ -o Compare_two_splines.exe Compare_two_splines.cpp -llapacke -llapack -lblas
	g++ -o Task_C.exe Task_C.cpp -llapacke -llapack -lblas
	g++ -o BSpline_arbitrary_order.exe BSpline_arbitrary_order.cpp
	g++ -o Task_D.exe Task_D.cpp -llapacke -llapack -lblas
	@echo "\n----All programs have been complied.----\n"

execute:
	./ppForm_linear.exe < ../data/ppForm_linear.in | python draw.py
	./BSpline_linear.exe < ../data/BSpline_linear.in | python draw.py
	./Task_A.exe | python draw_TaskA.py
	./Compare_two_splines.exe < ../data/Compare_two_splines.in | python draw.py
	./Task_C.exe | python draw_TaskC.py
	./BSpline_arbitrary_order.exe < ../data/BSpline_arbitrary_order.in | python draw.py
	./Task_D.exe > ../data/Task_D.out
	@echo "\n----All programs have been executed.----\n"

clean:
	rm -f *.exe
