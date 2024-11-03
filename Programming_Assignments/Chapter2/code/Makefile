run: compile execute

compile:
	g++ Task_B.cpp -o Task_B.exe -Wall
	g++ Task_C.cpp -o Task_C.exe -Wall
	g++ Task_D.cpp -o Task_D.exe -Wall
	g++ Task_E.cpp -o Task_E.exe -Wall
	g++ Task_F.cpp -o Task_F.exe -Wall

execute:
	./Task_B.exe | python Task_B_draw.py
	./Task_C.exe | python Task_C_draw.py
	./Task_D.exe < ../data/Task_D_data.txt > ../results/Task_D/Task_D_output.txt
	./Task_E.exe < ../data/Task_E_data.txt
	python Task_E_draw.py < ../results/Task_E/Task_E_scatter.txt
	./Task_F.exe
	python Task_F_draw.py < ../results/Task_F/scatter_plot_10.txt
	python Task_F_draw.py < ../results/Task_F/scatter_plot_40.txt
	python Task_F_draw.py < ../results/Task_F/scatter_plot_160.txt

clean:
	rm -f *.exe