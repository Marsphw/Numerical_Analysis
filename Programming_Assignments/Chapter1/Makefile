all: run report

run:
	$(MAKE) -C code run

report:
	$(MAKE) -C tex report

clean:
	$(MAKE) -C code clean
	$(MAKE) -C tex clean
	rm -f result/*.txt