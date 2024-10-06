all: report

report:
	xelatex Report#1.tex
	bibtex Report#1
	xelatex Report#1.tex
	xelatex Report#1.tex

clean:
	rm -f *.pdf *-blx.bib *.aux *.bbl *.blg *.out *.run.xml *.log