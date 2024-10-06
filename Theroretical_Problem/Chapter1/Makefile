all: report

report:
	xelatex Report.tex
	bibtex Report
	xelatex Report.tex
	xelatex Report.tex

clean:
	rm -f *.pdf *-blx.bib *.aux *.bbl *.blg *.out *.run.xml *.log