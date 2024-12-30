all: report

report: report_design report_report

report_design:
	pdflatex design
	pdflatex design

report_report:
	pdflatex report
	pdflatex report

clean:
	rm -f *.pdf *.bbl *.blg *.toc *.out *.run.xml *.log *.aux *.bib
