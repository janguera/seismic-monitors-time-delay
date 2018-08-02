report.pdf: report.tex refs.bib ecmi-logo.png
	pdflatex report.tex
	bibtex report
	pdflatex report.tex


