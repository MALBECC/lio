all: finalstep

pdf0:
	pdflatex main.tex -halt-on-error
pdf1: pdf0
	biber main
pdf2: pdf1
	pdflatex main.tex
pdf3: pdf2
	pdflatex main.tex
finalstep: pdf3
	mv main.pdf ../Manual.pdf

clean:
	rm main.aux  main.log  main.out  main.toc main.aux.blg  main.bbl  main.bcf  main.blg  main.pdf  main.run.xml main.tex.blg ../Manual.pdf
