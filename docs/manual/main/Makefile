all: finalstep

pdf1:
	pdflatex main.tex -halt-on-error
pdf2: pdf1
	pdflatex main.tex
pdf3: pdf2
	pdflatex main.tex
finalstep: pdf3
	mv main.pdf ../Manual.pdf

clean:
	rm main.aux  main.log  main.out  main.toc ../Manual.pdf
