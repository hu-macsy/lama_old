
matrix = SparseMatrix.png 

all: ${matrix}

.SUFFIXES: .fig .pdf
.SUFFIXES: .fig .png
.SUFFIXES: .eps .pdf

.fig.pdf:
	fig2dev -L pdf $*.fig > $*.pdf 

.fig.png:
	fig2dev -L png $*.fig > $*.png

.eps.pdf:
	epstopdf $*.eps

clean:
	rm -rf ${matrix}
