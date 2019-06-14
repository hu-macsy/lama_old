
IMAGES = matrix_distribution_row matrix_distribution_col sparse_vector dense_vector \
         single_io independent_io collective_io file_io

IMAGES_SVG = $(IMAGES:=.svg)
IMAGES_PDF = $(IMAGES:=.pdf)

all: $(IMAGES_SVG) $(IMAGES_PDF) 

clean:
	rm -rf $(IMAGES_SVG) $(IMAGES_PDF) *.bak

.SUFFIXES: .svg .fig
.SUFFIXES: .pdf .fig

%.svg: %.fig
	fig2dev -Lsvg  < $< > $@

%.pdf: %.fig
	fig2dev -Lpdf  < $< > $@
