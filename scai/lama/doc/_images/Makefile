
IMAGES = matrix_distribution_row matrix_distribution_col

IMAGES_SVG = $(IMAGES:=.svg)

all: $(IMAGES_SVG)

clean:
	rm -rf $(IMAGES_SVG) *.bak

.SUFFIXES: .svg .fig

%.svg: %.fig
	fig2dev -Lsvg  < $< > $@
