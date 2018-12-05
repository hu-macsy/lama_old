
IMAGES = splitting block_distribution cyclic2_distribution genblock_distribution general_distribution grid_distribution shift

IMAGES_SVG = $(IMAGES:=.svg)

all: $(IMAGES_SVG)

clean:
	rm -rf $(IMAGES_SVG) *.bak

.SUFFIXES: .svg .fig

%.svg: %.fig
	fig2dev -Lsvg  < $< > $@
