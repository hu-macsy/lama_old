# This makefile is used to generate pdf and svg images from the fig files

IMAGES = bucket_sort gather

IMAGES_SVG = $(IMAGES:=.svg)
IMAGES_PDF = $(IMAGES:=.pdf)

all: $(IMAGES_PDF) $(IMAGES_SVG)

clean:
	rm -rf *.bak

distclean:
	rm -rf $(IMAGES_SVG) $(IMAGES_PDF) *.bak

.SUFFIXES: .svg .fig
.SUFFIXES: .pdf .fig

%.svg: %.fig
	fig2dev -Lsvg  < $< > $@

%.pdf: %.fig
	fig2dev -Lpdf  < $< > $@
