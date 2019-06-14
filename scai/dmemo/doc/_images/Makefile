# This makefile is used to generate pdf and svg images from the fig files

IMAGES = splitting block_distribution cyclic2_distribution genblock_distribution \
         general_distribution grid_distribution shift \
         compute_owners compute_owners_general \
         communication_plan exchange_by_plan \
         global_exchange global_exchange_plan \
         global_addressing_plan \
         global_gather global_scatter global_redistribute \
         redistribution redistribute_plan \
         halo halo_plan halo_update halo_reduce \
         collective_io independent_io

IMAGES_SVG = $(IMAGES:=.svg)
IMAGES_PDF = $(IMAGES:=.pdf)

all: $(IMAGES_PDF) $(IMAGES_SVG)

clean:
	rm -rf $(IMAGES_SVG) $(IMAGES_PDF) *.bak

.SUFFIXES: .svg .fig
.SUFFIXES: .pdf .fig

%.svg: %.fig
	fig2dev -Lsvg  < $< > $@

%.pdf: %.fig
	fig2dev -Lpdf  < $< > $@
