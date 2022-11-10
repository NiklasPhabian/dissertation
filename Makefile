.DEFAULT_GOAL := all
SOURCES := C0_Introduction.md C1_OCCUR.md C2_Software.md C3_Spires.md

PDFS = $(SOURCES:.md=.pdf) 
LATEX = $(SOURCES:.md=.tex)


%.pdf : %.md	
	pandoc --pdf-engine=xelatex --toc --toc-depth=5 --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --listings $^ -o $@
	#pandoc --pdf-engine=xelatex --toc --toc-depth=5 --citeproc --lua-filter context/pandoc-gls.lua --listings $^ -o $@
	#pandoc --from markdown --to latex --filter pandoc-citeproc -o output.pdf --lua-filter context/pandoc-gls.lua $^


%.tex: %.md
	pandoc --write=latex --biblatex --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --listings $^ -o full/chapters/$@	
	

pdflatex: 
	rsync -rupE images/ full/images/
	make -C full/

all : $(PDFS) $(LATEX) pdflatex

clean:
	rm -fv $(PDFS) 
