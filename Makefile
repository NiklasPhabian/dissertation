SOURCES := C1_OCCUR.md C2_Software.md C3_Spires.md

PDFS = $(SOURCES:.md=.pdf) 


%.pdf : %.md	
	pandoc --pdf-engine=xelatex --toc --toc-depth=5 --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --listings $^ -o $@
	#pandoc --pdf-engine=xelatex --toc --toc-depth=5 --citeproc --lua-filter context/pandoc-gls.lua --listings $^ -o $@
	#pandoc --from markdown --to latex --filter pandoc-citeproc -o output.pdf --lua-filter context/pandoc-gls.lua $^
	


all : $(PDFS) 	

 
clean:
	rm -fv $(PDFS) 	
