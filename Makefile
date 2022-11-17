.DEFAULT_GOAL := all

MDDIR=mds
CHAPTERDIR=chapters
STANDALONEDIR=chapters_standalone
DISSERTATIONDIR=dissertation

MDS =$(wildcard $(MDDIR)/*.md)
TEXCHAPTERS=$(patsubst $(MDDIR)/%.md,$(CHAPTERDIR)/%.tex,$(MDS))
TEXSTANDALONES=$(patsubst $(MDDIR)/%.md,$(STANDALONEDIR)/%.tex,$(MDS))
STANDALONEPDFS=$(patsubst $(STANDALONEDIR)/%.tex,$(STANDALONEDIR)/%.pdf,$(TEXSTANDALONES))

.PHONY: tex pdfs clean 

chapters: $(TEXCHAPTERS)
texstandalone: $(TEXSTANDALONES)
standalone: $(STANDALONEPDFS) 
dissertation.pdf: $(TEXCHAPTERS) tex/* dissertation.tex
	pdflatex -shell-escape -output-directory $(DISSERTATIONDIR) $(basename $@).tex | awk 'BEGIN{IGNORECASE = 1}/warning|!/,/^$$/;' 	
	makeglossaries -d $(DISSERTATIONDIR) $(basename $(notdir $@))
	./biber-linux_x86_64 $(DISSERTATIONDIR)/$(basename $@)
	pdflatex -shell-escape -output-directory $(DISSERTATIONDIR) $(basename $@).tex | awk 'BEGIN{IGNORECASE = 1}/warning|!/,/^$$/;'
	mv $(DISSERTATIONDIR)/dissertation.pdf dissertation.pdf


$(STANDALONEDIR)/%.pdf: $(STANDALONEDIR)/%.tex
	pdflatex -shell-escape -output-directory $(STANDALONEDIR) $^ | awk 'BEGIN{IGNORECASE = 1}/warning|!/,/^$$/;' 	
	makeglossaries -d $(STANDALONEDIR) $(basename $(notdir $^))
	./biber-linux_x86_64  $(basename $^)
	pdflatex -shell-escape -output-directory $(STANDALONEDIR) $^ | awk 'BEGIN{IGNORECASE = 1}/warning|!/,/^$$/;' 			
	
	
$(STANDALONEDIR)/%.tex: $(MDDIR)/%.md	
	pandoc -s --to=latex --no-highligh --biblatex --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --lua-filter context/minted.lua $^ -o $@
	#pandoc -s --to=latex --no-highligh --biblatex --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --listings $^ -o $@
	

$(CHAPTERDIR)/%.tex: $(MDDIR)/%.md
	pandoc --to=latex --biblatex --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --lua-filter context/minted.lua $^ -o $@
	#pandoc --to=latex --biblatex --filter pandoc-citeproc --lua-filter context/pandoc-gls.lua --lua-filter --listings $^ -o $@

cleanstandalone:
	rm -fv $(STANDALONEDIR)/*

cleandissertation:
	
	rm -fv $(DISSERTATIONDIR)/dissertation.[a-oq-su-z]* 
	rm -fv $(DISSERTATIONDIR)/dissertation.tex.*
	rm -fv $(DISSERTATIONDIR)/dissertation.toc
	rm -fv $(DISSERTATIONDIR)/dissertation.pdf
	rm -fvr $(DISSERTATIONDIR)/_minted-dissertation
	
test: 
	@echo $(TEXSTANDALONES)
	
