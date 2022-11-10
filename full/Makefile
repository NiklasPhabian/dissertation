basename=dissertation

$(basename).pdf: $(basename).tex chapters/* tex/*
	pdflatex        	$(basename).tex
	makeglossaries  	$(basename)
	./biber-linux_x86_64      $(basename)
	pdflatex        	$(basename).tex
	pdflatex        	$(basename).tex

clean:
	rm -fv $(basename).[a-oq-su-z]* $(basename).tex.* $(basename).toc $(basename).pdf

realclean: clean
	rm -fv $(basename).pdf
