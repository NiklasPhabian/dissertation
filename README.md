# Dissertation

# Dependencies
- pandoc
- pandoc-citeproc
- texlive-latex-extra
- texlive-science
- texlive-bibtex-extra
- make
- librsvg2-bin might be required to include the SVGs; at least it made the error: "check that rsvg-convert is in path." go away ...

## Glossaries
https://github.com/tomncooper/pandoc-gls


## For tex
```bash
sudo apt install texlive-bibtex-extra texlive-science texlive-latex-extra biber
```
    
## Biber
Using lates built of biber (in full/). Needs libsombok-dev installed

To build new, follow https://github.com/plk/biber



    
# Build 

    make

## Display warnings only

    pdflatex --interaction=nonstopmode <filename> | awk 'BEGIN{IGNORECASE = 1}/warning|!/,/^$/;'

    

# Convert tex to md
    
## tex > md

    pandoc
    
## Replace latex citations with md citations

    %s/\\cite{\([^}]*\)}/[@\1]/g 
    

## Replace images

    %s/<img src="\([^"]*\)".*alt="\([^"]*\)" \/>/![\2](\1)/g 


## Slash FU

    %s/\\\([\[\]]\)/\1/g 

    
