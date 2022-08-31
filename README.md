# Dissertation

# Dependencies
- librsvg2-bin might be required to include the SVGs; at least it made the error: "check that rsvg-convert is in path." go away ...

# Glossaries
https://github.com/tomncooper/pandoc-gls


## For tex

sudo apt install texlive-bibtex-extra texlive-science texlive-latex-extra biber
    
# Build with latex
    
    biber occur
    makeglossaries occur.glo
    pdflatex occur.tex
    
    
## tex > md

    pandoc
    
## Replace latex citations with md citations

    %s/\\cite{\([^}]*\)}/[@\1]/g 
    

## Replace images

    %s/<img src="\([^"]*\)".*alt="\([^"]*\)" \/>/![\2](\1)/g 


## Slash FU

    %s/\\\([\[\]]\)/\1/g 

    
