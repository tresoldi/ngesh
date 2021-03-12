pandoc -V graphics="true" -V geometry:margin=1in --pdf-engine=xelatex --bibliography=paper.bib --csl=apa.csl --template latex.template --filter pandoc-citeproc -o paper.tex paper.md
latexmk -xelatex paper.tex

