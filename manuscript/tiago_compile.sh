#!/bin/bash
pandoc -V graphics="true" -V geometry:margin=1in --pdf-engine=xelatex --bibliography=paper.bib -o paper.pdf paper.md
