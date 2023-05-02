.PHONY: main clean FORCE

main: poster.pdf

poster.pdf: FORCE
	latexmk  -shell-escape -pdflatex='lualatex -f -interaction nonstopmode' -pdf poster.tex

clean:
	latexmk -pdf -C
