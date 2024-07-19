all: lecture_0.pdf lecture_5.pdf lecture_6.pdf lecture_7.pdf lecture_8.pdf lecture_13.pdf

%.pdf: %.md
	pandoc -t beamer -V theme:Boadilla --pdf-engine=lualatex -V fontsize=10pt -V 'mainfont:NewComputerModernSans10' -fmarkdown-implicit_figures -H header.tex -s $< -o $@

.PHONY: clean

clean:
	rm -f *.pdf

.PHONY: clean
