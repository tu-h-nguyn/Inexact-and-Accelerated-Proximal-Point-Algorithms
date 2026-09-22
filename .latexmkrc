# Cau hinh latexmk cho repo nay.
$pdf_mode      = 1;          # pdflatex
$bibtex_use    = 2;          # chay bibtex, don dep .bbl khi `latexmk -C`
$max_repeat    = 5;          # du vong lap cho muc luc + trich dan on dinh
$out_dir       = '.';
$clean_ext     = 'bbl nav snm vrb run.xml synctex.gz';
