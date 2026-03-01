# snec_mag.gp — Plot g,r band light curves from bolometric-correction magnitudes
magfile = "Data/magnitudes.dat"

BASE_FONT  = "Helvetica,20"
TITLE_FONT = ",26"
LABEL_FONT = ",22"
KEY_FONT   = ",20"
TIC_FONT   = ",20"

set terminal qt size 1100,750 enhanced font BASE_FONT title "SNEC g,r Light Curves"
set title  "SNEC g,r Light Curves (bolometric correction)" font TITLE_FONT
set xlabel "Time since explosion [days]" font LABEL_FONT
set ylabel "Magnitude" font LABEL_FONT
set key font KEY_FONT
set tics font TIC_FONT
set border lw 3
set grid lw 3

# astronomical convention: brighter = up
set yrange [] reverse

# Column mapping for magnitudes.dat:
# time(s)  Teff  PTF_R  u   g   r   i ...
# 1        2     3      4   5   6   7

set style line 1 lt 1 lc rgb "#008800" pt 7 ps 2.4 lw 4   # g: dark green
set style line 2 lt 1 lc rgb "#cc0000" pt 7 ps 2.4 lw 4   # r: dark red

plot \
    magfile using ($1/86400.0):5 with linespoints ls 1 title "g", \
    magfile using ($1/86400.0):6 with linespoints ls 2 title "r"

pause mouse close
