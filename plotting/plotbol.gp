# snec_lbol.gp
lumfile = "Data/lum_observed.dat"

BASE_FONT  = "Helvetica,20"
TITLE_FONT = ",26"
LABEL_FONT = ",22"
KEY_FONT   = ",20"
TIC_FONT   = ",20"

set terminal qt size 1100,750 enhanced font BASE_FONT title "SNEC Bolometric LC"

set title  "SNEC Bolometric Light Curve" font TITLE_FONT
set xlabel "Time since explosion [days]" font LABEL_FONT
set ylabel "Bolometric luminosity [erg/s]" font LABEL_FONT
set key font KEY_FONT
set tics font TIC_FONT
set border lw 3
set grid lw 3

# Normal y direction:
set yrange [*:*] noreverse

# Often useful:
set logscale y

plot \
    lumfile using ($1/86400.0):2 with lines lw 4 lc rgb "#0000cc" title "L_{bol}"

pause mouse close

