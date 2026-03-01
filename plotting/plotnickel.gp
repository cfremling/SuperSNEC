# plot_ni_profiles.gp
# Visualize initial 56Ni distribution from SNEC

ni_file   = "Data/Ni_init_frac.dat"
mass_file = "Data/mass_initial.dat"
rad_file  = "Data/rad_initial.dat"
vel_file  = "Data/vel.xg"

# Physical constants
MSUN = 1.989e33     # [g]
RSUN = 6.957e10     # [cm]

# Target time for velocity slice (2 days after explosion)
T_TARGET_DAYS = 2.0
T_TARGET_SEC  = T_TARGET_DAYS * 86400.0

set terminal qt size 1200,900 enhanced font "Helvetica,20"
set grid lw 3
set tics font ",20"
set border lw 3

KEY_FONT = ",20"
LABEL_FONT = ",22"
TITLE_FONT = ",26"

set key font KEY_FONT

set multiplot layout 2,2 title "Initial 56Ni distribution (SNEC)" font TITLE_FONT

# --------------------------------------------------------------------
# 1) X(Ni) vs zone index
# Ni_init_frac.dat: col1 = zone index, col2 = X(Ni)
# --------------------------------------------------------------------
set xlabel "Zone index" font LABEL_FONT
set ylabel "X(^{56}Ni)" font LABEL_FONT
set title  "X(^{56}Ni) vs zone index" font TITLE_FONT

plot ni_file using 1:2 with lines lw 4 title "X(^{56}Ni)"

# --------------------------------------------------------------------
# 2) X(Ni) vs enclosed mass [M_sun]
# mass_initial.dat: col1 = zone index, col2 = mass [g]
# Ni_init_frac.dat: col1 = zone index, col2 = X(Ni)
# After paste: col2 = mass, col4 = X(Ni)
# --------------------------------------------------------------------
set xlabel "Enclosed mass [M_{sun}]" font LABEL_FONT
set ylabel "X(^{56}Ni)" font LABEL_FONT
set title  "X(^{56}Ni) vs mass" font TITLE_FONT

plot "< paste Data/mass_initial.dat Data/Ni_init_frac.dat" \
     using ($2/MSUN):4 with lines lw 4 title "X(^{56}Ni)(m)"

# --------------------------------------------------------------------
# 3) X(Ni) vs radius [R_sun]
# rad_initial.dat: col1 = zone index, col2 = radius [cm]
# Ni_init_frac.dat: col1 = zone index, col2 = X(Ni)
# After paste: col2 = radius, col4 = X(Ni)
# --------------------------------------------------------------------
set xlabel "Radius [R_{sun}]" font LABEL_FONT
set ylabel "X(^{56}Ni)" font LABEL_FONT
set title  "X(^{56}Ni) vs radius" font TITLE_FONT

plot "< paste Data/rad_initial.dat Data/Ni_init_frac.dat" \
     using ($2/RSUN):4 with lines lw 4 title "X(^{56}Ni)(r)"

# --------------------------------------------------------------------
# 4) X(Ni) vs velocity [km/s] at ~2 days
#
# We use the AWK helper to:
#   - read n = #zones from Ni_init_frac.dat
#   - scan vel.xg and pick the block whose Time is closest to T_TARGET_SEC
#   - output mass, velocity for that block
# Then we paste those with Ni_init_frac.dat and plot v (km/s) vs X(Ni).
#
# After paste: col1 = mass, col2 = vel [cm/s], col3 = zone, col4 = X(Ni)
# --------------------------------------------------------------------
set xlabel sprintf("Velocity [km/s] (t \\approx %.1f d)", T_TARGET_DAYS) font LABEL_FONT
set ylabel "X(^{56}Ni)" font LABEL_FONT
set title  "X(^{56}Ni) vs velocity" font TITLE_FONT

cmd_v2d = sprintf("< awk -v ttarget=%g -f plotting/vel_at_time.awk %s %s | paste - %s", \
                  T_TARGET_SEC, ni_file, vel_file, ni_file)

plot cmd_v2d using ($2/1e5):4 with lines lw 4 title "X(^{56}Ni)(v, t\\~2d)"

unset multiplot

pause mouse close
