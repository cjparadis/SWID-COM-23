#plots outputs from SWID COM 23

# =================================================================================================
# Plot 1 Injection Phase Concentration
# =================================================================================================

set terminal pngcairo
set output 'plot_conc_s1.png'
set title "Injection"
set xlabel "x (ft)"
set ylabel "C (mg/L)"
set key below

plot "swid_out_s1.conc" using 1:2 title "sim1",\
"swid_out_s1.conc" using 1:3 title "sim2",\
"swid_out_s1.conc" using 1:4 title "sim3",\
"swid_out_s1.conc" using 1:5 title "ana1",\
"swid_out_s1.conc" using 1:6 title "ana2",\
"swid_out_s1.conc" using 1:7 title "ana3"

# =================================================================================================	 
# Plot 2 Drift Phase Concentration
# =================================================================================================

set terminal pngcairo
set output 'plot_conc_s2.png'
set title "Drift"
set xlabel "t (days)"
set ylabel "C (mg/L)"
set key below

plot "swid_out_s2.conc" using 1:2 title "sim1",\
"swid_out_s2.conc" using 1:3 title "sim2",\
"swid_out_s2.conc" using 1:4 title "sim3",\
#"swid_hoss_s2.conc" using 1:($2*500) title "hoss1",\
#"swid_hoss_s2.conc" using 1:($3*500) title "hoss2"

# =================================================================================================
# Plot 3 Injection Phase Mass Budget
# =================================================================================================

set terminal pngcairo
set output 'plot_disc_s1.png'
set title "Injection"
set xlabel "t (hours)"
set ylabel "Disc. (%)"
set key below

plot "swid_out_s1.disc" using 1:2 title "sim1",\
"swid_out_s1.disc" using 1:4 title "sim2",\
"swid_out_s1.disc" using 1:6 title "sim3"

# =================================================================================================	 
# Plot 4 Drift Phase Mass Budget
# =================================================================================================

set terminal pngcairo
set output 'plot_disc_s2.png'
set title "Drift"
set xlabel "t (days)"
set ylabel "Disc. (%)"
set key below

plot "swid_out_s2.disc" using 1:2 title "sim1",\
"swid_out_s2.disc" using 1:4 title "sim2",\
"swid_out_s2.disc" using 1:6 title "sim3"