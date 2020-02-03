#set encoding iso_8859_1
reset

# set title "" 0.0,3.0
#set key title 5,5

set xlabel "E [keV]"
#################
set xrange [1:1E6]
#################
#set autoscale x
#set xzeroaxis lt -1 lw 2
#set xtics nomirror 1
set mxtics 10
set format x "10^{%L}"
set logscale x


######################################

set ylabel "Cross section [cm2/g]"
##################
set yrange [0.001:100]
##################
#set autoscale y
#set yzeroaxis lt -1 lw 2
#set ytics .1
set mytics 10
set format y "10^{%L}"
set logscale y


#set size ratio -1
######################################
#set object circle at 1,1 radius 0.05
#set label "P"  at 1.15,0.7
#set object circle at -3,9 radius 0.05
#set label "P'"  at -3.4,8.5
#set object circle at -1,-3 radius 0.05
#set label "Q"  at -0.8,-3.2
#set arrow 1 from 0.523 , 0.5 to 0.523 , 0 nohead                      
######################################

set grid
#set sample 500
#set key 4.8,-1

#set terminal postscript eps enhanced color size 10,7  "Arial" 40
#set terminal pdf

set terminal svg size 800,400
#set terminal wxt persist
file = "SynSpec_material_CdTe.txt"

# PHOTON       SCATTERING         PHOTO-    PAIR PRODUCTION    TOTAL ATTENUATION
# ENERGY    COHERENT INCOHER.    ELECTRIC    IN        IN       WITH     WITHOUT
#                               ABSORPTION NUCLEAR  ELECTRON  COHERENT  COHERENT
#                                           FIELD     FIELD    SCATT.    SCATT.
#   (MeV)    (cm2/g)   (cm2/g)   (cm2/g)   (cm2/g)   (cm2/g)   (cm2/g)   (cm2/g)

set output "./interaction.svg"
plot file using ($1*1000):3 title "incoherent" w lines lw 2,\
     file using ($1*1000):4 title "photoelectric" w lines lw 2,\
     file using ($1*1000):($3+$4) title "photoelectric + incoherent" w lines lw 2,\
     file using ($1*1000):($5+$6) title "pair" w lines lw 2,\
     file using ($1*1000):7 title "total" w lines lw 2
