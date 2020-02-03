#set encoding iso_8859_1
reset

# set title "" 0.0,3.0
#set key title 5,5

 set xlabel "E [keV]"
#################
#set xrange [0:100]
#################
#set autoscale x
#set xzeroaxis lt -1 lw 2
#set xtics nomirror 1
# set mxtics 10 
# set format x "10^{%L}"
#set logscale x


######################################

set ylabel "Intensity []"
##################
#set yrange [0.000000001:1]
##################
#set autoscale y
#set yzeroaxis lt -1 lw 2
#set ytics .1
# set mytics 10
#set format y "10^{%L}"
#set logscale y


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
#set term wxt persist
#set term qt 
set output "./SynSpec_out_lin.svg"
#fit f(x) 's_data.dat' via a, w
plot "SynSpec_out.txt" using 1:2 title "" w lines lw 2

#set yrange [1E-10:2]
set mytics 10
set format y "10^{%L}"
set logscale y
set yrange [1E-2:1]
set output "./SynSpec_out_log.svg"
plot "SynSpec_out.txt" using 1:2 notitle w lines  lw 2 lt 1,\
     "SynSpec_out.txt" using 1:2 notitle w points ps 1 lt 1


#set terminal qt persist
#plot "SynSpec_out.txt" using 1:2 title "" w lines lw 2

