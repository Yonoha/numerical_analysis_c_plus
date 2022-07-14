set title 'boost::odeint Lotka Volterra'

set xl 'time'
set yl 'x and y'

set xr [0:5.1]
set yr [0.3:1.2]

set grid

plot 'data_boost_odeint_2dim_lotka_volterra.txt' u 1:2 w l title 'x',\
'data_boost_odeint_2dim_lotka_volterra.txt' u 1:3 w l title 'y'

set terminal pngcairo
set output 'boost_odeint_2dim_lotka_volterra.png'
replot 
