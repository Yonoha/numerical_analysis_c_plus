set title '1dim Burgers eq by fem'

set term gif animate optimize delay 5 size 700, 480
set output 'fem_1dim_burgers.gif'

set xl 'x'
set yl 'u'

set grid

XRANGE = 1.0
YRANGE = 1.0

do for [i = 0:100]{
    filename = sprintf('data_fem_1dim_burgers_%d.txt', i)
    plot_title = sprintf('time number = %d', i)
    set title plot_title
    plot [0:XRANGE] [0:YRANGE] filename w l 
}

set out
set terminal wxt enhanced
