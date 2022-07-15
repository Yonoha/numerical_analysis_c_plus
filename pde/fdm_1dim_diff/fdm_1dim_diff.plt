set title '1dim Burgers eq by fdm'

set term gif animate optimize delay 10 size 700, 480
set output 'fdm_1dim_nonste_diff.gif'

set xl 'x'
set yl 'u'

set grid

XRANGE = 1.0
URANGE = 1.0

do for [i = 0:100]{
    filename = sprintf('data_fdm_1dim_nonste_diff_%d.txt', i)
    plot_title = sprintf('time number = %d', i)
    set title plot_title
    plot [0:XRANGE] [0:URANGE] filename w l 
}

set out
set terminal wxt enhanced
