set title 'non-steady advection eq'

set term gif animate optimize delay 5 size 700, 480
set output 'fdm_2dim_nonste_advec.gif'

set xl 'x'
set yl 'y'
set zl 'u'

XRANGE = 1.0
YRANGE = 1.0
LOWCBR = -1.0
UPPCBR = 1.0

set grid
set isosamples 50
set palette rgbformula 22, 13, -31

do for [i = 0:100]{
    filename = sprintf('data_fdm_2dim_nonste_advec_%d.txt', i)
    plot_title = sprintf('time number = %d', i)
    set title plot_title
    splot [0:XRANGE] [0:YRANGE] [LOWCBR:UPPCBR] filename u 1:2:3 palette with pm3d
}

set out
set terminal wxt enhanced
