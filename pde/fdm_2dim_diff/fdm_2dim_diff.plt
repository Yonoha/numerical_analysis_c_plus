set title 'non-steady diffusion eq'

set term gif animate optimize delay 10 size 700, 480
set output 'fdm_2dim_nonste_diff.gif'

set xl 'x'
set yl 'y'
set zl 'u'

set cbr [-1:1]

set grid
set isosamples 50
set palette rgbformula 22, 13, -31

XRANGE = 2
YRANGE = 3
LOWCBR = -1
UPPCBR = 1

do for [i = 0:100]{
    filename = sprintf('data_nonste_diff_%d.txt', i)
    plot_title = sprintf('time number = %d', i)
    set title plot_title
    splot [0:XRANGE] [0:YRANGE] [LOWCBR:UPPCBR] filename using 1:2:3 palette with pm3d
}

set out
set terminal wxt enhanced
