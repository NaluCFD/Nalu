      set   autoscale                        # scale axes automatically
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set xlabel "z/d_o"
      set ylabel "Q/Q^o"
      set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 16 
      set output "CenterlineData.ps"
      set size ratio 0.5

      set key on right top
#      set key width -3

      set style line 1 lt 2
      set style line 2 lt 4

      normalize(x) = x/20.0

      plot   "CenterlineData.dat" using 8:2 title 'Z/Z^o' with line ls 1,\
             "CenterlineData.dat" using 8:(normalize($5)) title 'Uz/Uz^o' with line ls 2