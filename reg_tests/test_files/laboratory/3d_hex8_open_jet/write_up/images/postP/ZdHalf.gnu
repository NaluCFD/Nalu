      set   autoscale                        # scale axes automatically
      set xrange [0.0:5.0]
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set xlabel "r/r_h"
      set ylabel "Uz/Uz^o"
      set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 16 
      set output "ZdHalf.ps"
      set size ratio 0.5

      set key on right top
#      set key width -3

      set style line 1 lt 2
      set style line 2 lt 4

      normalizeUFive(x) = x/15.037
      normalizeRFive(x) = x/0.65861

      normalizeUTen(x) = x/10.337
      normalizeRTen(x) = x/0.9860591

      normalizeUTwenty(x) = x/6.2604
      normalizeRTwenty(x) = x/1.612608

      plot   "ZdFive.dat" using (normalizeRFive($6)):(normalizeUFive($5)) every 5 title 'z/d_o = 5' with linespoints ls 2,\
             "ZdTen.dat" using (normalizeRTen($6)):(normalizeUTen($5)) every 5 title 'z/d_o = 10' with linespoints ls 6,\
             "ZdTwenty.dat" using (normalizeRTwenty($6)):(normalizeUTwenty($5)) every 5 title 'z/d_o = 20' with linespoints ls 8