# Benchmarking test for the C library
First we make datasets:

    for j in 1 2 3 4 5 10 20 30 40 100 200 300 400 1000; do
      i=0
      while [ $i -le $j ]; do
        cat data/test.fa >> /tmp/test_$j.fa
        i=$((i + 1))
      done
    done

The tests are run like this:

    for i in 1 2 3 4 5 10 20 30 40 100 200 300 400 1000; do
      echo -n "$i "
      /usr/bin/time -f "%U" tssv /tmp/test_$i.fa data/library.csv -r /dev/null
    done

Save the files in `old.dat` and `new.dat` respectively.

## Speedup estimation
We just use some shell magic here:

    IFS="
    "
    for i in `paste -d ' ' old.dat new.dat`; do
      echo -n "$(echo $i | cut -f 1 -d ' ') "
      echo 3k $(echo $i | cut -f 2 -d ' ') $(echo $i | cut -f 4 -d ' ') / p | \
        dc
    done > speed.dat

## Plotting
We use `gnuplot` for visualisation:

The raw run times:

    set terminal postscript color
    set output "benchmark.eps"
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2
    set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 2
    set xlabel "dataset size"
    set ylabel "run time"
    plot "old.dat" with linespoints ls 1, "new.dat" with linespoints ls 2

    set terminal postscript color
    set output "speed.eps"
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2
    set xlabel "dataset size"
    set ylabel "speedup"
    plot "speed.dat" with linespoints ls 1

## Results
Here we see the plot for the raw run times:
![][benchmark.png]

And the plot for the speedup:
![][speed.png]
