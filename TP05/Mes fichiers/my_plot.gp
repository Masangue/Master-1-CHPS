data = ARG1
name = ARG2
xlbl = ARG3
ylbl = ARG4
ttl = ARG5

# call gnuplot -c 'data' 'methode'

in = sprintf("saved_data/%s_%s.dat",data, name)
out = sprintf("plots/plot_%s_%s.jpg",data, name)

set terminal jpeg
set output out
print name
set title ttl
set xlabel xlbl
set ylabel ylbl

set style data lines
plot in with lines 

ttl = sprintf("%s (log)", ttl)
set title ttl
set logscale y
out = sprintf("plots/plot_%s_%s_log.jpg",data, name)
set output out
plot in with lines 