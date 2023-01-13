data = ARG1
name = ARG2
datab = ARG3
nameb = ARG4
xlbl = ARG5
ylbl = ARG6

# call gnuplot -c 'data' 'methode'

in = sprintf("saved_data/diff_%s_%s_and_%s_%s.dat",data, name, datab, nameb)
out = sprintf("plots/plot_compare_%s_%s_and_%s_%s.jpg",data, name, datab, nameb)

set terminal jpeg
set output out
print name
ttl = sprintf("Différence entre %s de %s et %s de %s",data, name, datab, nameb)
set title ttl
set xlabel xlbl
set ylabel ylbl

set style data lines
plot in u ($1-$2)
