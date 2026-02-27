set term png
set output "../plots/Temperature.png"

set xlabel "Time (r.u.)"
set ylabel "Temperature (r.u.)"

unset key
t_ref=1.8781302170283807
plot '../out/energies_T.dat' u 1:5 w lp t 'Temp',  t_ref t 'T_ref'

set output '../plots/Energies.png'

set ylabel 'Energies (r.u.)'
set key outside top horizontal box center
dades='../out/energies_T.dat'

plot dades u 1:2 t 'Upot' w lp, dades u 1:3 t 'Kinetic' w lp, dades u 1:4 t 'Total' w lp

set output '../plots/rdf.png'

set xlabel 'r (r.u.)'
set ylabel 'g(r)'
unset key
dades='../out/gr_ArAr.dat'

plot dades u 1:2 w lp
