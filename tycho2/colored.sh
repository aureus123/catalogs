
echo "Add annotations with colored attributes"

mv cross_tyc2_south.csv cross_tyc2_south.bak
cp cross_tyc2_south_plain.csv cross_tyc2_south.csv
echo "  1) Adding plain stars"
plotann.py $1.wcs $1.png $1.tmp1.png --no-grid --no-const --no-ngc --no-bright --tycho2cat=tycho2.kd --tcolor green
cp cross_tyc2_south_dpl.csv cross_tyc2_south.csv
echo "  2) Adding double stars"
plotann.py $1.wcs $1.tmp1.png $1.tmp2.png --no-grid --no-const --no-ngc --no-bright --tycho2cat=tycho2.kd --tcolor magenta
cp cross_tyc2_south_color.csv cross_tyc2_south.csv
echo "  3) Adding color stars"
plotann.py $1.wcs $1.tmp2.png $1.ann.png --no-grid --no-const --no-ngc --no-bright --tycho2cat=tycho2.kd --tcolor red
mv cross_tyc2_south.bak cross_tyc2_south.csv
rm -f $1.tmp1.png $1.tmp2.png
