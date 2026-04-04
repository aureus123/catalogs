
echo "Process a FITS generated from DWARF 3 sensor and estimate variable star"
echo "Usage: estimate.sh fitsfile varstar controlstar firstseqstar secondseqstar controlvmag firstseqvmag secondseqvmag"
echo "Examples:"
echo " estimate.sh YCen 292721 263030 263043 263015 8.15 7.68 8.75"
echo " estimate.sh ROct 376550 376400 376362 376602 8.71 7.95 9.41"
echo

echo "*******************"
echo "* FITS extraction *"
echo "*******************"
echo

# Extract channel from FITS file (channel 1 = Green)
../wcstools-3.9.7/bin/imextract -v -o $1-TG 1 $1.fits

# Solve-plate green channel
solve-field --new-fits none --match none --solved none --rdls none --corr none --scale-units arcsecperpix --scale-low 2.76 --scale-high 2.79 $1-TG.fits --overwrite
tablist "$1-TG.axy[col flux;background]" >Flux.txt
wcs-xy2rd -w $1-TG.wcs -i $1-TG.axy -o $1-TG.radec
tablist "$1-TG.radec[col ra;dec]" >Coord.txt

echo "******************"
echo "* Binning method *"
echo "******************"
echo

cd ..

# Cross-match with PPM catalog by binning method
./cross_txt --txt estim

echo "****************************"
echo "* Variable star estimation *"
echo "****************************"
echo

# Estimation of variable star from PPM ensemble
./cross_txt --var estim $2

# Same use, but with control star
./cross_txt --var estim $2 $3

# Same use, but with control star and a sequence star
./cross_txt --var estim $2 $3 $4

# Same use, but with control star and interpolate two sequence stars
./cross_txt --var estim $2 $3 $4 $5

# Same use, but with custom magnitudes
./cross_txt --var estim $2 $3 $4 $5 $6 $7 $8

cd estim
