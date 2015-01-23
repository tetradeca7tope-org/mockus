# shell script for converting eps to pdf

for f in *.eps
do
  echo $f
  ps2pdf -dEPSCrop $f
done
