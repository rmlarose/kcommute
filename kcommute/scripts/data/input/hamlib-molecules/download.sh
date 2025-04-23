# This script downloads and unzips molecules from HamLib.
# This is required to run before running any scripts which use HamLib molecules
# (e.g., scripts/hamlib-molecules/molecules.ipynb).
url="https://portal.nersc.gov/cfs/m888/dcamps/hamlib/chemistry/electronic/standard/"
molecules="O2 B2 BeH BH CH HF C2 OH N2 Li2 NaLi"

for molecule in $molecules
do
  echo "Downloading $url$molecule.zip"
  curl -s "$url$molecule.zip" -L -o "$molecule.zip" > /dev/null
  unzip -o -qq "$molecule.zip"
  rm "$molecule.zip"
done
