# This script downloads and unzips Hydrogen chains from HamLib.
# This is required to run before running any scripts which use HamLib Hydrogen chains
# (e.g., scripts/hchains/hchains.py).
url="https://portal.nersc.gov/cfs/m888/dcamps/hamlib/chemistry/electronic/hydrogen_data/"
sizes="16 20 24 28 32"

for size in $sizes
do
  fname="${url}H${size}_linear/ES_H${size}_linear_R0.5_sto-6g"
  if [ $size == 16 ]; then
    fname="${fname}_ham"
  fi
  echo "Downloading" "${fname}.zip"
  curl -s "${fname}.zip" -L -o "$size.zip" > /dev/null
  unzip -o -qq "${size}.zip"
  rm "${size}.zip"
done
