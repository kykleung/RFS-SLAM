#! /bin/bash
# Stich images together for movie
# Keith Leung 2013

mkdir combined
for f in `ls ??????.png`
do
    echo "Progress: ${f:0:6}"
    f1="${f:0:6}_cp.png"
    f2="${f:0:6}_cp2.png"
    fout="combined/${f:0:6}_combined.png"
    convert $f1 $f2 $f +append $fout
done

avconv -r 20 -i combined/%06d_combined.png -b 30M compare.mp4
rm -rf combined


