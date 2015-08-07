#! /bin/bash
# Generate video from images
# Keith Leung 2013
 
if [ $# -ne 1 ] 
then
    echo "Usage: generateVideo.sh DATA_DIR"
    echo "There must be a folder in DATA_DIR named animation"
    echo "Files in DATA_DIR/animation must be numbered as %06d.png"
    echo "avconv must be installed for video rendering"
    exit
fi
echo "Improving image contrast"
for f in `ls $1/animation/??????.png`
do
    echo "Processing: $f"
    fShort="${f:${#f}-10:${#f}}"
    convert $1/animation/$fShort -contrast $1/animation/_$fShort
done

echo "Generating video from files in: $1"
avconv -r 20 -i $1/animation/_%06d.png -pass 1 -b 30M -f mp4 -an -y /dev/null
avconv -r 20 -i $1/animation/_%06d.png -pass 2 -b 30M -y $1/video.mp4
rm $1/animation/_*;


 