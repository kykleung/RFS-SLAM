#!/bin/sh
# script to copy the headers to all the source files and header files

for f in ../README; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    sed -i '1,/^$/ d' $f 
    echo "Removing License Header in $f"
  fi
  cat ../License $f > $f.new
  mv $f.new $f
  echo "New license Header copied to $f" 
done  

for f in ../include/*.hpp; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    sed -i '1,/^$/ d' $f 
    echo "Removing License Header in $f"
  fi
  cat ../License $f > $f.new
  mv $f.new $f
  echo "New License Header copied to $f" 
done  

for f in ../src/*.cpp; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    sed -i '1,/^$/ d' $f 
    echo "Removing License Header in $f"
  fi
  cat ../License $f > $f.new
  mv $f.new $f
  echo "New License Header copied to $f" 
done  

for f in ../cfg/*; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    sed -i '1,/^$/ d' $f 
    echo "Removing License Header in $f"
  fi
  cat ../License $f > $f.new
  mv $f.new $f
  echo "New License Header copied to $f" 
done  