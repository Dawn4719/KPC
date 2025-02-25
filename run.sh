#!/bin/bash

cmake --build build

for i in {5..9}
do
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  echo ./build/DegenCol 4 $i ../Dataset/dblpZJ4.csv
  timeout 7200 ./build/DegenCol 4 $i ../Dataset/dblpZJ4.csv || echo "!!!! Timeout."
done

for i in {5..9}
do
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  echo ./build/DegenCol 4 $i ../Dataset/youtubeZJ4.csv
  timeout 7200 ./build/DegenCol 4 $i ../Dataset/youtubeZJ4.csv || echo "!!!! Timeout."
done
#
#for i in {6..10}
#do
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/DegenCol 5 $i ../Dataset/dblpXL5.csv
#  timeout 7200 ./build/DegenCol 5 $i ../Dataset/dblpXL5.csv || echo "!!!! Timeout."
#done
#
#for i in {6..10}
#do
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/DegenCol 5 $i ../Dataset/youtubeXL5.csv
#  timeout 7200 ./build/DegenCol 5 $i ../Dataset/youtubeXL5.csv || echo "!!!! Timeout."
#done
#
#for i in {7..11}
#do
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/DegenCol 6 $i ../Dataset/dblpXL6.csv
#  timeout 7200 ./build/DegenCol 6 $i ../Dataset/dblpXL6.csv || echo "!!!! Timeout."
#done
#
#for i in {7..11}
#do
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/DegenCol 6 $i ../Dataset/youtubeXL6.csv
#  timeout 7200 ./build/DegenCol 6 $i ../Dataset/youtubeXL6.csv || echo "!!!! Timeout."
#done