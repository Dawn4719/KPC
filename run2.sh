#!/bin/bash

#cmake --build build

timeout 3600 ./build/NewP2 6 1 ../Dataset/KPC_Dataset/SDXL6.csv 1 2
timeout 3600 ./build/NewP2 6 1 ../Dataset/KPC_Dataset/SDXL6.csv 1 3
timeout 3600 ./build/NewP2 6 1 ../Dataset/KPC_Dataset/SDXL6.csv 2 3
timeout 3600 ./build/NewP2 6 1 ../Dataset/KPC_Dataset/SDXL6.csv 2 4


#for data in congressZJ4
#do
#  for i in {5..9}
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP1 4 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP1 4 $i ../Dataset/KPC_Dataset/${data}.csv || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done

#for data in CTXL5 youtubeXL5
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NULL 5 8 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NULL 5 8 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL6 youtubeXL6
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NULL 6 9 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NULL 6 9 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
## R1
#
#for data in CTZJ4 youtubeZJ4
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/R1 4 7 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/R1 4 7 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL5 youtubeXL5
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/R1 5 8 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/R1 5 8 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL6 youtubeXL6
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/R1 6 9 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/R1 6 9 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
## R1R2
#
#for data in CTZJ4 youtubeZJ4
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/R1R2 4 7 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/R1R2 4 7 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL5 youtubeXL5
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/R1R2 5 8 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/R1R2 5 8 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL6 youtubeXL6
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/R1R2 6 9 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/R1R2 6 9 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done

# NOIDORDER

#for data in CTZJ4 youtubeZJ4
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NOIDORDER 4 7 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NOIDORDER 4 7 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL5 youtubeXL5
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NOIDORDER 5 8 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NOIDORDER 5 8 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL6 youtubeXL6
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NOIDORDER 6 9 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NOIDORDER 6 9 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done

# NOALL

#for data in CTZJ4 youtubeZJ4
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NOALL 4 7 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NOALL 4 7 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL5 youtubeXL5
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NOALL 5 8 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NOALL 5 8 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done
#
#for data in CTXL6 youtubeXL6
#do
#  da=1
#  time3=$(date "+%Y-%m-%d %H:%M:%S")
#  echo $time3
#  echo ./build/NOALL 6 9 ../Dataset/KPC_Dataset/${data}.csv
#  timeout 3600 ./build/NOALL 6 9 ../Dataset/KPC_Dataset/${data}.csv || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#done