#!/bin/bash


##cmake --build build
##./build/R2
##./build/R3
#
#for data in congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
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
#
#for data in congressXL5 CTXL5 SDXL5 youtubeXL5
#do
#  for i in {6..10}
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP1 5 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP1 5 $i ../Dataset/KPC_Dataset/${data}.csv || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done



#for data in congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv 1 2 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv 1 3 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv 2 3 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 4 $i ../Dataset/KPC_Dataset/${data}.csv 2 4 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL5 CTXL5 SDXL5 youtubeXL5
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv 1 2 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL5 CTXL5 SDXL5 youtubeXL5
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv 1 3 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL5 CTXL5 SDXL5 youtubeXL5
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv 2 3 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL5 CTXL5 SDXL5 youtubeXL5
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 5 $i ../Dataset/KPC_Dataset/${data}.csv 2 4 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done



#for data in congressXL6 CTXL6 SDXL6 youtubeXL6
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv 1 2 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL6 CTXL6 SDXL6 youtubeXL6
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv 1 3 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL6 CTXL6 SDXL6 youtubeXL6
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv 2 3 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressXL6 CTXL6 SDXL6 youtubeXL6
#do
#  for i in 1
#  do
#    da=1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    echo ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv
#    timeout 3600 ./build/NewP2 6 $i ../Dataset/KPC_Dataset/${data}.csv 2 4 || da=2
#    if [ $da -eq 2 ]; then
#      break
#    fi
#  done
#done
#
#for data in congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
#do
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#    echo $time3
#    ./build/NewP3 4 8 ../Dataset/KPC_Dataset/${data}.csv 1 1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#        echo $time3
#    ./build/NewP3 4 8 ../Dataset/KPC_Dataset/${data}.csv 2 2
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#        echo $time3
#    ./build/NewP3 4 8 ../Dataset/KPC_Dataset/${data}.csv 3 3
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#        echo $time3
#    ./build/NewP3 4 8 ../Dataset/KPC_Dataset/${data}.csv 4 4
#done
#
#for data in congressXL5 CTXL5 SDXL5 youtubeXL5
#do
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#      echo $time3
#    ./build/NewP3 5 10 ../Dataset/KPC_Dataset/${data}.csv 1 1
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#        echo $time3
#    ./build/NewP3 5 10 ../Dataset/KPC_Dataset/${data}.csv 2 2
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#        echo $time3
#    ./build/NewP3 5 10 ../Dataset/KPC_Dataset/${data}.csv 3 3
#    time3=$(date "+%Y-%m-%d %H:%M:%S")
#        echo $time3
#    ./build/NewP3 5 10 ../Dataset/KPC_Dataset/${data}.csv 4 4
#done

#for i in 1
#do
#  da=1
#  timeout 3600 ./build/base2Problem2 5 10 ../Dataset/KPC_Dataset/congressXL5.csv 1 2 || da=2
#  if [ $da -eq 2 ]; then
#    break
#  fi
#  timeout 3600 ./build/base2Problem2 5 15 ../Dataset/KPC_Dataset/congressXL5.csv 1 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem2 5 15 ../Dataset/KPC_Dataset/congressXL5.csv 2 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem2 5 20 ../Dataset/KPC_Dataset/congressXL5.csv 2 4 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#done
#for i in 1
#do
#  da=1
#  timeout 3600 ./build/base2Problem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 1 1 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 2 2 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 3 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 4 4 || da=2
#  if [ $da -eq 2 ]; then
#        break
#      fi
#done
#
#for i in 1
#do
#  da=1
#  timeout 3600 ./build/mainProblem2 5 10 ../Dataset/KPC_Dataset/congressXL5.csv 1 2 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem2 5 15 ../Dataset/KPC_Dataset/congressXL5.csv 1 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem2 5 15 ../Dataset/KPC_Dataset/congressXL5.csv 2 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem2 5 20 ../Dataset/KPC_Dataset/congressXL5.csv 2 4 || da=2
#  if [ $da -eq 2 ]; then
#        break
#      fi
#done
#
#for i in 1
#do
#  da=1
#  timeout 3600 ./build/mainProblem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 1 1 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 2 2 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 3 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem3 5 50 ../Dataset/KPC_Dataset/congressXL5.csv 4 4 || da=2
#  if [ $da -eq 2 ]; then
#        break
#      fi
#done
#
#
#for i in 1
#do
#  da=1
#  timeout 3600 ./build/base2Problem2 6 12 ../Dataset/KPC_Dataset/congressXL6.csv 1 2 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem2 6 18 ../Dataset/KPC_Dataset/congressXL6.csv 1 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem2 6 18 ../Dataset/KPC_Dataset/congressXL6.csv 2 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/base2Problem2 6 24 ../Dataset/KPC_Dataset/congressXL6.csv 2 4 || da=2
#  if [ $da -eq 2 ]; then
#        break
#      fi
#done
#
#for i in 1
#do
#  da=1
#  timeout 3600 ./build/mainProblem2 6 12 ../Dataset/KPC_Dataset/congressXL6.csv 1 2 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem2 6 18 ../Dataset/KPC_Dataset/congressXL6.csv 1 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem2 6 18 ../Dataset/KPC_Dataset/congressXL6.csv 2 3 || da=2
#  if [ $da -eq 2 ]; then
#      break
#    fi
#  timeout 3600 ./build/mainProblem2 6 24 ../Dataset/KPC_Dataset/congressXL6.csv 2 4 || da=2
#  if [ $da -eq 2 ]; then
#        break
#      fi
#done



for data in congressXL6
do
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/base2Problem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 1 1 || da=2
  if [ $da -eq 2 ]; then
      break
    fi
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/base2Problem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 2 2 || da=2
  if [ $da -eq 2 ]; then
      break
    fi
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/base2Problem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 3 3 || da=2
  if [ $da -eq 2 ]; then
      break
    fi
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/base2Problem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 4 4 || da=2
  if [ $da -eq 2 ]; then
      break
    fi
done

for data in congressXL6
do
  da=1
      time3=$(date "+%Y-%m-%d %H:%M:%S")
    echo $time3
    timeout 3600 ./build/mainProblem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 1 1 || da=2
    if [ $da -eq 2 ]; then
      break
    fi

    time3=$(date "+%Y-%m-%d %H:%M:%S")
    echo $time3
    timeout 3600 ./build/mainProblem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 2 2 || da=2
    if [ $da -eq 2 ]; then
      break
    fi

    time3=$(date "+%Y-%m-%d %H:%M:%S")
    echo $time3
    timeout 3600 ./build/mainProblem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 3 3 || da=2
    if [ $da -eq 2 ]; then
      break
    fi

    time3=$(date "+%Y-%m-%d %H:%M:%S")
    echo $time3
    timeout 3600 ./build/mainProblem3 6 50 ../Dataset/KPC_Dataset/${data}.csv 4 4 || da=2
    if [ $da -eq 2 ]; then
      break
    fi
done

for data in SDXL6
do
  da=1
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/NewP3 6 12 ../Dataset/KPC_Dataset/${data}.csv 2 2 || da=2
  if [ $da -eq 2 ]; then
    break
  fi
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/NewP3 6 12 ../Dataset/KPC_Dataset/${data}.csv 3 3 || da=2
  if [ $da -eq 2 ]; then
    break
  fi
  time3=$(date "+%Y-%m-%d %H:%M:%S")
  echo $time3
  timeout 3600 ./build/NewP3 6 12 ../Dataset/KPC_Dataset/${data}.csv 4 4 || da=2
  if [ $da -eq 2 ]; then
    break
  fi
done