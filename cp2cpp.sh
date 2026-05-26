#!/bin/sh

for file in $(find . -path "./*.cu"); do
    echo $file | sed -e "s/\.cu/\.cpp/g" > au.txt
    # echo "$file" | sed -e "s/CI18/CI19/" -e "s/\./tmp/" > au.txt
    fileNew=`cat au.txt`
    cp $file $fileNew
done

rm au.txt
