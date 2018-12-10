#!/bin/bash

out=$1
organisms=($(cat $out/mash/*sorted.txt | awk '{ if( $4==0 ) print $1 }' | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "_" $2 }' | sed 's/.-//g' | sort | uniq ))
header="Sample"

for organism in ${organisms[@]}
do
  header=$(echo "$header\t$organism" )
done
echo -e "$header" > $out/mash/mash_results.txt

mash_results=($(ls $out/mash/*sorted.txt | sed 's!.*/!!' | cut -d "_" -f 1 | cut -d '.' -f 1 | sort | uniq ))
for mash_result in ${mash_results[@]}
do
  mash_line=$mash_result
  for organism in ${organisms[@]}
  do
    number=$(cat $out/mash/$mash_result*sorted.txt | awk '{ if( $4==0 ) print $1 }' | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "_" $2 }' | sed 's/.-//g' | sort | uniq -c | grep $organism | awk '{ print $1 }' )
    if [ -z "$number" ]
    then
      number="0"
    fi
    mash_line=$(echo "$mash_line\t$number" )
  done
echo -e "$mash_line" >> $out/mash/mash_results.txt
done
