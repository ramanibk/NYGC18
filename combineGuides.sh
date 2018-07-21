#!/bin/bash
#counter=1
#while [ $counter -lt 23 ]; do
#    python searchGuides.py "../chroms/chr"$counter".fa" $counter &
#    let counter=counter+1
#done

#wait

#counter=1
#while [ $counter -lt 23 ]; do
#    wc -l "../GUIDES/guideschr"$counter".txt" >> counts.txt
#    let counter=counter+1
#done

#python searchGuides.py "../chroms/chrX.fa" "X" &
#python searchGuides.py "../chroms/chrY.fa" "Y" &
#wait
#wc -l "../GUIDES/guideschrX.txt" >> counts.txt
#wc -l "../GUIDES/guideschrY.txt" >> counts.txt

counter=1
while [ $counter -lt 23 ]; do
    tail -n +2 "../chroms/chr"$counter.fa >> "chr"$counter".txt"
    let counter=counter+1
done
tail -n +2 "../chroms/chrX.fa" >> "chrX.txt"
tail -n +2 "../chroms/chrY.fa" >> "chrY.txt"
