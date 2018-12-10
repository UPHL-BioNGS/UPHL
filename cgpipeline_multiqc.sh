#!/bin/bash

out=$1

grep avgReadLength $out/cg-pipeline/*.out.txt | sort | uniq | head -n 1 | cut -f 2- -d ':' > $out/cg-pipeline/cg-pipeline-summary.txt
grep -v \"avgReadLength\" $out/cg-pipeline/*.out.txt | cut -f 2- -d ':' | sort | uniq >> $out/cg-pipeline/cg-pipeline-summary.txt
