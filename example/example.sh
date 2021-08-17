#!/bin/bash
echo "*** Running test analysis ***"
echo
perl ../bcpi.pl --snp snp.lst --f1 GF_1.fq.gz --f2 GF_2.fq.gz --m1 GM_1.fq.gz --m2 GM_2.fq.gz --c1 GC_1.fq.gz --c2 GC_2.fq.gz --out ./test_output --w 6
