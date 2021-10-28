#!/bin/sh

INPUTBEDPE=${1}

awk 'BEGIN {OFS="\t"} { print $1,$2,$3,"amp"NR"fwd",".","+"; print $4,$5,$6,"amp"NR"rev",".","-"}' $INPUTBEDPE | sort -k1,1 -k2n,3n
