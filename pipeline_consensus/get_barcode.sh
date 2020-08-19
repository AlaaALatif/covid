#!/bin/bash

f=$1
echo -e $f"\t"$(awk '$3 > 5' $f | cut -f 1 | sed "s/_rc//g" | sort | uniq | grep -v "unknown\|forward")"\t"$(awk '$3 > 5' $f | cut -f 2 | sed "s/_rc//g" | sort | uniq | grep -v "unknown\|reverse")
