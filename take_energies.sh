#!/bin/bash
## This script takes the energy calculated by AMBER single point

number_line=$(grep -n "ENERGY" $1 | cut -f1 -d: | head -n 1)
energy_line=$((number_line+1))

awk -v "line=$energy_line" 'NR==line' $1 | awk '{print $2}'

exit 0
