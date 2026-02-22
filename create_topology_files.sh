#!/bin/bash

#!! This script shall be run within the correct directory!!

base_file=$(basename $1)
input=${base_file%.*}


last_proper_torsion=$(grep -n "IMPROPER" ${input}.frcmod | cut -d: -f1)

file=$(head -n $last_proper_torsion ${input}.frcmod) 	#; As (dihedral) torsions are by default printed before the improper torsions then simply taking the file up to the 
							#; section title "IMPROPER" won't generate any conflict

#; In this version only 2 Torsions are accepted to re-parametrize
## Torsion 1
atom_1=$2
atom_2=$3
atom_3=$4
atom_4=$5

t_order1=$(echo "$file" | grep -n "$atom_1-$atom_2-$atom_3-$atom_4" | cut -d: -f1)
t_order2=$(echo "$file" | grep -n "$atom_4-$atom_3-$atom_2-$atom_1" | cut -d: -f1)
if [[ -n $t_order1 && -n $t_order2 && $t_order1 != $t_order2 ]]; then
	echo "More than one naming order for the torsion. Aborting..."
	exit 1
	#! How to deal with this? 
else
	[[ -z $t_order1 ]] && torsion=$t_order2
	[[ -z $t_order2 ]] && torsion=$t_order1
fi
#; Variable $torsion is now the line within the .frcmod file were the torsion is located
torsion=$(tr '\n' ' ' <<< $torsion)
#; Torsion is now a 1-line string --> Easier to "_play_with_" e.g.: s/ /PATTERN/ ; or #. 

## set torsion(s) to zero ;; Maybe creating an array with the torsion would be a great idea to run the loops ;; maybe not 
if (( $(wc -w <<< $torsion) >= 1 )); then
	deleting_lines=$(echo "${torsion#${torsion%% *} }" | sed 's/ /d;/g') ## All lines except the first one are deleted
	sed -e "${deleting_lines}" ${input}.frcmod > shut-dihe-1_${input}.frcmod
	torsion=${torsion%% *} #; Keep only the first number
	## If there is only 1 number inside the variable well, it is already the first number
elif (( $(wc -w <<< $torsion) != 1 )); then
	echo "Something broke... if only I knew what!" 
	exit 1
fi
#?echo "Torsion at \"$torsion\""
awk -v torsion="$torsion" '
	NR == torsion { $3 = "0.00"}
	{ print }
' OFS='\t' FS='[ \t]+' shut-dihe-1_${input}.frcmod > tmp.t && { sed -e "${torsion}s/#.*$//" tmp.t > shut-dihe-1_${input}.frcmod ; rm tmp.t ; }

shift 4
## Torsion 2
atom_1=$2
atom_2=$3
atom_3=$4
atom_4=$5

t_order1=$(echo "$file" | grep -n "$atom_1-$atom_2-$atom_3-$atom_4" | cut -d: -f1)
t_order2=$(echo "$file" | grep -n "$atom_4-$atom_3-$atom_2-$atom_1" | cut -d: -f1)
if [[ -n $t_order1 && -n $t_order2 && $t_order1 != $t_order2 ]]; then
	echo "More than one naming order for the torsion. Aborting..."
	exit 1
else
	[[ -z $t_order1 ]] && torsion=$t_order2
	[[ -z $t_order2 ]] && torsion=$t_order1
fi
torsion=$(tr '\n' ' ' <<< $torsion)

if (( $(wc -w <<< $torsion) >= 1 )); then
	deleting_lines=$(echo "${torsion#${torsion%% *} }" | sed 's/ /d;/g') 
	sed -e "${deleting_lines}" ${input}.frcmod > shut-dihe-2_${input}.frcmod
	torsion=${torsion%% *}
elif (( $(wc -w <<< $torsion) != 1 )); then
	echo "Something broke... if only I knew what!" 
	exit 1
fi
#?echo "Torsion at \"$torsion\""
awk -v torsion="$torsion" '
	NR == torsion { $3 = "0.00"}
	{ print }
' OFS='\t' FS='[ \t]+' shut-dihe-2_${input}.frcmod > tmp.t && { sed -e "${torsion}s/#.*$//" tmp.t > shut-dihe-2_${input}.frcmod ; rm tmp.t ; }




exit 0
