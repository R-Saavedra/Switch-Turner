#!/bin/bash

INSTALLATION_PATH="/home/rsaavedra/Photoswitches_others/protocol_reparametrization"

gaff_topology=$1
torsion=$2
shift 2

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fitter)    fit_mode="$2"
		shift 
		;;
        --reference) ref_file="$2"
	       	shift 
		;;
        --mm)        mm_file="$2" 
		shift 
		;;
        --angle)     angle_file="$2" 
		shift 
		;;
    esac
    shift
done


if [ ${fit_mode} == "normal" ]; then
	python ${INSTALLATION_PATH}/fitter.py $mm_file $ref_file $angle_file ${torsion}_results.dat > fit_parameters_${torsion}.gprm
elif [ ${fit_mode} == "hyper" ]; then
	python ${INSTALLATION_PATH}/hyper_fitter.py $mm_file $ref_file $angle_file ${torsion}_results.dat > fit_parameters_${torsion}.gprm
else
	echo "Fitting mode unmatched... performing normal fitting."
	python ${INSTALLATION_PATH}/fitter.py $mm_file $ref_file $angle_file ${torsion}_results.dat > fit_parameters_${torsion}.gprm
fi

## 
last_proper_torsion=$(grep -n "IMPROPER" ${gaff_topology} | cut -d: -f1)
file=$(head -n $last_proper_torsion $gaff_topology)

atom_1=$(sed -n '1p' ${torsion}.txt)
atom_2=$(sed -n '2p' ${torsion}.txt)
atom_3=$(sed -n '3p' ${torsion}.txt)
atom_4=$(sed -n '4p' ${torsion}.txt)

torsion_name_order_1=$(grep "$atom_1-$atom_2-$atom_3-$atom_4" <<< ${file} | head -n 1)
torsion_name_order_2=$(grep "$atom_4-$atom_3-$atom_2-$atom_1" <<< ${file} | head -n 1)
if [[ -n $torsion_name_order_1 && -n $torsion_name_order_2 ]]; then
        echo "More than one naming order for the torsion. Aborting..."
        exit 1
else
        [[ -n $torsion_name_order_1 ]] && { template_line=$torsion_name_order_1 ; torsion_atoms="$atom_1-$atom_2-$atom_3-$atom_4" ; } 
        [[ -n $torsion_name_order_2 ]] && { template_line=$torsion_name_order_2 ; torsion_atoms="$atom_4-$atom_3-$atom_2-$atom_1" ; }
fi

n_torsions=$(awk '{print $2}' <<< ${template_line})

n_terms=$(grep "Terms" fit_parameters_${torsion}.gprm | cut -d ':' -f2 | sed 's/,//; s/\[//g; s/\]//g' | wc -w)
for i in $(seq 1 ${n_terms}); do
	## 
	term=$(grep "Terms" fit_parameters_${torsion}.gprm | cut -d ':' -f2 | sed 's/,//g; s/\[//g; s/\]//g' | awk -v "i=$i" '{print $i}')
	phase=$(grep "Shift" fit_parameters_${torsion}.gprm | cut -d ':' -f2 | sed 's/,//g; s/\[//g; s/\]//g' | awk -v "i=$i" '{print $i}')
	raw_coeff=$(grep "Coefficients" fit_parameters_${torsion}.gprm | cut -d ':' -f2 | sed 's/,//g; s/\[//g; s/\]//g' | awk -v "i=$i" '{print $i}')
	coeff=$(awk -v c="$raw_coeff" -v n="$n_torsions" 'BEGIN {printf "%.8f", c/n}')
	#;> coeff=$(echo "scale=8; $coeff / $n_torsions" | bc)

	if [ $i -eq $n_terms ]; then
		printf "%-15s %2d %.8f %10.3f %5.1f      #%d/%d Fitted %s Functions\n" \
			"${torsion_atoms}" "${n_torsions}" "${coeff}" "${phase}" "${term}" "${i}" "${n_terms}" "${torsion}" >> ${torsion}_refined.amber.gprm
	else
		printf "%-15s %2d %.8f %10.3f %5.1f      #%d/%d Fitted %s Functions\n" \
			"${torsion_atoms}" "${n_torsions}" "${coeff}" "${phase}" "-${term}" "${i}" "${n_terms}" "${torsion}" >> ${torsion}_refined.amber.gprm
	fi
	##
done

torsion_lines=$(grep -n "$torsion_atoms" $gaff_topology | cut -d: -f1)
deleting_lines=$(echo "${torsion_lines}d" | sed 's/ /d;/g') ## All lines except the first one are deleted
sed -i -e "${deleting_lines}" $gaff_topology #?> tmp_top.tmp && { mv tmp_top.tmp ${gaff_topology} ; }

first_line=$(echo "${torsion_lines}" | head -n 1)
sed -i "${first_line}r ${torsion}_refined.amber.gprm" ${gaff_topology}


exit 0
