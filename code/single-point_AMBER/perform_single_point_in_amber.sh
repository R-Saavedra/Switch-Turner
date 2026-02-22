#!/bin/bash

## This script calculates the singlepoint energy of a PDB file
#; This will need as STD-IN-1 the PDB File
#; Right now the forcefield and minimization files are the specific ones from this folder
#; If the system has not the standard AMBER names then these must be manually changed


#! Flags to use:
usage() { 
cat << EOF
Usage:: bash perform_single_point_in_amber.sh [-d STORING_FILE.DAT] COORDINATES.PDB CHARGES.MOL2 PARAMETERS.FRCMOD
	-d <file name>	Store energy in a file. If file already exists append energy 
			at the end. If no file name is provided the flag will be ignored.
	-l <name>	Ligand Name, default is "LIG"

	-h		Displays this message.
EOF
exit 1
}

while getopts ":d:l:" opt; do
	case "${opt}" in 
		d)
			store_file="$OPTARG"
			;;
		l)
			residue_name="$OPTARG"
			;;
		*)
			usage
			exit 0
			;;
	esac
done
shift $((OPTIND - 1))


[[ -z ${residue_name} ]] && residue_name="LIG"
current_name=$(awk 'NR==3 {print $4}' $1)
if [ "${current_name}" != "${residue_name}" ]; then
	sed -i "s/${current_name}/${residue_name}/g" $1
fi
sed -i "/^CONECT/d" $1
sed -i "/^MASTER/d" $1
[[ $(tail -n 1 $1) == "END" ]] || { echo "PDB File does not finnish in _END_ ;; Aborting..." ; exit 1 ; }

${AMBERHOME}/bin/antechamber -i "${1%.*}.pdb" -fi pdb -o "${1%.*}.mol2" -fo mol2 > /dev/null
${AMBERHOME}/bin/antechamber -i "${1%.*}.mol2" -fi mol2 -o "${1%.*}_amber.pdb" -fo pdb > /dev/null

#; Run tleap
tleap -f - > /dev/null 2>&1 << EOF 
source leaprc.gaff2
loadamberparams ${3%.*}.frcmod
LIG = loadmol2 ${1%.*}.mol2
mol = loadpdb ${1%.*}_amber.pdb
set mol box {50.0 50.0 50.0}
saveamberparm mol ${1%.*}.prmtop ${1%.*}.inpcrd
quit 
EOF

[[ -f "leap.log" ]] && mv leap.log $(dirname "${1}")/${1%.*}_leap.log

#; The cpptraj command converts PDB into rst7 ; it uses a forcefield ofc
cpptraj -p ${1%.*}.prmtop -y ${1%.*}.pdb -x ${1%.*}.rst7 > /dev/null 2>&1

#; The pmemd command runs a minimization calculation with 0 steps, so basically runs a single point calculation
pmemd -O -i $(dirname -- "${BASH_SOURCE[0]}")/single_point.in -o ${1%.*}.out -c ${1%.*}.inpcrd -p ${1%.*}.prmtop -r ${1%.*}.rst7 -ref ${1%.*}.inpcrd > /dev/null 2>&1

[[ -f "mdinfo" ]] && rm -- mdinfo

number_line=$(grep -n "ENERGY" ${1%.*}.out | cut -f1 -d: | head -n 1)
energy_line=$((number_line+1))

if [[ -n "${store_file}" ]]; then
	[[ -f "${store_file}" ]] || touch "${store_file}"
	awk -v "line=$energy_line" 'NR==line {print $2}' ${1%.*}.out >> ${store_file}
else
	awk -v "line=$energy_line" 'NR==line {print $2}' ${1%.*}.out 
fi
	


exit 0
