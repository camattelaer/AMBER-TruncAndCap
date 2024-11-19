#!/bin/bash
set -x
debug="True"
### 1) generate crude truncated receptor for residue selection
# MODIFY: change the location of the input files (parm, reference and trajin lines)
#LigandRes="<LIGAND_RES>"
complex_original_pdb="complex_solv.pdb"
#Radius=<RADIUS>
residue_file="resNo_list.txt"
remove_single_residues="False"
regex_numeric='^[0-9]+$'

# clean pdb

grep -E '^ATOM|^HETATM' $complex_original_pdb > complex_orig.pdb

## extract cpptraj residue list

cpptraj -p complex_orig.pdb --resmask \* > cpptraj_reslist.txt

if [[ $LigandRes =~ $regex_numeric ]] ; then
	LigandResName=$(awk -v ligandres_awk=$LigandRes '$6==ligandres_awk {print $2}' cpptraj_reslist.txt)
else
	LigandResName=$LigandRes
fi

## test if pdb is configured with or without chainid -> we use this when truncating from predefined list

column5=$(head -1 complex_orig.pdb | awk '{print $5}')

if [[ $column5 =~ $regex_numeric ]] ; then
  chainID=' '
  starting_residue=$(head -1 complex_orig.pdb | awk '{print $5}')
  last_residue=$(awk '{print $4 " " $6}' complex_orig.pdb  | sort -u -n -k2,2 | grep -E '^GLY|^ALA|^VAL|^LEU|^ILE|^THR|^SER|^MET|^CYS|^PRO|^PHE|^TYR|^TRP|^HIS|^LYS|^ARG|^ASP|^GLU|^ASN|^GLN' | awk 'END {print $2}')
else
  chainID=$column5
  starting_residue=$(head -1 complex_orig.pdb | awk '{print $6}')
  last_residue=$(awk '{print $4 " "  $5 " " $6}' complex_orig.pdb  | sort -u -n -k3,3 | grep -E '^GLY|^ALA|^VAL|^LEU|^ILE|^THR|^SER|^MET|^CYS|^PRO|^PHE|^TYR|^TRP|^HIS|^LYS|^ARG|^ASP|^GLU|^ASN|^GLN' | awk 'END {print $3}')
fi

# calculate differences between ending and starting residues
awk 'BEGIN { OFS = FS } NR == 1 { last = $2; print;  next }{ $3 = $1 - last; last = $2  }1' ${residue_file} > resNo_diff.txt

# get lines where difference between end and start < 2 
awk '{if ($3 > 0 && $3 <= 2 && $3 != "") { print 1 } else { print 0 } }' resNo_diff.txt > lines_to_patch.txt

# patch residue list so gaps smaller than 2 residues are bridged by the deleted residue
# this should also work with subsequent gaps,  eg 183 185 and 187 in the example 
if [ -f resNo_final.txt ] ; then
rm resNo_final.txt
fi

lines=$(cat lines_to_patch.txt | wc -l )

for ((i=2 ; i<=$lines ; i+=1)); do

lineNo=$((i-1))
bool=$(awk -v var=$i 'NR==var { print $1} ' lines_to_patch.txt) 
if [ $bool == 0 ] ; then
	awk -v var=$((lineNo)) 'NR==var { print $0 }' ${residue_file} >> resNo_final.txt
elif [ $bool == 1 ] ; then
	startRes_curr=$(awk -v var=$((lineNo)) 'NR==var { print $1 }' ${residue_file})
	bool_curr=1
	while [ $bool_curr == 1 ] ; do
		lineNo=$((lineNo+1))
		bool_curr=$(awk -v var=$((lineNo)) 'NR==var { print $1 }' lines_to_patch.txt)
	done
	endRes_curr=$(awk -v var=$((lineNo-1)) 'NR==var { print $2 }' ${residue_file})
	echo "$startRes_curr $endRes_curr" >> resNo_final.txt
	i=$lineNo
fi

done

# check if last line has been patched
last_resNo_diff=$(tail -n 1 resNo_diff.txt | awk '{print $3}')

if [[ $last_resNo_diff > 2 ]] ; then
tail -n 1 ${residue_file} >> resNo_final.txt
fi

#use final filter to remove single residues
if [ $remove_single_residues == "True" ] ; then
	awk '!($1==$2)'  resNo_final.txt > resNo_final_tmp.txt && mv resNo_final_tmp.txt resNo_final.txt
fi

## 3.2) create STRIP string for residues within distance threshold

while IFS=" " read -r startRes endRes ; do
if [ $startRes == $endRes ] ; then
startRes=$(awk -v ligandres_awk=$startRes '$6==ligandres_awk {print $1}' cpptraj_reslist.txt)
Res="$startRes,$Res"
else
startRes=$(awk -v ligandres_awk=$startRes '$6==ligandres_awk {print $1}' cpptraj_reslist.txt)
endRes=$(awk -v ligandres_awk=$endRes '$6==ligandres_awk {print $1}' cpptraj_reslist.txt)
Res="$startRes-$endRes,$Res"
fi
done < <(tac resNo_final.txt)

Res=$(echo $Res | sed s/,$// )

## 3.5) combine the necessary STRIP strings for a cpptraj input file and generate extended truncated pdb

# MODIFY: change (paths for) input files
cat <<eof> trunc_cap.in
parm complex_orig.pdb
reference complex_orig.pdb
trajin complex_orig.pdb
#strip !(:$Res|:$LigandResName)
strip !(:$Res)
strip :WAT,HOH
strip :Na+
strip :Cl-
trajout trunc_cap.pdb
run
quit
eof

cpptraj -i trunc_cap.in -o trunc_cap.out

## 3.6) generate pymol script to add ACE and NME caps

cat <<eof> pymol-caps.py
from pymol import editor, cmd
cmd.load("trunc_cap.pdb")
# select & remove all non A altlocs
cmd.remove( "not alt ''+A" )
# reset the PDB information
cmd.alter("all", "alt=''")
eof

while IFS=" " read -r startRes endRes ; do
echo "editor.attach_amino_acid(\"resi $startRes and name N\", \"ace\")" >> pymol-caps.py
echo "editor.attach_amino_acid(\"resi $endRes and name C\", \"nme\")" >> pymol-caps.py
done < resNo_final.txt

echo "cmd.set(\"pdb_use_ter_records\",0)" >> pymol-caps.py
echo "cmd.save(\"trunc_cap_pymol.pdb\")" >> pymol-caps.py

source activate CompChem

python pymol-caps.py

### replace atom names for ACE and NME residues compatible with LEaP

sed -i "s|1HH3 ACE| H1  ACE|g" trunc_cap_pymol.pdb
sed -i "s|2HH3 ACE| H2  ACE|g" trunc_cap_pymol.pdb
sed -i "s|3HH3 ACE| H3  ACE|g" trunc_cap_pymol.pdb

sed -i "s| CH3 NME| C   NME|g" trunc_cap_pymol.pdb
sed -i "s|1HH3 NME| H1  NME|g" trunc_cap_pymol.pdb
sed -i "s|2HH3 NME| H2  NME|g" trunc_cap_pymol.pdb
sed -i "s|3HH3 NME| H3  NME|g" trunc_cap_pymol.pdb

#### 4) remove unnecessary files

if [ $debug == "False" ] ; then
 rm complex_orig.pdb cpptraj_reslist.txt trunc.in trunc.out trunc.pdb resNo_diff.txt lines_to_patch.txt trunc_cap.in trunc_cap.out # resNo_10A.txt resNo_diff.txt resNo_final.txt lines_to_patch.txt trunc.in trunc.out trunc.pdb
fi
