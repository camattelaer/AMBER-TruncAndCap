#!/bin/bash

### 1) generate crude truncated receptor for residue selection
# MODIFY: LigandRes should have the residue number of the ligand in the PDB file
# MODIFY: change the location of the input files (parm, reference and trajin lines)
LigandRes=1

cat <<eof> trunc.in
parm complex_solv.prmtop
reference complex_solv.pdb
trajin complex_solv.pdb
strip !(:$LigandRes<:8)
strip :$LigandRes
strip :WAT
strip :Na+
strip :Cl-
trajout trunc.pdb
run
quit
eof

cpptraj -i trunc.in -o trunc.out

### 2) create the list of protein residues present in the truncated version

grep 'ATOM' trunc.pdb | awk '{print $5}' | uniq | awk 'NR==1{first=$1;last=$1;next} $1 == last+1 {last=$1;next} {print first,last;first=$1;last=first} END{print first,last}' > resNo.txt

# grep : remove CRYST or TER lines
# awk '{print $5}' extracts residue numbers
# uniq : collaps redundant residue numbers (due to one residue number per atom)
# awk : 
#	NR==1{first=$1;last=$1;next}
#    	On the first line, initialize the variables first and last and skip to next line.

#	$1 == last+1 {last=$1;next}
#	If this line continues in the sequence from the last, update last and jump to the next line.

#	print first,last;first=$1;last=first
#	If we get here, we have a break in the sequence. Print out the range for the last sequence and reinitialize the variables for a new sequence.
#	END{print first,last}

# 	After we get to the end of the file, print the final sequence.

### 3) we need to patch the residue list to add extra residues which will then be mutated into Me Caps

## 3.1) clean initial residue list (concatenate sequences where only 1 residue separates the peptides, eg 183 and 185)

# calculate differences between ending and starting residues
awk 'BEGIN { OFS = FS } NR == 1 { last = $2; print;  next }{ $3 = $1 - last; last = $2  }1' resNo.txt > resNo_diff.txt

# get lines where difference between end and start < 2 
awk '{if ($3 <= 2 && $3 != "") { print 1 } else { print 0 } }' resNo_diff.txt > lines_to_patch.txt

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
	awk -v var=$((lineNo)) 'NR==var { print $0 }' resNo.txt >> resNo_final.txt
elif [ $bool == 1 ] ; then
	startRes_curr=$(awk -v var=$((lineNo)) 'NR==var { print $1 }' resNo.txt)
	bool_curr=1
	while [ $bool_curr == 1 ] ; do
		lineNo=$((lineNo+1))
		bool_curr=$(awk -v var=$((lineNo)) 'NR==var { print $1 }' lines_to_patch.txt)
	done
	endRes_curr=$(awk -v var=$((lineNo-1)) 'NR==var { print $2 }' resNo.txt)
	echo "$startRes_curr $endRes_curr" >> resNo_final.txt
	i=$lineNo
fi

done

tail -n 1 resNo.txt >> resNo_final.txt

## 3.2) create STRIP string for residues within distance threshold

while IFS=" " read -r startRes endRes ; do
if [ $startRes == $endRes ] ; then
Res="$startRes,$Res"
else
Res="$startRes-$endRes,$Res"
fi
done < <(tac resNo_final.txt)

Res=$(echo $Res | sed s/,$// )

## 3.3) create STRIP string for Nterm residues

startRes=$(awk '{print $1}' resNo_final.txt | tac)

for res in $startRes ; do

NtermRes=$((res-1))
AceRes="$NtermRes,$AceRes"

done

AceRes=$(echo $AceRes | sed s/,$// )

## 3.4) create STRIP string for Cterm residues

endRes=$(awk '{print $2}' resNo_final.txt | tac)

for res in $endRes ; do

CtermRes=$((res+1))
NmeRes="$CtermRes,$NmeRes"

done

NmeRes=$(echo $NmeRes | sed s/,$// )

## 3.5) combine the necessary STRIP strings for a cpptraj input file and generate extended truncated pdb

# MODIFY: change (paths for) input files
cat <<eof> trunc_cap.in
parm complex_solv.prmtop
reference complex_solv.pdb
trajin complex_solv.pdb
strip !(:$AceRes@C,O|:$NmeRes@N,H|:$Res|:$LigandRes)
strip :WAT
strip :Na+
strip :Cl-
trajout trunc_cap.pdb
run
quit
eof

cpptraj -i trunc_cap.in -o trunc_cap.out

## 3.6) rename Nterm residues to ACE

for res in $startRes ; do

NtermRes=$((res-1))

oldResName=$(awk -v var=$NtermRes '$5 == var {print $4}' trunc_cap.pdb | tail -n 1)

if [ $NtermRes -lt 100 ] ; then
sed -i "s/$oldResName    $NtermRes/ACE    $NtermRes/g" trunc_cap.pdb
elif [ $NtermRes -lt 1000 ] ; then
sed -i "s/$oldResName   $NtermRes/ACE   $NtermRes/g" trunc_cap.pdb
else
sed -i "s/$oldResName  $NtermRes/ACE  $NtermRes/g" trunc_cap.pdb
fi

done

## 3.7) rename Cterm residues to NME

for res in $endRes ; do

CtermRes=$((res+1))

oldResName=$(awk -v var=$CtermRes '$5 == var {print $4}' trunc_cap.pdb | tail -n 1)

if [ $CtermRes -lt 100 ] ; then
sed -i "s/$oldResName    $CtermRes/NME    $CtermRes/g" trunc_cap.pdb
elif [ $CtermRes -lt 1000 ] ; then
sed -i "s/$oldResName   $CtermRes/NME   $CtermRes/g" trunc_cap.pdb
else
sed -i "s/$oldResName  $CtermRes/NME  $CtermRes/g" trunc_cap.pdb
fi
done

### 4) remove unnecessary files

rm resNo.txt resNo_diff.txt resNo_final.txt lines_to_patch.txt trunc.in trunc.out trunc.pdb
