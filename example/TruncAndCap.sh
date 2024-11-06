#!/bin/bash
set -x
debug="false"
### 1) generate crude truncated receptor for residue selection
# MODIFY: change the location of the input files (parm, reference and trajin lines)
LigandRes="1"
complex_original_pdb="3uo4.pdb"
Radius=12
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

# generate truncated protein


# generate truncated protein

cat <<eof> trunc.in
parm $complex_original_pdb
reference $complex_original_pdb
trajin $complex_original_pdb
strip !(:$LigandResName<:$Radius)
strip :$LigandResName
strip :WAT,HOH
strip :Na+
strip :Cl-
trajout trunc.pdb
run
quit
eof

cpptraj -i trunc.in -o trunc.out

### 2) create the list of protein residues present in the truncated version

if [[ $chainID == ' ' ]] ; then
	grep 'ATOM' trunc.pdb | awk '{print $5}' | uniq | awk 'NR==1{first=$1;last=$1;next} $1 == last+1 {last=$1;next} {print first,last;first=$1;last=first} END{print first,last}' > resNo.txt
else
	grep 'ATOM' trunc.pdb | awk '{print $6}' | uniq | awk 'NR==1{first=$1;last=$1;next} $1 == last+1 {last=$1;next} {print first,last;first=$1;last=first} END{print first,last}' > resNo.txt
fi
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

# check if last line has been patched
last_resNo_diff=$(tail -n 1 resNo_diff.txt | awk '{print $3}')

if [[ $last_resNo_diff > 2 ]] ; then
tail -n 1 resNo.txt >> resNo_final.txt
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

## 3.3) create STRIP string for Nterm residues

startRes=$(awk '{print $1}' resNo_final.txt | tac)

for res in $startRes ; do

NtermRes=$((res-1))
NtermRes=$(awk -v ligandres_awk=$NtermRes '$6==ligandres_awk {print $1}' cpptraj_reslist.txt)
# do a check for residue name, if startRes is ACE, there is no additional NtermRes to be included
res_name=$(awk -v res_awk=$res '$6==res_awk {print $2}' cpptraj_reslist.txt)
if [[ $res_name == "ACE" ]] ; then
AceRes="$AceRes"
else
AceRes="$NtermRes,$AceRes"
fi

done

AceRes=$(echo $AceRes | sed s/,$// )

## 3.4) create STRIP string for Cterm residues

endRes=$(awk '{print $2}' resNo_final.txt | tac)

for res in $endRes ; do

CtermRes=$((res+1))
CtermRes=$(awk -v ligandres_awk=$CtermRes '$6==ligandres_awk {print $1}' cpptraj_reslist.txt)
# do a check for residue name, if endRes is NME, there is no additional CtermRes to be included
res_name=$(awk -v res_awk=$res '$6==res_awk {print $2}' cpptraj_reslist.txt)
if [[ $res_name == "NME" ]] ; then
NmeRes="$NmeRes"
else
NmeRes="$CtermRes,$NmeRes"
fi

done

NmeRes=$(echo $NmeRes | sed s/,$// )

## 3.5) combine the necessary STRIP strings for a cpptraj input file and generate extended truncated pdb

# MODIFY: change (paths for) input files
cat <<eof> trunc_cap.in
parm $complex_original_pdb
reference $complex_original_pdb
trajin $complex_original_pdb
strip !(:$AceRes@C,O,CA|:$NmeRes@N,H,CA|:$Res|:$LigandResName)
strip :WAT,HOH
strip :Na+
strip :Cl-
trajout trunc_cap.pdb
run
quit
eof

cpptraj -i trunc_cap.in -o trunc_cap.out

## 3.6) rename Nterm residues to ACE

AceRes_list=$(echo "${AceRes//,/ }")


for NtermRes in $AceRes_list ; do

#NtermRes=$((res-1))

if [[ -z "${chainID// /}" ]] ; then
oldResName=$(awk -v var=$NtermRes '$5 == var {print $4}' trunc_cap.pdb | tail -n 1)
else
oldResName=$(awk -v var=$NtermRes '$6 == var {print $4}' trunc_cap.pdb | tail -n 1)
fi

if [ $NtermRes -lt 10 ] ; then
sed -i "s/$oldResName $chainID   $NtermRes/ACE $chainID   $NtermRes/g" trunc_cap.pdb
elif [ $NtermRes -lt 100 ] ; then
sed -i "s/$oldResName $chainID  $NtermRes/ACE $chainID  $NtermRes/g" trunc_cap.pdb
elif [ $NtermRes -lt 1000 ] ; then
sed -i "s/$oldResName $chainID $NtermRes/ACE $chainID $NtermRes/g" trunc_cap.pdb
else
sed -i "s/$oldResName ${chainID}$NtermRes/ACE ${chainID}$NtermRes/g" trunc_cap.pdb
fi


done

## 3.7) rename Cterm residues to NME

NmeRes_list=$(echo "${NmeRes//,/ }")

for CtermRes in $NmeRes_list ; do

#CtermRes=$((res+1))

if [[ -z "${chainID// /}" ]] ; then
oldResName=$(awk -v var=$CtermRes '$5 == var {print $4}' trunc_cap.pdb | tail -n 1)
else
oldResName=$(awk -v var=$CtermRes '$6 == var {print $4}' trunc_cap.pdb | tail -n 1)
fi

if [ $CtermRes -lt 10 ] ; then
sed -i "s/$oldResName $chainID   $CtermRes/NME $chainID   $CtermRes/g" trunc_cap.pdb
elif [ $CtermRes -lt 100 ] ; then
sed -i "s/$oldResName $chainID  $CtermRes/NME $chainID  $CtermRes/g" trunc_cap.pdb
elif [ $CtermRes -lt 1000 ] ; then
sed -i "s/$oldResName $chainID $CtermRes/NME $chainID $CtermRes/g" trunc_cap.pdb
else
sed -i "s/$oldResName ${chainID}$CtermRes/NME ${chainID}$CtermRes/g" trunc_cap.pdb
fi
done

### 3.8) rename CA for ACE and NME residues

for res in ${AceRes//,/ } ; do
if [ $res -lt 10 ] ; then
sed -i "s/CA  ACE     $res/CH3 ACE     $res/g" trunc_cap.pdb
elif [ $res -lt 100 ] ; then
sed -i "s/CA  ACE    $res/CH3 ACE    $res/g" trunc_cap.pdb
elif [ $res -lt 1000 ] ; then
sed -i "s/CA  ACE   $res/CH3 ACE   $res/g" trunc_cap.pdb
else
sed -i "s/CA  ACE  $res/CH3 ACE  $res/g" trunc_cap.pdb
fi
done

for res in ${NmeRes//,/ } ; do
if [ $res -lt 10 ] ; then
sed -i "s/CA  NME     $res/C   NME     $res/g" trunc_cap.pdb
elif [ $res -lt 100 ] ; then
sed -i "s/CA  NME    $res/C   NME    $res/g" trunc_cap.pdb
elif [ $res -lt 1000 ] ; then
sed -i "s/CA  NME   $res/C   NME   $res/g" trunc_cap.pdb
else
sed -i "s/CA  NME  $res/C   NME  $res/g" trunc_cap.pdb
fi
done

#### 4) remove unnecessary files

if [ $debug == "False" ] ; then
 rm complex_orig.pdb cpptraj_reslist.txt trunc.in trunc.out trunc.pdb resNo.txt resNo_diff.txt lines_to_patch.txt trunc_cap.in trunc_cap.out # resNo.txt resNo_diff.txt resNo_final.txt lines_to_patch.txt trunc.in trunc.out trunc.pdb
fi
