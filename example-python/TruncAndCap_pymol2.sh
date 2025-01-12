#!/bin/bash

ligand_name="X01"
radius1=12
pdb="complex_solv.pdb"
## run truncation

cat <<eof> trunc-and-cap.pymol
from pymol import editor,cmd,stored
import numpy
cmd.load("${pdb}")
# select & remove all non A altlocs
cmd.remove( "not alt ''+A" )
# reset the PDB information
cmd.alter("all", "alt=''")
# select zone to remove
cmd.select("${radius1}A_model","(byres (all within ${radius1} of resn ${ligand_name})) and (name ca)")
#
stored.list=[]
cmd.iterate("${radius1}A_model", "stored.list.append(resi)")
## residue numbers are returned as strings, they need to be converted to integers
list_of_residues=[int(x) for x in stored.list]
# now we determine if sequences are divides by 2 residues (too short for both Me caps)
list1=numpy.array(list_of_residues[:-1])
list2=numpy.array(list_of_residues[1:])
indices=numpy.where(list2 == list1+2)[0]
for index in indices:
    a=(list2[index]+list1[index])/2
    list_of_residues.append(int(a.astype(numpy.int64)))
list_of_residues.sort()
# select residues in the new residue list
cmd.select("${radius1}A_model", "resi " + ",".join(str(res) for res in list_of_residues))
# remove residues outside of radius
cmd.remove("not(${radius1}A_model)")
# identify consecutive spans
ranges=sum((list(t) for t in zip(list_of_residues, list_of_residues[1:]) if t[0]+1 != t[1]), [])
starts_and_ends=list_of_residues[0:1] + ranges + list_of_residues[-1:]
starts=starts_and_ends[::2]
ends=starts_and_ends[1::2]
## write ACE residues
for i in starts :
    editor.attach_amino_acid("resi " + str(i) + " and name N", "ace")
## write NME residues
for i in ends :
    editor.attach_amino_acid("resi " + str(i) + " and name C", "nme")
# save pdb
cmd.set("pdb_use_ter_records",0)
cmd.save("receptor_capped_${radius1}A.pdb")
eof

source activate pymol

python trunc-and-cap.pymol

conda deactivate

## replace atom names for ACE and NME residues compatible with LEaP

sed -i "s|1HH3 ACE| H1  ACE|g" receptor_capped_${radius1}A.pdb
sed -i "s|2HH3 ACE| H2  ACE|g" receptor_capped_${radius1}A.pdb
sed -i "s|3HH3 ACE| H3  ACE|g" receptor_capped_${radius1}A.pdb

sed -i "s| CH3 NME| C   NME|g" receptor_capped_${radius1}A.pdb
sed -i "s|1HH3 NME| H1  NME|g" receptor_capped_${radius1}A.pdb
sed -i "s|2HH3 NME| H2  NME|g" receptor_capped_${radius1}A.pdb
sed -i "s|3HH3 NME| H3  NME|g" receptor_capped_${radius1}A.pdb
