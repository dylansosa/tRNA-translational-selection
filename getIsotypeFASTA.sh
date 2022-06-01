for i in *.highConfidence.fa  ; do

name=$(ls $i | rev | cut -f 4-  -d '.' | rev)
# store assembly name 

cat isotypes.txt | while read -r aa ; do

cat $i | egrep $aa -A1 --no-group-separator > $name.$aa.isotype.fa
# extract sequences per isotype 

spp=$(echo $name | cut -f 3 -d '_')
# store species name

sed -i "s/^>/>$spp./" $name.$aa.isotype.fa ; done ;done
# add species name to fasta header
