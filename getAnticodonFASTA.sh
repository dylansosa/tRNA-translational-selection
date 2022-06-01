for i in *.highConfidence.fa  ; do

name=$(ls $i | rev | cut -f 4-  -d '.' | rev)


cat anticodonsPerIsotype.txt | while read -r aa ; do

cat $i | egrep $aa -A1 --no-group-separator > $name.$aa.anticodonPerIsotype.fa ; done ; done

# get for all isotypes in a loop
# make files for each isotype anticodon 
