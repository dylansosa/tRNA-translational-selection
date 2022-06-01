for i in *.trnascan.gff  ; do
cat $i | egrep 'pseudo|Undet' | cut -f 5 > $i.pseudogeneCoordinates

name=$(ls $i | rev | cut -f 2- -d '.' | rev)
egrep -wvf $i.pseudogeneCoordinates $name.bed > $i.highConfidence.bed ; done

