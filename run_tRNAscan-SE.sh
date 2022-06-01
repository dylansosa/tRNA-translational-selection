for f in ../*_genomic.fna ; do

name=$(ls $f | cut -f 1-4 -d '_' | cut -f 2 -d '/')

#tRNAscan-SE -E --max -f -b -j -a -d -o $name.trnascan $f ; done
tRNAscan-SE -E -d -o ./$name.trnascan -a ./$name.trnascan.fa -b ./$name.trnascan.bed -f ./$name.trnascan.struct -j ./$name.trnascan.gff $f ; done
