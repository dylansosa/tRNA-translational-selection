for b in *.highConfidence.bed ; do
name=$(ls $b | rev | cut -f 4- -d '.' | rev)
assembly=$(ls $b | rev | cut -f 5- -d '.' | rev)
	bedtools getfasta -fi ../"$assembly"_genomic.fna -bed $b -fo $name.highConfidence.fa -name ; done
