#!/usr/bin/zsh
usage() {
	echo "Usage:"
	echo "\tgenef_pipeline.sh -g <genome> -q <query> -u <up> -d <down> -t <threads> -e <evalue>.\n"
	echo "The gffread, bedtools, blast and exonerate must be in ENV_PATH."
}
# Step1 read_arguments
threads=8
set -- $(getopt -q g:q:e:t:u:d:h "$@")
while [ -n "$1" ]
do
	case "$1" in
		-g) genome=`eval echo $2`
			shift;;
		-q) homology=`eval echo $2`
			shift;;
		-e) evalue=`eval echo $2`
			shift;;
		-t) threads=`eval echo $2`
			shift;;
		-u) up=`eval echo $2`
			shift;;
		-d) down=`eval echo $2`
			shift;;
		-h) usage
			exit;;
		--)shift
			break;;
		*) echo "$1 is not an option"
			exit;;
	esac
	shift
done
outname=`basename $genome`
outdir=`dirname $genome`
outprefix="`dirname $genome`/${outname%.fna}"

# Step2 blast
makeblastdb -dbtype nucl -in $genome -out $outprefix -parse_seqids
tblastn -db $outprefix -num_threads $threads -out ${outprefix}".blast" -outfmt 6 -evalue $evalue -query $homology

# Step3 filter overlapping region
python $GENEFPATH/main/region_tools.py -u $up -d $down -f -i $outprefix".blast" -g $genome -o $outprefix"_filter.blast"

# Step4 transform tabular into bed format
gawk 'BEGIN{FS="\t";OFS="\t"}{if ($9<$10){print $2,$9-1,$10,".",".","+"}else{print $2,$10-1,$9,".",".","-"}}'  $outprefix"_filter.blast" > $outprefix"_filter.bed"

# Step5 Use bedtools get region
bedtools getfasta -s -fi $genome -bed $outprefix"_filter.bed" -fo $outprefix"_candidate.fna"

# Step6 Split region and exonerate
python $GENEFPATH/main/split_genome.py -g $outprefix"_candidate.fna" -n $threads
ls $outdir | 
	grep -E ".*fna_[0-9]+$" |
	while read id
	do
		nohup exonerate -E -q $homology -t $id --model protein2genome:bestfit --querytype protein --targettype dna --showvulgar no --softmaskquery no --softmasktarget yes --showalignment no --showtargetgff yes --showcigar no --score 300 --bestn 0 --verbose 0 1> ${id%.fna*}.gff_${id##*_} 2> ${id%.fna*}.log_${id##*_} &
	done
wait
ls ${outdir}/*.gff_* | gawk -f $GENEFPATH/main/gff_sorter.awk | while read id 
do
	 cat $outdir/$id >> $outprefix".gff"
done
gawk -f $GENEFPATH/main/fix_gff.awk $outprefix".gff" | tr -s "     " > $outprefix"_fix.gff"  
gffread $outprefix"_fix.gff" -g $outprefix"_candidate.fna" -y $outprefix"_candidate.cds"
# Step6.1 Extract id
grep ">" $outprefix"_candidate.cds" | awk -F">" '{print $2}' > $outprefix"_candidate_id.txt"
# Step6.2 Transform into tabular
gawk -F":" '{print $3"\t"$1":"$2"\t"$4"\t.\t.\t.\t.\t.\t"$5"\t"$6"\t.\t"$7}' $outprefix"_candidate_id.txt" > $outprefix"_candidate_id.tabular"
# Step6.3 Sort and filter by scores
python $GENEFPATH/main/region_tools.py -s -f -i $outprefix"_candidate_id.tabular" -o $outprefix"_maxscore_id.tabular"
# Step6.3 Transform into idlist
gawk -F"\t" 'BEGIN{OFS=":"}{print $2,$1,$3,$9,$10,$12}' $outprefix"_maxscore_id.tabular" > $outprefix"_maxscore_id.txt"
# Step6.4 Get final cds
python $GENEFPATH/main/extract_by_id.py -f $outprefix"_candidate.cds" -i $outprefix"_maxscore_id.txt" > $outprefix"_maxscore.cds"
sed '/^>/!s/\./X/g' $outprefix"_maxscore.cds" > $outprefix"_maxscore_replaceX.cds"
