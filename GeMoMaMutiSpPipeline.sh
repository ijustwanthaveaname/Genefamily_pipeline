#!/usr/bin/zsh
# -*- coding: utf-8 -*-

# Get options from command line
usage() {
	echo "Usage\n\tzsh /pathToScript/ -t <TargetGenome> -q <DirToQuerySpecies> -o <OutMergedDir> -n <threads> -e <evalue>"
}
while {getopts t:q:o:n:e:h arg} {
	case $arg {
		(h)
		usage
		exit
		;;
		(t)
		if [[ ! -f $OPTARG ]] {
			echo "The target genome file does not exist"
			exit
		} else {
			target=$OPTARG
		}
		;;
		(q)
		if [[ ! -d $OPTARG ]] {
			echo "The directory you specified does not exist!"
			exit
		} else {
			query=$OPTARG
		}
		;;
		(e)
		if [[ ! -n $OPTARG ]] {
			echo "Evalue not be set, default value is 100"
			evalue=100
		} else {
			evalue=$OPTARG
		}
		;;
		(o)
		if [[ ! -d $OPTARG ]] {
			echo "The output directory does not exist and will be created."
			mkdir -p $OPTARG
		}
		outdir=$OPTARG
		;;
		(n)
		threads=$OPTARG
		;;
		(?)
		echo "Please Check Your options!"
		usage
		exit
		;;
	}
}

# Start main pipeline 
# Step.1 Get query regions and gff from reference genome and gff
ls $query | while read sp
do
	echo "Proccessing $sp ..."
	if [[ -d $query/$sp ]] {
		for file (`ls $query/$sp`) {
			if [[ $file =~ ".*genomic\.fna$" ]] {
				refgenome=$query/$sp/$file
			} elif [[ $file =~ ".*query\.faa$" ]] {
				queryprotein=$query/$sp/$file
			} elif [[ $file =~ ".*genomic\.gff$" ]] {
				refgff=$query/$sp/$file
			} 
		}
		if [[ -f $query/$sp/query.gff && -f $query/$sp/queryregion.fna ]] {
			echo "The query.gff and queryregion.fna are existing."
		} else {
			python $GENEFPATH/main/parsegff.py -g $refgenome -h $queryprotein -a $refgff -b $query/$sp/query.gff -f $query/$sp/queryregion.fna && echo "Successfully get queryGFF and queryRegion."
		}
		wait 
		echo "Start GeMoMa pipeline for $sp ..."
		zsh $GENEFPATH/GeMoMaOneSpPipeline.sh tblastn $target $query/$sp/query.gff $query/$sp/queryregion.fna $threads $query/$sp $evalue && echo "Successfully get predicted region and protein of $sp."
	}
	wait
done
wait

# Step.2 Merge gff and exclude overlap records in gff
echo "Start merging gff files..."
for rec (final_annotation.gff predicted_proteins.fasta protocol_GeMoMaPipeline.txt) {
	find $query -name $rec | while read fp;do mv $fp $outdir/${${fp:h}:t}_${fp:t};done
}
zsh $GENEFPATH/main/mergegff.sh -d $outdir && echo "Successfully merge gff files!"
python $GENEFPATH/main/excludeOverlapInSortedGff.py -g $outdir/merged.gff -o $outdir/filtermerged.gff

# Step.3 Extract cds from target genome
echo "Start extracting cds from target genome..."
gffread $outdir/filtermerged.gff -g $target -y $outdir/predicted_cds.faa && echo "Successfully extract cds(AA)!"
gffread $outdir/filtermerged.gff -g $target -x $outdir/predicted_cds.fna && echo "Successfully extract cds(DNA)!"
echo "All done!"
