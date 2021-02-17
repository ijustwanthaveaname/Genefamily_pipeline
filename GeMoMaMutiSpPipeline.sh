#!/usr/bin/zsh
# -*- coding: utf-8 -*-

# Get options from command line
usage() {
	echo "Usage\n\tzsh /pathToScript/ -t <TargetGenome> -d <DirToSpecies> -o <OutMergedDir> -n <threads>"
}
while {getopts t:d:o:n:h arg} {
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
		(d)
		if [[ ! -d $OPTARG ]] {
			echo "The directory you specified does not exist!"
			exit
		} else {
			workdir=$OPTARG
		}
		;;
		(o)
		if [[ ! -d $OPTARG ]] {
			echo "The output directory does not exist and will be created."
			mkdir $OPTARG
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
for sp (`ls $workdir`) {
	if [[ -d $workdir/$sp ]] {
		for file (`ls $workdir/$sp`) {
			if [[ $file =~ ".*\.fna$" ]] {
				refgenome=$workdir/$sp/$file
			} elif [[ $file =~ ".*\.faa$" ]] {
				queryprotein=$workdir/$sp/$file
			} elif [[ $file =~ ".*\.gff$" ]] {
				refgff=$workdir/$sp/$file
			} else {
				echo "The directory that you specified is not contain correct directories!"
				exit
			}
		}
		echo "Proccessing $sp ..."
		python $GENEFPATH/main/parsegff.py -g $refgenome -h $queryprotein -a $refgff -b ${queryprotein%.faa}.gff -f ${queryprotein%.faa}GenomeRegion.fna && echo "Successfully get queryGFF and queryRegion."
		wait 
		echo "Start GeMoMa pipeline for $sp ..."
		zsh $GENEFPATH/GeMoMaOneSpPipeline.sh tblastn $target ${queryprotein%.faa}.gff ${queryprotein%.faa}GenomeRegion.fna $threads $workdir/$sp && echo "Successfully get predicted region and protein of $sp."
	}
	wait
}
wait

# Step.2 Merge gff and exclude overlap records in gff
echo "Start merging gff files..."
find $workdir -name final_annotation.gff | while read gff;do mv $gff $outdir/${${gff:h}:t}_${gff:t};done
zsh $GENEFPATH/main/mergegff.sh -d $outdir && echo "Successfully merge gff files!"

# Step.3 Extract cds from target genome
echo "Start extracting cds from target genome..."
gffread $outdir/filtermerged.gff -g $target -y $outdir/predicted_cds.faa && echo "Successfully extract cds!\nAll done!" 
