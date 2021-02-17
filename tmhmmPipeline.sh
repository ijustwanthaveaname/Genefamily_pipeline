#!/usr/bin/zsh
### --------------------------------------
# Run tmhmm with multiple processes and extract the results for transmembrane segments larger than a certain number.
### --------------------------------------

# Get options from command line
usage() {
	echo "Usage:\n\tzsh /pathToScript/tmhmmPipeline.sh -p <protein> -n <threshold> -t <threads> "
}

while {getopts hp:o:t:n: arg} {
	case $arg {
		(h)
		usage
		exit
		;;
		(p)
		protein=$OPTARG
		if [[ ! -e $protein ]] {
			echo "Please input correct path to protein!"
			exit
		}
		;;
		(t)
		if [[ -n "`echo $OPTARG | sed 's/[0-9]//g'`" ]] {
			"You must use int to specify the threads!"
			exit
		} else {
			threads=$OPTARG	
		}
		;;
		(n)
		if [[ -n "`echo $OPTARG | sed 's/[0-9]//g'`" ]] {
			"You must use int to specify the counts of transmembrane!"
			exit
		} else {
			trans=$OPTARG	
		}
		;;
		(?)
		echo "You must use -p to specify protein file and -o to specify output file."
		usage
		exit
		;;
	}
}

# Start tmhmm pipeline
# Step1 Split cds
python $GENEFPATH/main/split_genome.py -g $protein -n $threads

# Step2 Running tmhmm
ls ${protein:h} |
	grep -E "${protein}_[0-9]+$" | 
	while read id
	do
		nohup tmhmm $id -short >tmhmm_result.txt_${id##*_} &
	done
wait
cat ${protein:h}/tmhmm_result.txt_* > ${protein:h}/tmhmm_result.txt && rm -rf ${protein:h}/tmhmm_result.txt_*  ${protein}_<->

# Step3 Get PredHel>=trans id
gawk -F"[\t=]" -v n=$trans '$9>=n{print $1}' ${protein:h}/tmhmm_result.txt > ${protein:h}/tmhmmGe${trans}_id.txt

# Step4 Get CDS by id
python $GENEFPATH/main/extract_by_id.py -f $protein -i ${protein:h}/tmhmmGe${trans}_id.txt > ${protein%.*}_tmhmm_Ge${trans}.cds
