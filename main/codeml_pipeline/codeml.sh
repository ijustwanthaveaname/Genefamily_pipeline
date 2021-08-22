#!/usr/bin/zsh
usage() {
	echo "Usage:"
	echo "\tcodeml.sh -a <alignment> -t <tree> -o <outdir>\n"
}
# Get options from command line
set -- $(getopt -q a:t:o:h "$@")
while [ -n "$1" ]
do
	case "$1" in
		-a) align=`eval echo $2`
			shift;;
		-t) tree=`eval echo $2`
			shift;;
		-o) outdir=`eval echo $2`
			shift;;
		-h) usage
			exit;;
		--)shift
			break;;
		*) echo "$1 is not an option"
			usage
			exit;;
	esac
	shift
done
# Change parameters in ctl
pd=$PWD;ls $GENEFPATH/main/codeml_pipeline/*.ctl | while read id; do if [[ ! -d $outdir/${${id%.ctl}:t} ]] { mkdir $outdir/${${id%.ctl}:t} };sed "s:<seq>:${align:a}:;s:<tree>:${tree:a}:;s:<out>:${outdir:a}/${${id%.ctl}:t}/${${id%.ctl}:t}:" $id > $outdir/${${id%.ctl}:t}/${${id}:t};cd $outdir/${${id%.ctl}:t};codeml ${id:t} > log.txt 2>&1 &; cd $pd; done && echo "Running codeml(s) successfully!"
