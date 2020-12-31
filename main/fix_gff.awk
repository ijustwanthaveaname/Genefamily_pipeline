#!/usr/bin/gawk -f 
BEGIN{FS="\t";OFS="\t"}$7=="+"{
	#if ($1~/\(\+\)/)
	#	$7="+"
	#else
	#	$7="-"
	if ($3=="gene"){
		scores = $6
		start=$4
		end=$5
		split($9, a, "[ ]")
		#if ($1~/\(\+\)/)
			n=++b[$1,a[5]]
		#else
		#	n=++c[$1,a[5]]
		$9=""
		id=$1":"a[5]":"n":"start":"end":"scores
		print $0,"gene_id "id}
	else if ($3=="cds"){
		$9=""
		print $0"transcript_id "id" ;gene_id "id" ;"}}
