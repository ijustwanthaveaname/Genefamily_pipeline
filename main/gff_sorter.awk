#!/usr/bin/gawk -f 
BEGIN{
	FS="gff_"
	OFS="gff_"
	PROCINFO["sorted_in"]="@val_num_asc"
}
{
	ARRAY[$0]=$2
}
END{
	for (i in ARRAY){
		print(i)
	}
}
