#!/usr/bin/perl -w
use strict;
my $sum_gc = 0;
my $sum_N = 0;
my $sum_all = 0;
while (<>){
	chomp;
	if ( ! /^>/){
		$sum_gc += tr/Gg/Gg/;
		$sum_gc += tr/Cc/Cc/;
		$sum_N += tr/Nn/Nn/;
		$sum_all += length;
	}
}
my $gc_content = $sum_gc / ($sum_all-$sum_N);
printf "%.4f\n", $gc_content;
