#!/usr/bin/env perl

use strict;
use warnings;

sub salva{
	# $lista[0] - pmid
	# $lista[1] - abstract
	my @lista = @_;
	open(my $OUTFH, ">", "abstracts/".$lista[0].".txt");
	$lista[1] =~ s/\t+/ /g;
	print $OUTFH $lista[1];
	close($OUTFH);
}

my $FH; 
open($FH, "<", "abstracts.txt");

my $pmid;
my $abstract = "";
my $abflag = 0;

while(my $line = <$FH>){
	chomp($line);
	if($line =~ /^PMID-\s+(\d+)/){
		if($abstract){
			&salva($pmid, $abstract);
		}
		$pmid = $1;
		$abstract = "";
	}
	if($line =~ /^AB\s*-\s*(.*)/){
		$abflag = 1;
		$abstract = $1." ";	
	} else {
		if($abflag){
			if($line =~ /^\s+(.*)/){
				$abstract .= $1." ";
			} else {
				$abflag = 0;
			}
		}
	}
}
&salva($pmid, $abstract);
