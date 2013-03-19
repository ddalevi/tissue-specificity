#!/usr/bin/perl -w
# This script takes n individual analyses (n files with inner scores)
# and combine them to calculate the total score.
# Syntax: ./combine-data.pl file1.txt file2.txt ... fileN.txt
  
use strict;

my %symbol = ();
my $outfile = "results.txt";

foreach my $e (@ARGV) 
{
    open( FILE, $e ) or die "Cannot open file: $e\n";
    while( <FILE> )
    { 
	my @p = split "\t\t\t", $_;
	$symbol{ $p[0] } .= $p[1];
    }
    close( FILE );
} 

open( OUT, ">$outfile" ) or die "Cannot open file: $outfile\n";
   
foreach my $gene (keys %symbol)
{
    my @set = split (/END/,$symbol{$gene});
    my %tissue=();
    foreach my $n( @set )
    {
	if ($n =~ /\t/)
	{ 
	    $n =~ s/\n//g;
	    $n =~ s/\r//g;
	    my @t=split (/\t/,$n);  #tt
	    for( my $i=0; $i<$#t; $i=$i+2 )
	    {
		$tissue{$t[$i]} += $t[($i+1)];
	    }
	}
    }
    
    foreach my $t ( keys %tissue )
    {	
	$tissue{ $t }= $tissue{ $t }/( $#set );
    }
    print OUT "$gene\t$#set\t";
    foreach my $ti (sort { $tissue {$b} <=> $tissue {$a}} (keys %tissue))
    {	
        print OUT "$ti\t$tissue{$ti}\t";
    }
    print OUT "\n";
}
close(OUT);

print "Results written to file: $outfile\n";
