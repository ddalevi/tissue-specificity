#!/bin/env ke_perl
#use List::Util qw(max);
#inside one dataset, combine probesets
my $file1 = shift;
my $file2 = shift;
my $outfile = shift;

open( FILE, $file1 ) or die "Cannot open file $file1\n"; 
open( UNIQ, $file2 ) or die "Cannot open file $file2\n";
open( OUT, ">$outfile") or die "Cannot open file output.txt\n";

my @probeset = ();
my @symbol = ();
my @nr=();
my @t1=(),
my @t2=();

my @uniq;
while(<UNIQ>){
    my @p = split(" ", $_);
    push (@uniq,$p[0]);
  }
while( <FILE> )
  {
    chomp;
	$_ =~ s/\r//g;
    my @p = split " ", $_;
	#unshift (@p,"nr");	#only for gse
	push( @probeset, $p[1]);
	push( @symbol, $p[2]);
	push( @nr, $p[3]);
	push( @t1, $p[4]);
	push( @t2, $p[5]);
  }
@map = @uniq;

for($i=0; $i<=$#uniq; $i++){
  for($j=0; $j<=$#symbol;$j++){

      if($uniq[$i] eq $symbol[$j]){
	      $map[$i] = $map[$i]."###".$j;   #ss
	  }     

 }
}
#print "@map\n";

foreach my $gene (@map){
  my %tissue=();
  my %score=();
  my @p = split (/###/,$gene);
  my @sorted_keys = ();
  my $hk = 0; #number of HK probesets per gene
  my $sp = 0; #number of probesets after deletion
  my $largest=0;
  my $largest2=0;
  for($i=1;$i<=$#p;$i++){
	  $tissue{$t1[$p[$i]]} += 1;
	  $tissue{$t2[$p[$i]]} += 1;
	  if($nr[$p[$i]]==0) {$hk=$hk+1;}
	  else{
	  $score{$t1[$p[$i]]} += 1/$nr[$p[$i]];
      $score{$t2[$p[$i]]} += 1/$nr[$p[$i]];
      }
  }
  if(exists $tissue{"NA"}){
   delete $tissue{"NA"};
  }
   if(exists $score{"NA"}){
   delete $score{"NA"};
  }

  if (keys %tissue){
    @sorted_keys = sort {$tissue{$b} <=> $tissue{$a}} (keys %tissue);
	
	if ($#sorted_keys>=0){$largest = $tissue{$sorted_keys[0]};}

    if ($#sorted_keys>=1){$largest2 = $tissue{$sorted_keys[1]};}
	}
	
    if ($largest>=$#p/2){$sp = $#p-$hk;} 
    else {$sp = $#p;} 
   

    if($largest2>$sp/2 ) {
      print OUT "$p[0]\t\t\t$sorted_keys[0]\t0.5\t$sorted_keys[1]\t0.5\tEND\n";
    }
	
    elsif($largest >= $sp/2 && $largest2 != $sp/2) {
      print OUT "$p[0]\t\t\t$sorted_keys[0]\t1\tEND\n";
      }
    else{   
     foreach $n (keys %score){
     $score{$n}="\t".$score{$n}/($sp)."\t";#tt
      }
     print OUT "$p[0]\t\t\t";#aa
     print OUT %score;
     print OUT "END\n";
    }

   
}
close(FILE);
close(UNIQ);
close(OUT);

