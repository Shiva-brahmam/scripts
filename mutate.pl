#!/usr/bin/perl
use strict;
print "enter a dna seq:\n";
chomp (my $dna=<STDIN>);
print "number of mutations to be done:";
chomp (my $n=<STDIN>);
my $i;
my $mutant;
srand (time |$$);
$mutant =mutate($dna);
print "the original dna sequence:\n$dna\n";
print "the mutated dna sequence :\n$mutant\n";
for ($i=0; $i<$n; $i++){
$mutant=mutate($mutant);
print "$mutant \n";
}
exit;
sub mutate{
my($dna)=@_;
my(@nucleotide)=('A','T','G','C');
my ($position)=randomposition($dna);
my($newbase)=randomnucleotide(@nucleotide);
substr($dna,$position,1,$newbase);
return $dna;
}
sub randomnucleotide{
my(@bases)=@_;
return $bases[rand@bases];
}
sub randomposition{
my($string)=@_;
return int rand length($string);
}


