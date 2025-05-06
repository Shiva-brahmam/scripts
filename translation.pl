#!/usr/bin/perl
use strict;
print "enter the dna sequence file name:\n";
chomp(my $file_name=<STDIN>);
my $dna=extraxt_dna_seq($file_name);
my $protein;
my $codon='';
for (my$i=0; $i<length($dna);$i+=3){
$codon = substr($dna,$i,3);
}
$protein =print_protein_seq($codon);
print "the protein sequence:\n\n$protein\n\n";

exit;
#first subroutine
sub extract_dna_Seq{
my($filename)=@_;
my $dnaseq='';
unless(open(DNA_SEQ,$filename)){
print "cannot open the file\n";
exit;
my @dna=<DNA_SEQ>;
foreach my$line (@dna){
if ($line =~ /^>/){
next;}
elsif ($line=~ /^\s*$/){
next;}
elsif ($line =~ /^\s*#/){
next;}
else {
$dnaseq.=$line;
}}
$dnaseq =~ s/\s//g ;
return $dnaseq;
}
#second subroutine
sub print_protein_seq{
my($codon)=@_;
my $protein='';
my(%aa)=(
  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
  TGT => "C", TGC => "C", TGA => "*", TGG => "W",
  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G",);
if (exists $aa{$codon}){
$protein.=$aa{$codon};
}
else{
print "bad codon!\n";
}
return $protein;
}
