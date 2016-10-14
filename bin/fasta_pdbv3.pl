#!/usr/bin/perl -w

use strict;
use FileHandle;
use Data::Dumper;

my $inFile = $ARGV[0];
my $outFile = $ARGV[1];
my $targetChain = $ARGV[2];
$targetChain = uc $targetChain;

my $hAAMap = {
"ALA" => "A",
"VAL" => "V",
"PHE" => "F",
"PRO" => "P",
"MET" => "M",
"ILE" => "I",
"LEU" => "L",
"ASP" => "D",
"GLU" => "E",
"LYS" => "K",
"ARG" => "R",
"SER" => "S",
"THR" => "T",
"TYR" => "Y",
"HIS" => "H",
"CYS" => "C",
"ASN" => "N",
"GLN" => "Q",
"TRP" => "W",
"GLY" => "G",
};

my $fhIn = new FileHandle($inFile,"r");
my $hRes = {};

#read in pdb data;
while(my $line = $fhIn->getline)
{
	if($line =~ /^ATOM/)
	{
	my $resNum = substr $line, 22,4;
	$resNum =~ s/\s+//;
	my $resName = substr $line, 17,3;
	#print $resNum." ".$resName."\n";
	my $chain = substr $line, 21, 1;
	$chain = uc $chain;
	
	if($chain =~ /$targetChain/)
	{
		$hRes->{$chain}{$resNum} = $resName;
	}
	}
}

#what's the numerical residue sequence we're going to print?
my $lowest = 0;
my $highest = 0;
foreach my $resNum (keys %{$hRes->{$targetChain}})
{
	if($lowest == 0)
	{
		$lowest = $resNum;
	}
	if($resNum < $lowest)
	{
		$lowest = $resNum;
	}
	if($resNum > $highest)
	{
		$highest = $resNum;
	}
}

#Now print the fasta file
my $fhOut = new FileHandle($outFile,"w");
print $fhOut ">".$inFile."| chain ".$targetChain."\n";

for(my $i=$lowest; $i <= $highest;$i++)
{
	if(exists $hRes->{$targetChain}{$i})
	{
		my $tlc = $hRes->{$targetChain}{$i};
		$tlc = uc $tlc;
		print $fhOut $hAAMap->{$tlc}
	}
	else
	{
		print $fhOut "-";
	}
}
