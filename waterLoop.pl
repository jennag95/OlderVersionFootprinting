#!/usr/bin/perl

use warnings; use strict; use Getopt::Std;
use vars qw($opt_n $opt_t $opt_l $opt_r $opt_q $opt_s $opt_g $opt_c);
getopts("n:t:l:r:q:s:g:c");

my $usage = "usage: $0 [options] -n <gene name> -t <conversion threshold> -l <minimum R-loop length> -r <path to reads> -q <minimum mapping quality for each read to be considered> -s <path to original gene sequence> -g <path to reference genome for indexing> 

Options:
-c: consider Cs in CpG context

Extra information:
-n: gene name must be in all caps (e.g. CALM3)
-t: percentage (0.0-1.0) of Cs that must be converted to be considered an R-loop
-l: minimum length in base pairs to be considered an R-loop
-r: must be path to reads file itself named pacbio.fastq
-q: must be phred score (e.g. score of 40 = 99.99% base call accuracy)
-s: must be path to sequence file itself
-g: must be path to .fa file containing reference genome used for indexing

-Must have Bisulfite_Genome directory in same directory as this executable.
";

die $usage if not ($opt_n) or not ($opt_t) or not ($opt_l) or not ($opt_r) or not ($opt_q) or not ($opt_s) or not ($opt_g);

#runs bismark (output file will be pacbio.fastq_bismark_bt2.sam) only if it hasn't been ran previously
unless(-e "./pacbio.fastq_bismark_bt2.sam")
{
	#allows user to enter coordinates of PCR primers to create indexes for genes. Genes are indexed +/- 50bp on either side. 
	my $geneIndexes = "geneIndexes.bed";
	open(my $indexFile, ">", $geneIndexes) or die "Could not open $geneIndexes: $!";

	my $user = "y";
	print "Gene indexes:\n";
	while($user eq "y")
	{	
		print "Enter chromosome number(e.g. chr2) followed by starting coordinate, ending coordinate, and gene name(in all caps). Each field should be separated by tab.\n";
		my $index = <STDIN>;
		my @indexes = split("\t", $index);
		$indexes[1] -= 50;
		$indexes[2] += 50;
		my $ind = join("\t", @indexes);
		print $indexFile "$ind"; 
		print "Would you like to input another index? (y or n)\n";
		$user = <STDIN>;
		chomp($user);
	}

	system("bedtools getfasta -fi $opt_g -bed geneIndexes.bed -fo geneIndexes.fa -name") or die "Failed to run bedtools: $!\n";

	system("bismark_genome_preparation --bowtie2 ./") or die "Failed to run bismark genome preparation: $!\n";

   system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 ./ $opt_r") == 0 or die "Failed to run bismark: $!\n";
}

#takes sequence of gene and splits into array of individual bases
open(SEQ, $opt_s) or die "Could not open $opt_s: $!";
my @seq = split("", <SEQ>);
close SEQ;

my $positive = $opt_n . "Positive.txt";
my $negative = $opt_n . "Negative.txt";
open(my $positiveReads, ">", $positive) or die "Could not open $positive: $!";
open(my $negativeReads, ">", $negative) or die "Could not open $negative: $!";

my $samFile = "pacbio.fastq_bismark_bt2.sam";
open(my $sam, $samFile) or die "Could not open $samFile: $!";

#loops through each read and writes high quality reads into two separate files <gene>Positive.txt and <gene>Negative.txt
while(my $line = <$sam>)
{
	chomp($line);
   my @fields = split("\t", $line);
	#discounts the first 4 lines of information at the top of the .sam file
   if(@fields < 6)
   {
      next
   }
	#discounts if mapping quality is not at least phred score specified
   if($fields[4] < $opt_q)
   {
      next
   }
	#discounts any read that is not the gene of interest, the proper length, or at the proper location in genome(accounting for indexing)
	if($fields[2] ne "$opt_n" || length($fields[9]) < (@seq) || $fields[3] < 45 || $fields[3] > (@seq+55))
   {
      next
   }
	#counts number of CT conversions (takes into account whether or not the user wants to include Cs in CpG context)
   my $CT = ($fields[13] =~ tr/xhu/xhu/);;
	if($opt_c)
	{
		my $CT = ($fields[13] =~ tr/zxhu/zxhu/);
	}
	#writes positive reads into <gene>Positive.txt and negative reads into <gene>Negative.txt
   if($fields[1] == 0 && $fields[2] eq "$opt_n")
   {
      print $positiveReads "$CT\t" . $line . "\n";
   }
   if($fields[1] == 16 && $fields[2] eq "$opt_n")
   {
      print $negativeReads "$CT\t" . $line . "\n";
   }
}

#sorts by number of CT conversions, takes top 1000, and removes the number of CT conversions in prepartion for methylation extractor
my $finalPositive = $opt_n . "PositiveFinal.txt";
system("sort -k1,1rn $positive | head -n 1000 | cut -f 2- > $finalPositive");

my $finalNegative = $opt_n . "NegativeFinal.txt";
system("sort -k1,1rn $negative | head -n 1000 | cut -f 2- > $finalNegative");

my $CPGpos = "CpG_context_" . $finalPositive;
my $CPGneg = "CpG_context_" . $finalNegative;
my $CHGpos = "CHG_context_" . $finalPositive;
my $CHGneg = "CHG_context_" . $finalNegative;
my $CHHpos = "CHH_context_" . $finalPositive;
my $CHHneg = "CHH_context_" . $finalNegative;

#runs bismark methylation extractor on top 1000 reads of each strand
my $bismarkOutput = "./" . $CPGpos;
unless(-e $bismarkOutput)
{
	system("bismark_methylation_extractor -s --comprehensive $finalPositive");
	system("bismark_methylation_extractor -s --comprehensive $finalNegative");
}

#pulls the CHH, CHG, and CpG sites together into methylationPos<gene>.txt and methylationNeg<gene>.txt (takes into account -c option)
my $methylationPos = "methylationPos" . $opt_n . ".txt";
my $methylationNeg = "methylationNeg" . $opt_n . ".txt";
$methylationPos = "methylationPos" . $opt_n . "CG.txt" if($opt_c);
$methylationNeg = "methylationNeg" . $opt_n . "CG.txt" if($opt_c);

if($opt_c)
{
	system("cat $CPGpos $CHGpos $CHHpos | sort -n > $methylationPos");
	system("cat $CPGneg $CHGneg $CHHneg | sort -n > $methylationNeg");
}
else
{
	system("cat $CHGpos $CHHpos | sort -n > $methylationPos");
	system("cat $CHGneg $CHHneg | sort -n > $methylationNeg");
}

my $filePos = $opt_n . "PosPositions.txt";
my $fileNeg = $opt_n . "NegPositions.txt";

open(FILEPOS, ">", $filePos) or die "Could not open $filePos: $!";
open(FILENEG, ">", $fileNeg) or die "Could not open $fileNeg: $!";

#gets the position of each conversion (conversions=1; not converted=0)
my $j = 50;
for(my $i=0; $i<2; $i++)
{
	my $analyzeFile = $methylationPos if($i==0);
	$analyzeFile = $methylationNeg if($i==1);
	open(my $analyze, $analyzeFile) or die "Could not open $analyzeFile: $!";
	my (@read) = ("0\t")x(@seq);
	while(my $line = <$analyze>)
	{
		chomp($line);
		next if $line  =~ /Bismark/;
		my @fields = split("\t", $line);
		if($read[0] ne "$fields[0]\t")
   	{
      	if(@read == (@seq) && $read[0] ne "0\t" && $read[0] ne "Bismark methylation extractor version v0.13.1\t")
      	{
				if($i==0)
				{
         		print FILEPOS "@read\n";
				}
				if($i==1)
				{
					print FILENEG "@read\n";
      		}
			}
     		@read = ();
      	@read = ("0\t")x(@seq);
      	$read[0] = "$fields[0]\t";
   	}
		if($read[0] eq "$fields[0]\t")
   	{
      	if($fields[4] eq "x" || $fields[4] eq "h" || $fields[4] eq "z" || $fields[4] eq "u" && $fields[3]>$j)
      	{
         	my $index = ($fields[3] - $j + 1);
         	$read[$index] = "1\t";
      	}
   	}
	}
}

#determines which regions of conversion are R-loops based on conversion threshold
#0=not converted (grey)
#1=converted (green)
#2=non-C (white)
#3=converted CpG (blue)
#4=non-converted CpG  (black)
#5=converted CpG in R-loop (purple)
#9=converted C in R-loop (red)

my $start = 1;
my $end = $opt_l;
my $conC = 0;
my $nonConC = 0;
my $conPer = 0;

my $finalPos = $opt_n . "Pos" . ($opt_t*100) . ".txt";
my $finalNeg = $opt_n . "Neg" . ($opt_t*100) . ".txt";
$finalPos = $opt_n . "Pos" . ($opt_t*100) . "CG.txt" if($opt_c);
$finalNeg = $opt_n . "Neg" . ($opt_t*100) . "CG.txt" if($opt_c);
open(FINALPOS, ">", $finalPos) or die "Could not open $finalPos: $!";
open(FINALNEG, ">", $finalNeg) or die "Could not open $finalNeg: $!";
for(my $i=0; $i<2; $i++)
{
	my $fileLast = $opt_n . "PosPositions.txt" if($i==0);
	$fileLast = $opt_n . "NegPositions.txt" if($i==1);
	my @lineFinal = `cat $fileLast`;
  	for (my $p =0 ; $p < @lineFinal; $p++) 
	{
		my $lineFinal = $lineFinal[$p];
		chomp($lineFinal);
		my @fields = split("\t", $lineFinal);
		for(my $z=0; $z<@fields; $z++)
		{
			$fields[$z] =~ s/ //g;
		}
		for(my $z=0; $z<@seq; $z++)
		{
			#positive reads
			if($i==0)
			{
				if($seq[$z] ne "C")
				{
					$fields[$z+1] = 2;
				}
				if($opt_c)
				{
					if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G" && $fields[$z+1] == 1)
      			{
         			$fields[$z+1] = 3;
      			}
      			if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G" && $fields[$z+1] == 0)
      			{
         			$fields[$z+1] = 4;
      			}	
				}
				else
				{
					if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G")
      			{
         			$fields[$z+1] = 2;
      			}
				}	
			}
			#negative reads
			if($i==1)
			{
				if($seq[$z] ne "G")
				{
					$fields[$z+1] = 2;
				}
				if($opt_c)
				{
					if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C" && $fields[$z+1] == 1)
               {
                  $fields[$z+1] = 3;
               }
               if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C" && $fields[$z+1] == 0)
               {
                  $fields[$z+1] = 4;
               }
				}
				else
				{
					if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C")
      			{
         			$fields[$z+1] = 2;
      			}
				}
			}
		}
		while($end<@fields)
		{
			for(my $temp = $start; $temp<$end; $temp++)
			{
				if($fields[$temp] == 1 || $fields[$temp] == 9 || $fields[$temp] == 3 || $fields[$temp] == 5)
				{
            	$conC++;
         	}
         	if($fields[$temp] == 0 || $fields[$temp] == 4)
         	{
            	$nonConC++;
         	}
			}
			if($conC != 0 || $nonConC != 0)
			{
				$conPer = (($conC)/($conC+$nonConC));
			}
			else
			{
				$conPer = 0;
			}
			if($conPer >= $opt_t)
			{
				for(my $k=$start; $k<=$end; $k++)
				{
					if($fields[$k] == 1)
					{
						$fields[$k] = 9;
					}
					if($fields[$k] == 3)
					{
						$fields[$k] = 5;
					}
				}
			}
			$start++;
			$end++;
			$conC = 0;
			$nonConC = 0;
		}
		my $newfield = join("\t", @fields);
		if($i==0)
		{
			print FINALPOS "$newfield\n";
		}
		if($i==1)
		{
			print FINALNEG "$newfield\n";
		}
		$start = 1;
   	$end = 100;
   	$conC = 0;
   	$nonConC = 0;
   	$conPer = 0;
	} 
}

#makes heatmaps for positive and negative strand
my $finalPosPDF = $opt_n . "Pos" . ($opt_t*100) . ".pdf";
my $finalNegPDF = $opt_n . "Neg" . ($opt_t*100) . ".pdf";
$finalPosPDF = $opt_n . "Pos" . ($opt_t*100) . "CG.pdf" if($opt_c);
$finalNegPDF = $opt_n . "Neg" . ($opt_t*100) . "CG.pdf" if($opt_c);

my $Rscript = "MakeHeatmap.R";
open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";
if($opt_c)
{
	print $out "
		library(\"GMD\")
		df = read.table(\"./$finalPos\", sep=\"\t\", row.names=1)
		pdf(\"$finalPosPDF\")
		heatmap.3(
			x=df,
			dendrogram=\"row\",
			Rowv=TRUE, Colv=FALSE,
			labRow=FALSE,labCol=FALSE,
			breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,9.5),
			color.FUN=function(x) c(\"grey\",\"green\",\"white\",\"blue\",\"black\",\"purple\",\"red\")
		)
		df2 = read.table(\"./$finalNeg\", sep=\"\t\", row.names=1)
		pdf(\"$finalNegPDF\")
      heatmap.3(
         x=df2,
         dendrogram=\"row\",
         Rowv=TRUE, Colv=FALSE,
         labRow=FALSE,labCol=FALSE,
         breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,9.5),
         color.FUN=function(x) c(\"grey\",\"green\",\"white\",\"blue\",\"black\",\"purple\",\"red\")
      )
		dev.off()
	";
}
else
{
   print $out "
      library(\"GMD\")
      df = read.table(\"./$finalPos\", sep=\"\t\", row.names=1)
      pdf(\"$finalPosPDF\")
      heatmap.3(
         x=df,
         dendrogram=\"row\",
         Rowv=TRUE, Colv=FALSE,
         labRow=FALSE,labCol=FALSE,
         breaks=c(-0.5,0.5,1.5,2.5,9.5),
         color.FUN=function(x) c(\"grey\",\"green\",\"white\",\"red\")
      )
		df2 = read.table(\"./$finalNeg\", sep=\"\t\", row.names=1)
      pdf(\"$finalNegPDF\")
      heatmap.3(
         x=df2,
         dendrogram=\"row\",
         Rowv=TRUE, Colv=FALSE,
         labRow=FALSE,labCol=FALSE,
         breaks=c(-0.5,0.5,1.5,2.5,9.5),
         color.FUN=function(x) c(\"grey\",\"green\",\"white\",\"red\")
      )
      dev.off()
   ";
}
close $out;
system("R --vanilla --no-save < $Rscript"); 
