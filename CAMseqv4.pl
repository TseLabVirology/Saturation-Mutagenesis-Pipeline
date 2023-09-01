use strict; use warnings;

my $help = $ARGV[0];
if ($help =~ /help/){print "\n ***CAMseqv4 Help*** \n  
usage: sequencingfile.fastq 5'flankingsequence 3'flankingsequence [R]
Where flankingsequences surround your region of interest \n
Optional command-line input: R for results as reverse complement, default is no revcomp 
Flanking sequences must run in same direction as sequence as input fastq file \n
V4 forces a specific length between motifs. To allow any length between motifs, using V3 code found in Tse et al. PNAS 2017 \n

Further questions can be directed to tselabvirology\@gmail.com \n"; die;
}

die "usage: sequencingfile.fastq 5'flankingsequence 3'flankingsequence \n help for more info \n" unless 
@ARGV > 2 ;

#open .fastq sequencing file and push lines into array 
my $file = $ARGV[0];
my $infilename = "".$ARGV[0]."";
my @lines = ();
my $line = 0;
  open (IN, "< $ARGV[0]") or die "Can't open $file for read: $!";
while (<IN>) {
      push (@lines, $_);
      }

#count the total reads and the number of reads that include flanking sequences
my $totallines = scalar(@lines);
my $readcount = 0;
my $motif = $ARGV[1];
my $motif2= $ARGV[2];
my $linecount=0;
my @seqswithflanks= ();
#the nofastq file will allow output with just seq
foreach (@lines) {$line = $_; $linecount= ($linecount + 1); 
		if ( $line =~ /$motif.......................................$motif2/) {$readcount = ($readcount+1)};
		if ( $line =~ /$motif(.......................................)$motif2/) {push (@seqswithflanks, $1); 
		}
		}

my @finalDNAseqs= ();

if(@ARGV > 3){my $revcompdecision = $ARGV[3];
if ($revcompdecision =~ /R/){print "Reads have been revcomped \n"; 
	foreach my $revcompDNA (@seqswithflanks){my $revcomp = reverse($revcompDNA);
		$revcomp =~ tr/ACGTacgt/TGCAtgca/;
		push (@finalDNAseqs, $revcomp);
}}else{print "R option selects reverse complementing of sample, default is no revcomp \n";
	foreach(@seqswithflanks){my $seqtopush=$_; push (@finalDNAseqs, $seqtopush);
}}}

if(@ARGV < 4){foreach(@seqswithflanks){my $seqtopush=$_; push (@finalDNAseqs, $seqtopush);
}}

#Compute basic statistics of .fastq file, number of unique seqs, and reads per unique seq
my $seqswithflanks = scalar(@seqswithflanks);
my %readcount2;
map($readcount2{$_}++, @finalDNAseqs);
my $uniquecount = keys (%readcount2);


#Create output file with basic statistics of sequencing file
my $samplename = $infilename;
$samplename =~ s/\.fastq//;
open (OUT, ">" . $samplename . "statistics.txt") or die "error creating 
" . $samplename . "statistics.txt \n";
my $sequences = ($totallines/4);
print OUT "Total number of lines: " . $totallines . "\n";
print OUT "Total sequencing lines: " . $sequences . "\n";
print OUT "Reads with " . $motif . " start motif: " . $readcount . "\n";
print OUT "Reads with both " . $motif . " start and " . $motif2 . " end motif: " 
. $seqswithflanks . "\n";
print OUT "Unique reads with both start and end motifs: " . $uniquecount 
. "\n";
if(@ARGV > 3){my $revcompdecision = $ARGV[3];
if ($revcompdecision =~ /R/){print OUT "Samples output in Reverse Complement \n"}};
close OUT or die "Cannot close " . $samplename . "statistics.txt \n";


#Create output file with list of nucleotide sequences & number of occurrences(high to low)
open (OUT, ">" . $samplename . "nucleotidelist.txt") or die "error creating 
" . $samplename . "nucleotidelist.txt \n";
print OUT join("\n",map($_."\t".$readcount2{$_}, sort
{$readcount2{$b}<=>$readcount2{$a}}
keys (%readcount2))) . "\n";
close OUT or die "Cannot close " . $samplename . "nucleotidelist.txt \n";

#Create output file with list of translated amino acid seqs & number of occurrences(high to low)
open (OUT, ">" . $samplename . "aminoacidlist.txt") or die "error creating 
" . $samplename . "aminoacidlist.txt \n";
sub codon2aa{
my($codon)=@_;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F',
'TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C',
'TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
'CCA'=>'P','CCG'=>'P', 'CCC'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R',
'CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T',
'ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAG'=>'K', 'AAA'=>'K',
'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V',
'GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G',
'GGG'=>'G','GGT'=>'G');
if(exists $g{$codon})
{
return $g{$codon};
}
else
{
$codon = lc $codon; }
}

#Translate DNA sequences, maintaining any extra out of frame nucleotides as nucleotides
my $protein;
my $codon;
my @translated = ();
my @alltranslated = ();
foreach my $DNA (@finalDNAseqs){for(my $i=0; $i<(length($DNA)-2);
$i+=3){my $AAnumber = (length($DNA) % 3); 
	if($i<(length($DNA)-4-$AAnumber)){
        $codon = substr($DNA, $i, 3);
        $protein=&codon2aa($codon); push(@translated, $protein);}else{
	if ($AAnumber==0){
        $codon = substr($DNA, $i, 3); $protein=&codon2aa($codon);
        push(@translated, $protein);
        push(@alltranslated, join("", @translated)); 
        @translated = (); next;}else{$codon=substr($DNA, $i, 3);
	$protein=&codon2aa($codon); push(@translated, $protein);
	$codon=substr($DNA, $i + 3, $AAnumber);
	$protein=&codon2aa($codon); push(@translated, $protein);

	push(@alltranslated, join("", @translated)); @translated=(); next; }
        }
	}
	}
my %readcount3;
map($readcount3{$_}++, @alltranslated);
my $uniquecount2 = keys (%readcount3);
print OUT  "Unique Amino Acid sequences: " . $uniquecount2 . "\n";
print  OUT join("\n",map($_."\t".$readcount3{$_}, sort 
{$readcount3{$b}<=>$readcount3{$a}}
keys (%readcount3))) . "\n";
close OUT or die "Cannot close " . $samplename . "aminoacidlist.txt \n";