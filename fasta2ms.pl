#!/usr/bin/perl -w
use strict;

#######################################################################################################################
##                                                                                                                   ##
##  USAGE: perl fasta2ms.pl <infile> <m|p> [<out_prefix>] [<0|1>]                                                    ##
##                                                                                                                   ##
##  infile = file in fasta format with aligned sequences.  Sequences must be on a single lines                       ##
##  m = output all sites with < 3 alleles, (m)onomorphic; p = only output sites with 2 alleles, (p)olymorphic)       ##
##  out_prefix = optional, prefix to creat uniquely named outfiles, default = out                                    ##
##  0 = limited output to screen (default), 1 = verbose                                                              ##
##                                                                                                                   ##
##  fasta2ms.pl is designed to take as input fasta formated sequence file and tranform into ms formated files with   ##
##  or without monomorphic sites.  Ideal for processing long sequences (full chromosomes) from a limited number of   ##
##  individuals.  Note that fasta format must have sequence on a single line.                                        ##
##                                                                                                                   ##
##  Code written by G. John Lazur & Leah K. McHale (2012)                                                            ##
##  Last updated: 06_03_2014                                                                                         ##
##  The Ohio State University, Dept of Horticulture & Crop Science                                                   ##
##                                                                                                                   ##
#######################################################################################################################


my @nucleotide = ();            # nucleotides present in each sequence at given $pos
my $pos = 0;                    # nucleotide position on the sequence
my $seqnum = 0;                 # count of sequences
my $totalseq = 0;               # total # of sequences
my %position = ();              # key = position on the sequence, value = @nucleotide at that position
my $position = ();
my $seq;                        # sequence line from fasta file
my %count = ();                 # key = nucleotide, used to determine $countuniq
my $countuniq;                  # count of unique alleles (nucleotide) at each position
my $seqlength;                  # total length of input sequence, all sequences should have same total length
my %dump = ();                  # key = nucleotide position, value = 0 or 1 (1= exclude position from outfile)
my $subseq = "";                # string of length $sub derived from chromsome sequence
my $temppos = 1;                # count of temp files being generated
my @tempfiles = ();             # array of temp file names
my $subseqlength = 0;           # length of subsequence as string is processed
my $sub = 250000;               # length of subsequences which full sequences are broken into; can be optimized
my $out = "out";                # prefix of outfile name; uniquely ID outfile & tempfiles to allow same folder runs
my $segsites = 0;               # count of sites with 2 alleles
my $verbose = 0;                # provide progress update as subsequences are processed (1)
my $monomorph = 1;              # include monomorphic sites (1), or only include polymorphic sites (0)

sub help {
    print "USAGE: 
	
    perl fasta2ms.pl <infile> <m|p> [<out_prefix>] [<0|1>]
	
    for help file:
    perl fasta2ms.pl h
	
    infile = file in fasta format with aligned sequences  
        Note: sequences must be on a single line
	
    m = output all sites with < 3 alleles, (m)onomorphic 
    p = only output sites with 2 alleles, (p)olymorphic)
	
    out_prefix = optional, prefix to creat uniquely named outfiles, default = out
	
    0 = output to screen supressed (default) 
    1 = verbose output to screen
	
    fasta2ms.pl is designed to take as input fasta formated sequence file and tranform into ms formated files with
    or without monomorphic sites.  Ideal for processing long sequences (full chromosomes) from a limited number of
    individuals.  Note that fasta format must have sequence on a single line.
	
    Code written by G. John Lazur & Leah K. McHale (2012)
    Last updated: 06_03_2014
    The Ohio State University, Dept of Horticulture & Crop Science\n";
	
	die "\n";
}

if ($ARGV[0] eq "h") {
    help();
} elsif (@ARGV < 2) {
    die "USAGE:\nChrmFasta2ms.pl USAGE: perl fasta2ms.pl <infile> <m|p> [<out_prefix>] [<0|1>]\nfor help file:\nperl",
        " fasta2ms.pl h\n\n";
}

if (exists $ARGV[2]) {
    $out = $ARGV[2];
}

if (exists $ARGV[3]) {
    $verbose = $ARGV[3];
}

if ($ARGV[1] eq "p") {
    $monomorph = 0;
} elsif ($ARGV[1] eq "m") {
    $monomorph = 1;
} else {
    die "USAGE:\nperl fasta2ms.pl <infile> <m|p> [<out_prefix>] [<0|1>]\nfor help file:\nperl fasta2ms.pl h\n\n";
}

open (FASTA, $ARGV[0]) or die "USAGE:\nUSAGE: perl fasta2ms.pl <infile> <m|p> [<out_prefix>] [<0|1>]\nfor help file:",
                              "\nperl fasta2ms.pl h\n\n";
while (<FASTA>) {
    chomp;
    # count fasta sequences
    if (/^>/) {;
        $seqnum++;
    } elsif ($seqnum == 1 && $_ !~ m/^\s*$/) {
        # break nucleotide sequence into temp files, each containing fragment size $sub bp
        $temppos = 1;
        $seq = $_;
        $seqlength = length $seq;
        if ($verbose == 1) {
            print "total sequence length is ", $seqlength, "\n";
        }
        while (length $seq > 0) {
            push(@tempfiles, $temppos); 
            open (TEMP, ">>$out.$temppos.temp");
            $subseq = substr($seq, 0, $sub, "");
            print TEMP $subseq, "\n";
            close (TEMP);
            if ($verbose == 1) {
                print "Position ", $temppos, " to ", $temppos + (length $subseq) - 1, " from sequence ",$seqnum," to ",
                      "temp file\n";
            }
            if (length $seq > 0) {
                $temppos = $temppos + $sub;
            }
        }
    } elsif ($seqnum > 1 && $_ !~ m/^\s*$/) {
        $temppos = 1;
        $seq = $_;
        # check that sequence length of $seqnum 2+ equals the $seqlength of the first sequence
        if (length $seq != $seqlength) {
            die $seqnum + 1, " is unexpected length";
        }
        while (length $seq > 0) {
            open (TEMP, ">>$out.$temppos.temp");
            $subseq = substr($seq, 0, $sub, "");
            print TEMP $subseq, "\n";
            close (TEMP);
            if ($verbose == 1) {
                print "Position ", $temppos, " to ", $temppos + (length $subseq) - 1, " from sequence ",$seqnum," to ", 
                      "temp file\n";
            }
            if (length $seq > 0) {
                $temppos = $temppos + $sub;
            }
        }
    }
}
close (FASTA);
# total number of sequences in input file
$totalseq = $seqnum;
# process each TEMP file
foreach $temppos (@tempfiles) {
    $seqnum = 1;
    open (TEMP, "$out.$temppos.temp");
    while (<TEMP>) {
        if (!/^\s*$/) {
            chomp $_;
            # set $seq equal to nucleotide sequence string
            $subseq = $_;
            $subseqlength = length ($subseq);
            $pos = 1;
            # steps for first nucleotide sequence ($seq)
            if ($seqnum == 1) {
                while ($pos <= $subseqlength) {
                    # Count current positon along $seq
                    # hash key is the position on the sequence, value is the @nucleotide at that position
                    @$position{$pos} = @nucleotide;
                    # remove the first nucleotide from $subseq and set as the $seqnum element of @nucleotide
                    $position{$pos}[$seqnum] = substr($subseq, 0, 1, "");
                    $dump{$pos} = 0;
                    # if nucleotide is equal to "N" remove that position in that key value pair from the hash
                    if ($position{$pos}[$seqnum] eq "N") {
                        $dump{$pos} = 1;
                    }
                    $pos++;
                }
                if ($verbose == 1) {
                    print "Read sequence ", $seqnum, " from ",$temppos," to ",$temppos + $subseqlength - 1, " bp into",
                          " hash of arrays\n";
                }
            } else {
                # steps for remaining sequences ($seqnum > 1). Same as for first sequence except array already exists.
                while ($pos <= $subseqlength) {
                    $position{$pos}[$seqnum] = substr($subseq, 0, 1, "");
                    if ($position{$pos}[$seqnum] eq "N") {
                        $dump{$pos} = 1;
                    }
                    $pos++;
                }
                if ($verbose == 1) {
                    print "Read sequence ", $seqnum, " from ",$temppos," to ",$temppos + $subseqlength - 1," bp into ",
                          "hash of arrays\n";
                }
            }
            $seqnum++;
        }
    }
    close (TEMP);
    unlink ("$out.$temppos.temp");
    # If number of unique elements per array is greater than 2, delete.
    # For p option, if number of unique elements = 1, delete 
    # To do this: go through each key/value pair, create a new hash with the @nucleotide elements as a key
    $pos = 1;
    while ($pos <= $subseqlength) {
        %count = ();
        $seqnum = 1;
        while ($seqnum <= $totalseq) {
            if (defined $position{$pos}[$seqnum]) {
                %count = (%count, $position{$pos}[$seqnum] => 1);
            }
            $seqnum++;
        }
        $countuniq = scalar keys %count;
        if ($countuniq > 2) {
            $dump{$pos} = 1;
        }
        if ($monomorph == 0 && $countuniq == 1) {
            $dump{$pos} = 1;
        }
        if ($countuniq == 2 && $dump{$pos} == 0) {
            $segsites = $segsites + 1;
        }
        # Translate nucleotides to 0s & 1s for non-"N" positions with 2 alleles (p option) or <= 2 alleles (m option).
        if ($dump{$pos} == 0) {
            if (exists $count{A}) {
                $seqnum = 1;
                while ($seqnum <= $totalseq) {
                    $position{$pos}[$seqnum] =~ tr/A[C|G|T]/01/;
                    $seqnum++;
                }
            } elsif (exists $count{C}) {
                $seqnum = 1;
                while ($seqnum <= $totalseq) {
                    $position{$pos}[$seqnum] =~ tr/C[G|T]/01/;
                    $seqnum++;
                }
            } elsif (exists $count{G}) {
                $seqnum = 1;
                while ($seqnum <= $totalseq) {
                    $position{$pos}[$seqnum] =~ tr/GT/01/;
                    $seqnum++;
                }
            } else {
                $seqnum = 1;
                while ($seqnum <= $totalseq) {
                    $position{$pos}[$seqnum] =~ tr/T/0/;
                    $seqnum++;
                }
            }
        }
        $pos++;
    }
    if ($verbose == 1) {
        print "Removed all positions containing more than 2 alleles within subsequence from ",$temppos," to ",$temppos 
        + $subseqlength - 1,"\n";
    }
    # print out hash keys and elements of array in tab delimited format
    open (OUTHEAD, ">>outhead.temp");
    foreach $pos (sort { $a <=> $b } keys %position) {
        if ($dump{$pos} == 0 && $pos <= $subseqlength) {
            my $fraction = ($pos + $temppos - 1)/$seqlength;
            print OUTHEAD sprintf("%e\t", $fraction);
        }
    }
    close (OUTHEAD);
    $seqnum = 1;
    while ($seqnum <= $totalseq) {
        open (OUTDATA, ">>$out.outdata.$seqnum.temp");
        $pos = 1;
        while ($pos <= $subseqlength) {
            if ($dump{$pos} == 0) {
                print OUTDATA $position{$pos}[$seqnum];
            }
            $pos++
        }
        $seqnum++;
    }
    close (OUTDATA);
    if ($verbose == 1) {
        print "Two allele data added to temp file for subsequence\n", $temppos, " to ",$temppos + $subseqlength - 1, 
              "\n";
    }
}
# concatenate processed sequences (OUTDATA) with header (OUTHEAD)
if ($verbose == 1) {
    print "concatenating final data\n";
}
open (OUT, ">>$out.ms");
print OUT "//\nsegsites: ", $segsites, "\npositions:\t";
open (OUTHEAD, "outhead.temp");
while (<OUTHEAD>) {
    print OUT;
}
close (OUTHEAD);
unlink ("outhead.temp");
print OUT "\n";
$seqnum = 1;
while ($seqnum <= $totalseq){
    open (OUTDATA, "$out.outdata.$seqnum.temp");
    while (<OUTDATA>) {
        print OUT;
    }
    close (OUTDATA);
    unlink ("$out.outdata.$seqnum.temp");
    print OUT "\n";
    $seqnum++;
}
close (OUT);
if ($verbose == 1) {
    print "ms file written to $out.ms\n";
}
