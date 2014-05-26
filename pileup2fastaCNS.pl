#!/usr/bin/perl -w

#######################################################################################################################
##                                                                                                                   ##
##  USAGE: perl pileup2fasta.pl -i <infile> -o outfile [-c <coverage threshold>]                                     ##
##  [-s <snp threshold>] [-r <reference base threshold>]                                                             ##
##                                                                                                                   ##
##  For help:  perl pileup2fasta.pl -h                                                                               ##
##                                                                                                                   ##
##  required input:                                                                                                  ##
##  i = infile = file in pileup formatted file (http://samtools.sourceforge.net/pileup.shtml)                        ##
##  o = outfile = file for output.  Will overwrite if file exists.                                                   ##
##                                                                                                                   ##
##  optional input:                                                                                                  ##
##  c = coverage threshold = minimum number of reads coverage at a site.  Must be an integer > 0. (default = 3)      ##
##  s = snp threshold = minimum fraction of reads with the variant (non-reference) nucleotide to record the variant  ##
##      in the consensus sequence.  Must be a fraction > 0.5 and <= 1. (default = 0.9)                               ##
##  r = reference base threshold = minimum fraction of reads with the reference nucleotide to record the reference   ##
##      in the consensus sequence.  Must be a fraction > 0.5 and <= 1. (default = 0.5)                               ##
##                                                                                                                   ##
##  pileup2fasta.pl is designed to take as input a pileup formatted file and transform into fasta formatted          ##
##  consensus sequence(s).  This script ignores insertions and deletions and only considers SNPs.  Where thresholds  ##
##  are not met, an 'N' is recorded in the consensus sequence.                                                       ##
##                                                                                                                   ##
##  Code written by G. John Lazur (2012)                                                                             ##
##  Last updated: 05_23_2014                                                                                         ##
##  The Ohio State University, Dept of Horticulture & Crop Science                                                   ##
##                                                                                                                   ##
#######################################################################################################################

use strict;
use Text::Wrap;
use Getopt::Std;
use vars qw/$opt_i $opt_o $opt_c $opt_s $opt_r $opt_h/;

my ($infile, $outfile, $cov_thr, $snp_thr, $ref_thr) = cmd();
$Text::Wrap::columns = 72;

# test for infile and outfile and open them
open(INFILE, $infile) or die("Error: Cannot open input file\n");
open(OUTFILE, ">$outfile") or die("Error: Cannot open output file\n");

my $entry_name; # sequence identifier for each line
my $pos = 0; # current position in sequence
my $prev_pos; # previous position of sequence that was processed
my $ref_seq; # nucleotide at current position
my $reads; # number of reads
my $snp; # nucleotides read at current position
my $cns_seq; # final sequence
my $checked_name = ''; # compared to $entry_name for processing files with multiple identifiers
my $count = 0; # used to handle initial case to set $checked_name to $entry_name

while(<INFILE>) {
    chomp;
    my @line = split /\t/;
    $prev_pos = $pos;
    $entry_name = $line[0];
    $pos = $line[1];
    $ref_seq = $line[2];
    $reads = $line[3];
    $snp = $line[4];

    # first time through while loop, ensure $checked_name = $entry_name
    if($count == 0) {
        $checked_name = $entry_name;
        $count++;
    }
    # if $checked_name and $entry_name differ, print header and $cns_seq and reset values
    if($checked_name ne $entry_name) {
        print OUTFILE ">$checked_name\n";
        print OUTFILE wrap '', '', $cns_seq, "\n";
        $cns_seq = '';
        $prev_pos = 0;
        $checked_name = $entry_name;
    }
    # if $prev_pos is not position immediately before $pos fill $cns_seq with 'N' until it is
    if($pos != ($prev_pos + 1)) {
        for(($prev_pos + 1) .. ($pos - 1)) {
            $cns_seq .= 'N';
        }
    }
    # remove insertions/deletions
    $snp =~ s/[\+\-]\d+//;
    $_ = $snp;
    # create %counts and fill with ACGT or 'ref' as keys and their respective occurences as values
    my %counts = (tr/aA//, 'A', tr/cC//, 'C', tr/gG//, 'G', tr/tT//, 'T', tr/.,//, 'ref');
    my $max = 0;
    # determine which key has the highest value in %counts
    grep{$max = ($_ > $max) ? $_ : $max} keys %counts;
    # determine which base to add to $cns_seq
    if($reads >= $cov_thr) {    
        if($max / $reads >= $snp_thr) {
            if($counts{$max} eq 'ref') {
                $cns_seq .= $ref_seq;
            } else {
                $cns_seq .= $counts{$max};
            }
        } else {
            $cns_seq .= 'N';
        }
    } elsif($max / $reads >= $ref_thr) {
        if($counts{$max} eq 'ref') {
            $cns_seq .= $ref_seq;
        } else {
            $cns_seq .= 'N';
        }
    } else {
        $cns_seq .= 'N';
    }
}
# print last header and $cns_seq
print OUTFILE ">$checked_name\n";
print OUTFILE wrap '', '', $cns_seq, "\n";

close INFILE;
close OUTFILE;

sub cmd {
    # create variables for command line processing and set default thresholds
    my %options = ();
    my @cmd_line;
    my @list;
    my $cmd_args;
    my $item;
    my $infile = '';
    my $outfile = '';
    my $cov_thr = 3;
    my $snp_thr = 0.90;
    my $ref_thr = 0.50;

    getopts('i:o:c:s:r:h', \%options);

    @list = keys %options;

    foreach $item (@list) {
        if ($item eq "i") {
            $infile = $options{$item};
        } elsif ($item eq "o") {
            $outfile  = $options{$item};
        } elsif ($item eq "c") {
            $cov_thr = $options{$item};
        } elsif ($item eq "s") {
            $snp_thr = $options{$item};
        } elsif ($item eq "r") {
            $ref_thr = $options{$item};
        } elsif ($item eq "h") {
            help();
        }
    }

    @cmd_line = ($infile, $outfile, $cov_thr, $snp_thr, $ref_thr);
    return @cmd_line;
}

sub help {
    print "    USAGE:

    perl pileup2fasta.pl -i PILEUP file (required) 

    -o OUTPUT filename (required) -c NUMBER (default = 3) 
    -s NUMBER (default = 0.90) -r NUMBER (default = 0.50)

    -c minimum base coverage
    -s threshold for SNPs
    -r reference sequence threshold
    ";

    die("\n");
}