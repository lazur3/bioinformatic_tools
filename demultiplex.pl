#!/usr/bin/perl -w

#######################################################################################
##                                                                                   ##
##  USAGE: perl demultiplex.pl <index> <infile1> <infile2> ... <infilex>             ##
##                                                                                   ##
##  index = list of indices                                                          ##
##                                                                                   ##
##  infile1 = fastq file (read 1) with index in header                               ##
##                                                                                   ##
##  infile2 = fastq file (read 2) with index in header                               ##
##                                                                                   ##
##  infilex = fastq file (read x) with index in header                               ##
##                                                                                   ##
##  demultiplex.pl demultiplexes (splits files acording to index).  It takes input   ##
##  from multiplexed illumina reads and a list of indices.  Output is a fastq file   ##
##  for each index and infile. Outfiles are named <index dot read file number>.      ##
##  For example, if a read in the third file after the index file in user input      ##
##  matches the index ACGT, the read would be found in the outfile "ACGT.3".         ##
##                                                                                   ##
##  example header:                                                                  ##
##      @HWI-ST1052:74:D10TLACXX:5:1101:2540:91596 1:N:0:CTTGGAA                     ##
##                                                                                   ##
##  index must follow the third ':' after the space in the header.                   ##
##                                                                                   ##
##  Code written by G. John Lazur & Leah K. McHale (2012)                            ##
##  Last updated: 06_03_2014                                                         ##
##  The Ohio State University, Dept of Horticulture & Crop Science                   ##
##                                                                                   ##
#######################################################################################

use strict;

# default hash to store each record
my %file = ("UNIQ" => "",
            "BAR" => "",
            "SEQ" => "",
            "PLUS" => "+",
            "QUAL" => "");
# array to hold the barcodes
my @barcodes = ();
# parallel array for filenames
my @filenames = ();
# shift @ARGV to give $infile the file containing list of barcodes
my $infile = shift(@ARGV);

# open file containing list of all barcodes
open INFILE, $infile or die "Unable to open barcode file.";
# read in barcode file line by line
while(<INFILE>) {
    chomp;
    # add barcode to array of filenames
    push @filenames, $_;
    # split on each character to make a character array for regular expressions
    my @line = split('');
    # push each character array into barcodes array
    push @barcodes, [@line];
}
close INFILE;

# $len is length of the barcodes used for generating regular expressions
my $len = length($filenames[0]);
# call regular expression generating subroutine
my $regex = expr();

#determine number of files being processed
my $num_files = scalar @ARGV;
#delete all files in current directory that will be created by this program (BARCODE.ARGVPOSITION)
foreach my $file (@filenames) {
    my $tempext = 1;
    while($tempext <= $num_files) {
        my $tempfile = $file . "." . $tempext;
        if(-e $tempfile) {
            unlink($tempfile);
        }
        $tempext++;
    }
}

#process the files
my $file_ext = 1;
while($file_ext <= scalar @ARGV) {
    open INFILE, $ARGV[$file_ext - 1] or die "Unable to open file number $file_ext";
    # there are four lines to each record, start with line 1
    my $line_number = 1;
    while(<INFILE>) {
        chomp;
        if($line_number == 1) {
            # first line, retrieve unique identifier and barcode and update line number
            my @line1 = split(' ');
            my @line2 = split (':', $line1[1]);
            $file{"UNIQ"} = $line1[0] . " " . $line2[0] . ":" . $line2[1] . ":" . $line2[2];
            $file{"BAR"} = $line2[scalar @line2 - 1];
            $line_number++;
        } elsif($line_number == 2) {
            # second line, retrieve sequence and update line number
            $file{"SEQ"} = $_;
            $line_number++;
        } elsif($line_number == 3) {
            # third line, retrieve nothing and update line number
            $line_number++;
        } elsif($line_number == 4) {
            # fourth line, retrieve quality and update line number
            $file{"QUAL"} = $_;
            my $sub = 0;
            # while there are still more barcodes to test the data against
            while($sub < scalar @barcodes) {
                # if current data member's barcode matches barcode or all but one character of barcode
                if($file{"BAR"} =~ m/$$regex[$sub]/) {
                    # obtain filename from @filenames position parallel to current @barcodes position
                    my $filewriter = $filenames[$sub] . "." . $file_ext;
                    # create or open outfile
                    open OUTFILE, ">>", $filewriter;
                    print OUTFILE $file{"UNIQ"}, "::", $file{"BAR"}, "\n", $file{"SEQ"}, "\n", $file{"PLUS"}, "\n", $file{"QUAL"}, "\n";
                    close OUTFILE;
                }
                # increment current subscript
                $sub++;
            }        
            # reset hash to default values
            $line_number = 1;
            %file = ("UNIQ" => "",
                     "BAR" => "",
                     "SEQ" => "",
                     "PLUS" => "+",
                     "QUAL" => "");
        }
    }
    # increment current file extension
    $file_ext++;
}

# subroutine to build regular expression based on length and number of barcodes
sub expr {
    # $line gathers characters to create a full line
    my $line = "";
    # $final_exp gathers lines to create full regular expression
    my $final_exp = "";
    # @regex gathers regular expression to achive full set for all barcodes
    my @regex = ();
    my $s = 0;
    # create a set of regular expressions for each barcode
    while($s < scalar @barcodes) {
        my $count1 = 0;
        my $count2 = 0;
        while($count1 < $len) {
            while($count2 < $len) {
                # each line begins with '^'
                $line .= "^";
                # slide the wildcard character 1 pos per line
                for(my $y = 0; $y < $len; $y++) {
                    if($y == $count2) {
                        $line .= ".";
                    } else {
                        $line .= "$barcodes[$s][$y]";
                    }
                }
                $count2++;
                # each line ends with '$'
                $line .= "\$";
                # include '|' to end each but the last line
                if($count2 < $len) {
                    $line .= '|' . "\n";
                }        
                $final_exp .= $line;
                $line = '';
            }
            $count1++;
        }
        # add completed line to expression and reset $final_exp
        push @regex, $final_exp;
        $final_exp = "";
        $s++;
    }
    # return reference to regular expression array
    return \@regex;
}
