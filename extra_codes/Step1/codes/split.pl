#!/usr/bin/perl

use strict; 
use warnings;

my $word_size = 9000;

$/="\n>"; # read sequence blocks
while (<>) {
    s/>//g;
    my ($id, @seq) = split (/\n/, $_);
    my $seq = join "", @seq;
    while  ($seq) {
        my $sub_seq = substr($seq, 0, $word_size);
        substr($seq, 0, $word_size) = '';
        print ">$id\n$sub_seq\n";
    }
}
