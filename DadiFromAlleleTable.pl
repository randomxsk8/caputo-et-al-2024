#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_a $opt_b $opt_o $opt_c);

# Usage
my $usage = "
DadiFromAlleleTable.pl - makes dadi input SFS file from allele table
Copyright (C) 2019 by Jacob A Tennessn 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl DadiFromAlleleTable.pl options
 required:
  -f  (path to) an allele table (numerical genotypes 0, 1, or 2)
  -a  comma-delimited list of sample freq positions [samples start at 1], first group
  -b  comma-delimited list of sample freq positions [samples start at 1], second group
  -o  (path to) output file
optional:
  -c  comma-delimited list of sample freq positions [samples start at 1], third group
";

#############

# command line processing.
getopts('f:a:b:o:c:');
die $usage unless ($opt_f);
die $usage unless ($opt_a);
die $usage unless ($opt_b);
die $usage unless ($opt_o);

my ($genotypes, $grouppos1, $grouppos2, $grouppos3, $outfile);

$genotypes	= $opt_f if $opt_f;
$grouppos1 = $opt_a if $opt_a;
$grouppos2 = $opt_b if $opt_b;
$grouppos3 = $opt_c if $opt_c;
$outfile = $opt_o if $opt_o;

my @g1f = split ",", $grouppos1;

my @g2f = split ",", $grouppos2;

my @g3f;

if (defined $grouppos3) {
  @g3f = split ",", $grouppos3;
}

my %sfs;

my $size1 = 2*scalar(@g1f);

my $size2 = 2*scalar(@g2f);

my $size3 = 2*scalar(@g3f);

my $total = $size1 + $size2 + $size3;

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[0] =~ /Site/) {
      next;
    }
    my $group1alleletotal = 0;
    my $good1samples = 0;
    foreach my $g1 (@g1f) {
      if ($data[$g1] =~ /\d/) {
        $group1alleletotal += $data[$g1];
        $good1samples +=1;
      }
    }
    my $group2alleletotal = 0;
    my $good2samples = 0;
    foreach my $g2 (@g2f) {
      if ($data[$g2] =~ /\d/) {
        $group2alleletotal += $data[$g2];
        $good2samples +=1;
      }
    }
    my $group3alleletotal = 0;
    my $good3samples = 0;
    foreach my $g3 (@g3f) {
      if ($data[$g3] =~ /\d/) {
        $group3alleletotal += $data[$g3];
        $good3samples +=1;
      }
    }
    if (2*($good1samples + $good2samples + $good3samples) >= $total) {
      if (2*($group1alleletotal + $group2alleletotal + $group3alleletotal) > $total) {
        $group1alleletotal = 2*$good1samples - $group1alleletotal;
        $group2alleletotal = 2*$good2samples - $group2alleletotal;
        $group3alleletotal = 2*$good3samples - $group3alleletotal;
      }
      my $combo = "$group1alleletotal\t$group2alleletotal";
      if (defined $grouppos3) {
        $combo = "$group1alleletotal\t$group2alleletotal\t$group3alleletotal";
      }
      if (defined $sfs{$combo}) {
        $sfs{$combo} +=1;
      } else {
        $sfs{$combo} = 1;
      }
    }
}

close (IN);

my @outsfs;

my @outmask;

for (my $a1 = 0; $a1 <= $size1; $a1++) {
  for (my $a2 = 0; $a2 <= $size2; $a2++) {
    if (defined $grouppos3) {
      for (my $a3 = 0; $a3 <= $size3; $a3++) {
        if ((($a1 + $a2 + $a3) == 0)||(($a1 + $a2 + $a3) >= $total)) {
          push @outsfs, 0;
          push @outmask, 1;
        } else {
          my $combo = "$a1\t$a2\t$a3";
          unless (defined $sfs{$combo}) {
            $sfs{$combo} = 0;
          }
          push @outsfs, $sfs{$combo};
          if (2*($a1 + $a2 + $a3) > $total) {
            push @outmask, 1;
          } else {
            push @outmask, 0;
          }
        }
      }
    } else {
      if ((($a1 + $a2) == 0)||(($a1 + $a2) >= $total)) {
        push @outsfs, 0;
        push @outmask, 1;
      } else {
        my $combo = "$a1\t$a2";
        unless (defined $sfs{$combo}) {
          $sfs{$combo} = 0;
        }
        push @outsfs, $sfs{$combo};
        if (2*($a1 + $a2) > $total) {
          push @outmask, 1;
        } else {
          push @outmask, 0;
        }
      }
    }
  }
}

my $outsfs = join " ", @outsfs;

my $outmask = join " ", @outmask;

my $outsize1 = $size1 + 1;

my $outsize2 = $size2 + 1;

my $outsize3 = $size3 + 1;

my $title = "$outsize1 $outsize2 folded";

if (defined $grouppos3) {
  $title = "$outsize1 $outsize2 $outsize3 folded";
}

unless ( open(OUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OUT "$title\n$outsfs\n$outmask";
close (OUT);


