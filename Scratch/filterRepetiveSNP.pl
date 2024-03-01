#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use IO::File;

my $fh     = new IO::File( $ARGV[0], "r" );
my $line   = $fh->getline();

while ($line) {
	chomp($line);
	if ($line) {
		my @tmpAr = split( /\t/, $line );
		my $seq   = $tmpAr[9];
		if ( $line =~ /X0:i:(\d+)/ ) {
			my $count = $1;
			if ( $line =~ /.*X1:i:(\d+)/ ) {
				$count += $1;
			}
			if ( $count <= 1 ) {
				print $line . "\n";
			}
		}
		else {
			print $line . "\n";
		}
	}
	$line = $fh->getline();
}
$fh->close();

