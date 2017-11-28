#!/usr/bin/perl

sub process {
	my ($gt,$t,$gl,$l,$g) = @_;
	if ($gt==$t && $gl==$l) {
		print "$g\n";
		return 1;
	}
	if ($gt!=$t) {
		return 0;
	}
	my ($tt,$ll,$acc,$n,$tot,@c) = split /\s+/, $g;
	$n = $c[$l];
	$acc = $n/$tot;
	print "$t $l $acc $n $tot ";
	print (join ' ', @c);
	print "\n";
}

open(F,'<','icpsr.txt') || die "could not open icpsr.txt";
open(G,'<',$ARGV[0]) || die "could not open $ARGV[0]";
$_ = <F>;

$g = <G>;
chomp $g;
($gt,$gl) = split /\s+/, $g;
while(<F>) {
	chomp;
	($l,$t) = split /\s+/, $_;
	$t += 365-10956;
	last if ($t>=$gt);
}
while (1) {
	while (process($gt,$t,$gl,$l,$g)) {
		$_ = <F>;
		($l,$t) = split /\s+/, $_;
		$t += 365-10956;
	}
	$g = <G>;
	last if (eof(G));
	chomp $g;
	($gt,$gl) = split /\s+/, $g;
}
