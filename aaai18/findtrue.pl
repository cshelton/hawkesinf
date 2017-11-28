#!/usr/bin/perl

($nsec,$pid,$nump) = @ARGV;

@fs = `ls cfg*.txt`;
for($i=$pid;$i<@fs;$i+=$nump) {
	$fn = $fs[$i];
	chomp $fn;
	open(F,"./timetst --config=$fn -a 1 -n 1 -r 1 -t $nsec -T $nsec --nopause --noplot --nowritedata |") || die "could not launch timetst";
	$v = "error";
	while(<F>) {
		if (/\(([0-9e+.-]+)\)/) {
			$v = $1;
		}
	}
	close(F);
	print "$fn $v\n";
}
