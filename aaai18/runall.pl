#!/usr/bin/perl

($nsec,$pid,$nump) = @ARGV;

@fs = `ls cfg*.txt`;
for($i=$pid;$i<@fs;$i+=$nump) {
	$fn = $fs[$i];
	chomp $fn;
	print("./timetst --config=$fn --nopause\n");
	system("./timetst --config=$fn --nopause");
}
