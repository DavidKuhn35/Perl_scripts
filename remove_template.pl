#!/usr/bin/perl -w
open (SNP, "$ARGV[0]") || die "Can't open XML infile";
	while (<SNP>)
{
	unless (m/^<template>|^<properties/)  	{	#If line doesn't have <template> or <properties> at the start
open (OUTFILE, ">>$ARGV[0]\_notemplate\_noproperties") || die "Can't open outfile";
print OUTFILE "$_";
			}
}
