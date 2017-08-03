#!/usr/bin/perl

####################################
# Batch-run CASSIS on SMIPS output #
####################################

use 5.14.2;
use strict;
use warnings;
use File::Basename;
use File::Which;
use Text::CSV;
use List::AllUtils qw(any);

if ( @ARGV < 4 or ( ( @ARGV - 1 ) % 3 != 0 ) or $ARGV[0] !~ m/^\d+$/)
{
	say "usage: perl smips2cassis.pl <CPUs> <SMIPS output> <annotation> <genome> [ <SMIPS ouptput> <annotation> <genome> … ]";
	exit;
}

# path to CASSIS
my $cassis = which("cassis");
$cassis = "/home/twolf/CASSIS/cassis.pl" if (!$cassis);  # location of CASSIS executeable, if not in path

# CASSIS parameters
my $meme = "1.0e+005";
my $fimo = "0.00006";
my $freq = "14";
my $gap  = "2";
my $cpus = shift;

# limit cluster prediction to certain anchor gene types
my @wanted_types = (
	"PKS",
	"NRPS",
	"DMATS",
	"NRPS-like",
	"NRPS-like*",
	"PKS-like",
	"PKS-like*",
	"NRPS-PKS hybrid",
	"NRPS-PKS hybrid-like*",
	);

# organisms to exclude
my %exclude = (
	# "./Ajellomyces_capsulatus_NAm1" => 1
);

my $current_organism = 0;
my $total_organism   = @ARGV / 3;

# input example: aspergillus_flavus_2_proteins.fasta.tsv.anchor_genes.csv
while ( my $input_file = shift )
{
	++$current_organism;

	# check input file aka SMIPS output
	if(-s $input_file)
	{
		say "input file from SMIPS: \"$input_file\"";
	}
	else
	{
		die "check input file \"$input_file\"";
	}

	# get working directory (== dir of input file aka SMIPS output)
	my $working_dir = ( fileparse("$input_file") )[1];
	next if ( $exclude{"$working_dir"} );    # current working dir (== organism) is in the exclude list --> skip it

	# my $annotation;
	# my $genome;
	# 
	# search for genome and annotation file in given dir
	# opendir( my $dir_handle, $working_dir ) or die $!;
	# while ( my $basename = readdir $dir_handle )
	# {
		# if ( !$annotation )
		# {
			# if ( "$basename" =~ m/(_genome_summary_per_gene(_sort)?\.csv)|((_sort)?\.gff3?$)|(_chromosomal_feature(_sort)?\.tab$)|(_Tabelle\.txt)/i )
			# {
				# $annotation = "$working_dir$basename";
				# next;
			# }
		# }
	# 
		# if ( !$genome )
		# {
			# if ( "$basename" =~ m/(_supercontigs\.fa)|(_chromosomes\.fa)|(_Unmasked_Assembly\.fasta)|(_contig\.fa)|(_scaffolds\.fasta)|(_AssemblyScaffolds.fasta)$/i )
			# {
				# $genome = "$working_dir$basename";
				# next;
			# }
		# }
	# }
	# closedir $dir_handle;
	# 
	# die "Could not find annotation file:" if ( !$annotation );
	# die "Could not find genome file:"     if ( !$genome );

	# check annotation file
	my $annotation = shift;
	if (-s $annotation)
	{
		say "annotation file: \"$annotation\"";
	}
	else
	{
		die "check annotation file \"$annotation\"";
	}

	# check genome file
	my $genome = shift;
	if(-s $genome)
	{
		say "genome file: \"$genome\"";
	}
	else
	{
		die "check genome file \"$genome\"";
	}

	# extract anchor, cluster and type from SMIPS output
	my $backbone_csv = Text::CSV->new(
		{
			# allow_loose_quotes => 1   # bad practice in csv format (blame the creator of the input file?)
			sep_char => "\t",
		}
	);

	my @cassis_main_params;
	open( my $backbone_fh, "<", $input_file ) or die $!;
	while ( my $row = $backbone_csv->getline($backbone_fh) )
	{
		next if ( index( $row->[0], "#" ) == 0 );    # skip comments

		my $backbone = $row->[0];
		my $cluster  = $row->[0];                    # don't know the name of the cluster --> call it like its anchor gene
		my $type     = $row->[1];

		push( @cassis_main_params, [ $backbone, $cluster, $type ] ) if ( any{ $type eq $_ } @wanted_types );
	}

	if ( !$backbone_csv->eof )
	{
		$backbone_csv->error_diag;
		die $backbone_csv->error_input;
	}
	close $backbone_fh;

	# run CASSIS with obtained parameters
	for ( 0 .. $#cassis_main_params )
	{
		my $current_backbone = $_ + 1;
		my $total_backbone   = @cassis_main_params;

		my @args = (
			"$cassis",
			"--dir", "$working_dir",
			"--annotation", "$annotation",
			"--genome", "$genome",
			"--meme", $meme,
			"--anchor", "$cassis_main_params[$_][0]",
			"--cluster", "$cassis_main_params[$_][1]",
			"--fimo", $fimo,
			"--frequency", $freq,
			"--prediction", $gap,
			"--verbose",
			"--num-cpus", $cpus
		);

		say "---";
		say "$current_organism/$total_organism $working_dir $current_backbone/$total_backbone $cassis_main_params[$_][0] $cassis_main_params[$_][2]";
		say join(" ", @args);
		say "---";

		system( "perl", @args );
		say "";
	}
}
