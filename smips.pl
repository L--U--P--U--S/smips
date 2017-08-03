#!/usr/bin/perl

##############################################################
# SMIPS - (S)econdary (M)etablolites from (I)nter(P)ro(S)can #
##############################################################

# Copyright (C) 2015 Leibniz Institute for Natural Product Research and
# Infection Biology -- Hans-Knoell-Institute (HKI)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Dieses Programm ist Freie Software: Sie können es unter den Bedingungen
# der GNU General Public License, wie von der Free Software Foundation,
# Version 3 der Lizenz oder (nach Ihrer Wahl) jeder neueren
# veröffentlichten Version, weiterverbreiten und/oder modifizieren.
#
# Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber
# OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
# Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
# Siehe die GNU General Public License für weitere Details.
#
# Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
# Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
#
# Contact: ekaterina.shelest@hki-jena.de, thomas.wolf@hki-jena.de

use 5.14.2;

use strict;
use warnings;

use Getopt::Long;
use List::AllUtils qw(any);
use Scalar::Util qw(looks_like_number);
use Sort::Naturally;
use IPC::System::Simple qw(system);    # improved behavior of Perl's "system()" call
use File::Which;                       # emulation of "which" command line tool
use File::Basename;
use Bio::SeqIO;
use Data::Dumper;
$Data::Dumper::Terse    = 1;
$Data::Dumper::Sortkeys = 1;
use Text::Autoformat;                  # wrap text automatically (especially helpful for usage text)

use Text::CSV;                         # check if Text::CSV_XS is also installed, for huge gain in reading/writing speed
my $tsv = Text::CSV->new(
	{
		allow_loose_quotes  => 1,      # bad practice in csv format (blame the creator of the input file?)
		allow_loose_escapes => 1,      # bad practice in csv format (blame the creator of the input file?)
		sep_char            => "\t",
		eol                 => "\n",
		binary              => 1
	}
);

######################
# Domain assignments #
######################
sub domain_assignments { }

my %iprs = (

	### PKS ###

	# "general" PKS domains
	IPR015083 => "xPKS",
	#   IPR011141 => "pks", (only in PIRSF, plant PKS type III)
	IPR013601 => "xPKS",

	IPR006765 => "Cyc",

	IPR014043 => "AT",    # also FAS beta subunit in Candida dubliniensis
	IPR020801 => "AT",
	IPR016035 => "AT",    # also FAS beta subunit in Candida dubliniensis
	IPR016036 => "AT",    # sounds like ACP/PP, but is small part of AT

	IPR014030 => "KS",
	IPR014031 => "KS",
	IPR020841 => "KS",
	IPR018201 => "KS",

	IPR013968 => "KR",
	IPR020842 => "KR",

	IPR020807 => "DH",

	IPR020843 => "ER",
	IPR011032 => "ER",
	IPR013154 => "ER",
	IPR013149 => "ER",

	IPR013217 => "MT",
	IPR013216 => "MT",
	IPR004033 => "MT",
	IPR020803 => "MT",
	IPR029063 => "MT",

	### NRPS ###

	# "general" NRPS domains
	IPR010060 => "xNRPS",
	IPR013624 => "xNRPS",
	IPR012728 => "xNRPS",

	IPR010071 => "A",
	IPR000873 => "A",
	IPR020845 => "A",
	IPR025110 => "A",

	IPR001242 => "C",

	### PKS or NRPS ###

	IPR009081 => "PP",
	IPR020806 => "PP",
	IPR006162 => "PP",
	IPR003231 => "PP",
	IPR008278 => "PP",
	IPR004568 => "PP",

	IPR001031 => "TE",
	IPR020802 => "TE",
	IPR010080 => "TE",

	### DMATS ###

	IPR017795 => "Prenyltransferase",
	IPR012148 => "Prenyltransferase",
	IPR017796 => "Prenyltransferase",
	IPR000537 => "Prenyltransferase",
);

my %domains = (
	xPKS => "PKS",
	Cyc  => "PKS",
	AT   => "PKS",
	KS   => "PKS",
	KR   => "PKS",
	DH   => "PKS",
	ER   => "PKS",
#	MT   => "PKS", # can also be found in NRPSs, hence MT does not necessarily result in PKS or hybrid anchor gene

	xNRPS => "NRPS",
	A     => "NRPS",
	C     => "NRPS",

	PP => "PKS/NRPS",
	TE => "PKS/NRPS",
	MT => "PKS/NRPS", # see MT

	Prenyltransferase => "DMATS",
);

sub usage()
{
	my $usage = autoformat(
		"Copyright (C) 2015 Leibniz Institute for Natural Product Research and Infection Biology -- Hans-Knoell-Institute (HKI)\n\n"
		  . "This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions. See the file COPYING, which you should have received along with this program, for details.\n\n"
		  . "     ##############################################################\n"
		  . "     # SMIPS - (S)econdary (M)etablolites from (I)nter(P)ro(S)can #\n"
		  . "     ##############################################################\n\n"
		  . "SMIPS is a tool to predict secondary metabolite (SM) anchor (or backbone) genes in protein sequences. These enzymes play a major role in synthesizing secondary metabolites and their genes are often clustered together with others participating in the same metabolic pathway. Common SM anchor genes are Polyketide synthases (PKS) and Non-ribosomal peptide synthetases (NRPS). The predictions are based on protein domain annotations provided by InterProScan.",
		{ justify => "full", all => 1 }
	  )
	  . "\n\n"
	  . "usage: smips.pl [ options ] [ <interproscan_file.[tsv|out|tab]> | <protein_sequences.fasta> ]\n\n"
	  . autoformat(
		"options:\n\n"
		  . "   --path-to-interproscan, -p  <path>\n\n"
		  . "        Path to InterProScan software, including the file name. Usually the file name is \"interproscan.sh\". If not specified, SMIPS uses the \"which\" command to figure out the location of InterProScan. Only needed if a FASTA file is given.\n\n"
		  . "\n"
		  . "SMIPS accepts one (many) InterProScan output file(s) in TSV or OUT format, or the corresponding JGI files in TAB format, as input. It analyses them for typical SM anchor/backbone genes and writes their domain structure, along with additionl information, to two files in CSV format, extending the file name(s) of the input file(s).\n\n"
		  . "As an alternative, SMIPS accepts protein sequence files in FASTA format. If InterProScan is installed on your machine, SMIPS first calls InterProScan with the given FASTA file(s) and than proceeds with the resulting TSV file(s).\n\n"
		  . "If you don't have a proper InterProScan output at hand, you may download and install InterProScan from",
		{ justify => "full", all => 1 }
	  )
	  . "\n"
	  . "   http://www.ebi.ac.uk/interpro/download.html\n" . "\n"
	  . autoformat(
		"and apply the following command-line to get the desired file (SMIPS will apply the very same command-line, if you run it with a FASTA file)",
		{ justify => "full", all => 1 }
	  )
	  . "\n"
	  . "   interproscan.sh -appl PrositePatterns,SuperFamily,PfamA,SMART,TIGRFAM -i <fasta file with protein sequences> -f tsv,svg -iprlookup\n"
	  . "\n"
	  . autoformat(
"InterProScan HTML output is broken, hence using SVG. However, it's only for graphical (webservice-like) output and optional - you may omit that option.\n\n"
		  . "These mappings inside the SMIPS' source code are used to assign protein domains and anchor genes to this or that category (change them, if you like):",
		{ justify => "full", all => 1 }
	  )
	  . "\n"
	  . "InterPro IPR ID => Anchor gene domain type (%iprs variable):\n"
	  . Dumper( \%iprs ) . "\n"
	  . "Anchor gene domain type => Anchor gene type (%domains variable):\n"
	  . Dumper( \%domains ) . "\n\n"
	  . "Contact: ekaterina.shelest\@hki-jena.de, thomas.wolf\@hki-jena.de";

	# print README to file
	open( my $readme, ">", "README" ) or die "Cannot write to file \"./README\".\n", $!;
	say $readme $usage;
	close $readme;

	# and show it on the screen
	if ( which("tty") )
	{
		if ( system( "tty", "-s" ) == 0 )    # check if stdin is a terminal
		{
			system( "less", "README" );
		}
		else
		{
			system( "cat", "README" );
		}
	}
	else                                     # no tty --> possibly running on windows
	{
		system("type README");
	}

	exit;
}

# read arguments from command line
my $help;
my $path_to_interproscan;

GetOptions(
	"help|h"                   => \$help,
	"path-to-interproscan|p=s" => \$path_to_interproscan
) or die "[ERROR] Please check your command line.\nStopped";

# if no input file given or asked for help
# --> print help information and domain assignments, and exit
usage() if ( !@ARGV or $help );

# if more than one input file --> print current file
my $print_file = 1 if ( @ARGV > 1 );

# run SMIPS for each given input file
while ( my $input_file = shift )
{
	say "\nWorking an \"$input_file\":" if ($print_file);

	# STEP 0a - guess input file type
	my %input_type;
	if ( $input_file =~ m/\.tsv$/ )
	{
		$input_type{interpro} = 1;
	}
	elsif ( $input_file =~ m/\.out$/ )    # same like .tsv from InterPro, just different file extension
	{
		$input_type{interpro} = 1;
	}
	elsif ( $input_file =~ m/\.tab$/ )
	{
		$input_type{jgi} = 1;
	}
	elsif ( $input_file =~ m/\.fasta$/ )
	{
		$input_type{fasta} = 1;
	}
	else
	{
		die
"[ERROR] Sorry, but SMIPS only supports InterProScan files in TSV or OUT format, JGI files in TAB format, or protein sequences in FASTA format.\nStopped";
	}

	# if input file == protein sequences
	if ( $input_type{fasta} )
	{
		# remove stars from sequences, otherwise InterProScan will complain
		say "Removing \'*\' characters from input file, if any. InterProScan doesn't like them.";
		remove_stars($input_file);

		# die if InterProScan not installed or path unkown
		if ( !$path_to_interproscan )
		{
			if ( !which("interproscan") )
			{
				die "[ERROR] Please install InterProScan first or tell me its location via --path-to-interproscan.\nStopped";
			}
			else
			{
				$path_to_interproscan = which("interproscan");
			}
		}

		# run interproscan on protein fasta input file
		say "Running InterProScan on protein sequences file \"$input_file\". This may take several hours. Please be patient …";
		say "\n>>>>> InterProScan BEGIN >>>>>";
		system( $path_to_interproscan, "-appl", "PrositePatterns,SuperFamily,PfamA,SMART,TIGRFAM", "-i", $input_file, "-f", "tsv,svg", "-iprlookup" );
		say "<<<<< InterProScan END <<<<<\n";

		# point to new interproscan output file
		$input_file .= ".tsv";
		$input_type{interpro} = 1;
		say "InterProScan finished with output file \"$input_file\".";
	}

	#	# STEP 0b - create backup and remove some common unescaped double-quotes
	#	my $count_substitutions;
	#	my @lines;
	#
	#	# read original input file
	#	open( my $original, "<", $input_file ) or die "Cannot read from InterProScan file \"$input_file\".\n", $!;
	#	while (<$original>)
	#	{
	#		$count_substitutions += ( $_ =~ s/"G-D-X-G"/G-D-X-G/ig );
	#		$count_substitutions += ( $_ =~ s/"G-D-S-L"/G-D-S-L/ig );
	#		$count_substitutions += ( $_ =~ s/Appr-1"-p/Appr-1-p/ig );
	#		$count_substitutions += ( $_ =~ s/"FY-rich"/FY-rich/ig );
	#		$count_substitutions += ( $_ =~ s/"Winged helix"/Winged helix/ig );
	#		$count_substitutions += ( $_ =~ s/"zincins"/zincins/ig );
	#		$count_substitutions += ( $_ =~ s/"zinc finger"/zinc finger/ig );
	#		push( @lines, $_ );    # save changed and unchanged lines to array
	#	}
	#	close $original;
	#
	#	if ($count_substitutions)
	#	{
	#		# rename original input file
	#		# write lines to new file with name of the original one
	#		rename( $input_file, "$input_file.bak" ) or die "Cannot rename \"$input_file\" to \"$input_file.bak\".\n", $!;
	#		open( my $modified, ">", $input_file ) or die "Cannot write to InterProScan file \"$input_file\".\n", $!;
	#		print $modified @lines;
	#		close $modified;
	#		say "[Warning] Removed $count_substitutions unescaped double quote(s) from input file \"$input_file\"";
	#	}

	# STEP 1 - collect all possible anchor genes
	sub all_backbones { }

	my %all_backbones;

	open( my $in, "<", $input_file ) or die "Cannot read from InterProScan file \"$input_file\".\n", $!;

	while ( my $row = $tsv->getline($in) )
	{
		next if ( "@$row" eq "" or index( $row->[0], "#" ) == 0 );    # comment | emtpy line --> next line

		# some checks to recognize faulty formatted input files
		if ( $input_type{interpro} )
		{
			die "Too little/many columns in line "
			  . $tsv->record_number()
			  . ". Maybe wrong input file format? SMIPS only supports InterProScan files with 11-15 columns.\nStopped"
			  if ( @$row < 11 or @$row > 15 );
			die "Start of protein domain \"$row->[6]\" in line " . $tsv->record_number() . " isn't numeric. Maybe wrong input file format?.\nStopped"
			  if ( !looks_like_number( $row->[6] ) );
			die "Stop of protein domain \"$row->[7]\" in line " . $tsv->record_number() . " isn't numeric. Maybe wrong input file format?.\nStopped"
			  if ( !looks_like_number( $row->[7] ) );
		}
		elsif ( $input_type{jgi} )
		{
			die "Too little/much columns in line "
			  . $tsv->record_number()
			  . ". Maybe wrong input file format? SMIPS only supports JGI files with 10 columns.\nStopped"
			  if ( @$row != 10 );
			die "Start of protein domain \"$row->[7]\" in line " . $tsv->record_number() . " isn't numeric. Maybe wrong input file format?.\nStopped"
			  if ( $row->[7] !~ m/^[0-9|,|-]+$/ );
			die "Stop of protein domain \"$row->[8]\" in line " . $tsv->record_number() . " isn't numeric. Maybe wrong input file format?.\nStopped"
			  if ( $row->[8] !~ m/^[0-9|,|-]+$/ );
			die "No IPR number \"$row->[1]\" in column 2, line " . $tsv->record_number() . ". Maybe wrong input file format?.\nStopped"
			  if ( $row->[1] !~ m/^IPR\d+/ );
		}

		# IPR assigned to row
		# and IPR corresponds to one of the collected anchor gene domains
		if ( $input_type{interpro} and $row->[11] and $iprs{ $row->[11] } )
		{
			my $id          = $row->[0];
			my $ipr         = $row->[11];
			my $domain      = $iprs{$ipr};
			my $domain_type = $domains{$domain};
			my $start       = $row->[6];
			my $stop        = $row->[7];
			my $annotation  = $row->[12];

			die "One of the necessary columns in line " . $tsv->record_number() . " is empty. Maybe wrong input file format?.\nStopped"
			  if ( any { $_ eq "" } ( $id, $ipr, $start, $stop, $annotation ) );

			# gene ID not present yet --> add new ID
			if ( !$all_backbones{$id} )
			{
				$all_backbones{$id} = {
					domains => {
						$domain => [
							{ start => $start, stop => $stop, ipr => { $ipr => 1 }, domain_type => $domain_type, annotation => { $annotation => 1 } }
						]
					}
				};
			}

			# gene ID present, but not for current domain --> add new domain to present ID
			elsif ( !$all_backbones{$id}{domains}{$domain} )
			{
				$all_backbones{$id}{domains}{$domain} =
				  [ { start => $start, stop => $stop, ipr => { $ipr => 1 }, domain_type => $domain_type, annotation => { $annotation => 1 } } ];
			}

			# gene ID and domain already present --> add another instance of that domain
			else
			{
				push(
					@{ $all_backbones{$id}{domains}{$domain} },
					{ start => $start, stop => $stop, ipr => { $ipr => 1 }, domain_type => $domain_type, annotation => { $annotation => 1 } }
				);
			}
		}
		elsif ( ( $input_type{jgi} and $row->[1] and $iprs{ $row->[1] } ) )
		{
			my $id          = $row->[0];
			my $ipr         = $row->[1];
			my $domain      = $iprs{$ipr};
			my $domain_type = $domains{$domain};
			my $start       = $row->[7];
			my $stop        = $row->[8];
			my $annotation  = $row->[2];

			die "One of the necessary columns in line " . $tsv->record_number() . " is empty. Maybe wrong input file format?.\nStopped"
			  if ( any { $_ eq "" } ( $id, $ipr, $start, $stop, $annotation ) );

			# remove trailing "," if present
			$start =~ s/,$//;
			$stop  =~ s/,$//;

			my @starts = split( ",", $start );
			my @stops  = split( ",", $stop );
			die "More/less domain start than stop positions in line " . $tsv->record_number() . ".\nStopped" if ( @starts != @stops );

			for ( 0 .. $#starts )
			{
				my $start = $starts[$_];
				my $stop  = $stops[$_];

				# possible for JGI files --> don't know how to treat --> skip
				if ( $start < 0 )
				{
					warn "[WARNING] Start position \"$start\" in line " . $tsv->record_number(), " is negative. Skipping domain.\n";
					next;
				}
				elsif ( $stop < 0 )
				{
					warn "[WARNING] Stop position \"$stop\" in line " . $tsv->record_number() . " is negative. Skipping domain.\n";
					next;
				}

				# gene ID not present yet --> add new ID
				if ( !$all_backbones{$id} )
				{
					$all_backbones{$id} = {
						domains => {
							$domain => [
								{
									start       => $start,
									stop        => $stop,
									ipr         => { $ipr => 1 },
									domain_type => $domain_type,
									annotation  => { $annotation => 1 }
								}
							]
						}
					};
				}

				# gene ID present, but not for current domain --> add new domain to present ID
				elsif ( !$all_backbones{$id}{domains}{$domain} )
				{
					$all_backbones{$id}{domains}{$domain} =
					  [ { start => $start, stop => $stop, ipr => { $ipr => 1 }, domain_type => $domain_type, annotation => { $annotation => 1 } } ];
				}

				# gene ID and domain already present --> add another instance of that domain
				else
				{
					push(
						@{ $all_backbones{$id}{domains}{$domain} },
						{ start => $start, stop => $stop, ipr => { $ipr => 1 }, domain_type => $domain_type, annotation => { $annotation => 1 } }
					);
				}
			}
		}
	}

	# if we didn't reach the end of the file, something went wrong
	# --> let's see what happened
	if ( !$tsv->eof )
	{
		$tsv->error_diag;
		die $tsv->error_input;
	}
	close $in;

	# print all anchor domains, found in InterProScan output, to file
	open( my $out_all_domains, ">", $input_file . ".anchor_gene_domains.csv" )
	  or die "Cannot write to csv output file \"" . $input_file . ".anchor_gene_domains.csv" . "\".\n", $!;
	$tsv->print( $out_all_domains, [ "# ID", "Domain", "Start", "Stop", "IPR", "Type", "Annotation" ] ) or die $tsv->error_diag();

	for my $id ( nsort( keys %all_backbones ) )
	{
		for my $domain ( sort keys $all_backbones{$id}{domains} )
		{
			for my $domain_occurrence ( @{ $all_backbones{$id}{domains}{$domain} } )
			{
				# only one key for ipr and domain_type at the moment
				$tsv->print(
					$out_all_domains,
					[
						$id,                        $domain,                        $domain_occurrence->{start},
						$domain_occurrence->{stop}, keys $domain_occurrence->{ipr}, $domain_occurrence->{domain_type},
						keys $domain_occurrence->{annotation}
					]
				) or die $tsv->error_diag();
			}
		}
	}

	close $out_all_domains;

	# STEP 2 - select all anchor genes with a minimal set of core domains
	# --> check for existence of core domains
	# --> determine anchor type
	# --> count number of their occurrences
	sub minimal_backbones { }

	my %minimal_backbones;
	my %counts;

	for ( keys %all_backbones )
	{
		my @domains = keys $all_backbones{$_}{domains};    # collect all domains present in the current anchor gene

		# pure PKS
		if (    exists $all_backbones{$_}{domains}{KS}
			and exists $all_backbones{$_}{domains}{AT}
			and exists $all_backbones{$_}{domains}{PP} )
		{
			# or NRPS-PKS hybrid
			if (   exists $all_backbones{$_}{domains}{A}
				or exists $all_backbones{$_}{domains}{C} )
			{
				$all_backbones{$_}{backbone_type} = "NRPS-PKS hybrid";
			}
			else
			{
				$all_backbones{$_}{backbone_type} = "PKS";
			}
		}

		# pure NRPS
		elsif ( exists $all_backbones{$_}{domains}{A}
			and exists $all_backbones{$_}{domains}{PP}
			and exists $all_backbones{$_}{domains}{C} )
		{
			# or NRPS-PKS hybrid
			if (   exists $all_backbones{$_}{domains}{KS}
				or exists $all_backbones{$_}{domains}{AT}
				or exists $all_backbones{$_}{domains}{Cyc}
				or exists $all_backbones{$_}{domains}{KR}
				or exists $all_backbones{$_}{domains}{DH}
				or exists $all_backbones{$_}{domains}{ER}
#				or exists $all_backbones{$_}{domains}{MT}
				)
			{
				$all_backbones{$_}{backbone_type} = "NRPS-PKS hybrid";
			}
			else
			{
				$all_backbones{$_}{backbone_type} = "NRPS";
			}
		}

		# KS-only
		elsif ( exists $all_backbones{$_}{domains}{KS} )
		{
			if ( any { $_ ne "KS" } @domains )
			{
				# or NRPS-PKS hybrid-like
				if (   exists $all_backbones{$_}{domains}{A}
					or exists $all_backbones{$_}{domains}{C} )
				{
					$all_backbones{$_}{backbone_type} = "NRPS-PKS hybrid-like";
				}
				# or PKS-like?
				else
				{
					$all_backbones{$_}{backbone_type} = "PKS-like";
				}
			}
			else
			{
				$all_backbones{$_}{backbone_type} = "KS-only";
			}
		}

		# C-only
		elsif ( exists $all_backbones{$_}{domains}{C} )
		{
			if ( any { $_ ne "C" } @domains )
			{
				# or NRPS-PKS hybrid-like
				if (   exists $all_backbones{$_}{domains}{KS}
					or exists $all_backbones{$_}{domains}{AT}
					or exists $all_backbones{$_}{domains}{Cyc}
					or exists $all_backbones{$_}{domains}{KR}
					or exists $all_backbones{$_}{domains}{DH}
					or exists $all_backbones{$_}{domains}{ER}
#					or exists $all_backbones{$_}{domains}{MT}
					)
				{
					$all_backbones{$_}{backbone_type} = "NRPS-PKS hybrid-like";
				}
				# or NRPS-like
				else
				{
					$all_backbones{$_}{backbone_type} = "NRPS-like";
				}
			}
			else
			{
				$all_backbones{$_}{backbone_type} = "C-only";
			}
		}

		# NRPS/PKS-like … (* similar to definition by SMURF) ################################################
		elsif ( @domains >= 2 )
		{
			# NRPS-PKS hybrid-like* --> at one PKS and one NRPS domain
#			if ( ( any { /^(KS|AT|Cyc|KR|DH|ER|MT)$/ } @domains ) and ( any { /^(A|C)$/ } @domains ) )
			if ( ( any { /^(KS|AT|Cyc|KR|DH|ER)$/ } @domains ) and ( any { /^(A|C)$/ } @domains ) )
			{
				$all_backbones{$_}{backbone_type} = "NRPS-PKS hybrid-like*";
			}
			# PKS-like --> at least two PKS domains
#			elsif ( ( grep { /^(KS|AT|Cyc|KR|DH|ER|MT|PP)$/ } @domains ) >= 2 )
			elsif ( ( grep { /^(KS|AT|Cyc|KR|DH|ER|PP)$/ } @domains ) >= 2 )
			{
				$all_backbones{$_}{backbone_type} = "PKS-like*";
			}
			# NRPS-like --> at least two NRPS domains
			elsif ( ( grep { /^(A|C|PP)$/ } @domains ) >= 2 )
			{
				$all_backbones{$_}{backbone_type} = "NRPS-like*";
			}
		}
		#####################################################################################################

		# trans-AT (exactly one AT domain which might belong to a NRPS that lacks the A domain)
		elsif ( @domains == 1 and $domains[0] eq "AT" )
		{
			$all_backbones{$_}{backbone_type} = "trans-AT";
		}

		# pure DMATS
		elsif ( exists $all_backbones{$_}{domains}{Prenyltransferase} )
		{
			# or DMATS-NRPS-PKS hybrid
			if (
				(
					   exists $all_backbones{$_}{domains}{KS}
					or exists $all_backbones{$_}{domains}{AT}
					or exists $all_backbones{$_}{domains}{Cyc}
					or exists $all_backbones{$_}{domains}{KR}
					or exists $all_backbones{$_}{domains}{DH}
					or exists $all_backbones{$_}{domains}{ER}
#					or exists $all_backbones{$_}{domains}{MT}
				)
				and (  exists $all_backbones{$_}{domains}{A}
					or exists $all_backbones{$_}{domains}{C} )
			  )
			{
				$all_backbones{$_}{backbone_type} = "DMATS-NRPS-PKS hybrid";
			}

			# or DMATS-PKS hybrid
			elsif (exists $all_backbones{$_}{domains}{KS}
				or exists $all_backbones{$_}{domains}{AT}
				or exists $all_backbones{$_}{domains}{Cyc}
				or exists $all_backbones{$_}{domains}{KR}
				or exists $all_backbones{$_}{domains}{DH}
				or exists $all_backbones{$_}{domains}{ER}
#				or exists $all_backbones{$_}{domains}{MT}
				)
			{
				$all_backbones{$_}{backbone_type} = "DMATS-PKS hybrid";
			}
			# or DMATS-NRPS hybrid
			elsif (exists $all_backbones{$_}{domains}{A}
				or exists $all_backbones{$_}{domains}{C} )
			{
				$all_backbones{$_}{backbone_type} = "DMATS-NRPS hybrid";
			}
			else
			{
				$all_backbones{$_}{backbone_type} = "DMATS";
			}
		}

		# current gene got an anchor type --> minimal set of core domains present
		if ( $all_backbones{$_}{backbone_type} )
		{
			$counts{ $all_backbones{$_}{backbone_type} }++;
			$minimal_backbones{$_} = $all_backbones{$_};
		}
	}

	# STEP 3 - sort and merge domains of "true" anchor genes
	# --> sort domains by position
	# --> merge domains which belong together
	sub sort_and_merge { }

	my %minimal_backbones_sorted;
	for my $id ( keys %minimal_backbones )
	{
		# make all instances of all domains of a gene sortable by their start and stop positions
		# and sort them by position, huh
		my @domains;
		for my $domain ( keys $minimal_backbones{$id}{domains} )
		{
			push( @domains,
				map { { start => $_->{start}, stop => $_->{stop}, name => $domain, content => $_ } } @{ $minimal_backbones{$id}{domains}{$domain} } );
		}
		@domains =
		  map { { name => $_->{name}, %{ $_->{content} } } } sort { $a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop} } @domains;

		# set anchor type and prepare domains array
		$minimal_backbones_sorted{$id} = { backbone_type => $minimal_backbones{$id}{backbone_type}, domains => [] };

		# merge all domains $j=2…n, which are of the same type as domain $i=1…n-1, into domain $i
		for ( my $i = 0 ; $i <= $#domains ; $i++ )
		{
			my $jump = $i;

			for ( my $j = $i + 1 ; $j <= $#domains ; $j++ )
			{
				# overlapping domains of same type
				# or adjecent domains of same type with completely different IPR(s)?
				# --> merge both into first domain
				if ( $domains[$i]{name} eq $domains[$j]{name}
					and ( $domains[$i]{stop} >= $domains[$j]{start} or !any { $domains[$j]{ipr}{$_} } keys $domains[$i]{ipr} ) )
				{
					my %merged_iprs        = ( %{ $domains[$i]{ipr} },        %{ $domains[$j]{ipr} } );
					my %merged_annotations = ( %{ $domains[$i]{annotation} }, %{ $domains[$j]{annotation} } );

					my $start;
					my $stop;
					( $domains[$i]{start} < $domains[$j]{start} ) ? ( $start = $domains[$i]{start} ) : ( $start = $domains[$j]{start} );
					( $domains[$i]{stop} > $domains[$j]{stop} )   ? ( $stop  = $domains[$i]{stop} )  : ( $stop  = $domains[$j]{stop} );

					$domains[$i] = {
						name        => $domains[$i]{name},
						ipr         => \%merged_iprs,
						start       => $start,
						stop        => $stop,
						domain_type => $domains[$i]{domain_type},
						annotation  => \%merged_annotations
					};
					$jump = $j;
				}
				else
				{
					$jump = $j - 1;
					last;    # end inner for loop if adjacent domains are not of same type (anymore)
				}
			}

			push( @{ $minimal_backbones_sorted{$id}{domains} }, $domains[$i] );    # save type and content of (maybe merged) domain
			$i = $jump;    # next $i will be current $j ($jump=$j-1) OR one position in advance of current $j ($jump=$j)
		}
	}

	# STEP 4 - print results file
	# --> basically, content of %minimal_backbones in csv format
	sub print { }

	open( my $out_backbones, ">", $input_file . ".anchor_genes.csv" )
	  or die "Cannot write to csv output file \"" . $input_file . ".anchor_genes.csv" . "\".\n", $!;

	# MS Excel is just stupid!!!
	# (1) Either uncomment to following line, to tell excel which column separator to use
	# (but this will add a non-comment otherwise senseless first line to the output file)
	# OR (2) change the file extension to ".txt".
	# Without those cheats, excel doesn't recognize the tabs and will show each line in a single column.
	# While applying one of those cheats, Excel will "see" now (why? nothing changed!) the tabs.
	#
	# say $out "sep=\t";

	# first, some statistics and explanations
	if ( any { /NRPS|PKS|DMATS/ and !/like/ } keys %counts )
	{
		say $out_backbones "## NRPS, PKS, and DMATS genes:";
		say $out_backbones "# $_ $counts{$_}" for ( sort grep { /NRPS|PKS|DMATS/ and !/like/ } keys %counts );
		say $out_backbones "#";
	}
	if ( any { /NRPS|PKS/ and /like/ and !/hybrid/ } keys %counts )
	{
		say $out_backbones "## NRPS- and PKS-like genes:";
		say $out_backbones "# $_ $counts{$_}" for ( sort grep { /NRPS|PKS/ and /like/ and !/hybrid/ } keys %counts );
		say $out_backbones "#";
	}
	if ( any { /only|trans/ } keys %counts )
	{
		say $out_backbones "## Single domains (genes where only AT, C, or KS are reported):";
		say $out_backbones "# $_ $counts{$_}" for ( sort grep { /only|trans/ } keys %counts );
		say $out_backbones "#";
	}

	say $out_backbones "# -only: incomplete anchor genes, with only KS or C domain(s)"          if ( any { /-only$/ } keys %counts );
	say $out_backbones "# -like: incomplete anchor genes, with at least one KS and/or C domain" if ( any { /-like$/ } keys %counts );
	say $out_backbones
	  "# -like*: incomplete anchor genes, with (at least) two typical PKS and/or NRPS domains (similar to definition by the SMURF tool)"
	  if ( any { /-like\*$/ } keys %counts );
	say $out_backbones "# trans-AT: possible trans-AT domain, which is used by an somewhere else located PKS (which is lacking the AT domain)"
	  if ( any { /^trans-AT$/ } keys %counts );
	say $out_backbones "#";
	say $out_backbones "# A: Adenylation, AMP-dependent synthetase/ligase"
	  if ( any { /^A$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# AT: Acyl transferase" if ( any { /^AT$/ } map  { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# C: Condensation"      if ( any { /^C$/ } map   { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# Cyc: Cyclase"         if ( any { /^Cyc$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# DH: Dehydratase"      if ( any { /^DH$/ } map  { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# ER: Enoylreductase, GroES-like, Alcohol dehydrogenase"
	  if ( any { /^ER$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# KR: Keto reductase"                         if ( any { /^KR$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# KS: Beta-ketoacyl synthase"                 if ( any { /^KS$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# MT: Methyltransferase"                      if ( any { /^MT$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# PP: Phosphopantetheine"                     if ( any { /^PP$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# TE: Thioesterase, Thioester reductase-like" if ( any { /^TE$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# xNRPS: typical NRPS domain with yet unkown function"
	  if ( any { /^xNRPS$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "# xPKS: typical PKS domain with yet unkown function"
	  if ( any { /^xPKS$/ } map { keys $_->{domains} } values %minimal_backbones );
	say $out_backbones "#";
	say $out_backbones "# \'X-Y\': depiction of distinct domains";
	say $out_backbones "# \'X:Y\': depiction of overlapping domains";
	say $out_backbones "#";

	# next, the table head
	say $out_backbones "## " . join( "\t", "Gene ID", "Anchor gene type", "Domain arrangement", "Domains" );

	# finally, the content of %minimal_backbones
	# NATURALLY sorted by the gene ID (or protein ID or whatever …)
	for my $id ( nsort( keys %minimal_backbones_sorted ) )
	{
		my @domain_details;
		my $domain_arrangement;
		my $backbone_type;

		# stringify the domain details
		for ( @{ $minimal_backbones_sorted{$id}{domains} } )
		{
			push(
				@domain_details,
				join( ", ",
					"name=$_->{name}",   "type=$_->{domain_type}",
					"start=$_->{start}", "stop=$_->{stop}",
					"ipr=" . join( "|", sort keys $_->{ipr} ), "annotation=" . join( " | ", sort keys $_->{annotation} ) )
			);
		}

		# gene/protein has more than one domain?
		if ( @{ $minimal_backbones_sorted{$id}{domains} } > 1 )
		{
			for ( 0 .. $#{ $minimal_backbones_sorted{$id}{domains} } - 1 )
			{
				my $delimiter = "-";

				# overlapping domains of different types?
				# WARNING: domains overlapping with more than one other domain will not be recognized and marked with ":"
				# --> only the first overlap will have a proper representation in the output
				if ( $minimal_backbones_sorted{$id}{domains}[$_]{stop} >= $minimal_backbones_sorted{$id}{domains}[ $_ + 1 ]{start} )
				{
					warn "[WARNING] Overlapping domains in protein \"$id\": "
					  . "<$minimal_backbones_sorted{$id}{domains}[$_]{start}"
					  . " [$minimal_backbones_sorted{$id}{domains}[$_]{name}] "
					  . "$minimal_backbones_sorted{$id}{domains}[$_]{stop}>" . " and "
					  . "<$minimal_backbones_sorted{$id}{domains}[$_+1]{start}"
					  . " [$minimal_backbones_sorted{$id}{domains}[$_+1]{name}] "
					  . "$minimal_backbones_sorted{$id}{domains}[$_+1]{stop}>\n";
					$delimiter = ":";
				}
				$domain_arrangement .= $minimal_backbones_sorted{$id}{domains}[$_]{name} . $delimiter;
			}
			$domain_arrangement .= $minimal_backbones_sorted{$id}{domains}[-1]{name};
		}

		# gene/protein has only one domain? --> quite easy domain "arrangement"
		else
		{
			$domain_arrangement = $minimal_backbones_sorted{$id}{domains}[0]{name};
		}

		# the anchor type is already there
		$backbone_type = $minimal_backbones_sorted{$id}{backbone_type};

		# print line to file
		$tsv->print( $out_backbones, [ $id, $backbone_type, $domain_arrangement, @domain_details ] ) or die $tsv->error_diag();
	}
	close $out_backbones;

	# also print short results to console
	if ( any { /NRPS|PKS|DMATS/ and !/like/ } keys %counts )
	{
		say "NRPS, PKS, and DMATS genes:";
		say "$_ $counts{$_}" for ( sort grep { /NRPS|PKS|DMATS/ and !/like/ } keys %counts );
	}
	if ( any { /NRPS|PKS/ and /like/ and !/hybrid/ } keys %counts )
	{
		say "NRPS- and PKS-like genes:";
		say "$_ $counts{$_}" for ( sort grep { /NRPS|PKS/ and /like/ and !/hybrid/ } keys %counts );
	}
	if ( any { /only|trans/ } keys %counts )
	{
		say "Single domains (genes where only AT, C, or KS are reported):";
		say "$_ $counts{$_}" for ( sort grep { /only|trans/ } keys %counts );
	}
}

sub remove_stars
{
	my $infile = shift;
	my @temp_seq;    # keep (changed) sequences in memory before writing to output

	my $in = Bio::SeqIO->new( -file => "<$infile" );
	while ( my $seq_obj = $in->next_seq() )
	{
		my $seq  = $seq_obj->seq();
		my $desc = $seq_obj->desc();

		# check if protein input file actually contains protein sequences only
		die "[ERROR] Sequence \"$desc\" is not a protein sequence. Please check your input file.\nStopped" if ( $seq_obj->alphabet() ne "protein" );

		# remove all asterisks (if any)
		$seq =~ s/\*//g;

		$seq_obj->seq($seq);
		$seq_obj->desc($desc);

		push( @temp_seq, $seq_obj );
	}

	my $out = Bio::SeqIO->new( -file => ">$infile" );    # overwrite infile
	$out->write_seq($_) for @temp_seq;
}
