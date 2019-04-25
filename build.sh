#!/bin/bash

# Perl modules
# sudo apt-get install --no-upgrade libperl-dev cpanminus
# sudo cpanm --skip-satisfied PAR::Packer LWP::UserAgent Bio::Seq Data::Dumper File::Basename File::Copy::Recursive File::Which Getopt::Long IPC::System::Simple List::AllUtils Parallel::ForkManager Roman Scalar::Util Sort::Naturally Sys::CpuAffinity Text::CSV Text::CSV_XS Text::CSV_PP

# build SMIPS binary
sed -i 's/usage: smips\.pl/usage: smips/' smips.pl
perl smips.pl > /dev/null || exit
pp -M Text::CSV_XS -M Text::CSV_PP -o smips smips.pl || exit

# pack SMIPS binary
rm "SMIPS Linux 64bit.7z"
7z a "SMIPS Linux 64bit.7z" smips README COPYING

# clean up
rm smips
sed -i 's/usage: smips/usage: smips\.pl/' smips.pl
