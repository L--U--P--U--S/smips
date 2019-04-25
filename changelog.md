# 2019-03-08 #
* [web] move SMIPS to a dedicated web server
* [web] update InterProScan to version 5.33-72.0
* [web] HTML output bug should be fixed in this version of InterProScan --> switch from SVG to HTML output

# 2018-03-13 #
* [web] tried to update the server version of InterProScan from 5.19 to 5.27
  * v5.27 does not work because of outdated glibc and libgnutls libraries
  * reverted to v5.19
  * need newer server OS to update InterProScan

# 2017-08-14 #
* Add citing information to readme

# 2017-08-03 #
* [web] move source code, binaries and changelog to github

# 2017-04-24 #
* change MT domain from beeing PKS specific to just beeing an anchor gene domain in general (PKS or NRPS)

# 2017-02-20 #
* [web] fix error resulting in empty input file for protein input option
* [web] properly remove old files from results dir

# 2016-12-23 #
* small adjustment to copyright text

# 2016-07-26 #
* improve error messages
* revise method to remove asterisks ('*') in protein FASTA files

# 2016-07-05 #
* [web] fix a bug with the SMIPS protein input option: uploading an input file with more than one dot in the filename resulted in no email beeing sent to the user

# 2016-04-01 #
* better error handling
* small improvements
