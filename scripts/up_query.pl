# Reference for uniprot query fields can be found here: http://www.uniprot.org/help/programmatic_access

use strict;
use warnings;
use LWP::UserAgent;

my $file = shift(@ARGV);
open(FILE, $file) || die "cannot open $file";
while(<FILE>) {
         chomp;
         my $query_term = $_;
         # query UniProt for the given identifier and retrieve tab-delimited format
         system ("wget -O output -q  \"http://www.uniprot.org/uniprot/?query=id:$query_term&format=tab&columns=id,reviewed,length,organism,organism-id,mass,database(dbSNP),comment(ALTERNATIVE%20PRODUCTS),annotation%20score,existence,protein%20names\"");
         open(WGET, "output") || die "cannot open output";
         while(<WGET>) {
                 print $_ if (/[A-Z]/ && !/Entry/); # do not print header
         }
         close(WGET);
}
close(FILE);
