use lib "./perl/";
use Summary;
use XML::Simple;
use XML::Writer;
use IO::File;
use List::Util qw[min max];
use Data::Dumper;

my @allreports=Summary::collectreports();
Summary::make_summary (@allreports);




Summary::get_git_hash("../src/version.inc");