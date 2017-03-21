
#!/usr/bin/env perl

package version;
use strict;

use vars '@ISA', '@EXPORT', '$NAME', '$VERSION', '$DATE', '$AUTHOR';
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(logCmdLine);


BEGIN {

    $SIG{__DIE__} = sub {
	                  my $message = shift;
			  print "$message";
			  open LOG, ">>output/logs/command_history.log";
			  print LOG localtime . " -- DIED WITH ERROR:\n$message\n";
			  close LOG;
			  exit("179");
                         }

}


sub logCmdLine() {

    my $command = shift;
    foreach(@_) { $command .= /\s/ ? " \'$_\'" : " $_"; }
    my $VERSION = `git -C $FindBin::Bin describe --always --dirty --tags`;
    chomp $VERSION;

    if (-d "output/logs") {
	open LOG, ">>output/logs/command_history.log" or die "Can't write to output/logs/command_history.log: $!\n\n";
	print LOG "\n";
	print LOG localtime . " -- SONAR $VERSION run with command:\n\t$command\n";
	close LOG;
	
	our $PRINT_LOGS = 1;
    } else {
	warn "SONAR log directory not found; command line and output will not be saved\n";
    }

}

END {

    if ( -defined $version::PRINT_LOGS ) {

	if ($? != 179) { #error code 179 indicates a die command which was already logged by $SIG{__DIE__} above
	    open LOG, ">>output/logs/command_history.log";
	    my $status =  $? == 0 ? "finished successfully" : "exited unexpectedly with error code $?";
	    print LOG localtime . " -- Program $status\n";
	    close LOG;
	}

    }
}

1;

