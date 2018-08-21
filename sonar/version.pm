
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
			  print STDERR "$message";
			  open LOG, ">>$version::LOGFILE";
			  print LOG localtime() . " -- DIED WITH ERROR:\n$message\n";
			  close LOG;
			  exit("179");
                         }

}


sub logCmdLine() {

    my $command = shift;
    foreach(@_) { $command .= /(\s|\*)/ ? " \'$_\'" : " $_"; }
    my $VERSION = `git --git-dir $FindBin::Bin/../../.git --work-tree=$FindBin::Bin/../ describe --always --dirty --tags`;
    chomp $VERSION;

    our $LOGFILE = "SONAR_command_history.log";
    if (-d "output/logs") {
	$LOGFILE = "output/logs/command_history.log";
    }

    my $LOGMESSAGE = "\n" . localtime() . " -- SONAR $VERSION run with command:\n\t$command\n";
    print STDERR $LOGMESSAGE;
    
    my $check = open LOG, ">>$LOGFILE";
    if ($check) {
	print LOG $LOGMESSAGE;
	close LOG;
	
	our $PRINT_LOGS = 1;
    } else {
	warn("Directory appears to be read-only; command line and output will not be saved\n");
    }

}

END {

    if ( -defined $version::PRINT_LOGS ) {

	if ($? != 179) { #error code 179 indicates a die command which was already logged by $SIG{__DIE__} above
	    open LOG, ">>$version::LOGFILE";
	    my $status =  $? == 0 ? "finished successfully" : "exited unexpectedly with error code $?";
	    print LOG localtime() . " -- Program $status\n";
	    close LOG;
	}

    }
}

1;

