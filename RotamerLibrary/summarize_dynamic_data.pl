#!/usr/bin/perl
BEGIN { push @INC, "/users/scouras/code/RotLib" }
use warnings;
use strict;
use Data::Dumper;
use Rotamer;
use Utility;




my ( $ROOT, $SUMM, @LIBRARY_ORDER, %LIBRARIES );

my $CONFIG_FILE = $ARGV[0] || "/users/scouras/code/RotLib/native_810.dynamic.dat";

if ( not -e $CONFIG_FILE ) { 
  die "Config file '$CONFIG_FILE' not found.";
}

my $CONFIG_CONTENTS = `cat $CONFIG_FILE`;
eval $CONFIG_CONTENTS;
if ( not -e $SUMM ) { mkdir $SUMM or die "Couldn't mkdir $SUMM. $!" }

#print Dumper ( \%LIBRARIES );

my %AA = %Utility::AA;

my %CONFIG = (
  rotamer_mode        => "dynamic",
  home_dir            => "/users/scouras",
  root_dir            => "/users/scouras/rotamers",
  overwrite_rotamers  => 0,
  overwrite_dihedrals => 0,
);



#============================================================= Files

my $FILE_IN_STATE_COUNTS  = "state_counts.dat";
my $FILE_OUT_ASS          = "$SUMM/dyn.assignable.dat";
my $FILE_OUT_SIG          = "$SUMM/dyn.significant.dat";

my $FILE_IN_WAIT_COUNTS   = "dynamics_values.mean_counts.dat";
my $FILE_OUT_WAIT_COUNTS  = "$SUMM/dyn.wait_counts.dat";


#======================================== Load Residue Libraries


my $ROT = RotLib::Rotamer -> Initialize ( \%CONFIG, \%AA );

#my @OFFICIAL_AAS = qw(LEU);
my @OFFICIAL_AAS = qw(ARG ASN ASP CYH CYS GLN GLU HID HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL);
my @ALL_AAS = ('ALA', 'GLY', @OFFICIAL_AAS);

my %RESIDUE_CONVERSION = (
  HID => 'HIS',
  HIE => 'HIS',
  CYH => 'CYS',
  CYS => 'CYH',
);


#================================== Load Residue Library Directories

my @LIB_PAIRS = ();
foreach my $i ( 0..$#LIBRARY_ORDER ) { 
  foreach my $j ( $i+1..$#LIBRARY_ORDER ) { 
    push @LIB_PAIRS, [$LIBRARY_ORDER[$i], $LIBRARY_ORDER[$j]];
  }
}

while ( my ( $LIB_NAME, $LIB ) = each %LIBRARIES ) {
  foreach my $res_orig ( @ALL_AAS ) { 
    my $res_virt = $res_orig;

    #=== Find a residue data directory, possibly using the virtual residue name
    my $res_dir = (grep { /$res_orig/ } @{$LIB->{'dir_globs'}})[0];
    if ( not defined $res_dir ) { 
      if ( exists $RESIDUE_CONVERSION{$res_orig} ) {
        $res_virt = $RESIDUE_CONVERSION{$res_orig};
        $res_dir = (grep { /$res_virt/ } @{$LIB->{'dir_globs'}})[0];
        if ( not defined $res_dir ) { 
          die "Couldn't find residue directory for library $LIB_NAME residue $res_orig $res_virt";
        }
      }
    }

    $LIB->{'virts'}{$res_orig} = $res_virt;
    $LIB->{'dirs' }{$res_orig} = $res_dir;

  }
}




my %MINOR_INTS  = ();
my %MINOR_NAMES = ();
my %MINOR_LOOKUP = ();

foreach my $res ( @OFFICIAL_AAS ) { 
  
  foreach my $major ( 1..$AA{$res}{'total_rotamers'} ) { 
    my @minors = $ROT->Get_Rotamer_Minors($res, $major-1);
    my @names  = $ROT->Get_Rotamer_Minor_Names($res, [ map{$_-1} @minors ] );
    $MINOR_INTS {$res}[$major-1] = \@minors;
    $MINOR_NAMES{$res}[$major-1] = \@names;
  }

  $MINOR_LOOKUP{$res} = $ROT->Get_Rotamer_Minor_Name_Options($res);
}


#  485  grep 0.5000 */state_counts.dat > ../summary/counts.assignable.dat
#  486  grep 0.1000 */state_counts.dat > ../summary/counts.significant.dat 
#  490  grep 0.1000 */state_counts.dat > ../summary/counts.significant.dynamic.dat 
#  491  grep 0.5000 */state_counts.dat > ../summary/counts.assignable.dynamic.dat 
#  500  cat */dynamics_values.mean_counts.dat > ../summary/dynamics_values.mean_counts.dat

####################################################################
#                                  ASSIGNABLE AND SIGNIFICANT COUNTS
####################################################################

sub Get_Assignability { 
  my $file = $_[0];
  my $THRESHOLD = "0.5000";
  my $line = `grep $THRESHOLD $file`;
  my ($t, $d0, $d1, $dx) = split /\s+/, $line;
  return $d1;
}

sub Get_Significant { 
  my $file = $_[0];
  my $THRESHOLD = "0.1000";
  my $line = `grep $THRESHOLD $file`;
  my ( $t, @data ) = split /\s+/, $line;
  my $sum = 0;
  map { $sum += ( $_ * $data[$_] ) } 0..$#data;
  return $sum;
}

open SC_SIG, ">$FILE_OUT_SIG" or die "Couldn't open dyn sig, '$FILE_OUT_SIG'. $!";
open SC_ASS, ">$FILE_OUT_ASS" or die "Couldn't open dyn ass, '$FILE_OUT_ASS'. $!";

print SC_SIG "RESIDUE\t" . (join "\t", @LIBRARY_ORDER) . "\n";
print SC_ASS "RESIDUE\t" . (join "\t", @LIBRARY_ORDER) . "\n";

foreach my $res ( @OFFICIAL_AAS ) { 

  print SC_SIG "$res";
  print SC_ASS "$res";

  foreach my $lib ( @LIBRARY_ORDER ) { 

    my $file = "$LIBRARIES{$lib}{'dirs'}{$res}/$FILE_IN_STATE_COUNTS";
    if ( not -e $file ) { 
      die "Couldn't find file, '$file'";
    }
  
    my $sig = Get_Significant   ( $file );
    my $ass = Get_Assignability ( $file );
    
    printf SC_ASS "\t%.4f", $ass;
    printf SC_SIG "\t%.4f", $sig;
  }

  print SC_ASS "\n";
  print SC_SIG "\n";
}

close SC_SIG or die "Couldn't close dyn sig, '$FILE_OUT_SIG'. $!";
close SC_ASS or die "Couldn't close dyn ass, '$FILE_OUT_ASS'. $!";





####################################################################
#                                                         WAIT TIMES
####################################################################

open WAITS, ">$FILE_OUT_WAIT_COUNTS" or die "Couldn't open wait counts, '$FILE_OUT_WAIT_COUNTS'. $!";

print WAITS "RESIDUE\tLIB\tANGLE\tWAIT\n";

foreach my $res ( @OFFICIAL_AAS ) { 
  my @angles = @{$AA{$res}{'rotamer_names'}};

  foreach my $lib ( @LIBRARY_ORDER ) { 

    my $file = "$LIBRARIES{$lib}{'dirs'}{$res}/$FILE_IN_WAIT_COUNTS";
    my $line = `grep mean_counts $file`;
    my ($mc, $res, @data) = split /\s+/, $line;

    foreach my $i ( 0..$#angles ) { 
      my $chi = $angles[$i];
      my $wait = $data[$i];
      printf WAITS "%s\t%s\t%s\t%.4f\n", $res, $lib, $chi, $wait;
    }
  }
}
close WAITS or die "Couldn't close wait counts, '$FILE_OUT_WAIT_COUNTS'. $!";



####################################################################
#                                             RESIDUE LEVEL DYNAMICS
####################################################################

if ( not -e "$SUMM/op" ) { 
  mkdir "$SUMM/op" or die "Couldn't mkdir $SUMM/op. $!";
}
foreach my $lib ( @LIBRARY_ORDER ) { 
  foreach my $res ( @ALL_AAS ) { 
    my $file_source = "$LIBRARIES{$lib}{'dirs'}{$res}/dynamics.residues.dat";
    my $file_target = "$SUMM/op/dynamics.$res.$lib.dat";
    `cp $file_source $file_target`;
  }
}

####################################################################
#                                                                   
####################################################################


####################################################################
#                                                                   
####################################################################


####################################################################
#                                                                   
####################################################################


