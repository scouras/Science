#!/usr/bin/perl

BEGIN { push @INC, "/users/scouras/code/RotLib",
                   "/users/scouras/lab/code/RotLib" }

use warnings;
use strict;
use Data::Dumper;
use Rotamer;
use Utility;


my ( $ROOT, $SUMM, @LIBRARY_ORDER_STATIC, @LIBRARY_ORDER_DYNAMIC, %LIBRARIES_STATIC, %LIBRARIES_DYNAMIC );

my $CONFIG_FILE = $ARGV[0] || "/users/scouras/code/RotLib/native_810.static.dat";

if ( not -e $CONFIG_FILE ) { 
  die "Config file '$CONFIG_FILE' not found.";
}

my $CONFIG_CONTENTS = `cat $CONFIG_FILE`;
eval $CONFIG_CONTENTS;
if ( not -e $SUMM ) { mkdir $SUMM or die "Couldn't mkdir $SUMM. $!" }

#print Dumper ( \%LIBRARIES_STATIC );

my %AA = %Utility::AA;

my %CONFIG = (
  rotamer_mode        => "static",
  home_dir            => "/users/scouras",
  #root_dir            => "/users/scouras/lab/rotamers",
  root_dir            => "/users/scouras/rotamers",
  overwrite_rotamers  => 0,
  overwrite_dihedrals => 0,
);


my $FILE_COUNTS           = "$SUMM/residue_counts.dat";
my $FILE_LIBRARY          = "$SUMM/library.dat";
my $FILE_MINOR            = "$SUMM/minors.dat";
my $FILE_DISPLACE_DETAIL  = "$SUMM/displacement.majors.dat";
my $FILE_DISPLACE_SUMMARY = "$SUMM/displacement.summary.dat";
my $FILE_DISPLACE_MINORS  = "$SUMM/displacement.minors.dat";
my $FILE_SIG              = "$SUMM/significant.dat";



#======================================== Load Residue Libraries


my $ROT = RotLib::Rotamer -> Initialize ( \%CONFIG, \%AA );

#my @OFFICIAL_AAS = qw(LEU);
my @OFFICIAL_AAS = qw(ARG ASN ASP CYH CYS GLN GLU HID HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL);

my %RESIDUE_CONVERSION = (
  HID => 'HIS',
  HIE => 'HIS',
  CYH => 'CYS',
  CYS => 'CYH',
);


#================================== Load Residue Library Directories

my @LIB_PAIRS = ();
foreach my $i ( 0..$#LIBRARY_ORDER_STATIC ) { 
  foreach my $j ( $i+1..$#LIBRARY_ORDER_STATIC ) { 
    push @LIB_PAIRS, [$LIBRARY_ORDER_STATIC[$i], $LIBRARY_ORDER_STATIC[$j]];
  }
}


while ( my ( $LIB_NAME, $LIB ) = each %LIBRARIES_STATIC ) {
  foreach my $res_orig ( @OFFICIAL_AAS ) { 
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


#die Dumper (\%LIBRARIES_STATIC);



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


####################################################################
#                                           SUMMARIZE RESIDUE COUNTS
####################################################################

my %RESIDUE_COUNTS = ();

while ( my ( $LIB_NAME, $LIB ) = each %LIBRARIES_STATIC ) {
  foreach my $res ( @OFFICIAL_AAS ) { 

    my $res_dir = $LIB->{'dirs'}{$res};
    my $file_state_counts = "$res_dir/state_counts.dat";
    if ( not -e $file_state_counts ) { 
      die "Couldn't find state counts file, $LIB_NAME  $res $res_dir '$file_state_counts'.";
    }
    my $line = `grep 'Total residues' $file_state_counts`;
    $line =~ /(\d+)\s*$/;
    my $count = $1;

    $LIB->{'count'}{$res} = $count;
  }
}

open COUNT, ">$FILE_COUNTS" or die "Couldn't open counts, '$FILE_COUNTS'. $!";
print COUNT (join "\t", "RES", @LIBRARY_ORDER_STATIC) . "\n";
foreach my $res ( @OFFICIAL_AAS ) { 
  print COUNT "$res\t" . (join "\t", map { $LIBRARIES_STATIC{$_}{'count'}{$res} } @LIBRARY_ORDER_STATIC) . "\n";
}
close COUNT or die "Couldn't close counts, '$FILE_COUNTS'. $!";



####################################################################
#                                                  READ ROTAMER DATA
####################################################################

my @SIG_CUTOFFS = (0.0, 0.001, 0.01, 0.05);
my %SIG = ();

while ( my ( $LIB_NAME, $LIB ) = each %LIBRARIES_STATIC ) {
  foreach my $res ( @OFFICIAL_AAS ) { 

    my $total_rotamers = 0;
    my @angle_names = @{$AA{$res}{'rotamer_names'}};

    my $res_dir = $LIB->{'dirs'}{$res};
    my $file_library = "$res_dir/library.new.dat";
    if ( not -e $file_library ) { 
      die "Couldn't find library, '$file_library'. $!";
    }

    my @data = map { [ split /\s+/, $_ ] } split "\n", `cat $file_library`;
    foreach my $datum ( @data[1..$#data] ) {
      my ($res_virt, $major, $min1, $min2, $min3, $min4, $p) = @{$datum}[0..6];

      #==== Record Major Populations
      $LIB->{'library'}{$res}[$major-1] = {
        major   => $major,
        minor1  => $min1,
        minor2  => $min2,
        minor3  => $min3,
        minor4  => $min4,
        p       => $p,
      };

      foreach my $i ( 0..$#SIG_CUTOFFS ) { 
        if ( $p > $SIG_CUTOFFS[$i] ) { 
          $SIG{$LIB_NAME}{$res}[$i]++;
        }
      }

      #==== Record Minor Populations
      my @minors = ($min1, $min2, $min3, $min4);
      foreach my $i (0..$#angle_names) { 
        $LIB->{'minors'}{$res}[$i][$minors[$i]-1] += $p;
      }


      $total_rotamers++;
      #print "$LIB_NAME\t$res\t$total_rotamers\t$AA{$res}{'total_rotamers'}\n";
    }
    if ( $total_rotamers != $AA{$res}{'total_rotamers'} ) { 
      #die "Wrong rotamer count in $LIB_NAME $res ($total_rotamers vs. $AA{$res}{'total_rotamers'}";
    }

  }
}


#===== Output Major Rotamer Populations

open LIB, ">$FILE_LIBRARY" or die "Couldn't open library, '$FILE_LIBRARY'. $!";
print LIB "RES\tMAJOR\tM1\tM2\tM3\tM4\t" . (join "\t", @LIBRARY_ORDER_STATIC) . "\n";
foreach my $res ( @OFFICIAL_AAS ) { 
  my $total_rotamers = $AA{$res}{'total_rotamers'};
  foreach my $major ( 1..$total_rotamers ) {
    my @names = @{$MINOR_NAMES{$res}[$major-1]};
    print LIB "$res\t$major";
    foreach my $i (0..3) { print LIB "\t" . (defined $names[$i] ? $names[$i] : '') }
    foreach my $l (@LIBRARY_ORDER_STATIC) { 
      my $p = $LIBRARIES_STATIC{$l}{'library'}{$res}[$major-1]{'p'};
      #print STDERR "$res\t$major\t$l\t$p\n";
      print LIB "\t$p";
    }
    print LIB "\n";
  }
}

#===== Output Minor Rotamer Populations

open MINOR, ">$FILE_MINOR" or die "Couldn't open minor, '$FILE_MINOR'. $!";
print MINOR "RES\tANGLE\tBIN\t" . (join "\t", @LIBRARY_ORDER_STATIC) . "\n";
foreach my $res ( @OFFICIAL_AAS ) { 

  my @angle_names = @{$AA{$res}{'rotamer_names'}};
  foreach my $a ( 0..$#angle_names ) { 
    my $angle = $angle_names[$a];
    my @minor_names = @{$MINOR_LOOKUP{$res}[$a]};;
    foreach my $m ( 0..$#minor_names ) { 
      my $minor = $minor_names[$m];
      print MINOR "$res\t$angle\t$minor";
      foreach my $lib ( @LIBRARY_ORDER_STATIC ) {
        printf MINOR "\t%.8f", $LIBRARIES_STATIC{$lib}{'minors'}{$res}[$a][$m];
      }
      print MINOR "\n";
    }
  }
}
close MINOR or die "Couldn't close minor, '$FILE_MINOR'. $!";



#===== Output significant rotamer counts

open SIG, ">$FILE_SIG" or die "Couldn't open sig file, '$FILE_SIG. $!";
print SIG "RES\tCUTOFF\t" . (join "\t", @LIBRARY_ORDER_STATIC) . "\n";
foreach my $i ( 0..$#SIG_CUTOFFS ) { 
  foreach my $res ( @OFFICIAL_AAS ) { 
    printf SIG "%s\t%.1f%%", $res, $SIG_CUTOFFS[$i]*100;
    foreach my $lib ( @LIBRARY_ORDER_STATIC ) { 
      printf SIG "\t%u", $SIG{$lib}{$res}[$i] || 0;
    }
    print SIG "\n";
  }
}
close SIG or die "Couldn't close sig file, '$FILE_SIG'. $!";




####################################################################
#                          POPULATION DISPLACEMENT BETWEEN LIBRARIES_STATIC
####################################################################

my %DISPLACEMENTS = ();
my %DISPLACEMENT_SUM = ();

open DISP, ">$FILE_DISPLACE_DETAIL"  or die "Coudln't open '$FILE_DISPLACE_DETAIL'. $!";
open SUMM, ">$FILE_DISPLACE_SUMMARY" or die "Coudln't open '$FILE_DISPLACE_SUMMARY'. $!";

print DISP "RES\tMAJOR\tM1\tM2\tM3\tM4\t" . (join "\t", map { "$_->[0]-$_->[1]" } @LIB_PAIRS) . "\n";
print SUMM "RES\t"                        . (join "\t", map { "$_->[0]-$_->[1]" } @LIB_PAIRS) . "\n";

foreach my $res ( @OFFICIAL_AAS ) { 
  my $total_rotamers = $AA{$res}{'total_rotamers'};
  print SUMM $res;

  foreach my $major ( 1..$total_rotamers ) {
    print DISP "$res\t$major\t" . (join "\t", map { $MINOR_NAMES{$res}[$major-1][$_] || '' } (0..3));
    #print STDERR "$res\t$major\t" . (join "\t", map { $MINOR_NAMES{$res}[$major-1][$_] || '' } (0..3)) . "\n";
    foreach my $pair ( @LIB_PAIRS ) { 
      my ($LIB1_NAME, $LIB2_NAME) = @$pair;

      #print STDERR "\t$LIB1_NAME\t$LIB2_NAME\n";

      my $d = $LIBRARIES_STATIC{$LIB2_NAME}{'library'}{$res}[$major-1]{'p'}
            - $LIBRARIES_STATIC{$LIB1_NAME}{'library'}{$res}[$major-1]{'p'};

      $DISPLACEMENTS{$LIB1_NAME}{$LIB2_NAME}{$res}[$major-1] = $d;
      $DISPLACEMENT_SUM{$LIB1_NAME}{$LIB2_NAME}{$res} += abs($d/2);

      printf DISP "\t%.8f", $d;

    }
    print DISP "\n";
  }

  foreach my $pair ( @LIB_PAIRS ) { 
    my ($LIB1_NAME, $LIB2_NAME) = @$pair;
    printf SUMM "\t%.8f", $DISPLACEMENT_SUM{$LIB1_NAME}{$LIB2_NAME}{$res};
  }
  print SUMM "\n";
}

close DISP or die "$!";
close SUMM or die "$!";


#exit;


####################################################################
#                             MINOR ROTAMER POPULATION DISPLACEMENTS
####################################################################

open MINOR, ">$FILE_DISPLACE_MINORS" or die "Couldn't open minor displacement, $FILE_DISPLACE_MINORS. $!";

print MINOR "RES\tANGLE\tBIN\t" . (join "\t", map { "$_->[0]-$_->[1]" } @LIB_PAIRS) . "\n";

foreach my $res ( @OFFICIAL_AAS ) { 
  
  my @angle_names = @{$AA{$res}{'rotamer_names'}};
  
  foreach my $a ( 0..$#angle_names ) { 
    my $angle = $angle_names[$a];
    my @minor_names = @{$MINOR_LOOKUP{$res}[$a]};;
  
    foreach my $m ( 0..$#minor_names ) { 
      my $minor = $minor_names[$m];
      print MINOR "$res\t$angle\t$minor";

      foreach my $pair ( @LIB_PAIRS ) { 
        my ($LIB1_NAME, $LIB2_NAME) = @$pair;

        my $d = $LIBRARIES_STATIC{$LIB1_NAME}{'minors'}{$res}[$a][$m]
              - $LIBRARIES_STATIC{$LIB2_NAME}{'minors'}{$res}[$a][$m];

        printf MINOR "\t%.8f", $d;
      }
      print MINOR "\n";
    }

  }
}















####################################################################
####################################################################
####################################################################
####################################################################
####################################################################





#============================================================= Files

my $FILE_IN_STATE_COUNTS  = "state_counts.dat";
my $FILE_OUT_DYN_ASS          = "$SUMM/dynamic.assignable.dat";
my $FILE_OUT_DYN_SIG          = "$SUMM/dynamic.significant.dat";
my $FILE_OUT_STA_ASS          = "$SUMM/static.assignable.dat";
my $FILE_OUT_STA_SIG          = "$SUMM/static.significant.dat";

my $FILE_IN_WAIT_COUNTS   = "dynamics_values.mean_counts.dat";
my $FILE_OUT_WAIT_COUNTS  = "$SUMM/dyn.wait_counts.dat";



#================================== Load Residue Library Directories
my @ALL_AAS = ('ALA', 'GLY', @OFFICIAL_AAS);
while ( my ( $LIB_NAME, $LIB ) = each %LIBRARIES_DYNAMIC ) {
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






####################################################################
#                                  ASSIGNABLE AND SIGNIFICANT COUNTS
####################################################################

sub Get_Assignability { 
  my $file = $_[0];
  my $THRESHOLD = "0.5000";
  my $line = `grep "^$THRESHOLD" $file`;
  print "LINE $file\n$line\n";
  my ($t, $d0, $d1, $dx) = split /\s+/, $line;
  return $d1;
}

sub Get_Significant { 
  my $file = $_[0];
  my $THRESHOLD = "0.1000";
  my $line = `grep "^$THRESHOLD" $file`;
  print "LINE $file\n$line\n";
  my ( $t, @data ) = split /\s+/, $line;
  my $sum = 0;
  map { $sum += ( $_ * $data[$_] ) } 0..$#data;
  return $sum;
}

open SC_SIG, ">$FILE_OUT_DYN_SIG" or die "Couldn't open dyn sig, '$FILE_OUT_DYN_SIG'. $!";
open SC_ASS, ">$FILE_OUT_DYN_ASS" or die "Couldn't open dyn ass, '$FILE_OUT_DYN_ASS'. $!";

print SC_SIG "RESIDUE\t" . (join "\t", @LIBRARY_ORDER_DYNAMIC) . "\n";
print SC_ASS "RESIDUE\t" . (join "\t", @LIBRARY_ORDER_DYNAMIC) . "\n";

foreach my $res ( @OFFICIAL_AAS ) { 

  print SC_SIG "$res";
  print SC_ASS "$res";

  foreach my $lib ( @LIBRARY_ORDER_DYNAMIC ) { 

    my $file = "$LIBRARIES_DYNAMIC{$lib}{'dirs'}{$res}/$FILE_IN_STATE_COUNTS";
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

close SC_SIG or die "Couldn't close dyn sig, '$FILE_OUT_DYN_SIG'. $!";
close SC_ASS or die "Couldn't close dyn ass, '$FILE_OUT_DYN_ASS'. $!";




open SC_SIG, ">$FILE_OUT_STA_SIG" or die "Couldn't open dyn sig, '$FILE_OUT_STA_SIG'. $!";
open SC_ASS, ">$FILE_OUT_STA_ASS" or die "Couldn't open dyn ass, '$FILE_OUT_STA_ASS'. $!";

print SC_SIG "RESIDUE\t" . (join "\t", @LIBRARY_ORDER_DYNAMIC) . "\n";
print SC_ASS "RESIDUE\t" . (join "\t", @LIBRARY_ORDER_DYNAMIC) . "\n";

foreach my $res ( @OFFICIAL_AAS ) { 

  print SC_SIG "$res";
  print SC_ASS "$res";

  foreach my $lib ( @LIBRARY_ORDER_DYNAMIC ) { 

    my $file = "$LIBRARIES_STATIC{$lib}{'dirs'}{$res}/$FILE_IN_STATE_COUNTS";
    if ( not -e $file ) { 
      die "Couldn't find file, '$file'";
    }
  
    my $ass = Get_Assignability ( $file );
    my $sig = Get_Significant   ( $file );
    
    printf SC_ASS "\t%.4f", $ass;
    printf SC_SIG "\t%.4f", $sig;
  }

  print SC_ASS "\n";
  print SC_SIG "\n";
}

close SC_SIG or die "Couldn't close dyn sig, '$FILE_OUT_STA_SIG'. $!";
close SC_ASS or die "Couldn't close dyn ass, '$FILE_OUT_STA_ASS'. $!";





####################################################################
#                                                         WAIT TIMES
####################################################################

open WAITS, ">$FILE_OUT_WAIT_COUNTS" or die "Couldn't open wait counts, '$FILE_OUT_WAIT_COUNTS'. $!";

print WAITS "RESIDUE\tLIB\tANGLE\tWAIT\n";

foreach my $res ( @OFFICIAL_AAS ) { 
  my @angles = @{$AA{$res}{'rotamer_names'}};

  foreach my $lib ( @LIBRARY_ORDER_DYNAMIC ) { 

    my $file = "$LIBRARIES_DYNAMIC{$lib}{'dirs'}{$res}/$FILE_IN_WAIT_COUNTS";
    my $line = `grep mean_counts $file`;
    my ($mc, $res, $major, @data) = split /\s+/, $line;

    printf WAITS "%s\t%s\t%s\t%.4f\n", $res, $lib, "Major", $major;
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
foreach my $lib ( @LIBRARY_ORDER_DYNAMIC ) { 
  foreach my $res ( @ALL_AAS ) { 
    my $file_source = "$LIBRARIES_DYNAMIC{$lib}{'dirs'}{$res}/dynamics.residues.dat";
    my $file_target = "$SUMM/op/dynamics.$res.$lib.dat";
    `cp $file_source $file_target`;
  }
}

















