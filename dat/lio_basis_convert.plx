#!/usr/bin/perl -w

# This script converts from a deMon2K format basis set file to
# the format expected by lio
#
# Usage: lio_basis_convert.plx <deMon2K file>
# The lio basis set file is written to stdout

%ELEMENT_MAP = (
    HYDROGEN => 1 ,
    HELIUM => 2 ,
    LITHIUM => 3 ,
    BERYLLIUM => 4 ,
    BORON => 5 ,
    CARBON => 6 ,
    NITROGEN => 7 ,
    OXYGEN => 8 ,
    FLUORINE => 9 ,
    NEON => 10 ,
    SODIUM => 11 ,
    MAGNESIUM => 12 ,
    ALUMINUM => 13 ,
    SILICON => 14 ,
    PHOSPHORUS => 15 ,
    PHOSPHOROUS => 15 ,
    SULFUR => 16 ,
    CHLORINE => 17 ,
    ARGON => 18 ,
    POTASSIUM => 19 ,
    CALCIUM => 20 ,
    SCANDIUM => 21 ,
    TITANIUM => 22 ,
    VANADIUM => 23 ,
    CHROMIUM => 24 ,
    MANGANESE => 25 ,
    IRON => 26 ,
    COBALT => 27 ,
    NICKEL => 28 ,
    COPPER => 29 ,
    ZINC => 30 ,
    GALLIUM => 31 ,
    GERMANIUM => 32 ,
    ARSENIC => 33 ,
    SELENIUM => 34 ,
    BROMINE => 35 ,
    KRYPTON => 36 ,
    RUBIDIUM => 37 ,
    STRONTIUM => 38 ,
    YTTRIUM => 39 ,
    ZIRCONIUM => 40 ,
    NIOBIUM => 41 ,
    MOLYBDENUM => 42 ,
    TECHNETIUM => 43 ,
    RUTHENIUM => 44 ,
    RHODIUM => 45 ,
    PALLADIUM => 46 ,
    SILVER => 47 ,
    CADMIUM => 48 ,
    INDIUM => 49 ,
    TIN => 50 ,
    ANTIMONY => 51 ,
    TELLURIUM => 52 ,
    IODINE => 53 ,
    XENON => 54 ,
    CESIUM => 55 ,
    BARIUM => 56 ,
    LANTHANUM => 57 ,
    CERIUM => 58 ,
    PRASEODYMIUM => 59 ,
    NEODYMIUM => 60 ,
    PROMETHIUM => 61 ,
    SAMARIUM => 62 ,
    EUROPIUM => 63 ,
    GADOLINIUM => 64 ,
    TERBIUM => 65 ,
    DYSPROSIUM => 66 ,
    HOLMIUM => 67 ,
    ERBIUM => 68 ,
    THULIUM => 69 ,
    YTTERBIUM => 70 ,
    LUTETIUM => 71 ,
    HAFNIUM => 72 ,
    TANTALUM => 73 ,
    TUNGSTEN => 74 ,
    RHENIUM => 75 ,
    OSMIUM => 76 ,
    IRIDIUM => 77 ,
    PLATINUM => 78 ,
    GOLD => 79 ,
    MERCURY => 80 ,
    THALLIUM => 81 ,
    LEAD => 82 ,
    BISMUTH => 83 ,
    POLONIUM => 84 ,
    ASTATINE => 85 ,
    RADON => 86 ,
    FRANCIUM => 87 ,
    RADIUM => 88 ,
    ACTINIUM => 89 ,
    THORIUM => 90 ,
    PROTACTINIUM => 91 ,
    URANIUM => 92 ,
    NEPTUNIUM => 93 ,
    PLUTONIUM => 94 ,
    AMERICIUM => 95 ,
    CURIUM => 96 ,
    BERKELIUM => 97 ,
    CALIFORNIUM => 98 ,
    EINSTEINIUM => 99 ,
    FERMIUM => 100 ,
    MENDELEVIUM => 101 ,
    NOBELIUM => 102 ,
    LAWRENCIUM => 103 ,
    RUTHERFORDIUM => 104 ,
    DUBNIUM => 105 ,
    SEABORGIUM => 106 ,
    BOHRIUM => 107 ,
    HASSIUM => 108 ,
    MEITNERIUM => 109 ,
    DARMSTADTIUM => 110 ,
    ROENTGENIUM => 111 ,
    COPERNICIUM => 112 ,
    UNUNTRIUM => 113 ,
    FLEROVIUM => 114 ,
    UNUNPENTIUM => 115 ,
    LIVERMORIUM => 116 ,
    UNUNSEPTIUM => 117 ,
    UNUNOCTIUM => 118
);

$cartesian = -1;
while(<>) {
    if (/^\s*#/ or /^\s*$/) {
        if (/This basis set uses/) {
            $cartesian = 1 if /cartesian/;
            $cartesian = 0 if /spherical/;
        }
        print;
        while(<>) {
            if (/This basis set uses/) {
                $cartesian = 1 if /cartesian/;
                $cartesian = 0 if /spherical/;
            }
            last if !(/^\s*#/ or /^\s*$/);
            print;
        }
    }

    die "There was an error reading the basis file - it was not the expected format\n" if not /O-([A-Z]+)/;

    chomp;
    s/O-([A-Z]+)\s+.*/$1/;
    my $atomic_number = $ELEMENT_MAP{$_};
    print "# $_\n";
    while(<>) {
        last if !(/^\s*#/ or /^\s*$/);
        print;
    }
    my @line = split;
    my $num_funcs = $line[0];

    my @a_vals = ();
    my @c_vals = ();
    my @prims = ();
    my @types = ();
    my $total_prims = 0;

    print "gaussian\n";
    for (1..$num_funcs) {
        while (<>) {
            last if !(/^\s*#/ or /^\s*$/);
            print;
        }
        @line = split;
        push @types, $line[1];
        push @prims, $line[2];
        my $prim = $line[2];
        $total_prims += $prim;
        for (1..$prim) {
            while (<>) {
                last if !(/^\s*#/ or /^\s*$/);
                print;
            }
            s/D/E/g;
            s/d/E/g;
            @line = split;
            push @a_vals,$line[0];
            push @c_vals,$line[1];
        }
    }

    print " $atomic_number $total_prims $num_funcs\n";
    print " ";
    for (@prims) { print "$_ "; }
    print "\n ";
    for (@types) { print "$_ "; }
    print "\n";
    for (1..$total_prims) {
        printf "%24.10E%24.10E\n",$a_vals[$_-1],$c_vals[$_-1];
    }
}
print "endbasis\n";

#if    ($cartesian == 1) { print "BASIS CONTAINS CARTESIAN COMPONENTS\n"; }
#elsif ($cartesian == 0) { print "BASIS CONTAINS SPHERICAL COMPONENTS\n"; }
#else                    { die "COULD NOT DETERMINE BASIS TYPE\n"; }
