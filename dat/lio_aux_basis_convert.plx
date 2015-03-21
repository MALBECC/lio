#!/usr/bin/perl -w

%ELEMENT_MAP = (
H => 1 ,
He => 2 ,
Li => 3 ,
Be => 4 ,
B => 5 ,
C => 6 ,
N => 7 ,
O => 8 ,
F => 9 ,
Ne => 10 ,
Na => 11 ,
Mg => 12 ,
Al => 13 ,
Si => 14 ,
P => 15 ,
S => 16 ,
Cl => 17 ,
Ar => 18 ,
K => 19 ,
Ca => 20 ,
Sc => 21 ,
Ti => 22 ,
V => 23 ,
Cr => 24 ,
Mn => 25 ,
Fe => 26 ,
Co => 27 ,
Ni => 28 ,
Cu => 29 ,
Zn => 30 ,
Ga => 31 ,
Ge => 32 ,
As => 33 ,
Se => 34 ,
Br => 35 ,
Kr => 36 ,
Rb => 37 ,
Sr => 38 ,
Y => 39 ,
Zr => 40 ,
Nb => 41 ,
Mo => 42 ,
Tc => 43 ,
Ru => 44 ,
Rh => 45 ,
Pd => 46 ,
Ag => 47 ,
Cd => 48 ,
In => 49 ,
Sn => 50 ,
Sb => 51 ,
Te => 52 ,
I => 53 ,
Xe => 54 ,
Cs => 55 ,
Ba => 56 ,
La => 57 ,
Ce => 58 ,
Pr => 59 ,
Nd => 60 ,
Pm => 61 ,
Sm => 62 ,
Eu => 63 ,
Gd => 64 ,
Tb => 65 ,
Dy => 66 ,
Ho => 67 ,
Er => 68 ,
Tm => 69 ,
Yb => 70 ,
Lu => 71 ,
Hf => 72 ,
Ta => 73 ,
W => 74 ,
Re => 75 ,
Os => 76 ,
Ir => 77 ,
Pt => 78 ,
Au => 79 ,
Hg => 80 ,
Tl => 81 ,
Pb => 82 ,
Bi => 83 ,
Po => 84 ,
At => 85 ,
Rn => 86 ,
Fr => 87 ,
Ra => 88 ,
Ac => 89 ,
Th => 90 ,
Pa => 91 ,
U => 92 ,
Np => 93 ,
Pu => 94 ,
Am => 95 ,
Cm => 96 ,
Bk => 97 ,
Cf => 98 ,
Es => 99 ,
Fm => 100 ,
Md => 101 ,
No => 102 ,
Lr => 103 ,
Rf => 104 ,
Db => 105 ,
Sg => 106 ,
Bh => 107 ,
Hs => 108 ,
Mt => 109 ,
Ds => 110 ,
Rg => 111 ,
Cn => 112 ,
Uut => 113 ,
Fl => 114 ,
Uup => 115 ,
Lv => 116 ,
Uus => 117 ,
Uuo => 118
);

$line = <>;
print "# $line";

$curr_el = -1;
$curr_el_name = "";
@a_vals = ();
@c_vals = ();
@p_a_vals = ();
@p_c_vals = ();
$s_funcs = 0;
$p_funcs = 0;
$d_funcs = 0;
$comment = "";

while (<>) {
    if (/^\s*#/ or /^\s*$/) {
        if (@a_vals > 0) {
            print $comment;
        }
        $comment = $_;
        next;
    }
    if (not /[A-Z]/) { die "Error 1 while parsing - incorrect format\n"; }

    if (/END/) { print "endbasis\n"; last; }

    chomp;
    s/^\s+//;
    s/\s+$//;

    my @line = split;
    if (@line !=2 ) { die "Error 2 while parsing - incorrect format\n"; }

    my $el = $line[0];
    my $el_num = $ELEMENT_MAP{$el};
    my $type = $line[1];

    if ($el_num != $curr_el and $curr_el != -1) {
        if (@p_a_vals != 0) {
            push @a_vals, @p_a_vals;
            push @c_vals, @p_c_vals;
            @p_a_vals = (); @p_c_vals = ();
        }
        print "# $curr_el_name\n";
        print "gaussian\n";
        my $total_funcs = $s_funcs + $p_funcs + $d_funcs;
        printf " %d %d %d\n",$curr_el,$total_funcs,$total_funcs;
        print " ";
        print "1 " for (1..$total_funcs);
        print "\n ";
        print "0 " for (1..$s_funcs);
        print "1 " for (1..$p_funcs);
        print "2 " for (1..$d_funcs);
        print "\n";
        for (1..$total_funcs) {
            printf "%24.10E%24.10E\n",$a_vals[$_-1],$c_vals[$_-1];
        }
        $s_funcs = 0; $p_funcs = 0; $d_funcs = 0;
        @a_vals = (); @c_vals = ();
    }
    $curr_el = $el_num;
    $curr_el_name = $el;

    my $next_line = <>;
    chomp $next_line;
    $next_line =~ s/^\s+//;
    $next_line =~ s/\s+$//;
    my @next_line = split /\s+/,$next_line;

    if (@next_line != 2 and @next_line != 3) { die "Error 3 while parsing - incorrect format\n"; }

    if ($type eq 'S') {
        if (@next_line != 2) { die "Error while parsing - incorrect format\n"; }
        $s_funcs++;
        push @a_vals, $next_line[0];
        push @c_vals, $next_line[1];
    }
    elsif ($type eq 'SP') {
        if (@next_line != 3) { die "Error while parsing - incorrect format\n"; }
        $s_funcs++; $p_funcs++;
        push @a_vals, $next_line[0];
        push @c_vals, $next_line[1];
        push @p_a_vals, $next_line[0];
        push @p_c_vals, $next_line[2];
    }
    elsif ($type eq 'P') {
        if (@next_line != 2) { die "Error while parsing - incorrect format\n"; }
        $p_funcs++;
        if (@p_a_vals != 0) {
            push @a_vals, @p_a_vals;
            push @c_vals, @p_c_vals;
            @p_a_vals = (); @p_c_vals = ();
        }
        push @a_vals, $next_line[0];
        push @c_vals, $next_line[1];
    }
    elsif ($type eq 'D') {
        if (@next_line != 2) { die "Error while parsing - incorrect format\n"; }
        $d_funcs++;
        if (@p_a_vals != 0) {
            push @a_vals, @p_a_vals;
            push @c_vals, @p_c_vals;
            @p_a_vals = (); @p_c_vals = ();
        }
        push @a_vals, $next_line[0];
        push @c_vals, $next_line[1];
    }
    else { die "Error while parsing - unexpected function type\n"; }
}
