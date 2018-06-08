#!/bin/env perl

# Takes as args a csv file and a glycoct file.
# Print on stdout the corresponding inp format for POLYS tool.

# use
# perl convertor1.0.pl phi_psi_all_disacc.csv HP_dp02_0001.glycoct > HP_dp02_0001.inp

use strict;
use warnings;

(@ARGV == 2) || die "E: This script needs 2 arg: a csv file with phi-psi value and the glycoct file.\n";
my ($csvFile, $glycoctFile) = @ARGV;


# @pile will contain information of our glycoct file
# each element of @pile is a hash (or residue in our glycoct file)
# all hashs have keys: 'type' (type of residue: basetype or substituent) and 'nom' (name of residue)
# somes hashs have keys: 'elem_suivant' (index of next residue in @pile which is a basetype),
#                        'liaison_elem_suivant' (bound with next residue which is a basetype),
#                        'liste_feuilles' (table of index of next residues in @pile which are substituents),
#                        'liste_liaisons_feuilles' (table of bound with next residues which are substituents),
#                        'nom_complet_polys' (polys name of basetype residue)
my @pile;

# %glycoct2polys holds for each residue in glycoct format, corresponding residue in polys format
# key is glycoct format, value is polys format
my %glycoct2polys = ("a-dglc-HEX-1:5" => "aDGlcp", 
		     "a-lido-HEX-1:5|6:a" => "aLIdopA",
		     "b-dgal-HEX-1:5" => "bDGalp",
		     "b-dglc-HEX-1:5" => "bDGlcp",
		     "b-dglc-HEX-1:5|6:a" => "bDGlcpA",
		     "amino" => "N",
		     "n-acetyl" => "NAc",
		     "n-sulfate" => "NS",
		     "sulfate" => "S");

########################################################
# 1st part:read csv file and glycoct file and store data
########################################################
# %phiPsiAngles will contain information of our csv file
# key is disaccharide in polys format, value is phi-psi angle value
my %phiPsiAngles;
&fillAnglesHash($csvFile, \%phiPsiAngles);

open(FILE_IN, $glycoctFile) || die "E: Cannot open glycoct file $glycoctFile\n";
# $section: "RES" or "LIN", says which section we're in
my $section;
while(my $line = <FILE_IN>){
    chomp($line);
    if($line eq "RES"){
	    ($section) && die "E: Found RES header but section already defined as $section , shouldn't happen!\n";
	    $section = "RES";
    }
    elsif($line eq "LIN"){
	    ($section eq "RES") || die "E: found LIN header but we were not in RES section? section was $section\n";
	    $section = "LIN";
    }
    elsif($line =~ m/^\d+[bs]:.+/){ # RES bloc line with basetype or substituent (eg: 1b:a-dglc-HEX-1:5)
	    ($section eq "RES") || die "E: found RES line but not in RES section:\n$line\n";
        &treatmentResiduLine(\@pile, $line);
    }
    elsif($line =~ m/^\d+:\d+[odhn]\(\d\+\d\)\d+[odhn]$/){ # LIN bloc line (eg: 1:1d(2+1)2n)
	    ($section eq "LIN") || die "E: found LIN line but not in LIN section:\n$line\n";
        &treatmentLinkageLine(\@pile, $line);
    }
    elsif($line =~ m/^\d+[ra]:.+/){ # RES bloc line with repeating unit or alternative unit (eg: 1r:r1)
        die "E: This program can't process repeating or alternative unit\n";
    }
    else{
	    die "E parsing glycoctFile $glycoctFile: unexpected line:\n$line\n";
    }
}
close FILE_IN;
################################################################################
# 2nd part:for each basetype node, creation of full POLYS name with substituents
################################################################################
foreach my $i (0..$#pile){
    if($pile[$i]{"type"} eq "b"){
        # case basetype
        &createFullName(\@pile, $i);
    }
}
###########################
# 3rd part:write polys file
###########################
print("PRIMARY\n");
while(scalar(@pile)){
    my $last = pop(@pile);
    if(${$last}{"type"} eq "b"){
        # case basetype
        my $printLine = " <".${$last}{"nom_complet_polys"}.">"; # eg: <bDGlcpA>
        if(defined ${$last}{"elem_suivant"}){
            # case our basetype has a next basetype
            my $connexion = ${$last}{"liaison_elem_suivant"};
            ($connexion =~ m/^(\d+)\+(\d+)$/) || die "E in connexion of ".${$last}{"nom_complet_polys"}.":\n$connexion\n";
            my $parentAttachmentPosition = $1;
            my $childAttachmentPosition = $2;
            my $indexNextElem = ${$last}{"elem_suivant"};
            $indexNextElem--; # because array starts at 0
            my $key = ${$last}{"nom_complet_polys"}." ".$childAttachmentPosition."-".$parentAttachmentPosition." ".$pile[$indexNextElem]{"nom_complet_polys"}; # eg: aLIdopA_2S_1C4 1-4 aDGlcpN_6S
            ($phiPsiAngles{$key}) || die "E: This key '".$key."' doesn't exist in our phi-psi value hash\n";
            $printLine = $printLine." ( ".$childAttachmentPosition.":".$parentAttachmentPosition."; ".$phiPsiAngles{$key}.")"; # eg: <bDGlcpA> ( 1:4; 40.00; -110.00)
        }
        print($printLine."\n");
    }
}
print("\nSTOP\n");

###########
# functions
###########

# take as args name of csv file with phi-psi values
# and a reference of a hash
# extract information in hash 
sub fillAnglesHash{
    (@_ == 2) || die "E: fillAnglesHash requires 2 arg\n";
    my ($csvFile, $ref_hash) = @_;

    open(FILE_IN, $csvFile) || die "E: Cannot open csv file $csvFile\n";
    {
	# skip header
	my $line = <FILE_IN>;
	chomp($line);
	($line eq "Disaccharide,phi; psi") || die "E in fillAnglesHash: unexpected header in csv: $line\n";
    }
    while(my $line = <FILE_IN>){
        chomp($line);
        if($line =~ m/^([ab][^,]+),(.+)$/){
            my $disaccharide = $1;
            my $phi_psi_value = $2;
            ${$ref_hash}{$disaccharide} = $phi_psi_value;
        }
    	else {
    	    die "E in fillAnglesHash: cannot parse line from csvFile $csvFile:\n$line\n";
	    }
    }
    close FILE_IN;
}

# take as args a reference of an array which is a stack
# and a string which is a RES bloc line
# extract information from line and put it in hash then put hash in stack
sub treatmentResiduLine{
    (@_ == 2) || die "E: treatmentResiduLine requires 2 args\n";
    my ($ref_pile, $line) = @_;
    
    # %struct will contain information of our residue
    my %struct;
    ($line =~ m/^\d+([bs]):(.+)/) || die "E in treatmentResiduLine: bad line\n$line\n";
    $struct{"type"} = $1;
    $struct{"nom"} = $2;

    push(@{$ref_pile}, \%struct);
}

# take as args a reference of an array which is a stack
# and a string which is a LIN bloc line
# extract information from line and modifie hash in stack
sub treatmentLinkageLine{
    (@_ == 2) || die "E: treatmentLinkageLine requires 2 args\n";
    my ($ref_pile, $line) = @_;

    ($line =~ m/^\d+:(\d+)[a-z]\((\d\+\d)\)(\d+)[a-z]$/) || die "E in treatmentLinkageLine: bad line\n$line\n";
    # extract 2 nodes id and bound value
    my $nodeA = $1;
    my $nodeB = $3;
    my $bound = $2;

    if(${$ref_pile}[$nodeA - 1]{"type"} eq "b" && ${$ref_pile}[$nodeB - 1]{"type"} eq "b"){ # -1 because array starts at 0
        # case A and B are basetypes
        ${$ref_pile}[$nodeB - 1]{"elem_suivant"} = $nodeA;
        ${$ref_pile}[$nodeB - 1]{"liaison_elem_suivant"} = $bound;
    }
    elsif(${$ref_pile}[$nodeA - 1]{"type"} eq "b" && ${$ref_pile}[$nodeB - 1]{"type"} eq "s"){
        # case A is basetype and B is substituent
        # because a basetype can have many substituents, I put information in arrays
        push(@{${$ref_pile}[$nodeA - 1]{"liste_feuilles"}}, $nodeB);
        push(@{${$ref_pile}[$nodeA - 1]{"liste_liaisons_feuilles"}}, $bound);
    }
    else{
        die "E: error in line $line\n";
    }
}

# take as args a reference of an array which is a stack
# and position of element in this stack
# create full polys name of a basetype with his substituents
sub createFullName{
    (@_ == 2) || die "E: createFullName requires 2 args\n";
    my ($ref_pile, $i) = @_;

    my %node = %{${$ref_pile}[$i]};

    ($glycoct2polys{$node{"nom"}}) || die "E: This key '".$node{"nom"}."' doesn't exist in our convert hash\n";
    my $basetypePolysName = $glycoct2polys{$node{"nom"}};
    if(defined $node{"liste_feuilles"}){
        # case basetype has substituent(s)
        foreach my $j (0..$#{$node{"liste_feuilles"}}){
            # for each substituent(s)
            my $indexLeaf = ${$node{"liste_feuilles"}}[$j];
            $indexLeaf--; # because array starts at 0
            my $leafPolysName = $glycoct2polys{${$ref_pile}[$indexLeaf]{"nom"}};
            if($leafPolysName eq "S"){
                my $connexion = $node{"liste_liaisons_feuilles"}[$j]; # eg: 6+1
                my ($parentAttachmentPosition) = ($connexion =~ m/^(\d+)\+\d+$/);
                if($leafPolysName =~ m/_/){
                    # case polys name already contains an '_'
                    $leafPolysName = $parentAttachmentPosition.$leafPolysName;
                }
                else{
                    # case polys name doesn't contains an '_' yet
                    $leafPolysName = "_".$parentAttachmentPosition.$leafPolysName;
                }
            }
            $basetypePolysName = $basetypePolysName.$leafPolysName;
        }
    }
    if($basetypePolysName =~ m/Ido/){
        # case idose: conformation 1C4 or 2S0 doesn't exist in GlycoCT, we use 1C4 conformation by default in polys
        $basetypePolysName = $basetypePolysName."_1C4";
    }
    ${$ref_pile}[$i]{"nom_complet_polys"} = $basetypePolysName;
}
