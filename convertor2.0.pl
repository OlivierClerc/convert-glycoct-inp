#!/bin/env perl

# Takes as args a csv file and a glycoct file.
# Print on stdout the corresponding inp format for POLYS tool.

# use
# perl convertor2.0.pl phi_psi_all_disacc.csv HP_dp02_0001.glycoct > HP_dp02_0001.inp

use strict;
use warnings;

(@ARGV == 2) || die "E: This script needs 2 arg: a csv file with phi-psi value and the glycoct file.\n";
my ($csvFile, $glycoctFile) = @ARGV;


# @pile will contain information of our glycoct file
# each element of @pile is a hash (or residue in our glycoct file)
# all hashs have keys: 'type' (type of residue: basetype, substituent or repeating unit) and 'nom' (name of residue)
# somes hashs have keys: 'elem_suivant' (index of next residue in @pile which is a basetype or repeating unit),
#                        'liaison_elem_suivant' (bond with next residue which is a basetype or repeating unit),
#                        'liste_feuilles' (table of index of next residues in @pile which are substituents),
#                        'liste_liaisons_feuilles' (table of bond with next residues which are substituents),
#                        'nom_complet_inp' (inp name of basetype residue),
# repeating unit hash have keys: 'liaison_rep' (bond with next repetition),
#                                'bornes' (lower and upper bound of repetition, '-1' value mean infinity),
#                                'indice_elem_debut' (index of first residue in repeating unit),
#                                'indice_elem_fin' (index of last residue in repeating unit),
#                                'sous_pile' (substructure as @pile but only with basetypes or substituents)
my @pile;

# %glycoct2inp holds for each residue in glycoct format, corresponding residue in inp format
# key is glycoct format, value is inp format
my %glycoct2inp = ("a-dglc-HEX-1:5" => "aDGlcp",
             "a-lara-HEX-1:5|4:d|6:a" => "aL4deoxyIdopA",
		     "a-lido-HEX-1:5|6:a" => "aLIdopA",
		     "b-dgal-HEX-1:5" => "bDGalp",
		     "b-dglc-HEX-1:5" => "bDGlcp",
             "b-dglc-HEX-1:5|6:a" => "bDGlcpA",
		     "b-dxyl-HEX-1:5|4:d|6:a" => "bD4deoxyGlcpA",
		     "amino" => "N",
		     "n-acetyl" => "NAc",
		     "n-sulfate" => "NS",
		     "sulfate" => "S",
             "methyl" => "Me"); # maybe delete because no case with methyl in POLYS

# in case of repeating unit, %rep will contain position of repeating unit in @pile
# key is id (number) of repeating unit, value is index in @pile
my %rep;

########################################################
# 1st part:read csv file and glycoct file and store data
########################################################
# %phiPsiAngles will contain information of our csv file
# key is disaccharide in inp format, value is phi-psi angle value
my %phiPsiAngles;
&fillAnglesHash($csvFile, \%phiPsiAngles);

open(FILE_IN, $glycoctFile) || die "E: Cannot open glycoct file $glycoctFile\n";
# $section: "RES" or "LIN", says which section we're in
my $section;
# $option: "REP" or "", says which case we're in
my $option;
# in case of repetition, $id_rep is the id (number) of the current repeating unit
my $id_rep;
while(my $line = <FILE_IN>){
    chomp($line);
    if($line eq "RES"){
	    ($section && $option ne "REP") && die "E: Found RES header but section already defined as $section , shouldn't happen!\n";
	    $section = "RES";
    }
    elsif($line eq "LIN"){
	    ($section eq "RES") || die "E: found LIN header but we were not in RES section? section was $section\n";
	    $section = "LIN";
    }
    elsif($line eq "REP"){
        ($section eq "LIN" || $section eq "RES") || die "E: found REP header not after a RES or LIN section, section was $section\n";
        $option = "REP";
    }
    elsif($line =~ m/^\d+[brs]:.+/){ # RES bloc line with basetype, substituent or repeating unit (eg: 1b:a-dglc-HEX-1:5)
	    ($section eq "RES") || die "E: found RES line but not in RES section:\n$line\n";
        if(! $option){
            &treatmentResiduLine(\@pile, $line);
            if($line =~ m/^(\d+)r:r(\d+)$/){ # eg: 3r:r1
                # case repeat
                $rep{$2} = $1;
            }
        }
        elsif($option eq "REP"){
            &treatmentResiduLine(\@pile, $line, $rep{$id_rep});
        }
    }
    elsif($line =~ m/^\d+:\d+[odhn]\(\d\+\d\)\d+[odhn]$/){ # LIN bloc line (eg: 1:1d(2+1)2n)
	    ($section eq "LIN") || die "E: found LIN line but not in LIN section:\n$line\n";
        if(! $option){
            &treatmentLinkageLine(\@pile, $line);
        }
        elsif($option eq "REP"){
            &treatmentLinkageLine(\@pile, $line, $rep{$id_rep});
        }
    }
    elsif($line =~ m/^REP(\d+):\d+[odhn]\(\d\+\d\)\d+[odhn]=(-1|\d+)-(-1|\d+)$/){ # REP bloc line (eg: REP1:7o(4+1)5d=1--1)
        $id_rep = $1;
        &treatmentRepeatLine(\@pile, $line, $rep{$id_rep});
    }
    elsif($line =~ m/^\d+a:.+/){ # RES bloc line with alternative unit (eg: 1a:a1)
        die "E: This program can't process alternative unit\n";
    }
    elsif($line =~ m/^UND/){ # UND bloc line
        die "E: This program can't process underdetermined structures\n";
    }
    elsif($line =~ m/^ISO/){ # ISO bloc line
        die "E: This program can't process ISO section\n";
    }
    elsif($line eq ""){
        # sometimes we can have empty lines at the end of the file
    }
    else{
	    die "E parsing glycoctFile $glycoctFile: unexpected line:\n$line\n";
    }
}
close FILE_IN;
##############################################################################
# 2nd part:for each basetype node, creation of full inp name with substituents
##############################################################################
foreach my $i (0..$#pile){
    if(($pile[$i]) && ($pile[$i]{"type"} eq "b")){
        # case basetype
        &createFullName(\@pile, $i);
    }
    elsif(($pile[$i]) && ($pile[$i]{"type"} eq "r")){
        # repeat case
        foreach my $j (0..$#{$pile[$i]{"sous_pile"}}){
            if(($pile[$i]{"sous_pile"}[$j]) && ($pile[$i]{"sous_pile"}[$j]{"type"} eq "b")){
                &createFullName(\@pile, $i, $j);
            }
        }
    }
}
#########################
# 3rd part:write inp file
#########################
print("PRIMARY\n");
while(scalar(@pile)){
    my $lastElemPile = pop(@pile);
    if((${$lastElemPile}{"type"}) && (${$lastElemPile}{"type"} eq "b")){
        # case basetype
        my $printLine = " <".${$lastElemPile}{"nom_complet_inp"}.">"; # eg: <bDGlcpA>
        if(${$lastElemPile}{"elem_suivant"}){
            # case our basetype has a next basetype or repeating unit
            my $indexNextElem = ${$lastElemPile}{"elem_suivant"};
            my $connexion = ${$lastElemPile}{"liaison_elem_suivant"}; # eg: 4+1
            ($connexion =~ m/^(\d+)\+(\d+)$/) || die "E in connexion of ".${$lastElemPile}{"nom_complet_inp"}.":\n$connexion\n";
            my $parentAttachmentPosition = $1;
            my $childAttachmentPosition = $2;
            my $key;
            if($pile[$indexNextElem]{"type"} eq "b"){
                # normal case
                $key = ${$lastElemPile}{"nom_complet_inp"}." ".$childAttachmentPosition."-".$parentAttachmentPosition." ".$pile[$indexNextElem]{"nom_complet_inp"}; # eg: aLIdopA_2S_1C4 1-4 aDGlcpN_6S
            }
            elsif($pile[$indexNextElem]{"type"} eq "r"){
                # repeat case
                my $indexNextElemInRepeatUnit = $pile[$indexNextElem]{"indice_elem_debut"};
                $key = ${$lastElemPile}{"nom_complet_inp"}." ".$childAttachmentPosition."-".$parentAttachmentPosition." ".$pile[$indexNextElem]{"sous_pile"}[$indexNextElemInRepeatUnit]{"nom_complet_inp"}; # eg: aLIdopA_2S_1C4 1-4 aDGlcpN_6S
            }
            ($phiPsiAngles{$key}) || die "E: This key '".$key."' doesn't exist in our phi-psi value hash\n";
            $printLine = $printLine." ( ".$childAttachmentPosition.":".$parentAttachmentPosition."; ".$phiPsiAngles{$key}.")"; # eg: <bDGlcpA> ( 1:4; 40.00; -110.00)
        }
        print($printLine."\n");
    }
    elsif((${$lastElemPile}{"type"}) && (${$lastElemPile}{"type"} eq "r")){
        # repeat case
        my $printLine = "["; # repeat start
        my @currentTab = @{@{$lastElemPile}{"sous_pile"}};
        while(scalar(@currentTab)){
            my $lastElemRepPile = pop(@currentTab);
            if((${$lastElemRepPile}{"type"}) && (${$lastElemRepPile}{"type"} eq "b")){
                # case basetype
                $printLine = $printLine." <".${$lastElemRepPile}{"nom_complet_inp"}.">"; # eg: <bDGlcpA>
                if(${$lastElemRepPile}{"elem_suivant"}){
                    # case our basetype has a next basetype
                    my $connexion = ${$lastElemRepPile}{"liaison_elem_suivant"}; # eg: 4+1
                    ($connexion =~ m/^(\d+)\+(\d+)$/) || die "E in connexion of ".${$lastElemRepPile}{"nom_complet_inp"}.":\n$connexion\n";
                    my $parentAttachmentPosition = $1;
                    my $childAttachmentPosition = $2;
                    my $indexNextElem = ${$lastElemRepPile}{"elem_suivant"};
                    my $key = ${$lastElemRepPile}{"nom_complet_inp"}." ".$childAttachmentPosition."-".$parentAttachmentPosition." ".$currentTab[$indexNextElem]{"nom_complet_inp"}; # eg: aLIdopA_2S_1C4 1-4 aDGlcpN_6S
                    ($phiPsiAngles{$key}) || die "E: This key '".$key."' doesn't exist in our phi-psi value hash\n";
                    $printLine = $printLine." ( ".$childAttachmentPosition.":".$parentAttachmentPosition."; ".$phiPsiAngles{$key}.")\n"; # eg: <bDGlcpA> ( 1:4; 40.00; -110.00)                    
                }
            }
        }
        $printLine = $printLine."\n]"; # repeat end
        # number of repeat
        my $numberOfRepeat = ${$lastElemPile}{"bornes"};
        $numberOfRepeat =~ m/^(-1|\d+)-(-1|\d+)$/;
        my $lowerRange = $1;
        my $upperRange = $2;
        if($lowerRange eq "-1"){ # range equal to -1 means infinity, we decided to replace by 1
            $lowerRange = 1;
        }
        $printLine = $printLine.$lowerRange; # we decided to keep only the lower range because there is only one figure in inp repeat
        # connexion repeat
        my $connexion = ${$lastElemPile}{"liaison_rep"};
        ($connexion =~ m/^(\d+)\+(\d+)$/) || die "E in connexion of ".${$lastElemPile}{"nom_complet_inp"}.":\n$connexion\n";
        my $parentAttachmentPosition = $1;
        my $childAttachmentPosition = $2;
        my $indexNextElem = ${$lastElemPile}{"indice_elem_debut"};
        my $indexPreviousElem = ${$lastElemPile}{"indice_elem_fin"};
        my $key = ${@{$lastElemPile}{"sous_pile"}}[$indexPreviousElem]{"nom_complet_inp"}." ".$childAttachmentPosition."-".$parentAttachmentPosition." ".${@{$lastElemPile}{"sous_pile"}}[$indexNextElem]{"nom_complet_inp"}; # eg: aDGlcpN_6S 1-4 aLIdopA_2S_1C4
        ($phiPsiAngles{$key}) || die "E: This key '".$key."' doesn't exist in our phi-psi value hash\n";
        $printLine = $printLine." ( ".$childAttachmentPosition.":".$parentAttachmentPosition."; ".$phiPsiAngles{$key}.")"; # eg: ]1 ( 1:4; 40.00; -110.00)
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
    my ($csvFile, $hashRef) = @_;

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
            ${$hashRef}{$disaccharide} = $phi_psi_value;
        }
    	else {
    	    die "E in fillAnglesHash: cannot parse line from csvFile $csvFile:\n$line\n";
	    }
    }
    close FILE_IN;
}

# take as args a reference of an array which is a stack,
# a string which is a RES bloc line and sometimes the position of
# a repeating unit in the stack
# extract information from line and put it in hash then put hash in stack
sub treatmentResiduLine{
    ((@_ == 2) || (@_ == 3)) || die "E: treatmentResiduLine requires 2 or 3 args\n";
    my ($stackRef, $line, $positionOfRepInStack) = @_;
    
    # %struct will contain information of our residue
    my %struct;
    ($line =~ m/^(\d+)([brs]):(.+)/) || die "E in treatmentResiduLine: bad line\n$line\n";
    my $residuId = $1;
    $struct{"type"} = $2;
    $struct{"nom"} = $3;

    if($positionOfRepInStack){
        # repeat case
        ${$stackRef}[$positionOfRepInStack]{"sous_pile"}[$residuId] = \%struct;
    }else{
        # normal case
        ${$stackRef}[$residuId] = \%struct;
    }
}

# take as args a reference of an array which is a stack,
# a string which is a LIN bloc line and sometimes the position of
# a repeating unit in the stack
# extract information from line and modifie hash in stack
sub treatmentLinkageLine{
    ((@_ == 2) || (@_ == 3)) || die "E: treatmentLinkageLine requires 2 or 3 args\n";
    my ($stackRef, $line, $positionOfRepInStack) = @_;

    ($line =~ m/^\d+:(\d+)[a-z]\((\d\+\d)\)(\d+)[a-z]$/) || die "E in treatmentLinkageLine: bad line\n$line\n";
    # extract 2 nodes id and bond value
    my $nodeA = $1;
    my $nodeB = $3;
    my $bond = $2;

    # $workingStack is an array corresponding to stack in reference or if it exists, the substack
    my $workingStack;
    if($positionOfRepInStack){
        # repeat case
        $workingStack = ${$stackRef}[$positionOfRepInStack]{"sous_pile"};
    }else{
        # normal case
        $workingStack = $stackRef;
    }
    if(${$workingStack}[$nodeA]{"type"} =~ m/^[br]{1}$/ && ${$workingStack}[$nodeB]{"type"} =~ m/^[br]{1}$/){
        # case A and B are basetypes or repeating unit
        ${$workingStack}[$nodeB]{"elem_suivant"} = $nodeA;
        ${$workingStack}[$nodeB]{"liaison_elem_suivant"} = $bond;
    }elsif(${$workingStack}[$nodeA]{"type"} eq "b" && ${$workingStack}[$nodeB]{"type"} eq "s"){
        # case A is basetype and B is substituent
        # because a basetype can have many substituents, I put information in arrays
        push(@{${$workingStack}[$nodeA]{"liste_feuilles"}}, $nodeB);
        push(@{${$workingStack}[$nodeA]{"liste_liaisons_feuilles"}}, $bond);
    }else{
        die "E: error in line $line\n";
    }
}

# take as args a reference of an array which is a stack,
# a string which is a REP bloc line and the position of a
# repeating unit in the stack
# extract information from line and modifie hash in stack
sub treatmentRepeatLine{
    (@_ == 3) || die "E: treatmentRepeatLine requires 3 args\n";
    my ($stackRef, $line, $positionOfRepInStack) = @_;

    ($line =~ m/^REP\d+:(\d+)[odhn]\((\d\+\d)\)(\d+)[odhn]=(-1|\d+)-(-1|\d+)$/) || die "E in treatmentRepeatLine: bad line\n$line\n";
    my $repetitionBond = $2;
    my $parentResiduNumber = $1;
    my $childResiduNumber = $3;
    my $range = $4."-".$5;

    ${$stackRef}[$positionOfRepInStack]{"liaison_rep"} = $repetitionBond;
    ${$stackRef}[$positionOfRepInStack]{"bornes"} = $range;

    ${$stackRef}[$positionOfRepInStack]{"indice_elem_debut"} = $parentResiduNumber;
    ${$stackRef}[$positionOfRepInStack]{"indice_elem_fin"} = $childResiduNumber;
}

# take as args a reference of an array which is a stack,
# a position of element in this stack and sometimes a position
# of element in substack
# create full inp name of a basetype with his substituents
sub createFullName{
    ((@_ == 2) || (@_ == 3)) || die "E: createFullName requires 2 or 3 args\n";
    my ($stackRef, $indexInStack, $indexInSubstack) = @_;

    # $index is the index in stack or when it exists, the index in substack
    my $index;
    # @currentTab is an array corresponding of the stack in reference or when it exists, the substack in the stack
    my @currentTab;
    if($indexInSubstack){
        # repeat case
        @currentTab = @{${$stackRef}[$indexInStack]{"sous_pile"}};
        $index = $indexInSubstack;
    }else{
        # normal case
        @currentTab = @{$stackRef};
        $index = $indexInStack;
    }
    my %node = %{$currentTab[$index]}; # for simplification
    ($glycoct2inp{$node{"nom"}}) || die "E: This key '".$node{"nom"}."' doesn't exist in our convert hash\n";
    my $basetypeInpName = $glycoct2inp{$node{"nom"}};
    if($node{"liste_feuilles"}){
        # case basetype has substituent(s)
        foreach my $k (0..$#{$node{"liste_feuilles"}}){
            # for each substituent(s)
            my $indexLeaf = ${$node{"liste_feuilles"}}[$k];
            ($glycoct2inp{$currentTab[$indexLeaf]{"nom"}}) || die "E: This key '".$currentTab[$indexLeaf]{"nom"}."' doesn't exist in our convert hash\n";
            my $leafInpName = $glycoct2inp{$currentTab[$indexLeaf]{"nom"}};
            if($leafInpName eq "S"){
                my $connexion = $node{"liste_liaisons_feuilles"}[$k]; # eg: 6+1
                my ($parentAttachmentPosition) = ($connexion =~ m/^(\d+)\+\d+$/);
                if($basetypeInpName =~ m/_/){
                    # case inp name already contains an '_'
                    $leafInpName = $parentAttachmentPosition.$leafInpName;
                }else{
                    # case inp name doesn't contains an '_' yet
                    $leafInpName = "_".$parentAttachmentPosition.$leafInpName;
                }
            }
            $basetypeInpName = $basetypeInpName.$leafInpName;
        }
    }
    if($basetypeInpName =~ m/aLIdo/){
        # case idose: conformation 1C4 or 2S0 doesn't exist in GlycoCT, we use 1C4 conformation by default in inp
        $basetypeInpName = $basetypeInpName."_1C4";
    }
    $currentTab[$index]{"nom_complet_inp"} = $basetypeInpName;
}
