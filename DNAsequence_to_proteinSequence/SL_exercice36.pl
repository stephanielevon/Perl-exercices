use strict;
use warnings;
# Nom du programme : SL_exercice36.pl
# Auteur : Stéphanie Levon
# But du programme : Traduire une séquence nucléotidique en une séquence protéique suivant les six cadres de lecture possible. 

# Demander à l'utilisateur la séquence nucléotidique au format Fasta à traduire en protéine
print "Entrez le nom de fichier a ouvrir : \n";
my$nom_fichier=<>;
chomp $nom_fichier; # Chomp  : retrait du saut de ligne

# Mettre en forme la séquence (cf. ci-dessous)
my$ADN=mise_en_forme($nom_fichier);

# Traduire la séquence selon les trois cadres du brin inchangé (cf. ci-dessous)
my$proteine_cadre1= traduire_cadre($ADN, 0);#cadre 1
my$proteine_cadre2= traduire_cadre($ADN, 1);#cadre 2
my$proteine_cadre3= traduire_cadre($ADN, 2);#cadre 3

# Obtenir le brin inverse
my$brin_inverse=$ADN;
$brin_inverse=reverse $brin_inverse;
$brin_inverse=~ tr/ACGTacgt/TGCAtgca/;

# Traduire la séquence selon les trois cadres du brin inverse
my$proteine_cadre4= traduire_cadre($brin_inverse, 0);#cadre 1
my$proteine_cadre5= traduire_cadre($brin_inverse, 1);#cadre 2
my$proteine_cadre6= traduire_cadre($brin_inverse, 2);#cadre 3

# Créer un fichier contenant les six cadres de lecture
new_fasta($proteine_cadre1,$proteine_cadre2,$proteine_cadre3,$proteine_cadre4,$proteine_cadre5,$proteine_cadre6);

################################################################################################################

### Créer une fonction qui permet la mise en forme de la séquence  
sub mise_en_forme {
    my($nom_fichier)=@_; # Récupèrer le nom du fichier à convertir
    unless (open(FIC_PROTEINE, $nom_fichier)){  # Ouvrir le fichier FASTA et transférer les données dans la liste
        print "Impossible d'ouvrir le fichier $nom_fichier!\n"; # Afficher un message d'erreur lorsque le fichier ne peut pas etre ouvert
        exit; 
        }

    my@proteine=<FIC_PROTEINE>; # Enregistrer le fichier de séquence ADN dans une liste appelée protéine
    close FIC_PROTEINE; # Fermer le descripteur de fichier

   
    my$ADN=''; # Creer la variable ADN vide
    foreach my$ligne(@proteine){ # Pour chaque ligne du contenu 
        if($ligne=~/^>/){ # Si la ligne commence par le caractère ">"
            next; # Passer à la ligne suivante
        } 
        
        elsif($ligne=~/^\s*$/ ){ # Si la ligne est une suite d'espace, du debut à la fin
            next; # Passer à la ligne suivante
        } 
        
        else{ # Si ce n'est pas le cas et donc si c'est une sequence ADN, ajouter chaque ligne dans la variable $ADN
            $ADN .=$ligne;
        } 
    }
    
    $ADN=~ s/\s//g; # Supprime les espaces eventuels de la sequence 
    return($ADN); 
      
}

### Créer une fonction qui permet de traduire la séquence en protéine
sub traduire_cadre {
    my($sequence,$debut)=@_; # Récupérer la séquence et son cadre de lecture
    my$taille=length($sequence); # Définire la taille de la séquence à traduire

    my$cadre=substr($sequence,$debut,$taille-$debut); # Paramètres de la fonction substr (séquence, intervalle sur lequel on effectue la traduction)

    my$proteine=""; # Créer la variable vide qui contiendra la protéine
    for (my$i=0;$i<(length($cadre)-2);$i+=3){ # Compteur de base est de 0, traduit jusqu'à ce que l'on arrive à taille de la séquence -2. Compteur de 3 en 3. 
        my $triplet=substr($cadre,$i,3); # Sous-chaine de caractère qui commence à i et d'une taille de 3
        $proteine .= code_genetique($triplet); # Appel le dictionnaire contenant le code génétique
    }
    return $proteine;
}

### Créer le dictionnaire contenant les acides aminés. A chaque triplet de nucléotide est associé l'acide aminé correspondant
sub code_genetique {
    my($codon)=@_;
    my %code_genetique =(
        "GCT"=>"A", # Alanine
        "GCC"=>"A", # Alanine
        "GCA"=>"A", # Alanine
        "GCG"=>"A", # Alanine
        "TTA"=>"L", # Leucine
        "TTG"=>"L", # Leucine
        "CTT"=>"L", # Leucine
        "CTC"=>"L", # Leucine
        "CTA"=>"L", # Leucine
        "CTG"=>"L", # Leucine
        "CGT"=>"R", # Arginine
        "CGC"=>"R", # Arginine
        "CGA"=>"R", # Arginine
        "CGG"=>"R", # Arginine
        "AGA"=>"R", # Arginine
        "AGG"=>"R", # Arginine
        "AAA"=>"K", # Lysine
        "AAG"=>"K", # Lysine
        "AAT"=>"N", # Asparagine
        "AAC"=>"N", # Asparagine
        "ATG"=>"M", # Méthionine
        "GAT"=>"D", # Acide aspartique
        "GAC"=>"D", # Acide aspartique
        "TTT"=>"F", # Phénylalanine
        "TTC"=>"F", # Phénylalanine
        "TGT"=>"C", # Cystéine
        "TGC"=>"C", # Cystéine
        "CCT"=>"P", # Proline
        "CCC"=>"P", # Proline
        "CCA"=>"P", # Proline
        "CCG"=>"P", # Proline
        "CAA"=>"Q", # Glutamine
        "CAG"=>"Q", # Glutamine
        "TCT"=>"S", # Sérine
        "TCC"=>"S", # Sérine
        "TCA"=>"S", # Sérine
        "TCG"=>"S", # Sérine 
        "AGT"=>"S", # Sérine
        "AGC"=>"S", # Sérine 
        "GAA"=>"E", # Acide glutamique
        "GAG"=>"E", # Acide glutamique 
        "ACT"=>"T", # Thréonine
        "ACC"=>"T", # Thréonine
        "ACA"=>"T", # Thréonine
        "ACG"=>"T", # Thréonine
        "GGT"=>"G", # Glycine
        "GGC"=>"G", # Glycine
        "GGA"=>"G", # Glycine
        "GGG"=>"G", # Glycine
        "TGG"=>"W", # Tryptophane
        "CAT"=>"H", # Histidine
        "CAC"=>"H", # Histidine
        "TAT"=>"Y", # Tyrosine
        "TAC"=>"Y", # Tyrosine
        "ATT"=>"I", # Isoleucine
        "ATC"=>"I", # Isoleucine
        "ATA"=>"I", # Isoleucine
        "GTT"=>"V", # Valine
        "GTC"=>"V", # Valine
        "GTA"=>"V", # Valine
        "GTG"=>"V", # Valine
        "TAG"=>"*", # Codon STOP
        "TGA"=>"*", # Codon STOP
        "TAA"=>"*"); # Codon STOP
    return $code_genetique{$codon};
}
    
### Créer une fonction qui, une fois la traduction effectuée affiche la sortie dans un fichier au format Fasta. 
sub new_fasta {

    my($proteine_cadre1,$proteine_cadre2,$proteine_cadre3,$proteine_cadre4,$proteine_cadre5,$proteine_cadre6)=@_; 
    my@nofasta=split(/\./, $nom_fichier); # Se sert du nom de fichier entré par l'utilisateur pour nommer le nouveau fichier
    my$rendu = $nofasta[0]."-sixcadres.fasta"; # Ajoute le suffixe "_sixcadres.fasta" au nom de fichier


    open(FSOR,">$rendu"); #On créer le fichier avec le nom correspondant
    
    my@cadre=($proteine_cadre1,$proteine_cadre2,$proteine_cadre3,$proteine_cadre4,$proteine_cadre5,$proteine_cadre6);
    for (my$i=0;$i<6;$i++){ # Boucle pour chaque cadre
        print FSOR ">------------Cadre ".($i+1)." --------------\n";
        for (my$j=0;$j<length($cadre[$i]);$j+=70) { # Boucle pour chaque ligne du cadre correspondant, afficher 70 acides aminés par ligne 
		my$sous_chaine=substr($cadre[$i],$j,70);
        	print FSOR "$sous_chaine\n";
        }
    }
        
    close FSOR; # Fermer le fichier

    print "Un fichier Fasta contenant les six cadres de lectures a été crée sous le nom \"$rendu\"\n"; 
   }
<>;
