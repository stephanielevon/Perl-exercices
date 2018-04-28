# Nom du programme : Exercie_2.pl
# Auteur : Stéphanie Levon
# But du programme : A partir d'un fichier VCF, représenter l'ensemble des SNPs répertoriés sous forme de camemberts; D'une part pour l'ensemble des chromosomes et d'autre part pour chaque chromosome. 
use strict;
use warnings;
use Statistics::R;

# Demander à l'utilisateur le fichier VCF à analyser
print "Entrez le nom de fichier VCF a ouvrir : \n";
my$nom_fichier=<>;
chomp $nom_fichier; # Retirer le saut de ligne
my%tableau=mise_en_forme($nom_fichier); # Créer un tableau associatif qui contient les données importantes du fichier vcf (cf. ci-dessous)

# Créer un fichier sous R et le nommer avec une extension ".pdf"
my@liste_nom_fichier=split(/\./, $nom_fichier); # Créer une liste avec le nom de fichier à ouvrir. Séparer le nom de fichier en deux sous-entités : avant et après le "." 
my $R = Statistics::R->new(); # Créer un nouveau fichier sous R grâce au module "statistics::R" 
$R->run(qq`pdf("$liste_nom_fichier[0].pdf")`); # Lancer la commande R qui va créer notre fichier .pdf contenant toutes les représentations graphiques
foreach my$elmt(sort keys %tableau) { # Pour chaque élèment du tableau contenant les données, trier les clés du tableau soit les chromosomes
    my@valeurs_R=(); # Valeur des différents SNPs
    my@nom_R=(); # Nom du changement de base (AT, AG, AC, etc.)
    my@couleurs_R=(); # Couleur associée au changement de base
    foreach my$second_elmt(sort keys %{$tableau{$elmt}}) { # Pour chaque élèment du second tableau

        push @valeurs_R, $tableau{$elmt}{$second_elmt}; # Ajouter les valeurs des SNP dans une liste
        push @nom_R, "\"".$second_elmt."\""; # Ajouter le nom des SNP associés aux valeurs, sans oublier les guillemets
        push @couleurs_R, "\"".changement_couleurs($second_elmt)."\""; # Faire la même chose pour la couleur associée à chaque SNP
    }
    # Enregistrer les élèments des trois listes dans leurs vecteurs correspondant pour R, en les séparant par une virgule
    my $valeurs="valeurs <- c(". join(",", @valeurs_R) .")" ; # Vecteur R contenant les valeurs
    my $noms="noms <- c(". join(",", @nom_R) .")" ; # Vecteur R contenant les noms associés aux valeurs
    my $couleurs="c(". join(",", @couleurs_R) .")" ; # Vecteur R contenant les couleurs associées aux noms
    # Lancer l'interaction avec R, charger les 3 vecteurs et on fait un PIE avec (fonction pour les camemberts)
    # Calculer le pourcentage de chaque SNP dans chaque chromosome
    $R->run(<<"CMD"
    	$valeurs 
    	$noms
    	pourcentage <- round(valeurs/sum(valeurs)*100) 
    	noms <- paste(noms, pourcentage)
    	noms <- paste(noms,"%",sep="")
    	pie($valeurs, labels = noms, main="$elmt", col = $couleurs)
CMD
);
}
# Quitter le second tableau
$R->run(q`dev.off()`); # Indiquer à R qu'on a fini de travailler dans le fichier .pdf
$R->stopR(); # Fermer la session R et l'interaction avec Perl
print("Votre fichier a bien ete crée : liste_nom_fichier[0].pdf \n"); # Informer l'utilisateur de la réation du fichier pdf contenant les représentations graphiques




######################################################################################################## 
sub mise_en_forme { # Fonction qui va parser le fichier .vcf et enregistrer les informations dans un tableau associatif
    my($nom_fichier)=@_; # Enregistrer le nom du fichier dans une variable
    unless (open(FIC_VCF, $nom_fichier)){  # Ouvrir le fichier .vcf
        print "Impossible d'ouvrir le fichier $nom_fichier!\n"; # Afficher un message d'erreur lorsque le fichier ne peut pas etre ouvert
        exit; 
        }
# Récuperer les lignes d'informations utiles
    my@fichier_vcf=<FIC_VCF>; # Enregistrer le fichier dans un tableau appelé VCF
    close FIC_VCF; # Fermer le descripteur de fichier

    foreach my$ligne(@fichier_vcf){ # Pour chaque ligne du tableau
        if($ligne=~/^#/){ # Si la ligne commence par le caractère "#"
            next; # Passer à la ligne suivante
        } 
        
        elsif($ligne=~/^\s*$/ ){ # Si la ligne est une suite d'espace, du debut à la fin
            next; # Passer à la ligne suivante
        } 
        
        else{ # Si ce n'est pas le cas et donc si c'est une information utile
            my @vcf = split("\t",$ligne); # Chaque élèment séparé par une tabulation se retrouve dans une colonne du tableau associatif        
         
            if ((length($vcf[3]) == 1)&&(length($vcf[4]) == 1)){ # Si la longueur de la chaine contenue en colonne 3 et 4 est de 1 (changement d'un seul nucléotide)
                my$SNP=uc($vcf[3].$vcf[4]); # Alors créer une variable "SNP" dont le nom contiendra deux caractère : l'ancienne base (colonne 3) et la nouvelle base (colonne 4)
                $tableau{$vcf[0]}{$SNP}++; # Associer cette variable SNP au chromosome sur lequel il appartient (colonne 0)
                $tableau{'Total'}{$SNP}++; # Ajouter cette variable SNP sans tenir compte de sa localisation pour établir le camembert total
            }
        }
    }
    return(%tableau); # Retourner les informations sous la forme d'un tableau associatif
      
}

# Fonction qui va attribuer une couleur invariable à chaque changement de base
sub changement_couleurs { 
    my($SNP)=@_; # Récupérer le SNP dans une variable
    my%couleur_SNP = ( # Créer un tableau associatif ou une couleur est associée à chaque SNP
        'TG' => 'coral','TC' => 'coral3','TA' => 'coral4', 
        'AT' => 'cadetblue1','AG' => 'cadetblue3','AC' => 'cadetblue4', 
        'CT' => 'darkseagreen1','CG' => 'darkseagreen3','CA' => 'darkseagreen4', 
        'GT' => 'darkgoldenrod1','GC' => 'darkgoldenrod3','GA' => 'darkgoldenrod4', 
     );
    return $couleur_SNP{$SNP}; 
}
<>; 