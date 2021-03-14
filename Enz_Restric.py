################################### choix 1

def Afficher_Dico_Enz():                                                    #fonction qui affiche le dictionnaire                                  
    print("")
    dico_enz={}                                                             #création du dictionnaire d'enzyme                
    lien=open("dico_enz.txt","r")                                           #ouverture du fichier texte en mode lecture (read)    
    for ligne in lien:                                                      #parcours du fichier ligne par ligne 
        (nom_enz,donnees)=ligne.split()                                     #mise en forme de la liste d'enz à afficher
        dico_enz[nom_enz]=donnees
        donnees = donnees.split(",")
        seq_siteR = str(donnees[0])
        index_coupe = int(donnees[1])

        seq1 = seq_siteR[0:index_coupe]
        seq2 = seq_siteR[index_coupe:len(seq_siteR)]
        affichage_siteR = seq1 + "/" + seq2

        print(nom_enz,"(",affichage_siteR,")")                              
    lien.close()                                                            #fermeture du fichier

#################################### choix 2

def Ajout_Enz():                                                            #fonction pour ajouter une enzyme au dictionnaire(fichier texte)
    lien=open("dico_enz.txt","a")                                           #ouverture du fichier texte en mode ajout (append)
    print("L'enzyme doit être de la forme: nomEnz (SiteR,indiceCoupure) ")  
    print("Exemple: AluI (AGCT,2) pour l'enzyme AluI(AG/CT)")
    enz=input("Entrez une nouvelle enzyme: ")
    new_enz="\n"+enz                                                        #permet d'écrire de la nouvelle enzyme sur la ligne suivante du fichier
    lien.write(new_enz)                                                     #ajout de la nouvelle enzyme au fichier
    lien.close()                                                            

#################################### choix 3 
def creation_dico():                                                        #fonction qui initialise le dictionnaire d'enzymes pour tout le programme
    dico_enz={}                     
    lien=open("dico_enz.txt","r")
    for ligne in lien:
        (nom_enz,donnees)=ligne.split()
        dico_enz[nom_enz]=donnees
    return(dico_enz)



def verif_Enz(enzyme,dico):                                                 #fonction pour vérifier si l'enzyme donnée appartient au dictionnaire
    for cle in dico.keys():                                                 #parcours les clés du dictionnaire
        if(enzyme==cle):                                                    #test si le nom de l'enzyme correspond à l'une des clés
            return True

    return False                                                            #si l'on finit la boucle for alors l'enzyme n'est pas dans le dictionnaire


def verif_seq(seq):                                                         #fonction qui vérifie si la séquence entrée est bien de l'ADN
    ARN= False                                                              #bouléenne pour vérifier si la séquence n'est pas de l'ARN
    t=0                                                                     #compteur pour le nombre de T
    for i in range(len(seq)):
        if seq[i] not in "ATGC":                                            #teste si le caractère i de la séquence n'est pas un nucléotide A,T G ou C
            if (seq[i]!= "U"):                                              #teste si le caractère i de la séquence n'est pas un U
                print("La séquence n'est pas constituée que de nucléotides : ATGC") 
                return(False)
            else:                                                           #si le caractère est un U, alors nous somme en présence d'ARN
                ARN= True
        else:
            if seq[i]=="T":                                                 #on vérifie si le caractère est un T pour le compter
                t=t+1
    if (ARN):
        if(t==0):                                                           #Si on a des U et qu'il n'y a pas de T la séquence est de l'ARN
            print("Attention, la séquence entrée est une séquence d'ARN")
            return(False)
        else:                                                               #Si on a des U et des T, il ya eu erreur de saisi (ni ADN,ni ARN)
            print("Erreur, il y a des T et des U dans votre séquence")
            return(False)
    else:                                                                   #Si uniquement ATCG , on a de l'ADN 
        return(True)      


def trouver_siteR(w,seq_siteR):                                             #fonction qui indique si on se trouve sur le site de restriction d'une enzyme
        if(w == seq_siteR):
            return True
        else:
            return False


def trie_liste_frag(liste):                                                 #fonction qui trie selon la taille (décroissant) des éléments comme migration 
    modif=True
    while(modif==True):                                                     #tant qu'on fait une modification on reparcourt la liste
        modif=False
        for i in range(len(liste)-1):
            if(len(liste[i+1])>len(liste[i])):                              #trie à bulle 
               liste[i],liste[i+1]=liste[i+1],liste[i]
               modif=True
               
    return(liste)                                                           #sorti du while => liste triée



def afficheResultDig(enz,liste):                                            #fonction pour afficher le résultat de digestion d'une ou plusieurs enzymes
    trie_liste_frag(liste)                                                  #trie des fragments par taille décroissante
    print("Le nombre de fragment généré par la/les enzyme(s) ",enz,"est ",len(liste))
    print("Voici la liste des fragments obtenues ainsi que leur taille respectives:")
    for i in range(len(liste)):                                             #parcour liste de fragment pour l'afficher
        print("fragment",i+1,":",liste[i],"  taille:",len(liste[i]))


def digestion(ENZ, seq,dico):                                               #fonction de digestion d'une séquence d'ADN par une seule enzyme
    var=dico[ENZ]                                                           #variable pour récupérer les informations du site de coupure de l'enzyme dans le dico (site et indice de coupure)
    var=var.split(",")                                                      #.split pour séparer les données récupéré du dico
    seq_siteR = var[0]                                                      
    index_coupe=int(var[1])
    dernier_coupure=0                                                       #donne l'indice de la dernière coupure par l'enzyme
    l_frag=[]                                                               #liste des fragments qu'on obtiendra après digestion
    for i in range(len(seq)-len(seq_siteR)+1):                              #incrémentation de 1 du début à fin séq ADN (on retire la taille du site de coupure pour ne pas être hors cadre de lecture à la fin)
        cadre=seq[i:i+len(seq_siteR)]                                       #cadre= fenêtre de lecture= taille du site de restriction
        if (trouver_siteR(cadre,seq_siteR)):                                
            l_frag.append(seq[dernier_coupure:i+index_coupe])               #mise à jour liste frag
            dernier_coupure=i+index_coupe                                   #mise à jour du dernier site de coupure
                
        #sortie boucle for
    l_frag.append(seq[dernier_coupure:])
    return l_frag

########################## choix 4

def digestion_multiple(l_enz,seq,dico):                                     #fonction pour gérer les digestions par plusieurs enzymes
    l_frag=[]                                                               #liste final de tous les fragments digérés
    for i in range(len(l_enz)):                                             #parcours les enzymes
        if (i==0):
            l_frag=digestion(l_enz[i],seq,dico)
        else:                               
            l_frag2=[]                                                      #liste intermidiaire de fragments généré par 1 enzyme
            for j in range(len(l_frag)):                                    #parcours les fragments générés par les digestions précédentes
                l_frag2+=digestion(l_enz[i],l_frag[j],dico)                 #ajout des nouveaux fragments digérés
            l_frag=l_frag2   
    return l_frag

##########################
def Menu():                                                                 #fonction principale qui va gérer les différents choix possible du programme
    utilisation=True
    print("Bienvenue sur notre plateforme!")                                #Présentation des fonctionnalités
    print("Les options sont:")
    print("1. Afficher le dictionnaire des enzymes disponibles")
    print("2. Ajouter une enzyme au dictionnaire")
    print("3. Digestion d'une séquence par l'enzyme choisie")
    print("4. Digestion d'une séquence par plusieurs enzymes")
    print("5. Quittez l'interface")
    print("Attention de ne pas taper une lettre sous peine d'erreur de programme!")
    while(utilisation):                                                     #Boucle while pour pouvoir faire plusieurs choix à la suite
        choix=int(input("\n Quel est votre choix: "))                       #Demande du choix de l'utilisateur

        if (choix==1):                                                      #choix 1 pour afficher les enzymes du dictionnaires
            Afficher_Dico_Enz()
            
        elif(choix==2):                                                     #choix 2 pour ajouter une enzyme au dictionnaire
            Ajout_Enz()
            print("Votre enzyme a bien été ajoutée ! ^^ ")
            
        elif(choix ==3):                                                    #choix 3 pour digérer une séquence par 1 enzyme
            dico_enz=creation_dico()                                        #création du dico
            seq=input("Entrez la séquence d'ADN à digérer: ")
            seq=seq.upper()                                                 #on met la sequence donné en majuscule pour la traiter 
            enz=input("Entrez le nom de l'enzyme: ")
            if(verif_Enz(enz,dico_enz) and verif_seq(seq)):                 #vérification si aucune erreur pour lancer la digestion
                l_frag=digestion(enz,seq,dico_enz)
                afficheResultDig(enz,l_frag)
                
            else:                                                           #si au moins une erreur, affichage de message personalisé selon l'erreur
                if(not(verif_Enz(enz,dico_enz))):
                    print("Erreur!!!!!!!!!!!! ",enz," n'est pas dans le dictionnaire!!!!!")
                    print("Vérifiez l'écriture ou ajoutez la à partir du menu")
                if(not(verif_seq(seq))):
                    print("Donnez une séquence d'ADN pour faire une digestion")

        elif (choix == 4):                                                  #choix 4
            dico_enz=creation_dico()
            seq=input("Entrez la séquence d'ADN à digérer: ")
            seq=seq.upper()
            n=int(input("Entrez le nb enzymes pour digestion multiple : ")) #Demande le nombre d'enzyme à entrer
            test_enz=True
            enz=[]                                                          #liste des enzyme donné par l'utilisateur
            bug=[]                                                          #liste des enzymes n'existant pas
            
            for i in range(n):                                              #Entrée des n enzymes
                new_enz= input("Entrez le nom des enzymes : ")
                enz.append(new_enz)
                if(verif_Enz(enz[i],dico_enz)==False):                      #vérification de la validité des enzymes
                    bug.append(enz[i])                                      
                    test_enz=False
    
            test_seq=verif_seq(seq)
    
            if(test_enz and test_seq):                                      #si aucune erreur, on procède à la digestion
                l_frag=digestion_multiple(enz,seq,dico_enz)
                afficheResultDig(enz,l_frag)

            else:                                                           #sinon, affichage de message d'erreur personnalisé selon l'erreur
                if(not(test_enz)):
                    print("!!! Erreur ",bug," n'est pas dans le dictionnaire!!!")
                    print("Vérifiez l'écriture ou ajoutez la à partir du menu")
                if(not(test_seq)):
                    print("Donnez une séquence d'ADN pour faire une digestion")
            
        elif(choix==5):                                                     #choix 5, arrêt de la boucle while et du programme
            utilisation=False
            print("Merci de votre visite, à bientôt !!! <3")
        else:                                                               #Valeur correspondant à aucun choix
            print("Entrée non conforme aux choix")
