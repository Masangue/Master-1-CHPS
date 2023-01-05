Exercice 3.  
1.  
Il n'y a pas de manière particulière de déclarer et allouer une matrice avec BLAS ou LAPACK.   

2.  
Cela signifie que les colonnes de la matrice sont stockées les unes à la suite des autres.  

3.  
La 'leading dimension' est le nombre d'éléments à sauter pour passer à la ligne/colonne suivante. Respectivement pour raw major et column major, il vaut généralment le nombre de colonnes et le nombre de lignes.  
Si A en raw major, on peut donner LDA = 2*nb_colonnes, on accède alors à une ligne sur 2.  

4.  
dgbmv fait le produit matrice-vecteur avec une matrice bande général en double précision.  
Elle implémenta la méthode 

5.  
dgbtrf fait la factorisation LU d'une matrice bande général en double précision.  
Elle implémenta la méthode des pivots partiels.

6.  
dgbtrs résout un système d'équations linéaires avec une matrice bande général en double précision. La matrice est une factorisation LU préalablement calculée par dgbtrf.  
Elle implémenta la méthode   

7.  
dgbsv resout un système d'équations linéaires avec une matrice bande général en double précision.
Elle utilise dgbtrf et dgbtrs pour calculer la factorisation LU et résoudre le système d'équations linéaires.  

8.  
Aucune Idée  

Exercice 4.  

Exercice 6.
2.
On compare la factorisation LU et a solutions obtenus par rapport au résultat exact et la version BLAS.
