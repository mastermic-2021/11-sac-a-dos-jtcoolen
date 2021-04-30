# Sac à dos

Le fichier ``input.txt`` contient une clef publique ``[N,[b1,...bn]]`` et
un chiffré ``m_1*b_1+...m_n*b_n mod N``. On demande d'afficher la suite
de 0 et de 1 ``m_1...m_n``.

## Format

Vous devez écrire un programme `main.gp` en Pari/GP dont l'exécution via
```
gp -fq < main.gp
```
affiche (uniquement) la solution cherchée.

taper `make check` permet de vérifier que votre solution est bonne.

Veillez à mettre des commentaires sur votre démarche et sa validité
dans le fichier ``main.gp``.

De manière alternative, vous pouvez écrire un programme `main.c` qui
fasse la même chose.

## Validation

Il reste à faire un commit et à envoyer votre solution pour l'enregistrer
```
git add main.gp
git commit -m 'ma solution'
git push
```

