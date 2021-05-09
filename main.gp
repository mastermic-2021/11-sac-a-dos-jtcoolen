\\ Comme vu en cours, on va appliquer l'algorithme LLL pour résoudre l'égalité du sac à dos
\\ et ainsi retrouver le clair.

\\ Le chiffré c étant pris modulo N, il existe un k entier naturel pour lequel sum_i m[i] * k[i] = c+kN (*).
\\ On va donc prendre comme base du réseau:
\\ (1, 0, ..., 0, B k[1])
\\ (0, 1, ..., 0, B k[2])
\\ ....
\\ (0, 0, ..., 1, B k[140])
\\ (-1/2, -1/2, ..., -1/2, -B (c + kN))
\\ pour B > sqrt(length(k))

\\ On vérifie que la solution (m[1],..., m[140]) in {0,1}^140 recherchée est bien 
\\ dans le réseau pour un certain k tel que (*) :
\\ [ Id,(-1/2); (B k[i])_i, -B (c+ kN) ] * (m[1], ..., m[140], 1)~
\\ = (m[1]-1/2, ..., m[140]-1/2, (sum_i B k[i] * m[i]) - B (c+kN))~
\\ = (m[1]-1/2, ..., m[140]-1/2, 0)~          (**)

\\ On note que la fonction LLL de PARI/GP n'est pas très bien implémentée, car elle ne trouve pas directement
\\ la solution. On va donc réduire le problème aux poids k[2], ..., k[140] en fixant le premier bit m[1], en
\\ espérant que LLL va trouver une solution.

printmessage(m) = print(concat([Str(x)|x<-m]));
\\ La solution doit être translatée de 1/2 d'après (**)
shiftvec(v) = res=[x+1/2|x<-v]; res[length(v)] -= 1/2; res;

[N,k,c] = readvec("input.txt");
\\ N : modulo ; k : clé publique ; c : chiffré

\\ On récupère la valeur du premier poids
k1 = k[1];
\\ On raccourci la clé en ôtant le premier poids
k_red = k[2..-1];

\\ La dimension de la matrice dont les colonnes forment la base du réseau est (n+1)x(n+1) où n=length(k)
l = length(k_red) + 1;

\\ la constante multiplicative B est prise > sqrt(n) pour "forcer" l'équation sac à dos
B = ceil(sqrt(l * 2^l));

\\ construction de la base du réseau
m = matid(l);
m[l, l] = - B * c;
for(i = 1, l - 1, m[l, i] = B * k_red[i]);
for(i = 1, l - 1, m[i, l] = -1/2); \\ cette constante "translate" le réseau pour privilégier les m[i] dans {0,1}, plutôt que {-1,0,1}

\\ teste que la colonne satisfait l'équation sac à dos avec le chiffré comme capacité (mod N)
check_eq(col) = ((col * k~) % N == c);

\\ vérifie que la colonne est faite de deux éléments répétés. On attend 0 et une constante lambda
check_col(c) = (length(Set(c)) == 2 && setsearch(Set(c), 0)); \\ On passe c dans un ensemble (Set), qui conserve les éléments distincts

\\ attaque "low density", possible car la densité est 0.94 < 0.9408
\\ cf. M. J. Coster, A. Joux, B. A. LaMacchia, A. M. Odlyzko, C. P. Schnorr, and J. Stern,
\\ “Improved low-density subset sum algorithms,” Computational Complexity, Vol. 2, pp. 111 – 128, 1992.
\\ (https://www.researchgate.net/publication/2802927_Improved_Low-Density_Subset_Sum_Algorithms.)
job() = {
  i = 0;
  while(1,
    \\ Cas où m[1] = 0
    m[l, l] = - B * (c + i * N);
    lll = m * qflll(m); \\ Base réduite du réseau
    for(j = 1, l,
      sol = shiftvec(lll[,j])[1..-2]; \\ On ne garde pas le dernier élément qui doit valoir 0 pour être solution
      s = concat([0],sol);
      if(check_col(sol) && check_eq(s), printmessage(s); return)); \\ Si la colonne est de la bonne forme et vérifie l'égalité sac à dos

    \\ Cas où m[1] = 1
    m[l, l] = - B * (c - k1 + i * N); \\ On retranche le poids b[1]
    lll = m * qflll(m);
    for(j = 1, l,
      sol = shiftvec(lll[,j])[1..-2];
      s = concat([1],sol);
      if(check_col(sol) && check_eq(s), printmessage(s); return));
    i++;
  );
}

job();

\\ génération de tests
/*
generate_superincreasing_sequence(N) = {
  my(s, seq, r);
  s = 0;
  seq = List();
  for(i = 1, N, r = random(200); s += (s+r); print("r=", r, " ; s = ", s); listput(seq, s));
  seq;
}

privk = Vec(generate_superincreasing_sequence(len));
N = nextprime(23000);
trapdoor = nextprime(234);
k = [(trapdoor*x)%N|x<-privk];
msg = vector(len, n, random(2));
c = msg*k~;
cmod = c % N;
\\ la solution est trouvée pour i = c\N, pourvu que la densité soit faible et avec un peu de chance
*/
