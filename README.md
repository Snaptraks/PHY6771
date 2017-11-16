# PHY6771
Projet commun de code de modèle d'atmosphère pour le cours PHY6771

Hello world!

Compiler le code et l'exécuter
Dans .cshrc:
```
alias ismx 'ifort -o \!:1:r.x \!* -W0 -O1 -warn none -save -lsm -lX11 -lm && \!:1:r.x'
```

Dans le terminal:
```
ismx <nomfichier>.f
```
ou 
```
ifort -o <nomfichier>.x <nomfichier>.f -W0 -O1 -warn none -save -lsm -lX11 -lm && <nomfichier>.x
```
