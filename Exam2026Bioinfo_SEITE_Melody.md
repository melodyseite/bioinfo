# Bioinfo exam 2026

## Exercise 1: Take inventory of the dataset

**Context:** in population genetics, if one individual has 10 times more reads than another, comparisons between them are biased. Before anything else you must know how many reads each individual has.

### 1. Look at the first record of `popA_ind01.fastq` with `head`. Identify the four lines.

```bash
e22206540@ed58:~/HOME/data$ head popA_ind01.fastq
```

```text
@SRRpopA_ind01.1 GRC13_0027_FC:4:1:1:7 length=50
TGCAGGCATTGCGACGGGCCATTAATGTCGACCCTCTTGCTGTATTGAAG
+SRRpopA_ind01.1 GRC13_0027_FC:4:1:1:7 length=50
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII@IIIIIIIIIIIIIII
@SRRpopA_ind01.2 GRC13_0027_FC:4:1:2:14 length=50
TGCATACTATTGAAAGATGAGGGGCAACCAGAATGATCTGCGTTAGCTAT
+SRRpopA_ind01.2 GRC13_0027_FC:4:1:2:14 length=50
IIIIIIIIIIIIIIII@IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRRpopA_ind01.3 GRC13_0027_FC:4:1:3:21 length=50
TGCAGACAGTGTTATTACGCCCGTCAGGAACCGCTGCGAATCATCTGGGC
```

Les 4 premières lignes du fichier sont l'identifiant de la séquence, avec des informations comme la longueur, la deuxième ligne correspond à la séquence de nucléotides, puis la troisième correspond au séparateur, et enfin la dernière correspond à la qualité de la lettre écrite par des symboles et qui pourra être lue à l'aide d'une table.

### 2. Count the lines of `popA_ind01.fastq`, then compute the number of reads using `$(( ... ))`. Store it in a variable `n` and display it with `echo`.

```bash
e22206540@ed58:~/HOME/data$ n=$((` wc -l < popA_ind01.fastq` / 4 ))
e22206540@ed58:~/HOME/data$ echo $n
300
```

### 3. Write a script `ex1_count_reads.sh` that takes any number of files as arguments and prints one line per file: `filename population number_of_reads`, separated by tabs.

```bash
e22206540@ed58:~/HOME/data$ touch ex1_count_reads.sh
e22206540@ed58:~/HOME/data$ nano ex1_count_reads.sh
```

Dans la page qui s'ouvre avec nano :

```bash
#!/bin/bash
for f in "$@"; do
n=$(( $(wc -l < "$f") / 4))
pop=$(basename "$f" | cut -d_ -f1)
printf "%s\t%s\t%s\n" "$f" "$pop" "$n"
done
```

`Ctrl+O` puis `Entrée` pour sauvegarder, puis `Ctrl+X` pour quitter nano.

Dans le terminal :

```bash
e22206540@ed58:~/HOME/data$ chmod a+x ex1_count_reads.sh
e22206540@ed58:~/HOME/data$ ./ex1_count_reads.sh popA_ind01.fastq
popA_ind01.fastq	popA	300
```

### 4. Run it on the files of population A only (`popA_*.fastq`), then on all files.

Pour popA :

```bash
e22206540@ed58:~/HOME/data$ ./ex1_count_reads.sh popA_*.fastq
popA_ind01.fastq	popA	300
popA_ind02.fastq	popA	250
popA_ind03.fastq	popA	200
```

Pour tous les fichiers :

```bash
e22206540@ed58:~/HOME/data$ ./ex1_count_reads.sh *.fastq
popA_ind01.fastq	popA	300
popA_ind02.fastq	popA	250
popA_ind03.fastq	popA	200
popB_ind01.fastq	popB	400
popB_ind02.fastq	popB	275
popB_ind03.fastq	popB	225
```

### 5. How many reads in total for population A? For population B? (Pipe `cat` into `wc -l`.)

```bash
e22206540@ed58:~/HOME/data$ echo $(( $(cat popA_*.fastq | wc -l) / 4 ))
750
e22206540@ed58:~/HOME/data$ echo $(( $(cat popB_*.fastq | wc -l) / 4 ))
900
```

### 6. Trap. Run `grep -c "^@" popA_ind01.fastq` and compare with your count. Why is it different? Try `grep -c "^@SRR" popA_ind01.fastq`. What do you conclude about counting reads with `grep`?

```bash
e22206540@ed58:~/HOME/data$ grep -c "^@" popA_ind01.fastq
307
```

Cela ne donne pas le même résultat et c'est possible parce que `@` peut aussi être un caractère de qualité, donc cette ligne compte aussi certaines lignes de qualité.

```bash
e22206540@ed58:~/HOME/data$ grep -c "^@SRR" popA_ind01.fastq
300
```

Avec cette ligne on ne garde que les lignes dont le header commence par `@SRR` (peu probable qu'une ligne de qualité commence par cela).

**Conclusion :** `grep` sur `^@` n'est pas fiable, `wc -l / 4` est plus sûr, ou à la limite un motif spécifique au header.

---

## Exercise 2: Which individuals look unusable?

**Context:** sequencers write `N` when a base could not be called. Individuals with many `N` are suspicious (poor DNA, failed library).

### 1. For `popA_ind01.fastq`, extract the sequence lines (awk one-liner above) and count how many contain no `N` (`grep -v`).

`awk 'NR%4==2'` garde la ligne 2 de chaque bloc de 4 (la séquence), `grep -v` exclut celles contenant `N` et `-c` compte le reste.

```bash
e22206540@ed58:~/HOME/data$ awk 'NR%4==2' popA_ind01.fastq | grep -vc "N"
278
```

### 2. Compute the percentage of reads without any `N` (integer arithmetic is fine).

Dans le script `exo2.3.sh` :

```bash
#!/bin/bash
f=$1
lignes=$(wc -l < "$f")
total=$((lignes / 4))
clean=$(awk 'NR%4==2' "$f" | grep -vc "N")
echo $((clean * 100 / total))
```

Dans le terminal :

```bash
e22206540@ed58:~/HOME/data$ ./exo2.3.sh popA_ind01.fastq
92
```

### 3. Write `ex2_quality_N.sh`: for every file given as argument, print `file n_reads n_reads_without_N percent`. Print a header line first.

Dans le script `ex2_quality_N.sh` :

```bash
#!/bin/bash
printf "file\tn_reads\tn_reads_without_N\tpercent\n"
for f in "$@"; do
lignes=$(wc -l < "$f")
total=$((lignes / 4))
clean=$(awk 'NR%4==2' "$f" | grep -vc "N")
pct=$((clean * 100 /total))
printf "%s\t%s\t%s\t%s\n" "$f" "$total" "$clean" "$pct"
done
```

### 4. Run it on all individuals. Which one would you remove from a population genetic analysis?

Dans le terminal :

```bash
e22206540@ed58:~/HOME/data$ chmod a+x  ex2_quality_N.sh
e22206540@ed58:~/HOME/data$ ./ex2_quality_N.sh *.fastq
```

```text
file	n_reads	n_reads_without_N	percent
popA_ind01.fastq	300	278	92
popA_ind02.fastq	250	217	86
popA_ind03.fastq	200	181	90
popB_ind01.fastq	400	364	91
popB_ind02.fastq	275	252	91
popB_ind03.fastq	225	58	25
```

Il serait intéressant d'enlever le dernier individu (`popB_ind03.fastq`) qui possède un nombre de reads sans `N` très faible (58), c'est aussi visible dans le % (25 %). L'ADN est sans doute de mauvaise qualité ou la librairie est ratée, le génotype est peu fiable.

---

## Exercise 3: Check the restriction site

**Context:** in RAD-seq, DNA is cut by a restriction enzyme, so nearly all reads should start with the same short motif. In the *Pedicularis* data of the course this motif is `TGCA`. Reads that do not start with it come from something else (contamination, adapters, a failed digest).

### 1. Count the reads of `popA_ind01.fastq` whose sequence starts with `TGCA`. (`^` anchors a pattern at the start of a line.)

```bash
e22206540@ed58:~/HOME/data$ awk 'NR%4==2' popA_ind01.fastq | grep -c "^TGCA"
284
```

### 2. Write `ex3_cutsite.sh` taking the motif as first argument and the files as the following arguments, for example: `./ex3_cutsite.sh TGCA popA_*.fastq` (hint: `motif=$1` then `shift`).

Dans le script `ex3_cutsite.sh` :

```bash
#!/bin/bash
motif=$1
shift
printf "file\tn_reads\tn_reads_motif\tpercent\n"
for f in "$@" ; do
lignes=$(wc -l < "$f")
total=$((lignes / 4))
ok=$(awk 'NR%4==2' "$f" |grep -c "^$motif")
pct=$((ok * 100 /total))
printf "%s\t%s\t%s\t%s\n" "$f" "$total" "$ok" "$pct"
done
```

Dans le terminal :

```bash
e22206540@ed58:~/HOME/data$ ./ex3_cutsite.sh TGCA *.fastq
```

```text
file	n_reads	n_reads_motif	percent
popA_ind01.fastq	300	284	94
popA_ind02.fastq	250	233	93
popA_ind03.fastq	200	194	97
popB_ind01.fastq	400	382	95
popB_ind02.fastq	275	260	94
popB_ind03.fastq	225	52	23
```

### 3. Output: `file motif n_reads_with_motif percent`.

Dans le script :

```bash
#!/bin/bash
motif=$1
shift
printf "file\tn_reads\tn_reads_motif\tpercent\n"
for f in "$@" ; do
lignes=$(wc -l < "$f")
total=$((lignes / 4))
ok=$(awk 'NR%4==2' "$f" |grep -c "^$motif")
pct=$((ok * 100 /total))
printf "%s\t%s\t%s\t%s\n" "$f" "$motif" "$ok" "$pct"
done
```

Dans le terminal :

```bash
e22206540@ed58:~/HOME/data$ ./ex3_cutsite.sh TGCA *.fastq
```

```text
file	n_reads	n_reads_motif	percent
popA_ind01.fastq	TGCA	284	94
popA_ind02.fastq	TGCA	233	93
popA_ind03.fastq	TGCA	194	97
popB_ind01.fastq	TGCA	382	95
popB_ind02.fastq	TGCA	260	94
popB_ind03.fastq	TGCA	52	23
```

### 4. Which individual behaves differently? Give two possible explanations. Is it the same individual as in exercise 2?

L'individu `popB_ind03.fastq` a un très faible nombre de reads avec motif et donc un très faible % (23 %). C'est le même individu que dans l'exercice 2.

Les explications possibles sont :

- une digestion ratée par l'enzyme, qui ne coupe pas au bon endroit, donc la plupart des reads ne commencent pas au site de restriction ;
- ou une contamination d'ADN d'une autre espèce.

---

## Exercise 4: A small test dataset, and archiving

**Context:** you never test a new pipeline on the full dataset. You build a small subsample first, and you archive and compress what you send to collaborators.

### 1. Create a folder `subsample`.

```bash
e22206540@ed58:~/HOME/data$ mkdir subsample
```

### 2. Write `ex4_subsample.sh` taking the number of reads to keep as `$1` and the files as the other arguments. It writes, in `subsample/`, the first n reads of each file. Careful: `head` counts lines, and a read is 4 lines. Compute the number of lines in the script.

Dans le script :

```bash
#!/bin/bash
n=$1
shift
lignes=$((n * 4))
mkdir -p subsample
for f in "$@" ; do
head -n "$lignes" "$f" > "subsample/$(basename "$f")"
done
```

Dans le terminal :

```bash
e22206540@ed58:~/HOME/data$ ./ex4_subsample.sh 100 *.fastq
e22206540@ed58:~/HOME/data/subsample$ wc -l *.fastq
```

```text
 400 popA_ind01.fastq
 400 popA_ind02.fastq
 400 popA_ind03.fastq
 400 popB_ind01.fastq
 400 popB_ind02.fastq
 400 popB_ind03.fastq
2400 total
```

### 3. Check your result: `wc -l subsample/*`. Does it match what you asked for?

```bash
e22206540@ed58:~/HOME/data$ wc -l subsample/*
```

```text
400 subsample/popA_ind01.fastq
400 subsample/popA_ind02.fastq
400 subsample/popA_ind03.fastq
400 subsample/popB_ind01.fastq
400 subsample/popB_ind02.fastq
400 subsample/popB_ind03.fastq
2400 total
```

Oui, cela correspond bien.

### 4. Make a compressed archive `subsample.tar.gz` of the folder with `tar`.

```bash
e22206540@ed58:~/HOME/data$ tar -czvf subsample.tar.gz subsample/
e22206540@ed58:~/HOME/data$ ls -lh subsample.tar.gz
```
