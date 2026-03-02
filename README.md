# Ilea_Watson_Advanced_Python_Portfolio
This is a portfolio of the python code that I completed during BISC 4503

# Sequence Objects Pt. 1
### Load package

```python
from Bio.Seq import Seq
```


```python
# Make a sequence object
my_seq = Seq("GATCG")
```


```python
# Number each nucleotide on the sequence
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
#We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
# We can search for objects by their locations in the sequence
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
# We can also do a .count function to count the occurances
Seq("AAAA").count("AA")
```




    2




```python
# The len() and .count functions can be used to messure the GC content
my_seq = Seq("GATCGATGGCATACGTCA")
```


```python
len(my_seq)
```




    18




```python
my_seq.count("G")
```




    5




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    50.0




```python
# Calculating GC content is built into biopython
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGATGGCATACGTCA")
```


```python
gc_fraction(my_seq) * 100
```




    50.0




```python
# Sequences can be sliced into multiple parts
my_seq[0::3]
```




    Seq('GCTCAT')




```python
my_seq[1::3]
```




    Seq('AGGACC')




```python
my_seq[2:3]
```




    Seq('T')




```python
# Negatives can be used to start from the other end
# This prints the sequence backwards
my_seq[::-1]
```




    Seq('ACTGCATACGGTAGCTAG')




```python
# Seq can be saved as strings
str(my_seq)
```




    'GATCGATGGCATACGTCA'




```python
# You can label or name the string
# This is labeling the string in a fasta format
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATCGATGGCATACGTCA
    



```python
# Adding sequences together
# Sequences will be added based on order
seq1 = Seq("ACTG")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACTGAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACTG')




```python
# Sequence Objects Pt. 2
# Spacer can be used to join contigs
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
# For case sensitive sequences the .upper() & .lower() functions can be used
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq = dna_seq.upper()
```


```python
# You can search for specific nucleotides in a sequence
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GTTACGCTTCAGGCAAATCGT")
```


```python
# the complement() function gives the complementary strand sequence
# the reverse_complement() function gives the complementary strand starting at the opposite end
my_seq.complement()
```




    Seq('CAATGCGAAGTCCGTTTAGCA')




```python
my_seq.reverse_complement()
```




    Seq('ACGATTTGCCTGAAGCGTAAC')




```python
#This can also be done for protein sequences
# Pyton will do it, but it does not make sense
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
# Lets create a DNA sequence
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# Now lets make the reverse complement of our DNA sequence, like transcription
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
# This creates the transcription of the coding_dna
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# This mimics biological transcription, but should give the same result
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# The .back_transcribe() fuction can be used to get the original DNA sequence
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# The .translate function translates the rna sequence to amino acids
# "*" indicate stop codons
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# Sequence objects Pt. 3
# You can specify which codon table to use
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
#This can also be found using the ncbi table number
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
# You can have the translation stop at the stop codon
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table =2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
# You can change the symbol that denotes the stop codon
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# This tells python the sequnce is a complete DNA sequence and changes the start codon from Val to Met
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# Lets look at codon tables
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# You can search the stable for specific codons
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# Sequence comparisons
# The name of the seq is equivalent to the sequence
seq = Seq("ACGT")
```


```python
"ACGT" == seq
```




    True




```python
seq == "ACGT"
```




    True




```python
# You can create a sequence with an unknown sequence of a specific length
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
# Sequence objects Pt. 4
# We can create a sequence with some known sequencing of a specific length
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# This will just show the length since the sequence is unknown
seq[1000:1020]
```




    Seq(None, length=20)




```python
# This will show the known sequence
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# This will give you the known sequnce and the length to the end
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
# You can "mutate" your sequence
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# .remove will remove the first instance of the sequence
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# This will flip the sequence
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# Make a new sequence makes it no longer mutable
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# This can be used instead of seq objects
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'





