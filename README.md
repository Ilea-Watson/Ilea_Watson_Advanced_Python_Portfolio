# Ilea_Watson_Advanced_Python_Portfolio
This is a portfolio of the python code that I completed during BISC 4503 in Winter Quarter 2025-2026.

# Sequence Objects
### Part 1


```python
# Load package
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


### Part 2

```python
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
# This can also be done for protein sequences
# Pyton will do it, but it does not make sense
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
# Let's create a DNA sequence
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


### Part 3

```python
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


### Part 4

```python
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

# Sequence Annotations
### Part 1

```python
# In terminal pip install biopython
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
# help(SeqRecord) to pull up the help file
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
# Create a seq record to store sequences with annotaions
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
# You can search for specific information or pull up the whole seq record
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
# Annotations can be added that do not show when you pull up the record
simple_seq_r.annotations["evidance"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidance"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
# Save fasta file as "NC_005816.fna" and upload fasta file to JupyterLab
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# You can only pull up information saved in the record
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []


### Part 2

```python
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
# Save file as "NC_005816.gb" and upload fasta file to JupyterLab
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
# You can mark the postition of something in the sequence
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
# You can show the specific start and end locations
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
# You can show just the number values
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
# You can put the locations together
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)


### Part 3

```python
from Bio.SeqRecord import SeqRecord
```


```python
# We're creating the record and putting it in fasta file format
record = SeqRecord(Seq("MMVQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
                      "GAGVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
                      "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILQVEAKLGLKEEKMRATREVLSEYGNM"
                      "SSAC"),
                              id="gi|14150838|gb|AAK54648.1|AF376133_1",
                  description="chalcone synthase [Cucmis sativus]",
        )
```


```python
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucmis sativus]
    MMVQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKN
    IEKSLKEAFTPLGISDWNSTFWIAHPGGPAILQVEAKLGLKEEKMRATREVLSEYGNMSS
    AC
    



```python
# The record doesn't have to be printed in a fasta file format
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucmis sativus]
    Number of features: 0
    Seq('MMVQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
# You can search for different things in the record using the following functions
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# Tells the length
len(record)
```




    9609




```python
# The number of features
len(record.features)
```




    41




```python
# Display a specific feature
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# We can subdivide the record
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
# The number of records will be perserved in the subdivison
len(sub_record.features)
```




    2




```python
# Python starts counting at 0, so there are 2 features in our subdivision, feature 0 and 1
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
# All information will remain the same except the location because it is looking at the location within the context of the subrecord
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# We can add to annotations of the record
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
# Gives the file id
sub_record.id
```




    'NC_005816.1'




```python
# In this case the id is similar to the name
sub_record.name
```




    'NC_005816'




```python
# We can edit the description to match the edits we've made to the record in the subdivision
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



### Part 4


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id ="Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0

# Sequence Input/Output
### Part 1


```python
# In terminal pip install biopython
from Bio import SeqIO
```


```python
# We can manipulate these long database files with multiple organisms for a single genes using the parse function
# This example uses the fasta file
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# We can do a similar thing for our other file type
# This example uses the genbank file
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# Let's extract things as a list
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
# We can create a list and extract things from that list
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
#Next is used to extract the next iteration in the sequence
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
#If we want to access a specific record in the sequence
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740

### Part 2

```python
# Review from pt.1
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
#This is giving us a readable version of the first record
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
# We can pull up the annotations
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
# Or we can just pull up the keys from the annoations
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
#This gives the annotation values without the keys
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
# We can pull out the information we need
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
#Source will typically give the common name and organism will give the scientific name
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
# Let's create an empty vector and parse all the species in the file
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# This is a different way of doing the same thing
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank" )
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
# Lets do the same thing for the fasta file
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
# It gives different results because the annotations are different in the fasta file
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
# Let's alter some data
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
# Let's rename this record
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
#This does not change the data in the original file
```


```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
# Let's add a description to this record
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
print(first_record.format("fasta"[:200]))
```

    >new_id Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    

### Part 3

```python
# This is how you write a file
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
    Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
        "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
        "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
        "SSAC",
),
    id="gi|14150838|gb|AAK54658.1|AF376133_1",
    description = "chalcone synthase [Cucumis sativus]",
)
```


```python
rec2 = SeqRecord(
Seq("YPDYYFRINREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
       "DMVVVEIPKLGKEAAVKAIKEWGQ",
   ),
    id = "gi | 13919613 |gb| AAK33142.1|",
    description = "chalcone synthase [Fragaria vesc subsp. bracteata]",
)
```


```python
rec3 = SeqRecord(
Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
       "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
       "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
       "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
       "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
       "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
       "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
   ),
    id="gi|13925890|gb|AAK49457.1|",
    description="chalcone synthase [Nicotiana tabacum]",
)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa", "fasta")
```




    3




```python
# Let's convert to a fasta file
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
# Let's do some reverse complements 
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
# This code annotates the reverse complements we just made
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
len(records)
```




    94




```python
# We can select for records based on the length and annotate those records
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
len(records)
```




    18




```python
# Let's put this all together and create a fasta file
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18
 
# Multiple Sequence Alignment
### Part 1

```python
# pip install biopython
```


```python
from Bio import AlignIO
```


```python
# Download https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05371_seed.sth
# Let's generate our alignment from the file
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Alignmnet length %i" % alignment.get_alignment_length())
```

    Alignmnet length 52



```python
# This is what AlignIO.read does
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# We can call python's built in function on the alignment object to show it's of a particular file format
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
# We can write an alignment file
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
align1 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
        SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
        SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
    ]
)

align2 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
        SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
        SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
    ]
)

align3 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
        SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
        SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
    ]
)
```


```python
my_alignment = [align1, align2, align3]
```


```python
print(my_alignment)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f096056b710>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f09701ddd50>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f0972b3be50>]



```python
from Bio import AlignIO
AlignIO.write(my_alignment, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
# The file has to be converted to a list to pull out individual alignments
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
# There are 3 alignments in the file, so this will print the last alignment
last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma

### Part 2

```python
from Bio import AlignIO
```


```python
# Convert the alignment file to a clustal file
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignmnents" % count)
```

    Converted 1 alignmnents



```python
# This is another way to do that
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignmnets" % count)
```

    Converted 1 alignmnets



```python
# We can also convert this way
# The sample names are truncated in this file
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
# This creates a file with the entire sample name
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
# We can number all the elements in the file
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
# We can apply this to our other files
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
# Let's do some slicing alignments
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# You can pull out specific segments
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
# This pulls out what is located in row 2 column 6
print(alignment[2, 6])
```

    T



```python
#This does the same thing
print(alignment[2].seq[6])
```

    T



```python
# This prints the content of column 6 in each row
print(alignment[:, 6])
```

    TTT---T



```python
# This selects a range of the first 3 rows and the first 6 characters of wach row
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
# This takes the first 6 characters of all of the rows
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
# This takes column 6-9 of all of the rows
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
# This shows all the rows and everthing after the 9th column
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# You can remove blocks of columns
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# This sorts based on the id in the file
edited.sort()
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49

### Part 3

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id="seq1"),
        SeqRecord(Seq("AAT-CTA"), id="seq2"),
        SeqRecord(Seq("CCTACT-"), id="seq3"),
        SeqRecord(Seq("TCTCCTC"), id="seq4")
    ])
```


```python
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
# This calculates all pairs of rows in the alignment and counts how many times they are paired together
substituations = alignment.substitutions
```


```python
print(substituations)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
# We can add another substituation matrix
m = substituations.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# Or you can change the order
m = substituations.select("ACTG")
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
# You can do alignments in the commandline or python
```


```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
# We can make a tree using a newick file
from Bio import Phylo
```


```python
tree = Phylo.read("opuntia.dnd", "newick")
```


```python
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    

### Part 4

```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
# We can change the match score
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
# Let's create out target and quey. The query is compared to the target.
target = "GAACT"
```


```python
query = "GAT"
```


```python
# This gives us the number of overlaps of the query to the target sequence
score = aligner.score(target, query)
```


```python
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
# Let's do another example
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    5.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
# This lets you know that algorithm used and the significance of the alignment
aligner.algorithm
```




    'Smith-Waterman'




```python
aligner.epsilon
```




    1e-06




```python
# Lets make an alignment object
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
alignment = alignments[0]
```


```python
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f0960485ad0>




```python
# By making it an object we can pull out desired information
alignment.score
```




    3.0




```python
alignment.target
```




    'GAACT'




```python
alignment.query
```




    'GAT'




```python
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
len(alignment)
```




    2




```python
alignment.shape
```




    (2, 5)




```python
# This changes the mode of the alignment
aligner.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.shape
```




    (2, 4)




```python
aligner.mode = "global"
```


```python
aligner = Align.PairwiseAligner(match_score = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    

### Part 5

```python
from Bio.Seq import reverse_complement
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
# Aligned to positive strand
aligner.score(target, query)
```




    4.0




```python
# Aligned to reverse complement strand
aligner.score(target, reverse_complement(query))
```




    0.0




```python
# Aligned to the reverse complement of the negative strand (which is the reverse complement) (like a double negative)
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
# Aligned to the negative strand (which is the reverse complement)
aligner.score(target, query, strand = "-")
```




    0.0




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# We can save this alignment to a different file format
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
# This will give the same result
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
# This alignment will not print because there is no alignment
alignments = aligner.align(target, query, strand = "-")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# This gives us the same score as before
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
# Mapping to the negative strand gives a worse score
aligner.score(target, reverse_complement(query), strand = "+")
```




    -2.0

# Blast

```python
# pip install biopython
```


```python
from Bio.Blast import NCBIWWW
```


```python
# Having an email connected to NCBI can let you know when things are done running
NCBIWWW.email = "irw003@latech.edu"
```


```python
# You can Blast using the NCBI number
# This opens the connection to the NCBI search
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
# Or you can download the fasta file and run it from your local files
#https://raw.githubusercontent.com/biopython/biopython/refs/heads/master/Doc/examples/m_cold.fasta
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta", format = "fasta")
```


```python
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
# This writes and reads the results from the NCBI search, then closes the connection
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
# We're going to save it as a XML file
from Bio.Blast import NCBIXML
```


```python
# Then reopen the connection to the NCBI search
result_handle = open("my_blast.xml")
```


```python
# This lets us read in the saved result handle
blast_record = NCBIXML.read(result_handle)
```


```python
# We can select for the e value
E_VALUE_THRESH = 0.04
```


```python
# This makes the results of the Blast nicer to look at
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.length)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: 1173
    length: 1173
    e value: 6.88175e-117
    ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTC...
    || ||||||||| |||| | |||| ||  |||| |||| | |||| ||| | |||| ||| ||| ||||| | ||...
    ACCGAAAATGGGCAGAGGAGTGAATTATATGGCAATGACACCTGAGCAACTAGCCGCGGCCAATTTGATCAACTC...
    ****ALIGNMENT****
    sequence: 752
    length: 752
    e value: 1.84663e-111
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||| |||  |||| | || ||||| |||||||| || ||||| |||| ||| ||| ||||||||||||||...
    AAAATGGGTAGACGAATGGATTATTTGGCGATGAAAACCGAGCAATTAGCCGCGGCCAATTTGATCGATTCCGAT...
    ****ALIGNMENT****
    sequence: 621
    length: 621
    e value: 4.9552e-106
    TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACA...
    ||||||||||||||||| ||| ||||  |||||||| |||| ||||  ||||| ||||| ||||| || ||    ...
    TTGGCCATGAAAACTGAGCAAATGGCGTTGGCTAATTTGATAGATTATGATATGAATGAACTTAAGATCGCTTTG...
    ****ALIGNMENT****
    sequence: 1020
    length: 1020
    e value: 1.72953e-105
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    ||||||||| |||  | |  | |||||||||||||||||||    ||||  |||  || ||||||| || |||| ...
    AATGGGGAG-GAA--GGATAATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATAA...
    ****ALIGNMENT****
    sequence: 1028
    length: 1028
    e value: 6.03667e-105
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  ||||| |||| |||||||| |   |||  |||| |  ||||  |||| |||...
    AAAATGGGGAGG---ATGGAGTTTTTGGCTATGAGAACTGATCCA---GCCACGGCTGACTTGATAAATTCTGAT...
    ****ALIGNMENT****
    sequence: 1020
    length: 1020
    e value: 7.35417e-104
    ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
    |||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
    ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
    ****ALIGNMENT****
    sequence: 1036
    length: 1036
    e value: 7.35417e-104
    ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
    |||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
    ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
    ****ALIGNMENT****
    sequence: 909
    length: 909
    e value: 8.95921e-103
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
    ||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
    AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
    ****ALIGNMENT****
    sequence: 1025
    length: 1025
    e value: 8.95921e-103
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
    ||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
    AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
    ****ALIGNMENT****
    sequence: 1052
    length: 1052
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 868
    length: 868
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 952
    length: 952
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 881
    length: 881
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 980
    length: 980
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 913
    length: 913
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 978
    length: 978
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 893
    length: 893
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 968
    length: 968
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: 905
    length: 905
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || |||||  ||| ||||...
    AAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: 1009
    length: 1009
    e value: 2.92878e-96
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    |||||||||||||    ||| | || ||||| ||||| ||||||||   ||||   ||| || | |||  ||| |...
    AGAAAATGGGGAGG---ATGGAATATTTGGCTATGAAGACTGATCA---GGCCACTGCTGATCTCATCACTTCTG...
    ****ALIGNMENT****
    sequence: 1071
    length: 1071
    e value: 1.24535e-94
    TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
    || |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
    TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
    ****ALIGNMENT****
    sequence: 1065
    length: 1065
    e value: 1.24535e-94
    TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
    || |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
    TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
    ****ALIGNMENT****
    sequence: 1069
    length: 1069
    e value: 4.3467e-94
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: 1072
    length: 1072
    e value: 4.3467e-94
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: 1091
    length: 1091
    e value: 4.3467e-94
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: 824
    length: 824
    e value: 1.51715e-93
    ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
    ||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
    ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
    ****ALIGNMENT****
    sequence: 844
    length: 844
    e value: 1.51715e-93
    ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
    ||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
    ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
    ****ALIGNMENT****
    sequence: 999
    length: 999
    e value: 5.29536e-93
    GAACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGAT...
    ||| ||||||||||||||   | | | || ||||| ||||| |||||||    ||  ||||   || |||||   ...
    GAAAAGAAAATGGGGAGA---ACGGATTATTTGGCGATGAAGACTGATCC---GGTTGTGGACGATTTGATCAGC...
    ****ALIGNMENT****
    sequence: 1015
    length: 1015
    e value: 6.45107e-92
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    |||||||||    ||||| || ||||| |||||||||||||    |||  ||| | || ||||||  || |||||...
    AATGGGGAGG---ATGAATTATTTGGCTATGAAAACTGATCC---GGCAATGGATGATTTGATCGGCTCTGATAT...
    ****ALIGNMENT****
    sequence: 1152
    length: 1152
    e value: 6.45107e-92
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    ||||||||||||   ||| | || ||||| ||||| || |||| |  |||  || | ||| ||||| |||| || ...
    AAAATGGGGAGA---ATGGATTATTTGGCGATGAAGACCGATCCAGCGGC--TGCCGAAT-TGATCAATTCGGAC...
    ****ALIGNMENT****
    sequence: 837
    length: 837
    e value: 2.25165e-91
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    |||||||||    ||| | | ||||||||||||||| ||||    |||||  |   || ||||| ||||||| ||...
    AATGGGGAGG---ATGGATTTCTTGGCCATGAAAACAGATCCCGCGGCCGCCG---ATTTGATCAATTCCGACAT...
    ****ALIGNMENT****
    sequence: 1110
    length: 1110
    e value: 2.74307e-90
    GCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGT...
    ||||  ||| | |||||||| || |||||||| ||||| || || ||  | | | | || || |   | || |  ...
    GCCGATGCTGAAATGATCGACTCGGATATCAACGAGCTGAAGATCGCGGCCAAGCGACTGATTAGCCACGCCACC...
    ****ALIGNMENT****
    sequence: 850
    length: 850
    e value: 9.57424e-90
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    |||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
    TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
    ****ALIGNMENT****
    sequence: 2364
    length: 2364
    e value: 1.16638e-88
    GAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
    ||||||||||||      | |||||||||||||||||||||||| |  ||| |  |   || |||||  ||| ||...
    GAAAATGGGGAGGTC---GGAGTACTTGGCCATGAAAACTGATCCAGCGGCTGCCG---ATTTGATCAGTTCGGA...
    ****ALIGNMENT****
    sequence: 787
    length: 787
    e value: 1.16638e-88
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    |||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
    TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
    ****ALIGNMENT****
    sequence: 1020
    length: 1020
    e value: 4.07107e-88
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||| ||       || |||||||||||||| || | |   |||  |||| || ||||  |||| ||||...
    AAATGGGGAGGGAT------TATTTGGCCATGAAAACGGACCCA---GCCACGGCTGATTTGATAAATTCTGATA...
    ****ALIGNMENT****
    sequence: 1132
    length: 1132
    e value: 4.07107e-88
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
    ****ALIGNMENT****
    sequence: 980
    length: 980
    e value: 4.07107e-88
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
    ****ALIGNMENT****
    sequence: 982
    length: 982
    e value: 4.95958e-87
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT-GGCCGTGGCTAATATGATCGATTCCGA...
    |||||||||||    ||| | ||  |||  ||||||||||| ||||  ||     ||| |  ||||  | || ||...
    AAAATGGGGAGG---ATGGATTATCTGGGAATGAAAACTGA-CAATCAGGTTACTGCTGAGGTGATTAACTCTGA...
    ****ALIGNMENT****
    sequence: 1018
    length: 1018
    e value: 1.73106e-86
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: 1178
    length: 1178
    e value: 1.73106e-86
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: 910
    length: 910
    e value: 6.042e-86
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    |||||||||||||  ||||| |||| ||||| ||||| || |||||    |||    |  |   |||  ||||||...
    AGAAAATGGGGAGGAAAATGGAGTATTTGGCAATGAAGACCGATCATCCAGCCCCAACCCAATCGATGAATTCCG...
    ****ALIGNMENT****
    sequence: 945
    length: 945
    e value: 2.10887e-85
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    ||||||||||    ||| ||||||||   ||||  ||||||||| |  |   |||  || ||||| | || ||| ...
    AAATGGGGAGG---ATGGAGTACTTGAAAATGAGTACTGATCAAGTACC---GGCCGATTTGATCAAGTCTGATC...
    ****ALIGNMENT****
    sequence: 937
    length: 937
    e value: 2.10887e-85
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   | |||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGGTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: 1175
    length: 1175
    e value: 2.10887e-85
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: 811
    length: 811
    e value: 2.56913e-84
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
    ****ALIGNMENT****
    sequence: 929
    length: 929
    e value: 2.56913e-84
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
    ****ALIGNMENT****
    sequence: 1126
    length: 1126
    e value: 2.56913e-84
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    |||||||||||||| |   | |||| ||||| |||||||||||||    ||||   ||  |  ||||| ||||||...
    AGAAAATGGGGAGAAA---GGAGTATTTGGCTATGAAAACTGATCC---GGCCACAGCCGAATTGATCAATTCCG...
    ****ALIGNMENT****
    sequence: 850
    length: 850
    e value: 8.96713e-84
    ATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATTACGGGT...
    |||| |||||||| ||||||||||| ||| | | ||||||  |  |  |||| |  | |||||||    | ||  ...
    ATTCTGATATCAACGAGCTTAAAATTGCAGCCACGAGGCTTCTTGAACATGCCACCAAGCTCGGTGGAAAGGGCC...
    ****ALIGNMENT****
    sequence: 1083
    length: 1083
    e value: 3.12984e-83
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    || ||||||||    ||| | ||  |||| |||||||||||  |   ||  |  ||| |  ||||| | || |||...
    AAGATGGGGAGG---ATGGATTATCTGGCAATGAAAACTGACGATCAGGTTGCTGCTGACTTGATCAACTCTGAT...


# Challenge

```python
# pip install Biopython
# For this challenge I chose the gene for telomerase reverse transcriptase (TERT)
# TERT contributes to cancer formation by extending telomeres shortened by rapid cell divison
```


```python
# Load in the needed packages and establish email
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "irw003@latech.edu"
```


```python
from Bio import SeqIO
```


```python
# I downloaded this fasta file for the human TERT gene from NCBI
# https://www.ncbi.nlm.nih.gov/nuccore/NM_198253.3?report=fasta
record = SeqIO.read("TERT.fasta", format = "fasta")
```


```python
# Confirm the file is correct
print(record)
```

    ID: NM_198253.3
    Name: NM_198253.3
    Description: NM_198253.3 Homo sapiens telomerase reverse transcriptase (TERT), transcript variant 1, mRNA
    Number of features: 0
    Seq('CTCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGA...AAA')



```python
# Run blast using the sequence from the created record
# For the purpose of this challenge I am comparing the human TERT gene to the Chimpanzee TERT gene
# entrez_query used to select for Chimp TERT gene sequence, remove clone and unverified sequences, and restrict search to RefSeq's
result_handle = NCBIWWW.qblast(
    "blastn",
    "nt",
    record.seq,
    entrez_query = "Pan troglodytes[Organism] NOT clone[Title] NOT unverified[Title] AND RefSeq[Filter]"
)
```


```python
# This writes and reads the results from the NCBI search, then closes the connection
with open("my_blast2.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
# Then reopen the connection to the NCBI search
result_handle = open("my_blast2.xml")
```


```python
# This lets us read in the saved result handle
blast_record = NCBIXML.read(result_handle)
```


```python
# We can select for the e value
E_VALUE_THRESH = 0.04
```


```python
# This makes the results of the Blast nicer to look at

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("Sequence:", alignment.title)
            print("length:", alignment.length)
            print("E Value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    Sequence: gi|2697811307|ref|XM_016952902.4| PREDICTED: Pan troglodytes telomerase reverse transcriptase (TERT), transcript variant X1, mRNA
    length: 4043
    E Value: 0.0
    CTCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CTCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCC...
    ****ALIGNMENT****
    Sequence: gi|2697811308|ref|XR_010157111.1| PREDICTED: Pan troglodytes telomerase reverse transcriptase (TERT), transcript variant X2, misc_RNA
    length: 4745
    E Value: 0.0
    TCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCCC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCCC...
    ****ALIGNMENT****
    Sequence: gi|2697811309|ref|XR_001717411.4| PREDICTED: Pan troglodytes telomerase reverse transcriptase (TERT), transcript variant X3, misc_RNA
    length: 4020
    E Value: 0.0
    TCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCCC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCCC...



```python
# For the human TERT gene compared to all 3 varients of the Chimp TERT gene, E-Value = 0.0
```

# Open CV
### Part 1


```python
# Import the packages we need
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# pip install --upgrade pip
# pip install opencv-python
import cv2
```


```python
# Load in the desired image
# I picked Spiderman because I like him 
img = cv2.imread("Spiderman.jpg")
```


```python
# If ran correctly it will show as this type
type(img)
```




    numpy.ndarray




```python
# It can run even without the image
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
# But it will show as this type
type(img_wrong)
```




    NoneType




```python
# This will show the image, but matplotlib will change the colors
# This is an issue with how CV views RBG images
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f6ac623b510>




![png](output_26_1.png)



```python
# This changes the order CV reads the image and fixes the colors
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7f6ac49cd810>




![png](output_9_1.png)



```python
# We're going to make the image grayscale
img_gray = cv2.imread("Spiderman.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (272, 185)




```python
# Like before the color isn't correct
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7f6ac493abd0>




![png](output_11_2.png)



```python
# This will convert the image to grayscale
# grayscale images have less information so converting images to grayscale can make them less bulky
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f6ac49209d0>




![png](output_12_1.png)



```python
# We can resize images
fix_img.shape
```




    (272, 185, 3)




```python
new_img = cv2.resize(fix_img, (1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f6ac3083990>




![png](output_14_2.png)



```python
new_img.shape
```




    (400, 1000, 3)




```python
# We can also change the image size by changing the ratio
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f6ac30677d0>




![png](output_17_2.png)



```python
new_img.shape
```




    (136, 92, 3)




```python
# We can flip images
# This will just flip the image down
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7f6ac2fc9f50>




![png](output_19_2.png)



```python
# This will flip the image down and backwards
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7f6ac2f2cc10>




![png](output_20_1.png)



```python
type(fix_img)
```




    numpy.ndarray




```python
# We can save the manipulated image
# This includes the color change, so when opened in JupyterLabs Spiderman will be the wrong color again
# But when opened using cv it will be the 'correct' color"
cv2.imwrite('Spiderman_fixed_image.jpg', flip_img)
```




    True




```python
# Spiderman as he should be
final_spiderman = cv2.imread("Spiderman_fixed_image.jpg")
plt.imshow(final_spiderman)
```




    <matplotlib.image.AxesImage at 0x7f6ac2f11d90>




![png](output_23_1.png)



### Part 2


```python
# This is the same thing we did before to fix the color of the image
img = cv2.imread("Spiderman.jpg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f6ac2e7ef50>




![png](output_26_1.png)



```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f6ac2e55ed0>




![png](output_28_1.png)



```python
# We can change the image to other file types
# HSV is hue saturation value
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f6ac2db5890>




![png](output_30_1.png)



```python
# This is another way to change the file color
# HLS is hue light saturation
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f6ac2d96650>




![png](output_32_1.png)



```python
# This is an example of a watermark
img1 = cv2.imread('Venom.jpg')
img2 = cv2.imread('Spiderman.jpg')
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f6ac2cf5290>




![png](output_34_1.png)



```python
# Let's fix the color on the image
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f6ac2ccefd0>




![png](output_36_1.png)



```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f6ac2c33a50>




![png](output_37_1.png)



```python
# Let's resize the images
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200,1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
# This overwrites the two images where they are equal sizes and transparency
# This is useful for things like microscopy when you need to merge images
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f6ac2c17e90>




![png](output_41_1.png)



```python
# alpha is how transparent the images are so we can adjust the transparency of the images
alpha = 0.8
beta = 0.2

blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f6ac2387810>




![png](output_43_1.png)



```python
# Let's try imposing the images over each other in a different way
# First we need to define the images, fix the color, and adjust the size of our smaller image
img1 = cv2.imread('Venom.jpg')
img2 = cv2.imread('Spiderman.jpg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (100,100))
```


```python
# Then we set the position of our smaller image over the larger image
large_img = img2
small_img = img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7f6ac2bb52d0>




![png](output_45_1.png)

### Part 3

```python
# https://github.com/worklifesg/Python-for-Computer-Vision-with-OpenCV-and-Deep-Learning
```


```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f6ac2b91a50>




![png](output_49_1.png)



```python
# Let's do some thresholding
img = cv2.imread('rainbow.jpg', 0)
```


```python
# This gives a gradient
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f6ac2af6d90>




![png](output_51_1.png)



```python
# Thresholding will make it just white and black
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
ret1
```




    127.0




```python
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f6ac2ae45d0>




![png](output_54_1.png)



```python
# This will change the limits
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f6ac2a42990>




![png](output_55_1.png)



```python
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f6ac2a25890>




![png](output_56_1.png)



```python
# Let's do a real world application of background removal
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f6ac29854d0>




![png](output_57_1.png)



```python
# This function cuts out other aspects we have to do
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```


![png](output_59_0.png)



```python
# This threshold does not work well for the purposes of this example
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_60_0.png)



```python
# Let's try another one
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_61_0.png)



```python
# This uses the adaptive threshold, which is a little clearer
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
show_pic(th2)
```


![png](output_63_0.png)



```python
# This layers two of the thresholds we made
blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_64_0.png)



```python
# Let's look at gaussian
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_65_0.png)

# Corner Detection

```python

# pip install --upgrade pip
# pip install opencv-python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
flat_chess = cv2.imread('ChessBoardGreen.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fb3d53c0b10>




![png](output_1_1.png)



```python
# For our purposes we're going to make the chess board grayscale
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fb3d4300f10>




![png](output_2_1.png)



```python
# Now let's load in our other chess board and fix the colors
real_chess = cv2.imread("ChessBoard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fb3d426b310>




![png](output_4_1.png)



```python
# We're also converting this to grayscale
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fb3d5467ed0>




![png](output_5_1.png)



```python
# Now we are going to do Harris corner detection
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
# When the results are greater than 1% of the max value it will be asigned the color red
flat_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fb3d4249550>




![png](output_7_1.png)



```python
# Now we're going to try it on the real chess board
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fb3d41ac7d0>




![png](output_8_1.png)



```python
# Now let's try Shi-Tomasi Corner Detection

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y), 3, (255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fb3d4198390>




![png](output_10_1.png)



```python
# Let's try it on the real chess board
corners = cv2.goodFeaturesToTrack(gray_real_chess, 64, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fb3d40f9b10>




![png](output_11_1.png)

# Edge Detection

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# Spider man is back for this one
img = cv2.imread("Spiderman.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fb3d406ca10>




![png](output_13_1.png)



```python
# This is using the median threshold values
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb3d4051910>




![png](output_14_1.png)



```python
# Let's clean up the image
med_value = np.median(img)
med_value
```




    250.0




```python
# Let's changes the lower and upper limits of the threshold
lower = int(max(0, 0.7*med_value))
upper = int(min(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb3ccd599d0>




![png](output_16_1.png)



```python
# Let's continue to adjust the threshold
edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb3ccd38ed0>




![png](output_17_1.png)



```python
# Let's look at bluring get more clarity
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb3ccc9db50>




![png](output_18_1.png)



```python
# We can also increase the kernel size for better clartiy, but for my image this is not necessary 
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb3ccc88510>




![png](output_19_1.png)

# Feature Matching 

```python
# pip uninstall opencv-python
# pip install opencv-contrib-python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display(img, cmap='gray'):
    fig = plt.figure(figsize = (12, 10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
# For this example I am using action figures from the Across the Spiderverse movie and detecting Miles Morales
Miles_Morales = cv2.imread("Miles_Morales.jpg", 0)
display(Miles_Morales)
```


![png](output_22_0.png)



```python
SpiderVerse = cv2.imread("Many_Spidermen.jpg", 0)
display(SpiderVerse)
```


![png](output_23_0.png)



```python
# The first method we're trying is brute force matching
# orb is looking for key points in the thing we're searching for
orb = cv2.ORB_create()

kp1, des1 = orb.detectAndCompute(Miles_Morales, mask = None)
kp2, des2 = orb.detectAndCompute(SpiderVerse, mask = None)
```


```python
# This is collecting our matches
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
```


```python
# Now let's sort our matches
matches = sorted(matches, key = lambda x:x.distance)
```


```python
Miles_Morales_Matches = cv2.drawMatches(Miles_Morales, kp1, SpiderVerse, kp2, matches [:25], None, flags = 2 )
```


```python
display(Miles_Morales_Matches)
```


![png](output_28_0.png)



```python
# As you can see the match isn't perfect, so we can refine it
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(Miles_Morales, None)
kp2, des2 = sift.detectAndCompute(SpiderVerse, None)
```


```python
# Then apply the sifters to the brute force
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
```


```python
# Now we will apply a ratio test
good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('Length of total matches:', len(matches))
print('Length of good matches:', len(good))
```

    Length of total matches: 7929
    Length of good matches: 47



```python
sift_matches = cv2.drawMatchesKnn(Miles_Morales, kp1, SpiderVerse, kp2, good, None, flags = 2)
display(sift_matches)
```


![png](output_34_0.png)



```python
# Let's try flanmatcher
# same as before sift through both images
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(Miles_Morales, None)
kp2, des2 = sift.detectAndCompute(SpiderVerse, None)
```


```python
# Then identify details in both images
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, tree = 5)
search_params = dict(checks = 50)
```


```python
# Apply a ratio test for good matches and compare the two images
flann = cv2.FlannBasedMatcher(index_params, search_params)

match = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
        
flann_matches = cv2.drawMatchesKnn(Miles_Morales, kp1, SpiderVerse, kp2, good, None, flags = 0)
display(flann_matches)
```


![png](output_37_0.png)



```python
# Let's add a mask to it
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(Miles_Morales, None)
kp2, des2 = sift.detectAndCompute(SpiderVerse, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees = 5)
search_params = dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
for i, (match1, match2) in enumerate(matches):
    if match1.distance < 0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0, 255, 0),
                  singlePointColor = (255, 0, 0),
                  matchesMask = matchesMask,
                  flags = 0)
```


```python
flann_matches = cv2.drawMatchesKnn(Miles_Morales, kp1, SpiderVerse, kp2, matches, None, **draw_params)

display(flann_matches)
```


![png](output_43_0.png)

# Object Detection

```python
# python -m pip install --upgrade pip
# pip install opencv_python
```


```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# We are going to detect objects in images
# I am using carnations for this example
# This is the training image
full = cv2.imread('Training.jpg')
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7fb3cc52e090>




![png](output_46_1.png)



```python
# This is the test image
test = cv2.imread('Testing.jpg')
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7fb3cc502950>




![png](output_47_1.png)



```python
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (230, 219, 3)
    Training image shape: (851, 1280, 3)



```python
# These are all out different methonds
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
for m in methods:
    
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
        
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255, 0, 0), 10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')
    
    plt.suptitle(m)
    
    plt.show()
    print('\n')
    print('\n')
```


![png](output_50_0.png)


    
    
    
    



![png](output_50_2.png)


    
    
    
    



![png](output_50_4.png)


    
    
    
    



![png](output_50_6.png)


    
    
    
    



![png](output_50_8.png)


    
    
    
    



![png](output_50_10.png)


    
    
    
    




