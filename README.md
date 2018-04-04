# blast_getLCA
## Introduction
The blast_getLCA algorithm parses a blast file (format 6), and assigns each read to the lowest common ancestor of the best hit(s) in the blast file.

### Example
```python blast_getLCA.py -t 0.9 test.taxid.blast```

The output is a tsv file, with one row for each read, and the following columns: 

- **qseqid**  The identifier of the query sequence 
- **LCA**   The lowest common ancestor of the read, with info on genus, family, order and class of the assigned taxon
- **rank** Rank of the lowest common ancestor
- **TaxIDs** TaxIDs of the the best matching reference sequences in the database used to infer the lowest common ancestor, separated by colon
- **Stats** Basic statistics of the assignment (separated by underscore): 
  - *Tothits* = total number of alignments for the read, 
  - *accepted-hits* = total number of best (equally good) hits for the read.
  - *Min-Nm* = Edit distance for the best hit(s)
  - *IDp* = Identity percentage for the best hit(s)
 - **qseqlength**  Query sequence length
 - **IDp** Identity percentage for the best hit(s)
 - **gapmm** Gaps, mismatches and length of the best alignment(s), separated by underscore
 - **drop** Information on whether the LCA have been dropped to genus or family level based on a low identity percent to the best hit(s)

## Options
```
Options:
  -h, --help            show this help message and exit
  -i WRONGTAX, --ignoretaxid=WRONGTAX
                        csv file of taxids to ignore, first column should
                        contain taxids
  -t THRESHOLD, --threshold=THRESHOLD
                        Ignores reads where the best alignment has less than a
                        certain percentage identity to the refenrence
                        [default=0.95]
```
## Download and Install

#### 1) Download and unpack the blast_getLCA master
```
wget https://github.com/frederikseersholm/blast_getLCA/archive/master.tar.gz
tar xpvf master.tar.gz
cd blast_getLCA-master
```
#### 2) Run the script prepare_getLCA_script.sh
This script downloads the NCBI taxonomy files nodes.dmp and names.dmp to the folder 'taxdump', and adds the paths to the files to the python script.
```
bash prepare_getLCA_script.sh
```
#### 3) Run the script on the supplied test blast-file, modified for the blast_getLCA algorithm
This script downloads the NCBI taxonomy files nodes.dmp and names.dmp to the folder 'taxdump', and adds the paths to the files to the python script.
```
python blast_getLCA.py test.taxid.blast
```

## Generate a suitable blast file for analysis

#### 1) Blast fasta file with custom output format 6
In order for the blast_getLCA script to work, fasta files should be blasted with the following output format:
  - outfmt "6 qseqid sacc sseqid pident qlen length mismatch gapopen gaps evalue bitscore nident"
```
blastn -outfmt "6 qseqid sacc sseqid pident qlen length mismatch gapopen gaps evalue bitscore nident" -max_target_seqs 100 -reward 1 -db $DB -query ${FILENAME}.fasta -out ${FILENAME}.blast
```
