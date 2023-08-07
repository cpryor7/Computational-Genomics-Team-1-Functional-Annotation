# Team1-FunctionalAnnotation

Jiyeong Choi, Asmita Lagwankar, Chloe Pryor, Hannah Snyder, Likitha Venkatesh, Jiahong Zang

# Overview
Functional annotations describe the biochemical and biological function of proteins. There are homology-based and *ab intio* tools which are evidence based and using intrinsic characteristic of genome sequence. Our strategy is to run five different tools with FASTA files from gene prediction group and merge the outputs together. The chosen tools are eggNog-mapper, CARD-RGI, SingalP, Deep TMHMM, and VFDB.

# Running Tools
## Homology-based Tools
### eggNog-mapper
eggNog-mapper is a tool for fast functional annotation tool for novel sequences. It uses precomputed orthologus groups (OGs) and phylogenies from eggNOG database. Input is FASTA file and output files are annotations in gff file format. More information about eggNog-mapper can be found [here](https://github.com/eggnogdb/eggnog-mapper)

```
# setting eggnog database directory
export EGGNOG_DATA_DIR=/home/team2/databases/eggnog

# general usage
python emapper.py -m diamond -i {input file} -o {output file} --output_dir {output file directory}

# sample command for the analysis
python emapper.py -m diamond -i /home/team1/annotation/eggnog_output/CGT1058.faa -o /home/team1/annotation/eggnog_output/CGT1058.faa --output_dir /home/team1/annotation/eggnog_output --cpu 5 --pfam_realign realign
```

### CARD-RGI
The Comprehensive Antibiotic Resistance Database (CARD) is a biological database that collects and organizes reference information on antimicrobial resistance genes, proteins and phenotypes. RGI is the application uses reference data from the CARD. More information about CARD-RGI can be found [here](https://github.com/arpcard/rgi)

Input file is FASTA file and output is json and txt file.


```
# setting CARD database directory
rgi load -i /path/to/card.json --local

# rgi search against antibiotics resistance gene
rgi main --input_sequence /path/to/*.faa --output_file /path/to/result --local --clean -t protein
```
### VFDB
Virulence Factor Database (VFDB) is a database that allows get virulence factors in bacterial pathogens. BLAST can be used to utilize this tool. More information about VFDB can be found [here](http://www.mgc.ac.cn/VFs/)

```
#making blast database
makeblastdb -in /home/team1/annotation/vfdb_database/VFDB_setA_pro.fas -parse_seqids -dbytype prot -out /home/team1/annotation/vfdb_output/vfdb_prot

#generating blast output for all .faa files
for i in CGT*.faa; do blastp -db vfdb_prot -query $i -out “${i}_out” -outfmt 6 -num_threads 4 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1; done 

#converting blast output into .gff files
for i in CGT*.faa_out; do blast2gff blastdb $i “${i%_*}.gff”; done
```

## *Ab Initio* Tools
### SignalP
SignalP uses a deep neural network-based method for improved SP prediction. It predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms.

Input file is FASTA file and outputs are raw cleavage site score(C-score), signal peptide score (S-score), and combined cleavage site score(Y-score). More information about SignalP can be found [here](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md)

```
# Copy the model files to the location at which the signalp module got installed.
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/

# generate outputs
signalp6 --fastafile /home/team1/annotation/signalp6_output/CGT1572.faa --organism other --output_dir /home/team1/annotation/signalp6_output/output/CGT1572.faa --format txt --mode fast --write_procs 8
```
 

### Deep TMHMM
Deep TMHMM uses a deep learning protein language model-based algorithm. Its encoder consists of three components; a pre-trained language model (ESM-1b), a bi-directional LSTM and a dense layer with drop-out. This tool can take a protein sequence as input and outputs the corresponding per-residue sequence of labels. More information about Deep TMHMM can be found [here](https://dtu.biolib.com/DeepTMHMM)

```
# general usage 
biolib run DTU/DeepTMHMM --fasta input.fasta
```
Due to technical problem, only few files were run with command line and others were run in GUI mode. 



# Merging Outputs
To merge the GFF files we first parsed through the original fasta files to collect the complete list of sequence_id's. Next the output from each of the tools was uploaded and looped through. A dictionary was created for each of the tools with the key of the dictionary being the sequence_id and the dictionary was filled with a list corresponding to the valuable information from each tool. Next, the list of sequence_id's from the original fasta file were looped through. For each tool dictionary, the keys were searched for the corresponding sequence_id. If it was present, the information was added to a master dictionary, which had a key of the fasta sequence_id. If the key was not found for a tool, a list of "." the same length was added to show it was null (correlating with GFF formatting). The values in this master list were printed to an output for each corresponding isolate titled {isolate}.gff in the folder located at /team1/home/annotated/output/merged_gff. 

The python script, merge.py, was used to merge all outputs together.

For merging on the unannotated GFF file generated from Gene Prediction group, the job could be done through ```whole_merge.py``` with inputting sequence_id, result from eggnog, card-rgi, vfdb, signalP, tmhmm and unannotated GFF file from gene prediction. 

```
#example command for one contig
#generate sequence_id
grep '>' /path/to/faa > /path/to/store/id_file

#merging on the unannotated GFF
python whole_merge.py /path/to/id /path/to/eggnog_output /path/to/card_rgi_output /path/to/vfdb_output /path/to/signalp6_output /path/to/tmhmm_output /path/to/unannotated_gff > /path/to/store/whole_merge.gff; done
```
