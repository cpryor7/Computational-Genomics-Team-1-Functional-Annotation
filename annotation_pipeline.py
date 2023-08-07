import glob
import subprocess
import pathlib
import os
import sys
import shutil


pred_faa = glob.glob('/home/team1/prediction/final_results_prodigal/faas/*.faa') # Location of predicted .faa files (to make directories)
ids = [x.split("/")[-1].split('.')[0] for x in pred_faa]


def eggnog():
    '''
    run eggnog-mapper with faa files from gene prediction group, generating annotated gff file
    this script assumes that eggnog-mapper is already installed
    '''
    
    os.system("cd /home/team1/annotation")
    os.system("mkdir eggnog_output")
    os.system("cd /home/team1/annotation/eggnog-mapper")
    os.system("export EGGNOG_DATA_DIR=/home/team2/databases/eggnog")
    os.system("for i in pred_faa; do python emapper.py -m diamond -i $i -o $i --output_dir /home/team1/annotation/eggnog_output --cpu 5 --pfam_realign realign; done")


def card_rgi():
    '''
    setting CARD database directory
    rgi load -i /path/to/card.json --local
    rgi search against antibiotics resistance gene
    rgi main --input_sequence /path/to/*.faa --output_file /path/to/result --local --clean -t protein
    '''
    p = subprocess.Popen(['cd', '/home/team1/prediction/final_results_prodigal/faas'])
    p = subprocess.Popen(['rgi','load','-i','/home/team1/annotation/card.json','--local'])	
    p=subprocess.Popen(['for','i','in','*.faa',';','do','rgi', 'main', '--input_sequence', '$i', '--output_file', '/home/team1/annotation/card_rgi/$i/result', '--local', '--clean', '-t','protein',';','done'])	

def signalp(ids):
    '''
    #installing signalp!
    pip install signalp-6-package/
    #Copy the model files to the location at which the signalp module got installed.
    SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
    cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
    '''
    for file in ids:
        p = subprocess.Popen(["signalp6", "--fastafile", "f{file}", "--organism","other","output_dir" , "/home/team1/annotation/signalp6_output/output1", "--format", "txt", "--mode", "fast", "--write_procs", "8"], stdout=2, stderr=subprocess.STDOUT)
def tmhmm(ids):
    for file in ids:
        os.system(f"biolib run DTU/DeepTMHMM --fasta {ids}.faa")


def blast():

    '''
    run BLAST with VFDB database, generating blast output and converting to gff files
    '''
    
    os.system("conda install -c bioconda blast")
    os.system("conda install -c bioconda mgkit")
    os.system("cd /home/team1/annotation/vfdb_output")
    os.system("makeblastdb -in /home/team1/annotation/vfdb_database/VFDB_setA_pro.fas -parse_seqids -dbtype prot -out /home/team1/annotation/vfdb_output/vfdb_prot")
    os.system("for i in CGT*.faa; do blastp -db vfdb_prot -query $i -out ${i}_out -outfmt 6 -num_threads 4 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1; done")
    os.system("for i in CGT*.faa_out; do blast2gff blastdb $i ${i%_*}.gff; done")


def merge(i): 
    '''
    merges all of the files in a directory using common IDs
    creates a new file in the output directory with the information from every tool for the IDs
    also creates a annotated fasta file for each
    '''

    output_file = open(f"/home/team1/annotation/output/merged_gff/{i}.gff", "w")
    merge_fasta = open(f"/home/team1/annotation/output/annotated_fasta/{i}_fasta", "w")
    track_list=[]
    gf_dic = {}

    fasta = open(f"/home/team1/prediction/final_results_prodigal/faas/{i}.faa", "r")
    fasta_dic = {}
    header_list = []
    sequence = []
    sequence_list = []
    count = 0
    for f in fasta:
        f_list = (f.split("\n"))
        index = 1
        head = ""
        if count == (len(f)+1):
            sequence_list.append(("").join(sequence))
        if "NODE" in f:
            sequence_list.append(("").join(sequence))
            header_list.append((f.split("#")[0]).replace(">","").strip())
            index = 0
            sequence = []
        if index > 0:
            sequence.append(f)
            index += 1
        count += 1

    gff_list =[]
    gff_dict = {}
    gff = open(f"/home/team1/prediction/final_results_prodigal/gffs/{i}.gff", "r")
    for g in gff:
        g_list = (g.split("\t"))
        for fi in g_list:
            if "NODE" in fi.split("#")[0]:
                gff_node = (fi.split("#")[0]).replace(">","")
                gff_list.append(gff_node)
                gff_dict[gff_node] = [g_list[2], "Start: " + g_list[3],"Stop: " +  g_list[4]]

    for x in range(0, len(gff_list)):
        current_header = header_list[x]
        gff_info = gff_dict[gff_list[x]]
        gf_dic[current_header] = gff_info
        try:
            fasta_convert =( header_list[x]+ " | " + gff_list[x] +"\n" +  sequence_list[x+ 1])
        except IndexError as e:
            fasta_convert = (header_list[x] + " | " + gff_list[x] + "\n")
        merge_fasta.write(fasta_convert + "\n")


    cr_final = {} #has a combination of the information from card_rgi and spaces if there isnt information for a header from the original fasta
    cr_nodes=[] #to store the ID from the card_rgi output file
    cr_dict={} # stores the ID as a key and the columns as values for the card_rgi output file
    cr = open(f"/home/team1/annotation/card_rgi/{i}.faa/result.txt", "r")
    for cr_line in cr:
        cr_lines = cr_line.split("\t")
        if ("NODE" in cr_lines[0]):
            cr_key = (cr_lines[0]).split("#")
            cr_nodes.append(cr_key[0].strip())
            cr_dict[cr_key[0].strip()]= ["CR-CutOff: " + cr_lines[5],"CR-Pass_bit: " + cr_lines[6], "CR-BestHitBit: " + cr_lines[7], "CR-BestHitARO: " + cr_lines[8],"CR-BestIdent:" + cr_lines[9], "CR-ARO: " + cr_lines[10],"CR-ModelType: " + cr_lines[11], "CR-SNPinARO: " + cr_lines[12],"CR-OtherSNP: " + cr_lines[13],"CR-DrugClass: " + cr_lines[14],"CR-ResMech: " + cr_lines[15],"CR-AMRFam: " + cr_lines[16],"CR-PredDNA: " +  cr_lines[17],"CR-%LenRefSeq: " + cr_lines[20],"CR-ID: " + cr_lines[21]]


    eggnog_final = {}
    eggnog_nodes = [] #stores the ID from eggnog domain/motif output
    eggnog_dict = {} #stores the ID as a key and the columns as values for the eggnog output file
    eggnog = open(f"/home/team1/annotation/eggnog_output/{i}.faa.emapper.annotations", "r")
    for eggnog_line in eggnog:
        eggnog_lines = eggnog_line.split("\t")
        if ("NODE" in eggnog_lines[0]):
            eggnog_key = (eggnog_lines[0].split("#"))
            eggnog_nodes.append(eggnog_key[0].strip())
            eggnog_dict[eggnog_key[0].strip()]= ["Eggnog: " + eggnog_lines[1],"Eggnog: " + eggnog_lines[2],"Eggnog: " + eggnog_lines[3]]


    tmhmm_final = {}
    tmhmm_nodes = []
    tmhmm_dict = {}
    tmhmm = open(f"/home/team1/annotation/tmhmm/finalresults/{i}.gff3", "r")
    for tmhmm_line in tmhmm:
        tmhmm_lines = tmhmm_line.split("\t")
        for t in tmhmm_lines:
            if ("Number of predicted TMRs" in t):
                num_tmr = t.split(":")[1]
                tmhmm_value = (t.split(":")[0])
                tmhmm_key = tmhmm_value.split()
                tmhmm_key = (tmhmm_key[1].strip())
                tmhmm_nodes.append(tmhmm_key)
                tmhmm_dict[tmhmm_key] = "NumPredTMR: " + num_tmr


    blast_final ={}
    blast_nodes = []
    blast_dict = {}
    blast =open(f"/home/team1/annotation/vfdb_output/{i}.faa.gff", "r")
    for blast_line in blast:
        blast_lines = blast_line.split("\t")
        if ("NODE" in blast_lines[0]):
            blast_nodes.append(blast_lines[0].strip())
            blast_dict[blast_lines[0].strip()]= ["BlastVirStart: " + blast_lines[3],"BlastVirStop: " + blast_lines[4]]


    signalp_final = {}
    signalp_nodes = [] #stores the ID from signalp
    signalp_dict = {} #stores the ID as a key and the columns as values for the signal p output file
    signalp = open(f"/home/team1/annotation/signalp6_output/output/{i}.faa/output.gff3", "r")
    for signalp_line in signalp:
        '''
        important information from signalp
        index 2: type of region
        index 3: predicted start signal peptide region
        index 4: predicted end signal peptide region

        '''
        signalp_lines = signalp_line.split("\t")
        #count = 0
        if ("NODE" in signalp_lines[0]):
            signalp_key = (signalp_lines[0].split("#"))
            signalp_nodes.append(signalp_key[0].strip())
            signalp_dict[signalp_key[0].strip()] = ["SigP-Region: " + signalp_lines[2], "SigP-RegionStart: " + signalp_lines[3],"SigP-RegionEnd: " + signalp_lines[4]]


    for id in header_list:
        if id in cr_dict.keys():
            cr_final[id] = cr_dict[id]
        else:
            cr_final[id] = [".",".",".",".",".",".",".",".",".",".", ".",".",".",".","."]
        if id in eggnog_dict.keys():
            eggnog_final[id] = eggnog_dict[id]
        else:
            eggnog_final[id] = [".",".","."]
        if id in tmhmm_dict.keys():
            tmhmm_final[id] = tmhmm_dict[id]
        else:
            tmhmm_final[id] = "."
        if id in blast_dict.keys():
                blast_final[id] = blast_dict[id]
        else:
            blast_final[id] = [".","."]
        if id in signalp_dict.keys():
            signalp_final[id] = signalp_dict[id]
        else:
            signalp_final[id] = [".",".","."]
        if(id in signalp_final.keys()) and (id in eggnog_final.keys()) and (id in cr_final.keys()) and (id in blast_final.keys()) and (id in tmhmm_final.keys()):
            merged = gf_dic[id] + blast_final[id] + cr_final[id] + eggnog_final[id] + signalp_final[id]
            output_file.write(id + "\t"  + "\t".join(merged)+ "\t" + tmhmm_final[id] + "\n")

    output_file.close()
    merge_fasta.close()


if __name__ == "__main__":
    eggnog()
    card_rgi()
    blast()
    signalp(ids)
    tmhmm()
    os.mkdir("/home/team1/annotation/output")
    os.mkdir("/home/team1/annotation/output/merged_gff") #location of merged_gff files
    os.mkdir("/home/team1/annotation/output/annotated_fasta") #location of annotated fasta files

    #merge the gff files!
    for i in ids:
        try: 
            merge(i)
        except FileNotFoundError as e:
           print(f"File {i} not found!")
           pass


