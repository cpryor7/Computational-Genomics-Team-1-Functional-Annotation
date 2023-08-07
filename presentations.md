# Background and Strategy Presentation
## Team 1 - 8:06 to 8:24 (18 min)
6 of 6 presenting
### General Topics to mention
  1. background, overview, and core concepts
  1. description of tasks, algorithms uses, qualitative comparisons
  1. proposed approach and analysis workflow
  1. task descriptions and delegations

### Results Expected (to plan)
1. predict various aspects of gene (protein) function
1. perform both ab initio and homology-based prediction as appropriate
1. aspects of function to predict – biochemical activity, molecular function,
(sub)cellular localization, domain and motif composition, higher level features
such as protein families or operons, enzymatic activity, virulence factors etc (note
that this list is not exhaustive)

### Deliverables Expected (to plan)
1. deliverable #1 – functional annotation pipeline on GitHub
1. deliverable #2 – TSV listing tool comparisons plus optional figure formats
1. deliverable #3 – GFF files for comparative genomics group
1. deliverable #4 – annotated FastA files for gene nucleotide and protein sequences

### Info Presented
- historically was done "by hand"
- 4 different tools and functions
  - 2 ab initio
  - 2 homology-based
1. SignalP
1. Deep TMHMM
1. eggNOG
1. CARD-RGI
- domain vs motif reviewed

### Comments and Suggestions
- I strongly encourage agreeing on a subset (maybe 5 smallest genome sizes) and set those for an earlier internal team deadline to allow downstream evaluations to happen. If you wait for all 50 genomes to be performed, the downstream work could run out of time if just 1 person was late.
- Team 1 needs to **add a homology-based virulence factor annotation**, so 5 instead of 4 tools to be used for annotation.
- Work delegation slide was good but revise based on my understanding of difficulty to even the workload out more:
    1. @lvenkatesh7 please be a secondary helper in the eggNOG work that @jchoi768 is the lead on
    1. @cpryor7 please identify an appropriate homology-based virulence factor annotation tool, or use the VFDB described [here](https://pubmed.ncbi.nlm.nih.gov/15608208/) for all genomes
    1. writing the pipeline delegations were not listed but each lead on a tool should be responsible for including their part 
