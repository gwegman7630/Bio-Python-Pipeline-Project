TOOLS NEEDED:
Python: https://www.python.org/
R: https://www.r-project.org/
kallisto: https://github.com/pachterlab/kallisto
BLAST: https://github.com/ncbi/blast_plus_docs
Bowtie2: https://github.com/BenLangmead/bowtie2
SPAdes: https://github.com/ablab/spades
SRA Toolkit: https://github.com/ncbi/sra-tools

DOWNLOADING FULL SAMPLE DATA:
This repo provides you with a subsample of the actual data used so that the runtime is produced.
The full sample data was retrieved via these terminal commands:
prefetch SRR5660030
fasterq-dump SRR5660030
prefetch SRR5660033
fasterq-dump SRR5660033
prefetch SRR5660044
fasterq-dump SRR5660044
prefetch SRR5660045
fasterq-dump SRR5660045
NOTE: The SRA toolkit is needed for downloading SRA files in this way. See tools needed section.

RUNNING WORKFLOW:

1. Clone the repo to your desired working directory, this will create a directory called Bio-Python-Pipeline-Project. Run in terminal:
  git clone https://github.com/gwegman7630/Bio-Python-Pipeline-Project
2. Now change your working directory by running in terminal:
  cd Bio-Python-Pipeline-Project/
4. Now run the snakemake workflow. In terminal run this command:
  snakemake -s done.smk --cores 4 --latency-wait 120 -p
5. Output file to view results is WegmanGabriel_PipelineReport.txt

MISC:

If for any reason you need to stop the workflow while its running hit CTRL + C
Then run this in terminal: 
  snakemake -s done.smk --cores 4 clean


