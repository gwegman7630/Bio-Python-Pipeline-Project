TOOLS NEEDED:
Python: https://www.python.org/
R: https://www.r-project.org/
kallisto: https://github.com/pachterlab/kallisto
BLAST: https://github.com/ncbi/blast_plus_docs
Bowtie2: https://github.com/BenLangmead/bowtie2
SPAdes: https://github.com/ablab/spades

RUNNING WORKFLOW:

1. Clone the repo to your desired working directory, this will create a directory called Bio-Python-Pipeline-Project
2. in terminal run this command:
snakemake -s Bio-Python-Pipeline-Project/done.smk --cores 4 --latency-wait 120 -p
3. Output file to view results is WegmanGabriel_PipelineReport.txt

MISC:

If for any reason you need to stop the workflow while its running hit CTRL + C
Then run this in terminal:


