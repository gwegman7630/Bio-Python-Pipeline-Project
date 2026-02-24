#STEP 2
#download full genome from NCBI
#use biopython to parse through the file and
#build a fasta file with each header as the protein ID

#things:
#if you define an input or output in the rule, you can use it in run

#rule all must include everything the code needs to run for input
#and output checks that the code output everythign it was supposed to, if a file doesnt exist that is in output, it errors
#COMMON ERRORS:
#NO spaces in between rule,input,output,shell
#ALWAYS indent input,output,shell, and EVERYTHING under that unless it doesnt give u an error (i think only for variables?)
rule all:
    input:
        'cds_count.txt',
        'fdr05_results.txt',
        'countbowtie2.txt',
        'catit.txt',
        'reference_transcriptome.fasta',
        'ncbi/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna',
        'PipelineReport.txt'
#apparently you dont put stuff in the directory in rule all, the input should be what snakemake is gonnabuild
#pull only the cds regions of the desired protein and make sure that the cds file is THERE
rule cds:
    output:
        "ncbi/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna"
    shell:
        """
        mkdir -p ncbi
        datasets download genome accession GCF_000845245.1 --include cds --filename ncbi/genome.zip
        unzip ncbi/genome.zip -d ncbi
        """
rule protein_ID: #this uses biopython with some fancy parsing to create a new fasta file where each record id is the protein ID yp
    input:
        fun ='ncbi/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna'
    output:
        fun = 'reference_transcriptome.fasta'

    params:
        window = 83 #website recommeneded this, not sure what it does

    run:
        from Bio import SeqIO #gonna use seqio to parse through each record and the string of each record id

        records = [] #this is my list that will store everything that i edit and i will use it to write to an outfile

        for record in SeqIO.parse(input.fun , 'fasta'): #DO NOT CODE HARD PATHS
            search = record.description #record.description is basically the header, i want to search for protein id in header, and then figure out how to parse through the string syntax

            if 'protein_id=' in search: #i know each description will have protein id, so i can use this if statement to verify and also begin parsin

                pid = search.split('protein_id=')[1].split(']')[0] #important that i break this down
            #by splitting the string before and after protein_id. 0 will index what came before, and 1 will index what came after, i want the protein id which comes after so index 1
            #i need to split again by the bracket which comes at the end of [protein_id=YP_0908] and then since i want everything before it i index 0 which gives me the code
            #super useful for the future

                record.id = pid #these two lines tell the fasta file? i guess, that the name and id will now be the variable which is the current protein ID i cut out
                record.name = pid
                record.description="" #this makes sure there is nothing else left in the record

                records.append(record) #now each record name and id is the protein ID, and each record has the nucloetides under it, so i can append it to a list for each record and write
        SeqIO.write(records, output.fun, 'fasta') #now write out a new file use seqio and specify that its fasta
#make an output file that states how many CDS sequences are in the HCMV genome
rule cds_count:
    input:
        fasta_in = 'reference_transcriptome.fasta'

    output:
        funny = 'cds_count.txt'

    shell: """
    echo "The HCMV Genome (GCF_000845245.1) has $(grep -c '>' {input.fasta_in}) CDS." > {output.funny} #echo will write/print something to a designated file, i use one > since there is noting in pipeline.txt rn, to not overwrite i do >>
    """
#grep -c will count every line containing a > (every fasta file starts with this)

#making the kallisto index so i can quantify TPM in each sample
rule kallisto:
    input:
        i = 'reference_transcriptome.fasta'
    output:
        o = 'index.idx'
    shell: """
        kallisto index -i {output.o} {input.i}
        """
#alright the fact that you cant have any spaces inbetween inputs out puts and shell is frustrating
#quantify TPM in each sample using kallisto, using the previous index built
rule quantify:
    input:
        I_I = 'SRR5660030_1.fastq',
        I_II = 'SRR5660030_2.fastq',
        II_I = 'SRR5660033_1.fastq',
        II_II = 'SRR5660033_2.fastq',
        III_I = 'SRR5660045_1.fastq',
        III_II = 'SRR5660045_2.fastq',
        IV_I = 'SRR5660044_1.fastq',
        IV_II = 'SRR5660044_2.fastq',
        index = 'index.idx'
    output:
        o1 = directory('results/SRR5660030'),#had to make them all directories using that command since snakemake was throwing a fit saying those didnr exist (because it created them as files, not directories)
        o2 = directory('results/SRR5660033'),
        o3 = directory('results/SRR5660045'),
        o4 = directory('results/SRR5660044')
    shell: """
    mkdir -p results
    kallisto quant -i {input.index} -o {output.o1} -b 30 -t 4 {input.I_I} {input.I_II}
    kallisto quant -i {input.index} -o {output.o2} -b 30 -t 4 {input.II_I} {input.II_II}
    kallisto quant -i {input.index} -o {output.o3} -b 30 -t 4 {input.III_I} {input.III_II}
    kallisto quant -i {input.index} -o {output.o4} -b 30 -t 4 {input.IV_I} {input.IV_II}
    """
#use 30 bootstraps and four threads (standard), TPM quantification is simple command line script
#snakemake wasnt making the results files from my specified output, so i had to do it manually in shell 

#this is gonna use sleuth and an R script to basically output differential expression between 2dpi and 6dpi for the samples
#the sleuth package is going to use sleuth_table.tsv (which gives proper headers condition etc.) to make the results table
rule sleuth:
    input:
        slth = 'sleuth.R',
        tab = 'sleuth_table.tsv',
        wait = directory('results/SRR5660030'),
        wai = directory('results/SRR5660033'),
        wa = directory('results/SRR5660044'),
        w = directory('results/SRR5660045')
    output:
        'fdr05_results.txt'
    shell: """
    Rscript {input.slth}
    """
#you can run R in shell (or at least on the server) using Rscript, which is helpful to run whatever sript you have in your working directory
#the same way i downloaded the OG genome i need to do again
#but i only downloaded cds last time, i dont wanna change stuff around cuz it works so i am just gonna download the full thing and pull from the HCMV file so it doesnt overwrite
rule genome:
    output: 
        'HCMV/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna'
    shell:"""
    mkdir -p HCMV
    datasets download genome accession GCF_000845245.1 --filename HCMV/genome.zip
    unzip HCMV/genome.zip -d HCMV
    """
#you can make dir in shell snakemake but you CANT cd in snakemake, why tho?
#i dont really know how else to make the bowtie two index using HCMV genome
rule bowtieindex:
    input:
        genome = 'HCMV/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna'
    output:
        q = 'HCMV_genome_index.1.bt2',
        l = 'HCMV_genome_index.2.bt2',
        k = 'HCMV_genome_index.3.bt2',
        v = 'HCMV_genome_index.4.bt2',
        f = 'HCMV_genome_index.rev.1.bt2',
        w = 'HCMV_genome_index.rev.2.bt2'
#i guess i need to put all these files in (which i checked if they generated by doing this locally) so snakemake can confirm these are downloaded before moving on
    shell: """
    bowtie2-build {input.genome} HCMV_genome_index
    """
#building the bowtie2 genome so i can map reads to genome
#now im gonna map reads to the genome
#i dont need to call all .bt2 files just name of index so snakemake will be like hey we need this .bt2 file, go to the rule that makes it, which will in turn generat all the other bt2 files i need
rule bowtie2:
    input:
        bindex = 'HCMV_genome_index.1.bt2',
        I_I = 'SRR5660030_1.fastq',
        I_II = 'SRR5660030_2.fastq',
        II_I = 'SRR5660033_1.fastq',
        II_II = 'SRR5660033_2.fastq',
        III_I = 'SRR5660045_1.fastq',
        III_II = 'SRR5660045_2.fastq',
        IV_I = 'SRR5660044_1.fastq',
        IV_II = 'SRR5660044_2.fastq'
    output:
        I = 'SRR5660030.sam',
        II = 'SRR5660033.sam',
        III = 'SRR5660044.sam',
        IV = 'SRR5660045.sam',
        I_Im = 'SRR5660030_mapped_1.fq',
        I_IIm = 'SRR5660030_mapped_2.fq',
        II_Im = 'SRR5660033_mapped_1.fq',
        II_IIm = 'SRR5660033_mapped_2.fq',
        III_Im = 'SRR5660045_mapped_1.fq',
        III_IIm = 'SRR5660045_mapped_2.fq',
        IV_Im =  'SRR5660044_mapped_1.fq',
        IV_IIm = 'SRR5660044_mapped_2.fq'
    shell:"""
    bowtie2 -x HCMV_genome_index -1 {input.I_I} -2 {input.I_II} -S {output.I} --al-conc SRR5660030_mapped_%.fq
    bowtie2 -x HCMV_genome_index -1 {input.II_I} -2 {input.II_II} -S {output.II} --al-conc SRR5660033_mapped_%.fq
    bowtie2 -x HCMV_genome_index -1 {input.III_I} -2 {input.III_II} -S {output.III} --al-conc SRR5660045_mapped_%.fq
    bowtie2 -x HCMV_genome_index -1 {input.IV_I} -2 {input.IV_II} -S {output.IV} --al-conc SRR5660044_mapped_%.fq
    """
#-al--conc will write reads that only mmapped to genome to a separate fastq file, dont think i need to do it for both forward and reverse using % but its easiest for now
rule countbowtie2:
    input:
        I = 'SRR5660030_1.fastq',
        II = 'SRR5660033_1.fastq',
        III = 'SRR5660044_1.fastq',
        IV = 'SRR5660045_1.fastq',
        Im = 'SRR5660030_mapped_1.fq',
        IIm = 'SRR5660033_mapped_1.fq',
        IIIm = 'SRR5660045_mapped_1.fq',
        IVm =  'SRR5660044_mapped_1.fq'
    output:
        e = 'countbowtie2.txt'
    shell:"""
    echo "Sample SRR5660030 had $(( $(wc -l < {input.I}) / 4 )) read pairs before and $(( $(wc -l < {input.Im}) / 4 )) read pairs after Bowtie2 filtering." > {output.e}
    echo "Sample SRR5660033 had $(( $(wc -l < {input.II}) / 4 )) read pairs before and $(( $(wc -l < {input.IIm}) / 4 )) read pairs after Bowtie2 filtering." >> {output.e}
    echo "Sample SRR5660045 had $(( $(wc -l < {input.III}) / 4 )) read pairs before and $(( $(wc -l < {input.IIIm}) / 4 )) read pairs after Bowtie2 filtering." >> {output.e}
    echo "Sample SRR5660044 had $(( $(wc -l < {input.IV}) / 4 )) read pairs before and $(( $(wc -l < {input.IVm}) / 4 )) read pairs after Bowtie2 filtering." >> {output.e}
    """
#so this is gonna take the original fq and then mapped fq (using bt2) and its gonna tell shell to output how many lines contain '<' (which is each sequence)
#echo will write this to the file, i had to use $(( $)) to tell shell what i want done in the terminal and not have it be in the quotations
#echo it all to a txt file so i can cat it to my pipeline report later since again snakemake won let me edit the same file over the course of the workflow
rule spades:
    input:
        I_Im = 'SRR5660030_mapped_1.fq',
        I_IIm = 'SRR5660030_mapped_2.fq',
        II_Im = 'SRR5660033_mapped_1.fq',
        II_IIm = 'SRR5660033_mapped_2.fq',
        III_Im = 'SRR5660045_mapped_1.fq',
        III_IIm = 'SRR5660045_mapped_2.fq',
        IV_Im =  'SRR5660044_mapped_1.fq',
        IV_IIm = 'SRR5660044_mapped_2.fq'
#i just want to check i the folder is made, kinda silly bc i should tell snakemake to check if every correct thing is in there but thats a lot
    output: #i need to specify the EXACT file path in output or else snakemake makes a fuss when i try to call it later
        "SRR5660030_assembly/contigs.fasta",
        "SRR5660033_assembly/contigs.fasta",
        "SRR5660045_assembly/contigs.fasta",
        "SRR5660044_assembly/contigs.fasta"
    shell:"""
    spades.py -k 127 -t 2 --only-assembler -1 {input.I_Im} -2 {input.I_IIm} -o SRR5660030_assembly/
    spades.py -k 127 -t 2 --only-assembler -1 {input.II_Im} -2 {input.II_IIm} -o SRR5660033_assembly/
    spades.py -k 127 -t 2 --only-assembler -1 {input.III_Im} -2 {input.III_IIm} -o SRR5660045_assembly/
    spades.py -k 127 -t 2 --only-assembler -1 {input.IV_Im} -2 {input.IV_IIm} -o SRR5660044_assembly/
    """
#standard run of spades, use 2 threads at 127 kmer length, in only assembler mode to cut down on time and resources
#assemble using only the ones that mapped to the genome to make it easier/run less long
rule longestcontig:
    input:
        I = 'SRR5660030_assembly/contigs.fasta',
        II = 'SRR5660033_assembly/contigs.fasta',
        III = 'SRR5660045_assembly/contigs.fasta',
        IV = 'SRR5660044_assembly/contigs.fasta'
    output:
        L1 = 'SRR5660030_longest_contig.fasta',
        L2 = 'SRR5660033_longest_contig.fasta',
        L3 = 'SRR5660045_longest_contig.fasta',
        L4 = 'SRR5660044_longest_contig.fasta'
    run:

        from Bio import SeqIO #basically im gonna find the longest contig in each assembly
        length = 0 #baseline length

        contig = None #this will eb reset and written to file once longest contig is found
        #doing four for loops is stupid but at the time this was how i could think of generating the outputs i need

        for record in SeqIO.parse(input.I, 'fasta'): #go through each contig in each assembly

            if len(record) > length: #if the length of the current record beats the length of previous maximum (it will always beat it first time because 0)

                length = len(record) #then new longest contig is the current iteration and keep moving through the entire file, once file is parsed through the longest contig will have beaten out all the others

                contig = record #after this, make this max length record the new contig (resets it for everyone found but doesnt matter)
        SeqIO.write(contig, output.L1, 'fasta') #now write it to the proper output, yes you call outputs in python like you do in shell with the output. 
        length = 0 #reset length for each time
        for record in SeqIO.parse(input.II, 'fasta'):

            if len(record) > length: #do this for tyhe rest of the mapped fasta files

                length = len(record)

                contig = record
        SeqIO.write(contig, output.L2, 'fasta')
        length = 0
        for record in SeqIO.parse(input.III, 'fasta'):

            if len(record) > length:

                length = len(record)

                contig = record
        SeqIO.write(contig, output.L3, 'fasta')
        length = 0
        for record in SeqIO.parse(input.IV, 'fasta'):

            if len(record) > length:

                length = len(record)

                contig = record
        SeqIO.write(contig, output.L4, 'fasta')

# i am gonna use nucleotide to nucleotide blast because the longes contig is in nucleotide form and im searching the nucleotide database
#grabbing the genome for Betaherpesvirinae
rule betavisgenome:
    output:
        'betaviz/ncbi_dataset/data/genomic.fna'
    shell:"""
    mkdir -p betaviz
    datasets download virus genome taxon Betaherpesvirinae --filename betaviz/betavis.zip
    unzip -o betaviz/betavis.zip -d betaviz
    """
#-o and -p overwrite the files if they already exist, since these arent integral files and only exist for the workflow idc if they get overwritten each run
rule betavizblast:
    input:
        i = 'betaviz/ncbi_dataset/data/genomic.fna'
    output:
        'betaviz/betaviz_db.nsq',
        'betaviz/betaviz_db.nin',
        'betaviz/betaviz_db.nhr' #so you cannot pecify outputs like betavizdb because it creates multiple files, so you need to actually write it out in the shell script and then tell snakemake that these r the files it needed to have generated
    shell:"""
    makeblastdb -in {input.i} -out betaviz/betaviz_db -title betavisdb -dbtype nucl
    """
#making my blast database to align the contigs to so i can get percent identity of this genome and other metrics
rule blasting:
    input:
        I = 'SRR5660030_longest_contig.fasta',
        II = 'SRR5660033_longest_contig.fasta',
        III = 'SRR5660045_longest_contig.fasta',
        db = 'betaviz/betaviz_db.nsq', #needs to flag that this file is made (if its made then the other ones are for sure), if its not made, snakemake will go to the rule that makes it
        IV = 'SRR5660044_longest_contig.fasta'
    output:
        R1 = 'betaviz_results1.csv',
        R2 = 'betaviz_results2.csv',
        R3 = 'betaviz_results3.csv',
        R4 = 'betaviz_results4.csv'

    shell:""" #same thing down here, blast needs the common names of the files, not output.I specified as a file, because its NOT a file
    blastn -query {input.I} -db betaviz/betaviz_db -out {output.R1} -outfmt "6 qseqid sacc pident length qstart qend sstart send bitscore evalue stitle"
    blastn -query {input.II} -db betaviz/betaviz_db -out {output.R2} -outfmt "6 qseqid sacc pident length qstart qend sstart send bitscore evalue stitle"
    blastn -query {input.III} -db betaviz/betaviz_db -out {output.R3} -outfmt "6 qseqid sacc pident length qstart qend sstart send bitscore evalue stitle"
    blastn -query {input.IV} -db betaviz/betaviz_db -out {output.R4} -outfmt "6 qseqid sacc pident length qstart qend sstart send bitscore evalue stitle"
    """ #yeah i have to just say the name of the database not as input.variable
#i had to use outfmt 6 to make it tab separated instead of comma because stitle has spaces after which makes more columns, so when parsing with pandas it throws error later
#taking my results and writing top 5 hits to a tab separated file so i can cat them latr
rule writingblasttofile:
    input:
        R1 = 'betaviz_results1.csv',
        R2 = 'betaviz_results2.csv',
        R3 = 'betaviz_results3.csv',
        R4 = 'betaviz_results4.csv'
    output:
        O1 = 'top5_1.tsv',
        O2 = 'top5_2.tsv',
        O3 = 'top5_3.tsv',
        O4 = 'top5_4.tsv'
    run:
        import pandas as pd

        myresults = [input.R1,input.R2,input.R3,input.R4] #i have to make an input and an output file list so i can call them and keep syntax with snakemake
        outputfiles = [output.O1,output.O2,output.O3,output.O4]
        headers = ["qseqid","sacc","pident","length","qstart","qend","sstart","send","bitscore","evalue","stitle"] #specified by assignment
        #outfmt -10 DOES NOT produce headers, so i have to make a list of the headers and then do cols = names to set them as headers in these tsv files
        for i in range(len(myresults)): #i ma iterating through the length so i can pull out both the first input file, grab top 5 hits, and then paste it into first output file using same index
            frameit = pd.read_csv(myresults[i], sep="\t", names =headers) #this reads in csv and sets it to variable also had to change to reading in tab delimited since i changed outfmt to 6
        
            grabtop = frameit.sort_values('bitscore', ascending = False) #sort the csv file by bitscore to get the top 5 hits and make ascending = False so tyhey are sorted largest to smallest

            top5 = grabtop.groupby('qseqid').head(5) #this is where it grabs the top 5 
            top5 = top5.drop(columns=['qseqid']) #now that ive sorted the top 5 by qseqid i dont need it anymore and assignment doesnt want me including it in output
            top5.to_csv(outputfiles[i], sep='\t', index = False) #make a new tsv file by noting sep='\t' and indexing proper output file
#index = False makes it so that column numbers dont show
#combinging all files into one txt output that i can echo to my pipeline report later
rule catit:
    input:
        I = 'top5_1.tsv',
        II = 'top5_2.tsv',
        III = 'top5_3.tsv',
        IV = 'top5_4.tsv'

    output:
        o = 'catit.txt'

    shell: """
    echo 'SRR5660030:' >> {output.o}
    cat {input.I} >> {output.o}

    echo 'SRR5660033:' >> {output.o}
    cat {input.II} >> {output.o}

    echo 'SRR5660044:' >> {output.o}
    cat {input.IV} >> {output.o}

    echo 'SRR5660045:' >> {output.o}
    cat {input.III} >> {output.o}
    """
#double >> appends to, not overwrites

#files to cat
#'cds_count.txt', 'fdr05_results.txt', 'countbowtie2.txt', 'catit.txt
#so i learned that you cant have snakemake writing to the same file in multiple rules without order specification, so i wrote to unique files and now im gonna cat them to pipeline.txt here
#these are the necessary files that pipeline report requires to be built, snakemake will traceback all rules needed to generate these files so pipelinereport can be made
rule letsdothis:
    input:
        i = 'cds_count.txt',
        ii = 'fdr05_results.txt',
        iii = 'countbowtie2.txt',
        iv = 'catit.txt'
    output:
        o ='PipelineReport.txt'
    shell:"""
    cat {input.i} > {output.o}
    cat {input.ii} >> {output.o}
    cat {input.iii} >> {output.o}
    cat {input.iv} >> {output.o}
    """
#use cat and >> so it doesnt overwrite
rule clean: #this is so i can do freshs start runs, will delete all intermediate and output files
    shell:"""
        rm -rf ncbi_dataset
        rm -f reference_transcriptome.fasta
        rm -f ncbi_dataset.zip
        rm -f PipelineReport.txt
        rm -f cds_count.txt
        rm -f fdr05_results.txt
        rm -f countbowtie2.txt
        rm -f catit.txt
        rm -f *.sam
        rm -f md5sum.txt
        rm -rf ncbi
        rm -rf betaviz
        rm -f index.idx
        rm -rf HCMV
        rm -f *.bt2
        rm -rf results
        rm -f betaviz*
        rm -f *_mapped*
        rm -rf *_assembly
        rm -f top*
        rm -f *longest*
        """
#so if my snakemake runs halfway and stops, i can run clean and the -f part will make it so that shell wont throw an error if the file doesnt exist
