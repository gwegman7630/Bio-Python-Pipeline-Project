library(sleuth) #load package

#read in the table describing samples and kallisto output
#assign to variable name stab

stab = read.table('sleuth_table.tsv',header=TRUE)


#initialize sleuth object using sleuth_prep function from sleuth library
so = sleuth_prep(stab)

#fit model comparing 2dpi to 6dpi
so = sleuth_fit(so, ~condition, 'full')

#fit reduced model to fit likelihood ration test

so = sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test for differential expression between conditions
so = sleuth_lrt(so, 'reduced', 'full')

#sidenote: see the help page for any function in R

#load the dplyr package for data.frame filtering

library(dplyr)

sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <=0.05) |> dplyr::arrange(pval)
#had to change gval to qval because sleuth was throwing errorss!
#print top 10 transcripts
head(sleuth_significant, n=10)

#write FDR < 0.05 transcripts to file
write.table(sleuth_significant[, c('target_id', 'test_stat', 'pval', 'qval')], file='fdr05_results.txt', quote=FALSE, row.names=FALSE)
#[,1 4] will only write the first four columns to file, which are the ones i want
#had to change it to a specific order, so i just used c(column names i want in order)