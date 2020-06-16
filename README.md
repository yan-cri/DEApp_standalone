# DE_tools
DE analysis function for edgeR, DESeq2, and Limma-Voom corresponding to the published manuscript at https://scfbm.biomedcentral.com/articles/10.1186/s13029-017-0063-4, the DEApp shiny application is accessible at https://gallery.shinyapps.io/DEApp/ or https://yanli.shinyapps.io/DEApp/

If you are going to use the scripts here for your DE analysis, please cite ["Yan Li and Jorge Andrade. DEApp: an interactive web interface for differential expression analysis of next generation sequence data. Source Code for Biology and Medicine, 12:2, February 2017."](https://scfbm.biomedcentral.com/articles/10.1186/s13029-017-0063-4).

## How To Use It

### 1. step1. download this repo to your local PC. Directions for downloading, please click top right green color 'Clone or download' button.

### 2. step2. go into the downloaded folder, you can see R script named as 'DE_analysis.R', which is used to run all 3 different DE analysis and generate overlapping Venn-diagram.

To run 'DE_analysis.R', you can use command 'Rscript DE_analysis.R' with additional 11 parameter flags defined as below: 

1. raw count input: /fullpathTo/rawCount.txt

2. metadata table input: /fullpathTo/metatable.txt

3. group 1 name: e.g. Control, consistent with metatable.txt group information

4. group 2 name: e.g. KO, consistent with metatable.txt group information

5. output directory: /fullpathTo/DEG_analysis_results/

6. group cutoff: 2

7. cpm cutoff: 3

8. logical, whether to save the results: True

9. logical, whether to use orignal p-value or FDR corrected p-value for DEGs identification: True

10. p-values / FRD adjusted p-value: 0.05

11. FC(Fold Change) cutoff: 1.5

The full command to run the test sample is 'Rscript DE_analysis.R $PWD/testData/pnas-count_singleFactor.txt $PWD/testData/pnas-count_singleFactor-meta.txt Control DHT $PWD/DEG_analysis_results 2 1 True True 0.05 1.5 > DE_analysis.Rout.log 2> DE_analysis.Rerr.log'

The log files are saved under the current directory with name 'DE_analysis.Rout.log' and 'DE_analysis.Rerr.log', where 'DE_analysis.Rout.log' includes the results summary, and 'DE_analysis.Rerr.log' includes the error message if program runs incorrectly.

## DE analysis results outputs

The entire DE analysis are saved inside direcotry named as 'DEG_analysis_results' (be aware to provide the full path where this results is saved at), where 5 folders and 1 txt file are saved inside this directory as following:

1. 'txt' file named as 'rmlow-logcpm-cpmCutoof_1-sampCutoff_2-17633.txt' is the logrithm transformed CPM (counts per million) expression value after low expression genes removal, here the file name could change due to different cpm and sample cutting off values usage.

2. 3 folders named as 'DEseq2-res-comp-XX_XX', 'edgeR-res-comp-XX_XX', and 'limma-voom-res-comp-XX_XX' have entire DE analysis results saved inside with respect to each DE analysis methods, here XX_XX indicateing group1 and group2 name for your own DE anlaysis. In each folder, there are 3 files saved insider, where '\*-full-GLMDisp\*' has entire DE analysis results, '\*-DEall\*-full.txt' and '\*-DEall\*-GeneName.txt' include the identified DE analyiss results and gene names respectively with corresponding specified FDR corrected p-value and FC thresholds.

3. folder named as 'ol-res' includes the overlapped venn-diagram of 3 different DE analysis methods and corresponding overlapped DE analysis results with respect to venn-diagram.

4. folder named as 'plot-res' includes the volcano plot of each DE analysis methods.

## Feedback

If you have any comments or suggestions for this tool, please contact Yan Li, Center for Research Informatics, the University of Chicago at e-mail yli22@bsd.uchicago.edu

