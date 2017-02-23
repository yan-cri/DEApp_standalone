# DE_tools
DE analysis function for edgeR, DESeq2, and Limma-Voom.

##How To Use It

Rscript /fullpathTo/DE_analysis_edgeR.R /fullpathTo/featureCount-mRNA-20362.txt /fullpathTo/metatable-fc-de.txt Control KO /fullpathTo/outputDir 2 3 True True /fullpathTo/DE-Analysis-AllComb.R > DE_analysis_edgeR.Rout1 2> DE_analysis_edgeR.Rerr1

Where, the following parameters are:

1. raw count input /fullpathTo/featureCount-mRNA-20362.txt

2. metadata table input /fullpathTo/metatable-fc-de.txt

3. group 1: Control

4. group 2: KO

5. output directory: /fullpathTo/outputDir

6. group cutoff: 2

7. cpm cutoff: 3

8. whether to save the results: True

9. whether to save the results2: True

10. full path to the source code: /fullpathTo/DE-Analysis-AllComb.R

The same for DESeq2 and Limma-Voom analysis.
