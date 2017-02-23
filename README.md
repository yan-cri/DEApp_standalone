# DE_tools
DE analysis function for edgeR, DESeq2, and Limma-Voom.

##How To Use It

Rscript DE_analysis.R /fullpathTo/rawCount.txt /fullpathTo/metatable.txt Control KO /fullpathTo/DEG_analysis_res/ 2 3 True True /fullpathTo/DE-Analysis-AllComb-Fn.R 0.05 1.5 > /fullpathTo/DE_analysis.Rout.log 2> /fullpathTo/DE_analysis.Rerr.log

Where, the following parameters are:

1. raw count input: /fullpathTo/rawCount.txt

2. metadata table input: /fullpathTo/metatable.txt

3. group 1 name: e.g. Control, consistent with metatable.txt group information

4. group 2 name: e.g. KO, consistent with metatable.txt group information

5. output directory: /fullpathTo/outputDir

6. group cutoff: 2

7. cpm cutoff: 3

8. logical, whether to save the results: True

9. logical, whether to use orignal p-value or FDR corrected p-value for DEGs identification: True

10. full path to the source code: /fullpathTo/DE-Analysis-AllComb-Fn.R

The log output - DE_analysis.Rout.log include the results summary, and log output DE_analysis.Rerr.log includes the error message if program runs incorrectly.
