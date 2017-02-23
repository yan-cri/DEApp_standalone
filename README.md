# DE_tools
DE analysis function for edgeR, DESeq2, and Limma-Voom.

##How To Use It

###To run DE_analysis.R

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

11. p-values / FRD adjusted p-value: 0.05

12. FC(Fold Change) cutoff: 1.5

The log output - DE_analysis.Rout.log include the results summary, and log output DE_analysis.Rerr.log includes the error message if program runs incorrectly.

###To run DE_analysis_visualization.R

Rscript DE_analysis_visualization.R /fullpathTo/DE_analysis_vis_Fn.R /fullpathTo/DEG_analysis_res group1 group2 2 3 0.05 1.5

Where, the following parameters are:

1. function code: /fullpathTo/DE_analysis_vis_Fn.R

2. DE analysis results directory in above running: /fullpathTo/DEG_analysis_res

3. group 1 name: e.g. Control, consistent with above running 

4. group 2 name: e.g. KO, consistent with above running 

5. group cutoff: 2

6. cpm cutoff: 3

7. p-values / FRD adjusted p-value: 0.05

8. FC(Fold Change) cutoff: 1.5





