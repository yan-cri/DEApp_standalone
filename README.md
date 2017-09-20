# DE_tools
DE analysis function for edgeR, DESeq2, and Limma-Voom corresponding to the published manuscript at https://scfbm.biomedcentral.com/articles/10.1186/s13029-017-0063-4, the DEApp shiny application is accessible at https://gallery.shinyapps.io/DEApp/ or https://yanli.shinyapps.io/DEApp/

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

The full command to run the test sample is 'Rscript DE_analysis.R testData/pnas-count_singleFactor.txt testData/pnas-count_singleFactor-meta.txt Control DHT $PWD/DEG_analysis_results 2 1 True True 0.05 1.5 > DE_analysis.Rout.log 2> DE_analysis.Rerr.log'

The log files are saved under the current directory with name 'DE_analysis.Rout.log' and 'DE_analysis.Rerr.log', where 'DE_analysis.Rout.log' includes the results summary, and 'DE_analysis.Rerr.log' includes the error message if program runs incorrectly.

## Feedback

If you have any comments or suggestions for this tool, please contact Yan Li, Center for Research Informatics, the University of Chicago at e-mail yli22@bsd.uchicago.edu

