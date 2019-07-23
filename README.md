# ABRF_iPRG_2015_SpC

## Spectral Counting (SpC) Analysis with R and edgeR
### Phil Wilmarth, OHSU
### October 1, 2018

The data is from the ABRF iPRG 2015 study and is described in this publication:

> Choi, M., Eren-Dogu, Z.F., Colangelo, C., Cottrell, J., Hoopmann, M.R., Kapp, E.A., Kim, S., Lam, H., Neubert, T.A., Palmblad, M. and Phinney, B.S., 2017. ABRF Proteome Informatics Research Group (iPRG) 2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC–MS/MS Experiments. Journal of proteome research, 16(2), pp.945-957.

> **Abstract:**
Detection of differentially abundant proteins in label-free quantitative shotgun liquid chromatography-tandem mass spectrometry (LC-MS/MS) experiments requires a series of computational steps that identify and quantify LC-MS features. It also requires statistical analyses that distinguish systematic changes in abundance between conditions from artifacts of biological and technical variation. The 2015 study of the Proteome Informatics Research Group (iPRG) of the Association of Biomolecular Resource Facilities (ABRF) aimed to evaluate the effects of the statistical analysis on the accuracy of the results. The study used LC-tandem mass spectra acquired from a controlled mixture, and made the data available to anonymous volunteer participants. The participants used methods of their choice to detect differentially abundant proteins, estimate the associated fold changes, and characterize the uncertainty of the results. The study found that multiple strategies (including the use of spectral counts versus peak intensities, and various software tools) could lead to accurate results, and that the performance was primarily determined by the analysts' expertise. This manuscript summarizes the outcome of the study, and provides representative examples of good computational and statistical practice. The data set generated as part of this study is publicly available.

Six proteins were prepared in 4 different abundance mixes and spiked into a yeast cell lysate background. Each of the 4 different spike-in experiments were analyzed in triplicate on a Q-Exactive instrument. Each sample was analyzed using a single 2-hour LC run.

The RAW files were downloaded and analyzed with Comet and the PAW pipeline. Although the study was really designed for feature intensity analysis using Skyline, spectral counting was done by a few groups in the study cohort. Spectral counting was also done here to demonstrate some considerations needed for that type of data. Differential expression statistical analysis was done in R using edgeR.

Successful proteomics data analyses often need more than good training, domain knowledge, and intelligence. They need common sense and a flair for the practical. Preparing the data for analysis from the (all too often) horrendous summary files takes skill and practice. There may be data that needs to excluded on the basis of measurement limitations. This is conceptually very different from excluding a specific sample. The analysis presented here touches on some of these aspects, which can also have applicability in other types of proteomics studies.

A direct link to the rendered notebook is [here](https://pwilmart.github.io/TMT_analysis_examples/ABRF_2015_edgeR.html)

---

**_File Key:_**
* ABRF_2015_edgeR.ipynb - _main Jupyter notebook_
* ABRF_2015_grouped_protein_summary_8.xlsx - _spreadsheet of grouped protein results_
* ABRF_2015_pipeline.log - _consolated log file of all PAW steps_
* JD_sample1-A_peptide_results_8.txt - _first of 12 detailed PSM files_
* ...
* JD_sample4_C_peptide_results_8.txt - _last of 12 detailed PSM files_
* PAW_grouped_proteins_with_stats.txt - _output file from notebook_
* PAW_protein_grouper.log - _log file from protein grouping step_
* PAW_results.log - _log file from protein inference step_
* README.md - _this file_
* extras_iPRG2015_both.fasta - _FASTA database used in Comet searches_
* grouped_peptide_summary_8.txt - _peptide summary for grouped proteins_
* grouped_protein_summary_8.txt - _protein summary after protein grouping_
* peptide_summary_8.txt - _peptide summary for inferred proteins_
* protein_summary_8.txt - _list of inferred proteins (group members are explicit)_
