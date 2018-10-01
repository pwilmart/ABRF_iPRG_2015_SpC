# ABRF_iPRG_2015_SPC

## Spectral Counting (SPC) Analysis with R and edgeR
### Phil Wilmarth, OHSU
### October 1, 2018
The data is from the ABRF iPRG 2015 study and is described in (Choi-2017). Six proteins were prepared in 4 different abundance mixes and spiked into a yeast cell lysate background. Each of the 4 different spike-in experiments were analyzed in triplicate on a Q-Exactive instrument. Each sample was analyzed using a single 2-hour LC run.

The RAW files were downloaded and analyzed with Comet and the PAW pipeline. Although the study was really designed for feature intensity analysis using Skyline, spectral counting was done by a few groups in the study cohort. Spectral counting was also done here to demonstrate some considerations needed for that type of data. Differential expression statistical analysis was done in R using edgeR.

Successful proteomics data analyses often need more than good training, domain knowledge, and intelligence. They need common sense and a flair for the practical. Preparing the data for analysis from the (all too often) horrendous summary files takes skill and practice. There may be data that needs to excluded on the basis of measurement limitations. This is conceptually very different from excluding a specific sample. The analysis presented here touches on some of these aspects, which can also have applicability in other types of proteomics studies.

> Choi, M., Eren-Dogu, Z.F., Colangelo, C., Cottrell, J., Hoopmann, M.R., Kapp, E.A., Kim, S., Lam, H., Neubert, T.A., Palmblad, M. and Phinney, B.S., 2017. ABRF Proteome Informatics Research Group (iPRG) 2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC–MS/MS Experiments. Journal of proteome research, 16(2), pp.945-957.
