# TFM
Bioinformatics thesis project

Analysis in lung cancer (LUAD-TCGA)
“Mechanistic analysis of signaling pathways under in silico EGFR inhibition and their relationship with somatic mutations in lung adenocarcinoma using HiPathia and dN/dS”

This project explores how *in silico* inhibition of **EGFR (Epidermal Growth Factor Receptor)** affects 
cell signaling pathways in **lung adenocarcinoma**, and how these effects relate to **somatic mutations** 
in patient samples. The analysis integrates **pathway activity modeling (HiPathia)** with 
**selection analysis (dN/dS)** to identify driver mutations and mechanistic links between mutations 
and pathway dysregulation.


**Abstract**

This study aimed to evaluate the differential activity of signaling pathways in lung adenocarcinoma (LUAD) and their relationship with somatic mutations after in silico silencing of EGFR. Transcriptomic and mutational data from TCGA-LUAD was analyzed determining pathway activity with HiPathia, and assessing selective pressures on genes using the dNdScv model. The FoxO signaling pathway emerged as a potential bypass to EGFR inhibition, suggesting its role in therapeutic resistance. Conversely, the oxytocin signaling pathway, mediated by CDKN1A/p21, showed significant and beneficial activity for disease outcomes. This result was reinforced by the identification of genes under positive selection, including EGFR, through dN/dS analysis. Notably, the L858R missense mutation provided an adaptive advantage, correlating directly with oxytocin CDKN1A pathway activation and suggesting improved response to EGFR silencing. The integration of functional and evolutionary computational approaches demonstrated effectiveness in directly linking somatic mutations with signaling pathways. These results highlight novel mechanisms of resistance and open new perspectives for the development of personalized therapeutic strategies in LUAD.


Attempt for AML consisted in (didn't choose this approach at the end but I leave in this repository relevant data for anyone interested):
First approach was using HiPathia to perform a mechanistic analysis on a cohort of more than 3000 RNAseq samples to study the effect of the inhibition of METTL3 in cellular pathways in AML.

Some opportunities/obstacles:

    -no paired data: AML is in blood so no healthy samples can be taken in parallel for the same patient, a cohort of AML patients and a cohort of healthy patients has to be taken separately. I used TARGET-AML cohort for AML and GTex for healthy.

    - METTL3 gene is not part of HiPathia so can't directly inhibitd it in silico and simulate the effect: my alternative is to use a new cohort were inhibition has been carried out in vivo https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94613 from study : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94613.



