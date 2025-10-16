# TFM Trabajo de Fin de Máster (Master's Thesis Project)
Bioinformatics thesis project

Analysis in lung cancer (LUAD-TCGA)
“Mechanistic analysis of signaling pathways under in silico EGFR inhibition and their relationship with somatic mutations in lung adenocarcinoma using HiPathia and dN/dS”

This project explores how *in silico* inhibition of **EGFR (Epidermal Growth Factor Receptor)** affects 
cell signaling pathways in **lung adenocarcinoma**, and how these effects relate to **somatic mutations** 
in patient samples. The analysis integrates **pathway activity modeling (HiPathia)** with 
**selection analysis (dN/dS)** to identify driver mutations and mechanistic links between mutations 
and pathway dysregulation.

# Citation
Cite this specific workflow as María Granados Presa (2025). If published, cite will be updated in this README file.

I do not recognise authorship of the tools used, please cite tools separately:
- Hipathia: as specified by Bioconductor DOI: 10.18129/B9.bioc.hipathia https://www.bioconductor.org/packages/release/bioc/html/hipathia.html
- dndscv: https://github.com/im3sanger/dndscv or DOI: 10.1016/j.cell.2017.09.042  (Martincorena et al., 2017).

**Abstract**
This study aimed to evaluate the differential activity of signaling pathways in lung adenocarcinoma (LUAD) and their relationship with somatic mutations after in silico silencing of EGFR. 
Transcriptomic and mutational data from TCGA-LUAD was analyzed determining pathway activity with HiPathia, and assessing selective pressures on genes using the dNdScv model. EGFR knockout leads to strong activation of the FoxO pathway with sex-dependant differences, although significant in both (female and male). Suggesting it as a potential bypass or compensatory mechanism to EGFR inhibition should need further approaches. Conversely, the oxytocin signaling pathway, mediated by CDKN1A/p21, showed a significant and indirect 
inhibition, being alarming for disease outcomes. This result was reinforced by the identification of genes under positive selection, including EGFR, through dN/dS analysis. Notably, the L858R missense mutation provided an evolutive disadvantage, correlating directly with oxytocin CDKN1A pathway silencing and suggesting worse response to EGFR-targeted treatments. The integration of functional and evolutionary computational approaches demonstrated effectiveness in directly linking somatic mutations with signaling pathways. These results provide a computational framework complementary for personalized therapeutic strategies in LUAD.


---------------------------------------------------------------------------------
Attempt for AML consisted in (didn't choose this approach at the end but I leave in this repository relevant data for anyone interested):
First approach was using HiPathia to perform a mechanistic analysis on a cohort of more than 3000 RNAseq samples to study the effect of the inhibition of METTL3 in cellular pathways in AML.

Some opportunities/obstacles:

    -no paired data: AML is in blood so no healthy samples can be taken in parallel for the same patient, a cohort of AML patients and a cohort of healthy patients has to be taken separately. I used TARGET-AML cohort for AML and GTex for healthy.

    - METTL3 gene is not part of HiPathia so can't directly inhibitd it in silico and simulate the effect: my alternative is to use a new cohort were inhibition has been carried out in vivo https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94613 from study : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94613.



