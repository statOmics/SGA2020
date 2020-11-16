---
layout: default
title: Statistical Genomics Analysis 2020 (SGA2020)
---

![IntroFig](./pages/figs/IntroFig.png)

### Course Description
High throughput 'omics studies generate ever larger datasets and, as a consequence, complex data interpretation challenges. This course focuses on the statistical concepts involved in preprocessing, quantification and differential analysis of high throughput omics data. Moreover, more advanced experimental designs and blocking will also be introduced. The core focus will be on shotgun proteomics and next generation sequencing. The course will rely exclusively on free and userfriendly opensource tools in R/Bioconductor. The course will provide a solid basis for beginners, but will also bring new perspectives to those already familiar with standard data analysis workflows for proteomics and next-generation sequencing applications.

### Target Audience
This course is oriented towards biologists and bioinformaticians with a particular interest in differential analysis for quantitative 'omics.

### Prerequisites
The prerequisites for the Statistical Genomics course are the successful completion of a basic course of statistics that covers topics on data exploration and descriptive statistics, statistical modeling, and inference: linear models, confidence intervals, t-tests, F-tests, anova, chi-squared test.

The basis concepts may be revisited in my online course [https://gtpb.github.io/PSLS20/](https://gtpb.github.io/PSLS20/) and in [https://statomics.github.io/statistiekCursusNotas/](https://gtpb.github.io/PSLS20/)

A primer to R and Data visualisation  in R can be found in:

- R Basics: [https://dodona.ugent.be/nl/courses/335/](https://dodona.ugent.be/nl/courses/335/)
- R Data Exploration: [https://dodona.ugent.be/nl/courses/345/](https://dodona.ugent.be/nl/courses/345/)



---

#### Topics

**Introduction**

  - Slides: [Intro](assets/intro.pdf)
  - Software: [Install and Launch Statistical Software](pages/software4stats.md)
  - Recap linear models: [case study](assets/recapGeneralLinearModel.html)
  - Entire analysis for KPNA2 gene: [KPNA2](assets/08-multipleRegression_KPNA2.html)

**Part I: Quantitative proteomics**

  - [Download Tutorial Data](https://github.com/statOmics/SGA2019/tree/data)


  1. Bioinformatics for proteomics
  - Slides: [Bioinformatics for Proteomics](assets/martens_proteomics_bioinformatics_20190923.pdf)
  - Students can sharpen their background knowledge on Mass Spectrometry, Proteomics & Bioinformatics for Proteomics
 here:[Mass Spectrometry and Bioinformatics for Proteomics](pages/techVideos.md)

 2. Identification
 - Slides:  [False Discovery Rate and Target Decoy Approach](assets/1_Identification_Evaluation_Target_Decoy_Approach.pdf)
 - Tutorial: [Evaluating Target Decoy Quality](pages/Identification.md), [example script identification](assets/identification.html),
 [All searches](assets/identification_all.html)

 3. Preprocessing & Analysis of Label Free Quantitative Proteomics Experiments with Simple Designs
 - Install Software: [Installation instructions msqrob2](pages/installMsqrob2.md)
 - Slides: [Preprocessing](assets/2_MSqRob_data_analysisI.pdf)
 - Tutorial: [preprocessing](pages/sdaMsqrobSimple.md)


 4. Statistical Inference & Analysis of Experiments with Factorial Designs
 - Slides: [Inference](assets/2_MSqRob_data_analysisII.pdf)
 - Tutorial: [Statistical Data Analysis with MSqRob for Factorial Designs](pages/sdaMsqrobDesign.md)

 5. Reading Material and Technical details
    - Paper: [Sticker et al. (2020) Robust summarization and inference in proteome-wide label-free quantification](https://www.biorxiv.org/content/10.1101/668863v1)
    - part of PhD dissertation: [Extensive Background on proteomics and proteomics data analysis](assets/backgroundProteomicsDataAnalysis.pdf)
    - [Inference upon summarization](assets/technicalDetailsProteomics.html)


6. Stagewise testing: Omnibus test and post hoc analysis: [slides](assets/stagewiseTesting.pdf)

7. Solutions
  - [cptac median](assets/cptac_median.html)
  - [cptac robust](assets/cptac.html)
  - [cptac maxLFQ](assets/cptac_maxLfQ.html)
  - [cancer 3 vs 3](assets/cancer2_3x3.html)
  - [cancer 6 vs 6](assets/cancer2_6x6.html)
  - [cancer 9 vs 9](assets/cancer2_9x9.html)
  - [mouse CRD](assets/mouseCRD2.html)
  - [mouse RCB](assets/mouseRCB2.html)
  - [mouse RCB wrong analysis](assets/mouseRCBwrongAnalysis.html)
  - [heart](assets/heartMainInteraction.html)
  - [heart StageR](assets/heartMainInteractionStageR.html)


---

**Part II: Next-generation sequencing**

  - [Download Tutorial Data](https://github.com/statOmics/SGA2020/tree/data-rnaseq)

  1. Introduction to transcriptomics with next generation sequencing

      - slides: [intro](assets/rnaseq1.pdf)
      - Background: [RNA sequencing data hitchhiker's guide to expression analysis](https://www.zora.uzh.ch/id/eprint/181231/1/Ann_Rev_Biomed_Data_Science_-_RNA_sequencing_data__hitchhiker_s__guide_to_expression_analysis.pdf)
      - tutorial

        - Mapping: [html](assets/elegansMappingCountTable.html)
        - Differential Analysis: [html](assets/elegans.html), which source of variability is not included in the analysis and how could we account for this? Try to adjust the script accordingly.  

        - Background for the airway example (count table on small fastQ files available in the Tutorial Data):
      [html](assets/airwayMappingCountTable.html)

  2. More Complex Designs

      - Researchers assessed the effect of spinal nerve ligation (SNL) on the transcriptome of rats. In this experiment, transcriptome profiling occurred at two weeks and two months after treatment, for both the SNL group and a control group. Two biological replicates are used for every treatment - time combination. The researchers are interested in early and late effects and in genes for which the effect changes over time. The data can be downloaded from the ReCount project website (http://bowtie-bio.sourceforge.net/recount/, dataset Hammer et al.). The following code can be used to download an R/Bioconductor expression set object.

      ```
      file <- "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData"
      load(url(file))
      hammer.eset
      ```

      - Paired-end sequencing was performed on primary cultures from parathyroid tumors of 4 patients at 2 time points over 3 conditions (control, treatment with diarylpropionitrile (DPN) and treatment with 4-hydroxytamoxifen (OHT)). DPN is a selective estrogen receptor agonist and OHT is a selective estrogen receptor modulator. One sample (patient 4, 24 hours, control) was omitted by the paper authors due to low quality. Data, the count table and information on the experiment is available at http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37211. It is not required to do the read mapping!

3. Technical details on transcriptomics with next generation sequencing. Generalized linear models are introduced in the slides and bulk RNA-seq tools via their corresponding papers

    - [slides on GLM](assets/rnaseq2.pdf)
    - [Poisson GLM and parameter estimation](assets/poissonIRWLS-implemented.html)
    - [edgeR: Negative Binomial](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378882/)
    - [DESeq2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/)
    - [voom](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053721/)
    - [edgeR: Quasi Negative Binomial](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.362.9634&rep=rep1&type=pdf)

4. Solutions

    - Airway Example: [GenomeIndex](assets/airwayGenomeIndex.html), [read mapping and count table](assets/airwayMappingCountTableCorr.html)

##### [Instructors](pages/instructors.md)
