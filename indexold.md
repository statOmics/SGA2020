---
layout: default
title: Statistical Genomics Analysis 2019 (SGA2019)
---

![IntroFig](./pages/figs/IntroFig.png)

### Course Description
High throughput 'omics studies generate ever larger datasets and, as a consequence, complex data interpretation challenges. This course focuses on the statistical concepts involved in preprocessing, quantification and differential analysis of high throughput omics data. Moreover, more advanced experimental designs and blocking will also be introduced. The core focus will be on shotgun proteomics and next generation sequencing. The course will rely exclusively on free and userfriendly opensource tools in R/Bioconductor. The course will provide a solid basis for beginners, but will also bring new perspectives to those already familiar with standard data analysis workflows for proteomics and next-generation sequencing applications.

### Target Audience
This course is oriented towards biologists and bioinformaticians with a particular interest in differential analysis for quantitative 'omics.

### Prerequisites
The prerequisites for the Statistical Genomics course are the successful completion of a basic course of statistics that covers topics on data exploration and descriptive statistics, statistical modeling, and inference: linear models, confidence intervals, t-tests, F-tests, anova, chi-squared test.

The basis concepts may be revisited in the free e-book Practical Regression and Anova using R of J. Faraway. The book and additional material is freely available on [http://www.maths.bath.ac.uk/~jjf23/book/](http://www.maths.bath.ac.uk/~jjf23/book/).

- Brief introduction to R: appendix C
- Linear models: Chapter 1-3, 7, 8.1-8.2, 12
- Anova: Chapter 16.1-16.2


---

#### Topics

**Introduction**

  - Slides: [Intro](assets/intro.pdf)
  - Software: [Install and Launch Statistical Software](pages/software4stats.md)

**Part I: Quantitative proteomics**

  - [Download Tutorial Data](https://github.com/statOmics/SGA2019/tree/data)


  1. Bioinformatics for proteomics
  - Slides: [Bioinformatics for Proteomics](assets/martens_proteomics_bioinformatics_20190923.pdf)
  - Students can sharpen their background knowledge on Mass Spectrometry, Proteomics & Bioinformatics for Proteomics
 here:[Mass Spectrometry and Bioinformatics for Proteomics](pages/techVideos.md)

  2. Identification
  - Slides:  [False Discovery Rate and Target Decoy Approach](assets/1_Identification_Evaluation_Target_Decoy_Approach.pdf)
  - Tutorial: [Evaluating Target Decoy Quality](pages/Identification.md)  

  3. Preprocessing & Analysis of Label Free Quantitative Proteomics Experiments with Simple Designs
  - Slides: [Preprocessing](assets/2_MSqRob_data_analysisI.pdf)
  - Tutorial: [preprocessing](pages/sdaMsqrobSimple.md)

  4. Statistical Inference & Analysis of Experiments with Factorial Designs
  - Slides: [Inference](assets/2_MSqRob_data_analysisII.pdf)
  - Tutorial: [Statistical Data Analysis with MSqRob for Factorial Designs](pages/sdaMsqrobDesign.md)

  5. Technical details of linear models
  - Slides: [Inference](assets/parameterEstimationInferenceLinearModels.pdf)
  - Example: [html](pages/robustRegression.nb.html), [rmd](https://raw.githubusercontent.com/statOmics/SGA2019/master/rmarkdownExamples/robustRegression.Rmd)

  6. Stagewise testing: Omnibus test and post hoc analysis: [slides](assets/stagewiseTesting.pdf)

  7. Homework 2: Analysis of the heart example from the tutorial page in 4. Do the analysis with MSqRob. For one protein we will do the analysis with the functions lm, rlm and matrix algebra. Use the rmarkdownfile below as a template.
  - Use as name: Namegroupmember1Namegroupmember2Namegroupmember3_SGA2019_Homework2.Rmd
  - [Homework2.Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/rmarkdownExamples/Homework2.Rmd)
  - The homework is due by Tuesday 12/11/2019.

  8. Solutions

  - Identification Uniprot:[html](assets/assessDecoysUniprot.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/assessDecoysUniprot.Rmd)
  - Identification Swissprot:[html](assets/assessDecoysSwissprot.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/assessDecoysSwissprot.Rmd)
  - francisella:[html](assets/francisella.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/francisella.Rmd)
  - cancer3x3: [html](assets/cancer_3x3.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/cancer_3x3.Rmd)
  - cancer3x3: [html](assets/cancer_9x9.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/cancer_9x9.Rmd)
  - mouse CRD: [html](assets/mouseCRD.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/mouseCRD.Rmd)
  - mouse RCB: [html](assets/mouseRCB.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/mouseRCB.Rmd)
  - mouse RCB wrong analysis (no blocking): [html](assets/mouseRCB_wrong_analysis.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/mouseRCB_wrong_analysis.Rmd)

**Part II: Next-generation sequencing**

  - [Download Tutorial Data](https://github.com/statOmics/SGA2019/tree/data-rnaseq)

  1. Introduction to transcriptomics with next generation sequencing

      - slides: [intro](assets/rnaseq1.pdf)
      - tutorial

        - Mapping: [html](pages/elegansMappingCountTable.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/elegansMappingCountTable.Rmd)
        - Differential Analysis: [html](pages/elegans.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/elegans.Rmd), which source of variability is not included in the analysis and how could we account for this? Try to adjust the script accordingly.  
        - Background for the airway example (count table on small fastQ files available in the Tutorial Data to be prepared by Monday November 4th 2019):
      [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/airwayMappingCountTable.Rmd)
        - Airway entire analysis: genome index [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/airwayGenomeIndex.Rmd), [html](assets/airwayGenomeIndex.html), read mapping and count table [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/airwayMappingCountTableCorr.Rmd),[html](assets/airwayMappingCountTableCorr.html), DE analysis [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/airway.Rmd), [html](assets/airway.html)

  2. More Complex Designs

      - Researchers assessed the effect of spinal nerve ligation (SNL) on the transcriptome of rats. In this experiment, transcriptome profiling occurred at two weeks and two months after treatment, for both the SNL group and a control group. Two biological replicates are used for every treatment - time combination. The researchers are interested in early and late effects and in genes for which the effect changes over time. The data can be downloaded from the ReCount project website (http://bowtie-bio.sourceforge.net/recount/, dataset Hammer et al.). The following code can be used to download an R/Bioconductor expression set object.

      ```
      file <- "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData"
      load(url(file))
      hammer.eset
      ```

      - Paired-end sequencing was performed on primary cultures from parathyroid tumors of 4 patients at 2 time points over 3 conditions (control, treatment with diarylpropionitrile (DPN) and treatment with 4-hydroxytamoxifen (OHT)). DPN is a selective estrogen receptor agonist and OHT is a selective estrogen receptor modulator. One sample (patient 4, 24 hours, control) was omitted by the paper authors due to low quality. Data, the count table and information on the experiment is available at http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37211. It is not required to do the read mapping!

      - Submit your RMD scripts via Ufora by Wednesday November 20th 2019  

  3. Technical details on transcriptomics with next generation sequencing. Generalized linear models are introduced in the slides and bulk RNA-seq tools via their corresponding papers

      - [slides on GLM](assets/rnaseq2.pdf)
      - Poisson GLM and parameter estimation: [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/poissonIRWLS-implemented.Rmd), [html](assets/poissonIRWLS-implemented.html)
      - [edgeR: Negative Binomial](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378882/)
      - [DESeq2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/)
      - [voom](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053721/)
      - [edgeR: Quasi Negative Binomial](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.362.9634&rep=rep1&type=pdf)

  4. DE analysis starting from transcript level counts

      - [Soneson et al. 2016](https://f1000research.com/articles/4-1521/v2)
      - Building index with salmon
      ```
      salmon index --gencode -t gencode.v32.transcripts.fa -i gencode.v32_salmon_index
      ```
      - Mapping one sample with salmon:
      ```
salmon quant -i gencode.v32_salmon_index -l A --gcBias -1 SRR1039508_subset_1.fastq -2 SRR1039508_subset_2.fastq --validateMappings -o quant/SRR1039508_subset_quant
```
      - [Intro to salmon](https://combine-lab.github.io/salmon/getting_started/)
      - airway with DESeq2: [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/airway_salmon_DESeq2.Rmd);[html](assets/airway_salmon_DESeq2.html)
      - airway with EdgeR: [Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/airway_salmon_edgeR.Rmd);[html](assets/airway_salmon_edgeR.html)


4. Solutions

    - hammer edgeR: [html](assets/hammer.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/hammer.Rmd)
    - hammer alternative coding: [html](assets/hammerAlternativeCoding.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/hammerAlternativeCoding.Rmd)
    - parathyroid (few comments in file, template for discussions): [html](assets/parathyroid.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/parathyroid.Rmd)
    - hammer QL: [html](assets/hammer_quasi.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/hammer_quasi.Rmd)
    - hammer DESeq2:
    [html](assets/hammer_DESEQ2.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/hammer_DESEQ2.Rmd)
    - hammer QL + stagewiseTesting:
    [html](assets/hammer_quasi_stagewise.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/hammer_quasi_stagewise.Rmd)
    - parathyroid Stagewise Testing(few comments in file, template for discussions): [html](assets/parathyroid_stagewise.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/parathyroid_stagewise.Rmd)
    - parathyroid Stagewise Testing with average contrast (few comments in file, template for discussions): [html](assets/parathyroid_stagewise_avgContrast.html),[Rmd](https://raw.githubusercontent.com/statOmics/SGA2019/gh-pages/assets/parathyroid_stagewise_avgContrast.Rmd)

5. Single Cell analysis
    - [Slides](assets/scRNA-seq.pdf)
    - Orchestrating Single Cell: [workshopVignette](http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/OSCABioc2019__OSCABioc2019/), [ebook](https://osca.bioconductor.org)
    - Analysis of multi-sample
multi-group scRNA-seq data:  [pseudoBulk](http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/muscWorkshop__vignette/)
    - Project work: See [ufora.ugent.be](ufora.ugent.be)





---

##### [Instructors](pages/instructors.md)
