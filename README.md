# Whole-genome sequencing of _Staphylococcus aureus_ clinical isolates

<h3>Software used in the pipeline.</h3>

**FastQC** is widely used and is robust, efficient, and versatile quality control software for a varied range of raw genetic data. It outputs a quality report which can be viewed to give an indication on how good the respective reads are. For more information, please visit the FastQC website.

**MultiQC** is used create a single report with interactive plots for multiple bioinformatics analyses across many samples. It reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control. For more information, please visit the MultiQC website.

**Trim galore** is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing) https://github.com/FelixKrueger/TrimGalore.

**Shovill** <br>
It is a pipeline for assembly of bacterial isolate genomes from Illumina paired-end reads. It uses SPAdes at its core, but alters the steps before and after the primary assembly step to get similar results in less time. Shovill also supports other assemblers like SKESA, Velvet and Megahit https://github.com/tseemann/shovill.

**QUAST** <br>
It is a quality assessment tool for evaluating and comparing genome assemblies by computing various metrics and works both with and without reference genomes. It produces many reports, summary tables and plots to help scientists in their research and in their publications. For more information, please visit the Quast website.

**Prokka** <br>
It is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files. It identifies features of interest in a set of genomic DNA sequences, and labelling them with useful information https://github.com/tseemann/prokka.

**Busco**

**MLST** Multilocus sequence typing (MLST) is an unambiguous procedure for characterising isolates of bacterial species using the sequences of internal fragments of (usually) seven house-keeping genes https://github.com/tseemann/mlst.

**spaTyper** is a computational method for finding spa types. The spa typing method is based on sequencing of the polymorphic X region of the protein A gene (spa), present in all strains of Staphylococcus aureus. https://github.com/HCGB-IGTP/spaTyper

**AgrVATE** is a tool for rapid identification of Staphylococcus aureus agr locus type and also reports possible variants in the agr operon https://github.com/VishnuRaghuram94/AgrVATE.

**Abricate** <br>
It is a tool for mass screening of contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB. For more information, please visit https://github.com/tseemann/abricate.

**ResFinder**  identifies acquired genes and/or finds chromosomal mutations mediating antimicrobial resistance in total or partial DNA sequence of bacteria https://bitbucket.org/genomicepidemiology/resfinder/src/master/.

**MOB-suite** is a software tools for clustering, reconstruction and typing of plasmids from draft assemblies (https://github.com/phac-nml/mob-suite).

**PhiSpy** is a tool that identifies prophages in Bacterial (and probably Archaeal) genomes https://github.com/linsalrob/PhiSpy. 


