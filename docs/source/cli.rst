Command-line interface
=======================


``run-pipeline``
----------------
Run the whole single-sample pipeline on the input WGS alignment file in BAM or CRAM format::

    mitopy run-pipeline [OPTIONS] BAM
    

.. note::
  The input alignment file should be coordinate-sorted and indexed, however if these prerequisities are not met, mitopy will coordinate-sort and index the input file.

.. list-table::
   :widths: 20 10 70
   :header-rows: 1
   :class: tight-table 

   * - Option
     - Deafult
     - Description
   * - ``--bai``
     - null
     - Alignment index file (BAI/CRAI). If not provided, it is assumed it resides in the same directory as input BAM.
   * - ``--mt-ref``
     - rCRS
     - Mitochondrial reference. By default, variants are detected and analyzed with respect to **rCRS** reference. We include **RSRS** as an optional mitochondrial reference.
   * - ``--reference-fa``
     - null
     - Reference genome FASTA file. **Only required when input is a CRAM file**.
   * - ``--contig-name``
     - null
     - Name of the mitochondrial contig in the alignment file. If not provided, it will be automatically detected.
   * - ``--out-dir`` ``-o``
     - BAM_DIR
     - Output directory. By default, results are outputed in the directory of input BAM file.
   * - ``--prefix`` ``-p``
     - BAM_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 
   * - ``--ncores`` ``-c``
     - 1
     - Number of cores.
   * - ``--tmp-dir``
     - tmp
     - Directory for intermediate files.
   * - ``--remove-tmp``
     - false
     - If true, remove intermediate files after the analysis is completed.
   * - ``--m2-extra-args``
     - ""
     - Extra arguments to pass onto Mutect2 variant caller.

Variant postprocessing options

.. list-table::
   :widths: 20 10 70
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--f-score-beta``
     - 1.0
     - F score beta. Specifies the relative weight of recall and precision for the filtering strategy.
   * - ``--vaf-treshold``
     - 0.0
     - Minimum variant allele fraction treshold. All sites with variant allele fraction below the treshold will be filtered.
   * - ``--blacklisted-sites``
     - null
     - Custom BED file containing blacklisted sites (the BED index file has to be present as well). If not specified, the default blacklist for chosen MT reference will be used. 
   * - ``--autosomal-coverage``
     - 0.0
     - Median autosomal coverage. Set to activate filter against erroneously mapped nuclear mitochondrial DNA segments (NuMTs). To estimate median autosomal coverage from WGS BAM, `Picard CollectWgsMetrics <https://gatk.broadinstitute.org/hc/en-us/articles/360036804671-CollectWgsMetrics-Picard->`_ can be used.
   * - ``--contamination-filter``
     - false
     - Contamination filter. If enabled, sample contamination level will be estimated using `haplocheck <https://mitoverse.readthedocs.io/haplocheck/haplocheck/>`_ and variants will be filtered (valid only for rCRS mitochondrial reference).
   * - ``--remove-non-pass``
     - true
     - Remove variants not passing the enabled filters from final VCF file.
   * - ``--normalize``
     - true
     - Split multi-allelic sites and left-align variant calls. 

Annotation options

.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--min-hom-treshold``
     - 0.95
     - Minimum homoplasmy level treshold. Annotate variants above this treshold as homoplasmic, otherwise heteroplasmic.
   * - ``--population-freqs``
     - true
     - Annotate variants with population frequencies from `gnomAD <https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna>`_ database.
   * - ``--conservation-scores``
     - true
     - Include conservation scores from `PhyloP100way <https://genome.ucsc.edu/cgi-bin/hgc?hgsid=784677241_vYLABfJrjxNKeDTusOROCSUBXtnK&c=chrM&l=0&r=16569&o=0&t=16569&g=phyloP100way&i=phyloP100way>`_ and `PhastCons100way <https://genome.ucsc.edu/cgi-bin/hgc?hgsid=916826631_g8XasCQqrg8t9dxczEQmzhNA9Nyc&c=chr12&l=53858048&r=53859044&o=53858048&t=53859044&g=phastCons100way&i=phastCons100way>`_ in annotations.
   * - ``--patho-predictions``
     - true
     - Annotate variants with in-silico pathogenicity predictions from `SIFT <https://sift.bii.a-star.edu.sg/>`_, `MitoTIP <https://www.mitomap.org/MITOMAP/MitoTipInfo>`_ and `PON-mt-tRNA <http://structure.bmc.lu.se/PON-mt-tRNA/datasets.html/>`_.
   * - ``--phenotype-annot``
     - true
     - Annotate variants with phenotype information from `MITOMAP <https://www.mitomap.org/MITOMAP>`_ and `ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`_ databases.
   * - ``--create-annotation-report``
     - true
     - Export annotated variants to human-readable CSV format.

Visualization options

.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--split-strands``
     - false
     -  Split H and L strand of mitochondrial genome in the visualization.
   * - ``--save-as-png``
     - false
     - Additionally, save plot as static PNG image.



``preprocess``
--------------
Prepare input WGS alignment file in BAM or CRAM format for the mitochondrial analysis::

    mitopy preprocess [OPTIONS] BAM
    

.. list-table::
   :widths: 20 10 70
   :header-rows: 1
   :class: tight-table 

   * - Option
     - Default
     - Description
   * - ``--bai``
     - null
     - Alignment index file (BAI/CRAI). If not provided, it is assumed it resides in the same directory as input BAM.
   * - ``--reference-fa``
     - null
     - Reference genome FASTA file. **Only required when input is a CRAM file**.
   * - ``--contig-name``
     - null
     - Name of the mitochondrial contig in the alignment file. If not provided, it will be automatically detected.
   * - ``--out-dir`` ``-o``
     - BAM_DIR
     - Output directory. By default, results are outputed in the directory of input BAM file.
   * - ``--prefix`` ``-p``
     - BAM_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 


``align``
----------

Align mitochondrial reads in unmapped BAM format  (`uBAM <https://gatk.broadinstitute.org/hc/en-us/articles/360035532132-uBAM-Unmapped-BAM-Format>`_) to canonical or shifted mitochondrial reference using `bwa-mem2 <https://github.com/bwa-mem2/bwa-mem2>`_::

    mitopy align [OPTIONS] UBAM

.. list-table::
   :widths: 25 5 70
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--mt-ref``
     - rCRS
     - Mitochondrial reference to align against. By default, the reads are aligned to **rCRS** mitochondrial reference. We include **RSRS** as an optional mitochondrial reference.
   * - ``--shifted``
     - false
     - Shifted mode. If enabled, the alignment is performed against shifted mitochondrial reference.
   * - ``--ncores`` ``-c``
     - 1
     - Number of cores.
   * - ``--out-dir`` ``-o``
     - UBAM_DIR
     - Output directory. By default, results are outputed in the directory of input UBAM file.
   * - ``--prefix`` ``-p``
     - UBAM_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 


``call`` 
---------
Call variants in non-control region (using canonical mitochondrial reference) or control region (using shifted mitochondrial reference) of mitochondrial genome using `Mutect2 <https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2>`_::

    mitopy call [OPTIONS] BAM


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--mt-ref``
     - rCRS
     - Mitochondrial reference. By default, variants are called against **rCRS** mitochondrial reference. We include **RSRS** as optional mitochondrial reference.
   * - ``--m2-extra-args``
     - ""
     - Extra arguments to pass onto Mutect2 variant caller.
   * - ``--shifted``
     - false
     - Shifted mode. If enabled, the variant are called against shifted mitochondrial reference (**control region**).
   * - ``--out-dir`` ``-o``
     - BAM_DIR
     - Output directory. By default, results are outputed in the directory of input BAM file.
   * - ``--prefix`` ``-p``
     - BAM_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 



``merge``
---------

Merge variant calls and stats from control and non-control region::

    mitopy merge [OPTIONS] VCF VCF_SHIFTED


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--stats``
     - null
     - File containing Mutect2 stats for VCF file. By default, it is assumed that stats file is located in the same directory as input VCF.
   * - ``--stats-shifted``
     - null
     - File containing Mutect2 stats for shifted VCF file. By default, it is assumed that shifted stats file is located in the same directory as input VCF_SHIFTED.
   * - ``--mt-ref``
     - rCRS
     - Mitochondrial reference used in variant calling process. By default, it is assumed that variants were called against **rCRS**.
   * - ``--out-dir`` ``-o``
     - VCF_DIR
     - Output directory. By default, results are outputed in the directory of input VCF file.
   * - ``--prefix`` ``-p``
     - VCF_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 


``postprocess`` 
----------------

Postprocess raw variant calls to remove potential false-positives by applying several filters and normalize the VCF file::

    mitopy postprocess [OPTIONS] VCF


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--stats``
     - null
     - File containing Mutect2 stats for VCF file. By default, it is assumed that stats file is located in the same directory as input VCF.
   * - ``--mt-ref``
     - rCRS
     - Mitochondrial reference used in variant calling process. By default, it is assumed that variants were called against **rCRS**.
   * - ``--f-score-beta``
     - 1.0
     - F score beta. Specifies the relative weight of recall and precision for the filtering strategy.
   * - ``--vaf-treshold``
     - 0.0
     - Minimum variant allele fraction treshold. All sites with variant allele fraction below the treshold will be filtered.
   * - ``--blacklisted-sites``
     - null
     - Custom BED file containing blacklisted sites (the BED index file has to be present as well). If not specified, the default blacklist for chosen MT reference will be used. 
   * - ``--autosomal-coverage``
     - 0.0
     - Median autosomal coverage. Set to activate filter against erroneously mapped nuclear mitochondrial DNA segments (NuMTs). To estimate median autosomal coverage from WGS BAM, `Picard CollectWgsMetrics <https://gatk.broadinstitute.org/hc/en-us/articles/360036804671-CollectWgsMetrics-Picard->`_ can be used.
   * - ``--contamination-filter``
     - false
     - Contamination filter. If enabled, sample contamination level will be estimated using `haplocheck <https://mitoverse.readthedocs.io/haplocheck/haplocheck/>`_ and variants will be filtered (valid only for rCRS mitochondrial reference).
   * - ``--remove-non-pass``
     - true
     - Remove variants not passing the enabled filters from final VCF file.
   * - ``--normalize``
     - true
     - Split multi-allelic sites and left-align variant calls. 
   * - ``--out-dir`` ``-o``
     - VCF_DIR
     - Output directory. By default, results are outputed in the directory of input VCF file.
   * - ``--prefix`` ``-p``
     - VCF_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 



``coverage`` 
------------

Calculate per-base coverage using `mosdepth <https://github.com/brentp/mosdepth>`_ . Coverage is combined from control (``SHIFTED_MT_BAM``) and non-control region (``MT_BAM``)::

    mitopy coverage [OPTIONS] MT_BAM SHIFTED_MT_BAM


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--mt-bai``
     - null
     - Index file for input BAM containing reads aligned against mitochondrial reference. By default, it is assumed index file is located in the same directory as MT_BAM.
   * - ``--shifted-mt-bai``
     - null
     - Index file for shifted input BAM containing reads aligned against shifted mitochondrial reference. By default, it is assumed index file is located in the same directory as SHIFTED_MT_BAM.
   * - ``--create-plot``
     - true
     - Create coverage plot.
   * - ``--out-dir`` ``-o``
     - BAM_DIR
     - Output directory. By default, results are outputed in the directory of input BAM file.
   * - ``--prefix`` ``-p``
     - BAM_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 


``annotate``
------------

Annotate mitochondrial variant calls with functional effects using `SnpEff <http://pcingola.github.io/SnpEff/snpeff/introduction/>`_ and optionally add other annotations::

    mitopy annotate [OPTIONS] VCF


.. note::
  The annotation stage assumes that the input VCF file is normalized (specifically multi-allelic sites are split). 
  
  The variants are annotated with respect to **rCRS** mitochondrial reference.


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--min-hom-treshold``
     - 0.95
     - Minimum homoplasmy level treshold. Annotate variants above this treshold as homoplasmic, otherwise heteroplasmic.
   * - ``--population-freqs``
     - true
     - Annotate variants with population frequencies from `gnomAD <https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna>`_ database.
   * - ``--conservation-scores``
     - true
     - Include conservation scores from `PhyloP100way <https://genome.ucsc.edu/cgi-bin/hgc?hgsid=784677241_vYLABfJrjxNKeDTusOROCSUBXtnK&c=chrM&l=0&r=16569&o=0&t=16569&g=phyloP100way&i=phyloP100way>`_ and `PhastCons100way <https://genome.ucsc.edu/cgi-bin/hgc?hgsid=916826631_g8XasCQqrg8t9dxczEQmzhNA9Nyc&c=chr12&l=53858048&r=53859044&o=53858048&t=53859044&g=phastCons100way&i=phastCons100way>`_ in annotations.
   * - ``--patho-predictions``
     - true
     - Annotate variants with in-silico pathogenicity predictions from `SIFT <https://sift.bii.a-star.edu.sg/>`_, `MitoTIP <https://www.mitomap.org/MITOMAP/MitoTipInfo>`_ and `PON-mt-tRNA <http://structure.bmc.lu.se/PON-mt-tRNA/datasets.html/>`_.
   * - ``--phenotype-annot``
     - true
     - Annotate variants with phenotype information from `MITOMAP <https://www.mitomap.org/MITOMAP>`_ and `ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`_ databases.
   * - ``--create-csv``
     - true
     - Export annotated variants to human-readable CSV format.
   * - ``--out-dir`` ``-o``
     - VCF_DIR
     - Output directory. By default, results are outputed in the directory of input VCF file.
   * - ``--prefix`` ``-p``
     - VCF_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 



``visualize``
-------------

Visualize variant calls::

    mitopy visualize [OPTIONS] VCF


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--coverage-csv``
     - null
     - CSV file with calculated per-base coverage. If provided, coverage will be included in the final visualization.
   * - ``--split-strands``
     - false
     -  Split H and L strand of mitochondrial genome in the visualization.
   * - ``--save-as-png``
     - false
     - Save plot as static PNG image.
   * - ``--out-dir`` ``-o``
     - VCF_DIR
     - Output directory. By default, results are outputed in the directory of input VCF file.
   * - ``--prefix`` ``-p``
     - VCF_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.


``identify-haplogroup``
-----------------------

Identify sample haplogroup using `haplogrep3 <https://haplogrep.readthedocs.io/en/latest/>`_::

    mitopy identify-haplogroup [OPTIONS] VCF


.. list-table::
   :widths: 25 10 65
   :header-rows: 1
   :class: tight-table  

   * - Option
     - Default
     - Description
   * - ``--mt-ref``
     - rCRS
     - Mitochondrial reference. By default, haplogroup is classified with respect to **rCRS** reference. We include **RSRS** as an optional mitochondrial reference.
   * - ``--out-dir`` ``-o``
     - VCF_DIR
     - Output directory. By default, results are outputed in the directory of input VCF file.
   * - ``--prefix`` ``-p``
     - VCF_BASENAME
     - Prefix for output files. By default, resulting files will be prefixed with the input's file basename.
   * - ``--verbose`` ``-v``
     - false
     - Verbosity. If true, logs generated by underlying tools will be recorded. 














