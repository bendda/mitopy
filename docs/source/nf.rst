Nextflow implementation
========================

:Source code: https://github.com/bendda/mitopy-nf

To run the mitopy pipeline at scale on multiple samples, we provide a Nextflow-based implementation of the pipeline, which can be seamlessly deployed to cloud-based environments or computing clusters.

By this implementation, we also demonstrate usage of mitopy as an underlying library for implementation of scalable pipelines in the workflow manager of choice. 

How to run
-----------

Prerequisities
***************
To run the pipeline, please install `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_ and `Docker <https://docs.docker.com/desktop/>`_ on your system.

Input parameters
*****************

.. list-table::
   :widths: 20 10 70
   :header-rows: 1
   :class: tight-table  

   * - Parameter
     - Default
     - Description
   * - ``--alignments``
     - null
     - Path to directory containing alignment files (``.bam``) and respective index files (``.bai``). Ex: ``path/to/bams/*.{bam,bai}``
   * - ``--mt-reference``
     - rCRS
     - Mitochondrial reference genome to use for analysis. By default, **rCRS** reference is used. **RSRS** reference is included as optional reference.
   * - ``--outdir``
     - ./outputs
     - Output directory

Run pipeline
*************

Run the pipeline on your data using following command (adjust alignments path)::

    # Run on bam alignment files located in TEST/ directory
    nextflow run bendda/mitopy-nf -r main -latest \
        --alignments 'TEST/*.{bam,bai}' \
        --outdir results


Outputs
********

.. code-block:: bash

    outputs/
    ├── alignments
    │   ├── sample.bam
    │   ├── sample.bam.bai
    │   ├── sample_shifted.bam
    │   ├── sample_shifted.bam.bai
    ├── annotation
    │   ├── sample_annotated.csv
    │   ├── sample_annotated.vcf
    ├── haplogroup_report
    │   ├── sample_haplogroup.txt
    ├── variant_calls
    │   ├── sample.vcf
    │   ├── sample.vcf.idx
    └── visualization
        ├── sample.html


