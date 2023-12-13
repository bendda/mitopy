# mitopy

mitopy is a command-line toolkit to perform identification and analysis of short mitochondrial variants. The mitochondrial variant detection is based on [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-). Functionality for analysis of detected variants includes annotation, visualization, haplogroup identification and coverage calculation. For futher information, please see [documentation](https://mitopy.readthedocs.io)


## Quickstart

### Prerequisities
The preferred way of running mitopy is through Docker. Please install Docker and Docker Compose on your system.

### Run mitopy

Clone the repository:

```
git clone https://github.com/bendda/mitopy
```

Start the Docker container in interactive mode:

```
cd mitopy
docker compose run mitopy
```

Run the whole mitopy pipeline on example data:

```
mitopy run-pipeline example_data/NA12878_20k_hg38.bam -o example_data/outputs
```

To run the analysis on your files, please add them to `data/` directory.

---
**NOTE**

It is important to add custom data required for analysis to `data/` directory, since the data/ directory is mounted to Docker container (for further information, see [Docker Volumes](https://docs.docker.com/storage/volumes/)).

---


### Run tests

To run a suite of mitopy tests, run following command (in the cloned mitopy directory):

```
docker compose run test
```