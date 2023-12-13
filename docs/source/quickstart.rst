Quickstart
==========

Prerequisities
--------------
The preferred way of running mitopy is through Docker.
Please install `Docker <https://docs.docker.com/desktop/>`_ and `Docker Compose <https://docs.docker.com/compose/>`_ on your system.


Run mitopy
-----------

Clone the repository::

    git clone https://github.com/bendda/mitopy

Start the Docker container in interactive mode::

    cd mitopy
    docker compose run mitopy

Run the whole mitopy pipeline on example data::

    mitopy run-pipeline example_data/NA12878_20k_hg38.bam -o example_data/outputs


To run the analysis on your files, please add them to ``data/`` directory. 

.. note::
    It is important to add custom data required for analysis to ``data/`` directory, since the ``data/`` directory is mounted to Docker container (for further information, see `Docker Volumes <https://docs.docker.com/storage/volumes/>`_).
    
   

Run tests
----------

To run a suite of mitopy tests, run following command (in the cloned mitopy directory)::

    docker compose run test

