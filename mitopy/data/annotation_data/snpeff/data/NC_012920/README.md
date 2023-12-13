Had to change chr name to ChrM

add this to config file:

NC_012920.genome : NC_012920
        NC_012920.chromosomes : chrM
        NC_012920.chrM.codonTable : Vertebrate_Mitochondrial


java -jar snpEff.jar  build -v -genbank NC_012920