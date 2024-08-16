# Fusion Transcript Benchmarking

Repo: git@github.com:fusiontranscripts/LR-FusionBenchmarking.git

## Benchmarking simulated fusion reads

- JAFFAL BadRead simulated divergent fusion reads:

  - [simulated_data/sim_jaffal](simulated_data/sim_jaffal)

- simulated hi fidelity reads

    - pacbio pbsim3 simualted reads: [simulated_data/our_pbsim3_sims/pbio_pbsim3_part5](simulated_data/our_pbsim3_sims/pbio_pbsim3_part5)

    - ONT pbsim3 simulated reads: [simulated_data/our_pbsim3_sims/ONT_pbsim3_part5](simulated_data/our_pbsim3_sims/ONT_pbsim3_part5)
    
- simulated reads for fusions between paralogous genes

    - [simulated_data/paralog_fusion_sim](simulated_data/paralog_fusion_sim)

## Benchmarking SeraCare Seraseq fusion standard

- [SeraCareFusions](SeraCareFusions)

## Benchmarking DepMap Nine Cell Lines (MAS-ISO-seq)

- [DepMap_Cell_Lines](DepMap_Cell_Lines)

## Benchmarking ONT direct RNA on three cell lines from the [SG-NEx project](https://registry.opendata.aws/sgnex/)

- [SGNex_ONT](SGNex_ONT)

## Installation notes
    
In each benchmarking analysis below, running the 'analyze...' script will run through all benchmarking.

In the case you find a Makefile, see the Makefile for options in benchmarking.

>Before running, install UpSetRbyFeature (see the top of the README there).
    
Prediction results from each method are organized in the corresponding prog_results/ directory.
    
>Workflows for running the prediction methods and analyses of the above benchmarking results are organized at [github.com/broadinstitute/CTAT-LRF-Paper](https://github.com/broadinstitute/CTAT-LRF-Paper)

