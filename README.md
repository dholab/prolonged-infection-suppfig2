## Overview

This workflow was created to efficiently process and visualize data presented in Supplemental Figure 2 of [Halfmann et al. 2022, _Evolution of a globally unique SARS-CoV-2 Spike E484T monoclonal antibody escape mutation in a persistently infected, immunocompromised individual_](https://www.medrxiv.org/content/10.1101/2022.04.11.22272784v1). Consensus sequences and metadata were pulled from GenBank and are downloaded automatically in this workflow (see below).

## Getting Started

To run this workflow, simply `git clone` it into your working directory of choice, like so:

```
git clone https://github.com/nrminor/prolonged-infection-suppfig2.git .
```

Once the workflow bundle is in place, first ensure that the workflow scripts are executable, like so:

```
chmod +x bin/*
```

Next, build the Docker image that contains the workflow's dependencies:

```
docker build --tag wisc-prolonged-infection:v1_0_1 config/
```

Note that to build the above docker container, you may need to increase the amount of memory allotted to Docker in the Docker Engine preferences.

### Nextflow Installation

This workflow uses the [NextFlow](https://www.nextflow.io/) workflow manager. We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. Install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

1. Run the following line in a directory where you'd like to install NextFlow, and run the following line of code:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, you are ready to proceed.

To run the workflow, simply change into the workflow directory and run the following in the BASH terminal:

```
nextflow run prolonged_infection_suppfig2.nf
```

If the workflow runs partway, but a computer outage or other issue interrupts its progress, no need to start over! Instead, run:

```
nextflow run prolonged_infection_suppfig2.nf -resume
```

The workflow's configurations (see below) tell NextFlow to plot the workflow and record run statistics. However, the plot the workflow, note that NextFlow requires the package GraphViz, which is easiest to install via the intructions on [GraphViz's website](https://graphviz.org/download/).

### Configuration

The following runtime parameters have been set for the whole workflow:

- `subsample_size` - the number of SARS-CoV-2 samples to compare with the persistent infection in the final plot. We chose 5000 as our default, but this can be changed to any number as shown in the code block below. _NOTE:_ If you choose a number other than our default, the workflow will need to gather a new subsample from GenBank. To prompt it to do so, simply delete the file `include_list.csv` from the `resources/` directory. This process of pulling from GenBank may take as long as 2 days.
- `min_date` - the earliest date to display in the plot; default is September 1st, 2020. Dates specified here must be in "YYYY-mm-dd" format.
- `max_date` - the latest date to display in the plot; default is February 1st, 2022. Dates specified here must be in "YYYY-mm-dd" format.
- `random_sample_seed` - A seed that determines which GenBank accessions are pseudo-randomly selected by the script `select_subsample.R`. Our default value for this seed is 14; changing this will result in a final figure that is not identical to the figure in our manuscript.
- `refseq` - path to the SARS-CoV-2 reference sequence (GenBank Accession MN9089473.3), which is in the workflow subdirectory `resources/`
- `results` - path to default workflow output directory
- `results_data_files` - path to a subdirectory of `results/` where data files like VCFs are placed.
- `visuals` - path to a subdirectory of `results/` where graphics are stored. This is where the final supplementary figure 2 PDF is placed.

These parameters can be altered in the command line with a double-dash flag, like so:

```
nextflow run prolonged_infection_workflow.nf --min_date "2020-07-01"
```

## Workflow summary

- First, in the process PULL_METADATA, the workflow looks for the pre-existing include file specified in the configuration file, `nextflow.config`. By default, this file is called `include_list.csv`. If this process does not find a pre-existing include file, it will pull metadata for all complete SARS-CoV-2 genomes in GenBank. To do so, it uses the excellent [NCBI datasets command-line interface](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/command-line/datasets/). This has generally taken 4 to 6 hours on a standard desktop computer.
- In the process REFORMAT_METADATA, the workflow either sees that the include file already exists, as in the process PULL_METADATA, or it converts the GenBank metadata from PULL_METADATA into a convenient tab-separated format. To do so, it uses the handy [NCBI dataformat tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/command-line/dataformat/).
- In the process SELECT_SUBSAMPLE, the workflow either a) sees that the include file already exists and streams it into the next process, or b) or takes the TSV-formatted metadata from REFORMAT_METADATA and generates a new include file with the R script `select_subsample.R`. Note that the workflow configuration file, `nextflow.config`, specifies 4 parameters for this script, as described above: `subsample_size`, `min_date`, `max_date`, and  `random_sample_seed`. These settings can all be changed to produce different results with different input data.
- Next, in the process PULL_FASTAs, the workflow splits the include file row by row, and pulls a FASTA file separately for the accession in each row. This makes it possible for multiple FASTAs to be pulled and then pushed to the next step at once.
- Each GenBank FASTA is then aligned to Wuhan-1 in the process SUBSAMPLE_ALIGNMENT. As in PULL_FASTAs, these alignments can take place in parallel to one another, making the workflow more efficient.
- In the process SUBSAMPLE_VARIANT_CALLING, we use the `callvariants.sh` script from [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) to identify mutations in each of the GenBank accessions.
- Finally, in the process SUPP_FIGURE_2_PLOTTING, all VCFs from SUBSAMPLE_VARIANT_CALLING are brought together with some pre-existing inputs and plotted with the script `SupplementalFigure2_global_roottotip_plot.R`. The final plot is placed in `results/visuals/`, and the associated plotting data are placed in `results/` as `{run date}.csv`. With that, the workflow is finished!

The steps described above are also visualized in the file [prolonged_infection_suppfig2_dag.png](https://github.com/nrminor/prolonged-infection-suppfig2/blob/332cff270dafc8213b8135fc21feba2e711f5ce4/prolonged_infection_suppfig2_dag.png).

### Bundled data

Bundled together with the workflow are:

- the SARS-CoV-2 Wuhan-1 sequence from GenBank Accession MN9089473.3 is in `resources/` and is called reference.fasta
- the file `data/patient_variant_counts.csv` comes from the FIGURE_2A_PLOTTING process in the main workflow, which is available at the [GitHub repository for this project](https://github.com/dholab/E484T-visualizations/tree/main).
- a list of GenBank accessions to include in the final plot called `include_list.csv`. If this file is absent, the workflow will pull the SARS-CoV-2 metadata from GenBank and generate the include list again. _Note_ that the R script in the `bin/` directory, `select_subsample.R`, uses a consistent seed 

### Output files

- `results/visuals/figsupp2_global_roottotip_plot.pdf` is the final plot, which can be edited in Illustrator or other vector-graphic software
- `results/{run date}.csv` is a table of all the GenBank accessions used in making the plot, along with when they were collected, how many mutations they have that separate them from Wuhan-1, and what pango lineage they are. This is essentially the dataset that is used to make Supplemental Figure 2.
- `results/data/*.vcf` are all the VCFs generated to make the figure. By default, these files are not included in the git repository, but they will be generated if you run the workflow yourself.

## For more information

This workflow was created by Nicholas R. Minor. To report any issues, please visit the [GitHub repository for this project](https://github.com/nrminor/prolonged-infection-suppfig2).
