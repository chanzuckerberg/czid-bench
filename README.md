# [CZ ID](https://czid.org/) &middot; [![GitHub license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/chanzuckerberg/czid-web/blob/master/LICENSE) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)


#### Infectious Disease Sequencing Platform
CZ ID is an unbiased global software platform that helps scientists identify pathogens in metagenomic sequencing data.

- **Discover** - Identify the pathogen landscape
- **Detect** - Monitor and review potential outbreaks
- **Decipher** - Find potential infecting organisms in large datasets

A collaborative open project of [Chan Zuckerberg Initiative](https://www.chanzuckerberg.com/) and [Chan Zuckerberg Biohub](https://czbiohub.org).

Check out our repositories:
- [czid-web](https://github.com/chanzuckerberg/czid-web) - Frontend portal
- [czid-workflows](https://github.com/chanzuckerberg/czid-workflows) - Bioinformatics workflows
- [czid-cli](https://github.com/chanzuckerberg/czid-cli) - Command line upload interface
- [czid-bench](https://github.com/chanzuckerberg/czid-bench) - Pipeline benchmarking tools (here)

# czid-bench
Benchmark generator for the [CZ ID Portal](https://czid.org).

So far just a thin wrapper around [InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/).

## setup
```
pip3 install git+https://github.com/chanzuckerberg/czid-bench.git --upgrade
```

## running
```
czid-bench-generate config_file.yaml
```

This produces zipped fastq files and config files to generate them. You can upload the fastq files to the [CZ ID Portal](https://czid.org) via [IDSEQ-CLI](https://github.com/chanzuckerberg/czid-cli).

## help
```
czid-bench-generate -h
```


## selecting organisms and chromosomes
Create a yaml file in the following format:
```
# A readable name for the benchmark
description: List of relevant genomes to use on standard benchmarks
# Number of reads per organism
reads_per_organism: 10000
# The sequencer model to emulate (determines the error model used by iss)
# Possible values: novaseq, miseq, hiseq
# It will generate one benchmark per specified model
models:
  - hiseq
abundance: uniform
genomes:
  - category: fungi
    organism: aspergillus_fumigatus
    lineage:
      - level: subspecies
        tax_id: 330879
      - level: species
        tax_id: 746128
      - level: genus
        tax_id: 5052
      - level: family
        tax_id: 1131492
    versioned_accession_ids:
      - NC_007194.1
      - NC_007195.1
      - NC_007196.1
      - NC_007197.1
      - NC_007198.1
      - NC_007199.1
      - NC_007200.1
      - NC_007201.1
    genome_assembly_url: https://www.ncbi.nlm.nih.gov/genome/18?genome_assembly_id=22576
```

See more examples in the examples folder.

## tweaking InSilicoSeq options
You can select different sets of [error models](https://insilicoseq.readthedocs.io/en/latest/iss/model.html).

The generated filenames will include the package version used to create it.

## interpreting the output
Each output file name reflects the params of its generation, like so:
```
norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v0.1.0__[R1, R2].fastq.gz
  -- number of organisms: 6
  -- number of accessions: 27
  -- distribution: uniform per organism
  -- error model: hiseq
  -- logical version: 4
```

We generate a summary file for each pair of fastqs, indicating read counts per organism,
and the average coverage of the organism's genome.  Each pair counts as 2 reads / 300 bases,
matching InSilicoSeq and CZ ID conventions.
```
READS  COVERAGE    LINEAGE                                          GENOME
----------------------------------------------------------------------------------------------------------------------
16656    215.3x    benchmark_lineage_0_37124_11019_11018            viruses__chikungunya__37124
16594      0.1x    benchmark_lineage_330879_746128_5052_1131492     fungi__aspergillus_fumigatus__330879
16564    352.1x    benchmark_lineage_0_463676_12059_12058           viruses__rhinovirus_c__463676
16074      0.5x    benchmark_lineage_1125630_573_570_543            bacteria__klebsiella_pneumoniae_HS11286__1125630
15078      0.8x    benchmark_lineage_93061_1280_1279_90964          bacteria__staphylococcus_aureus__93061
14894      0.1x    benchmark_lineage_36329_5833_5820_1639119        protista__plasmodium_falciuparum__36329
```
We alter the read id's generated by ISS to satisfy the input requirements of tools like [STAR](https://github.com/alexdobin/STAR).
This requires stripping those `_1` and `_2` pair indicators from all read id's, so that both reads in a pair have the exact same
read ID.  Then, each read ID gets a serial number, and a tag identifying the taxonomic lineage of the organism
that read was sourced from, like so:
```
@NC_016845.1_503__benchmark_lineage_1125630_573_570_543__s0000001169
```
This is helpful in tracking reads through complex bioinformatic pipelines and
scoring results.  We assume the pipelines would not cheat by inspecting those tags.

An even more detailed summary, including all ISS options, is generated in json format.

## For CZ ID developers: automated testing of CZ ID Portal

Just upload an output folder to `s3://czid-bench/<next-number>` and add
an entry for it to `s3://czid-bench/config.json` to specify frequency and environments in which that test should run.

## scoring a CZ ID Portal Run

After a benchmark sample has completed running through the CZ ID Portal, the QC pass rate and recall per benchmark organism can be scored by running, e.g.,
```
czid-bench-score <project_id> <sample_id> <pipeline_version:major.minor>
```
which produces JSON formatted output like so
```
{
  "per_rank": {
    "family": {
      "NT": {
        "543": {
          "total_reads": 10000,
          "post_qc_reads": 8476,
          "recall_per_read": {
            "count": 8461,
            "value": 0.9982302973100519
          }
        },
        ...
        "accuracy": {
          "count": 80137,
          "value": 0.8820803522289489
        },
        "total_simulated_taxa": 12,
        "total_correctly_identified_taxa": 11,
        "total_identified_taxa": 539,
        "recall": 0.9166666666666666,
        "precision": 0.02040816326530612,
        "f1-score": 0.03992740471869328
        "aupr": 0.9751017478206347,
        "l1_norm": 0.8389712437238702,
        "l2_norm": 0.07556827265305112
      },
      "NR": {
        "543": {
          "total_reads": 10000,
          "post_qc_reads": 8476,
          "recall_per_read": {
            "count": 7951,
            "value": 0.9380604058518169
          }
        },
        ...
      },
      "concordance": {
        "11018": {
          "count": 16048,
          "value": 1.9154929577464788
        },
        ...
    },
    "genus": {
      "NT": {
        "570": {
          ...
```

### Local files

For users who lack direct access to S3, scoring also works on a local download of sample results.  However, you must organize any locally downloaded files in versioned subfolders to match the S3 structure illustrated in the example above. Use the option `-p <local_path>` or `--local-path <local_path>` to use the local folder instead.

### Comparison to ground truth

Users can also compare any sample against a provided ground truth file. This file should be a TSV file with the following fields (without headers):
```
<taxon_id>	<absolute_abundance>	<relative_abundance>	<rank>	<taxon_name>
```

e.g

```
366648	100000.00000	0.01746	species	Xanthomonas fuscans
1685	100000.00000	0.01746	species	Bifidobacterium breve
486	100000.00000	0.01746	species	Neisseria lactamica
2751	100000.00000	0.01746	species	Carnobacterium maltaromaticum
28123	100000.00000	0.01746	species	Porphyromonas asaccharolytica
118562	100000.00000	0.01746	species	Arthrospira platensis
...
```

To compare against a ground truth run the scoring script with the following options:

```
czid-bench-score <project_id> <sample_id> <pipeline_version:major.minor> -t <truth_file_1.tsv> <truth_file_2.tsv> ...
```

## help
```
czid-bench-score -h
```

## Contributing

This project adheres to the Contributor Covenant code of conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

## Security

Please disclose security issues responsibly by contacting security@chanzuckerberg.com.
