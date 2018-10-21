# idseq-bench
Benchmark generator for the [IDSeq Portal](https://idseq.net).

So far just a thin wrapper around [InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/).

# setup
```
pip3 install InSilicoSeq
pip3 install ncbi-acc-download
```

# running
```
python3 ../generate.py
```

This produces zipped fastq files that you can upload to the [IDSeq Portal](https://idseq.net) via [IDSEQ-CLI](https://github.com/chanzuckerberg/idseq-cli).

# selecting organisms and chromosomes
Add/modify entries like this in `params.py`
```
    Genome("fungi", "aspergillus_fumigatus",
           [("subspecies", 330879), ("species", 746128), ("genus", 5052), ("family", 1131492)],
           ["NC_007194.1", "NC_007195.1", "NC_007196.1", "NC_007197.1", "NC_007198.1", "NC_007199.1",
            "NC_007200.1", "NC_007201.1"],
           "https://www.ncbi.nlm.nih.gov/genome/18?genome_assembly_id=22576"),
```

# tweaking InSilicoSeq options
Edit [params.py](params.py) or [main.py](main.py) as desired, e.g., to select a different set of [error models](https://insilicoseq.readthedocs.io/en/latest/iss/model.html).

If the output data changes, we expect the output file name to change as well.  It's always a good idea to increment
LOGICAL_VERSION (whole numbers only) after making changes to the code.

If your code or paramter changes are uncommitted, the program will refuse to run.   You can force it to run with the `--irreproducible` flag, but irreproducible outputs cannot be used for automated testing of the IDSeq Portal.

# interpreting the output
Each output file name reflects the params of its generation, like so:
```
norg_6__nacc_27__uniform_weight_per_organism__hiseq_reads__v4__[R1, R2].fastq.gz
  -- number of organisms: 6
  -- number of accessions: 27
  -- distribution: uniform per organism
  -- error model: hiseq
  -- logical version: 4
```
TODO:  Random generator seed control.

We generate a summary file for each pair of fastqs, indicating read counts per organism,
and the average coverage of the organism's genome.  Each pair counts as 2 reads / 300 bases,
matching InSilicoSeq and IDSeq conventions.
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

# automated testing of IDSeq Portal

Just upload an output folder to `s3://idseq-bench/<next-number>` and add
an entry for it to `s3://idseq-bench/config.json` to specify frequency and environments in which that test should run.
