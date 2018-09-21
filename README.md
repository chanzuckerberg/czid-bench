# idseq-bench
Infectious disease benchmarking.   So far just a thin wrapper around InSilicoSeq.

# setup
```
pip3 install InSilicoSeq
pip3 install ncbi-acc-download
```

# selecting organisms and chromosomes
Add entries like this one to `params.py`
```
    Genome("fungi", "aspergillus_fumigatus",
           [("subspecies", 330879), ("species", 746128), ("genus", 5052), ("family", 1131492)],
           ["NC_007194.1", "NC_007195.1", "NC_007196.1", "NC_007197.1", "NC_007198.1", "NC_007199.1",
            "NC_007200.1", "NC_007201.1"],
           "https://www.ncbi.nlm.nih.gov/genome/18?genome_assembly_id=22576"),
```

# tweaking InSilicoSeq options
Edit the `iss` command in `main()`.  Increment LOGICAL_VERSION (whole numbers only);  that's reflected in the output names.  If possible, include all options in the output name, because uncommitted versions are prone to clashing.

# running
python3 main.py

This produces fastq files that you can gzip and upload via IDSEQ-CLI to project Benchmarking.
