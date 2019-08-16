#!/usr/bin/env python3
#
# IDseq Benchmark Generator.
#
# Prerequsites:
#
#    pip3 install InSilicoSeq
#    pip3 install ncbi-acc-download
#
# Usage:
#
#    Customize params.py, then run main.py, then use idseq-cli to upload output to idseq.
#
import os
import re
import sys
import json
import argparse
from multiprocessing import cpu_count
from collections import defaultdict
import yaml
from .util import remove_safely, check_call, smart_open, ProgressTracker
from .genome import Genome
from . import __version__


class BenchmarkConfigError(Exception):
    def __init__(self, errors):
        super().__init__()
        self.errors = errors

    def __str__(self):
        errors_str = '\t' + '\n\t'.join(self.errors)
        return f"Benchmark config contains errors:\n{errors_str}"


class ISSRunContext:
    """Execution context for each InSilicoSeq run.  Encapsulates all temp/intermediate/output files, but not the logic.
    """

    iss_version = None

    def __init__(self, tmp_prefix, output_prefix):
        print(f"GENERATING {output_prefix}")
        if ISSRunContext.iss_version == None:
            ISSRunContext.iss_version = get_iss_version()
        self.idseq_bench_command = f"{os.path.basename(sys.argv[0])} {' '.join(sys.argv[1:])}"
        self.subdir = output_prefix
        self.output_files = self.paired_files(output_prefix, "fastq.gz")
        self.tmp_files = self.paired_files(tmp_prefix, "fastq")
        self.summary_file_txt = f"{self.subdir}/summary.txt"
        self.metadata_file_json = f"{self.subdir}/metadata.json"
        self.abundance_file = f"{self.subdir}/{tmp_prefix}_abundance.txt"
        self.genomes_file = f"{self.subdir}/{tmp_prefix}_genomes.fasta"
        self.clean_slate()

    def clean_slate(self):
        os.makedirs(self.subdir, exist_ok=True)
        for f in self.tmp_files + self.output_files + [self.summary_file_txt, self.metadata_file_json, self.abundance_file, self.genomes_file]:
            remove_safely(f)

    def cleanup(self):
        for f in self.tmp_files + [self.abundance_file, self.genomes_file]:
            remove_safely(f)

    def paired_files(self, prefix, suffix):
        return [f"{self.subdir}/{prefix}_{r}.{suffix}" for r in ["R1", "R2"]]


def benchmark_lineage_tag(g):
    return f"benchmark_lineage_{g.subspecies_taxid}_{g.species_taxid}_{g.genus_taxid}_{g.family_taxid}"


def extract_accession_id(raw_line):
    try:
        return re.search(r'^(@.+\.\d+)', raw_line).group(1)
    except AttributeError:
        return None


def augment_and_count_read_header(line, line_number):
    iss_read_id = extract_accession_id(line)
    serial_number = "s{:010d}".format(line_number // 4)
    g = Genome.by_accid[iss_read_id[1:]]
    benchmark_lineage = benchmark_lineage_tag(g)
    return f"{iss_read_id}__{benchmark_lineage}__{serial_number}\n", g.key


def annotate_and_count_reads(input_fastq, output_fastq, counters, accumulators):
    """Annotate read IDs by appending the consecutive read counter, after stripping
    that _1 or _2  paired-end indicator appended by ISS.  Required to run correctly
    through STAR.  Also annotate with lineage information, that can later be used
    to score idseq accuracy."""
    with smart_open(input_fastq, "rb") as input_f, smart_open(output_fastq, "w") as output_f:
        line_number = 1
        try:
            line = input_f.readline()
            while line:
                # The FASTQ format specifies that each read consists of 4 lines,
                # the first of which begins with @ followed by read ID.
                line = line.decode('utf-8')
                assert line[0] == "@", f"fastq format requires every 4th line to start with @"
                augmented_read_header, g_key = augment_and_count_read_header(line, line_number)
                counters[g_key] += 1
                output_f.write(augmented_read_header.encode('utf-8'))
                for i in range(4):
                    line = input_f.readline()
                    line_number += 1
                    if i == 0:
                        accumulators[g_key] += (len(line) - 1)
                    if i < 3:
                        output_f.write(line)
        except Exception as _:
            print(f"Error parsing line {line_number} in {input_fastq}.")
            raise


def output_summary_counters(rc, iss_command, counters, accumulators, **extra_metadata):
    contents = []
    total_reads = 0
    with smart_open(rc.summary_file_txt, "w") as mft:
        headers = "READ_COUNT\tREAD_SIZE\tCOVERAGE\tLINEAGE\tGENOME\n"
        mft.write(headers)
        print("")
        print(headers)
        for g_key, read_count in sorted(counters.items(), key=lambda pair: pair[1], reverse=True):
            g = Genome.all[g_key]
            benchmark_lineage = benchmark_lineage_tag(g)
            coverage = accumulators[g.key] / g.size
            read_size = accumulators[g.key] / read_count
            summary_line = f"{read_count}\t{read_size}\t{coverage:3.1f}x\t{benchmark_lineage}\t{g.key}\n"
            print(summary_line)
            mft.write(summary_line)
            contents.append({
                'genome': g.key,
                'benchmark_lineage': benchmark_lineage,
                'coverage': coverage,
                'read_count': read_count,
                'read_size' : read_size,
                'genome_size': g.size
            })
            total_reads += read_count

    metadata = extra_metadata or {}
    metadata.update({
        "fastqs": [os.path.basename(f) for f in rc.output_files],
        "iss_command": iss_command,
        "iss_version": rc.iss_version,
        "idseq_bench_command": rc.idseq_bench_command,
        "idseq_bench_version": __version__,
        "prefix": os.path.basename(rc.subdir),
        "verified_total_reads": total_reads,
        "verified_contents": contents
    })

    with smart_open(rc.metadata_file_json, "w") as mdf:
        mdf.write(json.dumps(metadata, indent=4))


def get_iss_version():
    remove_safely("iss_version.txt")
    check_call("iss --version > iss_version.txt")
    with open("iss_version.txt") as ivf:
        iss_version = ivf.readline().strip()
    assert iss_version.startswith("iss version ")
    iss_version = iss_version[len("iss version "):]
    remove_safely("iss_version.txt")
    return iss_version


def strictly_above(va, vb):
    for a, b in zip(va.split("."), vb.split(".")):
        if int(a) > int(b):
            return True
        if int(a) < int(b):
            return False
    return False


def run_iss(rc, iss_command, **extra_metadata):
    check_call(iss_command)
    counters = defaultdict(int)
    accumulators = defaultdict(int)
    for tmp_fastq, output_fastq in zip(rc.tmp_files, rc.output_files):
        annotate_and_count_reads(tmp_fastq, output_fastq, counters, accumulators)
    output_summary_counters(rc, iss_command, counters, accumulators, **extra_metadata)
    rc.cleanup()


def uniform_abundance_per_organism(genomes, abundance_file):
    """Ensures each organism has equal abundance in the mix, despite the varying numbers of chromosomes and accessions across organisms.
    """
    with open(abundance_file, "w") as af:
        sum_of_weights = 0
        for g in genomes:
            vaccid_weight = 1.0 / (len(g.versioned_accession_ids) * len(genomes))
            for vaccid in g.versioned_accession_ids:
                sum_of_weights += vaccid_weight
                af.write(f"{vaccid} {vaccid_weight}\n")
    assert -0.0005 < sum_of_weights - 1.0 < 0.0005, f"{sum_of_weights} != 1.0"


def concatenate_fasta(genomes, genomes_file):
    genome_fastas = " ".join(g.filename for g in genomes)
    command = f"cat {genome_fastas} > {genomes_file}"
    check_call(command)


def run_iss_multiplexed(genomes, num_reads, model, tmp_prefix, num_cpus, abundance='uniform', description=None):
    """Runs iss with all genomes multiplexed in the same file.

    Returns:
        str -- The output path where generated benchmark is stored.
    """
    print("Initializing and downloading genomes...")
    num_organisms = len(genomes)
    num_accessions = sum(len(g.versioned_accession_ids) for g in genomes)
    output_prefix_multiplexed = f"norg_{num_organisms}__nacc_{num_accessions}__{abundance}_weight_per_organism__{model}_reads__v{__version__}"
    rc = ISSRunContext(tmp_prefix, output_prefix_multiplexed)
    uniform_abundance_per_organism(genomes, rc.abundance_file)
    concatenate_fasta(genomes, rc.genomes_file)
    iss_multiplexed_command = f"iss generate --n_reads {num_reads} --genomes {rc.genomes_file} --model {model} --abundance_file {rc.abundance_file} --gc_bias --output {rc.subdir}/{tmp_prefix} --cpus {num_cpus}"
    run_iss(rc, iss_multiplexed_command, description=description)
    return rc.subdir


def create_benchmark(benchmark_config):
    # genomes, models, reads_per_organism, benchmark_config, abundance='uniform', description=None):


    print("Generating IDseq benchmark data...")
    genomes = initialize_genomes(benchmark_config['genomes'])
    num_cpus = cpu_count()
    num_reads = benchmark_config['reads_per_organism'] * len(genomes)
    overall_reads = num_reads * len(benchmark_config['models'])
    Genome.fetch_all()
    pid = os.getpid()
    tmp_prefix = f"tmp_{pid}"
    ticker = ProgressTracker(target=overall_reads)
    for model in benchmark_config['models']:
        # Then generate a multiplexed benchmark.
        output_path = run_iss_multiplexed(genomes, num_reads, model, tmp_prefix, num_cpus, benchmark_config['abundance'], description=benchmark_config['description'])
        ticker.advance(num_reads)

        with open(f'{output_path}/{output_path}.yaml', 'w') as output_file:
            yaml.dump(benchmark_config, output_file, default_flow_style=False)


def initialize_genomes(genome_configs):
    print(genome_configs)
    return [
        Genome(**genome_config)
        for genome_config in genome_configs
    ]

def parse_and_validate_config(config_file):
    benchmark_config = yaml.safe_load(config_file)

    # Add basic validation for required fields
    # TODO: deep and optional fields validation
    errors = [
        f"Config file must contain a field '{field}' of type <{field_type}>"
        for field, field_type in [
            ('genomes', list),
            ('models', list),
            ('reads_per_organism', int)
        ]
        if not field in benchmark_config or not isinstance(benchmark_config[field], field_type)
    ]

    if len(errors) > 0:
        raise BenchmarkConfigError(errors)

    return benchmark_config

def valid_path(path):
    if os.path.isdir(path):
        return path
    os.mkdir(path)
    return path

def main():
    parser = argparse.ArgumentParser(description="Generate MGS benchmark fastq files.")
    parser.add_argument(dest='config_file', help="Benchmark configuration file.", type=argparse.FileType('r', encoding='UTF-8'))
    parser.add_argument(
        '-d', '--downloads-dir', dest='downloads_dir',
        help="Folder to store downloaded genomes (if genomes already exist in this folder they will not be downloaded).", type=valid_path)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    args = parser.parse_args()

    # default downloads is a static variable in Genome class
    if args.downloads_dir:
        Genome.downloads_dir = args.downloads_dir

    try:
        benchmark_config = parse_and_validate_config(args.config_file)
    except BenchmarkConfigError as e:
        print(e)
    else:
        create_benchmark(benchmark_config)


if __name__ == "__main__":
    main()
