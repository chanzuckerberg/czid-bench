#!/usr/bin/env python3
#
# IDSeq Benchmark Scorer.
#
# After running an idseq-bench sample through the IDSeq Portal,
# score the portal output as follows:
#
#     python3 score.py s3://idseq-samples-prod/samples/16/8848/results/2.8
#
# This computes the QC and recall rate for each sample organism.
#
import os
import sys
import json
import re
import subprocess
from fnmatch import fnmatch
from multiprocessing import cpu_count
from collections import defaultdict
from util import remove_safely, check_call, smart_open, check_output, ProgressTracker


# TODO: Use NamedTuple instead.
# List ranks in the same order as benchmark_lineage tags.
TAXID_RANKS = ["subspecies", "species", "genus", "family"]


def glob_sample_data(sample, version):
    "Enumerate pipeline stage data to be counted."
    return {
        "input_fastq":
            smart_glob(f"{sample}/fastqs/*.fastq*", expected=2),
        "post_qc_fasta":
            smart_glob(f"{sample}/results/{version}/gsnap_filter_[1,2].fa*", expected=2),
        "output_fasta":
            smart_glob(f"{sample}/postprocess/{version}/taxid_annot.fasta", expected=1)
    }


def smart_glob(pattern, expected, ls_memory={}):  # pylint: disable=dangerous-default-value
    pdir, file_pattern = pattern.rsplit("/", 1)
    def match_pattern(filename):
        return fnmatch(filename, file_pattern)
    matching_files = list(filter(match_pattern, smart_ls(pdir, memory=ls_memory)))
    actual = len(matching_files)
    assert expected == actual, f"Expected {expected} files for {pattern}, found {actual}"
    return [f"{pdir}/{mf}" for mf in sorted(matching_files)]


def smart_ls(pdir, missing_ok=True, memory=None):
    "Return a list of files in dir.  Cached."
    result = memory.get(pdir) if memory else None
    if result == None:
        try:
            if pdir.startswith("s3"):
                s3_dir = pdir
                if not s3_dir.endswith("/"):
                    s3_dir += "/"
                output = check_output(["aws", "s3", "ls", s3_dir])
                rows = output.strip().split('\n')
                result = [r.split()[-1] for r in rows]
            else:
                output = check_output(["ls", pdir])
                result = output.strip().split('\n')
        except Exception as e:
            msg = f"Could not read directory: {pdir}"
            if missing_ok and isinstance(e, subprocess.CalledProcessError):
                print(f"INFO: {msg}")
                result = []
            else:
                print(f"ERROR: {msg}")
                raise
        if memory != None:
            memory[pdir] = result
    return result


def parse_result_dir(sample_result_dir):
    "Deconstruct '<sample_path>/results/2.8' into ('<sample_path>', '2.8')"
    try:
        assert "/results/" in sample_result_dir
        before, after = sample_result_dir.split("/results/")
        version = after.split("/")[0]
        for version_part in version.split("."):
            assert str(int(version_part)) == version_part
        return before, version
    except:
        print(f"ERROR: Expected something like '.../results/5.1', got {sample_result_dir}")
        raise


def benchmark_lineage_from_header(header_line):
    matches = re.search(r'__benchmark_lineage_\d+_\d+_\d+_\d+__', header_line)
    assert matches, "Non-benchmark reads present"
    benchmark_lineage = matches.group(0)[2:-2]
    return benchmark_lineage


def benchmark_lineage_to_taxid_strs(benchmark_lineage):
    return tuple(taxid_str for taxid_str in benchmark_lineage.split("_")[2:])


def count_fastq(input_fastq):
    assert ".fastq" in input_fastq
    accumulators = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    with smart_open(input_fastq, "r") as input_f:
        line_number = 1
        try:
            line = input_f.readline()
            while line:
                # The FASTQ format specifies that each read consists of 4 lines,
                # the first of which begins with @ followed by read ID.
                line = line.decode('utf-8')
                assert line[0] == "@", f"fastq format requires every 4th line to start with @"
                benchmark_lineage = benchmark_lineage_from_header(line)
                taxid_strings = benchmark_lineage_to_taxid_strs(benchmark_lineage)
                for taxid_rank, taxid_str in zip(TAXID_RANKS, taxid_strings):
                    accumulators[benchmark_lineage][taxid_rank][taxid_str] += 1
                for _ in range(4):
                    line = input_f.readline()
                    line_number += 1
        except Exception as _:
            print(f"Error parsing line {line_number} in {input_fastq}.")
            raise
    return accumulators


def increment(accumulators, delta):
    for benchmark_lineage in delta.keys():
        for taxid_rank in delta[benchmark_lineage].keys():
            for taxid_str in delta[benchmark_lineage][taxid_rank].keys():
                accumulators[benchmark_lineage][taxid_rank][taxid_str] += delta[benchmark_lineage][taxid_rank][taxid_str]


def main(args):
    assert len(args) == 2, "Sample dir argument is required.  See usage."
    sample_result_dir = args[1]
    print(f"Scoring IDSEQ benchmark output {sample_result_dir}")
    sample, version = parse_result_dir(sample_result_dir)
    sample_data = glob_sample_data(sample, version)
    print(json.dumps(sample_data, indent=4))
    r1, r2 = sample_data["input_fastq"]
    input_counts = count_fastq(r1)
    increment(input_counts, count_fastq(r2))
    print(json.dumps(input_counts, indent=4))


if __name__ == "__main__":
    main(sys.argv)
