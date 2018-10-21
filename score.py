#!/usr/bin/env python3
#
# IDSeq Benchmark Scorer.
#
# After running an idseq-bench sample through the IDSeq Portal,
# score the portal output as follows:
#
#     python3 score.py s3://idseq-samples/16/8848/results/2.8
#
# This computes the QC and recall rate for each sample organism.
#
import os
import sys
import json
from fnmatch import fnmatch
from multiprocessing import cpu_count
from collections import defaultdict
from util import remove_safely, check_call, smart_open, check_output, ProgressTracker


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


def smart_glob(pattern, expected):
    dir, file_pattern = pattern.rsplit("/", 1)
    def match_pattern(filename):
        return fnmatch(filename, file_pattern)
    matching_files = list(filter(match_pattern, cached_ls(dir)))
    actual = len(matching_files)
    assert expected == actual, f"Expected {expected} files for {pattern}, found {actual}"
    return [f"{dir}/{mf}" for mf in sorted(matching_files)]


def cached_ls(dir, cache={}):
    "Return a list of files in dir.  Cached."
    if dir not in cache:
        if dir.startswith("s3"):
            s3_dir = dir
            if not s3_dir.endswith("/"):
                s3_dir += "/"
            output = check_output(["aws", "s3", "ls", s3_dir])
            rows = output.strip().split('\n')
            cache[dir] = [r.split()[-1] for r in rows]
        else:
            output = check_output(["ls", dir])
            cache[dir] = output.strip().split('\n')
    return cache[dir]


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


def main(args):
    assert len(args) == 2, "Sample dir argument is required.  See usage."
    sample_result_dir = args[1]
    print(f"Scoring IDSEQ benchmark output {sample_result_dir}")
    sample, version = parse_result_dir(sample_result_dir)
    sample_data = glob_sample_data(sample, version)
    print(json.dumps(sample_data, indent=4))


if __name__ == "__main__":
    main(sys.argv)
