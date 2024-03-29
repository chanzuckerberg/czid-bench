
import argparse
import json
import re
import sys
from os.path import exists
from .scoring.czid import score_benchmark, score_sample
from .scoring.truth import extract_truth

def pipeline_version(pipeline_version_str):
  matched = re.match(r"\d+\.\d+$", pipeline_version_str)
  if not matched:
    raise ValueError()
  return pipeline_version_str

def output_file(output_file_path):
  print(output_file_path)
  if output_file_path and not output_file_path.endswith('.json') and not output_file_path == '-':
    raise argparse.ArgumentTypeError("Output file does not terminate in '.json'")

  if exists(output_file_path):
    overwrite = input(f"File {output_file_path} already exists. Overwrite (y/N)? ")
    if overwrite != 'y':
      raise argparse.ArgumentTypeError("Choose another output file.")

  return argparse.FileType('w')(output_file_path)

def main():
  parser = argparse.ArgumentParser(
    description="Score benchmark samples based on answer key (against key encode in fastq "
                "read IDs for IDseq Bench samples or against a provided answer key).")
  parser.add_argument(dest='project_id', type=int, help='ID of the project')
  parser.add_argument(dest='sample_id', type=int, help='ID of the sample')
  parser.add_argument(dest='pipeline_version', type=pipeline_version, help='Pipeline version to retrieve: <major>.<minor>')
  parser.add_argument(
    '-o', '--output-file', type=output_file, default=sys.stdout, dest='output_file',
    help='Path of output file to write the results to (if not set results are written to stdout)')
  parser.add_argument(
    '-t', '--truth_files', type=argparse.FileType('r'), dest='truth_files', nargs='*',
    help='Space separated list of truth files. Truth files are TSV files with the fields <tax_id, absolute abundance, relative abundance, rank, tax_name>'
  )
  parser.add_argument(
    '-e', '--env', type=str, dest='env',
    help='String name of environment in which sample exists; options are: prod (default), staging'
  )
  parser.add_argument(
    '-p', '--local-path', type=str, dest='local_path',
    help='Root of local path where IDseq files are stored. Must follow the same structure to store files as IDseq'
  )
  parser.add_argument(
    '--mono-aupr', dest='mono_aupr', action='store_true',
    help='Computed adjusted AUPR by forcing the precision/recall curve to be monotonic decreasing'
  )
  args = parser.parse_args()

  stats_json = None
  if args.truth_files:
    truth_taxa = extract_truth(args.truth_files)
    stats_json = score_sample(args.project_id, args.sample_id, args.pipeline_version, truth_taxa, env=args.env, local_path=args.local_path, force_monotonic=args.mono_aupr)
  else:
    stats_json = score_benchmark(args.project_id, args.sample_id, args.pipeline_version, env=args.env, local_path=args.local_path, force_monotonic=args.mono_aupr)

  if stats_json:
    args.output_file.write(json.dumps(stats_json, indent=2))
  else:
    print("[ERROR] Unable to score sample")


if __name__ == "__main__":
  main()
