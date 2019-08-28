
import argparse
import json
import re
import sys
from os.path import exists
from .scoring.idseq import score

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
  parser.add_argument(dest='project_id', help='ID of the project', type=int)
  parser.add_argument(dest='sample_id', help='ID of the sample', type=int)
  parser.add_argument(dest='pipeline_version', help='Pipeline version to retrieve: <major>.<minor>', type=pipeline_version)
  parser.add_argument(
    '-o', '--output-file',
    dest='output_file', help='Path of output file to write the results to (if not set results are written to stdout)', type=output_file, default=sys.stdout)
  args = parser.parse_args()

  stats_json = score(args.project_id, args.sample_id, args.pipeline_version)

  args.output_file.write(json.dumps(stats_json, indent=2))

if __name__ == "__main__":
  main()
