import re
from collections import defaultdict

REGEX_TRUTH_FILE_ENTRY = r"(?P<tax_id>\d+)\t(?P<abs_abundance>\d+(?:\.\d*)?)\t(?P<rel_abundance>0\.\d+)?\t(?P<rank>[a-zA-Z]+)\t(?P<tax_name>.+)(?:\n|$)"

def parse_truth_file(truth_file):
  """Parses a TSV truth file with the following format.

  Arguments:
      filename {str} -- Path of the filename to parse
  """
  # Note: we ignore the relative abundances because it might have rounding errors
  # (e.g. there were entries like 0.000000)
  for line in truth_file:
    matches = re.match(REGEX_TRUTH_FILE_ENTRY, line)
    yield {
      'tax_id': int(matches.group('tax_id')),
      'rank': matches.group('rank'),
      'abs_abundance': int(float(matches.group('abs_abundance'))),
      'tax_name': matches.group('tax_name')
    }

def extract_truth(truth_files):
  truth_taxa = defaultdict(lambda: [])
  for truth_file in truth_files:
    for entry in parse_truth_file(truth_file):
      truth_taxa[entry['rank']].append(entry)

  return truth_taxa
