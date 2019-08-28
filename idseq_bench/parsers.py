import re

FAST_FILE_TYPE = r"\.(?:fast|f)(?P<type>q|a)(?:\.|$)"

def extract_accession_id(raw_line):
  try:
    return re.search(r'^[@>](.+?\.\d+)', raw_line).group(1)
  except AttributeError:
    return None

def extract_fast_file_type_from_path(file_path):
  matches = re.findall(FAST_FILE_TYPE, file_path)
  if not matches:
    raise AttributeError(f"Unable to detect fast file type for '{file_path}'")
  return matches[-1]
