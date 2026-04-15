"""
STAFF BAM: Base-4 XOR Analysis
Thin wrapper -- imports the original tool's run() function.
"""
import sys
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent.parent / "tools"
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

from staff_base4 import run, compare  # noqa: E402


def run_and_save(bam_path, output_path=None, max_reads=None):
    import json
    results = run(bam_path, max_reads=max_reads)
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
    return results
