# dp_shim.py
import sys, os, shlex, runpy

dp_path = os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin", "diagonal_partition.py")
if not os.path.isfile(dp_path):
    dp_path = "./scripts/diagonal_partition.py"


for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    args = shlex.split(line, posix=True)
    # Emulate: python diagonal_partition.py -1 <args...>
    old_argv = sys.argv
    try:
        sys.argv = ["diagonal_partition.py", "-1", *args]
        try:
            runpy.run_path(dp_path, run_name="__main__")
        except SystemExit as e:
            # diagonal_partition.py may call sys.exit(); treat 0 as success
            if int(getattr(e, "code", 0) or 0) != 0:
                raise
    finally:
        sys.argv = old_argv
