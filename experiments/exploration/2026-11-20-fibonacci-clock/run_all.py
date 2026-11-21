from pathlib import Path
from subprocess import run
root=Path(__file__).resolve().parent
scripts=(
    "five_anyon_exact.py",
    "a4_double_e8.py",
    "fibonacci_history_kms.py",
    "binary_two_clocks.py",
    "braid_readout_chain.py",
    "affine_su5_e8.py",
    "cuntz_oq_no_go.py",
)
for s in scripts:
    print("\n#####",s,"#####",flush=True)
    run(["python",str(root/s)],check=True)
print("\nALL EXACT CHECKS PASSED")
