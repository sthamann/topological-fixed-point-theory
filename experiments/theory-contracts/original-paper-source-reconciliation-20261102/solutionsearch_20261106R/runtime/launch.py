"""Extend the existing field search with constructive source-program search."""
from pathlib import Path
import os
import runpy
import sys

HERE = Path(__file__).resolve().parents[1]


def main():
    args = sys.argv[1:]
    if '--mode' in args and args[args.index('--mode')+1] == 'solve':
        i = args.index('--mode')
        args = args[:i]+args[i+2:]
        sys.path.insert(0, str(HERE))
        from run_search import main as solve
        return solve(args)
    repo = Path(os.environ.get('TFPT_REPO_ROOT', '/Users/stefanhamann/Projekte/tfpt-theoryv4'))
    q = repo/'experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/launch.py'
    sys.path.insert(0, str(q.parent))
    runpy.run_path(str(q), run_name='__main__')


if __name__ == '__main__':
    raise SystemExit(main())
