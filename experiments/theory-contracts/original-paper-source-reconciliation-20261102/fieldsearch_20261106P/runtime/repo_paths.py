"""Locate the actual TFPT repository without a machine-specific default."""
from pathlib import Path
import os

ENGINE_RELATIVE = Path('experiments/theory-contracts/seam-constructive-search-20261015/gpu_search.py')

def _is_repo(path):
    return (path/'verification/v252_full_finite_triple.py').is_file() and (path/ENGINE_RELATIVE).is_file()

def repo_root():
    explicit=os.environ.get('TFPT_REPO_ROOT')
    if explicit:
        root=Path(explicit).expanduser().resolve()
        if not _is_repo(root):
            raise ValueError('TFPT_REPO_ROOT does not identify the original TFPT repository')
        return root
    for start in (Path(__file__).resolve().parent,Path.cwd().resolve()):
        for root in (start,*start.parents):
            if _is_repo(root):return root
    raise ValueError('Run the installed package inside its original TFPT repository or set TFPT_REPO_ROOT')

def original_path(raw,root=None):
    path=Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else ((root or repo_root())/path).resolve()
