"""Construct, calculate and refute source PROGRAMS using original native code.

OpenAI Sol proposals are bounded typed JSON ASTs. No model Python, commands,
module imports or self-attested proof results are executed.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import os
import signal
import sys
import time
import urllib.request

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from catalogue import inputs, operations, starting_programs
from completion import OriginalCompletion
from search.synthesis import Synthesizer, evidence_json, ast_hash, digest
from source_recovery.native_constructors import NativeInputMissing


def write(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(evidence_json(value), indent=2, ensure_ascii=False)+'\n')
    temp.replace(path)


def evidence(repo):
    from knowledge_context import build_knowledge_context
    return build_knowledge_context(repo)


def model_backend(repo, model='gpt-6.1-sol', log_dir=None):
    from sol_backend import SolBackend
    return SolBackend(log_dir=log_dir, model=model)


def flat_program(ast):
    instructions = []
    def visit(node):
        args = [visit(a) for a in node['args']]
        instructions.append({'op': node['op'], 'id': node['id'], 'args': args})
        return len(instructions)-1
    visit(ast)
    return instructions


def decode_program(program):
    if not isinstance(program, list) or not 1 <= len(program) <= 64:
        raise ValueError('program must be a nonempty bounded instruction list')
    values = []
    for i, node in enumerate(program):
        if not isinstance(node, dict) or set(node) != {'op', 'id', 'args'}:
            raise ValueError('instructions contain only op, id, args')
        refs = node['args']
        if not isinstance(refs, list) or any(type(j) is not int or not 0 <= j < i for j in refs):
            raise ValueError('instruction arguments must reference earlier instructions')
        values.append({'op': node['op'], 'id': node['id'], 'args': [values[j] for j in refs]})
    return values[-1]


def source_signature(value):
    if not isinstance(value, dict) or 'matrix' not in value or 'dictionary' not in value:
        return None
    # Whole symbolic source, NOT only its four zero-source jets. Higher source
    # powers and time dependence must not collapse to the same signature.
    d = value['dictionary']
    payload = {'schema':'whole-source-family-v2', 'kernel':value['matrix'],
               'source_symbols':value['source_symbols'],
               'P':d['P'], 'cap':d['cap'], 'reference':d['reference_scalar'],
               'outside':d['outer_scalar'],
               **{k:value[k] for k in ('time_symbols','mark_operators','domain_operators')if k in value}}
    return digest(evidence_json(payload))


def main(argv=None):
    from campaign import main as run_campaign
    return run_campaign(argv)


if __name__ == '__main__':
    raise SystemExit(main())
