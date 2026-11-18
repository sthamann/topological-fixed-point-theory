"""Static AST inventory of every file named by the TFPT script registry.

Parsing is not execution: extracted check strings can describe assumptions,
and neither their truth nor their runtime coverage is established here.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tokenize


def source_segment(text, node):
    return ast.get_source_segment(text, node) or ''


class Extractor(ast.NodeVisitor):
    def __init__(self, text):
        self.text, self.scope = text, []
        self.imports, self.assignments, self.checks = [], [], []

    def record(self, node, kind):
        return {'kind': kind, 'line': node.lineno,
                'end_line': getattr(node, 'end_lineno', node.lineno),
                'scope': '.'.join(self.scope) or '<module>',
                'source': source_segment(self.text, node)}

    def visit_FunctionDef(self, node):
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_AsyncFunctionDef = visit_FunctionDef
    visit_ClassDef = visit_FunctionDef

    def visit_Import(self, node):
        self.imports.append(self.record(node, 'import'))

    visit_ImportFrom = visit_Import

    def visit_Assert(self, node):
        self.checks.append(self.record(node, 'assert'))
        self.generic_visit(node)

    def visit_Call(self, node):
        f = node.func
        name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else ''
        if name in {'check', 'require', 'need', 'verify', 'ensure'} or name.startswith('assert_'):
            self.checks.append(self.record(node, 'call:'+name))
        self.generic_visit(node)

    def visit_Assign(self, node):
        try:
            ast.literal_eval(node.value)
        except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
            pass
        else:
            self.assignments.append(self.record(node, 'literal_assignment'))
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        if node.value is not None:
            self.visit_Assign(node)


def inventory(root):
    if sys.flags.optimize >= 2:
        raise RuntimeError('Run without -OO: Python AST parsing can discard source docstrings')
    root = Path(root)
    registry = root/'verification/script_registry.csv'
    with registry.open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    records = []
    for row in rows:
        name = row['script']
        record = {'registry': row, 'source_path': None, 'status': 'unread'}
        if not re.fullmatch(r'[A-Za-z0-9_]+', name):
            record.update(status='invalid_registry_name'); records.append(record); continue
        path = root/'verification'/(name+'.py')
        record['source_path'] = str(path.relative_to(root))
        if not path.is_file():
            record['status'] = 'missing'; records.append(record); continue
        try:
            data = path.read_bytes()
            record['sha256'] = hashlib.sha256(data).hexdigest()
            encoding, _ = tokenize.detect_encoding(io.BytesIO(data).readline)
            text = data.decode(encoding)
            tree = ast.parse(text, filename=str(path))
            extractor = Extractor(text)
            extractor.visit(tree)
            doc = ast.get_docstring(tree, clean=False)
            doc_node = tree.body[0] if doc is not None else None
            record.update(status='parsed', docstring=doc,
                          docstring_line=doc_node.lineno if doc_node else None,
                          docstring_source=source_segment(text, doc_node) if doc_node else None,
                          checks=extractor.checks, imports=extractor.imports,
                          literal_assignments=extractor.assignments)
        except (OSError, UnicodeError, SyntaxError, ValueError, MemoryError, RecursionError) as error:
            record.update(status='parse_error', error=f'{type(error).__name__}: {error}')
        records.append(record)
    counts = {status: sum(x['status'] == status for x in records)
              for status in ('parsed', 'missing', 'parse_error', 'invalid_registry_name')}
    return {'schema': 'tfpt.registry-static-ast.v1', 'registry_entries': len(rows),
            'registry_sha256': hashlib.sha256(registry.read_bytes()).hexdigest(),
            'counts': counts, 'records': records,
            'scope': 'Static source extraction only; no verification script execution, proof validation, or claim promotion.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = inventory(args.root)
    with args.out.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
    print(json.dumps({'registry_entries': result['registry_entries'], 'counts': result['counts'],
                      'checks_extracted': sum(len(x.get('checks', [])) for x in result['records']),
                      'imports_extracted': sum(len(x.get('imports', [])) for x in result['records']),
                      'literal_assignments_extracted': sum(len(x.get('literal_assignments', [])) for x in result['records'])}))


if __name__ == '__main__':
    main()
