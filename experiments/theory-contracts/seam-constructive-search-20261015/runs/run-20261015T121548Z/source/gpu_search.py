#!/usr/bin/env python3
"""Checkpointed constructive genetic search; finite diagnostics never imply proof.

Plugin: sample(rng,n), mutate(rng,parent), canonicalize(genotype),
evaluate(candidates,device). Optional valid(genotype), crossover(rng,a,b).
Evaluation rows contain objectives (finite floats, all minimized), feasible,
descriptor (a discrete diversity cell), and JSON details. No model inference.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import random
import sqlite3
import sys
import time

import numpy as np


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def atomic_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as stream:
        stream.write(encode(value))
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def load_problem(path):
    path = Path(path).resolve()
    spec = importlib.util.spec_from_file_location("constructive_problem", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    for method in ("sample", "mutate", "canonicalize", "evaluate"):
        if not callable(getattr(module, method, None)):
            raise ValueError(f"Plugin missing {method}")
    return module


def manifest(problem):
    paths = [Path(__file__), Path(problem.__file__)]
    paths += [Path(p) if Path(p).is_absolute() else Path(problem.__file__).parent / p
              for p in getattr(problem, "SOURCE_FILES", [])]
    files = [{"name": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
             for p in paths]
    # Absolute local/remote deployment paths are deliberately not hashed.
    return {"files": files, "sha256": hashlib.sha256(encode(files).encode()).hexdigest()}


class Archive:
    """One WAL database shared by independent islands, keyed by source + genotype."""
    def __init__(self, path, source, owner):
        self.source, self.owner = source, owner
        self.db = sqlite3.connect(path, timeout=60)
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA busy_timeout=60000")
        self.db.execute("""CREATE TABLE IF NOT EXISTS candidates
            (source TEXT, key TEXT, genotype TEXT, result TEXT, owner TEXT,
             created REAL, PRIMARY KEY(source,key))""")
        # An exclusive island lock is held by the caller. Only its abandoned
        # reservations may be reclaimed; other running islands are untouched.
        self.db.execute("DELETE FROM candidates WHERE source=? AND owner=? AND result IS NULL",
                        (source, owner))
        self.db.commit()

    def key(self, genotype):
        return hashlib.sha256(encode(genotype).encode()).hexdigest()

    def reserve(self, genotype):
        key = self.key(genotype)
        with self.db:
            cur = self.db.execute("INSERT OR IGNORE INTO candidates VALUES (?,?,?,?,?,?)",
                                  (self.source, key, encode(genotype), None, self.owner, time.time()))
        return key if cur.rowcount else None

    def finish(self, rows):
        with self.db:
            for row in rows:
                self.db.execute("UPDATE candidates SET result=? WHERE source=? AND key=? AND owner=?",
                                (encode(row["result"]), self.source, row["key"], self.owner))

    def release_pending(self):
        with self.db:
            self.db.execute("DELETE FROM candidates WHERE source=? AND owner=? AND result IS NULL",
                            (self.source, self.owner))

    def migrants(self, rng, limit):
        # Pull a bounded recent pool, then preserve objective and cell diversity.
        rows = self.db.execute("SELECT key,genotype,result FROM candidates WHERE source=? "
                               "AND result IS NOT NULL ORDER BY created DESC LIMIT ?",
                               (self.source, max(limit * 32, 256))).fetchall()
        return select([{"key": k, "genotype": json.loads(g), "result": json.loads(r)}
                       for k, g, r in rows], limit, rng)

    def counts(self):
        evaluated, pending = self.db.execute(
            "SELECT COUNT(result),COUNT(*)-COUNT(result) FROM candidates WHERE source=?",
            (self.source,)).fetchone()
        return {"evaluated": evaluated, "pending": pending}


def dominates(a, b):
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def pareto_indices(scores, eligible):
    """Exact finite-vector dominance, blocked to bound the comparison memory."""
    indices = np.flatnonzero(eligible)
    pool = scores[indices]
    front = []
    for offset in range(0, len(indices), 256):
        block = pool[offset:offset + 256]
        no_worse = np.ones((len(pool), len(block)), dtype=bool)
        strictly_better = np.zeros_like(no_worse)
        for dimension in range(scores.shape[1]):
            a, b = pool[:, dimension, None], block[None, :, dimension]
            no_worse &= a <= b
            strictly_better |= a < b
        undominated = ~(no_worse & strictly_better).any(axis=0)
        front.extend(indices[offset:offset + len(block)][undominated].tolist())
    return front


def select(records, count, rng):
    """Feasibility, Pareto front, then archive cells and objective extremes."""
    unique = {row["key"]: row for row in records}
    rows = list(unique.values())
    if len(rows) <= count:
        return rows
    scores = np.asarray([r["result"]["objectives"] for r in rows], dtype=float)
    feasible = np.asarray([r["result"]["feasible"] for r in rows], dtype=bool)
    best_class = feasible if feasible.any() else np.ones(len(rows), dtype=bool)
    front = pareto_indices(scores, best_class)
    span = np.ptp(scores, axis=0)
    norm = (scores - np.min(scores, axis=0)) / np.where(span > 0, span, 1)
    weights = rng.dirichlet(np.ones(scores.shape[1]))
    scalar = norm @ weights
    # Every discrete cell retains its best feasible member when capacity permits.
    order = sorted(range(len(rows)), key=lambda i: (not feasible[i], scalar[i], rows[i]["key"]))
    cell_best = {}
    for i in order:
        cell_best.setdefault(encode(rows[i]["result"]["descriptor"]), i)
    extremes = [min(range(len(rows)), key=lambda i: (not feasible[i], scores[i, j]))
                for j in range(scores.shape[1])]
    chosen = []
    # Reserve half the population for diverse cells, then Pareto candidates.
    diversity = list(cell_best.values())
    rng.shuffle(diversity)
    for i in extremes + diversity[:max(1, count // 2)] + sorted(front, key=lambda i: scalar[i]) + order:
        if i not in chosen:
            chosen.append(i)
        if len(chosen) == count:
            break
    return [rows[i] for i in chosen]


def validate_result(result, dimensions=None):
    objectives = result.get("objectives")
    if not isinstance(objectives, (list, tuple)) or not objectives:
        raise ValueError("Each evaluation needs a nonempty objectives vector")
    objectives = [float(x) for x in objectives]
    if not all(math.isfinite(x) for x in objectives):
        raise ValueError("Nonfinite objectives cannot enter the archive")
    if dimensions is not None and len(objectives) != dimensions:
        raise ValueError("Objective vector dimensions changed")
    if type(result.get("feasible")) is not bool:
        raise ValueError("feasible must be an explicit boolean")
    if "descriptor" not in result:
        raise ValueError("Missing diversity descriptor")
    result = {**result, "objectives": objectives}
    encode(result)
    return result


def torch_setup(device, seed):
    if device == "cpu":
        return None
    import torch
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError(f"Requested GPU unavailable: {device}")
    torch.cuda.set_device(device)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    return torch


def save_checkpoint(path, config, source, rng, generation, population, total, reason, torch):
    state = {"schema": 1, "config": config, "source": source, "rng": rng.bit_generator.state,
             "python_rng": random.getstate(), "generation": generation, "population": population,
             "total_evaluated_here": total, "reason": reason, "time": time.time()}
    if torch:
        state["torch_cpu_rng"] = torch.get_rng_state().tolist()
        state["torch_cuda_rng"] = [s.tolist() for s in torch.cuda.get_rng_state_all()]
    atomic_json(path, state)


def nested_tuple(value):
    return tuple(nested_tuple(x) for x in value) if isinstance(value, list) else value


def run(args, problem=None):
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    lock = (out / f"island-{args.island}.lock").open("w")
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        problem = problem or load_problem(args.problem)
        result = _run(args, out, problem)
        if args.exhaustive_after and result["reason"] in ("generation_limit", "no_novel_candidates"):
            result["classification"] = classify(args, out, problem)
        return result
    finally:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


def classify(args, out, problem):
    """Classify the declared finite shard once, retaining the genetic evidence."""
    if not callable(getattr(problem, "enumerate_candidates", None)):
        raise ValueError("Exhaustive classification requires problem.enumerate_candidates(shard,shards)")
    shard, shards = args.island, args.exhaustive_shards
    if shards < 1 or not 0 <= shard < shards:
        raise ValueError("Exhaustive island must be in range(exhaustive_shards)")
    source = json.loads((out / f"island-{shard}.json").read_text())["source"]
    archive = Archive(out / "archive.sqlite3", source["sha256"], str(shard))
    archive.db.execute("""CREATE TABLE IF NOT EXISTS classification
        (source TEXT, shards INTEGER, key TEXT, shard INTEGER, feasible INTEGER, result_hash TEXT,
         PRIMARY KEY(source,shards,key))""")
    archive.db.commit()
    expected = getattr(problem, "ENUMERATION_SIZE", None)
    digest, seen, survivors, deferred = hashlib.sha256(), set(), [], []
    count = new_evaluations = cached_count = 0
    started = time.monotonic()
    report_path = out / f"island-{shard}-classification.json"
    status_path = out / f"island-{shard}-status.json"
    torch = torch_setup(args.device, args.seed)

    def record(rows):
        nonlocal count
        with archive.db:
            for row in rows:
                result_hash=hashlib.sha256(encode({k:row["result"][k] for k in ("objectives","feasible")}).encode()).hexdigest()
                archive.db.execute("INSERT OR IGNORE INTO classification VALUES (?,?,?,?,?,?)",
                                   (archive.source, shards, row["key"], shard, row["result"]["feasible"],result_hash))
                owner = archive.db.execute("SELECT shard FROM classification WHERE source=? AND shards=? AND key=?",
                                          (archive.source, shards, row["key"])).fetchone()[0]
                if owner != shard:
                    raise ValueError("Enumeration shards overlap: " + row["key"])
                if row["result"]["feasible"]:
                    survivors.append(row)
                count += 1

    def progress(stage, publish_status=True):
        report = {"stage": stage, "shard": shard, "shards": shards, "source": source,
                  "classification_count": count, "enumerated_count": len(seen),
                  "classification_expected": expected,
                  "classification_feasible": len(survivors),
                  "classification_new_evaluations": new_evaluations, "cached_count": cached_count,
                  "enumeration_key_digest": digest.hexdigest(), "survivors": survivors,
                  "elapsed_seconds": time.monotonic()-started,
                  "scope": "complete finite search filter coverage; exact proof and physical selection separate"}
        atomic_json(report_path, report)
        if publish_status:
            status = json.loads(status_path.read_text())
            atomic_json(status_path, {**status, **{k: v for k, v in report.items()
                                                  if k not in ("source", "survivors")},
                                     "reason": "running" if stage=="classifying" else stage})
        return report

    try:
        progress("classifying")
        iterator = iter(problem.enumerate_candidates(shard, shards))
        batches = 0
        dimensions = None
        while True:
            if (out / "STOP").exists() or (out / f"STOP-{shard}").exists():
                return progress("classification_stopped")
            batch = list(itertools.islice(iterator, args.batch))
            if not batch:
                break
            rows = []
            for genotype in batch:
                genotype = problem.canonicalize(genotype)
                if hasattr(problem, "valid") and not problem.valid(genotype):
                    raise ValueError("Enumeration emitted an invalid canonical genotype")
                key = archive.key(genotype)
                if key in seen:
                    raise ValueError("Enumeration repeated a canonical genotype: " + key)
                seen.add(key)
                digest.update((key+"\n").encode())
                rows.append({"key": key, "genotype": genotype})
            keys = [r["key"] for r in rows]
            placeholders = ",".join("?" for _ in keys)
            existing = dict(archive.db.execute(
                f"SELECT key,result FROM candidates WHERE source=? AND key IN ({placeholders})",
                [archive.source]+keys).fetchall())
            pending, complete = [], []
            with archive.db:
                for row in rows:
                    value = existing.get(row["key"])
                    if value is not None:
                        row["result"] = validate_result(json.loads(value), dimensions)
                        dimensions = len(row["result"]["objectives"])
                        complete.append(row)
                        cached_count += 1
                        continue
                    cur = archive.db.execute("INSERT OR IGNORE INTO candidates VALUES (?,?,?,?,?,?)",
                        (archive.source,row["key"],encode(row["genotype"]),None,archive.owner,time.time()))
                    (pending if cur.rowcount else deferred).append(row)
            if pending:
                results = problem.evaluate([r["genotype"] for r in pending], args.device)
                if torch:
                    torch.cuda.synchronize()
                if len(results) != len(pending):
                    raise ValueError("Evaluator returned wrong exhaustive row count")
                for row,result in zip(pending,results):
                    row["result"] = validate_result(result,dimensions)
                    dimensions = len(row["result"]["objectives"])
                archive.finish(pending)
                new_evaluations += len(pending)
                complete += pending
            record(complete)
            batches += 1
            if batches % 32 == 0:
                progress("classifying")
        # An island can enter classification while the other finishes its final
        # genetic batch. Reuse its results rather than evaluate those genes twice.
        deadline = time.monotonic()+10
        for row in deferred:
            while True:
                value = archive.db.execute("SELECT result FROM candidates WHERE source=? AND key=?",
                                          (archive.source,row["key"])).fetchone()
                if value and value[0] is not None:
                    row["result"] = validate_result(json.loads(value[0]),dimensions)
                    record([row])
                    cached_count += 1
                    break
                if time.monotonic() >= deadline:
                    raise RuntimeError("Other island has unfinished candidate reservations; resume to reuse completed results")
                time.sleep(.05)
        report = progress("classification_complete",publish_status=False)
        # Serialize the small union report; per-island reports remain independent.
        with (out / "classification-union.lock").open("w") as lock:
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
            reports = [json.loads(p.read_text()) for p in out.glob("island-*-classification.json")]
            reports = [r for r in reports if r["source"]==source and r["shards"]==shards]
            counts = archive.db.execute("SELECT shard,COUNT(*),SUM(feasible) FROM classification "
                                       "WHERE source=? AND shards=? GROUP BY shard",(archive.source,shards)).fetchall()
            union_count = sum(r[1] for r in counts)
            complete = len({r["shard"] for r in reports if r["stage"]=="classification_complete"})==shards
            union_digest = hashlib.sha256()
            result_digest = hashlib.sha256()
            for key,result_hash in archive.db.execute("SELECT key,result_hash FROM classification WHERE source=? AND shards=? ORDER BY key",
                                                      (archive.source,shards)):
                union_digest.update((key+"\n").encode())
                result_digest.update((key+" "+result_hash+"\n").encode())
            union = {"complete": complete and (expected is None or union_count==expected),
                     "source": source, "shards": shards, "expected_total": expected,
                     "classified_unique": union_count, "feasible_unique":sum(r[2] for r in counts),
                     "shard_counts": counts, "canonical_key_digest": union_digest.hexdigest(),
                     "classification_result_digest": result_digest.hexdigest(),
                     "result_digest_rule":"SHA256 of sorted key SPACE SHA256(canonical objectives+feasible JSON) NEWLINE",
                     "survivors":[s for r in reports for s in r["survivors"]],
                     "scope":"finite filter classification; no automatic theorem or physical completion"}
            atomic_json(out / "classification-union.json",union)
            if complete and expected is not None and union_count != expected:
                raise ValueError("Enumeration union cardinality disagrees with ENUMERATION_SIZE")
        return progress("classification_complete")
    except BaseException:
        archive.release_pending()
        progress("classification_error")
        raise
    finally:
        archive.db.close()


def _run(args, out, problem=None):
    problem = problem or load_problem(args.problem)
    source = manifest(problem)
    if args.source_manifest:
        packet = Path(args.source_manifest)
        source["external_packet_sha256"] = hashlib.sha256(packet.read_bytes()).hexdigest()
        source["sha256"] = hashlib.sha256(encode(source).encode()).hexdigest()
    run_label = getattr(problem, "RUN_LABEL", "research_candidate_search")
    config = {"island": args.island, "seed": args.seed, "population": args.population,
              "run_label": run_label,
              "batch": args.batch, "migration_every": args.migration_every}
    if args.population < 2 or args.batch < 1 or args.generations < 0 or args.migration_every < 1:
        raise ValueError("population >=2, batch/migration_every >=1, generations >=0 required")
    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)
    torch = torch_setup(args.device, args.seed)
    checkpoint = out / f"island-{args.island}.json"
    owner = str(args.island)
    population, generation, total = [], 0, 0
    if checkpoint.exists():
        state = json.loads(checkpoint.read_text())
        if state["config"] != config or state["source"] != source:
            raise ValueError("Resume rejected: source or search configuration changed; use a new output directory")
        rng.bit_generator.state = state["rng"]
        random.setstate(nested_tuple(state["python_rng"]))
        population, generation, total = state["population"], state["generation"], state["total_evaluated_here"]
        if torch and "torch_cpu_rng" in state:
            torch.set_rng_state(torch.tensor(state["torch_cpu_rng"], dtype=torch.uint8))
            torch.cuda.set_rng_state_all([torch.tensor(s, dtype=torch.uint8) for s in state["torch_cuda_rng"]])
    archive = Archive(out / "archive.sqlite3", source["sha256"], owner)
    log_path = out / f"island-{args.island}.jsonl"
    started = time.monotonic()
    summary = {"island": args.island, "generation": generation,
               "total_evaluated_here": total, "device": args.device, "run_label": run_label}
    reason = "generation_limit"
    try:
        while generation < args.generations:
            if args.hours is not None and time.monotonic() - started >= args.hours * 3600:
                reason = "duration_limit"
                break
            if (out / "STOP").exists() or (out / f"STOP-{args.island}").exists():
                reason = "stop_sentinel"
                break
            pending = []
            attempts = 0
            while len(pending) < args.batch and attempts < args.batch * 40:
                attempts += 1
                if not population or rng.random() < 0.15:
                    genotype = problem.sample(rng, 1)[0]
                else:
                    a = population[int(rng.integers(len(population)))]["genotype"]
                    if hasattr(problem, "crossover") and rng.random() < 0.25:
                        b = population[int(rng.integers(len(population)))]["genotype"]
                        genotype = problem.crossover(rng, copy.deepcopy(a), copy.deepcopy(b))
                    else:
                        genotype = problem.mutate(rng, copy.deepcopy(a))
                genotype = problem.canonicalize(genotype)
                if hasattr(problem, "valid") and not problem.valid(genotype):
                    continue
                key = archive.reserve(genotype)
                if key:
                    pending.append({"key": key, "genotype": genotype})
            if not pending:
                reason = "no_novel_candidates"
                break
            if torch:
                torch.cuda.synchronize()
                torch.cuda.reset_peak_memory_stats()
                gpu_start = torch.cuda.Event(enable_timing=True)
                gpu_end = torch.cuda.Event(enable_timing=True)
                gpu_start.record()
            start = time.perf_counter()
            results = problem.evaluate([p["genotype"] for p in pending], args.device)
            if torch:
                gpu_end.record()
                torch.cuda.synchronize()
            elapsed = time.perf_counter() - start
            if len(results) != len(pending):
                raise ValueError("Evaluator returned wrong number of rows")
            dimensions = len(population[0]["result"]["objectives"]) if population else None
            for row, result in zip(pending, results):
                row["result"] = validate_result(result, dimensions)
                dimensions = len(row["result"]["objectives"])
            archive.finish(pending)
            total += len(pending)
            generation += 1
            candidates = population + pending
            if generation % args.migration_every == 0:
                candidates += archive.migrants(rng, args.population)
            population = select(candidates, args.population, rng)
            summary = {"island": args.island, "generation": generation, "new_evaluations": len(pending),
                       "run_label": run_label,
                       "total_evaluated_here": total, "device": args.device,
                       "evaluation_seconds": elapsed, "archive": archive.counts(),
                       "invocation_elapsed_seconds": time.monotonic() - started,
                       "feasible_population": sum(r["result"]["feasible"] for r in population),
                       "diversity_cells": len({encode(r["result"]["descriptor"]) for r in population}),
                       "objective_minima": np.min([r["result"]["objectives"] for r in population], axis=0).tolist(),
                       "finite_diagnostic_only": True, "time": time.time()}
            if torch:
                properties = torch.cuda.get_device_properties(torch.cuda.current_device())
                summary.update(gpu_name=torch.cuda.get_device_name(),
                               gpu_uuid=str(getattr(properties, "uuid", "unavailable")),
                               cuda_event_milliseconds=gpu_start.elapsed_time(gpu_end),
                               peak_allocated_bytes=torch.cuda.max_memory_allocated())
            with log_path.open("a") as stream:
                stream.write(encode(summary) + "\n")
            print(encode(summary), flush=True)
            save_checkpoint(checkpoint, config, source, rng, generation, population, total, "running", torch)
            atomic_json(out / f"island-{args.island}-status.json", {**summary, "reason": "running"})
            atomic_json(out / f"island-{args.island}-archive.json",
                        {"source": source, "generation": generation, "candidates": population, "run_label": run_label,
                         "scope": "finite numerical candidates; no exact or physical success implied"})
    except BaseException:
        archive.release_pending()
        archive.db.close()
        # Keep the previous complete generation and RNG state, not a half batch.
        raise
    save_checkpoint(checkpoint, config, source, rng, generation, population, total, reason, torch)
    atomic_json(out / f"island-{args.island}-status.json",
                {**summary, "reason": "running" if args.exhaustive_after and reason in
                 ("generation_limit","no_novel_candidates") else reason, "archive": archive.counts()})
    atomic_json(out / f"island-{args.island}-archive.json",
                {"source": source, "generation": generation, "candidates": population, "run_label": run_label,
                 "scope": "finite numerical candidates; no exact or physical success implied"})
    archive.db.close()
    return {"generation": generation, "total_evaluated_here": total, "reason": reason,
            "checkpoint": str(checkpoint)}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--problem", default=str(Path(__file__).with_name("problem.py")))
    p.add_argument("--output", "--run-dir", required=True)
    p.add_argument("--island", type=int, default=0)
    p.add_argument("--seed", type=int, default=731)
    p.add_argument("--device", default="cuda:0")
    p.add_argument("--population", type=int, default=128)
    p.add_argument("--batch", "--batch-size", type=int, default=128)
    p.add_argument("--generations", type=int, default=100)
    p.add_argument("--migration-every", type=int, default=5)
    p.add_argument("--hours", type=float, default=None, help="Wall time for this invocation, checked between batches")
    p.add_argument("--source-manifest", help="Optional Atlas/source packet; its exact bytes bind the resume state")
    p.add_argument("--exhaustive-after", action="store_true", help="After the genetic generation limit, classify the finite shard")
    p.add_argument("--exhaustive-shards", type=int, default=1, help="Disjoint enumeration shards; shard index is --island")
    return p


if __name__ == "__main__":
    print(encode(run(parser().parse_args())), flush=True)
