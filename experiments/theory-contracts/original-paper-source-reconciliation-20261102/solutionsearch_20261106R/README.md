Aktueller Stand: `SOURCE-CANDIDATE-STEERING-20261106R5.md`. Der Quellenkandidatenlauf enthält jetzt feste Zielslots, direkte Atlasgegenbefunde und gesonderte Zielfortschrittsbewertung. Alte endliche Kampagnen bleiben historisch.

# TFPT Originalquellenforschung

Der aktuelle automatische Lauf verwendet `research_runtime/launch.py`, sowohl für Herleitung als auch für unabhängige Gegenprüfung ausschließlich mit `gpt-6.1-sol`. Er rechnet, kritisiert und verfolgt die daraus entstehende nächste Forschungsfrage fortlaufend. Der ganze Verlauf sowie tatsächliche Werkzeugergebnisse bleiben erhalten.

Siehe `AUTONOMOUS-ORIGINAL-SOURCE-RESEARCH-20261106R4.md` für aktuellen Aufruf, Wiederaufnahme, Quellenkontext und wissenschaftlichen Geltungsbereich. Der Monitor wählt standardmäßig `sol-original-source-research`.

Die folgende Dokumentation beschreibt ausschließlich die ältere endliche Hypothesensuche. Sie wird als historische Vergleichsmethode erhalten. Ihr Standardstart endet bei fehlendem physikalischem Auswertungspfad mit `NO_COMPLETION_ROUTE`.

---

# Constructive common-source program search

This extends the existing TFPT field search with `runtime/launch.py --mode solve`.
Candidates are executable weighted even-kernel programs, including the scalar
body, higher interactions, whole Berezin sewing and marked source derivatives.
They are not a list of precomputed four values. The original native constructors
are called to calculate the proposed whole kernel and its four values.

The present search is **a constructive hypothesis search**, not a completed
automatic TFPT proof solver. Its finite retained interfaces, coefficient ring,
cost bounds, chart choices and outer dictionaries remain explicit. The active
grammar does not yet construct the missing original physical identity/time/
continuum proofs. The default physical verifier has no trusted checkers closing
those proofs. Thus a full TFPT hit is unreachable in this fixed grammar alone.
Model proposals additionally name the missing derivation and propose an origin
argument for separate original-source review. These arguments are not admitted
as proofs just because a model supplies them. The user's stronger completion
objective remains open.

## Run

From the TFPT repository, replace `PACKAGE` by this installed package path:

```sh
experiments/tfpt-discovery/.venv/bin/python -B PACKAGE/runtime/launch.py --mode solve --repo "$PWD" --out RUN_DIRECTORY --rounds 24 --target tfpt_complete
```

The proposal backend is exclusively **OpenAI Sol `gpt-6.1-sol`**, through the
existing authenticated Codex session. The user explicitly approved the selected
internal TFPT package and actual candidate counterfeedback on 6 November 2026.
There is no Ollama or Qwen fallback. The Ollama application and server were
stopped at the user's request. Credentials are not read by the search client.

Each round requests up to eight model proposals and evaluates a direct batch
of up to 64 whole programs. The typed registry executes only listed native
operations; model commands and self-attested PASS flags are never executed.
The systematic grammar changes action terms, U/Y/UY insertions, whole-kernel
corrections, coupled vertices and ordered/reflected whole Berezin sewing.
Whole programs, unique source calculations, repeated programs and model calls
are recorded separately. The first 256 constructions were statically checked;
16 selected programs were actually executed and their four values independently
recalculated. Those tests establish finite algebra and software behavior.

`--without-model` evaluates the direct systematic construction batches alone.
`--prepare-sol-payload` writes an exact local preview without calling OpenAI.
Numerical integration and diagnostic modes still delegate to unchanged Q.
Expression cost and proposal probabilities are never physical source weights.

The read-only live monitor runs from `ui/server.py --runs-root RUNS --port 8768
--default-run openai-sol-source-search`. It listens only on localhost. Selection,
filters, freezing the view and downloading a snapshot do not stop the search.

## Read results

* `worker.json`: actual process, model identity, progress and capability limit.
* `results.jsonl`: executed ASTs, exact four values, failed tests and open proofs.
* `model-rounds.jsonl`: model proposals, proposed origin arguments and actual feedback.
* `status.json`: scoped results and declared grammar coverage.
* `best-unresolved.json`: current smallest unresolved construction.
* `original-evidence.json`: preserved Atlas revision and relevant scoped results.

Partial results are retained. One missing input does not abort the search.
Accepted scoped hits require every declared obligation to pass the trusted
verifier. Existing original Q checks receive actual current arrays and data;
an attached independent packet must also match the candidate's word embedding,
all four forms and values, P, cap and absolute reference factors.

Resume keys bind actual inputs, code, closures, verifier state and explicit
dependencies. A changed verifier or reference cannot reuse an old successful
result. Checkpoints never convert hypotheses or finite controls into physical
proofs.

## Checks

The native constructor has eight exact controls and four targeted metadata
controls. The synthesis core has 19 method controls on the old finite v330
source, including an actual scoped success and a corrupted candidate failure.
These are software/finite-algebra checks, not TFPT completion evidence. The
independent review and the fixed cache/coherence counterexamples are in `review/`.

No physical ledger, paper proof status or empirical prediction was promoted.
