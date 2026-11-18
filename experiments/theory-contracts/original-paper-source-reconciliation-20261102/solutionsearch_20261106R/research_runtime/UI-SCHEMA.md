# Original source research runtime interface

Launch: `python research_runtime/launch.py --repo CANONICAL_TFPT --out RUN --rounds 0`.
Zero rounds means unbounded. Repeating the command resumes completed proposal and
review calls from RUN. No model is invoked by importing these modules or controls.

The installed CLI rejects simultaneous `--sandbox` and `default_permissions`.
The runtime uses its concrete named permissions profile with restricted scientific
reads, only the current scratch writable, and network disabled. This provides
the requested workspace writing behavior with narrower reads. Model calls use
`--ignore-user-config`, `--ignore-rules`, `--ephemeral`, OpenAI provider and
`gpt-6.1-sol`; approval is never escalated. Existing CLI login remains unchanged.
Explanatory response fields and the initial task are German. Exact source
quotations and technical identifiers retain their original spelling.
An optional `RUN/operator-notes.txt` is advertised by path before each new role
call. It is read as operator steering within existing scientific scope and
permissions, creates no gate, and does not interrupt a running call. A running
controller loaded before this feature needs an ordinary later resume to advertise
the notes automatically; existing work should be preserved until then.
Calls have no mandatory elapsed time limit. The backend accepts an optional
timeout, where null or zero leaves a healthy process running. The profile
explicitly denies `/opt/homebrew/bin/ollama` and its resolved target. Scientific
instructions also prohibit Ollama, Qwen and local model servers, including
indirect launches through older scripts.

`RUN/worker.json` is atomically replaced during real events:

* `kind`: `original_source_research`
* `pid`, `status`: running, stopped, completed, error
* `started_utc`, `last_update_utc`, `requested_rounds`, `finished_rounds`
* `phase`: `{id,label,description,current,total,started_at}`. Current and total
  are actual completed calls/requested calls, with total null for unbounded work.
* `active_question`: exact current research task
* `backend_identity`: `{provider:'openai_codex',model:'gpt-6.1-sol',receipt_verified:bool}`
* `backend_receipt`: real most recent completed CLI response receipt, initially null.
  A CLI completion confirms a response on the pinned requested model, not a
  server supplied independent model attestation or physical proof.
* `knowledge_context`: manifest including revision, full local context/history paths
* `trial_summary`: `{rounds,model_calls,tool_events,completed_commands,failed_commands}`
* `active_call`: role, call_dir, tool counts, latest event type; optional error.
* `physics_admission`: always `NOT_ADMITTED_BY_RUNTIME`.

`RUN/research-rounds.jsonl` has one durable row per completed proposal and review
pair, with `round_id` (`research-000001`), `question`, `started_utc`, `completed_utc`,
`proposal`, `review` (full structured model responses), `proposal_call`,
`review_call` (raw log/artifact paths, usage, receipts and tool evidence),
`actual_result` (host observed command exits and evidence paths),
`proposed_status`, `critic_verdict`, `next_task`, `physics_admission`.
Critic claims remain AI interpretations even when both models agree.

`RUN/rounds/ROUND/round-state.json` persists after each role for resume.
`RUN/rounds/ROUND/{proposer,critic}-ATTEMPT/` contains `scratch/`,
`stdout.jsonl`, `stderr.log`, `last-message.json`, `cli-call.json`, `metrics.json`,
and host extracted `tool-events.jsonl`. A terminated partial attempt is retained;
resume makes a fresh attempt and never represents its incomplete output as done.
All prior rounds and finite programs remain accessible without a last N window.
An OS termination preserves the partial attempt and stops only the runtime's
own process group. Provider errors are saved without trying another model.
`RUN/context/VERSION/` contains full knowledge context, full previous finite trial
history, source index, review corrections, and a compact manifest. Context is
reloaded when Atlas or live plan bytes change. Sources are read only for calls.

No percentages, new scientific ledger claims, or automatic physics PASS exist.

Targeted controls: `python3 -B research_runtime/check_controls.py`. These use an
explicitly named fake control backend and do not contact a model. The local
configuration control accepted the installed CLI configuration without a model
call. It does not attest that a scientific calculation has run.
