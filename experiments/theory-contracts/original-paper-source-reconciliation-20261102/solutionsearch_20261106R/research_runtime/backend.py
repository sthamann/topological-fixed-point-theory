"""Fixed Sol CLI with live, host-observed tool evidence and restricted reads."""
from datetime import datetime, timezone
from pathlib import Path
import json
import os
import selectors
import signal
import subprocess
import time

from context import write_json, digest
from prompts import response_schema
from sol_backend import DEFAULT_CODEX, MODEL, redact_credentials


def utc():
    return datetime.now(timezone.utc).isoformat()


class ProviderError(RuntimeError):
    pass


def command(folder, manifest, *, model=MODEL):
    if model != 'gpt-6.1-sol':
        raise ValueError('Only gpt-6.1-sol is authorized; no fallback exists')
    scratch = folder / 'scratch'
    # The installed CLI's named profile has workspace-write behavior with
    # restricted reads. --sandbox and default_permissions are mutually exclusive.
    fs = {':minimal': 'read', str(folder.parent.parent.parent): 'read', str(scratch): 'write'}
    roots = list(map(Path, manifest['scientific_readonly_roots']))
    for path in manifest['scientific_readonly_roots']:
        fs[path] = 'read'
    for path in manifest['exact_source_files']:
        if not any(Path(path).is_relative_to(root) for root in roots):
            fs[path] = 'read'
    for path in manifest['scientific_readonly_roots']:
        for p in Path(path).rglob('*'):
            if p.name in {'.env', '.git', '.codex', '.ssh', 'auth.json', 'credentials.json'} or p.name.startswith('.env.'):
                fs[str(p)] = 'deny'
    # Python/Lean installation libraries only, never their user configuration.
    for p in ['/opt/homebrew/bin', '/opt/homebrew/opt', '/opt/homebrew/Cellar', '/opt/homebrew/lib',
              '/opt/homebrew/Frameworks', '/Library/Developer/CommandLineTools', '/Applications/Xcode.app/Contents/Developer',
              str(Path.home() / '.elan/bin'), str(Path.home() / '.elan/toolchains')]:
        if Path(p).exists():
            fs[p] = 'read'
    # Explicit user stop: even trusted legacy scripts must not start a local
    # model server through an otherwise readable package installation.
    ollama = Path('/opt/homebrew/bin/ollama')
    fs[str(ollama)] = 'deny'
    fs[str(ollama.resolve())] = 'deny'
    profile = '{filesystem={' + ','.join(json.dumps(k) + '=' + json.dumps(v) for k, v in fs.items()) + '},network={enabled=false}}'
    overrides = {
        'default_permissions': 'tfpt_research', 'permissions.tfpt_research': None,
        'approval_policy': 'never', 'model_provider': 'openai', 'mcp_servers': {},
        'web_search': 'disabled', 'project_doc_max_bytes': 0,
        'features.apps': False, 'features.hooks': False, 'features.memories': False,
        'features.multi_agent': False, 'features.multi_agent_v2': False,
        'shell_environment_policy.inherit': 'none',
        'shell_environment_policy.set': {'PATH': str(Path(manifest['python']).parent) + ':/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:' + str(Path.home() / '.elan/bin'),
            'TMPDIR': str(scratch / 'tmp'), 'PYTHONDONTWRITEBYTECODE': '1',
            'XDG_CACHE_HOME': str(scratch / 'cache')},
        'allow_login_shell': False,
    }
    argv = [DEFAULT_CODEX, '--no-daemon', '-a', 'never', 'exec', '--model', MODEL,
            '--ignore-user-config', '--ignore-rules', '--ephemeral', '--skip-git-repo-check',
            '--json', '-C', str(scratch), '--output-schema', str(folder / 'response-schema.json'),
            '--output-last-message', str(folder / 'last-message.json')]
    for key, value in overrides.items():
        argv += ['-c', key + '=' + (profile if value is None else
                 toml_value(value))]
    return argv + ['-'], {'filesystem': fs, 'network': {'enabled': False}}


def toml_value(value):
    if isinstance(value, dict):
        return '{' + ','.join(json.dumps(k) + '=' + toml_value(v) for k, v in value.items()) + '}'
    return json.dumps(value)


class Backend:
    def __init__(self, stop, update, timeout=None):
        self.stop, self.update, self.timeout = stop, update, timeout
        self.child = None

    def terminate(self):
        if self.child and self.child.poll() is None:
            os.killpg(self.child.pid, signal.SIGTERM)

    def __call__(self, role, folder, manifest, instruction):
        folder.mkdir(parents=True)
        scratch = folder / 'scratch'
        scratch.mkdir(); (scratch / 'tmp').mkdir()
        schema = response_schema(role)
        write_json(folder / 'response-schema.json', schema)
        (folder / 'prompt.txt').write_text(instruction)
        argv, policy = command(folder, manifest)
        write_json(folder / 'cli-call.json', {'argv': argv, 'role': role,
            'provider': 'openai_codex', 'model': MODEL, 'sandbox_policy': policy,
            'prompt_sha256': digest(folder / 'prompt.txt'), 'started_at': utc()})
        counters = {'tool_events': 0, 'completed_commands': 0, 'failed_commands': 0}
        usage, completed, commands, seen_tools = {}, False, [], set()
        started = time.monotonic()
        # CLI receives existing login through its normal session. No auth file
        # is opened by this host, and the child tool profile cannot read it.
        env = {k: v for k, v in os.environ.items() if k in {
            'HOME', 'PATH', 'USER', 'LOGNAME', 'SHELL', 'TMPDIR', 'LANG', 'LC_ALL'}}
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        try:
            self.child = subprocess.Popen(argv, cwd=scratch, stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                bufsize=1, env=env, start_new_session=True)
            self.child.stdin.write(instruction); self.child.stdin.close()
            selector = selectors.DefaultSelector()
            selector.register(self.child.stdout, selectors.EVENT_READ, 'stdout')
            selector.register(self.child.stderr, selectors.EVENT_READ, 'stderr')
            with (folder / 'stdout.jsonl').open('w') as output, (folder / 'stderr.log').open('w') as errors, (folder / 'tool-events.jsonl').open('w') as tool_log:
                stopped_at = None
                while selector.get_map():
                    if self.stop() or (self.timeout and time.monotonic() - started > self.timeout):
                        if stopped_at is None:
                            self.terminate(); stopped_at = time.monotonic()
                        elif time.monotonic() - stopped_at > 5 and self.child.poll() is None:
                            os.killpg(self.child.pid, signal.SIGKILL)
                    for key, _ in selector.select(0.5):
                        line = key.fileobj.readline()
                        if not line:
                            selector.unregister(key.fileobj); continue
                        line = redact_credentials(line)
                        stream = output if key.data == 'stdout' else errors
                        stream.write(line); stream.flush()
                        if key.data != 'stdout':
                            continue
                        try:
                            event = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        item = event.get('item') or {}
                        typ = item.get('type')
                        if event.get('type') == 'turn.completed':
                            completed = True; usage = event.get('usage') or {}
                        if event.get('type') in {'model.rerouted', 'model_reroute'} or (event.get('from_model') and event.get('to_model')):
                            raise ProviderError('Provider attempted model reroute; no fallback allowed')
                        if typ and typ not in {'agent_message', 'reasoning', 'plan'}:
                            tool_log.write(json.dumps(event, ensure_ascii=False) + '\n'); tool_log.flush()
                            tool_id = item.get('id') or json.dumps(item, sort_keys=True)
                            if tool_id not in seen_tools:
                                seen_tools.add(tool_id); counters['tool_events'] += 1
                            if typ == 'command_execution' and event.get('type') == 'item.completed':
                                counters['completed_commands'] += 1
                                if item.get('exit_code') not in (0, None):
                                    counters['failed_commands'] += 1
                                commands.append({'command': item.get('command'), 'exit_code': item.get('exit_code'),
                                    'output': item.get('aggregated_output'), 'evidence': str(folder / 'tool-events.jsonl')})
                        self.update({'role': role, 'call_dir': str(folder), **counters,
                                     'latest_event_type': event.get('type'), 'last_event_utc': utc()})
            code = self.child.wait()
            metrics = {'role': role, 'call_dir': str(folder), 'exit_code': code,
                'completed': completed, 'elapsed_seconds': time.monotonic() - started,
                'usage': usage, **counters, 'commands': commands}
            write_json(folder / 'metrics.json', metrics)
            if self.stop():
                raise InterruptedError('Stopped; partial call preserved')
            if code or not completed or not (folder / 'last-message.json').exists():
                raise ProviderError('OpenAI Sol call failed; inspect ' + str(folder))
            raw = (folder / 'last-message.json').read_text()
            safe = redact_credentials(raw)
            (folder / 'last-message.json').write_text(safe)
            if safe != raw:
                raise ProviderError('Credential-shaped output redacted and rejected')
            result = json.loads(raw)
            if not isinstance(result, dict) or set(result) != set(schema['properties']):
                raise ProviderError('Structured response does not match required schema')
            for key, spec in schema['properties'].items():
                value = result[key]
                if (spec['type'] == 'string' and not isinstance(value, str)) or (spec['type'] == 'array' and
                    (not isinstance(value, list) or any(not isinstance(x, str) for x in value))) or ('enum' in spec and value not in spec['enum']):
                    raise ProviderError('Invalid structured field: ' + key)
            receipt = {'provider': 'openai_codex', 'model': MODEL, 'verified': True,
                'verification_kind': 'pinned_cli_completed_response', 'server_model_attested': False,
                'call_dir': str(folder), 'exit_code': code, 'completed_at': utc()}
            metrics['receipt'] = receipt
            write_json(folder / 'metrics.json', metrics)
            return result, metrics
        except BaseException as exc:
            self.terminate()
            write_json(folder / 'failure.json', {'type': type(exc).__name__, 'error': str(exc), 'at': utc()})
            raise
        finally:
            if self.child and self.child.poll() is None:
                self.terminate()
                try:
                    self.child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(self.child.pid, signal.SIGKILL); self.child.wait()
            self.child = None
