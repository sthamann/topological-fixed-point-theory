"""Bounded API checks using a temporary copy of actual run records."""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from urllib.request import urlopen,Request
from urllib.error import HTTPError
import errno,hashlib,json,shutil,threading
import server

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'runs/constructive-local-search-v3'

def main():
    checks={}
    with TemporaryDirectory(prefix='tfpt-readonly-ui-')as temp:
        base=Path(temp);runs=base/'runs';run=runs/'actual-copy';run.mkdir(parents=True)
        for name in ('worker.json','status.json','launch.json','results.jsonl','model-rounds.jsonl'):
            if (SOURCE/name).exists():shutil.copy2(SOURCE/name,run/name)
        worker=json.loads((run/'worker.json').read_text());worker.update(status='RUNNING',pid=999999999)
        (run/'worker.json').write_text(json.dumps(worker))
        cache=server.FileCache();original=cache.read(run/'results.jsonl',True)
        with (run/'results.jsonl').open('ab')as f:f.write(b'{"unfinished":"\xc3')
        partial=cache.read(run/'results.jsonl',True)
        checks['incomplete_UTF8_JSONL_tail_preserves_complete_records']=partial['data']==original['data']and partial['trailing_ignored']
        checks['unchanged_files_reuse_mtime_cache']=cache.read(run/'results.jsonl',True)is partial
        with patch.object(server.os,'kill',side_effect=OSError(errno.EPERM,'permission')):
            checks['EPERM_is_unknown_not_dead']=server.process_alive(worker['pid'])is None
        with patch.object(server.os,'kill',side_effect=OSError(errno.ESRCH,'missing')):
            checks['ESRCH_is_dead']=server.process_alive(worker['pid'])is False
        with patch.object(server,'process_alive',return_value=False):
            state=server.Snapshots(runs,'actual-copy').state()
        checks['exited_worker_not_reported_running']=state['process']['status']=='stopped'and state['process']['alive']is False
        checks['real_four_values_and_proof_status_separated']=all(c['has_computed_values']and c['status']=='PARTIAL'for c in state['candidates'])and state['summary']['counts']['scoped_hits']==0
        checks['model_rounds_fragments_candidates_counted_separately']=state['summary']['counts']['program_fragments']>state['summary']['counts']['whole_candidates']and state['summary']['counts']['model_rounds']>=1
        for event in state['events']:
            if event['kind']=='repeat':checks['repeat_narrative_from_actual_log']=True
        actual_goal=next(r for r in original['data']if r.get('output_type')=='WholeSourceCandidate')
        scripted=server.round_event({'evaluation':actual_goal,'proposal':{'hypothesis':'<script>literal model text</script>'}},999)
        checks['model_text_remains_plain_data']=scripted['hypothesis']=='<script>literal model text</script>'and scripted['kind']=='kernelcomputed'
        static=base/'static';static.mkdir();(static/'index.html').write_text('<!doctype html><title>UI API control</title>')
        (static/'style.css').write_text('body{}');(static/'app.js').write_text('"use strict";')
        (static/'fonts').mkdir();(static/'fonts/manrope.ttf').write_bytes(b'test-only-font-route')
        outside=base/'outside';outside.mkdir();(outside/'worker.json').write_text('{"secret":"not-served"}')
        (runs/'linked-outside').symlink_to(outside,target_is_directory=True)
        before={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in run.iterdir()if p.is_file()}
        http=server.ReadOnlyServer(('127.0.0.1',0),runs,'actual-copy',static)
        thread=threading.Thread(target=http.serve_forever,daemon=True);thread.start()
        origin='http://127.0.0.1:'+str(http.server_port)
        def request(path,method='GET',headers=None):
            try:
                with urlopen(Request(origin+path,method=method,headers=headers or{}),timeout=3)as response:
                    return response.status,response.headers,response.read()
            except HTTPError as e:return e.code,e.headers,e.read()
        try:
            code,headers,body=request('/api/state?run=actual-copy');actual=json.loads(body)
            checks['actual_state_endpoint']=code==200 and len(actual['candidates'])==len(state['candidates'])and headers['Content-Type'].startswith('application/json')
            checks['runs_endpoint_bounded_to_allowed_root']=request('/api/runs')[0]==200 and 'linked-outside'not in {r['id']for r in json.loads(request('/api/runs')[2])['runs']}
            checks['static_files_and_font_whitelisted']=request('/')[0]==200 and request('/fonts/manrope.ttf')[1]['Content-Type']=='font/ttf'
            checks['arbitrary_file_and_listing_refused']=all(request(path)[0]==404 for path in ('/server.py','/fonts/','/../worker.json','/results.jsonl'))
            checks['path_traversal_and_symlink_refused']=request('/api/state?run=..%2Foutside')[0]==400 and request('/api/state?run=linked-outside')[0]==404
            checks['mutating_methods_refused']=request('/api/state',method='POST')[0]==405
            checks['foreign_host_refused']=request('/api/state',headers={'Host':'foreign.example'})[0]==403
            after={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in run.iterdir()if p.is_file()}
            checks['API_does_not_mutate_search_files']=before==after
        finally:http.shutdown();http.server_close();thread.join(timeout=2)
    result={'passed':all(checks.values()),'checks':checks,'scope':'Read-only API method controls over a temporary copy of actual search run data. No search process was modified.',
            'server_sha256':hashlib.sha256((HERE/'server.py').read_bytes()).hexdigest()}
    (HERE/'server-controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));assert result['passed']
if __name__=='__main__':main()
