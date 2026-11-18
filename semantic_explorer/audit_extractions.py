#!/usr/bin/env python3
"""Check annotation references against source packets; this is not a proof checker."""
import argparse
import json
from pathlib import Path
import time
from extract import Corpus,HERE,REPO,atomic,now,validate,sha
from translate import claims,fingerprint


def semantic_source_refs(value):
    refs=set()
    if isinstance(value,dict):
        source_refs=value.get('source_refs')
        if isinstance(source_refs,list):refs.update(x for x in source_refs if isinstance(x,str))
        for child in value.values():refs.update(semantic_source_refs(child))
    elif isinstance(value,list):
        for child in value:refs.update(semantic_source_refs(child))
    return refs


def manual_issue_still_applies(issue,semantic,claim_id):
    """Do not replay a stale self-reference finding after its source link was repaired."""
    if issue.get('code')!='self_upstream':return True
    return any(item.get('origin_claim_id')==claim_id for item in semantic.get('premises',[])) or \
        any(item.get('upstream_claim_id')==claim_id for item in semantic.get('input_audit',[])) or \
        any(claim_id in item.get('upstream_claim_ids',[]) for item in semantic.get('derivations',[])) or \
        any(item.get('target_claim_id')==claim_id for item in semantic.get('relationships',[]))


def source_freshness_signature(node,source_index,repo_root=REPO):
    """Cheap follow-loop key; exact SHA checks still happen when validation runs."""
    attrs=node.get('attrs',{})
    refs=semantic_source_refs(attrs.get('semantic',{}))|set(x for x in attrs.get('source_refs',[]) if isinstance(x,str))
    state=[]
    for sid in sorted(refs):
        source=source_index.get(sid)
        if not source:
            state.append((sid,'missing'))
            continue
        path=repo_root/source.get('path','')
        try:
            st=path.stat()
            stat=(st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns)
        except OSError:
            stat='unavailable'
        state.append((sid,source.get('path'),source.get('sha256'),stat))
    return sha(json.dumps(state,sort_keys=True,separators=(',',':')).encode())


def load_validation_indexes():
    """Load sources and typed Lean declarations from the same canonical graph."""
    path=HERE/'data/knowledge.json'
    try:graph=json.loads(path.read_text(encoding='utf-8'))
    except (OSError,ValueError):return {},{}
    sources={source['id']:source for source in graph.get('sources',[])
             if isinstance(source,dict) and isinstance(source.get('id'),str)}
    lean={node['id']:node for node in graph.get('nodes',[])
          if isinstance(node,dict) and isinstance(node.get('id'),str)
          and node.get('type') in {'lean_theorem','lean_declaration'}}
    return sources,lean


def load_source_index():
    return load_validation_indexes()[0]


def source_manifest_stamp():
    stamp=[]
    path=HERE/'data/knowledge.json'
    try:
        st=path.stat();stamp.append((str(path),st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns))
    except OSError:stamp.append((str(path),'unavailable'))
    return tuple(stamp)


def verified_proposal_sources(proposal, repo_root=REPO):
    """Return proposal source records only when their paths, hashes and ranges are valid."""
    hashes = {item.get('source_id'): item for item in proposal.get('source_hashes', [])
              if isinstance(item, dict) and isinstance(item.get('source_id'), str)}
    records = {}
    for source in proposal.get('sources', []):
        if not isinstance(source, dict):
            continue
        sid, rel, digest = source.get('source_id'), source.get('path'), source.get('sha256')
        start, end = source.get('line_start'), source.get('line_end')
        hashed = hashes.get(sid, {})
        if (not isinstance(sid, str) or not isinstance(rel, str) or not isinstance(digest, str)
                or not isinstance(start, int) or isinstance(start, bool) or start < 1
                or not isinstance(end, int) or isinstance(end, bool) or end < start
                or hashed.get('path') != rel or hashed.get('sha256') != digest):
            continue
        path = (repo_root / rel).resolve()
        if not path.is_relative_to(repo_root.resolve()) or not path.is_file():
            continue
        try:
            raw = path.read_bytes()
            line_count = len(raw.decode('utf-8').splitlines())
        except (OSError, UnicodeDecodeError):
            continue
        if sha(raw) != digest or end > line_count:
            continue
        record = {'id': sid, 'path': rel, 'sha256': digest,
                  'line_start': start, 'line_end': end}
        if isinstance(source.get('page'), int) and not isinstance(source.get('page'), bool):
            record['page'] = source['page']
        records[sid] = record
    return records


def proposal_quality_findings(proposal, registered_sources):
    """Make proposal issue refs resolvable from the source records in quality.json."""
    issues, sources = [], {}
    for finding in proposal.get('review', {}).get('issues', []):
        if not isinstance(finding, dict):
            continue
        cited={e.get('source_id') for e in finding.get('evidence', [])
               if isinstance(e, dict) and isinstance(e.get('source_id'), str)}
        if not cited.issubset(registered_sources):
            issues.append({'code':'second_pass_source_unregistered',
                'message':'Ein Quellenverweis des zweiten Prüfungsangebots besitzt keine aktuell registrierten Pfad-, Hash- und Zeilenmetadaten.',
                'field_path':finding.get('field_path',''),
                'source_refs':sorted(cited-registered_sources.keys()),
                'proposal_only':True,'reviewed':False})
            continue
        for sid in cited:sources[sid]=registered_sources[sid]
        issues.append({'code':'second_pass_proposal','message':finding['finding'],
            'field_path':finding['field_path'],'source_refs':sorted(cited),
            'proposal_only':True,'reviewed':False})
    return issues, sources


def annotation_source_issues(semantic, source_index, repo_root=REPO):
    """Check refs against the full hashed graph registry, not a bounded packet assembly."""
    issues=[]
    for sid in sorted(semantic_source_refs(semantic)):
        source=source_index.get(sid)
        if not source:
            issues.append({'code':'source_ref_unregistered','message':'Semantischer Quellenverweis ist nicht im aktuellen Knowledge-Graph-Quellenregister vorhanden.','source_refs':[sid]})
            continue
        rel=source.get('path'); digest=source.get('sha256')
        if not isinstance(rel,str) or not isinstance(digest,str):
            issues.append({'code':'source_ref_unregistered','message':'Quellenregistereintrag besitzt keinen Pfad oder SHA-256.','source_refs':[sid]})
            continue
        path=(repo_root/rel).resolve()
        if not path.is_relative_to(repo_root.resolve()) or not path.is_file():
            issues.append({'code':'source_ref_stale','message':'Registrierte Quelle ist nicht mehr lesbar unter ihrem Repository-Pfad.','source_refs':[sid]})
            continue
        try: current=sha(path.read_bytes())
        except OSError:
            current=''
        if current!=digest:
            issues.append({'code':'source_ref_stale','message':'Dateiinhalt stimmt nicht mehr mit dem SHA-256 des Quellenregisters überein.','source_refs':[sid]})
    return issues


def validation_source_records(packet, semantic, source_index, repo_root=REPO, file_cache=None):
    """Supply cited, hash-current graph excerpts so claim refs are checked against their evidence."""
    file_cache=file_cache if file_cache is not None else {}
    records={s['id']:s for s in packet.get('sources',[]) if isinstance(s,dict) and isinstance(s.get('id'),str)}
    for sid in semantic_source_refs(semantic):
        if sid in records:continue
        source=source_index.get(sid)
        if not source:continue
        rel=source.get('path'); digest=source.get('sha256')
        if not isinstance(rel,str) or not isinstance(digest,str):continue
        path=(repo_root/rel).resolve()
        if not path.is_relative_to(repo_root.resolve()) or not path.is_file():continue
        try:
            if rel not in file_cache:
                raw=path.read_bytes()
                if sha(raw)!=digest:continue
                file_cache[rel]=raw.decode('utf-8').splitlines()
            lines=file_cache[rel]
        except (OSError,UnicodeDecodeError):continue
        start=max(1,int(source.get('line_start',1) or 1)); end=max(start,int(source.get('line_end',start) or start))
        if end>len(lines):continue
        records[sid]={'id':sid,'path':rel,'sha256':digest,'text':'\n'.join(lines[start-1:end])}
    return list(records.values())


def validation_lean_candidates(packet, semantic, lean_index, source_index, repo_root=REPO):
    """Add only cited, typed declarations whose exact owning source is current.

    This is a local audit supplement. The model's bounded extraction packet and
    its validation/schema limits remain unchanged.
    """
    candidates={item['theorem_id']:item for item in packet.get('lean_candidates',[])}
    for link in semantic.get('lean_links',[]):
        tid=link.get('theorem_id')
        if tid in candidates:continue
        node=lean_index.get(tid)
        if not node or node.get('type') not in {'lean_theorem','lean_declaration'}:continue
        attrs=node.get('attrs',{});rel=attrs.get('path');line=attrs.get('line')
        if (not isinstance(rel,str) or not rel.endswith('.lean')
                or not isinstance(line,int) or isinstance(line,bool) or line<1):continue
        refs=attrs.get('source_refs',[])
        if not refs:continue
        verified=True
        for sid in refs:
            source=source_index.get(sid,{})
            start,end=source.get('line_start'),source.get('line_end')
            if (source.get('path')!=rel or not isinstance(source.get('sha256'),str)
                    or not isinstance(start,int) or isinstance(start,bool)
                    or not isinstance(end,int) or isinstance(end,bool)
                    or not 1<=start<=line<=end):
                verified=False;break
            path=(repo_root/rel).resolve()
            if not path.is_relative_to(repo_root.resolve()) or not path.is_file():
                verified=False;break
            try:
                raw=path.read_bytes()
                if sha(raw)!=source['sha256'] or end>len(raw.decode('utf-8').splitlines()):
                    verified=False;break
            except (OSError,UnicodeDecodeError):
                verified=False;break
        if not verified:continue
        candidates[tid]={**attrs,'theorem_id':tid,'label':node.get('label',''),
                         'proposition_surface':node.get('summary','')}
    return list(candidates.values())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--follow',action='store_true');a=p.parse_args()
    corpus=Corpus(REPO);records={};checked={};manual_path=HERE/'data/source-reviews.json';source_index={};lean_index={};catalog_stamp=None
    validation_file_cache={}
    while True:
        current_catalog_stamp=source_manifest_stamp()
        if current_catalog_stamp!=catalog_stamp:
            source_index,lean_index=load_validation_indexes();catalog_stamp=current_catalog_stamp
            checked.clear();validation_file_cache.clear()
        manual=json.loads(manual_path.read_text()) if manual_path.exists() else {}
        proposals={};proposal_sources={}
        for path in (HERE/'data/second-pass-proposals').glob('*.json'):
            proposal=json.loads(path.read_text())
            review=proposal.get('review',{})
            if proposal.get('proposal_only') and review.get('claim_id'):
                proposals[review['claim_id']]=proposal
                proposal_sources[review['claim_id']]=verified_proposal_sources(proposal)
        nodes=claims()
        # Rebuild source links on every cycle: the per-claim check cache can skip
        # unchanged findings, while quality.json is replaced wholesale each pass.
        quality_sources={}
        by_claim={n['id']:n for n in nodes}
        for claim_id,proposal in proposals.items():
            node=by_claim.get(claim_id)
            if not node:continue
            semantic=node.get('attrs',{}).get('semantic',{})
            ann_sha=sha(json.dumps(semantic,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode())
            if proposal.get('source_annotation_sha256')!=ann_sha:continue
            fresh=all((REPO/s['path']).is_file() and sha((REPO/s['path']).read_bytes())==s['sha256'] for s in proposal.get('source_hashes',[]))
            if fresh:
                _,registered=proposal_quality_findings(proposal,proposal_sources.get(claim_id,{}))
                quality_sources.update(registered)
        for n in nodes:
            sem=n['attrs']['semantic'];digest=fingerprint(sem);cid=n['id'].removeprefix('claim:')
            annotation_sha=sha(json.dumps(sem,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode())
            source_sha=source_freshness_signature(n,source_index)
            proposal=proposals.get(n['id'],{})
            check_key=annotation_sha+source_sha+sha(json.dumps([manual.get(n['id'],{}),proposal],sort_keys=True).encode())
            if checked.get(cid)==check_key:continue
            packet=corpus.packet(corpus.by_id[cid]);issues=[]
            refs=semantic_source_refs(sem)
            # Annotation IDs are validated against the full hashed graph registry;
            # the packet is an intentionally bounded, claim-specific source sample.
            validation_packet={**packet,
                'sources':validation_source_records(packet,sem,source_index,file_cache=validation_file_cache),
                'lean_candidates':validation_lean_candidates(packet,sem,lean_index,source_index)}
            try:validate(sem,validation_packet)
            except (ValueError,KeyError,TypeError) as e:issues.append({'code':'annotation_integrity','message':str(e)})
            issues.extend(annotation_source_issues(sem,source_index))
            review=manual.get(n['id'],{})
            if review.get('source_semantic_sha256')==digest:
                issues.extend(issue for issue in review.get('issues',[]) if manual_issue_still_applies(issue,sem,cid))
            if proposal.get('source_annotation_sha256')==annotation_sha:
                fresh=all((REPO/s['path']).is_file() and sha((REPO/s['path']).read_bytes())==s['sha256'] for s in proposal.get('source_hashes',[]))
                if fresh:
                    proposal_issues,registered=proposal_quality_findings(proposal,proposal_sources.get(n['id'],{}))
                    issues.extend(proposal_issues)
                    quality_sources.update(registered)
                else:issues.append({'code':'second_pass_sources_changed','message':'Quellen des zweiten KI-Prüfvorschlags haben sich geändert; erneut prüfen.'})
            quality={'state':'needs_review' if issues else 'structurally_checked','source_semantic_sha256':digest,'annotation_sha256':annotation_sha,'checked_at':now(),'issues':issues,'annotation_source_ids':sorted(refs),'packet_source_count':len(packet.get('sources',[])),'packet_source_sha256':sha(json.dumps(sorted(s['id'] for s in packet.get('sources',[])),separators=(',',':')).encode()),'scope':'References are checked against the full hashed graph source registry. Claim-specific packet coverage is recorded separately; it is not treated as annotation corruption. Structural checks do not certify general semantic correctness or physical proof.'}
            records[n['id']]={'id':n['id'],'type':'claim','label':n['label'],'summary':'','attrs':{'extraction_quality':quality}}
            checked[cid]=check_key
        atomic(HERE/'data/quality.json',{'meta':{'generated_at':now(),'kind':'annotation_quality','claims_checked':len(records),'needs_review':sum(n['attrs']['extraction_quality']['state']=='needs_review' for n in records.values()),'notice':'Structural validation and a named source sample; no full semantic certification.'},'nodes':list(records.values()),'edges':[],'sources':list(quality_sources.values())})
        print(json.dumps({'checked':len(records),'needs_review':sum(n['attrs']['extraction_quality']['state']=='needs_review' for n in records.values())}),flush=True)
        if not a.follow:break
        if not (HERE/'data/.extract.lock').exists():
            time.sleep(15)
            if not (HERE/'data/.extract.lock').exists():break
        time.sleep(20)

if __name__=='__main__':main()
