#!/usr/bin/env python3
"""Build a deterministic, source-grounded supplement for TFPT semantic graph extraction.

This is a structural inventory only. It does not infer theorem validity, proof status,
or physical interpretation from prose or keyword overlap.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable

try:
    from . import source_inventory
    from .json_transport import json_transport_value
except ImportError:  # direct `python semantic_explorer/build_corpus.py`
    import source_inventory
    from json_transport import json_transport_value

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "semantic_explorer/data/corpus.json"
MAX_ARTIFACT_BYTES = 262_144
ALLOWED_EXPERIMENT_FILES = {"README.md", "contract_index.json", "PROOF.txt", "PROOF.md", "RESULTS.md", "results.json", "result.json"}
SKIP_DIRS = {"__pycache__", ".git", "node_modules", "outputs", "logs", "_archive"}
ROOT_TEX = [
    "introduction.tex", "tfpt_1_architecture_e8.tex",
    "tfpt_2_standard_model.tex", "tfpt_3_e8_audit_bootstrap.tex",
    "tfpt_4_frontier.tex", "tfpt_5_redteam.tex",
    "tfpt_research_contracts.tex", "tfpt_horizon_readouts.tex",
    "origin_theory.tex", "tfpt_safeguards.tex",
    "tfpt_prime_front.tex", "note_e8_gaussian_code.tex",
    "note_hilbert_polya_truncations.tex",
]
STANDALONE_TEX = [
    "TFPT-Kernessenz-Deutsch-2026-10-15-v2.tex",
    "parity_toeplitz_classification.tex", "big_picture_2026-08-02_de.tex",
    "multilevel_toeplitz_lowerbound.tex", "note_w1_suzuki_identification.tex",
    "note_w3_detector_structure.tex", "TFPT_Universalraum_Sitzungsdokumentation_2026-09-20_v2.tex",
    "tfpt_universalraum_einfach_2026-09-12.tex", "TFPT_Universalraum_Sitzungsdokumentation_2026-09-20.tex",
    "tfpt_universalraum_synthesis_2026-09-12.tex", "big_picture_2026-08-01_de.tex",
    # Proof papers linked by the public summaries; the numeric twins do not
    # replace their volume and thermodynamic-limit arguments.
    "articles/2026-08-30/det16_mirror_gap_theorem_en.tex",
    "articles/2026-08-30/tfpt_thermodynamic_dynamics_en.tex",
]
INPUT_RE = re.compile(r"\\(?:input|include)\s*\{([^}]+)\}")
SECTION_RE = re.compile(r"\\(chapter|section|subsection|subsubsection|paragraph|subparagraph)\*?\s*(?:\[[^]]*\])?\s*\{([^}]*)\}")
LABEL_RE = re.compile(r"\\label\s*\{([^}]+)\}")
REF_RE = re.compile(r"\\(?:ref|eqref|autoref|cref|Cref|pageref)\s*\{([^}]+)\}")
CLAIM_RE = re.compile(r"(?<![A-Za-z0-9_])(?:TFPT\.)?[A-Z][A-Z0-9_-]*(?:\.[A-Z0-9_-]+){1,4}(?![A-Za-z0-9_])")
ENV_RE = re.compile(r"\\begin\s*\{([A-Za-z*]+)\}")
FORMULA_ENVS = {"equation", "equation*", "align", "align*", "gather", "gather*", "multline", "multline*", "flalign", "flalign*", "displaymath"}
THEOREM_ENVS = {"theorem", "lemma", "proposition", "corollary", "definition", "axiom", "conjecture", "claim", "remark", "example", "proof"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def sid(path: str, digest: str) -> str:
    return "source:" + hashlib.sha256(f"{path}|{digest}".encode()).hexdigest()[:20]


def nid(kind: str, source_id: str, tag: str) -> str:
    return f"{kind}:" + hashlib.sha256(f"{source_id}|{tag}".encode()).hexdigest()[:20]


def clean_tex_lines(text: str) -> list[str]:
    # Remove unescaped TeX comments while keeping line positions and rough contexts.
    out=[]
    for line in text.splitlines():
        cut=len(line)
        for i,ch in enumerate(line):
            if ch == "%":
                n=0; j=i-1
                while j>=0 and line[j]=="\\": n+=1; j-=1
                if n%2==0: cut=i; break
        out.append(line[:cut])
    return out


def read_source(path: Path, kind: str, lang: str) -> tuple[dict[str, Any], str] | None:
    try: data=path.read_bytes()
    except OSError: return None
    text=data.decode("utf-8", errors="replace")
    digest=sha(data); p=rel(path); source_id=sid(p,digest)
    return ({"id":source_id,"path":p,"sha256":digest,"line_count":len(text.splitlines()),"language":lang,"kind":kind,"bytes":len(data)}, text)


def detect_lang(path: str, text: str) -> str:
    name=Path(path).name.lower()
    if re.search(r"(?:_de|deutsch|sitzungsdokumentation|einfach)",name): return "de"
    de=len(re.findall(r"\b(?:und|der|die|das|nicht|ist|wird|für|eine|dieser|diese)\b",text,re.I))
    en=len(re.findall(r"\b(?:the|and|not|is|are|will|for|this|that|with)\b",text,re.I))
    return "de" if de>en*1.4 else "en" if en>de*1.4 else "mixed"


def add_node(nodes, node): nodes.append(node)
def add_edge(edges, src, dst, typ, attrs=None): edges.append({"src":src,"dst":dst,"type":typ,"attrs":attrs or {}})

def validate_unique_ids(rows, label):
    seen=set()
    for row in rows:
        node_id=row.get("id")
        if not isinstance(node_id,str) or not node_id.strip():
            raise ValueError(f"{label} contains an empty or non-string id")
        if node_id in seen:
            raise ValueError(f"{label} contains duplicate id: {node_id}")
        seen.add(node_id)


def text_context(lines: list[str], start: int, end: int, max_chars: int=1600) -> str:
    s="\n".join(lines[max(0,start-1):min(len(lines),end)])
    s=re.sub(r"\s+", " ", s).strip()
    return s[:max_chars]


def tex_files(roots: Iterable[str]) -> tuple[list[Path], list[dict[str,str]]]:
    found=[]; gaps=[]; seen=set(); stack=[ROOT/p for p in roots]
    while stack:
        p=stack.pop(0)
        if p in seen: continue
        seen.add(p)
        if not p.is_file(): gaps.append({"path":str(p),"reason":"missing_tex_root_or_include"}); continue
        found.append(p)
        try: text=p.read_text(encoding="utf-8",errors="replace")
        except OSError as e: gaps.append({"path":rel(p),"reason":"read_error:"+type(e).__name__}); continue
        for inc in INPUT_RE.findall("\n".join(clean_tex_lines(text))):
            raw_target=Path(inc.strip())
            if raw_target.suffix=="": raw_target=raw_target.with_suffix(".tex")
            # TeX resolves inputs from the compilation working directory; fall back to
            # the included file's directory for standalone source fragments.
            target=ROOT/raw_target
            if not target.is_file(): target=p.parent/raw_target
            try: target.resolve().relative_to(ROOT)
            except ValueError: gaps.append({"path":rel(p),"reason":"include_outside_repo:"+inc}); continue
            if target not in seen: stack.append(target)
    return sorted(found,key=lambda p:rel(p)),gaps


def parse_tex(path: Path, sources, nodes, edges, refs, labels, claims_valid, topic_keywords):
    got=read_source(path,"tex","mixed")
    if not got: return
    source,text=got; sources.append(source); source_id=source["id"]; p=source["path"]
    source["language"]=detect_lang(p,text)
    lines=clean_tex_lines(text)
    paper_id="paper:"+hashlib.sha256(p.encode()).hexdigest()[:20]
    add_node(nodes,{"id":paper_id,"type":"paper_source","label":Path(p).stem,"summary":"LaTeX-Quelle; strukturelle Extraktion aus Originaltext, keine inhaltliche Validierung.","attrs":{"source_id":source_id,"path":p,"language":source["language"],"root_role":"reachable_tex"}})
    add_edge(edges,paper_id,source_id,"has_source")
    # Keep section spans and raw context for later formula association.
    secs=[]
    for i,line in enumerate(lines):
        m=SECTION_RE.search(line)
        if m: secs.append((i+1,m.group(1),m.group(2).strip()))
    sec_nodes=[]
    for j,(ln,level,title) in enumerate(secs):
        end=(secs[j+1][0]-1) if j+1<len(secs) else len(lines)
        section_id=nid("section",source_id,f"{ln}:{level}:{title}")
        sec_nodes.append((ln,end,section_id,title))
        add_node(nodes,{"id":section_id,"type":"paper_section","label":title or f"(untitled {level})","summary":"LaTeX section heading; semantic status is structural only.","attrs":{"source_id":source_id,"path":p,"line_start":ln,"line_end":end,"level":level,"language":source["language"],"reviewed":False}})
        add_edge(edges,paper_id,section_id,"contains_section")
    def section_at(ln):
        matches=[s for s in sec_nodes if s[0]<=ln<=s[1]]
        return matches[-1][2] if matches else paper_id
    for i,line in enumerate(lines,1):
        for label in LABEL_RE.findall(line):
            label_id="label:"+hashlib.sha256(label.encode()).hexdigest()[:20]
            label_node_id=nid("labelref",source_id,f"{i}:{label}")
            labels.setdefault(label,[]).append(label_node_id)
            add_node(nodes,{"id":label_node_id,"type":"paper_label","label":label,"summary":"Explizites LaTeX-Label; Zielauflösung erfolgt corpusweit.","attrs":{"source_id":source_id,"path":p,"line_start":i,"line_end":i,"label_id":label_id}})
            add_edge(edges,section_at(i),label_node_id,"defines_label")
        for m in REF_RE.finditer(line):
            for target in m.group(1).split(","):
                refs.append((source_id,p,i,target.strip(),section_at(i)))
        exact_ids=set(CLAIM_RE.findall(line))
        for claim in sorted(exact_ids):
            claim=claim.removeprefix("TFPT.") if claim.startswith("TFPT.") and claim.removeprefix("TFPT.") in claims_valid else claim
            mention=(source_id,p,i,claim,section_at(i))
            if claim in claims_valid: refs.append((source_id,p,i,"CLAIM:"+claim,section_at(i)))
            else: refs.append((source_id,p,i,"CANDIDATE:"+claim,section_at(i)))
        # Preserve complete theorem/formula environments, including their source wording.
        em=ENV_RE.search(line)
        if em and em.group(1) in THEOREM_ENVS|FORMULA_ENVS:
            env=em.group(1); close=re.compile(r"\\end\s*\{\s*"+re.escape(env)+r"\s*\}")
            end=i
            while end<len(lines) and not close.search(lines[end-1]): end+=1
            ctx=text_context(lines,i,end)
            typ="formula_context" if env in FORMULA_ENVS else "theorem_environment"
            eid=nid(typ,source_id,f"{i}:{env}")
            add_node(nodes,{"id":eid,"type":typ,"label":f"{env} at {Path(p).name}:{i}","summary":"Raw source context; no theorem/proof validity inferred.","attrs":{"source_id":source_id,"path":p,"line_start":i,"line_end":end,"environment":env,"context":ctx,"language":source["language"],"epistemic":"source_text_only","reviewed":False}})
            add_edge(edges,section_at(i),eid,"contains_source_context")
            if end>=len(lines) and not close.search(lines[-1] if lines else ""):
                source.setdefault("coverage_gaps",[]).append({"line_start":i,"reason":"unterminated_environment","environment":env})
    return paper_id


def integrate_reviewed_aspects(aspect_doc: dict[str, Any], sources: list[dict[str, Any]], nodes: list[dict[str, Any]], edges: list[dict[str, Any]], valid_claims: set[str], topics: set[str]) -> None:
    """Add curated source-reviewed aspect mappings; associations never assert proof."""
    by_path={s["path"]:s for s in sources}
    sections={n["attrs"].get("path"):[] for n in nodes if n["type"]=="paper_section"}
    for n in nodes:
        if n["type"]=="paper_section": sections.setdefault(n["attrs"].get("path"),[]).append(n)
    source_nodes={n["id"] for n in nodes if n["type"] in ("tex_source","experiment_source")}
    records_by_source={}
    for n in nodes:
        sid0=n.get("attrs",{}).get("source_id")
        if sid0: records_by_source.setdefault(sid0,[]).append(n)
    used_ids=set()
    for item in aspect_doc.get("aspects",[]):
        aid=item.get("id"); kind=item.get("type")
        if not isinstance(aid,str) or not aid.startswith("aspect:") or aid in used_ids: raise ValueError(f"invalid or duplicate aspect id: {aid}")
        used_ids.add(aid)
        if kind not in ("theory_aspect","research_question"): raise ValueError(f"unsupported aspect type: {kind}")
        if any(c not in valid_claims for c in item.get("claim_ids",[])): raise ValueError(f"unknown ledger claim in {aid}")
        if any(t not in topics for t in item.get("topic_ids",[])): raise ValueError(f"unknown theory topic in {aid}")
        refs=[]
        for ref in item.get("source_refs",[]):
            p=ref.get("path"); src=by_path.get(p)
            if not src: raise ValueError(f"aspect source missing from inventory: {p}")
            lo=max(1,int(ref.get("line_start",1))); hi=max(lo,int(ref.get("line_end",lo)))
            ref2={**ref,"source_id":src["id"],"line_start":lo,"line_end":min(hi,src.get("line_count") or hi)}
            target_sections=sections.get(p,[])
            sec=next((n for n in target_sections if n["attrs"].get("line_start")==lo),None)
            if sec is None: sec=next((n for n in target_sections if n["attrs"].get("line_start",0)<=lo<=n["attrs"].get("line_end",0)),None)
            if sec: ref2["section_id"]=sec["id"]
            try:
                source_path=(ROOT/p).resolve(); source_path.relative_to(ROOT.resolve())
                raw=source_path.read_text(encoding="utf-8",errors="replace").splitlines()
                ref2["source_excerpt"]=re.sub(r"\s+"," "," ".join(raw[lo-1:min(hi,len(raw))])).strip()[:1200]
            except (OSError,ValueError): ref2["source_excerpt"]=""
            refs.append(ref2)
        attrs={k:v for k,v in item.items() if k not in ("id","type","label","summary_de","summary_en","source_refs")}
        attrs.update({"summary_de":item.get("summary_de",""),"summary_en":item.get("summary_en",""),"source_refs":refs,"extraction_method":"source_reviewed_mapping","reviewed_scope":aspect_doc.get("reviewed_scope"),"reviewed":True})
        add_node(nodes,{"id":aid,"type":kind,"label":item["label"],"summary":item.get("summary_de") or item.get("summary_en") or item["label"],"attrs":attrs})
        for ref in refs:
            add_edge(edges,aid,ref["source_id"],"grounded_in_source",{"path":ref["path"],"line_start":ref["line_start"],"line_end":ref["line_end"],"extraction_method":"source_reviewed_mapping"})
            if ref.get("section_id"): add_edge(edges,aid,ref["section_id"],"grounded_in_section",{"path":ref["path"],"line_start":ref["line_start"],"line_end":ref["line_end"]})
            for record in records_by_source.get(ref["source_id"],[]):
                if record["type"]=="experiment_record": add_edge(edges,aid,record["id"],"grounded_in_experiment_record",{"path":ref["path"]})
        for tid in item.get("topic_ids",[]): add_edge(edges,aid,tid,"about_topic",{"extraction_method":"source_reviewed_mapping","reviewed":True})
        for claim in item.get("claim_ids",[]):
            add_edge(edges,aid,"claim:"+claim,"associated_claim",{"mapping":"relevance_only","supports":False,"proves":False,"status_changed":False,"rationale":item.get("mapping_rationale","")})


def resolve_declared_script(reference: str, artifact: Path, experiment: Path) -> Path | None:
    """Resolve explicit local paths before registered verification fallbacks.

    Never merge two checkers by basename or let a source path escape the repo.
    """
    raw=Path(reference)
    variants=[raw] if raw.suffix else [raw,Path(str(raw)+".py")]
    bases=[artifact.parent,experiment,ROOT/"verification"]
    if len(raw.parts)>1 and raw.parts[0] not in {'.','..'}:bases.insert(0,ROOT)
    for base in bases:
        for variant in variants:
            candidate=(base/variant).resolve()
            if not candidate.is_relative_to(ROOT.resolve()):continue
            if candidate.is_file():return candidate
    return None


def add_declared_script(path: Path,sources: list,nodes: list) -> str:
    """A script's canonical path is its identity; indexing never records a run."""
    relative=rel(path)
    # Registered verification nodes are already supplied by the base extractor.
    registry=ROOT/"verification/script_registry.csv"
    if path.parent.resolve()==(ROOT/"verification").resolve() and registry.is_file():
        with registry.open(encoding="utf-8-sig",newline="") as stream:
            registered={Path(row.get("script","")).stem for row in csv.DictReader(stream)}
        if path.stem in registered:return "script:"+path.stem
    node_id="experiment_script:"+hashlib.sha256(relative.encode()).hexdigest()[:20]
    if not any(n["id"]==node_id for n in nodes):
        got=read_source(path,"python","en")
        if got is None:raise ValueError("Declared script is unreadable: "+relative)
        source,_text=got
        if not any(item["id"]==source["id"] for item in sources):sources.append(source)
        add_node(nodes,{"id":node_id,"type":"experiment_script","label":relative,
            "summary":"Explicitly declared exploratory checker source; no execution or proof status inferred.",
            "attrs":{"path":relative,"sha256":source["sha256"],"source_id":source["id"],
                "source_role":"exploratory","reviewed":False,"executed":False}})
    return node_id


def resolve_contract_references(sources: list,nodes: list,edges: list,gaps: list) -> None:
    """Join declared IDs to the exact owning source record, never a label guess."""
    source_paths={source["id"]:source["path"] for source in sources}
    available={node["id"] for node in nodes}
    kept=[]
    for edge in edges:
        if edge["type"]!="references_contract" or edge["dst"] in available:
            kept.append(edge);continue
        path=source_paths.get(edge.get("attrs",{}).get("source_id"),"")
        parts=Path(path).parts
        canonical="contract:"+parts[2] if len(parts)>3 and parts[:2]==("experiments","theory-contracts") else None
        if canonical in available:
            kept.append({**edge,"dst":canonical,"attrs":{**edge.get("attrs",{}),
                "declared_contract_id":edge["dst"].removeprefix("contract:"),
                "identity_method":"exact_source_owning_contract_directory","semantic":False,
                "notice":"The declared identifier is recorded by this exact contract source; no equivalence of scientific claims inferred."}})
        else:
            gaps.append({"path":path,"reason":"contract_reference_unresolved","reference":edge["dst"],
                "source_id":edge.get("attrs",{}).get("source_id"),"scope":"No exact source-owned contract node; no graph edge emitted."})
    edges[:]=kept


def main() -> dict[str,Any]:
    sources=[]; nodes=[]; edges=[]; gaps=[]; refs=[]; labels={}
    ledger_path=ROOT/"verification/status_ledger.csv"
    valid_claims=set()
    if ledger_path.exists():
        with ledger_path.open(encoding="utf-8-sig",newline="") as f: valid_claims={r["claim_id"].strip() for r in csv.DictReader(f) if r.get("claim_id")}
    orient_path=ROOT/"semantic_explorer/orientation.json"
    topics=[]
    try:
        orient=json.loads(orient_path.read_text(encoding="utf8")); topics=orient.get("topics",[])
    except Exception: topics=[]
    topic_keywords={t["id"]:[str(k).lower() for k in t.get("keywords",[]) if len(str(k))>=4] for t in topics}

    roots=[p for p in ROOT_TEX+STANDALONE_TEX if (ROOT/p).is_file()]
    # docset and its explicit root list are themselves provenance-bearing input.
    if (ROOT/"tex-artefacts/tfpt_docset.tex").is_file(): roots.append("tex-artefacts/tfpt_docset.tex")
    tex,g=tex_files(roots); gaps.extend(g)
    for p in tex: parse_tex(p,sources,nodes,edges,refs,labels,valid_claims,topic_keywords)
    for src_id,p,ln,target,owner in refs:
        if target.startswith("CLAIM:"):
            c=target[6:]; add_edge(edges,owner,"claim:"+c,"explicit_claim_mention",{"source_id":src_id,"path":p,"line_start":ln,"line_end":ln,"semantic":False,"reviewed":False})
        elif target.startswith("CANDIDATE:"):
            c=target[10:]; rid=nid("claim_candidate",src_id,f"{ln}:{c}")
            add_node(nodes,{"id":rid,"type":"claim_reference_candidate","label":c,"summary":"Unverified claim-like token; unresolved and not linked to a claim node.","attrs":{"source_id":src_id,"path":p,"line_start":ln,"line_end":ln,"reviewed":False}})
            add_edge(edges,owner,rid,"candidate_reference",{"semantic":False,"reviewed":False})
        elif target in labels:
            add_edge(edges,owner,labels[target][0],"references_label",{"label":target,"source_id":src_id,"path":p,"line_start":ln,"line_end":ln})
        else:
            rid=nid("label_candidate",src_id,f"{ln}:{target}")
            add_node(nodes,{"id":rid,"type":"unresolved_label_reference","label":target,"summary":"Unresolved explicit LaTeX reference.","attrs":{"source_id":src_id,"path":p,"line_start":ln,"line_end":ln}})
            add_edge(edges,owner,rid,"unresolved_reference")

    # Explicit topic tokens from paper section/title/context are keyword candidates only.
    topics_by_id={t["id"]:t for t in topics}
    topic_id_re=re.compile(r"\btopic:[a-z0-9][a-z0-9_-]*\b")
    topic_ids=set(topics_by_id)
    for n in nodes:
        if n["type"] not in ("paper_source","paper_section","formula_context","theorem_environment","experiment_record"): continue
        hay=(n["label"]+" "+n["summary"]+" "+str(n["attrs"].get("context",n["attrs"].get("excerpt","")))).lower()
        raw=n["label"]+" "+n["summary"]+" "+str(n["attrs"].get("context",n["attrs"].get("excerpt","")))
        for tid in sorted(set(topic_id_re.findall(raw)) & topic_ids):
            add_edge(edges,n["id"],tid,"explicit_topic_reference",{"semantic":False,"reviewed":False,"method":"exact_topic_id"})
        for tid,words in topic_keywords.items():
            matched=sorted({w for w in words if w in hay})
            if matched:
                add_edge(edges,n["id"],tid,"candidate_topic_match",{"matched_terms":matched[:12],"method":"keyword_candidate","semantic":False,"reviewed":False})

    # Experiment inventory: whitelist named evidence artifacts only; no scripts/logs executed or ingested.
    exp_root=ROOT/"experiments"; exp_dirs=[]
    if exp_root.exists():
        for p in exp_root.iterdir():
            if p.is_dir() and p.name not in SKIP_DIRS: exp_dirs.append(p)
    contract_root=exp_root/"theory-contracts"
    contract_dirs=sorted([p for p in contract_root.iterdir() if p.is_dir() and p.name not in SKIP_DIRS],key=lambda p:p.name) if contract_root.exists() else []
    indexed_contracts=set()
    experiment_node_ids={}
    for d in sorted(exp_dirs,key=lambda p:p.name):
        eid="experiment:"+hashlib.sha256(rel(d).encode()).hexdigest()[:20]
        experiment_node_ids[d.resolve()]=eid
        add_node(nodes,{"id":eid,"type":"experiment","label":d.name,"summary":"Experiment directory inventory; artifacts are source records, not independent evidence.","attrs":{"path":rel(d),"reviewed":False}})
    contract_node_ids={}
    if contract_root.exists():
        parent_id=experiment_node_ids.get(contract_root.resolve())
        for cd in contract_dirs:
            cid="contract_project:"+hashlib.sha256(rel(cd).encode()).hexdigest()[:20]
            contract_node_ids[cd.resolve()]=cid
            add_node(nodes,{"id":cid,"type":"experiment","label":cd.name,"summary":"Individual theory-contract directory inventory, including directories without a contract_index.json.","attrs":{"path":rel(cd),"parent":"experiments/theory-contracts","reviewed":False,"contract_index_present":(cd/"contract_index.json").is_file()}})
            if parent_id: add_edge(edges,parent_id,cid,"contains_experiment_directory",{"path":rel(cd)})
    for d in sorted(exp_dirs,key=lambda p:p.name):
        eid=experiment_node_ids[d.resolve()]
        # Exact path whitelist; prune ignored directories and never read arbitrary logs or keys.
        files=[]
        for base,dirs,names in __import__("os").walk(d):
            dirs[:]=[x for x in dirs if x not in SKIP_DIRS and not x.startswith(".")]
            for name in names:
                if name in ALLOWED_EXPERIMENT_FILES: files.append(Path(base)/name)
            contract_dir=Path(base).resolve()
            if contract_dir in contract_node_ids and "contract_index.json" in names:
                indexed=source_inventory._read_json(contract_dir/"contract_index.json")
                files.extend(path for _field,path in source_inventory._declared_contract_artifacts(ROOT,contract_dir,indexed))
                manifest_path,manifest_entries=source_inventory._declared_source_manifest(ROOT,contract_dir,indexed)
                if manifest_path: files.append(manifest_path)
                files.extend(ROOT/entry["path"] for entry in manifest_entries)
        # Declared proofs/checkers/manifests may already have entered through
        # the regular allowlist walk. Process a canonical path once so the
        # same artifact cannot create duplicate source-backed graph nodes.
        files=list({fp.resolve():fp.resolve() for fp in files}.values())
        for fp in sorted(files,key=lambda p:rel(p)):
            try: size=fp.stat().st_size
            except OSError: continue
            p=rel(fp)
            if size>MAX_ARTIFACT_BYTES:
                raw=fp.read_bytes(); source_id=sid(p,sha(raw))
                # A manifest may list a large TeX file already inventoried as
                # a paper source. Do not add a second, weaker oversize record
                # with the same content-addressed ID.
                if any(existing.get("id")==source_id for existing in sources):
                    continue
                s={"id":source_id,"path":p,"sha256":sha(raw),"line_count":None,"language":"unknown","kind":"oversize_experiment_artifact","bytes":size}; sources.append(s)
                gaps.append({"path":p,"reason":"recognized_artifact_over_size_limit","bytes":size,"limit":MAX_ARTIFACT_BYTES}); continue
            kind=fp.suffix.lstrip(".").lower() or "text"; got=read_source(fp,kind,"mixed")
            if not got: continue
            source,text=got
            # A declared experiment dependency may also be an already indexed
            # paper. Keep both uses, but retain one provenance record and its
            # richer paper metadata, just as in the oversize branch above.
            existing_source=next((s for s in sources if s["id"]==source["id"]),None)
            if existing_source is None:
                source["language"]=detect_lang(p,text)
                sources.append(source)
            else:
                source=existing_source
            owner_id=eid
            if contract_root.resolve() in fp.resolve().parents:
                try: child=next(x for x in contract_dirs if x.resolve() in fp.resolve().parents)
                except StopIteration: child=None
                if child is not None: owner_id=contract_node_ids[child.resolve()]
            rid="experiment_record:"+hashlib.sha256(f"{owner_id}|{p}".encode()).hexdigest()[:20]
            label=fp.name
            summary="Bounded experiment source record; status remains as stated in source and is not promoted by indexing."
            attrs={"source_id":source["id"],"path":p,"line_start":1,"line_end":source["line_count"],"language":source["language"],"record_kind":fp.name,"reviewed":False}
            if fp.name in ("contract_index.json","results.json","result.json"):
                try:
                    payload=json.loads(text)
                    contract=payload.get("contract") if isinstance(payload,dict) else None
                    if isinstance(contract,str): attrs["contract_id"]=contract; label=contract
                    for k in ("question","verdict","verdict_enum","claims_tested","scripts","checker","scope_summary_pass","no_t1_t8_closure"):
                        if isinstance(payload,dict) and k in payload: attrs[k]=payload[k]
                    if fp.name=="contract_index.json" and contract: indexed_contracts.add(contract)
                except (ValueError,TypeError): attrs["parse_status"]="invalid_json"
            else:
                # Markdown/text outline only: retain headings and bounded opening excerpt.
                headings=[line.lstrip("# ").strip() for line in text.splitlines() if line.startswith("#")][:30]
                attrs["headings"]=headings; attrs["excerpt"]=" ".join(text.split())[:1200]
            add_node(nodes,{"id":rid,"type":"experiment_record","label":label,"summary":summary,"attrs":attrs})
            add_edge(edges,owner_id,rid,"has_record",{"source_id":source["id"]})
            if fp.name=="contract_index.json" and isinstance(attrs.get("question"),str):
                q=attrs["question"]
                qid="contract_question:"+hashlib.sha256(f"{source['id']}|{q}".encode()).hexdigest()[:20]
                qline=next((i for i,line in enumerate(text.splitlines(),1) if '"question"' in line),1)
                add_node(nodes,{"id":qid,"type":"contract_question","label":q[:240],"summary":"Question field copied from contract_index; no answer or proof status inferred.","attrs":{"source_id":source["id"],"path":p,"line_start":qline,"line_end":qline,"contract_id":attrs.get("contract_id"),"verdict_as_stated":attrs.get("verdict"),"reviewed":False}})
                add_edge(edges,rid,qid,"records_question",{"source_id":source["id"]})
            if isinstance(attrs.get("claims_tested"),list):
                for c in attrs["claims_tested"]:
                    if isinstance(c,str) and c in valid_claims: add_edge(edges,rid,"claim:"+c,"explicit_claim_reference",{"source_id":source["id"],"semantic":False,"reviewed":False})
                    elif isinstance(c,str): gaps.append({"path":p,"reason":"claim_id_not_in_ledger","claim_id":c})
            # Exact topic IDs are safe references; term overlap remains only a candidate edge.
            raw_record=json.dumps(attrs,ensure_ascii=False)
            known_topic_ids={t["id"] for t in topics}
            for tid in sorted(set(re.findall(r"\btopic:[a-z0-9][a-z0-9_-]*\b",raw_record)) & known_topic_ids):
                add_edge(edges,rid,tid,"explicit_topic_reference",{"source_id":source["id"],"semantic":False,"reviewed":False,"method":"exact_topic_id"})
            checker=attrs.get("checker")
            script_values=[]
            if isinstance(checker,str): script_values.append(checker)
            if isinstance(attrs.get("scripts"),list): script_values.extend(x for x in attrs["scripts"] if isinstance(x,str))
            for script_ref in script_values:
                match=resolve_declared_script(script_ref,fp,d)
                if match:
                    script_id=add_declared_script(match,sources,nodes)
                    add_edge(edges,rid,script_id,"references_testscript",{"source_id":source["id"],"path":rel(match),"exact_file":True,"semantic":False})
                else: gaps.append({"path":p,"reason":"testscript_reference_unresolved","reference":script_ref})
            cid=attrs.get("contract_id")
            if isinstance(cid,str):
                add_edge(edges,eid,"contract:"+cid,"references_contract",{"source_id":source["id"],"exact_id":True})
    expected_contracts={d.name for d in contract_dirs}
    gaps.extend({"path":rel(contract_root/name),"reason":"missing_contract_index"} for name in sorted(expected_contracts-indexed_contracts))

    # Keep source-backed historical records in the Atlas supplement without
    # reading the generated legacy theory_graph.json. The importer returns an
    # in-memory inventory only; every edge comes from an explicit source ID.
    inventory=source_inventory.build(ROOT,valid_claims)
    source_ids={s["id"] for s in sources}
    for source in inventory["sources"]:
        if source["id"] not in source_ids:
            sources.append(source)
            source_ids.add(source["id"])
    imported_types={"contract","rh","experiment","scorecard","failure"}
    for record in inventory["nodes"]:
        if record["type"] in {"research_plan","research_plan_document"}:
            add_node(nodes,record)
            continue
        if record["type"] not in imported_types:
            continue
        attrs=dict(record.get("attrs",{}))
        attrs.setdefault("attempt_kind",record["type"])
        attrs.setdefault("source_type",record["type"])
        attrs.setdefault("notice","Historical source record; not an independent proof of linked claims.")
        add_node(nodes,{**record,"type":"attempt","attrs":attrs})
    for edge_record in inventory["edges"]:
        add_edge(edges,edge_record["src"],edge_record["dst"],edge_record["type"],edge_record.get("attrs",{}))

    resolve_contract_references(sources,nodes,edges,gaps)

    # Give every provenance source an explicit graph node so local edges are closed.
    for source in sources:
        stype="tex_source" if source["kind"]=="tex" else "experiment_source"
        add_node(nodes,{"id":source["id"],"type":stype,"label":Path(source["path"]).name,"summary":"Source artifact with content hash and line-count provenance.","attrs":{k:source[k] for k in ("path","sha256","line_count","language","kind","bytes") if k in source}})
    for n in nodes:
        source_id=n.get("attrs",{}).get("source_id")
        if source_id in {s["id"] for s in sources} and n["id"]!=source_id:
            add_edge(edges,n["id"],source_id,"has_source_provenance",{"source_id":source_id})
    aspects_path=ROOT/"semantic_explorer/corpus_aspects.json"
    aspect_doc=json.loads(aspects_path.read_text(encoding="utf-8")) if aspects_path.is_file() else {"aspects":[]}
    integrate_reviewed_aspects(aspect_doc,sources,nodes,edges,valid_claims,{t["id"] for t in topics})

    # Describe coverage without changing the existing knowledge graph.
    corpus={"version":1,"purpose":"Deterministic structural source supplement plus curated source-reviewed mappings; not a proof/status authority.","aspect_review_method":"source_reviewed_mapping; claim edges mean relevance only and never support/prove or change ledger status.",
      "scope":{"latex_roots":ROOT_TEX,"standalone_tex":STANDALONE_TEX,"experiment_root":"experiments/","experiment_artifact_whitelist":sorted(ALLOWED_EXPERIMENT_FILES),"max_artifact_bytes":MAX_ARTIFACT_BYTES,"excluded":{"directories":sorted(SKIP_DIRS),"files":"logs, scripts, hidden files, keys, and unlisted artifacts"},"languages":sorted({s["language"] for s in sources}),"language_detection":"filename and small stopword heuristic; informational only"},
      "coverage":{"tex_sources":sum(1 for s in sources if s["kind"]=="tex"),"reviewed_aspects":len(aspect_doc.get("aspects",[])),"experiment_directories":len(exp_dirs),"contract_directories":len(contract_dirs),"contract_indexes":len(indexed_contracts),"experiment_records":sum(n["type"]=="experiment_record" for n in nodes),"source_inventory_records":sum(n["type"]=="attempt" for n in nodes),"source_inventory_edges":len(inventory["edges"]),"nodes":len(nodes),"edges":len(edges),"gaps":gaps},
      "sources":sorted(sources,key=lambda s:s["path"]),"nodes":sorted(nodes,key=lambda n:n["id"]),"edges":sorted(edges,key=lambda e:(e["src"],e["type"],e["dst"],json.dumps(e["attrs"],sort_keys=True)))}
    # Validate before replacing the existing corpus so a malformed import
    # cannot take down an otherwise usable Atlas snapshot.
    validate_unique_ids(corpus["sources"],"source inventory")
    validate_unique_ids(corpus["nodes"],"corpus node inventory")
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(json_transport_value(corpus),ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    return corpus

if __name__=="__main__":
    c=main(); summary={k:v for k,v in c["coverage"].items() if k!="gaps"}; summary["gap_count"]=len(c["coverage"]["gaps"]); print(json.dumps(summary,ensure_ascii=False,sort_keys=True))
