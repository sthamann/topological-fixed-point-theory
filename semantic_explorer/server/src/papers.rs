//! Bounded paper navigation inside the one Atlas graph. No scientific verdicts.
use crate::service::{KnowledgeGraph, Node};
use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs,
    path::{Component, Path},
};

pub const TOOLS: &[&str] = &[
    "list_papers",
    "paper_outline",
    "paper_section",
    "search_papers",
    "paper_evidence",
    "paper_coverage",
];
const NOTICE: &str = "Paper navigation and generated AI annotations are unreviewed; citations and matching source hashes do not establish a proof or a fresh kernel/checker run.";

#[derive(Default, Deserialize, Serialize)]
pub struct PaperArgs {
    pub query: Option<String>,
    pub scope: Option<String>,
    pub document_id: Option<String>,
    pub section_id: Option<String>,
    pub id: Option<String>,
    pub evidence_kind: Option<String>,
    pub evidence_type: Option<String>,
    pub offset: Option<usize>,
    pub limit: Option<usize>,
    pub max_chars: Option<usize>,
}

pub fn tool_specs() -> Vec<Value> {
    TOOLS.iter().map(|name| {
        let required = match *name { "paper_outline"=>vec!["document_id"], "paper_section"=>vec!["section_id"], "search_papers"=>vec!["query"], "paper_evidence"=>vec!["id"], _=>vec![] };
        json!({"name":name,"description":format!("{}: bounded paper navigation in the sole Atlas. Stable registry IDs, source freshness and unreviewed AI state remain explicit. Lists/search default to public; scope=local_private or all opts into local private content. max_chars bounds the serialized result JSON; MCP framing is additional. Use read_source for original passages. {}",name,NOTICE),"inputSchema":{"type":"object","properties":{"query":{"type":"string","maxLength":2000},"scope":{"type":"string","enum":["public","local_private","private","all","active","archive","exploration","historical","support"],"default":"public"},"document_id":{"type":"string"},"section_id":{"type":"string"},"id":{"type":"string"},"evidence_kind":{"type":"string"},"evidence_type":{"type":"string"},"offset":{"type":"integer","minimum":0,"default":0},"limit":{"type":"integer","minimum":1,"maximum":100,"default":if *name=="paper_outline"{30}else{20}},"max_chars":{"type":"integer","minimum":2000,"maximum":60000,"default":16000}},"required":required}})
    }).collect()
}

pub fn is_paper(n: &Node) -> bool {
    matches!(n.kind.as_str(), "paper_document" | "paper_section")
        && n.attrs.get("paper_layer") == Some(&Value::Bool(true))
}
pub fn source_current(n: &Node) -> bool {
    n.attrs.get("source_current") != Some(&Value::Bool(false))
}

pub fn summary_current(n: &Node, source_current: bool) -> bool {
    source_current
        && n.attrs.get("context_current") != Some(&Value::Bool(false))
        && n.attrs.get("summary_state").and_then(Value::as_str) == Some("complete")
}

pub fn searchable(n: &Node, source_current: bool) -> String {
    if !is_paper(n) {
        return String::new();
    }
    let mut text = format!("{} {}", n.id, n.label);
    if source_current {
        if let Some(body) = n.attrs.get("source_text").and_then(Value::as_str) {
            text.push(' ');
            text.push_str(body);
        }
    }
    if summary_current(n, source_current) {
        text.push(' ');
        text.push_str(&n.summary);
        if let Some(semantic) = n.attrs.get("semantic") {
            text.push(' ');
            text.push_str(&semantic.to_string());
        }
    }
    text
}

/// Paper placeholders may link to scientific nodes, but never replace them.
pub fn merge(g: &mut KnowledgeGraph, path: &Path) -> Result<()> {
    if !path.exists() {
        g.meta["papers"] = json!({"state":"not_imported","notice":"Run the paper importer; missing paper coverage is not a scientific verdict."});
        return Ok(());
    }
    if fs::metadata(path)?.len() > 256 * 1024 * 1024 {
        bail!("paper overlay exceeds 256 MiB safe read limit");
    }
    let overlay: KnowledgeGraph =
        serde_json::from_slice(&fs::read(path)?).context("invalid paper overlay JSON")?;
    let mut known: HashSet<String> = g.nodes.iter().map(|n| n.id.clone()).collect();
    let mut seen = HashSet::new();
    let mut preserved = 0;
    for mut n in overlay.nodes {
        n.attrs["reviewed"] = json!(false);
        if n.id.is_empty() || !seen.insert(n.id.clone()) {
            bail!("duplicate or empty paper overlay node ID");
        }
        if !matches!(
            n.kind.as_str(),
            "paper_document" | "paper_section" | "paper_evidence"
        ) {
            bail!("invalid paper overlay node type: {}", n.kind);
        }
        if (n.kind == "paper_document" && !n.id.starts_with("document:"))
            || (n.kind == "paper_section" && !n.id.starts_with("section:"))
        {
            bail!("invalid stable paper ID: {}", n.id);
        }
        if matches!(n.kind.as_str(), "paper_document" | "paper_section") {
            n.attrs["paper_layer"] = json!(true);
        }
        if is_paper(&n) && n.attrs.get("embedding").is_some() {
            let text: String = format!("{}\n{}", n.label, n.summary)
                .chars()
                .take(12000)
                .collect();
            let fingerprint = format!("{:x}", Sha256::digest(text.as_bytes()));
            if n.attrs.get("embedding_text_sha256").and_then(Value::as_str)
                != Some(fingerprint.as_str())
                || n.attrs.get("embedding_model").and_then(Value::as_str)
                    != Some("text-embedding-3-small")
                || n.attrs.get("embedding_dimensions").and_then(Value::as_u64) != Some(256)
                || crate::service::embedding(&n).is_none()
            {
                n.attrs.as_object_mut().unwrap().remove("embedding");
                n.attrs["embedding_state"] = json!("invalid_or_stale_fingerprint");
            }
        }
        if known.insert(n.id.clone()) {
            g.nodes.push(n);
        } else {
            preserved += 1;
        }
    }
    let mut edges: HashSet<_> = g
        .edges
        .iter()
        .map(|e| (e.src.clone(), e.dst.clone(), e.kind.clone()))
        .collect();
    for mut e in overlay.edges {
        e.attrs["reviewed"] = json!(false);
        if !known.contains(&e.src) || !known.contains(&e.dst) {
            bail!("paper edge has an unknown endpoint");
        }
        if edges.insert((e.src.clone(), e.dst.clone(), e.kind.clone())) {
            g.edges.push(e);
        }
    }
    let mut sources: HashMap<String, _> = g
        .sources
        .iter()
        .map(|s| (s.id.clone(), s.clone()))
        .collect();
    for s in overlay.sources {
        if let Some(old) = sources.get(&s.id) {
            if old.path != s.path
                || old.sha256 != s.sha256
                || old.line_start != s.line_start
                || old.line_end != s.line_end
            {
                bail!("paper source ID collision: {}", s.id);
            }
        } else {
            sources.insert(s.id.clone(), s.clone());
            g.sources.push(s);
        }
    }
    let errors = overlay.meta.get("errors").and_then(Value::as_array);
    let import = json!({"coverage":overlay.meta.get("coverage"),"model":overlay.meta.get("model"),"prompt_version":overlay.meta.get("prompt_version"),"include_local":overlay.meta.get("include_local"),"registry_fingerprints":overlay.meta.get("registry_fingerprints"),"errors_count":errors.map_or(0,Vec::len),"error_examples":errors.into_iter().flatten().take(3).map(|v|clipped(&v.to_string(),500)).collect::<Vec<_>>()});
    g.meta["papers"] = json!({"state":"loaded","import":import,"existing_scientific_nodes_preserved":preserved,"reviewed":false,"notice":NOTICE});
    Ok(())
}

#[derive(Default)]
struct FileData {
    bytes: Vec<u8>,
    sha256: String,
}
#[derive(Default)]
struct Freshness {
    files: HashMap<String, std::result::Result<FileData, String>>,
    nodes: HashMap<String, usize>,
    sources: HashMap<String, usize>,
}
impl Freshness {
    fn new(g: &KnowledgeGraph) -> Self {
        Self {
            files: HashMap::new(),
            nodes: g
                .nodes
                .iter()
                .enumerate()
                .map(|(i, n)| (n.id.clone(), i))
                .collect(),
            sources: g
                .sources
                .iter()
                .enumerate()
                .map(|(i, s)| (s.id.clone(), i))
                .collect(),
        }
    }
    fn file(&mut self, root: &Path, path: &str) -> &std::result::Result<FileData, String> {
        self.files.entry(path.to_owned()).or_insert_with(|| {
            let raw = Path::new(path);
            if raw.components().any(|c| matches!(c, Component::ParentDir)) {
                return Err("unsafe source path".into());
            }
            let full = if raw.is_absolute() {
                raw.to_owned()
            } else {
                root.join(raw)
            };
            let full = fs::canonicalize(full).map_err(|_| "source missing".to_owned())?;
            if !full.starts_with(root) || crate::service::sensitive_path(&full) {
                return Err("unsafe source path".into());
            }
            if fs::metadata(&full).map_err(|e| e.to_string())?.len() > 8 * 1024 * 1024 {
                return Err("source exceeds safe read limit".into());
            }
            fs::read(full)
                .map(|bytes| {
                    let sha256 = format!("{:x}", Sha256::digest(&bytes));
                    FileData { bytes, sha256 }
                })
                .map_err(|e| e.to_string())
        })
    }
    fn check(&mut self, root: &Path, g: &KnowledgeGraph, n: &Node) -> Value {
        if n.kind == "paper_document" {
            if n.attrs.get("source_current") == Some(&Value::Bool(false)) {
                return json!({"current":false,"state":"import_source_error","error":n.attrs.get("source_error")});
            }
            for id in n
                .attrs
                .get("source_refs")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
            {
                let Some(source) = id
                    .as_str()
                    .and_then(|id| self.sources.get(id))
                    .and_then(|i| g.sources.get(*i))
                else {
                    return json!({"current":false,"state":"unregistered_document_source"});
                };
                let bytes = match self.file(root, &source.path) {
                    Ok(b) => b,
                    Err(e) => return json!({"current":false,"state":e}),
                };
                if source.sha256.is_empty() || bytes.sha256 != source.sha256 {
                    return json!({"current":false,"state":"document_source_hash_changed"});
                }
            }
            let ids = outline_ids(g, n);
            if ids.is_empty() {
                let registered = n
                    .attrs
                    .get("source_refs")
                    .and_then(Value::as_array)
                    .is_some_and(|r| !r.is_empty());
                return json!({"current":registered,"state":if registered{"document_sources_hash_verified"}else{"no_section_sources"}});
            }
            for id in &ids {
                let Some(section) = self
                    .nodes
                    .get(id)
                    .and_then(|i| g.nodes.get(*i))
                    .filter(|s| s.kind == "paper_section")
                else {
                    return json!({"current":false,"state":"section_missing"});
                };
                let state = self.check(root, g, section);
                if state["current"] != true {
                    return state;
                }
            }
            return json!({"current":true,"state":"section_sources_hash_verified","sections_checked":ids.len()});
        }
        let refs = n
            .attrs
            .get("source_refs")
            .and_then(Value::as_array)
            .cloned()
            .unwrap_or_default();
        let mut checked = 0;
        for id in &refs {
            let Some(s) = id
                .as_str()
                .and_then(|id| self.sources.get(id))
                .and_then(|i| g.sources.get(*i))
            else {
                return json!({"current":false,"state":"unregistered_source"});
            };
            let bytes = match self.file(root, &s.path) {
                Ok(b) => b,
                Err(e) => return json!({"current":false,"state":e}),
            };
            if s.sha256.is_empty() || bytes.sha256 != s.sha256 {
                return json!({"current":false,"state":"source_hash_changed"});
            }
            checked += 1;
        }
        if let Some(source) = n.attrs.get("source").filter(|s| s.is_object()) {
            let path = source.get("path").and_then(Value::as_str).unwrap_or("");
            let bytes = match self.file(root, path) {
                Ok(b) => b,
                Err(e) => return json!({"current":false,"state":e}),
            };
            let expected = source.get("sha256").and_then(Value::as_str).unwrap_or("");
            if expected.is_empty() || bytes.sha256 != expected {
                return json!({"current":false,"state":"source_hash_changed"});
            }
            if let (Some(start), Some(end), Some(pin)) = (
                source.get("byte_start").and_then(Value::as_u64),
                source.get("byte_end").and_then(Value::as_u64),
                source.get("span_sha256").and_then(Value::as_str),
            ) {
                let Some(span) = bytes.bytes.get(start as usize..end as usize) else {
                    return json!({"current":false,"state":"invalid_source_span"});
                };
                if format!("{:x}", Sha256::digest(span)) != pin
                    || n.attrs
                        .get("source_text")
                        .and_then(Value::as_str)
                        .is_some_and(|t| t.as_bytes() != span)
                {
                    return json!({"current":false,"state":"section_span_changed"});
                }
            }
            checked += 1;
        }
        if checked == 0 {
            return json!({"current":false,"state":"source_not_registered"});
        }
        json!({"current":true,"state":"hash_verified","references_checked":checked})
    }
}

pub fn refresh_sources(g: &mut KnowledgeGraph, root: &Path) {
    let mut cache = Freshness::new(g);
    let states: Vec<_> = g
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| is_paper(n))
        .map(|(i, n)| (i, cache.check(root, g, n)))
        .collect();
    let document_states: HashMap<String, Value> = states
        .iter()
        .filter(|(i, _)| g.nodes[*i].kind == "paper_document")
        .map(|(i, state)| (g.nodes[*i].id.clone(), state.clone()))
        .collect();
    for (i, state) in states {
        let attrs = &mut g.nodes[i].attrs;
        let failures: Vec<Value> = attrs
            .get("context_dependencies")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|d| {
                let id = d
                    .get("id")
                    .or(d.get("document_id"))
                    .and_then(Value::as_str)
                    .unwrap_or("");
                let state = document_states
                    .get(id)
                    .cloned()
                    .unwrap_or_else(|| json!({"current":false,"state":"document_missing"}));
                (state["current"] != true).then(|| json!({"document_id":id,"freshness":state}))
            })
            .collect();
        let context_state = if failures.is_empty() {
            "current"
        } else if failures.iter().any(|v| changed_state(&v["freshness"])) {
            "stale"
        } else {
            "incomplete"
        };
        let recorded = attrs
            .get("summary_state")
            .cloned()
            .unwrap_or(json!("pending"));
        attrs["summary_state_as_recorded"] = recorded.clone();
        attrs["source_freshness"] = state.clone();
        attrs["source_current"] = state["current"].clone();
        attrs["context_current"] = json!(context_state == "current");
        attrs["context_complete"] = json!(context_state != "incomplete");
        attrs["context_state"] = json!(context_state);
        if state["current"] != true {
            attrs["source_error"] = state["state"].clone();
        }
        if !failures.is_empty() {
            attrs["context_error"] = json!(failures.iter().take(3).collect::<Vec<_>>());
            attrs["context_errors_omitted"] = json!(failures.len().saturating_sub(3));
        }
        if recorded == "complete" {
            if changed_state(&state) || context_state == "stale" {
                attrs["summary_state"] = json!("stale");
            } else if state["current"] != true {
                attrs["summary_state"] = json!("source_error");
            }
        }
    }
}
fn changed_state(state: &Value) -> bool {
    state["state"]
        .as_str()
        .is_some_and(|s| s.contains("changed"))
}

fn observed_freshness(n: &Node, g: &KnowledgeGraph, root: &Path, fresh: &mut Freshness) -> Value {
    n.attrs
        .get("source_freshness")
        .cloned()
        .unwrap_or_else(|| fresh.check(root, g, n))
}
fn query_tokens(query: &str) -> Vec<String> {
    const STOP: &[&str] = &[
        "der", "die", "das", "den", "dem", "des", "ein", "eine", "einer", "eines", "einem",
        "einen", "und", "oder", "ist", "sind", "war", "wird", "werden", "wie", "was", "welche",
        "welcher", "welches", "wo", "mit", "von", "für", "auf", "aus", "zu", "zur", "zum", "im",
        "in", "an", "am", "als", "auch", "es", "sich", "man", "nicht", "the", "and", "or", "is",
        "are", "was", "were", "how", "what", "which", "where", "with", "from", "for", "of", "to",
        "on", "in", "a", "an", "does", "do", "it", "this", "that",
    ];
    let mut seen = HashSet::new();
    query
        .to_lowercase()
        .split(|c: char| !c.is_alphanumeric() && c != '_' && c != '.' && c != ':')
        .filter(|s| s.chars().count() > 1 && !STOP.contains(s))
        .filter(|s| seen.insert((*s).to_owned()))
        .map(str::to_owned)
        .collect()
}

// Paper bodies are often English even when a German query or an unavailable
// annotation sends the search through the lexical path.
fn lexical_aliases(term: &str) -> &'static [&'static str] {
    match term {
        "cliffordquelle" => &["clifford"],
        "transferkern" => &["transfer kernel"],
        "hintergrund" => &["background"],
        "gekoppelt" | "gekoppelte" => &["coupled"],
        "ausschluss" => &["exclusion", "excluded"],
        "residuenmatrix" => &["residue matrix"],
        "quarkmassen" => &["quark masses"],
        "skala" => &["scale"],
        "masselücke" => &["mass gap"],
        "brückenmaß" | "brückenmaßes" => &["bridge measure"],
        "clusterbildung" => &["clustering"],
        "hamiltonkasten" => &["finite box", "hamiltonian"],
        "vollständigkeit" => &["completeness"],
        "thermodynamischer" | "thermodynamischen" | "thermodynamische" => &["thermodynamic"],
        "kontinuumslimes" | "kontinuumsgrenze" => &["continuum limit", "continuum"],
        "infrarot" => &["infrared", "ir universality"],
        _ => &[],
    }
}

fn clipped(s: &str, cap: usize) -> String {
    s.chars().take(cap).collect()
}
fn row(n: &Node, g: &KnowledgeGraph, root: &Path, fresh: &mut Freshness) -> Value {
    let freshness = observed_freshness(n, g, root, fresh);
    let current = freshness["current"] == true;
    let recorded = n
        .attrs
        .get("summary_state")
        .and_then(Value::as_str)
        .unwrap_or("pending");
    let summary = summary_current(n, current).then(|| clipped(&n.summary, 650));
    json!({"id":n.id,"type":n.kind,"label":clipped(&n.label,240),"summary":summary,"field_omissions":{"label_chars":n.label.chars().count().saturating_sub(240),"summary_chars":summary.as_ref().map_or(0, |_| n.summary.chars().count().saturating_sub(650))},
        "scope":n.attrs.get("scope"),"publication_scope":n.attrs.get("publication_scope").or(n.attrs.get("scope")),
        "content_scope":n.attrs.get("content_scope"),"source_role":n.attrs.get("source_role"),"document_role":n.attrs.get("document_role"),
        "source_refs":n.attrs.get("source_refs").and_then(Value::as_array).map(|r|r.iter().take(8).cloned().collect::<Vec<_>>()),"source_refs_omitted":n.attrs.get("source_refs").and_then(Value::as_array).map_or(0,|r|r.len().saturating_sub(8)),"source":n.attrs.get("source"),"source_freshness":freshness,"source_current":current,"source_error":if current{Value::Null}else{freshness["state"].clone()},
        "summary_state":if !current && recorded=="complete"{"stale"}else{recorded},"summary_state_as_recorded":n.attrs.get("summary_state_as_recorded").cloned().unwrap_or(json!(recorded)),
        "summary_current":summary_current(n, current),"context_current":n.attrs.get("context_current").cloned().unwrap_or(json!(true)),"context_error":n.attrs.get("context_error"),"context_state":n.attrs.get("context_state"),"context_complete":n.attrs.get("context_complete"),"reviewed":false})
}
fn matches_scope(n: &Node, scope: &str) -> bool {
    let publication = n
        .attrs
        .get("publication_scope")
        .or(n.attrs.get("scope"))
        .and_then(Value::as_str)
        .unwrap_or("public");
    match scope {
        "all" => true,
        "public" => publication == "public",
        "local_private" | "private" => publication == "local_private",
        other => {
            publication == "public"
                && n.attrs.get("content_scope").and_then(Value::as_str) == Some(other)
        }
    }
}
fn outline_ids(g: &KnowledgeGraph, n: &Node) -> Vec<String> {
    if let Some(ids) = n.attrs.get("ordered_section_ids").and_then(Value::as_array) {
        return ids
            .iter()
            .filter_map(Value::as_str)
            .map(str::to_owned)
            .collect();
    }
    g.edges
        .iter()
        .filter(|e| e.src == n.id && e.kind == "contains_section")
        .map(|e| e.dst.clone())
        .collect()
}
fn required<'a>(a: &'a Value, key: &str) -> Result<&'a str> {
    a.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| anyhow!("missing string argument {key}"))
}
fn number(a: &Value, key: &str, default: usize) -> Result<usize> {
    match a.get(key) {
        Some(Value::Null) | None => Ok(default),
        Some(x) => x
            .as_u64()
            .map(|n| n as usize)
            .ok_or_else(|| anyhow!("invalid {key}")),
    }
}
fn sized(v: &Value) -> usize {
    serde_json::to_string(v).unwrap_or_default().chars().count()
}
fn finalize(mut v: Value, cap: usize) -> Result<Value> {
    v["budget"] = json!({"max_chars":cap,"actual_chars":0,"truncated":false});
    for _ in 0..4 {
        let size = sized(&v);
        v["budget"]["actual_chars"] = json!(size);
    }
    if sized(&v) > cap {
        bail!("max_chars too small for paper response metadata; use a larger budget");
    }
    Ok(v)
}
fn page(
    tool: &str,
    items: Vec<Value>,
    offset: usize,
    limit: usize,
    cap: usize,
    extra: Value,
) -> Result<Value> {
    let total = items.len();
    let mut result = json!({"tool":tool,"items":[],"total":total,"offset":offset,"returned":0,"omitted":total.saturating_sub(offset),"next_offset":Value::Null,"reviewed":false,"notice":NOTICE,"context":extra});
    let mut chosen = Vec::new();
    for item in items.into_iter().skip(offset).take(limit) {
        chosen.push(item);
        result["items"] = json!(chosen);
        if sized(&result) + 250 > cap {
            chosen.pop();
            break;
        }
    }
    if chosen.is_empty() && offset < total {
        bail!("max_chars too small for the next paper item; increase max_chars");
    }
    let next = offset.saturating_add(chosen.len());
    result["items"] = json!(chosen);
    result["returned"] = json!(chosen.len());
    result["omitted"] = json!(total.saturating_sub(next));
    result["next_offset"] = if next < total {
        json!(next)
    } else {
        Value::Null
    };
    let mut result = finalize(result, cap)?;
    result["budget"]["truncated"] = json!(next < total);
    finalize_with_budget(result, cap)
}
fn finalize_with_budget(mut v: Value, cap: usize) -> Result<Value> {
    for _ in 0..4 {
        v["budget"]["actual_chars"] = json!(sized(&v));
    }
    if sized(&v) > cap {
        bail!("max_chars exceeded by paper metadata");
    }
    Ok(v)
}

pub fn coverage(g: &KnowledgeGraph, root: &Path) -> Value {
    let mut fresh = Freshness::new(g);
    let mut scopes = BTreeMap::<String, usize>::new();
    let mut states = BTreeMap::<String, usize>::new();
    let (
        mut docs,
        mut sections,
        mut current,
        mut stale,
        mut current_chunks,
        mut context_stale,
        mut context_incomplete,
    ) = (0, 0, 0, 0, 0, 0, 0);
    for n in g.nodes.iter().filter(|n| is_paper(n)) {
        if n.kind == "paper_document" {
            docs += 1;
            continue;
        }
        sections += 1;
        *scopes
            .entry(
                n.attrs
                    .get("scope")
                    .and_then(Value::as_str)
                    .unwrap_or("unknown")
                    .into(),
            )
            .or_default() += 1;
        let is_current = observed_freshness(n, g, root, &mut fresh)["current"] == true;
        if is_current {
            current += 1
        } else {
            stale += 1
        };
        let context_current = n.attrs.get("context_current") != Some(&Value::Bool(false));
        if n.attrs.get("context_state").and_then(Value::as_str) == Some("stale") {
            context_stale += 1;
        }
        if n.attrs.get("context_state").and_then(Value::as_str) == Some("incomplete") {
            context_incomplete += 1;
        }
        if is_current && context_current {
            current_chunks += n
                .attrs
                .pointer("/coverage/current_chunks")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize;
        }
        let state = n
            .attrs
            .get("summary_state")
            .and_then(Value::as_str)
            .unwrap_or("pending");
        *states
            .entry(
                if !is_current && state == "complete" {
                    "stale"
                } else {
                    state
                }
                .into(),
            )
            .or_default() += 1;
    }
    let imported = g.meta.pointer("/papers/import/coverage");
    let expected = imported
        .and_then(|v| v.pointer("/expected_registry_coverage/sections"))
        .cloned()
        .unwrap_or(Value::Null);
    let complete = *states.get("complete").unwrap_or(&0);
    json!({"documents":docs,"sections":sections,"imported_sections":sections,"expected_sections":expected,
        "physical_sources":imported.and_then(|v|v.get("physical_sources")),"public_sections":scopes.get("public").copied().unwrap_or(0),"private_sections":scopes.get("local_private").copied().unwrap_or(0),
        "sections_by_scope":scopes,"summary_states":states,"complete":complete,"pending":states.get("pending").copied().unwrap_or(0),"error":states.get("error").copied().unwrap_or(0),"stale":states.get("stale").copied().unwrap_or(0),"source_error":states.get("source_error").copied().unwrap_or(0),
        "current_chunks":current_chunks,"context_stale_sections":context_stale,"context_incomplete_sections":context_incomplete,"total_chunks":imported.and_then(|v|v.get("total_chunks")),
        "semantic_complete":sections>0&&complete==sections&&stale==0&&context_incomplete==0&&context_stale==0&&expected==sections,"source_current_sections":current,"source_stale_or_missing_sections":stale,"reviewed_sections":0,"complete_semantic_review":false,
        "evidence_entities":imported.and_then(|v|v.get("evidence_entities")),"evidence_source_errors":imported.and_then(|v|v.get("evidence_source_errors")),"document_source_errors":imported.and_then(|v|v.get("document_source_errors")),"science_unverified":true,
        "import_errors_count":g.meta.pointer("/papers/import/errors_count").and_then(Value::as_u64).unwrap_or(0),
        "embedded_paper_nodes":g.nodes.iter().filter(|n|is_paper(n)&&crate::service::embedding(n).is_some()).count(),"embedding_model":"text-embedding-3-small","embedding_dimensions":256,"state":g.meta.pointer("/papers/state"),"notice":NOTICE})
}

pub fn query(
    g: &KnowledgeGraph,
    root: &Path,
    tool: &str,
    a: &Value,
    vector: Option<&[f32]>,
) -> Result<Value> {
    let cap = number(a, "max_chars", 16000)?;
    if !(2000..=60000).contains(&cap) {
        bail!("max_chars must be between 2000 and 60000");
    }
    let offset = number(a, "offset", 0)?;
    let limit = number(a, "limit", if tool == "paper_outline" { 30 } else { 20 })?;
    if !(1..=100).contains(&limit) {
        bail!("limit must be between 1 and 100");
    }
    let scope = a.get("scope").and_then(Value::as_str).unwrap_or("public");
    if ![
        "public",
        "local_private",
        "private",
        "all",
        "active",
        "archive",
        "exploration",
        "historical",
        "support",
    ]
    .contains(&scope)
    {
        bail!("invalid paper scope");
    }
    let q = a.get("query").and_then(Value::as_str).unwrap_or("");
    if q.chars().count() > 2000 {
        bail!("query exceeds 2000 characters");
    }
    let mut fresh = Freshness::new(g);
    let by_id: HashMap<_, _> = g.nodes.iter().map(|n| (n.id.as_str(), n)).collect();
    match tool {
        "paper_coverage" => finalize(coverage(g, root), cap),
        "list_papers" | "search_papers" => {
            if tool == "search_papers" && q.trim().is_empty() {
                bail!("search_papers requires query");
            }
            let tokens = query_tokens(q);
            let document = a.get("document_id").and_then(Value::as_str);
            let allowed_ids = document.and_then(|id| {
                by_id
                    .get(id)
                    .filter(|n| n.kind == "paper_document" && is_paper(n))
                    .map(|n| outline_ids(g, n))
            });
            if document.is_some() && allowed_ids.is_none() {
                bail!("paper document not found");
            }
            let evidence = a.get("evidence_kind").and_then(Value::as_str);
            let mut scored = Vec::new();
            let mut vector_candidates = 0usize;
            let mut eligible_nodes = 0usize;
            for n in g.nodes.iter().filter(|n| {
                if tool == "list_papers" {
                    n.kind == "paper_document" && is_paper(n)
                } else {
                    is_paper(n)
                }
            }) {
                if !matches_scope(n, scope)
                    || (document.is_some()
                        && Some(n.id.as_str()) != document
                        && !allowed_ids.as_ref().is_some_and(|ids| ids.contains(&n.id)))
                {
                    continue;
                }
                if evidence.is_some_and(|kind| {
                    !g.edges.iter().any(|e| {
                        e.src == n.id
                            && (e.kind == kind
                                || by_id.get(e.dst.as_str()).is_some_and(|other| {
                                    other.kind == kind
                                        || other
                                            .attrs
                                            .get("evidence_type")
                                            .or(other.attrs.get("original_type"))
                                            .and_then(Value::as_str)
                                            == Some(kind)
                                }))
                    })
                }) {
                    continue;
                }
                eligible_nodes += 1;
                let current = observed_freshness(n, g, root, &mut fresh)["current"] == true;
                let text = searchable(n, current).to_lowercase();
                let heading = n.label.to_lowercase();
                let lexical: f32 = tokens.iter().map(|term| {
                    let aliases = lexical_aliases(term);
                    let matches = |field: &str| field.contains(term.as_str())
                        || aliases.iter().any(|alias| field.contains(alias));
                    if matches(&heading) { 2.0 } else if matches(&text) { 1.0 } else { 0.0 }
                }).sum();
                let embedding = if summary_current(n, current) {
                    crate::service::embedding(n)
                } else {
                    None
                };
                if embedding.is_some() {
                    vector_candidates += 1;
                }
                let semantic = vector
                    .zip(embedding.as_deref())
                    .map(|(v, e)| crate::service::cosine(v, e))
                    .unwrap_or(0.0);
                if !tokens.is_empty() && lexical == 0.0 && semantic <= 0.15 {
                    continue;
                }
                let score = if tokens.is_empty() {
                    0.0
                } else {
                    2.0 * lexical / tokens.len() as f32
                } + semantic.max(0.0) * 3.0;
                let historical = matches!(
                    n.attrs.get("content_scope").and_then(Value::as_str),
                    Some("archive" | "historical")
                );
                let mut item = row(n, g, root, &mut fresh);
                item["score"] = json!(score);
                scored.push((score, !current, historical, n.id.clone(), item));
            }
            scored.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)).then(x.3.cmp(&y.3)));
            page(
                tool,
                scored.into_iter().map(|x| x.4).collect(),
                offset,
                limit,
                cap,
                json!({"query":q,"scope":scope,"search_mode":if vector.is_some()&&vector_candidates>0{"hybrid"}else{"lexical_full_text"},"paper_vectors_available":vector_candidates,"semantic_candidates_scored":if vector.is_some(){vector_candidates}else{0},"eligible_nodes":eligible_nodes,"embedding_coverage_partial":vector_candidates<eligible_nodes,"full_source_text_searched":true}),
            )
        }
        "paper_outline" => {
            let id = required(a, "document_id")?;
            let n = by_id
                .get(id)
                .filter(|n| n.kind == "paper_document" && is_paper(n))
                .ok_or_else(|| anyhow!("paper document not found: {id}"))?;
            let items = outline_ids(g, n)
                .iter()
                .filter_map(|id| by_id.get(id.as_str()))
                .map(|n| row(n, g, root, &mut fresh))
                .collect();
            page(
                tool,
                items,
                offset,
                limit,
                cap,
                json!({"document":row(n,g,root,&mut fresh)}),
            )
        }
        "paper_section" => {
            let id = required(a, "section_id")?;
            let n = by_id
                .get(id)
                .filter(|n| n.kind == "paper_section" && is_paper(n))
                .ok_or_else(|| anyhow!("paper section not found: {id}"))?;
            let mut result = json!({"section":row(n,g,root,&mut fresh),"semantic_as_recorded":n.attrs.get("semantic"),"coverage":n.attrs.get("coverage"),"document_ids":n.attrs.get("document_ids"),"reviewed":false,"notice":NOTICE,"source_read_tool":"read_source(source_id, line_start?, line_end?)","omitted_fields":[]});
            if sized(&result) + 150 > cap {
                result["semantic_as_recorded"] = json!({"summary":clipped(&n.summary,cap.saturating_sub(sized(&result["section"])+1000))});
                result["omitted_fields"] =
                    json!(["semantic detail beyond summary; inspect the registered source"]);
            }
            let truncated = result["omitted_fields"]
                .as_array()
                .is_some_and(|a| !a.is_empty());
            let mut result = finalize(result, cap)?;
            result["budget"]["truncated"] = json!(truncated);
            finalize_with_budget(result, cap)
        }
        "paper_evidence" => {
            let id = required(a, "id")?;
            let n = by_id
                .get(id)
                .ok_or_else(|| anyhow!("paper/evidence ID not found: {id}"))?;
            if matches!(n.kind.as_str(), "paper_document" | "paper_section") && !is_paper(n) {
                bail!(
                    "legacy structural paper ID; use get_dossier or an imported registry paper ID"
                );
            }
            let mut ids = HashSet::from([id.to_owned()]);
            if n.kind == "paper_document" {
                ids.extend(outline_ids(g, n));
            }
            let filter = a.get("evidence_type").and_then(Value::as_str);
            let mut items = Vec::new();
            for e in &g.edges {
                if e.kind == "contains_section" {
                    continue;
                }
                let target = if ids.contains(&e.src) {
                    &e.dst
                } else if ids.contains(&e.dst) {
                    &e.src
                } else {
                    continue;
                };
                let Some(other) = by_id.get(target.as_str()) else {
                    continue;
                };
                if filter.is_some_and(|f| {
                    e.kind != f
                        && other.kind != f
                        && other
                            .attrs
                            .get("evidence_type")
                            .or(other.attrs.get("original_type"))
                            .and_then(Value::as_str)
                            != Some(f)
                }) {
                    continue;
                }
                items.push(json!({"relation":{"src":e.src,"dst":e.dst,"type":e.kind,"match_basis":e.attrs.get("match_basis"),"confidence":e.attrs.get("confidence"),"scope":e.attrs.get("scope"),"reviewed":false},"evidence":{"id":other.id,"type":other.kind,"label":clipped(&other.label,240),"path":other.attrs.get("path"),"evidence_type":other.attrs.get("evidence_type").or(other.attrs.get("original_type")),"status_as_recorded":other.attrs.get("status_as_recorded"),"scope":other.attrs.get("publication_scope").or(other.attrs.get("scope")),"source_current_as_recorded":other.attrs.get("source_current"),"source_freshness":fresh.check(root,g,other)}}));
            }
            page(
                tool,
                items,
                offset,
                limit,
                cap,
                json!({"id":id,"scope":n.attrs.get("publication_scope").or(n.attrs.get("scope")),"content_scope":n.attrs.get("content_scope"),"paper_layer":is_paper(n),"evidence_type":filter,"notice":"Literal citation/navigation relations and AI candidates retain their separate stored types; none establishes a proof."}),
            )
        }
        _ => bail!("unknown paper tool: {tool}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::service::{GraphStore, SearchRequest};
    use std::{
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };
    struct Fixture {
        root: PathBuf,
        graph: KnowledgeGraph,
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }
    impl Fixture {
        fn new() -> Self {
            static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            let root = std::env::temp_dir().join(format!(
                "tfpt-paper-test-{}-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            fs::create_dir_all(&root).unwrap();
            let root = fs::canonicalize(root).unwrap();
            let text = format!(
                "\\section{{Original}}\n{}\ntrailing_covariant_formula\n",
                "ordinary source text ".repeat(1400)
            );
            fs::write(root.join("paper.tex"), &text).unwrap();
            let sha = format!("{:x}", Sha256::digest(text.as_bytes()));
            let source = json!({"path":"paper.tex","sha256":sha,"byte_start":0,"byte_end":text.len(),"span_sha256":sha,"line_start":1,"line_end":3});
            let section = |id: &str, scope: &str, content: &str| json!({"id":id,"type":"paper_section","label":"Native clock geometry","summary":"An intrinsic clock connection","attrs":{"scope":scope,"publication_scope":scope,"content_scope":content,"summary_state":"complete","source_current":true,"reviewed":false,"source":source,"source_refs":["paper-source:fixture"],"source_text":text,"semantic":{"summary":"Semantic clock relation","concepts":["geometric transport"]},"coverage":{"chunks_total":2,"current_chunks":2}}});
            let doc = |id: &str, sid: &str, scope: &str| json!({"id":id,"type":"paper_document","label":id,"summary":"Document context","attrs":{"scope":scope,"content_scope":"active","source_current":true,"summary_state":"complete","ordered_section_ids":[sid],"source_refs":["paper-source:fixture"]}});
            let mut graph:KnowledgeGraph=serde_json::from_value(json!({"meta":{"papers":{"state":"loaded","import":{"coverage":{"expected_registry_coverage":{"sections":3},"physical_sources":1,"current_chunks":6,"total_chunks":6}}}},"nodes":[section("section:active","public","active"),section("section:archive","public","archive"),section("section:private","local_private","private"),doc("document:active","section:active","public"),doc("document:private","section:private","local_private"),{"id":"claim:CHECK","type":"paper_evidence","label":"Original statement","attrs":{"evidence_type":"claim_ledger_row"}}],"edges":[{"src":"section:active","dst":"claim:CHECK","type":"literal_citation","attrs":{"confidence":"source_explicit","match_basis":"literal source ID"}},{"src":"section:active","dst":"claim:CHECK","type":"ai_candidate","attrs":{"confidence":"lexical","match_basis":"candidate only"}}],"sources":[{"id":"paper-source:fixture","path":"paper.tex","sha256":sha,"line_start":1,"line_end":3}]})).unwrap();
            for n in &mut graph.nodes {
                if matches!(n.kind.as_str(), "paper_document" | "paper_section") {
                    n.attrs["paper_layer"] = json!(true);
                }
            }
            Self { root, graph }
        }
        fn query(&self, tool: &str, args: Value) -> Value {
            query(&self.graph, &self.root, tool, &args, None).unwrap()
        }
        fn write_overlay(&self) -> PathBuf {
            let mut overlay = self.graph.clone();
            overlay.meta = json!({"coverage":{"expected_registry_coverage":{"sections":3}}});
            fs::write(
                self.root.join("papers.json"),
                serde_json::to_vec(&overlay).unwrap(),
            )
            .unwrap();
            let data = self.root.join("knowledge.json");
            fs::write(&data,serde_json::to_vec(&json!({"meta":{},"nodes":[{"id":"claim:CHECK","type":"claim","label":"Authoritative claim","attrs":{"ledger":{"status":"E"}}}],"edges":[],"sources":[]})).unwrap()).unwrap();
            data
        }
    }
    #[test]
    fn full_tail_semantics_scopes_and_tie_ranking() {
        let f = Fixture::new();
        let r = f.query(
            "search_papers",
            json!({"query":"trailing_covariant_formula","limit":10}),
        );
        assert_eq!(r["total"], 2);
        assert_eq!(r["items"][0]["id"], "section:active");
        assert_eq!(r["items"][0]["source_current"], true);
        assert_eq!(
            f.query("search_papers", json!({"query":"geometric transport"}))["total"],
            2
        );
        assert_eq!(
            f.query(
                "search_papers",
                json!({"query":"trailing_covariant_formula","scope":"local_private"})
            )["items"][0]["id"],
            "section:private"
        );
        assert_eq!(f.query("list_papers", json!({}))["total"], 1);
        assert_eq!(f.query("list_papers", json!({"scope":"all"}))["total"], 2);
        assert_eq!(
            query_tokens("Wie ist die geometrische Clock und ihre Wirkung?"),
            vec!["geometrische", "clock", "ihre", "wirkung"]
        );
    }
    #[test]
    fn bounded_pages_outline_evidence_and_large_section() {
        let mut f = Fixture::new();
        let first = f.query(
            "search_papers",
            json!({"query":"clock","scope":"all","limit":1,"max_chars":3000}),
        );
        assert_eq!(first["next_offset"], 1);
        assert_eq!(first["budget"]["actual_chars"], sized(&first));
        assert!(sized(&first) <= 3000);
        assert_eq!(first["budget"]["truncated"], true);
        let second = f.query(
            "search_papers",
            json!({"query":"clock","scope":"all","limit":1,"offset":1}),
        );
        assert_ne!(first["items"][0]["id"], second["items"][0]["id"]);
        assert_eq!(
            f.query("paper_outline", json!({"document_id":"document:active"}))["items"][0]["id"],
            "section:active"
        );
        let evidence = f.query("paper_evidence", json!({"id":"section:active"}));
        assert_eq!(evidence["total"], 2);
        assert_eq!(
            f.query(
                "paper_evidence",
                json!({"id":"section:active","evidence_type":"ai_candidate"})
            )["total"],
            1
        );
        f.graph.nodes[0].attrs["semantic"]["summary"] = json!("long quote ".repeat(10000));
        let section = f.query(
            "paper_section",
            json!({"section_id":"section:active","max_chars":3000}),
        );
        assert!(sized(&section) <= 3000);
        assert_eq!(section["budget"]["actual_chars"], sized(&section));
        assert_eq!(section["budget"]["truncated"], true);
    }
    #[test]
    fn coverage_keeps_evidence_source_warnings_separate_from_section_completion() {
        let mut f = Fixture::new();
        f.graph.meta["papers"]["import"]["coverage"]["evidence_entities"] = json!(6026);
        f.graph.meta["papers"]["import"]["coverage"]["evidence_source_errors"] = json!(11);
        f.graph.meta["papers"]["import"]["coverage"]["document_source_errors"] = json!(0);
        f.graph.meta["papers"]["import"]["evidence_errors"] =
            json!([{"private_error_blob":"not a coverage payload"}]);
        let report = f.query("paper_coverage", json!({}));
        assert_eq!(report["evidence_entities"], 6026);
        assert_eq!(report["evidence_source_errors"], 11);
        assert_eq!(report["document_source_errors"], 0);
        assert_eq!(report["science_unverified"], true);
        assert_eq!(report["semantic_complete"], true);
        assert_eq!(report["complete_semantic_review"], false);
        assert!(!report.to_string().contains("private_error_blob"));
    }

    #[test]
    fn legacy_structural_paper_nodes_do_not_inflate_imported_coverage() {
        let mut f = Fixture::new();
        f.graph.nodes.push(Node {
            id: "section:legacy".into(),
            kind: "paper_section".into(),
            label: "Legacy structural source".into(),
            summary: "clock".into(),
            attrs: json!({"scope":"public"}),
        });
        assert_eq!(f.query("paper_coverage", json!({}))["imported_sections"], 3);
        let r = f.query("search_papers", json!({"query":"clock","scope":"all"}));
        assert!(!r["items"]
            .as_array()
            .unwrap()
            .iter()
            .any(|n| n["id"] == "section:legacy"));
        assert!(query(
            &f.graph,
            &f.root,
            "paper_section",
            &json!({"section_id":"section:legacy"}),
            None
        )
        .is_err());
    }
    #[test]
    fn german_query_finds_current_english_body_without_annotation() {
        let mut f = Fixture::new();
        f.graph.nodes[0].label = "Clifford transfer kernel".into();
        f.graph.nodes[0].attrs["summary_state"] = json!("stale");
        let hits = f.query("search_papers", json!({"query":"Welchen Transferkern hat die Cliffordquelle?"}));
        assert_eq!(hits["items"][0]["id"], "section:active");
        assert_eq!(hits["items"][0]["summary_current"], false);
    }

    #[test]
    fn stale_paper_annotations_do_not_rank_or_display_as_current() {
        let mut f = Fixture::new();
        f.graph.nodes[0].summary = "expiredkeyword".into();
        f.graph.nodes[0].attrs["semantic"] = json!({"summary":"expiredkeyword"});
        f.graph.nodes[0].attrs["summary_state"] = json!("stale");
        let stale = f.query("search_papers", json!({"query":"expiredkeyword"}));
        assert!(stale["items"].as_array().unwrap().is_empty());
        let original = f.query("search_papers", json!({"query":"trailing_covariant_formula"}));
        let active = original["items"].as_array().unwrap().iter()
            .find(|item| item["id"] == "section:active").unwrap();
        assert_eq!(active["summary"], Value::Null);
        assert_eq!(active["source_current"], true);
        let direct = f.query("paper_section", json!({"section_id":"section:active"}));
        assert_eq!(direct["semantic_as_recorded"]["summary"], "expiredkeyword");
        fs::write(f.root.join("paper.tex"), "updated text").unwrap();
        let old_body = f.query("search_papers", json!({"query":"trailing_covariant_formula"}));
        assert!(old_body["items"].as_array().unwrap().is_empty());
    }

    #[test]
    fn freshness_span_tampering_and_private_source_security() {
        let mut f = Fixture::new();
        f.graph.nodes[0].attrs["source_text"] = json!("fabricated");
        assert_eq!(
            f.query("paper_section", json!({"section_id":"section:active"}))["section"]
                ["source_current"],
            false
        );
        f.graph.nodes[0].attrs["source"]["path"] = json!("../secret");
        assert_eq!(
            f.query("paper_section", json!({"section_id":"section:active"}))["section"]
                ["source_freshness"]["state"],
            "unsafe source path"
        );
        fs::write(f.root.join("paper.tex"), "changed source").unwrap();
        let r = f.query("paper_section", json!({"section_id":"section:archive"}));
        assert_eq!(r["section"]["summary_state"], "stale");
        assert_eq!(r["section"]["summary_current"], false);
    }
    #[tokio::test]
    async fn merged_graph_reload_existing_search_read_source_and_export() {
        let f = Fixture::new();
        let data = f.write_overlay();
        let store = GraphStore::new(data, f.root.clone()).unwrap();
        let exported = store.export_graph().await.unwrap();
        let claim = exported
            .nodes
            .iter()
            .find(|n| n.id == "claim:CHECK")
            .unwrap();
        assert_eq!(claim.kind, "claim");
        assert_eq!(claim.attrs["ledger"]["status"], "E");
        assert!(exported.nodes.iter().any(|n| n.id == "section:active"));
        let result = store
            .search(SearchRequest {
                query: "trailing_covariant_formula".into(),
                kind: Some("paper_section".into()),
                limit: Some(10),
            })
            .await
            .unwrap();
        assert_eq!(result["total"], 3);
        assert_eq!(result["results"][0]["source_current"], true);
        assert!(store
            .source("paper-source:fixture", Some(1), Some(1))
            .await
            .unwrap()["text"]
            .as_str()
            .unwrap()
            .contains("Original"));
        let mut overlay = f.graph.clone();
        overlay.nodes[0].label = "Changed paper title after reload".into();
        fs::write(
            f.root.join("papers.json"),
            serde_json::to_vec(&overlay).unwrap(),
        )
        .unwrap();
        let r = store
            .paper_tool("paper_section", &json!({"section_id":"section:active"}))
            .await
            .unwrap();
        assert_eq!(r["section"]["label"], "Changed paper title after reload");
        fs::write(f.root.join("paper.tex"), "mutated").unwrap();
        assert!(store
            .source("paper-source:fixture", None, None)
            .await
            .is_err());
        let graph = store.export_graph().await.unwrap();
        assert_eq!(
            graph
                .nodes
                .iter()
                .find(|n| n.id == "section:active")
                .unwrap()
                .attrs["summary_state"],
            "stale"
        );
    }
    #[test]
    fn hybrid_ranking_and_stale_embedding_fingerprint() {
        let mut f = Fixture::new();
        let mut vector = vec![0.0; 256];
        vector[0] = 1.0;
        f.graph.nodes[0].attrs["embedding"] = json!(vector);
        let r = query(
            &f.graph,
            &f.root,
            "search_papers",
            &json!({"query":"unrelated natural language","limit":1}),
            Some(&vector),
        )
        .unwrap();
        assert_eq!(r["context"]["search_mode"], "hybrid");
        assert_eq!(r["items"][0]["id"], "section:active");
        f.write_overlay();
        let mut g = KnowledgeGraph {
            meta: json!({}),
            ..Default::default()
        };
        merge(&mut g, &f.root.join("papers.json")).unwrap();
        assert!(g
            .nodes
            .iter()
            .find(|n| n.id == "section:active")
            .unwrap()
            .attrs
            .get("embedding")
            .is_none());
    }
}

#[cfg(test)]
mod context_tests {
    use super::*;
    #[test]
    fn missing_declared_context_keeps_available_body_completion_without_claiming_fresh_context() {
        let dir = std::env::temp_dir().join(format!(
            "tfpt-paper-incomplete-context-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let dir = fs::canonicalize(dir).unwrap();
        fs::write(dir.join("body.tex"), "available original body").unwrap();
        let body = format!("{:x}", Sha256::digest(b"available original body"));
        let mut g: KnowledgeGraph = serde_json::from_value(json!({
            "meta": {}, "edges": [],
            "nodes": [
                {"id":"section:available", "type":"paper_section", "attrs":{
                    "paper_layer":true, "summary_state":"complete", "source_refs":["source:body"],
                    "context_dependencies":[{"id":"document:historical"}]}},
                {"id":"document:historical", "type":"paper_document", "attrs":{
                    "paper_layer":true, "source_current":false,
                    "source_error":"missing declared original input: missing-original.tex",
                    "ordered_section_ids":["section:available"], "source_refs":["source:body"]}}
            ],
            "sources":[{"id":"source:body", "path":"body.tex", "sha256":body}]
        }))
        .unwrap();
        refresh_sources(&mut g, &dir);
        let section = &g.nodes[0].attrs;
        assert_eq!(section["source_current"], true);
        assert_eq!(section["summary_state"], "complete");
        assert_eq!(section["summary_state_as_recorded"], "complete");
        assert_eq!(section["context_current"], false);
        assert_eq!(section["context_complete"], false);
        assert_eq!(section["context_state"], "incomplete");
        assert!(section["context_error"]
            .to_string()
            .contains("missing-original.tex"));
        let report = coverage(&g, &dir);
        assert_eq!(report["source_current_sections"], 1);
        assert_eq!(report["context_stale_sections"], 0);
        assert_eq!(report["context_incomplete_sections"], 1);
        assert_eq!(report["semantic_complete"], false);
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn dependency_change_stales_context_even_when_section_bytes_are_current() {
        let dir = std::env::temp_dir().join(format!("tfpt-paper-context-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let dir = fs::canonicalize(dir).unwrap();
        fs::write(dir.join("body.tex"), "body").unwrap();
        fs::write(dir.join("preamble.tex"), "original macros").unwrap();
        let body = format!("{:x}", Sha256::digest(b"body"));
        let preamble = format!("{:x}", Sha256::digest(b"original macros"));
        let mut g:KnowledgeGraph=serde_json::from_value(json!({"meta":{},"nodes":[{"id":"section:test","type":"paper_section","attrs":{"summary_state":"complete","source_refs":["source:body"],"context_dependencies":[{"id":"document:test"}]}},{"id":"document:test","type":"paper_document","attrs":{"ordered_section_ids":["section:test"],"source_refs":["source:preamble"]}}],"edges":[],"sources":[{"id":"source:body","path":"body.tex","sha256":body},{"id":"source:preamble","path":"preamble.tex","sha256":preamble}]})).unwrap();
        for n in &mut g.nodes {
            n.attrs["paper_layer"] = json!(true);
        }
        fs::write(dir.join("preamble.tex"), "changed macros").unwrap();
        refresh_sources(&mut g, &dir);
        assert_eq!(g.nodes[0].attrs["source_current"], true);
        assert_eq!(g.nodes[0].attrs["context_current"], false);
        assert_eq!(g.nodes[0].attrs["summary_state"], "stale");
        let _ = fs::remove_dir_all(dir);
    }
}
