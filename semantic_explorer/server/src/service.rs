use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{HashMap, HashSet, VecDeque},
    fs,
    path::{Component, Path, PathBuf},
    sync::{Mutex, RwLock as StdRwLock},
    time::SystemTime,
};
use tokio::sync::RwLock;

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct KnowledgeGraph {
    #[serde(default)]
    pub meta: Value,
    #[serde(default)]
    pub nodes: Vec<Node>,
    #[serde(default)]
    pub edges: Vec<Edge>,
    #[serde(default)]
    pub sources: Vec<Source>,
}
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct Node {
    pub id: String,
    #[serde(rename = "type", default)]
    pub kind: String,
    #[serde(default)]
    pub label: String,
    #[serde(default)]
    pub summary: String,
    #[serde(default)]
    pub attrs: Value,
}
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct Edge {
    #[serde(default)]
    pub id: String,
    pub src: String,
    pub dst: String,
    #[serde(rename = "type", default)]
    pub kind: String,
    #[serde(default)]
    pub attrs: Value,
}
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct Source {
    pub id: String,
    #[serde(default)]
    pub path: String,
    #[serde(default)]
    pub sha256: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub line_start: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub line_end: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub line_count: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub page: Option<usize>,
    #[serde(default)]
    pub text: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SearchRequest {
    #[serde(default)]
    pub query: String,
    #[serde(rename = "type")]
    pub kind: Option<String>,
    pub limit: Option<usize>,
}

pub struct GraphStore {
    data_path: PathBuf,
    context_path: PathBuf,
    translations_path: PathBuf,
    quality_path: PathBuf,
    corpus_path: PathBuf,
    embeddings_cache_path: PathBuf,
    orientation_path: PathBuf,
    repo_root: PathBuf,
    graph: RwLock<KnowledgeGraph>,
    stamp: Mutex<(
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
        Option<(SystemTime, u64)>,
    )>,
    embedding_cache: StdRwLock<HashMap<String, Vec<f32>>>,
}

impl GraphStore {
    pub fn new(data_path: PathBuf, repo_root: PathBuf) -> Result<Self> {
        let data_path = absolute_path(data_path)?;
        let context_path = data_path
            .parent()
            .unwrap_or(Path::new("."))
            .join("context.json");
        let translations_path = data_path
            .parent()
            .unwrap_or(Path::new("."))
            .join("translations.json");
        let quality_path = data_path
            .parent()
            .unwrap_or(Path::new("."))
            .join("quality.json");
        let corpus_path = data_path
            .parent()
            .unwrap_or(Path::new("."))
            .join("corpus.json");
        let embeddings_cache_path = data_path
            .parent()
            .unwrap_or(Path::new("."))
            .join("embeddings-cache.json");
        let orientation_path = data_path
            .parent()
            .and_then(Path::parent)
            .unwrap_or(Path::new("."))
            .join("orientation.json");
        let repo_root = fs::canonicalize(&repo_root)
            .with_context(|| format!("repository root unavailable: {}", repo_root.display()))?;
        if !data_path.exists() {
            bail!("knowledge graph unavailable: {}", data_path.display());
        }
        let graph = load_merged_graph(
            &data_path,
            &context_path,
            &translations_path,
            &quality_path,
            &corpus_path,
            &embeddings_cache_path,
            &orientation_path,
        )?;
        let stamp = (
            file_stamp(&data_path),
            file_stamp(&context_path),
            file_stamp(&translations_path),
            file_stamp(&quality_path),
            file_stamp(&corpus_path),
            file_stamp(&embeddings_cache_path),
            file_stamp(&orientation_path),
            file_stamp(&data_path.with_file_name("papers.json")),
        );
        Ok(Self {
            data_path,
            context_path,
            translations_path,
            quality_path,
            corpus_path,
            embeddings_cache_path,
            orientation_path,
            repo_root,
            graph: RwLock::new(graph),
            stamp: Mutex::new(stamp),
            embedding_cache: StdRwLock::new(HashMap::new()),
        })
    }

    async fn snapshot(&self) -> Result<KnowledgeGraph> {
        let modified = (
            file_stamp(&self.data_path),
            file_stamp(&self.context_path),
            file_stamp(&self.translations_path),
            file_stamp(&self.quality_path),
            file_stamp(&self.corpus_path),
            file_stamp(&self.embeddings_cache_path),
            file_stamp(&self.orientation_path),
            file_stamp(&self.data_path.with_file_name("papers.json")),
        );
        let current = self
            .stamp
            .lock()
            .map_err(|_| anyhow!("graph reload lock poisoned"))?
            .clone();
        if modified != current {
            let mut graph = self.graph.write().await;
            let current = self
                .stamp
                .lock()
                .map_err(|_| anyhow!("graph reload lock poisoned"))?
                .clone();
            let modified = (
                file_stamp(&self.data_path),
                file_stamp(&self.context_path),
                file_stamp(&self.translations_path),
                file_stamp(&self.quality_path),
                file_stamp(&self.corpus_path),
                file_stamp(&self.embeddings_cache_path),
                file_stamp(&self.orientation_path),
                file_stamp(&self.data_path.with_file_name("papers.json")),
            );
            if modified != current {
                *graph = load_merged_graph(
                    &self.data_path,
                    &self.context_path,
                    &self.translations_path,
                    &self.quality_path,
                    &self.corpus_path,
                    &self.embeddings_cache_path,
                    &self.orientation_path,
                )?;
                *self
                    .stamp
                    .lock()
                    .map_err(|_| anyhow!("graph reload lock poisoned"))? = modified;
            }
        }
        let mut graph = self.graph.read().await.clone();
        crate::papers::refresh_sources(&mut graph, &self.repo_root);
        match crate::research::load(&self.research_path()) {
            Ok(program) => crate::research::augment(&mut graph, &program, &self.repo_root),
            Err(error) => {
                graph.meta["research_program"] = json!({"state":"invalid","error":error.to_string(),"notice":"Research decisions unavailable; original corpus remains accessible. No decision should be inferred from this failure."})
            }
        }
        Ok(graph)
    }

    fn research_path(&self) -> PathBuf {
        self.data_path
            .parent()
            .and_then(Path::parent)
            .unwrap_or(Path::new("."))
            .join("research/program.json")
    }
    pub async fn research_state(
        &self,
        query: &str,
        lang: &str,
        max_chars: usize,
        offset: usize,
        detail_offset: usize,
        expected_revision: Option<&str>,
    ) -> Result<Value> {
        let g = self.snapshot().await?;
        let p = crate::research::load(&self.research_path())?;
        if expected_revision.is_some_and(|expected| expected != p.revision) {
            bail!(
                "research revision changed; restart pagination from detail_offset=0 (expected {}, current {})",
                expected_revision.unwrap(),
                p.revision
            );
        }
        crate::research::state(
            &p,
            &g,
            &self.repo_root,
            query,
            lang,
            max_chars,
            offset,
            detail_offset,
        )
    }
    pub async fn research_record(&self, id: &str, max_chars: usize) -> Result<Value> {
        crate::research::editable_record(
            &crate::research::load(&self.research_path())?,
            id,
            max_chars,
        )
    }
    pub async fn assess_research_plan(
        &self,
        problem: &str,
        strategy: &str,
        premises: &[String],
        lang: &str,
    ) -> Result<Value> {
        let p = crate::research::load(&self.research_path())?;
        crate::research::assess(&p, &self.repo_root, problem, strategy, premises, lang)
    }
    pub async fn record_research_update(
        &self,
        expected: &str,
        problem: Value,
        premises: Vec<crate::research::Premise>,
        reason: &str,
    ) -> Result<Value> {
        let g = self.snapshot().await?;
        crate::research::update(
            &self.research_path(),
            &self.repo_root,
            &g,
            expected,
            problem,
            premises,
            reason,
        )
    }

    pub async fn paper_tool(&self, tool: &str, args: &Value) -> Result<Value> {
        let g = self.snapshot().await?;
        let query = args.get("query").and_then(Value::as_str).unwrap_or("");
        if query.chars().count() > 2000 {
            bail!("query exceeds 2000 characters");
        }
        let vector = if tool == "search_papers"
            && !query.trim().is_empty()
            && g.nodes.iter().any(|n| {
                crate::papers::is_paper(n)
                    && n.attrs.get("source_current") == Some(&Value::Bool(true))
                    && n.attrs.get("context_current") != Some(&Value::Bool(false))
                    && embedding(n).is_some()
            }) {
            self.query_embedding(query).await.ok()
        } else {
            None
        };
        crate::papers::query(&g, &self.repo_root, tool, args, vector.as_deref())
    }

    pub async fn status(&self) -> Result<Value> {
        let g = self.snapshot().await?;
        let enriched_claims: Vec<_> = g
            .nodes
            .iter()
            .filter(|n| {
                n.kind == "claim"
                    && n.attrs.get("enrichment_state").and_then(Value::as_str) == Some("enriched")
            })
            .collect();
        let (
            mut hardcoded,
            mut free_parameter,
            mut upstream_result,
            mut declared_axiom,
            mut external_data,
            mut derived_in_script,
            mut unknown,
            mut other_origin,
        ) = (
            0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize,
        );
        let (mut conditional, mut circularity_risk, mut missing_derivation) =
            (0usize, 0usize, 0usize);
        let mut proof_gaps = 0usize;
        let mut assessments: HashMap<String, usize> = HashMap::new();
        for claim in g.nodes.iter().filter(|n| n.kind == "claim") {
            let category = claim
                .attrs
                .pointer("/assessment/category")
                .and_then(Value::as_str)
                .unwrap_or("unknown");
            *assessments.entry(category.to_owned()).or_default() += 1;
        }
        let claims_total = g.nodes.iter().filter(|n| n.kind == "claim").count();
        let claims_enriched = g
            .nodes
            .iter()
            .filter(|n| {
                n.kind == "claim"
                    && n.attrs.get("enrichment_state").and_then(Value::as_str) == Some("enriched")
            })
            .count();
        let claims_pending = g
            .nodes
            .iter()
            .filter(|n| {
                n.kind == "claim"
                    && n.attrs.get("enrichment_state").and_then(Value::as_str) == Some("pending")
            })
            .count();
        let claims_error = g
            .nodes
            .iter()
            .filter(|n| {
                n.kind == "claim"
                    && n.attrs.get("enrichment_state").and_then(Value::as_str) == Some("error")
            })
            .count();
        let mut meta = g.meta.clone();
        meta["paper_coverage"] = crate::papers::coverage(&g, &self.repo_root);
        let current_claim_states: HashMap<&str, &str> = g
            .nodes
            .iter()
            .filter(|node| node.kind == "claim")
            .filter_map(|node| {
                node.attrs
                    .get("enrichment_state")
                    .and_then(Value::as_str)
                    .map(|state| (node.id.as_str(), state))
            })
            .collect();
        let mut active_errors = 0usize;
        let mut superseded_errors = 0usize;
        let mut unlinked_errors = 0usize;
        let mut error_examples = Vec::new();
        let recorded_errors = meta
            .get("errors")
            .and_then(Value::as_array)
            .cloned()
            .unwrap_or_default();
        for error in &recorded_errors {
            let raw_id = error.get("claim_id").and_then(Value::as_str).unwrap_or("");
            let claim_id = if raw_id.starts_with("claim:") {
                raw_id.to_owned()
            } else {
                format!("claim:{raw_id}")
            };
            match current_claim_states.get(claim_id.as_str()).copied() {
                Some("enriched") => superseded_errors += 1,
                Some("error") | Some("pending") => {
                    active_errors += 1;
                    if error_examples.len() < 10 {
                        error_examples.push(json!({
                            "claim_id":claim_id,
                            "state":current_claim_states.get(claim_id.as_str()),
                            "error":clip(&value_text(error.get("error")),220)
                        }));
                    }
                }
                Some(_) => superseded_errors += 1,
                None => unlinked_errors += 1,
            }
        }
        if let Some(meta_object) = meta.as_object_mut() {
            meta_object.remove("errors");
            meta_object.insert(
                "error_summary".to_owned(),
                json!({
                    "recorded_total":recorded_errors.len(),
                    "unresolved_current_claims":active_errors,
                    "superseded_by_current_enrichment":superseded_errors,
                    "unlinked_records":unlinked_errors,
                    "examples":error_examples,
                    "notice":"Raw extractor errors are omitted. Entries for claims now marked enriched are counted as superseded metadata, not current failures."
                }),
            );
            if let Some(integrity) = g.meta.get("graph_integrity") {
                meta_object.insert("graph_integrity".into(), json!({
                    "unresolved_semantic_links":integrity.get("unresolved_semantic_links").cloned().unwrap_or(Value::Null),
                    "issue_count":integrity.get("issues").and_then(Value::as_array).map(Vec::len).unwrap_or(0),
                    "notice":integrity.get("notice").cloned().unwrap_or_else(||json!("Some referenced semantic links may not resolve to graph nodes."))
                }));
            }
        }
        if let Some(coverage) = meta.get_mut("coverage").and_then(Value::as_object_mut) {
            coverage.insert("claims_total".into(), json!(claims_total));
            coverage.insert("claims_enriched".into(), json!(claims_enriched));
            coverage.insert("claims_pending".into(), json!(claims_pending));
            coverage.insert("claims_error".into(), json!(claims_error));
        }
        for claim in &enriched_claims {
            let semantic = claim.attrs.get("semantic").unwrap_or(&Value::Null);
            if let Some(entries) = semantic.get("input_audit").and_then(Value::as_array) {
                for entry in entries {
                    match entry.get("origin").and_then(Value::as_str).unwrap_or("") {
                        "hardcoded" => hardcoded += 1,
                        "free_parameter" => free_parameter += 1,
                        "upstream_result" => upstream_result += 1,
                        "declared_axiom" => declared_axiom += 1,
                        "external_data" => external_data += 1,
                        "derived_in_script" => derived_in_script += 1,
                        "unknown" => unknown += 1,
                        _ => other_origin += 1,
                    }
                    match entry.get("impact") {
                        Some(Value::String(v)) => count_impact(
                            v,
                            &mut conditional,
                            &mut circularity_risk,
                            &mut missing_derivation,
                        ),
                        Some(Value::Array(vs)) => {
                            for v in vs {
                                if let Some(v) = v.as_str() {
                                    count_impact(
                                        v,
                                        &mut conditional,
                                        &mut circularity_risk,
                                        &mut missing_derivation,
                                    );
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
            proof_gaps += value_count(semantic.get("proof_gaps"));
        }
        Ok(
            json!({"meta":meta,"counts":{"nodes":g.nodes.len(),"edges":g.edges.len()},"source_count":g.sources.len(),"source_freshness":"checked_when_read","assessment_counts":assessments,"audit":{"enriched_claims":enriched_claims.len(),"input_audit_entries":{"total":hardcoded+free_parameter+upstream_result+declared_axiom+external_data+derived_in_script+unknown+other_origin,"by_origin":{"hardcoded":hardcoded,"free_parameter":free_parameter,"upstream_result":upstream_result,"declared_axiom":declared_axiom,"external_data":external_data,"derived_in_script":derived_in_script,"unknown":unknown,"other":other_origin}},"by_impact":{"conditional":conditional,"circularity_risk":circularity_risk,"missing_derivation":missing_derivation},"proof_gaps":proof_gaps}}),
        )
    }

    /// Export the current merged Atlas graph for legacy compatibility readers.
    /// This is a read-only snapshot; it is never written back to disk here.
    pub async fn export_graph(&self) -> Result<KnowledgeGraph> {
        self.snapshot().await
    }
    pub async fn audit_inputs(&self, id: &str) -> Result<Value> {
        const MAX_RELATED_CLAIMS: usize = 20;
        let g = self.snapshot().await?;
        let node = g
            .nodes
            .iter()
            .find(|n| n.id == id)
            .ok_or_else(|| anyhow!("node not found: {id}"))?;
        let semantic = node
            .attrs
            .get("semantic")
            .cloned()
            .unwrap_or_else(|| json!({}));
        match node.kind.as_str() {
            "claim" => Ok(
                json!({"node":public_node(node),"input_audit":semantic.get("input_audit").cloned().unwrap_or_else(||json!([])),"proof_gaps":semantic.get("proof_gaps").cloned().unwrap_or_else(||json!([]))}),
            ),
            "script" => {
                let candidates = node
                    .attrs
                    .pointer("/static_inputs/syntactic_input_candidates")
                    .or_else(|| semantic.get("syntactic_input_candidates"))
                    .or_else(|| {
                        semantic
                            .get("script_analysis")
                            .and_then(|v| v.get("syntactic_input_candidates"))
                    })
                    .cloned()
                    .unwrap_or_else(|| json!([]));
                let by_id: HashMap<_, _> = g.nodes.iter().map(|n| (n.id.as_str(), n)).collect();
                let mut seen = HashSet::from([id.to_owned()]);
                let mut frontier = vec![id.to_owned()];
                let mut claims = Vec::new();
                for _ in 0..3 {
                    if frontier.is_empty() {
                        break;
                    }
                    let mut next_frontier = Vec::new();
                    for current in frontier {
                        for edge in &g.edges {
                            let next = if edge.src == current {
                                Some(edge.dst.as_str())
                            } else if edge.dst == current {
                                Some(edge.src.as_str())
                            } else {
                                None
                            };
                            let Some(next) = next else { continue };
                            if !seen.insert(next.to_owned()) {
                                continue;
                            }
                            if let Some(candidate) = by_id.get(next) {
                                if candidate.kind == "claim" {
                                    claims.push(json!({"id":candidate.id,"label":candidate.label,"enrichment_state":candidate.attrs.get("enrichment_state").cloned().unwrap_or(Value::Null),"input_audit":candidate.attrs.pointer("/semantic/input_audit").cloned().unwrap_or_else(||json!([])),"proof_gaps":candidate.attrs.pointer("/semantic/proof_gaps").cloned().unwrap_or_else(||json!([]))}));
                                    if claims.len() >= MAX_RELATED_CLAIMS {
                                        break;
                                    }
                                } else {
                                    next_frontier.push(next.to_owned());
                                }
                            }
                        }
                        if claims.len() >= MAX_RELATED_CLAIMS {
                            break;
                        }
                    }
                    frontier = next_frontier;
                }
                Ok(
                    json!({"node":public_node(node),"syntactic_input_candidates":candidates,"related_claims":claims,"interpretation_note":"Syntactic constants and input candidates are signals for review; they are not proof gaps by themselves."}),
                )
            }
            _ => bail!("audit_inputs supports claim and script nodes"),
        }
    }
    pub async fn topics(&self) -> Result<Value> {
        let g = self.snapshot().await?;
        let topics: Vec<_> = g
            .nodes
            .iter()
            .filter(|n| n.kind == "topic")
            .map(public_node)
            .collect();
        Ok(json!({"topics":topics}))
    }
    pub async fn theory_overview(&self) -> Result<Value> {
        let map = match self.read_theory_map().await {
            Ok(map) => map,
            Err(error)
                if error
                    .to_string()
                    .contains("Theory map has not been generated yet") =>
            {
                return Ok(
                    json!({"available":false,"path":"theory_map.json","notice":"Theory map has not been generated yet.", "evidence_areas":crate::overview::areas("de", true)}),
                )
            }
            Err(error) => return Err(error),
        };
        Ok(compact_theory_map(&map, 20))
    }
    pub async fn theory_map(&self) -> Result<Value> {
        self.read_theory_map().await
    }
    async fn read_theory_map(&self) -> Result<Value> {
        let path = self
            .data_path
            .parent()
            .and_then(Path::parent)
            .unwrap_or(Path::new("."))
            .join("theory_map.json");
        let bytes = match fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                bail!("Theory map has not been generated yet")
            }
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("cannot read theory map: {}", path.display()))
            }
        };
        if bytes.len() > 64 * 1024 * 1024 {
            bail!("theory map exceeds 64 MiB safe read limit");
        }
        let map: Value = serde_json::from_slice(&bytes).context("invalid theory_map.json")?;
        if !map.is_object()
            || !map.get("nodes").is_some_and(Value::is_array)
            || !map.get("edges").is_some_and(Value::is_array)
        {
            bail!("theory map must contain node and edge arrays");
        }
        Ok(map)
    }
    pub async fn research_brief(&self, id: &str, max_chars: usize, offset: usize) -> Result<Value> {
        let max_chars = bounded_char_limit(max_chars)?;
        let graph = self.snapshot().await?;
        let node = graph
            .nodes
            .iter()
            .find(|node| node.id == id)
            .ok_or_else(|| anyhow!("node not found: {id}"))?;
        if !matches!(
            node.kind.as_str(),
            "topic" | "claim" | "goal" | "attempt" | "contract"
        ) {
            bail!("research_brief requires a topic, claim, goal, attempt, or contract id");
        }
        let by_id: HashMap<&str, &Node> = graph
            .nodes
            .iter()
            .map(|node| (node.id.as_str(), node))
            .collect();
        if matches!(node.kind.as_str(), "attempt" | "contract") {
            let mut claim_ids = HashSet::new();
            for edge in graph.edges.iter().filter(|edge| edge.kind == "tests") {
                let other = if edge.src == id {
                    Some(edge.dst.as_str())
                } else if edge.dst == id {
                    Some(edge.src.as_str())
                } else {
                    None
                };
                if let Some(claim) = other
                    .and_then(|other| by_id.get(other))
                    .filter(|claim| claim.kind == "claim")
                {
                    claim_ids.insert(claim.id.clone());
                }
            }
            for container in [Some(&node.attrs), node.attrs.get("original_contract")] {
                for key in ["claims_tested", "claim_ids"] {
                    if let Some(values) = container
                        .and_then(|attrs| attrs.get(key))
                        .and_then(Value::as_array)
                    {
                        claim_ids.extend(
                            values
                                .iter()
                                .filter_map(Value::as_str)
                                .filter_map(|value| normalize_claim_ref(value, &by_id)),
                        );
                    }
                }
            }
            let mut claims: Vec<&Node> = claim_ids
                .iter()
                .filter_map(|claim_id| by_id.get(claim_id.as_str()).copied())
                .collect();
            claims.sort_by(|a, b| a.id.cmp(&b.id));
            let linked_claims: Vec<Value> = claims.iter().take(6).map(|claim| json!({
                "id":claim.id,
                "label":clip(&claim.label,100),
                "assessment":assessment_category(claim).unwrap_or("unknown"),
                "scope":assessment_scope(claim).map(|scope|clip(scope,120)),
                "enrichment_state":claim.attrs.get("enrichment_state"),
                "summary":clip(claim.attrs.pointer("/semantic/overview_de").and_then(Value::as_str).unwrap_or(&claim.summary),150),
                "source_refs":source_path_refs(&graph,node_source_ids(claim),2)
            })).collect();
            let mut card = attempt_card(node, &node_source_ids(node));
            remove_duplicate_findings(&mut card);
            let own_source_ids = node_source_ids(node);
            card["source_refs"] = Value::Array(source_path_refs(&graph, own_source_ids.clone(), 5));
            card["source_refs_omitted"] = json!(own_source_ids
                .len()
                .saturating_sub(card["source_refs"].as_array().map(Vec::len).unwrap_or(0)));
            let typed_related_records = direct_typed_related_records(node, &graph, &by_id);
            let total_typed_records = typed_related_records.len();
            let mut direct_connections: Vec<Value> = graph.edges.iter()
                .filter(|edge| edge.src == id || edge.dst == id)
                .filter_map(|edge| {
                    let other = if edge.src == id { &edge.dst } else { &edge.src };
                    let target = by_id.get(other.as_str())?;
                    let mut connection = json!({"src":edge.src,"dst":edge.dst,"edge_type":edge.kind,"target_type":target.kind});
                    if let Some(reason) = edge.attrs.get("reason").and_then(Value::as_str) {
                        connection["reason"] = json!(clip(reason, 260));
                    }
                    if let Some(conditions) = edge.attrs.get("conditions").and_then(Value::as_array) {
                        connection["conditions"] = Value::Array(conditions.iter().take(2).map(|condition| json!(clip(&value_text(Some(condition)), 120))).collect());
                    }
                    if let Some(source_refs) = edge.attrs.get("source_refs").and_then(Value::as_array) {
                        connection["source_ids"] = Value::Array(source_refs.iter().filter_map(Value::as_str).take(3).map(|source| json!(source)).collect());
                    }
                    Some(connection)
                }).collect();
            direct_connections.sort_by(|a, b| {
                a["edge_type"]
                    .as_str()
                    .cmp(&b["edge_type"].as_str())
                    .then_with(|| a["dst"].as_str().cmp(&b["dst"].as_str()))
            });
            let total_connections = direct_connections.len();
            direct_connections.truncate(8);
            let mut source_ids = own_source_ids;
            source_ids.sort();
            source_ids.dedup();
            let total_sources = source_ids.len();
            source_ids.truncate(5);
            let what_we_have = attempt_what_we_have(node, &card);
            let mut result = json!({
                "id":node.id,
                "kind":node.kind,
                "question":clip(&node.label,220),
                "verdict":attempt_verdict(node),
                "scope":if attempt_scope(node).is_empty(){Value::Null}else{json!(clip(&attempt_scope(node),300))},
                "scope_recorded":!attempt_scope(node).is_empty(),
                "attempt_record":card,
                "what_we_have":what_we_have,
                "typed_related_records":typed_related_records,
                "linked_claims":{"items":linked_claims,"total":claims.len(),"omitted":claims.len().saturating_sub(6),"membership":"Only direct tests edges and explicit claims_tested/claim_ids are included; no topic expansion."},
                "counts":{"claims":claims.len(),"attempts":1,"typed_related_records":total_typed_records,"recorded_connections":total_connections,"source_ids":total_sources},
                "proof_gap_coverage":{"recorded_gaps_from_enriched_claims":Value::Null,"claims_without_enriched_gap_audit":Value::Null,"complete":Value::Null,"scope":if claims.is_empty(){"not_applicable_no_claim_membership"}else{"not_a_full_claim_audit"}},
                "recorded_connections":direct_connections,
                "source_refs":source_path_refs(&graph,source_ids.iter().cloned(),5),
                "source_ids":source_ids.clone(),
                "source_ids_omitted":total_sources.saturating_sub(5),
                "notice":"Exact attempt/contract evidence and direct typed links only. Stored test results and relations are not a proof of a linked physical claim.",
                "max_chars":max_chars,
                "truncated":false
            });
            compact_attempt_research_brief(&mut result);
            return Ok(fit_char_budget(result, max_chars));
        }
        let related_claims_total = if node.kind == "topic" {
            let (all_topic_claims, _, _, _) = topic_members(&graph, id, &by_id, true);
            all_topic_claims
                .iter()
                .filter_map(|claim_id| by_id.get(claim_id.as_str()).copied())
                .filter(|claim| claim_topic_role(claim, id) == Some("related"))
                .count()
        } else {
            0
        };
        let (claim_ids, attempt_ids, structure_ids, goal_ids) = if node.kind == "topic" {
            topic_members(&graph, id, &by_id, false)
        } else if matches!(node.kind.as_str(), "attempt" | "contract") {
            let mut claims = HashSet::new();
            for edge in &graph.edges {
                if edge.kind != "tests" {
                    continue;
                }
                let other = if edge.src == id {
                    Some(edge.dst.as_str())
                } else if edge.dst == id {
                    Some(edge.src.as_str())
                } else {
                    None
                };
                if let Some(claim_id) =
                    other.filter(|other| by_id.get(*other).is_some_and(|n| n.kind == "claim"))
                {
                    claims.insert(claim_id.to_owned());
                }
            }
            for key in ["claims_tested", "claim_ids"] {
                if let Some(values) = node.attrs.get(key).and_then(Value::as_array).or_else(|| {
                    node.attrs
                        .pointer(&format!("/original_contract/{key}"))
                        .and_then(Value::as_array)
                }) {
                    claims.extend(
                        values
                            .iter()
                            .filter_map(Value::as_str)
                            .filter_map(|value| normalize_claim_ref(value, &by_id)),
                    );
                }
            }
            let (structures, goals) = linked_structures_and_goals_typed(&graph, id, &by_id);
            (claims, HashSet::from([id.to_owned()]), structures, goals)
        } else if node.kind == "goal" {
            let claim_ids = goal_claim_ids(&graph, node, &by_id);
            let mut attempts = HashSet::new();
            for edge in &graph.edges {
                if edge.kind == "tests" {
                    let pair = if claim_ids.contains(&edge.dst) {
                        Some((&edge.src, &edge.dst))
                    } else if claim_ids.contains(&edge.src) {
                        Some((&edge.dst, &edge.src))
                    } else {
                        None
                    };
                    if let Some((attempt_id, _)) = pair {
                        if by_id
                            .get(attempt_id.as_str())
                            .is_some_and(|other| other.kind == "attempt")
                        {
                            attempts.insert(attempt_id.clone());
                        }
                    }
                } else if edge.kind == "leaves_open" && (edge.src == id || edge.dst == id) {
                    let attempt_id = if edge.src == id { &edge.dst } else { &edge.src };
                    if by_id
                        .get(attempt_id.as_str())
                        .is_some_and(|other| other.kind == "attempt")
                    {
                        attempts.insert(attempt_id.clone());
                    }
                }
            }
            let (structures, _) = linked_structures_and_goals(&graph, id, &by_id);
            let goals = HashSet::from([id.to_owned()]);
            (claim_ids, attempts, structures, goals)
        } else {
            let claim_ids = HashSet::from([id.to_owned()]);
            let mut attempts = HashSet::new();
            for edge in &graph.edges {
                if edge.kind == "tests" {
                    let other = if edge.src == id {
                        Some(edge.dst.as_str())
                    } else if edge.dst == id {
                        Some(edge.src.as_str())
                    } else {
                        None
                    };
                    if other
                        .and_then(|other| by_id.get(other))
                        .is_some_and(|other| other.kind == "attempt")
                    {
                        if let Some(other_id) = other {
                            attempts.insert(other_id.to_owned());
                        }
                    }
                }
            }
            let (structures, goals) = linked_structures_and_goals(&graph, id, &by_id);
            (claim_ids, attempts, structures, goals)
        };

        let mut claims: Vec<&Node> = claim_ids
            .iter()
            .filter_map(|id| by_id.get(id.as_str()).copied())
            .collect();
        let mut attempts: Vec<&Node> = attempt_ids
            .iter()
            .filter_map(|id| by_id.get(id.as_str()).copied())
            .collect();
        let mut structures: Vec<&Node> = structure_ids
            .iter()
            .filter_map(|id| by_id.get(id.as_str()).copied())
            .collect();
        let mut goals: Vec<&Node> = goal_ids
            .iter()
            .filter_map(|id| by_id.get(id.as_str()).copied())
            .collect();
        claims.sort_by(|a, b| a.id.cmp(&b.id));
        if node.kind == "topic" {
            attempts.sort_by(|a, b| {
                attempt_topic_score(b, node, &claim_ids, &graph)
                    .cmp(&attempt_topic_score(a, node, &claim_ids, &graph))
                    .then_with(|| a.id.cmp(&b.id))
            });
        } else {
            attempts.sort_by(|a, b| a.id.cmp(&b.id));
        }
        structures.sort_by(|a, b| a.id.cmp(&b.id));
        goals.sort_by(|a, b| a.id.cmp(&b.id));
        let mut category_counts: HashMap<String, usize> = HashMap::new();
        let mut state_counts: HashMap<String, usize> = HashMap::new();
        let mut result_samples: HashMap<String, Vec<Value>> = HashMap::new();
        let mut what_we_have = Vec::new();
        let mut what_we_seek = if node.kind == "goal" {
            goal_required_steps(node)
                .into_iter()
                .take(4)
                .map(|requirement| {
                    json!({"goal_id":node.id,"requirement":clip(&requirement,150),"source_ids":node_source_ids(node).into_iter().take(2).collect::<Vec<_>>(),"origin":"explicit goal required step"})
                })
                .collect()
        } else {
            Vec::new()
        };
        let mut hidden_input_risks = Vec::new();
        let mut source_ids = HashSet::new();
        let mut source_count_missing = 0usize;
        let mut source_ids_total = HashSet::new();
        let mut proof_gap_count = 0usize;
        let mut claims_without_enriched_gap_audit = 0usize;
        let mut claims_without_current_annotation_sources = 0usize;
        let mut claims_without_assessment = 0usize;
        let mut claims_with_risk_sample = HashSet::new();
        let mut static_candidate_count = 0usize;
        let mut examples_of_static_candidates = Vec::new();

        for claim in &claims {
            if claim
                .attrs
                .pointer("/cache_freshness/sources_current")
                .and_then(Value::as_bool)
                != Some(true)
            {
                claims_without_current_annotation_sources += 1;
            }
            let category = assessment_category(claim).unwrap_or("unknown").to_owned();
            *category_counts.entry(category.clone()).or_default() += 1;
            *state_counts
                .entry(
                    claim
                        .attrs
                        .get("enrichment_state")
                        .and_then(Value::as_str)
                        .unwrap_or("unknown")
                        .to_owned(),
                )
                .or_default() += 1;
            if assessment_category(claim).is_none() {
                claims_without_assessment += 1;
            }
            collect_node_source_ids(claim, &mut source_ids_total);
            if let Some(results) = claim
                .attrs
                .pointer("/semantic/results")
                .and_then(Value::as_array)
            {
                for result in results.iter().take(6) {
                    let refs = source_ref_ids(result);
                    for source in &refs {
                        source_ids.insert(source.clone());
                    }
                    let sample = json!({
                        "claim_id":claim.id,
                        "claim_label":clip(&claim.label, 110),
                        "category":category,
                        "assessment_scope":clip(assessment_scope(claim).unwrap_or(""), 130),
                        "result_id":result.get("id").cloned().unwrap_or(Value::Null),
                        "epistemic":result.get("epistemic").cloned().unwrap_or(Value::Null),
                        "annotation_freshness":claim.attrs.get("cache_freshness").cloned().unwrap_or(Value::Null),
                        "statement":clip(&value_text(result.get("statement")), 230),
                        "scope":clip(&value_text(result.get("scope")), 130),
                        "source_ids":refs.into_iter().take(3).collect::<Vec<_>>()
                    });
                    let samples = result_samples.entry(category.clone()).or_default();
                    if samples.len() < 10 {
                        samples.push(sample);
                    }
                }
            }
            let gaps = claim
                .attrs
                .pointer("/semantic/proof_gaps")
                .and_then(Value::as_array)
                .cloned()
                .unwrap_or_default();
            if claim.attrs.get("enrichment_state").and_then(Value::as_str) == Some("enriched") {
                proof_gap_count += gaps.len();
            } else {
                claims_without_enriched_gap_audit += 1;
            }
            // One gap per claim gives topic briefs wider claim coverage instead
            // of letting the first long dossiers consume the requirement list.
            for gap in gaps.iter().take(1) {
                let refs = source_ref_ids(gap);
                for source in &refs {
                    source_ids.insert(source.clone());
                }
                let requirement = json!({"claim_id":claim.id,"claim_label":clip(&claim.label, 85),"requirement":clip(&value_text(gap.get("missing_step").or_else(||gap.get("label"))), 170),"why_checks_insufficient":clip(&value_text(gap.get("why_checks_insufficient")), 100),"source_ids":refs.into_iter().take(2).collect::<Vec<_>>(),"origin":"semantic.proof_gaps","annotation_freshness":claim.attrs.get("cache_freshness").cloned().unwrap_or(Value::Null)});
                if what_we_seek.len() < 4 {
                    what_we_seek.push(requirement.clone());
                }
            }
            let text = claim
                .attrs
                .pointer("/semantic/overview_de")
                .and_then(Value::as_str)
                .or_else(|| (!claim.summary.is_empty()).then_some(claim.summary.as_str()))
                .unwrap_or("");
            if !text.is_empty() && what_we_have.len() < 2 {
                let statement_limit = if claim.id == node.id { 1_200 } else { 140 };
                what_we_have.push(json!({"claim_id":claim.id,"category":category,"scope":clip(assessment_scope(claim).unwrap_or(""),140),"statement":clip(text,statement_limit),"annotation_freshness":claim.attrs.get("cache_freshness").cloned().unwrap_or(Value::Null)}));
            }
            if let Some(inputs) = claim
                .attrs
                .pointer("/semantic/input_audit")
                .and_then(Value::as_array)
            {
                for input in inputs {
                    let origin = input
                        .get("origin")
                        .and_then(Value::as_str)
                        .unwrap_or("unknown");
                    let impact_risk = input
                        .get("impact")
                        .is_some_and(|impact| impact_matches_risk(impact));
                    let origin_risk = matches!(
                        origin,
                        "hardcoded" | "free_parameter" | "unknown" | "declared_axiom"
                    );
                    if (impact_risk || origin_risk)
                        && hidden_input_risks.len() < 3
                        && claims_with_risk_sample.insert(claim.id.clone())
                    {
                        let refs = source_ref_ids(input);
                        for source in &refs {
                            source_ids.insert(source.clone());
                        }
                        hidden_input_risks.push(json!({"claim_id":claim.id,"name":input.get("name").cloned().unwrap_or(Value::Null),"origin":origin,"impact":input.get("impact").cloned().unwrap_or(Value::Null),"value":clip(&value_text(input.get("value")),60),"used_for":clip(&value_text(input.get("used_for")),65),"source_ids":refs.into_iter().take(2).collect::<Vec<_>>(),"note":"Recorded input audit only; independently unverified."}));
                    }
                }
            }
        }
        // Static script candidates are sampled only as review signals, never
        // represented as proof gaps.
        let mut relevant_scripts = HashSet::new();
        for edge in &graph.edges {
            if edge.kind.contains("script") || edge.kind == "uses_script" || edge.kind == "checks" {
                if claim_ids.contains(&edge.src) {
                    relevant_scripts.insert(edge.dst.clone());
                }
                if claim_ids.contains(&edge.dst) {
                    relevant_scripts.insert(edge.src.clone());
                }
            }
        }
        let mut relevant_scripts: Vec<String> = relevant_scripts.into_iter().collect();
        relevant_scripts.sort();
        for script_id in relevant_scripts {
            if let Some(script) = by_id
                .get(script_id.as_str())
                .filter(|node| node.kind == "script")
            {
                if let Some(candidates) = script
                    .attrs
                    .pointer("/static_inputs/syntactic_input_candidates")
                    .and_then(Value::as_array)
                {
                    static_candidate_count += candidates.len();
                    for candidate in candidates.iter().take(2) {
                        examples_of_static_candidates.push(json!({"script_id":script.id,"candidate":clip(&value_text(candidate.get("name")),100),"line":candidate.get("line"),"notice":"Syntactic candidate only; not a proof gap."}));
                    }
                }
            }
        }
        let mut scoped_results = Vec::new();
        let categories = [
            "recorded_established",
            "conditional",
            "mixed",
            "open",
            "excluded",
            "assumed",
            "historical",
            "methodological_rule",
            "unknown",
        ];
        for category in categories {
            if let Some(samples) = result_samples.get_mut(category) {
                samples.sort_by(|a, b| {
                    a["claim_id"]
                        .as_str()
                        .cmp(&b["claim_id"].as_str())
                        .then_with(|| a["result_id"].as_str().cmp(&b["result_id"].as_str()))
                });
                let cap = if matches!(category, "mixed" | "open") {
                    2
                } else {
                    1
                };
                scoped_results.extend(samples.iter().take(cap).cloned());
            }
        }
        scoped_results.truncate(5);

        let attempt_nodes = attempts;
        let attempts_total = attempt_nodes.len();
        let typed_related_records = if matches!(node.kind.as_str(), "attempt" | "contract") {
            direct_typed_related_records(node, &graph, &by_id)
        } else {
            Vec::new()
        };
        let attempt_cards: Vec<Value> = attempt_nodes
            .iter()
            .skip(offset)
            .take(4)
            .map(|attempt| {
                collect_node_source_ids(attempt, &mut source_ids_total);
                let refs = node_source_ids(attempt);
                for source in &refs {
                    source_ids.insert(source.clone());
                }
                let mut card = attempt_card(attempt, &refs);
                remove_duplicate_findings(&mut card);
                card["source_refs"] = Value::Array(source_path_refs(&graph, refs.into_iter(), 4));
                card
            })
            .collect();
        let mut known_negatives = Vec::new();
        for claim in &claims {
            if assessment_category(claim) == Some("excluded") {
                let refs = node_source_ids(claim);
                for source in &refs {
                    source_ids.insert(source.clone());
                }
                if known_negatives.len() < 3 {
                    known_negatives.push(json!({"id":claim.id,"label":clip(&claim.label,100),"status":"excluded","scope":clip(assessment_scope(claim).unwrap_or(""),150),"source_ids":refs.into_iter().take(3).collect::<Vec<_>>(),"notice":"Recorded scoped exclusion; this does not generalize beyond its stated scope."}));
                }
            }
        }
        for attempt in &attempt_nodes {
            if let Some(verdict) = attempt_verdict(attempt) {
                if verdict_is_negative(&verdict) {
                    let refs = node_source_ids(attempt);
                    for source in &refs {
                        source_ids.insert(source.clone());
                    }
                    if known_negatives.len() < 3 {
                        known_negatives.push(json!({"id":attempt.id,"label":clip(&attempt.label,90),"verdict":verdict,"scope":clip(&attempt_scope(attempt),140),"source_ids":refs.into_iter().take(3).collect::<Vec<_>>(),"notice":"Recorded attempt verdict only; not proof of the linked claim."}));
                    }
                }
            }
        }
        let mut recorded_scoped_rule_outs = Vec::new();
        for attempt in &attempt_nodes {
            let rule_outs = attempt_rule_outs(attempt);
            if !rule_outs.is_empty() && recorded_scoped_rule_outs.len() < 1 {
                recorded_scoped_rule_outs.push(json!({
                    "id":attempt.id,
                    "label":clip(&attempt.label,120),
                    "verdict":attempt_verdict(attempt),
                    "scope":clip(&attempt_scope(attempt),140),
                    "rules_out":rule_outs.into_iter().take(2).map(|text|clip(&text,140)).collect::<Vec<_>>(),
                    "source_ids":node_source_ids(attempt).into_iter().take(3).collect::<Vec<_>>(),
                    "category":"documented_attempt_rule_out",
                    "notice":"Reported attempt conclusion within its recorded scope; not proof of a linked claim."
                }));
            }
        }
        for goal in &goals {
            if node.kind == "goal" && goal.id == node.id {
                continue;
            }
            collect_node_source_ids(goal, &mut source_ids_total);
            let refs = node_source_ids(goal);
            for source in &refs {
                source_ids.insert(source.clone());
            }
            let requirement = goal_requirement(goal);
            if !requirement.is_empty() {
                if what_we_seek.len() < 4 {
                    what_we_seek.push(json!({"goal_id":goal.id,"requirement":clip(&requirement,150),"source_ids":refs.into_iter().take(2).collect::<Vec<_>>(),"origin":"explicit linked goal"}));
                }
            }
        }
        for structure in &structures {
            let refs = node_source_ids(structure);
            for source in &refs {
                source_ids.insert(source.clone());
            }
        }
        let graph_source_ids: HashSet<&str> = graph
            .sources
            .iter()
            .map(|source| source.id.as_str())
            .collect();
        for source in &source_ids_total {
            if !graph_source_ids.contains(source.as_str()) {
                source_count_missing += 1;
            }
        }
        let mut all_sources: Vec<String> = source_ids.into_iter().collect();
        all_sources.sort();
        let source_ids_omitted = all_sources.len().saturating_sub(8);
        let sources: Vec<String> = all_sources.into_iter().take(8).collect();
        let theory = self
            .theory_overview()
            .await
            .unwrap_or_else(|_| json!({"available":false,"patterns":[]}));
        let topic_tokens = tokenize(&format!("{} {}", node.label, node.summary));
        let mut patterns = theory_patterns(&theory);
        let mut topic_related_ids = HashSet::from([node.id.clone()]);
        if node.kind == "topic" {
            let (member_ids, _, structure_ids, goal_ids) = topic_members(&graph, id, &by_id, false);
            topic_related_ids.extend(member_ids);
            topic_related_ids.extend(structure_ids);
            topic_related_ids.extend(goal_ids);
        }
        for edge in graph
            .edges
            .iter()
            .filter(|edge| edge.src == id || edge.dst == id)
        {
            topic_related_ids.insert(if edge.src == id {
                edge.dst.clone()
            } else {
                edge.src.clone()
            });
        }
        let mut topic_keywords = topic_tokens.clone();
        topic_keywords.extend(tokenize(&value_text(node.attrs.get("keywords"))));
        if let Some(orientation) = node.attrs.get("topic_orientation") {
            for key in [
                "role_de",
                "have_de",
                "seek_de",
                "connection_de",
                "next_bridge_de",
            ] {
                topic_keywords.extend(tokenize(&value_text(orientation.get(key))));
            }
        }
        topic_keywords.sort();
        topic_keywords.dedup();
        patterns.sort_by(|a, b| {
            let score = |pattern: &Value| {
                let related_score = pattern
                    .get("related_ids")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(Value::as_str)
                    .filter(|related| topic_related_ids.contains(*related))
                    .count()
                    * 20;
                let text = format!(
                    "{} {} {}",
                    value_text(pattern.get("label")),
                    value_text(pattern.get("summary")),
                    value_text(pattern.get("scope"))
                )
                .to_lowercase();
                related_score
                    + topic_keywords
                        .iter()
                        .filter(|token| text.contains(token.as_str()))
                        .count()
            };
            score(b)
                .cmp(&score(a))
                .then_with(|| a["id"].as_str().cmp(&b["id"].as_str()))
        });
        patterns.truncate(2);
        let patterns: Vec<Value> = patterns
            .iter()
            .map(|pattern| {
                json!({
                    "id":pattern.get("id"),
                    "type":pattern.get("type"),
                    "kind":pattern.get("kind"),
                    "label":clip(&value_text(pattern.get("label")),120),
                    "summary":clip(&value_text(pattern.get("summary")),200),
                    "source_ids":pattern.get("source_ids").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                    "related_ids":pattern.get("related_ids").and_then(Value::as_array).into_iter().flatten().take(4).cloned().collect::<Vec<_>>(),
                    "conditions":pattern.get("conditions").and_then(Value::as_array).into_iter().flatten().take(1).cloned().collect::<Vec<_>>(),
                    "counterevidence":pattern.get("counterevidence").and_then(Value::as_array).into_iter().flatten().take(1).cloned().collect::<Vec<_>>()
                })
            })
            .collect();
        let status = node
            .attrs
            .get("ledger")
            .and_then(|ledger| ledger.get("status"))
            .cloned()
            .unwrap_or(Value::Null);
        let proof_gap_value =
            if matches!(node.kind.as_str(), "attempt" | "contract") && claims.is_empty() {
                Value::Null
            } else if node.kind == "claim"
                && node.attrs.get("enrichment_state").and_then(Value::as_str) != Some("enriched")
            {
                Value::Null
            } else {
                json!(proof_gap_count)
            };
        let mut brief = json!({
            "id":node.id,
            "kind":node.kind,
            "enrichment_state":node.attrs.get("enrichment_state").cloned().unwrap_or(Value::Null),
        "publication_scope":node.attrs.get("publication_scope"),"content_scope":node.attrs.get("content_scope"),"summary_state":node.attrs.get("summary_state"),"source_current":node.attrs.get("source_current"),"source_freshness":node.attrs.get("source_freshness"),"reviewed":node.attrs.get("reviewed"),
            "question":clip(&node.label,220),
            "what_we_have":what_we_have,
            "topic_orientation":compact_topic_orientation(node.attrs.get("topic_orientation")),
            "what_we_seek":{"recorded_proof_gaps":proof_gap_value,"linked_goals":goals.len(),"goal_required_steps":if node.kind == "goal" { goal_required_steps(node).len() } else { 0 },"note":"Recorded gap counts cover enriched claims only; claims without enriched gap annotations are reported separately. See next_testable_requirements for sourced examples."},
            "counts":{"claims":claims.len(),"related_claims_total":related_claims_total,"claims_by_assessment":category_counts,"claims_by_enrichment_state":state_counts,"attempts":attempts_total,"structures":structures.len(),"goals":goals.len(),"recorded_proof_gaps":proof_gap_value,"claims_without_enriched_gap_audit":claims_without_enriched_gap_audit,"syntactic_input_candidates":static_candidate_count},
            "proof_gap_coverage":{"recorded_gaps_from_enriched_claims":proof_gap_value,"claims_without_enriched_gap_audit":claims_without_enriched_gap_audit,"claims_without_current_annotation_sources":claims_without_current_annotation_sources,"complete":if matches!(node.kind.as_str(), "attempt" | "contract") && claims.is_empty() { Value::Null } else { json!(claims_without_enriched_gap_audit == 0 && claims_without_current_annotation_sources == 0) },"scope":if matches!(node.kind.as_str(), "attempt" | "contract") && claims.is_empty() { "not_applicable_no_claim_membership" } else { "claim_annotations_in_scope" },"notice":"Historical extracted gaps describe the supplied annotation packet. Stale or unverified annotations are not a current research verdict; consult research_state and the full paper sources."},
            "membership_notice":"Claims marked related by topic_roles are excluded from topic briefs; claims without a topic role remain included.",
            "claim_status":status,
            "goal_status":if node.kind == "goal" { node.attrs.get("status").cloned().unwrap_or(Value::Null) } else { Value::Null },
            "scope":assessment_scope(node).map(|scope|clip(scope,300)),
            "scoped_results":scoped_results,
            "attempts":{"items":attempt_cards,"total":attempts_total,"offset":offset,"limit":5,"next_offset":offset+attempt_cards.len(),"truncated":offset+attempt_cards.len()<attempts_total,"notice":"Historical records describe their documented test scope; they do not prove a linked physical claim."},
            "attempt_record":if matches!(node.kind.as_str(), "attempt" | "contract") { let refs=node_source_ids(node); let mut card=attempt_card(node,&refs); remove_duplicate_findings(&mut card); card["source_refs"]=Value::Array(source_path_refs(&graph,refs.into_iter(),8)); card } else { Value::Null },
            "typed_related_records":typed_related_records,
            "known_negatives":known_negatives,
            "recorded_scoped_rule_outs":recorded_scoped_rule_outs,
            "hidden_input_risks":hidden_input_risks,
            "syntactic_input_sample":examples_of_static_candidates,
            "patterns":patterns,
            "source_refs":source_path_refs(&graph, sources.iter().cloned(), 8),
            "source_ids":sources,
            "source_ids_omitted":source_ids_omitted,
            "missing_coverage":{"claims_without_assessment":claims_without_assessment,"source_ids_missing_from_graph":source_count_missing,"source_refs_for_drilldown_only":true,"corpus_completeness":"not inferred from this brief"},
            "next_testable_requirements":what_we_seek,
            "notice":"Deterministic evidence summary. Recorded statuses, tests, source mentions, and extracted semantics are not upgraded into proof or physical conclusions.",
            "max_chars":max_chars
        });
        if let Some(record) = brief.get_mut("attempt_record") {
            if record.is_object() {
                let total = node_source_ids(node).len();
                let shown = record
                    .get("source_refs")
                    .and_then(Value::as_array)
                    .map(Vec::len)
                    .unwrap_or(0);
                record["source_refs_omitted"] = json!(total.saturating_sub(shown));
            }
        }
        compact_research_brief(&mut brief);
        Ok(fit_char_budget(brief, max_chars))
    }
    /// Build a small, deterministic research handoff from recorded graph data.
    /// The packet separates stored links and candidates from evidence; it never
    /// upgrades either into a derivation.
    pub async fn research_context(
        &self,
        id: &str,
        max_chars: usize,
        offset: usize,
    ) -> Result<Value> {
        let max_chars = bounded_char_limit(max_chars)?;
        let graph = self.snapshot().await?;
        let by_id: HashMap<&str, &Node> = graph
            .nodes
            .iter()
            .map(|node| (node.id.as_str(), node))
            .collect();
        let root = by_id
            .get(id)
            .copied()
            .ok_or_else(|| anyhow!("node not found: {id}"))?;
        if !matches!(
            root.kind.as_str(),
            "topic" | "claim" | "goal" | "attempt" | "contract"
        ) {
            bail!("research_context requires a topic, claim, goal, attempt, or contract id");
        }

        let (topic_claims, mut topic_attempt_ids, _topic_structures) = if root.kind == "topic" {
            let (claims, attempts, structures, _) = topic_members(&graph, id, &by_id, false);
            (claims, attempts, structures)
        } else if root.kind == "goal" {
            (
                goal_claim_ids(&graph, root, &by_id),
                HashSet::new(),
                HashSet::new(),
            )
        } else if matches!(root.kind.as_str(), "attempt" | "contract") {
            let mut claims = HashSet::new();
            for edge in &graph.edges {
                if edge.kind != "tests" {
                    continue;
                }
                let other = if edge.src == root.id {
                    Some(edge.dst.as_str())
                } else if edge.dst == root.id {
                    Some(edge.src.as_str())
                } else {
                    None
                };
                if let Some(claim_id) =
                    other.filter(|other| by_id.get(*other).is_some_and(|n| n.kind == "claim"))
                {
                    claims.insert(claim_id.to_owned());
                }
            }
            for key in ["claims_tested", "claim_ids"] {
                for container in [Some(&root.attrs), root.attrs.get("original_contract")] {
                    if let Some(values) = container
                        .and_then(|attrs| attrs.get(key))
                        .and_then(Value::as_array)
                    {
                        claims.extend(
                            values
                                .iter()
                                .filter_map(Value::as_str)
                                .filter_map(|value| normalize_claim_ref(value, &by_id)),
                        );
                    }
                }
            }
            let (structures, _goals) = linked_structures_and_goals_typed(&graph, id, &by_id);
            (claims, HashSet::from([root.id.clone()]), structures)
        } else {
            (
                HashSet::from([root.id.clone()]),
                HashSet::new(),
                HashSet::new(),
            )
        };
        // Keep primary/related claim membership scoped as requested, but do
        // not let that filter hide attempts which have a stored topic edge.
        // Some imported `about` edges are explicitly marked as inherited from
        // claim-topic membership; those remain visible with that provenance.
        let topic_attempt_link_basis = if root.kind == "topic" {
            topic_linked_attempts(&graph, id, &by_id)
        } else {
            HashMap::new()
        };
        topic_attempt_ids.extend(topic_attempt_link_basis.keys().cloned());
        if root.kind == "goal" {
            for edge in &graph.edges {
                if edge.kind != "leaves_open" || (edge.src != root.id && edge.dst != root.id) {
                    continue;
                }
                let attempt_id = if edge.src == root.id {
                    &edge.dst
                } else {
                    &edge.src
                };
                if by_id
                    .get(attempt_id.as_str())
                    .is_some_and(|node| node.kind == "attempt")
                {
                    topic_attempt_ids.insert(attempt_id.clone());
                }
            }
        }
        let mut claim_ids: Vec<String> = topic_claims.into_iter().collect();
        claim_ids.sort();
        let total_claims = claim_ids.len();
        let (related_claims_total, related_claims): (usize, Vec<&Node>) = if root.kind == "topic" {
            let (all_claims, _, _, _) = topic_members(&graph, id, &by_id, true);
            let mut related: Vec<&Node> = all_claims
                .iter()
                .filter_map(|claim_id| by_id.get(claim_id.as_str()).copied())
                .filter(|claim| claim_topic_role(claim, id) == Some("related"))
                .collect();
            related.sort_by(|a, b| {
                topic_related_claim_score(b, root)
                    .cmp(&topic_related_claim_score(a, root))
                    .then_with(|| a.id.cmp(&b.id))
            });
            let total = related.len();
            (total, related.into_iter().take(3).collect())
        } else {
            (0, Vec::new())
        };
        let all_claim_set: HashSet<String> = claim_ids.iter().cloned().collect();
        // Keep the cross-claim sample bounded; the full population remains
        // visible in coverage counts and can be paged with offset.
        let page_ids: Vec<String> = claim_ids.iter().skip(offset).take(40).cloned().collect();
        let page_set: HashSet<String> = page_ids.iter().cloned().collect();
        let claims: Vec<&Node> = page_ids
            .iter()
            .filter_map(|claim_id| by_id.get(claim_id.as_str()).copied())
            .collect();

        let mut seed_members: HashMap<String, Vec<String>> = HashMap::new();
        let mut object_members: HashMap<String, Vec<(String, String)>> = HashMap::new();
        let mut object_records = Vec::new();
        let mut assumption_records = Vec::new();
        let mut assumption_records_total = 0usize;
        let mut enrichment_counts: HashMap<String, usize> = HashMap::new();
        let mut status_counts: HashMap<String, usize> = HashMap::new();
        let mut shared_sources = HashSet::new();
        shared_sources.extend(node_source_ids(root));
        for claim in &related_claims {
            shared_sources.extend(node_source_ids(claim));
        }
        let mut sampled_gaps = 0usize;
        let mut premise_records = Vec::new();
        let mut premise_records_total = 0usize;
        for claim in &claims {
            let enrichment = claim
                .attrs
                .get("enrichment_state")
                .and_then(Value::as_str)
                .unwrap_or("unknown");
            *enrichment_counts.entry(enrichment.to_owned()).or_default() += 1;
            let status = assessment_category(claim).unwrap_or("unknown");
            *status_counts.entry(status.to_owned()).or_default() += 1;
            shared_sources.extend(node_source_ids(claim));
            let semantic = claim.attrs.get("semantic").unwrap_or(&Value::Null);
            if let Some(gaps) = semantic.get("proof_gaps").and_then(Value::as_array) {
                sampled_gaps += gaps.len();
            }
            if let Some(premises) = semantic.get("premises").and_then(Value::as_array) {
                premise_records_total += premises.len();
                for premise in premises
                    .iter()
                    .take(6usize.saturating_sub(premise_records.len()))
                {
                    let text = premise
                        .get("text")
                        .or_else(|| premise.get("statement"))
                        .or_else(|| premise.get("label"));
                    premise_records.push(json!({
                        "claim_id":claim.id,
                        "id":premise.get("id"),
                        "kind":premise.get("kind"),
                        "text":clip(&value_text(text),180),
                        "source_ids":source_ref_ids(premise).into_iter().take(3).collect::<Vec<_>>(),
                        "notice":"Stored premise annotation; source support and derivational sufficiency require inspection."
                    }));
                }
            }
            for key in ["seed_ids", "seeds"] {
                if let Some(seeds) = semantic.get(key).and_then(Value::as_array) {
                    for seed in seeds {
                        let seed_id = seed
                            .as_str()
                            .or_else(|| seed.get("id").and_then(Value::as_str));
                        if let Some(seed_id) = seed_id.filter(|value| !value.is_empty()) {
                            seed_members
                                .entry(seed_id.to_owned())
                                .or_default()
                                .push(claim.id.clone());
                        }
                    }
                }
            }
            if let Some(objects) = semantic.get("objects").and_then(Value::as_array) {
                for object in objects.iter().take(12) {
                    let identity_scope = object
                        .get("identity_scope")
                        .and_then(Value::as_str)
                        .unwrap_or("unspecified");
                    let object_id = object
                        .get("id")
                        .and_then(Value::as_str)
                        .filter(|value| !value.is_empty());
                    if object_records.len() < 4 {
                        object_records.push(json!({
                            "claim_id":claim.id,
                            "id":object_id,
                            "label":clip(&value_text(object.get("name").or_else(|| object.get("label"))),100),
                            "kind":object.get("kind"),
                            "identity_scope":identity_scope,
                            "source_ids":source_ref_ids(object)
                        }));
                    }
                    // A shared object is recorded only where extraction supplies
                    // a stable ID and explicitly says it is not claim-local.
                    if identity_scope != "claim_local" {
                        if let Some(object_id) = object_id {
                            object_members
                                .entry(object_id.to_owned())
                                .or_default()
                                .push((claim.id.clone(), identity_scope.to_owned()));
                        }
                    }
                }
            }
            if let Some(inputs) = semantic.get("input_audit").and_then(Value::as_array) {
                let mut ranked_inputs: Vec<&Value> = inputs
                    .iter()
                    .filter(|input| input_is_research_relevant(input))
                    .collect();
                ranked_inputs.sort_by_key(|input| input_importance(input));
                for input in ranked_inputs {
                    assumption_records_total += 1;
                    let origin = input
                        .get("origin")
                        .and_then(Value::as_str)
                        .unwrap_or("unknown");
                    let impact = input.get("impact").cloned().unwrap_or(Value::Null);
                    assumption_records.push(json!({
                        "claim_id":claim.id,
                        "name":clip(&value_text(input.get("name")),80),
                        "origin":origin,
                        "impact":impact,
                        "used_for":clip(&value_text(input.get("used_for")),100),
                        "source_ids":source_ref_ids(input).into_iter().take(2).collect::<Vec<_>>(),
                        "priority":input_importance(input).0
                    }));
                }
            }
        }
        let recorded_shared_seeds: Vec<Value> = seed_members
            .into_iter()
            .filter(|(_, members)| {
                let unique: HashSet<&str> = members.iter().map(String::as_str).collect();
                unique.len() > 1
            })
            .map(|(seed_id, members)| {
                let mut members = members;
                members.sort();
                members.dedup();
                json!({"seed_id":seed_id,"claim_ids":members,"relation":"same_recorded_seed_id"})
            })
            .collect();
        let mut recorded_shared_seeds = recorded_shared_seeds;
        recorded_shared_seeds.sort_by(|a, b| a["seed_id"].as_str().cmp(&b["seed_id"].as_str()));
        let mut recorded_shared_objects: Vec<Value> = object_members
            .into_iter()
            .filter_map(|(object_id, members)| {
                let unique: HashSet<&str> = members.iter().map(|(claim_id, _)| claim_id.as_str()).collect();
                if unique.len() < 2 { return None; }
                let mut claims: Vec<String> = unique.into_iter().map(str::to_owned).collect();
                claims.sort();
                let mut scopes: Vec<String> = members.into_iter().map(|(_,scope)|scope).collect();
                scopes.sort(); scopes.dedup();
                Some(json!({"object_id":object_id,"claim_ids":claims,"identity_scopes":scopes,"relation":"same_explicit_object_id","notice":"Only explicit non-claim-local IDs are grouped; matching labels or symbols are not identity evidence."}))
            })
            .collect();
        recorded_shared_objects
            .sort_by(|a, b| a["object_id"].as_str().cmp(&b["object_id"].as_str()));
        assumption_records.sort_by(|a, b| {
            a["priority"]
                .as_u64()
                .cmp(&b["priority"].as_u64())
                .then_with(|| a["claim_id"].as_str().cmp(&b["claim_id"].as_str()))
                .then_with(|| a["name"].as_str().cmp(&b["name"].as_str()))
        });
        assumption_records.truncate(4);
        for record in &mut assumption_records {
            if let Some(attrs) = record.as_object_mut() {
                attrs.remove("priority");
            }
        }
        let mut structure_members = HashMap::<String, Vec<String>>::new();
        for edge in &graph.edges {
            if !matches!(
                edge.kind.as_str(),
                "uses_structure"
                    | "about_structure"
                    | "implements"
                    | "related_structure"
                    | "documents_structure"
            ) {
                continue;
            }
            let link = if by_id
                .get(edge.src.as_str())
                .is_some_and(|node| node.kind == "claim")
                && by_id
                    .get(edge.dst.as_str())
                    .is_some_and(|node| node.kind == "structure")
            {
                Some((edge.src.as_str(), edge.dst.as_str()))
            } else if by_id
                .get(edge.dst.as_str())
                .is_some_and(|node| node.kind == "claim")
                && by_id
                    .get(edge.src.as_str())
                    .is_some_and(|node| node.kind == "structure")
            {
                Some((edge.dst.as_str(), edge.src.as_str()))
            } else {
                None
            };
            if let Some((claim_id, structure_id)) =
                link.filter(|(claim_id, _)| page_set.contains(*claim_id))
            {
                structure_members
                    .entry(structure_id.to_owned())
                    .or_default()
                    .push(claim_id.to_owned());
            }
        }
        let mut recorded_shared_structures: Vec<Value> = structure_members
            .into_iter()
            .filter_map(|(structure_id, mut members)| {
                members.sort(); members.dedup();
                if members.len() < 2 { return None; }
                let structure = by_id.get(structure_id.as_str())?;
                Some(json!({"structure_id":structure_id,"label":clip(&structure.label,100),"claim_count":members.len(),"claim_ids":members.iter().take(8).cloned().collect::<Vec<_>>(),"claim_ids_omitted":members.len().saturating_sub(8),"relation":"same_explicit_structure_node","source_ids":node_source_ids(structure).into_iter().take(2).collect::<Vec<_>>(),"notice":"Explicit structure links are a recorded correlation, not evidence that linked claims are proven or physically identified."}))
            })
            .collect();
        recorded_shared_structures
            .sort_by(|a, b| a["structure_id"].as_str().cmp(&b["structure_id"].as_str()));

        let mut recorded_connections = Vec::new();
        let mut relevant_ids: HashSet<String> =
            if matches!(root.kind.as_str(), "attempt" | "contract") {
                HashSet::new()
            } else {
                page_ids.iter().cloned().collect()
            };
        relevant_ids.insert(root.id.clone());
        relevant_ids.extend(related_claims.iter().map(|claim| claim.id.clone()));
        for edge in &graph.edges {
            let touches = relevant_ids.contains(&edge.src) || relevant_ids.contains(&edge.dst);
            if !touches {
                continue;
            }
            let (Some(src), Some(dst)) =
                (by_id.get(edge.src.as_str()), by_id.get(edge.dst.as_str()))
            else {
                continue;
            };
            if !matches!(
                src.kind.as_str(),
                "claim"
                    | "topic"
                    | "structure"
                    | "goal"
                    | "attempt"
                    | "script"
                    | "derivation"
                    | "result"
                    | "object"
                    | "assumption"
            ) && !matches!(
                dst.kind.as_str(),
                "claim"
                    | "topic"
                    | "structure"
                    | "goal"
                    | "attempt"
                    | "script"
                    | "derivation"
                    | "result"
                    | "object"
                    | "assumption"
            ) {
                continue;
            }
            let mut connection = json!({"src":edge.src,"dst":edge.dst,"edge_type":edge.kind,"src_type":src.kind,"dst_type":dst.kind});
            if let Some(reason) = edge.attrs.get("reason").and_then(Value::as_str) {
                connection["reason"] = json!(clip(reason, 260));
            }
            if let Some(conditions) = edge.attrs.get("conditions").and_then(Value::as_array) {
                connection["conditions"] = Value::Array(
                    conditions
                        .iter()
                        .take(2)
                        .map(|condition| json!(clip(&value_text(Some(condition)), 120)))
                        .collect(),
                );
            }
            if let Some(source_refs) = edge.attrs.get("source_refs").and_then(Value::as_array) {
                connection["source_ids"] = Value::Array(
                    source_refs
                        .iter()
                        .filter_map(Value::as_str)
                        .take(3)
                        .map(|source| json!(source))
                        .collect(),
                );
            }
            recorded_connections.push(connection);
        }
        let edge_priority = |edge: &Value| match edge["edge_type"].as_str().unwrap_or("") {
            "requires"
            | "required_together"
            | "derives"
            | "has_derivation"
            | "has_result"
            | "uses_claim_under_scope"
            | "implements"
            | "checks" => 0,
            // Direct claim mentions may carry a scoped later result. Keep the
            // recorded relation explicit, but ahead of broad seed correlations.
            "mentions" => 1,
            "uses_input" | "originates_from" | "uses_result" | "describes_object" => 2,
            "related_structure" | "documents_structure" => 3,
            "shares_seed" | "shared_seed" => 4,
            "references" | "about" => 5,
            _ => 4,
        };
        recorded_connections.sort_by(|a, b| {
            edge_priority(a)
                .cmp(&edge_priority(b))
                .then_with(|| a["src"].as_str().cmp(&b["src"].as_str()))
                .then_with(|| a["edge_type"].as_str().cmp(&b["edge_type"].as_str()))
                .then_with(|| a["dst"].as_str().cmp(&b["dst"].as_str()))
        });
        let total_connections = recorded_connections.len();
        recorded_connections.truncate(12);

        let mut attempt_ids = std::mem::take(&mut topic_attempt_ids);
        if matches!(root.kind.as_str(), "attempt" | "contract") {
            attempt_ids.insert(root.id.clone());
        }
        for edge in &graph.edges {
            if matches!(root.kind.as_str(), "attempt" | "contract") {
                break;
            }
            if edge.kind == "tests" && all_claim_set.contains(&edge.dst) {
                if by_id
                    .get(edge.src.as_str())
                    .is_some_and(|node| node.kind == "attempt")
                {
                    attempt_ids.insert(edge.src.clone());
                }
            }
            if edge.kind == "tests" && all_claim_set.contains(&edge.src) {
                if by_id
                    .get(edge.dst.as_str())
                    .is_some_and(|node| node.kind == "attempt")
                {
                    attempt_ids.insert(edge.dst.clone());
                }
            }
        }
        let mut attempts: Vec<&Node> = attempt_ids
            .iter()
            .filter_map(|attempt_id| by_id.get(attempt_id.as_str()).copied())
            .collect();
        attempts.sort_by(|a, b| {
            attempt_topic_score(b, root, &all_claim_set, &graph)
                .cmp(&attempt_topic_score(a, root, &all_claim_set, &graph))
                .then_with(|| a.id.cmp(&b.id))
        });
        let total_attempts = attempts.len();
        let prior_attempts: Vec<Value> = attempts
            .iter()
            .skip(offset)
            .take(3)
            .map(|attempt| {
                let mut card = attempt_card(attempt, &node_source_ids(attempt));
                remove_duplicate_findings(&mut card);
                let refs = node_source_ids(attempt);
                if matches!(root.kind.as_str(), "attempt" | "contract") {
                    card["source_refs_omitted"] = json!(refs.len());
                    card.as_object_mut()
                        .map(|object| object.remove("source_refs"));
                } else {
                    card["source_refs"] = Value::Array(source_path_refs(&graph, refs.clone(), 4));
                    card["source_refs_omitted"] = json!(refs
                        .len()
                        .saturating_sub(card["source_refs"].as_array().map(Vec::len).unwrap_or(0)));
                }
                let basis = if matches!(root.kind.as_str(), "attempt" | "contract") {
                    "root_attempt_id"
                } else {
                    topic_attempt_link_basis
                        .get(&attempt.id)
                        .map(String::as_str)
                        .unwrap_or("tests_primary_topic_claim")
                };
                card["topic_link_basis"] = json!(basis);
                card
            })
            .collect();
        let typed_related_records = if matches!(root.kind.as_str(), "attempt" | "contract") {
            direct_typed_related_records(root, &graph, &by_id)
        } else {
            Vec::new()
        };
        let mut scoped_negatives = Vec::new();
        for attempt in &attempts {
            let Some(verdict) = attempt_verdict(attempt) else {
                continue;
            };
            let rules_out = attempt_rule_outs(attempt);
            if !verdict_is_negative(&verdict) && rules_out.is_empty() {
                continue;
            }
            scoped_negatives.push(json!({"id":attempt.id,"verdict":verdict,"scope":clip(&attempt_scope(attempt),140),"rules_out":rules_out.into_iter().take(2).map(|value|clip(&value,120)).collect::<Vec<_>>(),"source_ids":node_source_ids(attempt).into_iter().take(3).collect::<Vec<_>>(),"notice":"Recorded within this attempt's stated scope; not a general exclusion or proof."}));
            if scoped_negatives.len() >= 3 {
                break;
            }
        }
        let mut pattern_relevant_ids: HashSet<String> = claim_ids.iter().cloned().collect();
        pattern_relevant_ids.insert(root.id.clone());
        let mut candidates: Vec<Value> = theory_patterns(&self.read_theory_map().await.unwrap_or_else(|_| json!({"patterns":[]})))
            .into_iter()
            .filter_map(|pattern| {
                let related: Vec<String> = pattern.get("related_ids").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str).map(str::to_owned).collect();
                let overlap = related.iter().filter(|related_id| pattern_relevant_ids.contains(*related_id)).count();
                (overlap > 0).then(|| json!({
                    "id":pattern.get("id"), "kind":pattern.get("kind").cloned().unwrap_or_else(||json!("authored_pattern")),
                    "direct_root_link":related.iter().any(|related_id| related_id == id),
                    "label":clip(&value_text(pattern.get("label_de").or_else(||pattern.get("label"))),100),
                    "label_en":clip(&value_text(pattern.get("label_en")),100),
                    "summary":clip(&value_text(pattern.get("summary_de").or_else(||pattern.get("summary"))),140),
                    "summary_en":clip(&value_text(pattern.get("summary_en")),140),
                    "reviewed":pattern.pointer("/provenance/reviewed").or_else(||pattern.get("reviewed")).unwrap_or(&Value::Null),
                    "conditions":pattern.get("conditions").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                    "counterevidence":pattern.get("counterevidence").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                    "source_ids":source_ref_ids(&pattern).into_iter().take(3).collect::<Vec<_>>(),
                    "related_ids":related.into_iter().take(8).collect::<Vec<_>>(),
                    "interpretation":"candidate_only", "overlap_count":overlap
                }))
            }).collect();
        candidates.sort_by(|a, b| {
            b["direct_root_link"]
                .as_bool()
                .cmp(&a["direct_root_link"].as_bool())
                .then_with(|| {
                    b["overlap_count"]
                        .as_u64()
                        .cmp(&a["overlap_count"].as_u64())
                })
                .then_with(|| a["id"].as_str().cmp(&b["id"].as_str()))
        });
        candidates.truncate(4);
        let mut explicit_seed_links: Vec<Value> = graph.edges.iter()
            .filter(|edge| matches!(edge.kind.as_str(), "shares_seed" | "shared_seed") && relevant_ids.contains(&edge.src) && relevant_ids.contains(&edge.dst))
            .map(|edge| json!({"src":edge.src,"dst":edge.dst,"edge_type":edge.kind,"relation":"explicit_recorded_seed_link","notice":"Stored relation only; not evidence of a proof or physical identity."}))
            .collect();
        explicit_seed_links.sort_by(|a, b| {
            a["src"]
                .as_str()
                .cmp(&b["src"].as_str())
                .then_with(|| a["dst"].as_str().cmp(&b["dst"].as_str()))
        });
        explicit_seed_links.truncate(8);
        let total_source_ids = shared_sources.len();
        let mut sources = Vec::new();
        if matches!(root.kind.as_str(), "attempt" | "contract") {
            let mut own_sources = node_source_ids(root);
            own_sources.sort();
            sources.extend(own_sources.into_iter().take(5));
        }
        let mut other_sources: Vec<String> = shared_sources.into_iter().collect();
        other_sources.sort();
        for source in other_sources {
            if sources.len()
                >= if matches!(root.kind.as_str(), "attempt" | "contract") {
                    5
                } else {
                    8
                }
            {
                break;
            }
            if !sources.contains(&source) {
                sources.push(source);
            }
        }
        let brief = self.research_brief(id, max_chars, offset).await?;
        let next_discriminators = brief
            .get("next_testable_requirements")
            .cloned()
            .unwrap_or_else(|| json!([]));
        let claim_status = compact_claim_status(
            brief
                .get("claim_status")
                .or_else(|| root.attrs.get("assessment"))
                .unwrap_or(&Value::Null),
        );
        let all_integrity_issues = graph
            .meta
            .pointer("/graph_integrity/issues")
            .and_then(Value::as_array);
        let relevant_integrity_issues: Vec<Value> = all_integrity_issues
            .into_iter()
            .flatten()
            .filter(|issue| {
                issue
                    .get("src")
                    .and_then(Value::as_str)
                    .is_some_and(|value| relevant_ids.contains(value))
                    || issue
                        .get("dst")
                        .and_then(Value::as_str)
                        .is_some_and(|value| relevant_ids.contains(value))
            })
            .take(5)
            .cloned()
            .collect();
        let packet = json!({
            "id":root.id,
            "kind":root.kind,
            "question":clip(&root.label,180),
            "topic_orientation":compact_topic_orientation(root.attrs.get("topic_orientation")),
            "claim_status":claim_status,
            "goal_status":if root.kind == "goal" { root.attrs.get("status").cloned().unwrap_or(Value::Null) } else { Value::Null },
            "scope":if matches!(root.kind.as_str(), "attempt" | "contract") { let scope=attempt_scope(root); (!scope.is_empty()).then(||clip(&scope,240)) } else { assessment_scope(root).map(|scope|clip(scope,240)) },
            "scope_recorded":if matches!(root.kind.as_str(), "attempt" | "contract") { !attempt_scope(root).is_empty() } else { assessment_scope(root).is_some() },
            "attempt_record":if matches!(root.kind.as_str(), "attempt" | "contract") { brief.get("attempt_record").cloned().unwrap_or(Value::Null) } else { Value::Null },
            "coverage":{"claims_total":total_claims,"claims_in_page":claims.len(),"claims_offset":offset,"claims_next_offset":offset+claims.len(),"claims_truncated":offset+claims.len()<total_claims,"related_claims_total":related_claims_total,"related_claims_in_packet":related_claims.len(),"related_claims_omitted":related_claims_total.saturating_sub(related_claims.len()),"enrichment_states_in_page":enrichment_counts,"assessment_categories_in_page":status_counts,"proof_gaps_recorded_in_page":sampled_gaps,"premise_records_total":premise_records_total,"premise_records_in_packet":premise_records.len(),"premise_records_omitted":premise_records_total.saturating_sub(premise_records.len()),"assumption_records_total":assumption_records_total,"assumption_records_in_packet":assumption_records.len(),"assumption_records_omitted":assumption_records_total.saturating_sub(assumption_records.len()),"recorded_connection_total":total_connections,"recorded_attempt_total":total_attempts,"source_ids_total":total_source_ids,"unresolved_semantic_links":graph.meta.pointer("/graph_integrity/unresolved_semantic_links").cloned().unwrap_or(Value::Null),"graph_integrity_issue_total":all_integrity_issues.map(Vec::len).unwrap_or(0),"corpus_completeness":"not inferred"},
            "relevant_graph_integrity_issues":relevant_integrity_issues,
            "what_we_have":brief.get("what_we_have").cloned().unwrap_or_else(||json!([])),
            "premise_records":premise_records,
            "related_claims":{"items":related_claims.iter().map(|claim|json!({"id":claim.id,"label":clip(&claim.label,120),"category":assessment_category(claim).unwrap_or("unknown"),"scope":assessment_scope(claim).map(|scope|clip(scope,140)),"summary":clip(claim.attrs.pointer("/semantic/overview_de").and_then(Value::as_str).unwrap_or(&claim.summary),180),"topic_role":"related","source_ids":node_source_ids(claim).into_iter().take(2).collect::<Vec<_>>()})).collect::<Vec<_>>(),"total":related_claims_total,"offset":0,"limit":3,"truncated":related_claims_total>related_claims.len(),"notice":"These are representative claims explicitly marked related to the topic. They are shown separately and excluded from primary claim counts and status totals."},
            "recorded_shared_inputs":{"seeds":recorded_shared_seeds,"seed_links":explicit_seed_links,"objects":recorded_shared_objects,"structures":recorded_shared_structures,"identity_notice":"Only repeated explicit seed IDs or non-claim-local object IDs are grouped. Equal labels, symbols, or numeric values do not establish a shared object."},
            "recorded_connections":recorded_connections,
            "typed_related_records":typed_related_records,
            "prior_attempts":{"items":prior_attempts,"total":total_attempts,"offset":offset,"next_offset":offset+3.min(total_attempts.saturating_sub(offset)),"truncated":offset+3<total_attempts,"notice":"Recorded attempts are historical observations and do not prove linked claims. Topic edges inherited from claim-topic membership are labeled as inherited links."},
            "scoped_negatives":scoped_negatives,
            "candidate_connections":candidates,
            "object_records":object_records,
            "assumption_records":assumption_records,
            "next_discriminators":next_discriminators,
            "source_refs":source_path_refs(&graph, sources.iter().cloned(), if matches!(root.kind.as_str(), "attempt" | "contract") { 5 } else { 8 }),
            "source_ids":sources,
            "source_ids_omitted":total_source_ids.saturating_sub(8),
            "notice":"Research orientation only. Stored edges, common inputs, annotations, and authored patterns are surfaced with their recorded scope; no new proof or physical interpretation is inferred.",
            "max_chars":max_chars,
            "truncated":false
        });
        let program = crate::research::load(&self.research_path())?;
        let state =
            crate::research::state(&program, &graph, &self.repo_root, id, "de", 4000, 0, 0)?;
        let guide = json!({"revision":state["revision"],"coverage":state["coverage"],"questions":state["problem_index"],"next_tool":"research_state","query":id,"notice":"Read the durable question state before proposing another experiment. This index preserves resolved subproblems and scoped exclusions independently of historical PARTIAL verdicts."});
        let guide_chars = serde_json::to_string(&guide)?.chars().count() + 40;
        if max_chars > guide_chars + 1000 {
            let mut result = fit_research_context_budget(packet, max_chars - guide_chars);
            result["research_state"] = guide;
            Ok(result)
        } else {
            Ok(fit_research_context_budget(packet, max_chars))
        }
    }

    pub async fn compare_aspects(&self, ids: &[String], max_chars: usize) -> Result<Value> {
        let max_chars = bounded_char_limit(max_chars)?;
        if !(2..=6).contains(&ids.len()) {
            bail!("compare_aspects needs 2 to 6 node IDs");
        }
        let graph = self.snapshot().await?;
        let by_id: HashMap<&str, &Node> = graph
            .nodes
            .iter()
            .map(|node| (node.id.as_str(), node))
            .collect();
        let mut items = Vec::new();
        let mut upstream_sets = Vec::new();
        let mut seed_sets = Vec::new();
        let mut source_sets = Vec::new();
        let mut graph_link_sets = Vec::new();
        let mut graph_link_signature_sets = Vec::new();
        let mut shared_graph_link_records = HashMap::<String, Value>::new();
        for id in ids {
            let node = by_id
                .get(id.as_str())
                .copied()
                .ok_or_else(|| anyhow!("node not found: {id}"))?;
            let semantic = node
                .attrs
                .get("semantic")
                .cloned()
                .unwrap_or_else(|| json!({}));
            let mut upstream = HashSet::new();
            if let Some(derivations) = semantic.get("derivations").and_then(Value::as_array) {
                for derivation in derivations {
                    for value in derivation
                        .get("upstream_claim_ids")
                        .and_then(Value::as_array)
                        .into_iter()
                        .flatten()
                        .filter_map(Value::as_str)
                    {
                        upstream.insert(value.to_owned());
                    }
                }
            }
            if let Some(inputs) = semantic.get("input_audit").and_then(Value::as_array) {
                for input in inputs {
                    if let Some(value) = input
                        .get("upstream_claim_id")
                        .and_then(Value::as_str)
                        .filter(|value| !value.is_empty())
                    {
                        upstream.insert(value.to_owned());
                    }
                }
            }
            let mut seeds = HashSet::new();
            for key in ["seed_ids", "seeds"] {
                if let Some(values) = semantic.get(key).and_then(Value::as_array) {
                    seeds.extend(values.iter().filter_map(Value::as_str).map(str::to_owned));
                }
            }
            if let Some(value) = node.attrs.get("seed").and_then(Value::as_str) {
                seeds.insert(value.to_owned());
            }
            let mut sources: HashSet<String> = node_source_ids(node).into_iter().collect();
            let mut graph_links = HashMap::<String, Value>::new();
            let mut graph_link_ids = HashSet::<String>::new();
            let mut graph_link_signatures = HashSet::<String>::new();
            for edge in graph
                .edges
                .iter()
                .filter(|edge| edge.src == *id || edge.dst == *id)
            {
                let (direction, other_id) = if edge.src == *id {
                    ("outbound", edge.dst.as_str())
                } else {
                    ("inbound", edge.src.as_str())
                };
                let Some(other) = by_id.get(other_id).copied() else {
                    continue;
                };
                let key = format!("{}:{}:{}", direction, edge.kind, other.id);
                graph_link_ids.insert(other.id.clone());
                let signature = format!("{}:{}:{}", direction, edge.kind, other.id);
                graph_link_signatures.insert(signature.clone());
                shared_graph_link_records
                    .entry(signature)
                    .or_insert_with(|| {
                        json!({
                            "direction":direction,
                            "edge_type":edge.kind,
                            "node_id":other.id,
                            "node_type":other.kind,
                            "label":clip(&other.label,120)
                        })
                    });
                graph_links.entry(key).or_insert_with(|| {
                    json!({
                        "direction":direction,
                        "edge_type":edge.kind,
                        "node_id":other.id,
                        "node_type":other.kind,
                        "label":clip(&other.label,120)
                    })
                });
            }
            let mut topic_contents = Value::Null;
            if node.kind == "topic" {
                let (member_ids, _, _, _) = topic_members(&graph, id, &by_id, false);
                let mut members: Vec<&Node> = member_ids
                    .iter()
                    .filter_map(|member_id| by_id.get(member_id.as_str()).copied())
                    .collect();
                members.sort_by(|a, b| a.id.cmp(&b.id));
                let mut category_counts = HashMap::<String, usize>::new();
                for member in &members {
                    let category = assessment_category(member).unwrap_or("unknown");
                    *category_counts.entry(category.to_owned()).or_default() += 1;
                    let member_semantic = member.attrs.get("semantic").unwrap_or(&Value::Null);
                    if let Some(derivations) =
                        member_semantic.get("derivations").and_then(Value::as_array)
                    {
                        for derivation in derivations {
                            for upstream_id in derivation
                                .get("upstream_claim_ids")
                                .and_then(Value::as_array)
                                .into_iter()
                                .flatten()
                                .filter_map(Value::as_str)
                            {
                                upstream.insert(upstream_id.to_owned());
                            }
                        }
                    }
                    for key in ["seed_ids", "seeds"] {
                        if let Some(values) = member_semantic.get(key).and_then(Value::as_array) {
                            seeds
                                .extend(values.iter().filter_map(Value::as_str).map(str::to_owned));
                        }
                    }
                    sources.extend(node_source_ids(member));
                }
                let sample: Vec<Value> = members
                    .iter()
                    .take(5)
                    .map(|member| {
                        json!({
                            "id":member.id,
                            "label":clip(&member.label,120),
                            "assessment":assessment_category(member).unwrap_or("unknown"),
                            "summary":clip(&member.summary,180)
                        })
                    })
                    .collect();
                topic_contents = json!({"claim_count":members.len(),"assessment_counts":category_counts,"claim_sample":sample,"sample_notice":"Representative member records only; use topic_overview for paged membership."});
            }
            let mut upstream_values: Vec<String> = upstream.iter().cloned().collect();
            upstream_values.sort();
            let mut seed_values: Vec<String> = seeds.iter().cloned().collect();
            seed_values.sort();
            let mut source_values: Vec<String> = sources.iter().cloned().collect();
            source_values.sort();
            let mut links: Vec<Value> = graph_links.into_values().collect();
            links.sort_by(|a, b| {
                a["node_id"]
                    .as_str()
                    .cmp(&b["node_id"].as_str())
                    .then_with(|| a["edge_type"].as_str().cmp(&b["edge_type"].as_str()))
                    .then_with(|| a["direction"].as_str().cmp(&b["direction"].as_str()))
            });
            links.truncate(20);
            let is_record = matches!(node.kind.as_str(), "attempt" | "contract");
            let mut ranked_inputs: Vec<&Value> = semantic
                .get("input_audit")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .filter(|input| input_is_research_relevant(input))
                .collect();
            ranked_inputs.sort_by_key(|input| input_importance(input));
            let (assumptions, assumptions_total, assumptions_omitted) = if is_record {
                let card = attempt_card(node, &node_source_ids(node));
                let recorded = card.pointer("/recorded_conclusions/assumptions");
                let values = recorded
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .take(8)
                    .map(|assumption| {
                        json!({"name":assumption,"origin":"recorded_attempt_contract","value":assumption,"used_for":"Recorded as an assumption; scope and source are not independently inferred.","impact":Value::Null,"source_ids":node_source_ids(node).into_iter().take(2).collect::<Vec<_>>()})
                    })
                    .collect::<Vec<_>>();
                let total = card
                    .pointer("/recorded_conclusions/omitted_counts/assumptions")
                    .and_then(Value::as_u64)
                    .unwrap_or(0) as usize
                    + values.len();
                let omitted = total.saturating_sub(values.len());
                (values, total, omitted)
            } else {
                let total = ranked_inputs.len();
                let values = ranked_inputs
                    .drain(..)
                    .take(5)
                    .map(|input| {
                        json!({
                            "name":input.get("name"),
                            "origin":input.get("origin"),
                            "value":clip(&value_text(input.get("value")),90),
                            "used_for":clip(&value_text(input.get("used_for")),140),
                            "impact":input.get("impact"),
                            "source_ids":source_ref_ids(input).into_iter().take(4).collect::<Vec<_>>()
                        })
                    })
                    .collect::<Vec<_>>();
                let omitted = total.saturating_sub(values.len());
                (values, total, omitted)
            };
            let record_card = is_record.then(|| attempt_card(node, &node_source_ids(node)));
            let derivations: Vec<Value> = semantic
                .get("derivations")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .take(4)
                .map(|derivation| {
                    json!({
                        "id":derivation.get("id"),
                        "label":derivation.get("label"),
                        "premise_ids":derivation.get("premise_ids"),
                        "upstream_claim_ids":derivation.get("upstream_claim_ids"),
                        "result_ids":derivation.get("result_ids"),
                        "argument":derivation.get("argument"),
                        "scope":derivation.get("scope")
                    })
                })
                .collect();
            let recorded_scope = if matches!(node.kind.as_str(), "attempt" | "contract") {
                let scope = attempt_scope(node);
                (!scope.is_empty()).then_some(scope)
            } else {
                assessment_scope(node).map(str::to_owned).or_else(|| {
                    node.attrs
                        .pointer("/ledger/canonical_status")
                        .and_then(Value::as_str)
                        .map(str::to_owned)
                })
            };
            let summary = semantic
                .get("overview_de")
                .and_then(Value::as_str)
                .filter(|value| !value.is_empty())
                .unwrap_or(&node.summary);
            upstream_sets.push(upstream);
            seed_sets.push(seeds);
            source_sets.push(sources);
            graph_link_sets.push(graph_link_ids);
            graph_link_signature_sets.push(graph_link_signatures);
            items.push(json!({
                "id":node.id,
                "type":node.kind,
                "label":clip(&node.label,160),
                "scope":recorded_scope,
                "summary":clip(summary,700),
                "assessment":node.attrs.get("assessment").cloned().unwrap_or(Value::Null),
                "ledger_status":node.attrs.pointer("/ledger/status").cloned().unwrap_or(Value::Null),
                "verdict":if is_record { attempt_verdict(node) } else { None },
                "verdict_raw":if is_record { node.attrs.get("original_contract").and_then(|contract|contract.get("verdict_raw")).or_else(||node.attrs.get("verdict_raw")).cloned().unwrap_or(Value::Null) } else { Value::Null },
                "enrichment_state":node.attrs.get("enrichment_state").cloned().unwrap_or(Value::Null),
                "upstream_claim_ids":upstream_values.into_iter().take(40).collect::<Vec<_>>(),
                "seeds":seed_values.into_iter().take(40).collect::<Vec<_>>(),
                "assumptions":assumptions,
                "assumptions_total":assumptions_total,
                "assumptions_omitted":assumptions_omitted,
                "recorded_conclusions":record_card.as_ref().and_then(|card|card.get("recorded_conclusions")).cloned().unwrap_or(Value::Null),
                "original_source_scope_correction":record_card.as_ref().and_then(|card|card.get("original_source_scope_correction")).cloned().unwrap_or(Value::Null),
                "derivations":derivations,
                "source_ids":source_values.into_iter().take(20).collect::<Vec<_>>(),
                "recorded_graph_links":links,
                "topic_contents":topic_contents
            }));
        }
        let shared_upstream = set_intersection_strings(&upstream_sets);
        let shared_seeds = set_intersection_strings(&seed_sets);
        let shared_graph_ids = set_intersection_strings(&graph_link_sets);
        let shared_graph_nodes: Vec<Value> = shared_graph_ids
            .iter()
            .filter_map(|id| by_id.get(id.as_str()).copied())
            .take(20)
            .map(|node| json!({"id":node.id,"type":node.kind,"label":clip(&node.label,140)}))
            .collect();
        let mut shared_graph_links: Vec<Value> = shared_graph_link_records
            .into_iter()
            .filter(|(signature, _)| {
                graph_link_signature_sets
                    .iter()
                    .all(|item_links| item_links.contains(signature))
            })
            .map(|(_, link)| link)
            .collect();
        shared_graph_links.sort_by(|a, b| {
            a["node_id"]
                .as_str()
                .cmp(&b["node_id"].as_str())
                .then_with(|| a["edge_type"].as_str().cmp(&b["edge_type"].as_str()))
                .then_with(|| a["direction"].as_str().cmp(&b["direction"].as_str()))
        });
        shared_graph_links.truncate(20);
        let brief = json!({
            "ids":ids,
            "items":items,
            "shared_explicit":{"upstream_claim_ids":shared_upstream,"seeds":shared_seeds,"recorded_graph_nodes":shared_graph_nodes,"recorded_graph_links":shared_graph_links},
            "divergent_assumptions":"Compare each item's recorded input list; different assumptions do not imply contradiction unless the source semantics establish comparable scope.",
            "status_comparison":"Statuses and scope labels are shown as recorded and may refer to different subclaims or conditions.",
            "shared_source_ids":set_intersection_strings(&source_sets),
            "notice":"Comparison aid only. Shared upstream references, seeds, source files, or status labels do not establish a logical conclusion or proof.",
            "max_chars":max_chars
        });
        Ok(fit_compare_budget(brief, max_chars))
    }
    pub async fn theme(
        &self,
        id: &str,
        limit: usize,
        offset: usize,
        include_related: bool,
    ) -> Result<Value> {
        if limit == 0 || limit > 100 {
            bail!("limit must be between 1 and 100");
        }
        let g = self.snapshot().await?;
        let topic = g
            .nodes
            .iter()
            .find(|n| n.id == id && n.kind == "topic")
            .ok_or_else(|| anyhow!("topic not found: {id}"))?;
        let by_id: HashMap<&str, &Node> = g.nodes.iter().map(|n| (n.id.as_str(), n)).collect();
        let mut members: HashSet<String> = HashSet::new();
        let mut topic_links: HashSet<String> = HashSet::from([id.to_owned()]);
        let topic_attempt_link_basis = topic_linked_attempts(&g, id, &by_id);
        for e in &g.edges {
            if e.kind == "theme_links" || e.kind == "theme_link" || e.kind == "about" {
                if e.src == id {
                    if include_related
                        || !is_inherited_attempt_topic_link(e, by_id.get(e.dst.as_str()).copied())
                    {
                        members.insert(e.dst.clone());
                    }
                }
                if e.dst == id {
                    if include_related
                        || !is_inherited_attempt_topic_link(e, by_id.get(e.src.as_str()).copied())
                    {
                        members.insert(e.src.clone());
                    }
                }
            }
            if (e.kind == "theme_links" || e.kind == "theme_link") && (e.src == id || e.dst == id) {
                topic_links.insert(if e.src == id {
                    e.dst.clone()
                } else {
                    e.src.clone()
                });
            }
        }
        for n in &g.nodes {
            if n.kind == "claim" {
                if let Some(ids) = n
                    .attrs
                    .pointer("/semantic/topic_ids")
                    .and_then(Value::as_array)
                {
                    if ids
                        .iter()
                        .any(|v| v.as_str().is_some_and(|s| topic_links.contains(s)))
                    {
                        members.insert(n.id.clone());
                    }
                }
            }
        }
        let related_claim_ids: HashSet<String> = members
            .iter()
            .filter_map(|member_id| by_id.get(member_id.as_str()).copied())
            .filter(|node| node.kind == "claim" && claim_topic_role(node, id) == Some("related"))
            .map(|node| node.id.clone())
            .collect();
        if !include_related {
            for claim_id in &related_claim_ids {
                members.remove(claim_id);
            }
        }
        // Historical topic-linked attempts remain visible even when related
        // claims are filtered from the claim lists. Keep their link basis so
        // direct topic records stay distinguishable from inherited membership.
        members.extend(topic_attempt_link_basis.keys().cloned());
        // Attempts enter a theme only through an explicit topic link or a tests edge to a member claim.
        for e in &g.edges {
            if e.kind == "tests" {
                let (attempt, claim) = if by_id
                    .get(e.src.as_str())
                    .is_some_and(|n| n.kind == "attempt")
                {
                    (&e.src, &e.dst)
                } else if by_id
                    .get(e.dst.as_str())
                    .is_some_and(|n| n.kind == "attempt")
                {
                    (&e.dst, &e.src)
                } else {
                    continue;
                };
                if members.contains(claim) {
                    members.insert(attempt.clone());
                }
            }
        }
        let categories = [
            "recorded_established",
            "conditional",
            "open",
            "mixed",
            "excluded",
            "assumed",
            "historical",
            "methodological_rule",
            "unknown",
        ];
        let mut grouped: HashMap<&str, Vec<&Node>> =
            categories.iter().map(|c| (*c, Vec::new())).collect();
        let mut attempts = Vec::new();
        let mut structures = Vec::new();
        let mut goals = Vec::new();
        let mut contracts = Vec::new();
        for n in &g.nodes {
            if !members.contains(&n.id) {
                continue;
            }
            match n.kind.as_str() {
                "attempt" => attempts.push(n),
                "structure" => structures.push(n),
                "goal" => goals.push(n),
                "contract" => contracts.push(n),
                "claim" => {
                    let category = n
                        .attrs
                        .pointer("/assessment/category")
                        .and_then(Value::as_str)
                        .unwrap_or("unknown");
                    grouped.entry(category).or_default().push(n);
                }
                _ => {}
            }
        }
        let cap_nodes = |nodes: &[&Node]| -> Vec<Value> {
            nodes
                .iter()
                .skip(offset)
                .take(limit)
                .map(|n| compact_theme_node(n, id))
                .collect()
        };
        let mut summary = serde_json::Map::new();
        for category in categories {
            summary.insert(
                category.to_owned(),
                json!(grouped.get(category).map(Vec::len).unwrap_or(0)),
            );
        }
        summary.insert("attempts".into(), json!(attempts.len()));
        summary.insert(
            "related_claims_total".into(),
            json!(related_claim_ids.len()),
        );
        let attempt_page: Vec<Value> = attempts
            .iter()
            .skip(offset)
            .take(limit)
            .map(|n| {
                let mut card = compact_theme_node(n, id);
                card["topic_link_basis"] = json!(topic_attempt_link_basis
                    .get(&n.id)
                    .map(String::as_str)
                    .unwrap_or("tests_primary_topic_claim"));
                card
            })
            .collect();
        let mut result = json!({"topic":compact_theme_node(topic,id),"summary":summary,"related_claims_total":related_claim_ids.len(),"include_related":include_related,"membership_notice":"Claims explicitly marked related are excluded by default; claims with no topic role remain included. Explicit and inherited topic-linked attempts remain visible with their link basis. Attempts linked only by tests edges to excluded related claims are omitted. Related membership is an organizational label, not a change to claim status.","node_projection":"Topic lists include display fields and evidence references; use the node dossier for full ledger, semantic annotations, and original contract text.","solved":cap_nodes(grouped.get("recorded_established").map(Vec::as_slice).unwrap_or(&[])),"conditional":cap_nodes(grouped.get("conditional").map(Vec::as_slice).unwrap_or(&[])),"open":cap_nodes(grouped.get("open").map(Vec::as_slice).unwrap_or(&[])),"mixed":cap_nodes(grouped.get("mixed").map(Vec::as_slice).unwrap_or(&[])),"excluded":cap_nodes(grouped.get("excluded").map(Vec::as_slice).unwrap_or(&[])),"assumed":cap_nodes(grouped.get("assumed").map(Vec::as_slice).unwrap_or(&[])),"historical":cap_nodes(grouped.get("historical").map(Vec::as_slice).unwrap_or(&[])),"methodological_rule":cap_nodes(grouped.get("methodological_rule").map(Vec::as_slice).unwrap_or(&[])),"unknown":cap_nodes(grouped.get("unknown").map(Vec::as_slice).unwrap_or(&[])),"attempts":attempt_page,"structures":cap_nodes(&structures),"goals":cap_nodes(&goals),"contracts":cap_nodes(&contracts),"totals":{"structures":structures.len(),"goals":goals.len(),"contracts":contracts.len()},"pagination":{}});
        // Keep per-category totals explicit even when a returned list is capped.
        for category in categories {
            result["totals"][category] = json!(grouped.get(category).map(Vec::len).unwrap_or(0));
            let total = grouped.get(category).map(Vec::len).unwrap_or(0);
            let result_key = if category == "recorded_established" {
                "solved"
            } else {
                category
            };
            result["pagination"][category] = json!({"offset":offset,"limit":limit,"total":total,"truncated":offset+result[result_key].as_array().map(Vec::len).unwrap_or(0)<total});
        }
        for (name, total) in [
            ("attempts", attempts.len()),
            ("structures", structures.len()),
            ("goals", goals.len()),
            ("contracts", contracts.len()),
        ] {
            result["pagination"][name] = json!({"offset":offset,"limit":limit,"total":total,"truncated":offset+result[name].as_array().map(Vec::len).unwrap_or(0)<total});
        }
        Ok(result)
    }
    pub async fn solution_chain(&self, id: &str, depth: usize) -> Result<Value> {
        if depth > 4 {
            bail!("depth must be at most 4");
        }
        const CAP: usize = 80;
        let g = self.snapshot().await?;
        let by_id: HashMap<&str, &Node> = g.nodes.iter().map(|n| (n.id.as_str(), n)).collect();
        let root = *by_id
            .get(id)
            .ok_or_else(|| anyhow!("node not found: {id}"))?;
        let research_root = root.kind.starts_with("research_");
        let proof_types: HashSet<&str> = [
            "requires",
            "required_together",
            "derives",
            "has_result",
            "has_derivation",
            "uses_claim_under_scope",
            "implements",
            "checks",
            "formalization_candidate",
            "has_finding",
            "has_task",
            "requires_task",
            "requires_premise",
            "resolves_in_scope",
            "evaluated_by",
        ]
        .into_iter()
        .collect();
        let allowed_edge = |e: &Edge| {
            proof_types.contains(e.kind.as_str())
                && !(research_root
                    && root.kind != "research_problem"
                    && matches!(e.kind.as_str(), "has_finding" | "has_task"))
        };
        let mut visited = HashSet::from([id.to_owned()]);
        let mut queue = VecDeque::from([(id.to_owned(), 0usize)]);
        let mut selected_edges = HashSet::new();
        let mut truncated = false;
        let mut truncation_reasons: HashSet<&str> = HashSet::new();
        let mut affected_derivations: HashSet<String> = HashSet::new();
        while let Some((current, level)) = queue.pop_front() {
            // A shared premise is an input, not a bridge to every task that mentions it.
            if research_root
                && current != id
                && by_id
                    .get(current.as_str())
                    .is_some_and(|n| n.kind == "research_premise")
            {
                continue;
            }
            if level >= depth {
                for e in &g.edges {
                    if !allowed_edge(e) {
                        continue;
                    }
                    let next = if e.src == current {
                        Some(e.dst.as_str())
                    } else if e.dst == current {
                        Some(e.src.as_str())
                    } else {
                        None
                    };
                    if let Some(next) = next.filter(|next| !visited.contains(*next)) {
                        if by_id.get(next).is_some_and(|n| n.kind == "topic") {
                            continue;
                        }
                        if by_id
                            .get(current.as_str())
                            .is_some_and(|n| n.kind == "derivation")
                        {
                            affected_derivations.insert(current.clone());
                        }
                        if by_id.get(next).is_some_and(|n| n.kind == "derivation") {
                            affected_derivations.insert(next.to_owned());
                        }
                        truncated = true;
                        truncation_reasons.insert("depth_limit");
                    }
                }
                continue;
            }
            for (i, e) in g.edges.iter().enumerate() {
                if !allowed_edge(e) {
                    continue;
                }
                let next = if e.src == current {
                    Some(e.dst.as_str())
                } else if e.dst == current {
                    Some(e.src.as_str())
                } else {
                    None
                };
                let Some(next) = next else { continue };
                if visited.contains(next) {
                    selected_edges.insert(i);
                    continue;
                }
                if by_id.get(next).is_some_and(|n| n.kind == "topic") {
                    continue;
                }
                if visited.len() >= CAP {
                    truncated = true;
                    truncation_reasons.insert("node_limit");
                    if by_id
                        .get(current.as_str())
                        .is_some_and(|n| n.kind == "derivation")
                    {
                        affected_derivations.insert(current.clone());
                    }
                    if by_id.get(next).is_some_and(|n| n.kind == "derivation") {
                        affected_derivations.insert(next.to_owned());
                    }
                    continue;
                }
                selected_edges.insert(i);
                visited.insert(next.to_owned());
                queue.push_back((next.to_owned(), level + 1));
            }
        }
        let nodes: Vec<Value> = g
            .nodes
            .iter()
            .filter(|n| visited.contains(&n.id))
            .map(chain_public_node)
            .collect();
        if selected_edges.len() > 300 {
            truncated = true;
            truncation_reasons.insert("edge_limit");
            for (_, edge) in g
                .edges
                .iter()
                .enumerate()
                .filter(|(i, _)| selected_edges.contains(i))
                .skip(300)
            {
                for endpoint in [&edge.src, &edge.dst] {
                    if by_id
                        .get(endpoint.as_str())
                        .is_some_and(|n| n.kind == "derivation")
                    {
                        affected_derivations.insert(endpoint.clone());
                    }
                }
            }
        }
        let edges: Vec<Edge> = g
            .edges
            .iter()
            .enumerate()
            .filter(|(i, _)| selected_edges.contains(i))
            .take(300)
            .map(|(_, e)| e.clone())
            .collect();
        let references: Vec<Edge> = g
            .edges
            .iter()
            .filter(|e| {
                e.kind == "references" && (visited.contains(&e.src) || visited.contains(&e.dst))
            })
            .take(200)
            .cloned()
            .collect();
        let affected:Vec<Value>=affected_derivations.iter().map(|id|json!({"id":id,"label":by_id.get(id.as_str()).map(|n|n.label.as_str()).unwrap_or(""),"all_premises_complete":false})).collect();
        let reasons: Vec<&str> = truncation_reasons.into_iter().collect();
        let reason = if truncated {
            Some("Traversal reached its node, edge, or depth bound; affected derivations may have omitted premises.")
        } else {
            None
        };
        Ok(
            json!({"root":chain_public_node(root),"depth":depth,"nodes":nodes,"edges":edges,"references":references,"truncated":truncated,"traversal_complete":!truncated,"reason":reason,"truncation_reasons":reasons,"all_premises_complete":if truncated { Value::Bool(false) } else { Value::Null },"affected_derivations":affected,"node_limit":CAP,"interpretation":{"edge_types":"Stored graph relationships; not proof certificates.","structural_edges":["requires","has_result","has_derivation"],"references_separate":true,"research_traversal_scope":if research_root {json!("Premises are terminal inputs; task/finding roots do not expand through problem ownership into unrelated branches.")}else{Value::Null},"traversal_completeness":"Whether this bounded graph traversal omitted reachable stored reasoning edges.","premise_completeness":"Not assessed; a complete traversal does not establish that all original premises were recorded, sufficient, or valid."}}),
        )
    }
    pub async fn list_attempts(
        &self,
        id: &str,
        status: Option<&str>,
        limit: usize,
        offset: usize,
    ) -> Result<Value> {
        if limit == 0 || limit > 100 {
            bail!("limit must be between 1 and 100");
        }
        let g = self.snapshot().await?;
        let node = g
            .nodes
            .iter()
            .find(|n| n.id == id)
            .ok_or_else(|| anyhow!("node not found: {id}"))?;
        let mut ids = HashSet::new();
        let mut topic_attempt_link_basis = HashMap::new();
        if node.kind == "topic" {
            let by_id: HashMap<&str, &Node> = g
                .nodes
                .iter()
                .map(|node| (node.id.as_str(), node))
                .collect();
            topic_attempt_link_basis = topic_linked_attempts(&g, id, &by_id);
            ids.extend(topic_attempt_link_basis.keys().cloned());
            let theme = self.theme(id, 100, 0, false).await?;
            if let Some(items) = theme["attempts"].as_array() {
                for item in items {
                    if let Some(s) = item["id"].as_str() {
                        ids.insert(s.to_owned());
                    }
                }
            }
            // theme's displayed page is capped; preserve all direct and inherited
            // topic links even when the caller requests a later attempt page.
            for e in &g.edges {
                if e.kind == "about" && (e.src == id || e.dst == id) {
                    let other = if e.src == id { &e.dst } else { &e.src };
                    let other_node = g.nodes.iter().find(|n| n.id == *other);
                    if other_node.is_some_and(|n| n.kind == "attempt")
                        && !is_inherited_attempt_topic_link(e, other_node)
                    {
                        ids.insert(other.clone());
                    }
                }
            }
            let claim_ids: HashSet<&str> = g
                .nodes
                .iter()
                .filter(|n| {
                    n.kind == "claim"
                        && claim_topic_role(n, id) != Some("related")
                        && (n
                            .attrs
                            .pointer("/semantic/topic_ids")
                            .and_then(Value::as_array)
                            .is_some_and(|a| a.iter().any(|v| v.as_str() == Some(id)))
                            || g.edges.iter().any(|e| {
                                e.kind == "about"
                                    && ((e.src == id && e.dst == n.id)
                                        || (e.dst == id && e.src == n.id))
                            }))
                })
                .map(|n| n.id.as_str())
                .collect();
            for e in &g.edges {
                if e.kind == "tests" {
                    let (a, c) = if g.nodes.iter().any(|n| n.id == e.src && n.kind == "attempt") {
                        (&e.src, &e.dst)
                    } else {
                        (&e.dst, &e.src)
                    };
                    if claim_ids.contains(c.as_str()) {
                        ids.insert(a.clone());
                    }
                }
            }
        } else if node.kind == "claim" {
            for e in &g.edges {
                if e.kind == "tests" {
                    let other = if e.src == id {
                        &e.dst
                    } else if e.dst == id {
                        &e.src
                    } else {
                        continue;
                    };
                    if g.nodes.iter().any(|n| {
                        n.id == *other && matches!(n.kind.as_str(), "attempt" | "contract")
                    }) {
                        ids.insert(other.clone());
                    }
                }
            }
        } else {
            bail!("list_attempts requires a topic or claim id");
        }
        let mut attempts: Vec<&Node> = g
            .nodes
            .iter()
            .filter(|n| {
                ids.contains(&n.id)
                    && matches!(n.kind.as_str(), "attempt" | "contract")
                    && status.map_or(true, |s| attempt_matches_status(n, s))
            })
            .collect();
        if node.kind == "claim" {
            let node_kinds: HashMap<&str, &str> = g
                .nodes
                .iter()
                .map(|candidate| (candidate.id.as_str(), candidate.kind.as_str()))
                .collect();
            let mut tested_claim_counts = HashMap::<&str, usize>::new();
            for edge in g.edges.iter().filter(|edge| edge.kind == "tests") {
                let pair = if node_kinds
                    .get(edge.src.as_str())
                    .is_some_and(|kind| matches!(*kind, "attempt" | "contract"))
                    && node_kinds
                        .get(edge.dst.as_str())
                        .is_some_and(|kind| *kind == "claim")
                {
                    Some((edge.src.as_str(), edge.dst.as_str()))
                } else if node_kinds
                    .get(edge.dst.as_str())
                    .is_some_and(|kind| matches!(*kind, "attempt" | "contract"))
                    && node_kinds
                        .get(edge.src.as_str())
                        .is_some_and(|kind| *kind == "claim")
                {
                    Some((edge.dst.as_str(), edge.src.as_str()))
                } else {
                    None
                };
                if let Some((attempt_id, _)) = pair {
                    *tested_claim_counts.entry(attempt_id).or_default() += 1;
                }
            }
            attempts.sort_by(|a, b| {
                tested_claim_counts
                    .get(a.id.as_str())
                    .unwrap_or(&0)
                    .cmp(tested_claim_counts.get(b.id.as_str()).unwrap_or(&0))
                    .then_with(|| a.id.cmp(&b.id))
            });
        }
        let total = attempts.len();
        let page: Vec<Value> = attempts
            .iter()
            .skip(offset)
            .take(limit)
            .map(|n| {
                if node.kind == "topic" {
                    let mut card = compact_theme_node(n, id);
                    card["topic_link_basis"] = json!(topic_attempt_link_basis
                        .get(&n.id)
                        .map(String::as_str)
                        .unwrap_or("tests_primary_topic_claim"));
                    card
                } else {
                    let mut card = compact_theme_node(n, id);
                    card["claim_link_basis"] = json!("Direct tests edge; association does not establish chronology, causality, or proof.");
                    card
                }
            })
            .collect();
        let mut response = json!({"id":id,"attempts":page,"total":total,"offset":offset,"limit":limit,"truncated":offset+page.len()<total,"status_filter":status});
        if node.kind == "claim" {
            response["claim_link_basis"] = json!("Only direct tests edges are listed; these records are not inferred chronological predecessors or proof of the claim.");
        }
        Ok(response)
    }
    pub async fn search(&self, req: SearchRequest) -> Result<Value> {
        let g = self.snapshot().await?;
        let limit = req.limit.unwrap_or(20);
        if limit == 0 || limit > 50 {
            bail!("limit must be between 1 and 50");
        }
        let query = req.query.trim();
        if query.len() > 2000 {
            bail!("query exceeds 2000 characters");
        }
        let lexical: Vec<String> = tokenize(query);
        let query_vector = if !lexical.is_empty() && g.nodes.iter().any(|n| embedding(n).is_some())
        {
            self.query_embedding(query).await.ok()
        } else {
            None
        };
        let search_mode = if query_vector.is_some() {
            "hybrid"
        } else {
            "lexical"
        };
        let query_lower = query.to_lowercase();
        let raw_claim_query = query_lower.strip_prefix("claim:").unwrap_or(&query_lower);
        let raw_contract_query = query_lower
            .strip_prefix("contract:")
            .or_else(|| query_lower.strip_prefix("attempt:"))
            .unwrap_or(&query_lower);
        let mut scored: Vec<(u8, f32, &Node)> = g
            .nodes
            .iter()
            .filter(|n| req.kind.as_deref().map_or(true, |k| k == n.kind))
            .filter_map(|n| {
                let node_id = n.id.to_lowercase();
                let raw_claim_id =
                    n.id.strip_prefix("claim:")
                        .map(str::to_lowercase)
                        .unwrap_or_default();
                let label = n.label.to_lowercase();
                let canonical_contract_name = n
                    .attrs
                    .get("name")
                    .or_else(|| n.attrs.get("contract_id"))
                    .and_then(Value::as_str)
                    .map(str::to_lowercase)
                    .or_else(|| {
                        n.id.strip_prefix("contract:")
                            .or_else(|| n.id.strip_prefix("attempt:"))
                            .map(str::to_lowercase)
                    });
                let exact_contract_id = matches!(n.kind.as_str(), "contract" | "attempt")
                    && canonical_contract_name
                        .as_deref()
                        .is_some_and(|name| name == raw_contract_query);
                let exact_id = node_id == query_lower
                    || (n.kind == "claim" && raw_claim_id == query_lower)
                    || (query_lower.starts_with("claim:")
                        && n.kind == "claim"
                        && raw_claim_id == raw_claim_query);
                let id_prefix = node_id.starts_with(&query_lower)
                    || (n.kind == "claim"
                        && raw_claim_id.starts_with(raw_claim_query)
                        && !raw_claim_query.is_empty());
                let label_match =
                    !lexical.is_empty() && lexical.iter().all(|token| label.contains(token));
                let priority = if exact_id || exact_contract_id {
                    0
                } else if id_prefix || label_match {
                    1
                } else {
                    2
                };
                let paper = crate::papers::is_paper(n);
                let paper_source_current = crate::papers::source_current(n);
                let paper_summary_current = crate::papers::summary_current(n, paper_source_current);
                let semantic = n.attrs.get("semantic");
                let translation = n.attrs.pointer("/translations/en");
                let title = [
                    semantic.and_then(|v| v.get("title_de")),
                    semantic.and_then(|v| v.get("title_en")),
                    n.attrs.get("label_en"),
                    translation.and_then(|v| v.get("title")),
                ]
                .into_iter()
                .flatten()
                .map(|v| value_text(Some(v)))
                .collect::<Vec<_>>()
                .join(" ");
                let overview = [
                    semantic.and_then(|v| v.get("overview_de")),
                    semantic.and_then(|v| v.get("overview_en")),
                    n.attrs.get("summary_en"),
                    translation.and_then(|v| v.get("overview")),
                    n.attrs.pointer("/topic_orientation/role_en"),
                ]
                .into_iter()
                .flatten()
                .map(|v| value_text(Some(v)))
                .collect::<Vec<_>>()
                .join(" ");
                let interpretation = [
                    semantic.and_then(|v| v.get("interpretation_de")),
                    semantic.and_then(|v| v.get("interpretation_en")),
                ]
                .into_iter()
                .flatten()
                .map(|v| value_text(Some(v)))
                .collect::<Vec<_>>()
                .join(" ");
                let keywords = value_text(n.attrs.get("keywords"));
                let ledger_text = [
                    n.attrs.pointer("/ledger/claim"),
                    n.attrs.pointer("/ledger/topic"),
                    n.attrs.pointer("/ledger/statement"),
                ]
                .into_iter()
                .flatten()
                .map(|value| value_text(Some(value)))
                .collect::<Vec<_>>()
                .join(" ");
                let contract_text = [
                    n.attrs.pointer("/original_contract/question"),
                    n.attrs.pointer("/original_contract/contract"),
                    n.attrs.get("question"),
                ]
                .into_iter()
                .flatten()
                .map(|value| value_text(Some(value)))
                .collect::<Vec<_>>()
                .join(" ");
                let id_text = n.id.to_lowercase();
                let label_text = n.label.to_lowercase();
                let summary_text = n.summary.to_lowercase();
                let summary_en_text = value_text(n.attrs.get("summary_en")).to_lowercase();
                let translated_overview_text =
                    value_text(n.attrs.pointer("/translations/en/overview")).to_lowercase();
                let title_text = title.to_lowercase();
                let overview_text = overview.to_lowercase();
                let interpretation_text = interpretation.to_lowercase();
                let keywords_text = keywords.to_lowercase();
                let ledger_text = ledger_text.to_lowercase();
                let contract_text = contract_text.to_lowercase();
                let paper_text = crate::papers::searchable(n, paper_source_current).to_lowercase();
                // Presence-based field scoring prevents repeated text inside a
                // large contract from overwhelming concise identity/context.
                let lexical_score: f32 = lexical
                    .iter()
                    .map(|token| {
                        if id_text.contains(token) {
                            9.0
                        } else if label_text.contains(token) {
                            7.0
                        } else if paper {
                            if paper_text.contains(token) { 1.5 } else { 0.0 }
                        } else if summary_text.contains(token)
                            || summary_en_text.contains(token)
                            || translated_overview_text.contains(token)
                        {
                            3.0
                        } else if ledger_text.contains(token) || contract_text.contains(token) {
                            2.5
                        } else if title_text.contains(token) || keywords_text.contains(token) {
                            2.0
                        } else if overview_text.contains(token)
                            || interpretation_text.contains(token)
                        {
                            1.5
                        } else {
                            0.0
                        }
                    })
                    .sum();
                let semantic_score = query_vector
                    .as_ref()
                    .filter(|_| !paper || paper_summary_current)
                    .and_then(|q| embedding(n).map(|e| cosine(q, &e).max(0.0) * 5.0))
                    .unwrap_or(0.0);
                let score = if lexical.is_empty() {
                    1.0
                } else {
                    lexical_score + semantic_score
                };
                (score > 0.0).then_some((priority, score, n))
            })
            .collect();
        scored.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then_with(|| b.1.total_cmp(&a.1))
                .then_with(|| a.2.id.cmp(&b.2.id))
        });
        let attempts_by_contract: HashMap<String, &Node> = g
            .nodes
            .iter()
            .filter(|node| node.kind == "attempt")
            .filter_map(|node| {
                let contract_id = node
                    .attrs
                    .get("name")
                    .and_then(Value::as_str)
                    .map(str::to_owned)
                    .or_else(|| node.id.strip_prefix("contract:").map(str::to_owned))?;
                Some((contract_id, node))
            })
            .collect();
        let question_parent_by_id: HashMap<&str, &Node> = g
            .nodes
            .iter()
            .filter(|node| node.kind == "contract_question")
            .filter_map(|question| {
                let contract_id = question.attrs.get("contract_id").and_then(Value::as_str)?;
                attempts_by_contract
                    .get(contract_id)
                    .map(|attempt| (question.id.as_str(), *attempt))
            })
            .collect();
        let scored_attempt_ids: HashSet<&str> = scored
            .iter()
            .filter(|(_, _, node)| node.kind == "attempt")
            .map(|(_, _, node)| node.id.as_str())
            .collect();
        let suppressed_question_ids: HashSet<&str> = question_parent_by_id
            .iter()
            .filter_map(|(question_id, parent)| {
                scored_attempt_ids
                    .contains(parent.id.as_str())
                    .then_some(*question_id)
            })
            .collect();
        let pattern_hits = if req
            .kind
            .as_deref()
            .is_none_or(|kind| kind == "theory_pattern")
        {
            self.matching_theory_patterns(query, limit).await
        } else {
            Vec::new()
        };
        let (results, total) = if req.kind.as_deref() == Some("theory_pattern") {
            let total = pattern_hits.len();
            (pattern_hits, total)
        } else {
            let mut combined: Vec<(u8, f32, String, Value)> = scored
                .into_iter()
                .filter(|(_, _, node)| !suppressed_question_ids.contains(node.id.as_str()))
                .map(|(priority, score, node)| {
                    let result = compact_search_result(
                        node,
                        score,
                        question_for_attempt(node, &question_parent_by_id, &g.nodes),
                    );
                    (priority, score, node.id.clone(), result)
                })
                .collect();
            combined.extend(pattern_hits.into_iter().map(|pattern| {
                let id = pattern["id"].as_str().unwrap_or_default().to_owned();
                let score = pattern["score"].as_f64().unwrap_or(0.0) as f32;
                // Authored navigation patterns participate in relevance ranking;
                // they must not displace an exact graph-node match by position.
                (2, score, id, pattern)
            }));
            combined.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then_with(|| b.1.total_cmp(&a.1))
                    .then_with(|| a.2.cmp(&b.2))
            });
            let total = combined.len();
            let results = combined
                .into_iter()
                .take(limit)
                .map(|(_, _, _, result)| result)
                .collect();
            (results, total)
        };
        Ok(
            json!({"query":query,"count":results.len(),"total":total,"limit":limit,"truncated":total>limit,"search_mode":search_mode,"results":results}),
        )
    }
    async fn matching_theory_patterns(&self, query: &str, limit: usize) -> Vec<Value> {
        let tokens = tokenize(query);
        if tokens.is_empty() {
            return Vec::new();
        }
        let Ok(map) = self.read_theory_map().await else {
            return Vec::new();
        };
        let mut matches: Vec<(usize, Value)> = map
            .get("patterns")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|pattern| {
                let text = [
                    pattern.get("label_de"),
                    pattern.get("label_en"),
                    pattern.get("summary_de"),
                    pattern.get("summary_en"),
                ]
                .into_iter()
                .flatten()
                .map(|field| value_text(Some(field)))
                .collect::<Vec<_>>()
                .join(" ")
                .to_lowercase();
                let matched: Vec<&str> = tokens
                    .iter()
                    .map(String::as_str)
                    .filter(|token| text.contains(*token))
                    .collect();
                if matched.is_empty() {
                    return None;
                }
                let source_ids: Vec<String> = pattern
                    .get("source_refs")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(|reference| {
                        reference
                            .get("source_id")
                            .and_then(Value::as_str)
                            .or_else(|| reference.as_str())
                    })
                    .take(8)
                    .map(str::to_owned)
                    .collect();
                let label = value_text(pattern.get("label_de"));
                let summary = value_text(pattern.get("summary_de"));
                let label_en = pattern.get("label_en").and_then(Value::as_str);
                let summary_en = pattern.get("summary_en").and_then(Value::as_str);
                let result = json!({
                    "id":pattern.get("id").and_then(Value::as_str).unwrap_or(""),
                    "type":"theory_pattern",
                    "kind":pattern.get("kind"),
                    "label":clip(&label,180),
                    "label_en":label_en.map(|text|clip(text,180)),
                    "summary":clip(&summary,349),
                    "summary_en":summary_en.map(|text|clip(text,349)),
                    "reviewed":pattern.get("reviewed").and_then(Value::as_bool).unwrap_or(false),
                    "match_terms":matched,
                    "match_reason":"Direct term match in an authored theory-map pattern; this record is a navigation aid, not a proof.",
                    "source_ids":source_ids,
                    "related_ids":pattern.get("related_ids").and_then(Value::as_array).into_iter().flatten().take(8).cloned().collect::<Vec<_>>(),
                    "conditions":pattern.get("conditions").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                    "counterevidence":pattern.get("counterevidence").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                    "score":matched.len() as f32 / tokens.len() as f32
                });
                Some((matched.len(), result))
            })
            .collect();
        matches.sort_by(|a, b| {
            b.0.cmp(&a.0)
                .then_with(|| a.1["id"].as_str().cmp(&b.1["id"].as_str()))
        });
        matches
            .into_iter()
            .take(limit.min(3))
            .map(|(_, result)| result)
            .collect()
    }
    pub async fn dossier(&self, id: &str) -> Result<Value> {
        let g = self.snapshot().await?;
        let Some(node) = g.nodes.iter().find(|n| n.id == id) else {
            if let Some(pattern) = self.theory_pattern_by_id(id).await? {
                let node = theory_pattern_node(&pattern)?;
                let by_id: HashMap<&str, &Node> = g
                    .nodes
                    .iter()
                    .map(|node| (node.id.as_str(), node))
                    .collect();
                let related_ids: Vec<String> = pattern
                    .get("related_ids")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(Value::as_str)
                    .filter(|related_id| by_id.contains_key(related_id))
                    .map(str::to_owned)
                    .collect();
                let related_nodes: Vec<Value> = related_ids
                    .iter()
                    .take(200)
                    .filter_map(|related_id| by_id.get(related_id.as_str()))
                    .map(|node| public_node(node))
                    .collect();
                let edges: Vec<Edge> = related_ids
                    .iter()
                    .map(|related_id| {
                        let related = by_id[related_id.as_str()];
                        Edge {
                            id: format!("{}:about:{}", node.id, related.id),
                            src: node.id.clone(),
                            dst: related.id.clone(),
                            kind: if related.kind == "topic" {
                                "about_topic".to_owned()
                            } else {
                                "related_to".to_owned()
                            },
                            attrs: json!({
                                "epistemic_role":"authored_navigation_reference",
                                "reviewed":pattern.get("reviewed").and_then(Value::as_bool).unwrap_or(false)
                            }),
                        }
                    })
                    .take(200)
                    .collect();
                let sources: Vec<Source> = pattern
                    .get("source_refs")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(theory_pattern_source)
                    .take(200)
                    .collect();
                return Ok(json!({
                    "node":public_node(&node),
                    "related_nodes":related_nodes,
                    "edges":edges,
                    "sources":sources,
                    "external_sources":pattern.get("source_refs").and_then(Value::as_array).into_iter().flatten().filter(|reference| reference.get("source_location_state").and_then(Value::as_str) == Some("moved_to_external_repository")).map(|reference| json!({"source_url":reference.get("source_url"),"repository":reference.get("external_repository"),"revision":reference.get("external_revision"),"path":reference.get("external_path"),"sha256":reference.get("external_sha256"),"line_start":reference.get("line_start"),"line_end":reference.get("line_end"),"reference_role":reference.get("reference_role"),"source_location_state":reference.get("source_location_state")})).collect::<Vec<_>>(),
                    "counts":{"related_nodes":related_ids.len(),"edges":related_ids.len(),"sources":sources.len()},
                    "truncated":related_ids.len()>200,
                    "notice":"Theory-map pattern and related records are navigation metadata; these references do not establish proof."
                }));
            }
            return Err(anyhow!("node not found: {id}"));
        };
        let all_edges: Vec<_> = g
            .edges
            .iter()
            .filter(|e| e.src == id || e.dst == id)
            .cloned()
            .collect();
        let related_ids: HashSet<_> = all_edges
            .iter()
            .map(|e| {
                if e.src == id {
                    e.dst.as_str()
                } else {
                    e.src.as_str()
                }
            })
            .collect();
        let all_related_nodes: Vec<_> = g
            .nodes
            .iter()
            .filter(|n| related_ids.contains(n.id.as_str()))
            .collect();
        let related_nodes: Vec<_> = all_related_nodes
            .iter()
            .take(200)
            .map(|n| public_node(n))
            .collect();
        let mut linked_sources: HashSet<String> = all_edges
            .iter()
            .filter_map(|e| {
                if e.src == id && is_source_type(&e.kind) {
                    Some(e.dst.clone())
                } else if e.dst == id && is_source_type(&e.kind) {
                    Some(e.src.clone())
                } else {
                    None
                }
            })
            .collect();
        for related in std::iter::once(node).chain(all_related_nodes.iter().copied()) {
            if let Some(refs) = related.attrs.get("source_refs").and_then(Value::as_array) {
                for r in refs {
                    if let Some(s) = r.as_str() {
                        linked_sources.insert(s.to_owned());
                    }
                }
            }
        }
        let id_sources: HashSet<String> = g
            .sources
            .iter()
            .filter(|s| s.id == id)
            .map(|s| s.id.clone())
            .collect();
        let sources: Vec<_> = g
            .sources
            .iter()
            .filter(|s| linked_sources.contains(&s.id) || id_sources.contains(&s.id))
            .take(200)
            .cloned()
            .collect();
        let edges: Vec<_> = all_edges.iter().take(200).cloned().collect();
        let source_count = g
            .sources
            .iter()
            .filter(|s| linked_sources.contains(&s.id) || id_sources.contains(&s.id))
            .count();
        let edge_count = all_edges.len();
        let related_count = all_related_nodes.len();
        let truncated = edge_count > 200 || related_count > 200 || source_count > 200;
        Ok(
            json!({"node":public_node(node),"related_nodes":related_nodes,"edges":edges,"sources":sources,"counts":{"related_nodes":related_count,"edges":edge_count,"sources":source_count},"truncated":truncated}),
        )
    }

    /// A small MCP-facing projection. The HTTP dossier remains the full record;
    /// clients can request the complete MCP form explicitly when needed.
    pub async fn dossier_compact(
        &self,
        id: &str,
        max_chars: usize,
        neighbor_limit: usize,
    ) -> Result<Value> {
        if !(1_000..=500_000).contains(&max_chars) {
            bail!("max_chars must be between 1000 and 500000");
        }
        if neighbor_limit == 0 || neighbor_limit > 50 {
            bail!("neighbor_limit must be between 1 and 50");
        }
        let dossier = self.dossier(id).await?;
        let compact = compact_dossier_projection(&dossier, max_chars, neighbor_limit);
        if json_char_count(&compact) > max_chars {
            bail!("compact dossier exceeds max_chars; request a smaller neighbor_limit or use the web dossier for full details");
        }
        Ok(compact)
    }
    pub async fn trace(&self, id: &str, direction: &str, depth: usize) -> Result<Value> {
        if !["inbound", "outbound", "both"].contains(&direction) {
            bail!("invalid direction");
        }
        if depth > 5 {
            bail!("depth must be at most 5");
        }
        let g = self.snapshot().await?;
        let by_id: HashMap<_, _> = g.nodes.iter().map(|n| (n.id.as_str(), n)).collect();
        let root = *by_id
            .get(id)
            .ok_or_else(|| anyhow!("node not found: {id}"))?;
        let mut visited = HashSet::from([id.to_owned()]);
        let mut queue = VecDeque::from([(id.to_owned(), 0usize)]);
        let mut edge_ids = HashSet::new();
        let mut truncated = false;
        while let Some((current, level)) = queue.pop_front() {
            if level >= depth {
                continue;
            }
            for (index, e) in g.edges.iter().enumerate() {
                let next = match direction {
                    "outbound" if e.src == current => Some(&e.dst),
                    "inbound" if e.dst == current => Some(&e.src),
                    "both" if e.src == current => Some(&e.dst),
                    "both" if e.dst == current => Some(&e.src),
                    _ => None,
                };
                if let Some(next) = next {
                    if visited.contains(next) {
                        edge_ids.insert(index);
                    } else if visited.len() < 150 {
                        edge_ids.insert(index);
                        visited.insert(next.clone());
                        if by_id
                            .get(next.as_str())
                            .is_some_and(|node| !is_context_node(&node.kind))
                        {
                            queue.push_back((next.clone(), level + 1));
                        }
                    } else {
                        truncated = true;
                    }
                }
            }
        }
        let nodes: Vec<_> = g
            .nodes
            .iter()
            .filter(|n| visited.contains(&n.id))
            .map(public_node)
            .collect();
        let edges: Vec<_> = g
            .edges
            .iter()
            .enumerate()
            .filter(|(index, _)| edge_ids.contains(index))
            .map(|(_, e)| e.clone())
            .collect();
        Ok(
            json!({"root":public_node(root),"direction":direction,"depth":depth,"nodes":nodes,"edges":edges,"truncated":truncated,"node_limit":150}),
        )
    }
    pub async fn source(
        &self,
        id: &str,
        requested_start: Option<usize>,
        requested_end: Option<usize>,
    ) -> Result<Value> {
        let g = self.snapshot().await?;
        // These handles describe the live execution plan, which can advance
        // independently of the cached scientific graph snapshot.
        let live_plan_sources = if id.starts_with("source:fundamental-plan:") {
            crate::research::current_plan_sources(&self.repo_root)
        } else {
            Vec::new()
        };
        let source_from_graph = g.sources.iter().find(|s| s.id == id);
        let source_from_pattern;
        let src = if let Some(source) = live_plan_sources.iter().find(|s| s.id == id) {
            source
        } else if let Some(source) = source_from_graph {
            source
        } else {
            source_from_pattern = self.theory_pattern_source_by_id(id).await?;
            source_from_pattern
                .as_ref()
                .ok_or_else(|| anyhow!("source not found: {id}"))?
        };
        let (text, line_start, line_end, extent) =
            self.source_text(src, requested_start, requested_end)?;
        Ok(
            json!({"id":src.id,"path":src.path,"line_start":line_start,"line_end":line_end,
            "page":src.page,"sha256":src.sha256,"text":text,"extent":extent}),
        )
    }
    async fn theory_pattern_by_id(&self, id: &str) -> Result<Option<Value>> {
        let Ok(map) = self.read_theory_map().await else {
            return Ok(None);
        };
        Ok(map
            .get("patterns")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .find(|pattern| pattern.get("id").and_then(Value::as_str) == Some(id))
            .cloned())
    }
    async fn theory_pattern_source_by_id(&self, id: &str) -> Result<Option<Source>> {
        let Ok(map) = self.read_theory_map().await else {
            return Ok(None);
        };
        Ok(map
            .get("patterns")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|pattern| pattern.get("source_refs").and_then(Value::as_array))
            .flatten()
            .filter_map(theory_pattern_source)
            .find(|source| source.id == id))
    }
    fn source_text(
        &self,
        source: &Source,
        requested_start: Option<usize>,
        requested_end: Option<usize>,
    ) -> Result<(String, Option<usize>, Option<usize>, Value)> {
        let raw = Path::new(&source.path);
        if source.path.trim().is_empty() {
            bail!("source path is missing");
        }
        if raw.components().any(|c| matches!(c, Component::ParentDir)) {
            bail!("unsafe source path");
        }
        let path = if raw.is_absolute() {
            raw.to_path_buf()
        } else {
            self.repo_root.join(raw)
        };
        let canonical = fs::canonicalize(&path)
            .with_context(|| format!("source file not found: {}", source.path))?;
        if !canonical.starts_with(&self.repo_root) || sensitive_path(&canonical) {
            bail!("unsafe source path");
        }
        let metadata = fs::metadata(&canonical)?;
        if metadata.len() > 8 * 1024 * 1024 {
            bail!("source file exceeds safe read limit");
        }
        let bytes = fs::read(&canonical)?;
        if source.sha256.is_empty() {
            bail!("source hash is missing");
        }
        let hash = format!("{:x}", Sha256::digest(&bytes));
        if !hash.eq_ignore_ascii_case(&source.sha256) {
            bail!("source hash changed; excerpt is stale");
        }
        if canonical
            .extension()
            .and_then(|v| v.to_str())
            .is_some_and(|v| v.eq_ignore_ascii_case("pdf"))
        {
            if requested_start.is_some() || requested_end.is_some() {
                bail!("line ranges are unavailable for stored PDF excerpts");
            }
            let text = source
                .text
                .as_ref()
                .context("PDF source has no stored text excerpt")?;
            let truncated = text.len() > 200_000;
            let returned = bound_text(text.clone());
            return Ok((
                returned.clone(),
                None,
                None,
                json!({
                "kind":"stored_pdf_excerpt","range_complete":false,"truncated":truncated,
                "registered_line_start":source.line_start,"registered_line_end":source.line_end,
                "returned_chars":returned.chars().count(),"stored_chars":text.chars().count(),
                "notice":"Stored PDF text preview; original page/line coverage is not certified."}),
            ));
        }
        let content = String::from_utf8(bytes).context("source file is not UTF-8")?;
        let lines: Vec<&str> = content.lines().collect();
        let actual_line_count = lines.len();
        if let Some(expected) = source.line_count {
            if expected != actual_line_count {
                bail!("source line count changed; reference is stale");
            }
        }
        let allowed_start = source.line_start.unwrap_or(1).max(1);
        let allowed_end = source
            .line_end
            .or(source.line_count)
            .unwrap_or(actual_line_count);
        if allowed_start > allowed_end || allowed_end > actual_line_count {
            bail!("registered source line range is invalid");
        }
        let start = requested_start.unwrap_or(allowed_start);
        let end = requested_end.unwrap_or_else(|| start.saturating_add(199).min(allowed_end));
        if start < allowed_start || end > allowed_end || end < start {
            bail!("requested source range is outside the registered extent");
        }
        if end - start >= 200 {
            bail!("source line window must not exceed 200 lines");
        }
        let excerpt = lines
            .iter()
            .skip(start - 1)
            .take(end - start + 1)
            .copied()
            .collect::<Vec<_>>()
            .join("\n");
        let truncated = excerpt.len() > 200_000;
        let mut byte_end = excerpt.len().min(200_000);
        while !excerpt.is_char_boundary(byte_end) {
            byte_end -= 1;
        }
        let returned_end = start + excerpt[..byte_end].lines().count().saturating_sub(1);
        let range_complete = !truncated;
        let next_start = if truncated {
            returned_end
        } else {
            end.saturating_add(1)
        };
        Ok((
            bound_text(excerpt),
            Some(start),
            Some(returned_end),
            json!({
                "kind":"original_lines","range_complete":range_complete,"truncated":truncated,
                "requested_line_start":start,"requested_line_end":end,
                "registered_line_start":allowed_start,"registered_line_end":allowed_end,
                "registered_extent_complete":range_complete && start==allowed_start && end==allowed_end,
                "last_line_complete":range_complete,
                "next_request":if next_start<=allowed_end {json!({"source_id":source.id,"line_start":next_start,
                    "line_end":next_start.saturating_add(199).min(allowed_end)})}else{Value::Null}
            }),
        ))
    }
    pub async fn graph_slice(&self, id: Option<&str>, limit: usize) -> Result<Value> {
        if limit == 0 || limit > 1000 {
            bail!("limit must be between 1 and 1000");
        }
        let g = self.snapshot().await?;
        let selected: HashSet<String> = if let Some(id) = id {
            let _ = g
                .nodes
                .iter()
                .find(|n| n.id == id)
                .ok_or_else(|| anyhow!("node not found: {id}"))?;
            let mut s = HashSet::from([id.to_owned()]);
            for e in &g.edges {
                if e.src == id {
                    s.insert(e.dst.clone());
                }
                if e.dst == id {
                    s.insert(e.src.clone());
                }
            }
            s
        } else {
            g.nodes.iter().take(limit).map(|n| n.id.clone()).collect()
        };
        let nodes: Vec<_> = g
            .nodes
            .iter()
            .filter(|n| selected.contains(&n.id))
            .take(limit)
            .map(public_node)
            .collect();
        let selected_ids: HashSet<_> = g
            .nodes
            .iter()
            .filter(|n| selected.contains(&n.id))
            .take(limit)
            .map(|n| n.id.as_str())
            .collect();
        let edges: Vec<_> = g
            .edges
            .iter()
            .filter(|e| {
                selected_ids.contains(e.src.as_str()) && selected_ids.contains(e.dst.as_str())
            })
            .take(limit * 2)
            .cloned()
            .collect();
        Ok(
            json!({"meta":g.meta,"nodes":nodes,"edges":edges,"sources":g.sources.iter().take(limit).cloned().collect::<Vec<_>>()}),
        )
    }
    async fn query_embedding(&self, query: &str) -> Result<Vec<f32>> {
        if let Some(v) = self
            .embedding_cache
            .read()
            .map_err(|_| anyhow!("embedding cache poisoned"))?
            .get(query)
            .cloned()
        {
            return Ok(v);
        }
        let key = api_key()?;
        let client = reqwest::Client::new();
        let response: Value = client
            .post("https://api.openai.com/v1/embeddings")
            .bearer_auth(key)
            .json(&json!({"model":"text-embedding-3-small","input":query,"dimensions":256}))
            .send()
            .await?
            .error_for_status()?
            .json()
            .await?;
        let values = response["data"][0]["embedding"]
            .as_array()
            .context("embedding response missing vector")?;
        let vector: Vec<f32> = values
            .iter()
            .map(|v| v.as_f64().unwrap_or(0.0) as f32)
            .collect();
        if vector.len() != 256 {
            bail!("embedding response has unexpected dimensions");
        }
        let mut cache = self
            .embedding_cache
            .write()
            .map_err(|_| anyhow!("embedding cache poisoned"))?;
        if cache.len() >= 128 {
            if let Some(key) = cache.keys().next().cloned() {
                cache.remove(&key);
            }
        }
        cache.insert(query.to_owned(), vector.clone());
        Ok(vector)
    }
}

fn parse_graph(bytes: &[u8]) -> Result<KnowledgeGraph> {
    if bytes.len() > 256 * 1024 * 1024 {
        bail!("knowledge graph exceeds 256 MiB safe read limit");
    }
    let g: KnowledgeGraph =
        serde_json::from_slice(bytes).context("invalid knowledge graph JSON")?;
    let mut ids = HashSet::new();
    for n in &g.nodes {
        if n.id.is_empty() || !ids.insert(n.id.as_str()) {
            bail!("graph has an empty or duplicate node id");
        }
    }
    for e in &g.edges {
        if !ids.contains(e.src.as_str()) || !ids.contains(e.dst.as_str()) {
            bail!("edge {} references missing endpoint", e.id);
        }
    }
    Ok(g)
}

fn file_stamp(path: &Path) -> Option<(SystemTime, u64)> {
    let metadata = fs::metadata(path).ok()?;
    Some((metadata.modified().ok()?, metadata.len()))
}

fn load_merged_graph(
    main_path: &Path,
    context_path: &Path,
    translations_path: &Path,
    quality_path: &Path,
    corpus_path: &Path,
    embeddings_cache_path: &Path,
    orientation_path: &Path,
) -> Result<KnowledgeGraph> {
    let bytes = fs::read(main_path)
        .with_context(|| format!("knowledge graph unavailable: {}", main_path.display()))?;
    let mut main = parse_graph(&bytes)?;
    let overlay: KnowledgeGraph = match fs::read(context_path) {
        Ok(overlay_bytes) => {
            if overlay_bytes.len() > 256 * 1024 * 1024 {
                bail!("context overlay exceeds 256 MiB safe read limit");
            }
            serde_json::from_slice(&overlay_bytes).context("invalid context overlay JSON")?
        }
        Err(_) => KnowledgeGraph::default(),
    };
    let main_claims: HashMap<&str, &Node> = main
        .nodes
        .iter()
        .filter(|n| n.kind == "claim")
        .map(|n| (n.id.as_str(), n))
        .collect();
    let eligible_priority_claims: HashSet<String> = overlay
        .nodes
        .iter()
        .filter(|n| n.kind == "claim")
        .filter_map(|n| {
            let annotation = n.attrs.get("priority_annotation")?;
            let priority_enriched =
                annotation.get("enrichment_state").and_then(Value::as_str) == Some("enriched");
            let main_claim = main_claims.get(n.id.as_str()).copied();
            let main_state = main_claim
                .and_then(|claim| claim.attrs.get("enrichment_state"))
                .and_then(Value::as_str);
            let replaces_current = annotation
                .get("replaces_semantic_sha256")
                .and_then(Value::as_str)
                .zip(main_claim.and_then(semantic_translation_fingerprint))
                .is_some_and(|(expected, actual)| expected == actual);
            (priority_enriched
                && (matches!(main_state, Some("pending") | Some("error"))
                    || (main_state == Some("enriched") && replaces_current)))
                .then(|| n.id.clone())
        })
        .collect();
    let mut node_index: HashMap<String, usize> = main
        .nodes
        .iter()
        .enumerate()
        .map(|(i, n)| (n.id.clone(), i))
        .collect();
    for node in overlay.nodes {
        if node
            .attrs
            .get("priority_owner")
            .and_then(Value::as_str)
            .is_some_and(|owner| !eligible_priority_claims.contains(owner))
        {
            continue;
        }
        if let Some(index) = node_index.get(&node.id).copied() {
            if let (Some(base), Some(extra)) = (
                main.nodes[index].attrs.as_object_mut(),
                node.attrs.as_object(),
            ) {
                let annotation = extra.get("priority_annotation").and_then(Value::as_object);
                let overlay_enriched = extra.get("enrichment_state").and_then(Value::as_str)
                    == Some("enriched")
                    || annotation
                        .and_then(|a| a.get("enrichment_state"))
                        .and_then(Value::as_str)
                        == Some("enriched");
                let may_promote = overlay_enriched && eligible_priority_claims.contains(&node.id);
                for (key, value) in extra {
                    if key == "ledger" || key == "priority_annotation" {
                        continue;
                    }
                    if matches!(key.as_str(), "semantic" | "provenance" | "enrichment_state")
                        && !may_promote
                    {
                        continue;
                    }
                    if key != "semantic" || may_promote {
                        base.insert(key.clone(), value.clone());
                    }
                }
                if may_promote {
                    if let Some(annotation) = annotation {
                        for key in [
                            "semantic",
                            "provenance",
                            "enrichment_state",
                            "source_refs",
                            "source_coverage",
                        ] {
                            if let Some(value) = annotation.get(key) {
                                base.insert(key.to_owned(), value.clone());
                            }
                        }
                    }
                }
                base.remove("priority_annotation");
            }
        } else {
            node_index.insert(node.id.clone(), main.nodes.len());
            main.nodes.push(node);
        }
    }
    let mut edge_keys: HashSet<(String, String, String)> = main
        .edges
        .iter()
        .map(|e| (e.src.clone(), e.dst.clone(), e.kind.clone()))
        .collect();
    for edge in overlay.edges {
        if edge
            .attrs
            .get("priority_owner")
            .and_then(Value::as_str)
            .is_some_and(|owner| !eligible_priority_claims.contains(owner))
        {
            continue;
        }
        // Keep the usable overlay subgraph; dangling overlay links cannot be represented safely.
        if !node_index.contains_key(&edge.src) || !node_index.contains_key(&edge.dst) {
            continue;
        }
        if edge_keys.insert((edge.src.clone(), edge.dst.clone(), edge.kind.clone())) {
            main.edges.push(edge);
        }
    }
    let source_ids: HashSet<String> = main.sources.iter().map(|s| s.id.clone()).collect();
    main.sources.extend(
        overlay
            .sources
            .into_iter()
            .filter(|s| !source_ids.contains(&s.id)),
    );

    // Corpus records can refer across the main/context/corpus graphs. Merge
    // nodes first and admit only edges whose endpoints are now available. When
    // an ID already exists, main/context remains authoritative; corpus attrs
    // only fill missing fields.
    let mut corpus_added = 0usize;
    let mut corpus_merged = 0usize;
    let mut corpus_edges_added = 0usize;
    let mut corpus_edges_skipped = 0usize;
    let mut corpus_sources_added = 0usize;
    let corpus_bytes = match fs::read(corpus_path) {
        Ok(bytes) => Some(bytes),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
        Err(error) => {
            return Err(error)
                .with_context(|| format!("cannot read corpus overlay: {}", corpus_path.display()))
        }
    };
    if let Some(corpus_bytes) = corpus_bytes {
        if corpus_bytes.len() > 256 * 1024 * 1024 {
            bail!("corpus overlay exceeds 256 MiB safe read limit");
        }
        let corpus: KnowledgeGraph =
            serde_json::from_slice(&corpus_bytes).context("invalid corpus overlay JSON")?;
        let mut corpus_ids = HashSet::new();
        for node in corpus.nodes {
            if node.id.is_empty() || !corpus_ids.insert(node.id.clone()) {
                bail!("corpus overlay has an empty or duplicate node id");
            }
            if let Some(index) = node_index.get(&node.id).copied() {
                corpus_merged += 1;
                if let Some(extras) = node.attrs.as_object() {
                    if !main.nodes[index].attrs.is_object() {
                        main.nodes[index].attrs = json!({});
                    }
                    let target = main.nodes[index].attrs.as_object_mut().unwrap();
                    for (key, value) in extras {
                        if matches!(key.as_str(), "semantic" | "ledger") {
                            continue;
                        }
                        target.entry(key.clone()).or_insert_with(|| value.clone());
                    }
                }
            } else {
                node_index.insert(node.id.clone(), main.nodes.len());
                main.nodes.push(node);
                corpus_added += 1;
            }
        }
        let mut edge_keys: HashSet<(String, String, String)> = main
            .edges
            .iter()
            .map(|edge| (edge.src.clone(), edge.dst.clone(), edge.kind.clone()))
            .collect();
        for edge in corpus.edges {
            if !node_index.contains_key(&edge.src) || !node_index.contains_key(&edge.dst) {
                corpus_edges_skipped += 1;
                continue;
            }
            if edge_keys.insert((edge.src.clone(), edge.dst.clone(), edge.kind.clone())) {
                main.edges.push(edge);
                corpus_edges_added += 1;
            }
        }
        let mut source_ids: HashSet<String> = main
            .sources
            .iter()
            .map(|source| source.id.clone())
            .collect();
        for source in corpus.sources {
            if source_ids.insert(source.id.clone()) {
                main.sources.push(source);
                corpus_sources_added += 1;
            }
        }
    }
    if !main.meta.is_object() {
        main.meta = json!({});
    }
    main.meta.as_object_mut().unwrap().insert(
        "corpus_merge".to_owned(),
        json!({
            "nodes_added": corpus_added,
            "nodes_merged": corpus_merged,
            "edges_added": corpus_edges_added,
            "edges_skipped_unknown_endpoint": corpus_edges_skipped,
            "sources_added": corpus_sources_added
        }),
    );

    // Translations are a presentation-only overlay. They are accepted only when
    // their fingerprint matches the exact German title/overview/interpretation
    // currently attached to the authoritative semantic annotation.
    if let Ok(translation_bytes) = fs::read(translations_path) {
        if translation_bytes.len() > 256 * 1024 * 1024 {
            bail!("translations overlay exceeds 256 MiB safe read limit");
        }
        let translations: KnowledgeGraph = serde_json::from_slice(&translation_bytes)
            .context("invalid translations overlay JSON")?;
        let main_nodes: HashMap<String, usize> = main
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, node)| node.kind == "claim")
            .map(|(index, node)| (node.id.clone(), index))
            .collect();
        for translated in translations
            .nodes
            .iter()
            .filter(|node| node.kind == "claim")
        {
            let Some(index) = main_nodes.get(&translated.id).copied() else {
                continue;
            };
            let Some(candidate) = translated.attrs.get("translations") else {
                continue;
            };
            let Some(expected) = candidate
                .get("source_semantic_sha256")
                .and_then(Value::as_str)
            else {
                continue;
            };
            let Some(actual) = semantic_translation_fingerprint(&main.nodes[index]) else {
                continue;
            };
            if expected == actual {
                if let Some(attrs) = main.nodes[index].attrs.as_object_mut() {
                    attrs.insert("translations".to_owned(), candidate.clone());
                }
            }
        }
    }
    if let Ok(quality_bytes) = fs::read(quality_path) {
        if quality_bytes.len() > 64 * 1024 * 1024 {
            bail!("quality overlay exceeds 64 MiB safe read limit");
        }
        let quality: KnowledgeGraph =
            serde_json::from_slice(&quality_bytes).context("invalid quality overlay JSON")?;
        let main_nodes: HashMap<String, usize> = main
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, node)| node.kind == "claim")
            .map(|(index, node)| (node.id.clone(), index))
            .collect();
        for issue_node in quality.nodes.iter().filter(|node| node.kind == "claim") {
            let Some(annotation) = issue_node
                .attrs
                .get("extraction_quality")
                .and_then(Value::as_object)
            else {
                continue;
            };
            if let Some(index) = main_nodes.get(&issue_node.id).copied() {
                let annotation_hash = annotation.get("annotation_sha256").and_then(Value::as_str);
                let summary_hash = annotation
                    .get("source_semantic_sha256")
                    .and_then(Value::as_str);
                let (expected, actual) = if let Some(expected) = annotation_hash {
                    (
                        Some(expected),
                        semantic_annotation_fingerprint(&main.nodes[index]),
                    )
                } else {
                    (
                        summary_hash,
                        semantic_translation_fingerprint(&main.nodes[index]),
                    )
                };
                let Some(expected) = expected else {
                    continue;
                };
                let mut value = Value::Object(annotation.clone());
                if actual.as_deref() != Some(expected) {
                    if let Some(fields) = value.as_object_mut() {
                        fields.insert("state".to_owned(), json!("outdated_review"));
                        if let Some(issues) = fields.remove("issues") {
                            fields.insert("stale_issues".to_owned(), issues);
                        }
                    }
                }
                if let Some(attrs) = main.nodes[index].attrs.as_object_mut() {
                    attrs.insert("extraction_quality".to_owned(), value);
                }
            }
        }
        let mut source_ids: HashSet<String> = main
            .sources
            .iter()
            .map(|source| source.id.clone())
            .collect();
        let mut quality_sources_added = 0usize;
        for source in quality.sources {
            if source_ids.insert(source.id.clone()) {
                main.sources.push(source);
                quality_sources_added += 1;
            }
        }
        if !main.meta.is_object() {
            main.meta = json!({});
        }
        main.meta.as_object_mut().unwrap().insert(
            "quality_merge".to_owned(),
            json!({"sources_added":quality_sources_added}),
        );
    }

    let mut orientation_fields_added = 0usize;
    if let Ok(bytes) = fs::read(orientation_path) {
        if bytes.len() > 64 * 1024 * 1024 {
            bail!("orientation overlay exceeds 64 MiB safe read limit");
        }
        let orientation: Value =
            serde_json::from_slice(&bytes).context("invalid orientation overlay JSON")?;
        let node_index: HashMap<String, usize> = main
            .nodes
            .iter()
            .enumerate()
            .map(|(index, node)| (node.id.clone(), index))
            .collect();
        for record in ["structures", "goals"].into_iter().flat_map(|group| {
            orientation
                .get(group)
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
        }) {
            let Some(id) = record.get("id").and_then(Value::as_str) else {
                continue;
            };
            let Some(index) = node_index.get(id).copied() else {
                continue;
            };
            if !main.nodes[index].attrs.is_object() {
                main.nodes[index].attrs = json!({});
            }
            let attrs = main.nodes[index].attrs.as_object_mut().unwrap();
            for field in ["label_en", "summary_en", "required_steps_en"] {
                if attrs.get(field).is_some_and(|value| !value.is_null()) {
                    continue;
                }
                let Some(value) = record.get(field).filter(|value| !value.is_null()) else {
                    continue;
                };
                attrs.insert(field.to_owned(), value.clone());
                orientation_fields_added += 1;
            }
        }
    }
    if !main.meta.is_object() {
        main.meta = json!({});
    }
    main.meta.as_object_mut().unwrap().insert(
        "orientation_overlay".to_owned(),
        json!({"english_fields_added":orientation_fields_added,"path":orientation_path.file_name().and_then(|name|name.to_str()).unwrap_or("orientation.json")}),
    );

    // A derivation that lists its own target claim as an upstream claim is a
    // self-reference, not an evidentiary dependency. Keep it visible but typed
    // so clients cannot mistake it for an independent proof path.
    let mut derivation_owners: HashMap<String, HashSet<String>> = HashMap::new();
    for edge in &main.edges {
        if edge.kind == "has_derivation" {
            derivation_owners
                .entry(edge.dst.clone())
                .or_default()
                .insert(edge.src.clone());
        }
    }
    // Older records encode the owner only in derivation:<claim-id>:<local-id>.
    for node in main.nodes.iter().filter(|node| node.kind == "derivation") {
        if let Some((owner, _)) = node
            .id
            .strip_prefix("derivation:")
            .and_then(|id| id.rsplit_once(':'))
        {
            derivation_owners
                .entry(node.id.clone())
                .or_insert_with(|| HashSet::from([format!("claim:{owner}")]));
        }
    }
    for edge in &mut main.edges {
        if edge.kind == "uses_claim_under_scope"
            && derivation_owners
                .get(&edge.dst)
                .is_some_and(|claims| claims.contains(&edge.src))
        {
            edge.kind = "self_reference".to_owned();
            if !edge.attrs.is_object() {
                edge.attrs = json!({});
            }
            edge.attrs.as_object_mut().unwrap().insert(
                "epistemic_role".to_owned(),
                json!("non_evidentiary_self_reference"),
            );
        }
    }
    let ids: HashSet<&str> = main.nodes.iter().map(|n| n.id.as_str()).collect();
    for edge in &main.edges {
        if !ids.contains(edge.src.as_str()) || !ids.contains(edge.dst.as_str()) {
            bail!("merged graph edge {} references missing endpoint", edge.id);
        }
    }
    crate::papers::merge(&mut main, &main_path.with_file_name("papers.json"))?;
    hydrate_embeddings_from_cache(&mut main, embeddings_cache_path);
    Ok(main)
}

fn hydrate_embeddings_from_cache(graph: &mut KnowledgeGraph, cache_path: &Path) {
    let mut cache_available = false;
    let mut cache_entries = 0usize;
    let mut valid_cache_entries = 0usize;
    let mut matched_nodes = 0usize;
    let mut existing_vectors = 0usize;
    let cache = match fs::read(cache_path) {
        Ok(bytes) if bytes.len() <= 256 * 1024 * 1024 => {
            let parsed = serde_json::from_slice::<HashMap<String, Vec<f32>>>(&bytes).ok();
            cache_available = parsed.is_some();
            parsed
        }
        _ => None,
    };
    if let Some(cache) = cache {
        cache_entries = cache.len();
        valid_cache_entries = cache
            .values()
            .filter(|vector| valid_embedding_vector(vector))
            .count();
        for node in &mut graph.nodes {
            if embedding(node).is_some() {
                existing_vectors += 1;
                continue;
            }
            let key = embedding_content_key(&node.label, &node.summary);
            let Some(vector) = cache
                .get(&key)
                .filter(|vector| valid_embedding_vector(vector))
            else {
                continue;
            };
            let Some(attrs) = node.attrs.as_object_mut() else {
                continue;
            };
            attrs.insert("embedding".to_owned(), json!(vector));
            attrs.insert(
                "embedding_model".to_owned(),
                json!("text-embedding-3-small"),
            );
            matched_nodes += 1;
        }
    }
    let vector_nodes = graph
        .nodes
        .iter()
        .filter(|node| embedding(node).is_some())
        .count();
    if !graph.meta.is_object() {
        graph.meta = json!({});
    }
    graph.meta.as_object_mut().unwrap().insert(
        "semantic_search".to_owned(),
        json!({
            "available":vector_nodes > 0,
            "model":"text-embedding-3-small",
            "cache_available":cache_available,
            "cache_entries":cache_entries,
            "valid_cache_entries":valid_cache_entries,
            "hydrated_nodes":matched_nodes,
            "existing_graph_vectors":existing_vectors,
            "vector_nodes":vector_nodes,
            "graph_nodes":graph.nodes.len(),
            "coverage":if graph.nodes.is_empty(){0.0}else{vector_nodes as f64 / graph.nodes.len() as f64},
            "coverage_scope":"Exact SHA-256 matches of stored node label + newline + summary. Cached vectors cover selected records only; coverage is partial.",
            "query_embedding":"On-demand embedding for explicit search requests; falls back to lexical ranking if unavailable."
        }),
    );
}

fn valid_embedding_vector(vector: &[f32]) -> bool {
    vector.len() == 256 && vector.iter().all(|value| value.is_finite())
}

fn embedding_content_key(label: &str, summary: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(format!("{label}\n{summary}").as_bytes());
    format!("{:x}", hasher.finalize())
}

fn semantic_translation_fingerprint(node: &Node) -> Option<String> {
    let semantic = node.attrs.get("semantic")?;
    let title = semantic.get("title_de")?.as_str()?;
    let overview = semantic.get("overview_de")?.as_str()?;
    let interpretation = semantic.get("interpretation_de")?.as_str()?;
    let mut hasher = Sha256::new();
    hasher.update(format!("{title}\n{overview}\n{interpretation}").as_bytes());
    Some(format!("{:x}", hasher.finalize()))
}
fn semantic_annotation_fingerprint(node: &Node) -> Option<String> {
    let semantic = node.attrs.get("semantic")?;
    let bytes = serde_json::to_vec(semantic).ok()?;
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    Some(format!("{:x}", hasher.finalize()))
}

fn bounded_char_limit(limit: usize) -> Result<usize> {
    if !(1000..=30_000).contains(&limit) {
        bail!("max_chars must be between 1000 and 30000");
    }
    Ok(limit)
}

fn json_char_count(value: &Value) -> usize {
    serde_json::to_string(value)
        .unwrap_or_default()
        .chars()
        .count()
}

fn compact_claim_status(value: &Value) -> Value {
    if let Some(text) = value.as_str() {
        return json!(clip(text, 180));
    }
    let Some(object) = value.as_object() else {
        return Value::Null;
    };
    let mut compact = serde_json::Map::new();
    for key in ["category", "label", "scope", "authority", "notice"] {
        if let Some(value) = object.get(key) {
            compact.insert(
                key.to_owned(),
                if value.is_string() {
                    json!(clip(&value_text(Some(value)), 140))
                } else {
                    value.clone()
                },
            );
        }
    }
    Value::Object(compact)
}

fn fit_compare_budget(mut value: Value, max_chars: usize) -> Value {
    if json_char_count(&value) <= max_chars {
        return value;
    }
    if let Some(links) = value
        .pointer_mut("/shared_explicit/recorded_graph_links")
        .and_then(Value::as_array_mut)
    {
        // When compacting, remove generic topic links before explicitly typed
        // shared-seed/structure links.
        links.sort_by_key(|link| {
            !link
                .get("edge_type")
                .and_then(Value::as_str)
                .is_some_and(|kind| kind.contains("seed") || kind.contains("shared"))
        });
    }
    for path in [
        "/shared_source_ids",
        "/shared_explicit/recorded_graph_nodes",
        "/shared_explicit/recorded_graph_links",
        "/shared_explicit/seeds",
        "/shared_explicit/upstream_claim_ids",
    ] {
        while json_char_count(&value) > max_chars {
            let Some(array) = value.pointer_mut(path).and_then(Value::as_array_mut) else {
                break;
            };
            if array.len() <= 1 {
                break;
            }
            array.pop();
        }
    }
    if json_char_count(&value) <= max_chars {
        value["truncated"] = json!(true);
        return value;
    }

    // At this point comparison detail is larger than the caller's budget.
    // Preserve every compared identity and its recorded status while shrinking
    // per-item detail; shared edges remain typed and distinguishable.
    let items: Vec<Value> = value
        .get("items")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .map(compact_compare_item)
        .collect();
    let shared = value.get("shared_explicit").unwrap_or(&Value::Null);
    let compact_nodes: Vec<Value> = shared
        .get("recorded_graph_nodes")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .take(4)
        .map(|entry| json!({"id":entry.get("id"),"type":entry.get("type")}))
        .collect();
    let compact_links: Vec<Value> = shared
        .get("recorded_graph_links")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .take(4)
        .map(|entry| json!({"edge_type":entry.get("edge_type"),"direction":entry.get("direction"),"node_id":entry.get("node_id")}))
        .collect();
    let compact = json!({
        "ids":value.get("ids").cloned().unwrap_or(Value::Null),
        "items":items.clone(),
        "shared_explicit":{
            "upstream_claim_ids":shared.get("upstream_claim_ids").and_then(Value::as_array).into_iter().flatten().take(4).cloned().collect::<Vec<_>>(),
            "seeds":shared.get("seeds").and_then(Value::as_array).into_iter().flatten().take(4).cloned().collect::<Vec<_>>(),
            "recorded_graph_nodes":compact_nodes,
            "recorded_graph_links":compact_links
        },
        "shared_source_ids":value.get("shared_source_ids").and_then(Value::as_array).into_iter().flatten().take(4).cloned().collect::<Vec<_>>(),
        "status_comparison":value.get("status_comparison").cloned().unwrap_or(Value::Null),
        "divergent_assumptions":clip(&value_text(value.get("divergent_assumptions")),100),
        "notice":"Budgeted comparison: identities, recorded statuses, and shared stored links are retained; common links do not establish proof.",
        "max_chars":max_chars,
        "truncated":true
    });
    if json_char_count(&compact) <= max_chars {
        compact
    } else {
        // The stated API minimum (1000 chars) is enough for all IDs and status
        // categories, even when topic labels and link detail have to be omitted.
        json!({
            "ids":value.get("ids").cloned().unwrap_or(Value::Null),
            "items":items.iter().map(|item| {
                let mut compact = serde_json::Map::new();
                for key in ["id", "type", "label", "assessment", "ledger_status", "verdict", "verdict_raw", "enrichment_state", "scope", "assumptions", "assumptions_total", "assumptions_omitted", "recorded_conclusions", "original_source_scope_correction", "topic_contents"] {
                    if let Some(value) = item.get(key).filter(|value| !value.is_null()) {
                        compact.insert(key.to_string(), value.clone());
                    }
                }
                Value::Object(compact)
            }).collect::<Vec<_>>(),
            "shared_explicit":{
                "upstream_claim_ids":shared.get("upstream_claim_ids").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                "seeds":shared.get("seeds").and_then(Value::as_array).into_iter().flatten().take(2).cloned().collect::<Vec<_>>(),
                "recorded_graph_nodes":compact_nodes.iter().take(2).cloned().collect::<Vec<_>>(),
                "recorded_graph_links":compact_links.iter().take(2).cloned().collect::<Vec<_>>()
            },
            "status_comparison":value.get("status_comparison").cloned().unwrap_or_else(|| json!([])),
            "notice":"Budgeted comparison. Shared stored links do not establish proof.",
            "max_chars":max_chars,
            "truncated":true
        })
    }
}

fn compact_compare_item(item: &Value) -> Value {
    let mut compact = serde_json::Map::new();
    for key in [
        "id",
        "type",
        "assessment",
        "ledger_status",
        "verdict",
        "verdict_raw",
        "enrichment_state",
        "assumptions_total",
        "assumptions_omitted",
        "recorded_conclusions",
        "original_source_scope_correction",
    ] {
        if let Some(value) = item.get(key).filter(|value| !value.is_null()) {
            compact.insert(key.to_string(), value.clone());
        }
    }
    if let Some(label) = item.get("label").and_then(Value::as_str) {
        compact.insert("label".to_string(), json!(clip(label, 100)));
    }
    if let Some(scope) = item.get("scope").and_then(Value::as_str) {
        compact.insert("scope".to_string(), json!(clip(scope, 240)));
    }
    if let Some(summary) = item.get("summary").and_then(Value::as_str) {
        compact.insert("summary".to_string(), json!(clip(summary, 600)));
    }
    for key in ["upstream_claim_ids", "seeds"] {
        if let Some(values) = item.get(key).and_then(Value::as_array) {
            compact.insert(
                key.to_string(),
                json!(values.iter().take(8).cloned().collect::<Vec<_>>()),
            );
        }
    }
    if let Some(inputs) = item.get("assumptions").and_then(Value::as_array) {
        let compact_inputs: Vec<Value> = inputs
            .iter()
            .take(8)
            .map(|input| {
                json!({
                    "name":input.get("name"),
                    "origin":input.get("origin"),
                    "value":input.get("value").and_then(Value::as_str).map(|s|clip(s,100)),
                    "used_for":input.get("used_for").and_then(Value::as_str).map(|s|clip(s,120)),
                    "impact":input.get("impact").cloned().unwrap_or(Value::Null),
                    "source_ids":input.get("source_ids").and_then(Value::as_array).into_iter().flatten().take(1).cloned().collect::<Vec<_>>()
                })
            })
            .collect();
        let assumptions_omitted = (item
            .get("assumptions_total")
            .and_then(Value::as_u64)
            .unwrap_or(inputs.len() as u64) as usize)
            .saturating_sub(compact_inputs.len());
        compact.insert("assumptions".to_string(), json!(compact_inputs));
        compact.insert(
            "assumptions_omitted".to_string(),
            json!(assumptions_omitted),
        );
        if let Some(value) = item.get("assumptions_total") {
            compact.insert("assumptions_total".to_string(), value.clone());
        }
    }
    if let Some(derivations) = item.get("derivations").and_then(Value::as_array) {
        let compact_derivations: Vec<Value> = derivations
            .iter()
            .take(3)
            .map(|derivation| {
                json!({
                    "id":derivation.get("id"),
                    "label":derivation.get("label").and_then(Value::as_str).map(|s|clip(s,100)),
                    "premise_ids":derivation.get("premise_ids").and_then(Value::as_array).into_iter().flatten().take(6).cloned().collect::<Vec<_>>(),
                    "upstream_claim_ids":derivation.get("upstream_claim_ids").and_then(Value::as_array).into_iter().flatten().take(6).cloned().collect::<Vec<_>>(),
                    "result_ids":derivation.get("result_ids").and_then(Value::as_array).into_iter().flatten().take(4).cloned().collect::<Vec<_>>(),
                    "argument":derivation.get("argument").and_then(Value::as_str).map(|s|clip(s,240)),
                    "scope":derivation.get("scope").and_then(Value::as_str).map(|s|clip(s,120))
                })
            })
            .collect();
        compact.insert("derivations".to_string(), json!(compact_derivations));
    }
    if let Some(source_ids) = item.get("source_ids").and_then(Value::as_array) {
        compact.insert(
            "source_ids".to_string(),
            json!(source_ids.iter().take(2).cloned().collect::<Vec<_>>()),
        );
    }
    if let Some(contents) = item.get("topic_contents").filter(|value| !value.is_null()) {
        let mut summary = serde_json::Map::new();
        for key in ["claim_count", "assessment_counts"] {
            if let Some(value) = contents.get(key).filter(|value| !value.is_null()) {
                summary.insert(key.to_string(), value.clone());
            }
        }
        if !summary.is_empty() {
            compact.insert("topic_contents".to_string(), Value::Object(summary));
        }
    }
    Value::Object(compact)
}

fn compact_brief_array_fields(items: Option<&mut Value>, fields: &[(&str, usize)]) {
    if let Some(items) = items.and_then(Value::as_array_mut) {
        for item in items {
            for (field, limit) in fields {
                if let Some(text) = item.get(*field).and_then(Value::as_str).map(str::to_owned) {
                    item[*field] = json!(clip(&text, *limit));
                }
            }
            if let Some(sources) = item.get_mut("source_ids").and_then(Value::as_array_mut) {
                sources.truncate(1);
            }
        }
    }
}

fn compact_research_brief(value: &mut Value) {
    let truncate = |value: &mut Value, path: &str, max: usize| {
        if let Some(items) = value.pointer_mut(path).and_then(Value::as_array_mut) {
            items.truncate(max);
        }
    };
    truncate(value, "/what_we_have", 1);
    let what_we_have_statement_limit = if value.get("kind").and_then(Value::as_str) == Some("claim")
    {
        1_200
    } else {
        110
    };
    compact_brief_array_fields(
        value.pointer_mut("/what_we_have"),
        &[("scope", 140), ("statement", what_we_have_statement_limit)],
    );
    truncate(value, "/scoped_results", 5);
    compact_brief_array_fields(
        value.pointer_mut("/scoped_results"),
        &[
            ("claim_label", 80),
            ("assessment_scope", 80),
            ("statement", 150),
            ("scope", 70),
        ],
    );
    truncate(value, "/next_testable_requirements", 3);
    compact_brief_array_fields(
        value.pointer_mut("/next_testable_requirements"),
        &[
            ("claim_label", 60),
            ("requirement", 100),
            ("why_checks_insufficient", 65),
        ],
    );
    truncate(value, "/attempts/items", 4);
    if let Some(attempts) = value.get_mut("attempts").and_then(Value::as_object_mut) {
        attempts.insert("limit".into(), json!(4));
        if let Some(items) = attempts.get_mut("items").and_then(Value::as_array_mut) {
            for item in items {
                for (field, limit) in [("label", 80), ("question_or_scope", 120)] {
                    if let Some(text) = item.get(field).and_then(Value::as_str).map(str::to_owned) {
                        item[field] = json!(clip(&text, limit));
                    }
                }
                if let Some(sources) = item.get_mut("source_ids").and_then(Value::as_array_mut) {
                    sources.truncate(1);
                }
                if let Some(conclusions) = item
                    .get_mut("recorded_conclusions")
                    .and_then(Value::as_object_mut)
                {
                    for values in conclusions.values_mut().filter_map(Value::as_array_mut) {
                        values.truncate(1);
                        if let Some(text) =
                            values.first().and_then(Value::as_str).map(str::to_owned)
                        {
                            values[0] = json!(clip(&text, 95));
                        }
                    }
                }
            }
        }
    }
    truncate(value, "/hidden_input_risks", 2);
    compact_brief_array_fields(
        value.pointer_mut("/hidden_input_risks"),
        &[("value", 45), ("used_for", 50)],
    );
    truncate(value, "/known_negatives", 2);
    compact_brief_array_fields(
        value.pointer_mut("/known_negatives"),
        &[("label", 70), ("scope", 90)],
    );
    truncate(value, "/recorded_scoped_rule_outs", 1);
    compact_brief_array_fields(
        value.pointer_mut("/recorded_scoped_rule_outs"),
        &[("label", 70), ("scope", 80)],
    );
    if let Some(items) = value
        .pointer_mut("/recorded_scoped_rule_outs")
        .and_then(Value::as_array_mut)
    {
        for item in items {
            item.as_object_mut().map(|object| object.remove("category"));
            if let Some(rules) = item.get_mut("rules_out").and_then(Value::as_array_mut) {
                rules.truncate(1);
                if let Some(text) = rules.first().and_then(Value::as_str).map(str::to_owned) {
                    rules[0] = json!(clip(&text, 100));
                }
            }
        }
    }
    truncate(value, "/patterns", 1);
    if let Some(patterns) = value.get_mut("patterns").and_then(Value::as_array_mut) {
        for pattern in patterns {
            for (field, limit) in [("label", 85), ("summary", 95)] {
                if let Some(text) = pattern
                    .get(field)
                    .and_then(Value::as_str)
                    .map(str::to_owned)
                {
                    pattern[field] = json!(clip(&text, limit));
                }
            }
            for field in ["source_ids", "related_ids", "conditions", "counterevidence"] {
                if let Some(items) = pattern.get_mut(field).and_then(Value::as_array_mut) {
                    items.truncate(if field == "related_ids" { 3 } else { 1 });
                    if matches!(field, "conditions" | "counterevidence") {
                        for item in items {
                            if let Some(object) = item.as_object_mut() {
                                for language in ["label_de", "label_en"] {
                                    if let Some(text) = object
                                        .get(language)
                                        .and_then(Value::as_str)
                                        .map(str::to_owned)
                                    {
                                        object.insert(language.into(), json!(clip(&text, 45)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if let Some(orientation) = value
        .get_mut("topic_orientation")
        .and_then(Value::as_object_mut)
    {
        for entry in orientation.values_mut() {
            if let Some(text) = entry.as_str().map(str::to_owned) {
                *entry = json!(clip(&text, 80));
            }
        }
    }
    if let Some(record) = value.get_mut("attempt_record") {
        *record = compact_root_attempt_record(record, 170, 8);
    }
    if let Some(records) = value
        .get_mut("typed_related_records")
        .and_then(Value::as_array_mut)
    {
        records.truncate(5);
        for record in records {
            for field in [
                "positive_results",
                "resolved_results",
                "assumptions",
                "kills",
                "leaves_open",
                "not_claimed",
                "source_refs",
            ] {
                if let Some(items) = record.get_mut(field).and_then(Value::as_array_mut) {
                    items.truncate(1);
                    for item in items.iter_mut() {
                        if let Some(text) = item.as_str().map(str::to_owned) {
                            *item = json!(clip(&text, 100));
                        }
                    }
                }
            }
        }
    }
    let mut source_ids_omitted = 0usize;
    if let Some(sources) = value.get_mut("source_ids").and_then(Value::as_array_mut) {
        let original_len = sources.len();
        sources.truncate(4);
        source_ids_omitted = original_len.saturating_sub(4);
    }
    let previous = value
        .get("source_ids_omitted")
        .and_then(Value::as_u64)
        .unwrap_or(0) as usize;
    value["source_ids_omitted"] = json!(previous + source_ids_omitted);
}

fn fit_char_budget(mut value: Value, max_chars: usize) -> Value {
    let buckets = [
        "/scoped_results",
        "/typed_related_records",
        "/what_we_have",
        "/what_we_seek",
        "/attempts/items",
        "/known_negatives",
        "/hidden_input_risks",
        "/syntactic_input_sample",
        "/patterns",
        "/recorded_scoped_rule_outs",
        "/source_ids",
        "/next_testable_requirements",
        "/items",
    ];
    let mut truncated = HashSet::new();
    while json_char_count(&value) > max_chars {
        let candidate = buckets
            .iter()
            .filter(|path| {
                !(matches!(
                    value.get("kind").and_then(Value::as_str),
                    Some("attempt" | "contract")
                ) && **path == "/typed_related_records")
            })
            .filter_map(|path| {
                value
                    .pointer(path)
                    .and_then(Value::as_array)
                    .filter(|array| array.len() > 1)
                    .map(|array| (*path, array.last().map(json_char_count).unwrap_or(0)))
            })
            .max_by_key(|(_, size)| *size);
        let Some((path, _)) = candidate else { break };
        if let Some(array) = value.pointer_mut(path).and_then(Value::as_array_mut) {
            array.pop();
            truncated.insert(path.to_owned());
        }
        if let Some(object) = value.as_object_mut() {
            let mut paths: Vec<_> = truncated.iter().cloned().collect();
            paths.sort();
            object.insert(
                "truncated_buckets".to_owned(),
                Value::Array(paths.into_iter().map(Value::String).collect()),
            );
        }
    }
    update_attempt_page(&mut value);
    if json_char_count(&value) <= max_chars {
        return value;
    }
    let full_counts = value.get("counts").unwrap_or(&Value::Null);
    let fallback_counts = json!({
        "claims":full_counts.get("claims").cloned().unwrap_or(Value::Null),
        "related_claims_total":full_counts.get("related_claims_total").cloned().unwrap_or(Value::Null),
        "attempts":full_counts.get("attempts").cloned().unwrap_or(Value::Null),
        "recorded_proof_gaps":full_counts.get("recorded_proof_gaps").cloned().unwrap_or(Value::Null),
        "claims_without_enriched_gap_audit":full_counts.get("claims_without_enriched_gap_audit").cloned().unwrap_or(Value::Null)
    });
    // A bounded claim brief must retain its own recorded conclusion before
    // historical attempts or related-record lists.  Without this, a large
    // graph packet could collapse to counts alone even when the root claim has
    // a useful, explicitly scoped overview.
    let root_statement = value
        .pointer("/what_we_have/0")
        .filter(|_| value.get("kind").and_then(Value::as_str) == Some("claim"))
        .map(|item| compact_context_value(item, 1_800, 1));
    let root_is_attempt_record = matches!(
        value.get("kind").and_then(Value::as_str),
        Some("attempt" | "contract")
    );
    let attempts_page = if root_is_attempt_record {
        json!({"offset":value.pointer("/attempts/offset").cloned().unwrap_or(json!(0)),"total":value.pointer("/attempts/total").cloned().unwrap_or(Value::Null),"next_offset":value.pointer("/attempts/offset").cloned().unwrap_or(json!(0)),"notice":"The root record is shown above; page direct linked attempts separately."})
    } else {
        Value::Null
    };
    let mut fallback = json!({
        "id":value.get("id").cloned().or_else(||value.get("ids").cloned()).unwrap_or(Value::Null),
        "kind":value.get("kind").cloned().unwrap_or(Value::Null),
        "question":clip(value.get("question").and_then(Value::as_str).unwrap_or("Research summary"),100),
        "claim_status":value.get("claim_status").cloned().unwrap_or(Value::Null),
        "enrichment_state":value.pointer("/what_we_have/enrichment_state").cloned().unwrap_or_else(|| value.get("enrichment_state").cloned().unwrap_or(Value::Null)),
        "counts":fallback_counts,
        "proof_gap_coverage":value.get("proof_gap_coverage").map(|coverage|json!({"recorded_gaps_from_enriched_claims":coverage.get("recorded_gaps_from_enriched_claims").cloned().unwrap_or(Value::Null),"claims_without_enriched_gap_audit":coverage.get("claims_without_enriched_gap_audit").cloned().unwrap_or(Value::Null),"complete":coverage.get("complete").cloned().unwrap_or(Value::Null),"scope":coverage.get("scope").cloned().unwrap_or(Value::Null),"claims_without_current_annotation_sources":coverage.get("claims_without_current_annotation_sources").cloned().unwrap_or(Value::Null),"notice":coverage.get("notice").cloned().unwrap_or(Value::Null)})).unwrap_or(Value::Null),
        "topic_orientation":value.get("topic_orientation").map(|orientation|json!({"role_de":clip(&value_text(orientation.get("role_de")),60),"have_de":clip(&value_text(orientation.get("have_de")),70),"seek_de":clip(&value_text(orientation.get("seek_de")),70),"next_bridge_de":clip(&value_text(orientation.get("next_bridge_de")),70)})).unwrap_or(Value::Null),
        "scope":value.get("scope").cloned().unwrap_or(Value::Null),
        "what_we_have":root_statement.map(|item|json!([item])).unwrap_or_else(||json!([])),
        "attempt_record":value.get("attempt_record").map(|item|if root_is_attempt_record { compact_root_attempt_record(item,180,2) } else { item.clone() }).unwrap_or(Value::Null),
        "typed_related_records":value.get("typed_related_records").and_then(Value::as_array).map(|items|Value::Array(items.iter().take(if root_is_attempt_record {0}else{4}).map(|item|compact_context_value(item,180,1)).collect())).unwrap_or_else(||json!([])),
        "attempts_page":attempts_page,
        "notice":"Summary compacted to fit requested budget; no inference was added.",
        "truncated":true
    });
    if matches!(
        value.get("kind").and_then(Value::as_str),
        Some("attempt" | "contract")
    ) {
        fallback
            .as_object_mut()
            .map(|object| object.remove("attempts_page"));
    } else if fallback.get("attempts_page").is_some_and(Value::is_null) {
        fallback
            .as_object_mut()
            .map(|object| object.remove("attempts_page"));
    }
    if json_char_count(&fallback) <= max_chars {
        fallback
    } else {
        if root_is_attempt_record {
            json!({
                "id":value.get("id"),
                "kind":value.get("kind"),
                "question":clip(value.get("question").and_then(Value::as_str).unwrap_or("Research record"),100),
                "scope":value.get("scope"),
                "attempt_record":value.get("attempt_record").map(|item|compact_root_attempt_record(item,120,1)).unwrap_or(Value::Null),
                "notice":"Root contract retained; linked records and secondary context were omitted for the requested budget.",
                "truncated":true
            })
        } else {
            json!({
                "id":value.get("id"),
                "kind":value.get("kind"),
                "enrichment_state":value.get("enrichment_state"),
                "counts":{"recorded_proof_gaps":value.pointer("/counts/recorded_proof_gaps")},
                "proof_gap_coverage":{
                    "complete":value.pointer("/proof_gap_coverage/complete"),
                    "claims_without_enriched_gap_audit":value.pointer("/proof_gap_coverage/claims_without_enriched_gap_audit"),
                    "claims_without_current_annotation_sources":value.pointer("/proof_gap_coverage/claims_without_current_annotation_sources")
                },
                "what_we_have":value.get("what_we_have").and_then(Value::as_array).and_then(|items|items.first()).map(|item|json!([{"claim_id":item.get("claim_id"),"statement":clip(&value_text(item.get("statement")),80)}])).unwrap_or_else(||json!([])),
                "notice":"Compacted. Historical annotation gaps are not current research verdicts; consult research_state.",
                "truncated":true
            })
        }
    }
}

fn compact_attempt_research_brief(value: &mut Value) {
    let related_fields = [
        "positive_results",
        "resolved_results",
        "conditional_results",
        "leaves_open",
        "scoped_exclusions",
        "not_claimed",
        "unclassified_findings",
    ];
    if let Some(records) = value
        .get_mut("typed_related_records")
        .and_then(Value::as_array_mut)
    {
        let total = records.len();
        for record in records.iter_mut() {
            for (field, max) in [
                ("label", 100),
                ("scope", 150),
                ("original_source_scope_correction", 240),
            ] {
                if let Some(text) = record.get(field).and_then(Value::as_str).map(str::to_owned) {
                    record[field] = json!(clip(&text, max));
                }
            }
            let mut omitted = serde_json::Map::new();
            for field in related_fields {
                if let Some(items) = record.get_mut(field).and_then(Value::as_array_mut) {
                    let count = items.len();
                    items.truncate(1);
                    if let Some(text) = items.first().and_then(Value::as_str).map(str::to_owned) {
                        items[0] = json!(clip(&text, 135));
                    }
                    if count > items.len() {
                        omitted.insert(field.to_owned(), json!(count - items.len()));
                    }
                }
            }
            if !omitted.is_empty() {
                record["omitted_counts"] = Value::Object(omitted);
            }
            if let Some(sources) = record.get_mut("source_refs").and_then(Value::as_array_mut) {
                let count = sources.len();
                sources.truncate(1);
                if count > sources.len() {
                    record["source_refs_omitted"] = json!(count - sources.len());
                }
            }
        }
        records.truncate(6);
        if total > records.len() {
            value["typed_related_records_omitted"] = json!(total - records.len());
        }
    }
    if let Some(items) = value
        .pointer_mut("/linked_claims/items")
        .and_then(Value::as_array_mut)
    {
        items.truncate(4);
        for item in items {
            for field in ["label", "summary"] {
                if let Some(text) = item.get(field).and_then(Value::as_str).map(str::to_owned) {
                    item[field] = json!(clip(&text, 120));
                }
            }
            if let Some(sources) = item.get_mut("source_refs").and_then(Value::as_array_mut) {
                sources.truncate(1);
            }
        }
    }
    if let Some(connections) = value
        .get_mut("recorded_connections")
        .and_then(Value::as_array_mut)
    {
        let total = connections.len();
        connections.truncate(6);
        if total > connections.len() {
            value["recorded_connections_omitted"] = json!(total - connections.len());
        }
    }
    if let Some(sources) = value.get_mut("source_refs").and_then(Value::as_array_mut) {
        let total = sources.len();
        sources.truncate(5);
        if total > sources.len() {
            value["source_refs_omitted"] = json!(total - sources.len());
        }
    }
    if let Some(sources) = value.get_mut("source_ids").and_then(Value::as_array_mut) {
        sources.truncate(8);
    }
}

fn compact_context_value(value: &Value, text_limit: usize, item_limit: usize) -> Value {
    match value {
        Value::String(text) => json!(clip(text, text_limit)),
        Value::Array(items) => Value::Array(
            items
                .iter()
                .take(item_limit)
                .map(|item| compact_context_value(item, text_limit, item_limit))
                .collect(),
        ),
        Value::Object(fields) => {
            let mut compact = serde_json::Map::new();
            for (key, item) in fields {
                let limit = match key.as_str() {
                    "id" | "source_id" | "path" | "type" | "kind" | "verdict" | "status"
                    | "category" | "origin" | "topic_link_basis" => 300,
                    _ => text_limit,
                };
                let child_item_limit = if key == "not_claimed" { 8 } else { item_limit };
                compact.insert(
                    key.clone(),
                    compact_context_value(item, limit, child_item_limit),
                );
            }
            Value::Object(compact)
        }
        _ => value.clone(),
    }
}

/// Keep the requested contract itself readable when a research packet must
/// fall back to a compact form. Every shortened outcome list carries forward
/// its previous omission count and adds the number removed here.
fn compact_root_attempt_record(value: &Value, text_limit: usize, source_limit: usize) -> Value {
    if !value.is_object() {
        return Value::Null;
    }
    let mut compact = serde_json::Map::new();
    for field in ["id", "label", "verdict", "question_or_scope", "notice"] {
        if let Some(item) = value.get(field) {
            let bounded = if let Some(text) = item.as_str() {
                json!(clip(text, if field == "id" { 300 } else { text_limit }))
            } else {
                item.clone()
            };
            compact.insert(field.to_owned(), bounded);
        }
    }
    if let Some(correction) = value
        .get("original_source_scope_correction")
        .and_then(Value::as_str)
    {
        compact.insert(
            "original_source_scope_correction".to_owned(),
            json!(clip(correction, text_limit.min(300))),
        );
    }

    let source_refs = value
        .get("source_refs")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let source_refs_omitted = value
        .get("source_refs_omitted")
        .and_then(Value::as_u64)
        .unwrap_or(0)
        .saturating_add(source_refs.len().saturating_sub(source_limit) as u64);
    compact.insert(
        "source_refs".to_owned(),
        Value::Array(
            source_refs
                .iter()
                .take(source_limit)
                .map(|source| compact_context_value(source, 120, 1))
                .collect(),
        ),
    );
    if source_refs_omitted > 0 {
        compact.insert("source_refs_omitted".to_owned(), json!(source_refs_omitted));
    }

    if let Some(conclusions) = value.get("recorded_conclusions").and_then(Value::as_object) {
        let mut compact_conclusions = serde_json::Map::new();
        let mut omitted_counts = conclusions
            .get("omitted_counts")
            .and_then(Value::as_object)
            .cloned()
            .unwrap_or_default();
        for (field, limit) in [
            ("positive_results", 6usize),
            ("resolved_results", 6),
            ("conditional_results", 4),
            ("assumptions", 6),
            ("rules_out", 2),
            ("missing_or_open", 3),
            ("not_claimed", 3),
            ("supporting_findings", 2),
            ("findings", 2),
        ] {
            let items = conclusions
                .get(field)
                .and_then(Value::as_array)
                .cloned()
                .unwrap_or_default();
            let retained: Vec<Value> = items
                .iter()
                .take(limit)
                .map(|item| match item.as_str() {
                    Some(text) => json!(clip(text, text_limit)),
                    None => compact_context_value(item, text_limit, 1),
                })
                .collect();
            if !retained.is_empty() {
                compact_conclusions.insert(field.to_owned(), Value::Array(retained.clone()));
            }
            let prior_omitted = omitted_counts
                .get(field)
                .and_then(Value::as_u64)
                .unwrap_or(0);
            let omitted =
                prior_omitted.saturating_add(items.len().saturating_sub(retained.len()) as u64);
            if omitted > 0 {
                omitted_counts.insert(field.to_owned(), json!(omitted));
            }
        }
        if let Some(recorded) = conclusions.get("missing_or_open_recorded") {
            compact_conclusions.insert("missing_or_open_recorded".to_owned(), recorded.clone());
        }
        if !omitted_counts.is_empty() {
            compact_conclusions.insert("omitted_counts".to_owned(), Value::Object(omitted_counts));
        }
        compact.insert(
            "recorded_conclusions".to_owned(),
            Value::Object(compact_conclusions),
        );
    }
    Value::Object(compact)
}

fn fit_research_context_budget(mut value: Value, max_chars: usize) -> Value {
    if json_char_count(&value) <= max_chars {
        return value;
    }
    let original = value.clone();
    let buckets = [
        "/recorded_connections",
        "/source_ids",
        "/related_claims/items",
        "/object_records",
        "/assumption_records",
        "/prior_attempts/items",
        "/scoped_negatives",
        "/recorded_shared_inputs/objects",
        "/recorded_shared_inputs/structures",
        "/recorded_shared_inputs/seed_links",
        "/recorded_shared_inputs/seeds",
        "/what_we_have",
        "/next_discriminators",
        "/typed_related_records",
        "/candidate_connections",
    ];
    let mut omitted = HashMap::<String, usize>::new();
    let what_we_have_minimum = if max_chars >= 16_000 {
        6
    } else if max_chars >= 10_000 {
        4
    } else {
        1
    };
    while json_char_count(&value) > max_chars {
        let candidate = buckets.iter().find_map(|path| {
            let minimum = match *path {
                "/prior_attempts/items"
                | "/typed_related_records"
                | "/related_claims/items"
                | "/candidate_connections"
                | "/scoped_negatives"
                | "/next_discriminators" => 1,
                "/what_we_have" => what_we_have_minimum,
                _ => 0,
            };
            value
                .pointer(path)
                .and_then(Value::as_array)
                .filter(|items| items.len() > minimum)
                .map(|items| (*path, items.len()))
        });
        let Some((path, _)) = candidate else { break };
        if let Some(items) = value.pointer_mut(path).and_then(Value::as_array_mut) {
            items.pop();
            *omitted.entry(path.to_owned()).or_default() += 1;
        }
        if let Some(n) = omitted.get("/prior_attempts/items") {
            let prior = value
                .get_mut("prior_attempts")
                .expect("prior attempts packet");
            prior["next_offset"] = json!(
                prior["offset"].as_u64().unwrap_or(0) as usize
                    + prior["items"].as_array().map(Vec::len).unwrap_or(0)
            );
            prior["truncated"] = json!(true);
            prior["items_omitted_for_budget"] = json!(n);
        }
        if let Some(coverage) = value.get_mut("coverage") {
            let mut truncated: Vec<Value> = omitted
                .iter()
                .map(|(path, count)| json!({"bucket":path,"omitted":count}))
                .collect();
            truncated.sort_by(|a, b| a["bucket"].as_str().cmp(&b["bucket"].as_str()));
            coverage["omissions_for_budget"] = Value::Array(truncated);
        }
        value["truncated"] = json!(true);
    }
    if json_char_count(&value) <= max_chars {
        return value;
    }
    let full_coverage = original.get("coverage").unwrap_or(&Value::Null);
    let prior = original.get("prior_attempts").unwrap_or(&Value::Null);
    let prior_items = prior
        .get("items")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let root_is_attempt_record = matches!(
        original.get("kind").and_then(Value::as_str),
        Some("attempt" | "contract")
    );
    let compact_prior: Vec<Value> = if root_is_attempt_record {
        Vec::new()
    } else {
        prior_items
            .iter()
            .take(1)
            .map(|item| compact_context_value(item, 120, 1))
            .collect()
    };
    let coverage = json!({
        "claims_total":full_coverage.get("claims_total"),
        "claims_in_page":full_coverage.get("claims_in_page"),
        "claims_next_offset":full_coverage.get("claims_next_offset"),
        "claims_truncated":full_coverage.get("claims_truncated"),
        "recorded_attempt_total":full_coverage.get("recorded_attempt_total"),
        "related_claims_total":full_coverage.get("related_claims_total"),
        "related_claims_in_packet":0,
        "related_claims_omitted":full_coverage.get("related_claims_total"),
        "premise_records_total":full_coverage.get("premise_records_total"),
        "premise_records_in_packet":0,
        "premise_records_omitted":full_coverage.get("premise_records_total"),
        "assumption_records_total":full_coverage.get("assumption_records_total"),
        "assumption_records_in_packet":0,
        "assumption_records_omitted":full_coverage.get("assumption_records_total"),
        "unresolved_semantic_links":full_coverage.get("unresolved_semantic_links"),
        "graph_integrity_issue_total":full_coverage.get("graph_integrity_issue_total"),
        "omissions_for_budget":[]
    });
    let optional_arrays: Vec<(&str, usize)> = if root_is_attempt_record {
        // The root contract card already carries its own results and source
        // paths. Do not spend its fallback budget repeating topic negatives,
        // links, or the same conclusions under what_we_have.
        vec![("source_refs", 2)]
    } else {
        vec![
            ("premise_records", 3),
            ("assumption_records", 3),
            ("what_we_have", what_we_have_minimum),
            ("next_discriminators", 2),
            ("source_refs", 3),
            ("scoped_negatives", 1),
            ("recorded_connections", 2),
            ("candidate_connections", 1),
        ]
    };
    let mut fallback = json!({
        "id":original.get("id"),
        "kind":original.get("kind"),
        "question":compact_context_value(original.get("question").unwrap_or(&Value::Null),120,1),
        "claim_status":original.get("claim_status"),
        "goal_status":original.get("goal_status"),
        "scope":compact_context_value(original.get("scope").unwrap_or(&Value::Null),180,1),
        "scope_recorded":original.get("scope_recorded"),
        "coverage":coverage,
        "attempt_record":original.get("attempt_record").map(|item|if root_is_attempt_record { compact_root_attempt_record(item,180,2) } else { compact_context_value(item,320,3) }).unwrap_or(Value::Null),
        "prior_attempts":{"items":compact_prior,"total":prior.get("total"),"offset":prior.get("offset"),"next_offset":prior.get("offset").and_then(Value::as_u64).unwrap_or(0)+compact_prior.len() as u64,"truncated":prior.get("total").and_then(Value::as_u64).unwrap_or(0)>compact_prior.len() as u64,"items_omitted_for_budget":prior_items.len().saturating_sub(compact_prior.len()),"notice":if root_is_attempt_record{"The root contract is shown under attempt_record; this page is for linked historical records."}else{"Historical records remain scoped and available through list_attempts."}},
        "notice":"Research context was compacted to preserve recorded status, premises, results, and source paths. Omitted records remain available through paged tools; no inference was added.",
        "max_chars":max_chars,
        "truncated":true
    });
    for (path, limit) in optional_arrays {
        let items: Vec<Value> = original
            .get(path)
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .take(limit)
            .map(|item| {
                let text_limit = if path == "what_we_have" { 1_200 } else { 120 };
                compact_context_value(item, text_limit, 2)
            })
            .collect();
        let total = original
            .get(path)
            .and_then(Value::as_array)
            .map_or(0, Vec::len);
        if !items.is_empty() {
            fallback[path] = Value::Array(items.clone());
        }
        if path == "premise_records" {
            fallback["coverage"]["premise_records_in_packet"] = json!(items.len());
            fallback["coverage"]["premise_records_omitted"] = json!(original
                .pointer("/coverage/premise_records_total")
                .and_then(Value::as_u64)
                .unwrap_or(total as u64)
                .saturating_sub(items.len() as u64));
        } else if path == "assumption_records" {
            fallback["coverage"]["assumption_records_in_packet"] = json!(items.len());
            fallback["coverage"]["assumption_records_omitted"] = json!(original
                .pointer("/coverage/assumption_records_total")
                .and_then(Value::as_u64)
                .unwrap_or(total as u64)
                .saturating_sub(items.len() as u64));
        }
        if total > items.len() {
            let mut omissions = fallback["coverage"]["omissions_for_budget"]
                .as_array()
                .cloned()
                .unwrap_or_default();
            omissions.push(json!({"bucket":path,"omitted":total-items.len()}));
            fallback["coverage"]["omissions_for_budget"] = Value::Array(omissions);
        }
    }
    fallback["source_ids_omitted"] = json!(
        original
            .get("source_ids")
            .and_then(Value::as_array)
            .map_or(0, |items| items.len().saturating_sub(3))
            + original
                .get("source_ids_omitted")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize
    );
    if json_char_count(&fallback) <= max_chars {
        fallback
    } else {
        let minimal = json!({
            "id":original.get("id"),
            "kind":original.get("kind"),
            "question":compact_context_value(original.get("question").unwrap_or(&Value::Null),80,1),
            "claim_status":original.get("claim_status"),
            "scope":compact_context_value(original.get("scope").unwrap_or(&Value::Null),100,1),
            "scope_recorded":original.get("scope_recorded"),
            "coverage":{"claims_total":full_coverage.get("claims_total"),"recorded_attempt_total":full_coverage.get("recorded_attempt_total"),"premise_records_total":full_coverage.get("premise_records_total"),"graph_integrity_issue_total":full_coverage.get("graph_integrity_issue_total")},
            "premise_records":original.get("premise_records").and_then(Value::as_array).map(|items|Value::Array(items.iter().take(1).map(|item|compact_context_value(item,if max_chars >= 10_000 { 400 } else { 120 },1)).collect())).unwrap_or_else(||json!([])),
            "attempt_record":original.get("attempt_record").map(|item|if root_is_attempt_record { compact_root_attempt_record(item,if max_chars >= 10_000 { 220 } else { 150 },2) } else { compact_context_value(item,if max_chars >= 10_000 { 400 } else { 120 },3) }).unwrap_or(Value::Null),
            "what_we_have":if root_is_attempt_record { json!([]) } else { original.get("what_we_have").and_then(Value::as_array).map(|items|Value::Array(items.iter().take(what_we_have_minimum).map(|item|compact_context_value(item,if max_chars >= 10_000 { 400 } else { 120 },1)).collect())).unwrap_or_else(||json!([])) },
            "notice":"Research context was reduced to the requested bound; status and any retained premise remain explicitly scoped.",
            "max_chars":max_chars,
            "truncated":true
        });
        if json_char_count(&minimal) <= max_chars {
            minimal
        } else {
            json!({"id":original.get("id"),"kind":original.get("kind"),"notice":"Research context was reduced to the requested character bound.","max_chars":max_chars,"truncated":true})
        }
    }
}

fn update_attempt_page(value: &mut Value) {
    let Some(attempts) = value.get_mut("attempts").and_then(Value::as_object_mut) else {
        return;
    };
    let offset = attempts.get("offset").and_then(Value::as_u64).unwrap_or(0);
    let total = attempts.get("total").and_then(Value::as_u64).unwrap_or(0);
    let consumed = attempts
        .get("items")
        .and_then(Value::as_array)
        .map_or(0, |items| items.len() as u64);
    let next_offset = offset.saturating_add(consumed);
    attempts.insert("next_offset".into(), json!(next_offset));
    attempts.insert("truncated".into(), json!(next_offset < total));
}

fn compact_topic_orientation(value: Option<&Value>) -> Value {
    let Some(value) = value.and_then(Value::as_object) else {
        return Value::Null;
    };
    let mut compact = serde_json::Map::new();
    for key in [
        "role_de",
        "role_en",
        "have_de",
        "have_en",
        "seek_de",
        "seek_en",
        "connection_de",
        "connection_en",
        "next_bridge_de",
        "next_bridge_en",
    ] {
        if let Some(text) = value.get(key).and_then(Value::as_str) {
            compact.insert(key.to_owned(), json!(clip(text, 130)));
        }
    }
    Value::Object(compact)
}

fn clip(value: &str, max_chars: usize) -> String {
    let mut chars = value.chars();
    let clipped: String = chars.by_ref().take(max_chars).collect();
    if chars.next().is_some() {
        format!("{clipped}…")
    } else {
        clipped
    }
}

fn value_text(value: Option<&Value>) -> String {
    match value {
        Some(Value::String(value)) => value.clone(),
        Some(Value::Null) | None => String::new(),
        Some(value) => value.to_string(),
    }
}

fn assessment_category(node: &Node) -> Option<&str> {
    node.attrs
        .pointer("/assessment/category")
        .and_then(Value::as_str)
}

fn assessment_scope(node: &Node) -> Option<&str> {
    node.attrs
        .pointer("/assessment/scope")
        .and_then(Value::as_str)
}

fn collect_source_refs(value: &Value, out: &mut HashSet<String>) {
    match value {
        Value::Object(object) => {
            if let Some(values) = object.get("source_refs").and_then(Value::as_array) {
                for reference in values {
                    if let Some(id) = reference
                        .as_str()
                        .or_else(|| reference.get("source_id").and_then(Value::as_str))
                    {
                        out.insert(id.to_owned());
                    }
                }
            }
            for value in object.values() {
                collect_source_refs(value, out);
            }
        }
        Value::Array(values) => {
            for value in values {
                collect_source_refs(value, out);
            }
        }
        _ => {}
    }
}

fn source_ref_ids(value: &Value) -> Vec<String> {
    let mut refs = HashSet::new();
    collect_source_refs(value, &mut refs);
    let mut refs: Vec<String> = refs.into_iter().collect();
    refs.sort();
    refs
}

fn input_importance(input: &Value) -> (u8, String) {
    let mut impacts: Vec<&str> = input
        .get("impact")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(Value::as_str)
        .collect();
    if let Some(impact) = input.get("impact").and_then(Value::as_str) {
        impacts.push(impact);
    }
    let origin = input
        .get("origin")
        .and_then(Value::as_str)
        .unwrap_or("unknown");
    let rank = if impacts.contains(&"circularity_risk") {
        0
    } else if impacts.contains(&"missing_derivation") {
        1
    } else if impacts.contains(&"conditional") {
        2
    } else {
        match origin {
            "unknown" => 3,
            "declared_axiom" => 4,
            "hardcoded" | "free_parameter" | "external_data" => 5,
            "upstream_result" | "derived_in_script" => 6,
            _ => 7,
        }
    };
    (rank, value_text(input.get("name")))
}

fn input_is_research_relevant(input: &Value) -> bool {
    let origin = input
        .get("origin")
        .and_then(Value::as_str)
        .unwrap_or("unknown");
    let impact = input.get("impact").and_then(Value::as_str).unwrap_or("");
    matches!(
        origin,
        "hardcoded" | "free_parameter" | "unknown" | "declared_axiom" | "external_data"
    ) || matches!(
        impact,
        "conditional" | "circularity_risk" | "missing_derivation"
    ) || input
        .get("impact")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(Value::as_str)
        .any(|impact| {
            matches!(
                impact,
                "conditional" | "circularity_risk" | "missing_derivation"
            )
        })
}

fn node_source_ids(node: &Node) -> Vec<String> {
    source_ref_ids(&node.attrs)
}

fn collect_node_source_ids(node: &Node, output: &mut HashSet<String>) {
    output.extend(node_source_ids(node));
}

fn topic_members(
    graph: &KnowledgeGraph,
    topic_id: &str,
    by_id: &HashMap<&str, &Node>,
    include_related: bool,
) -> (
    HashSet<String>,
    HashSet<String>,
    HashSet<String>,
    HashSet<String>,
) {
    let mut claims = HashSet::new();
    let mut attempts = HashSet::new();
    let mut structures = HashSet::new();
    let mut goals = HashSet::new();
    for edge in &graph.edges {
        let other = if edge.src == topic_id {
            Some(edge.dst.as_str())
        } else if edge.dst == topic_id {
            Some(edge.src.as_str())
        } else {
            None
        };
        let Some(other) = other else { continue };
        if let Some(node) = by_id.get(other) {
            match node.kind.as_str() {
                "claim"
                    if matches!(edge.kind.as_str(), "about" | "theme_link" | "theme_links")
                        && (include_related
                            || claim_topic_role(node, topic_id) != Some("related")) =>
                {
                    claims.insert(node.id.clone());
                }
                "attempt"
                    if matches!(edge.kind.as_str(), "about" | "theme_link" | "theme_links")
                        && (include_related
                            || !is_inherited_attempt_topic_link(edge, Some(node))) =>
                {
                    attempts.insert(node.id.clone());
                }
                "structure" => {
                    structures.insert(node.id.clone());
                }
                "goal" => {
                    goals.insert(node.id.clone());
                }
                _ => {}
            }
        }
    }
    for node in graph.nodes.iter().filter(|node| node.kind == "claim") {
        if (include_related || claim_topic_role(node, topic_id) != Some("related"))
            && node
                .attrs
                .pointer("/semantic/topic_ids")
                .and_then(Value::as_array)
                .is_some_and(|topics| topics.iter().any(|value| value.as_str() == Some(topic_id)))
        {
            claims.insert(node.id.clone());
        }
    }
    for edge in &graph.edges {
        if edge.kind != "tests" {
            continue;
        }
        let (attempt, claim) = match (by_id.get(edge.src.as_str()), by_id.get(edge.dst.as_str())) {
            (Some(src), Some(dst)) if src.kind == "attempt" && dst.kind == "claim" => (src, dst),
            (Some(src), Some(dst)) if dst.kind == "attempt" && src.kind == "claim" => (dst, src),
            _ => continue,
        };
        if claims.contains(&claim.id) {
            attempts.insert(attempt.id.clone());
        }
    }
    (claims, attempts, structures, goals)
}

fn is_inherited_attempt_topic_link(edge: &Edge, node: Option<&Node>) -> bool {
    edge.kind == "about"
        && node.is_some_and(|node| node.kind == "attempt")
        && edge.attrs.get("method").and_then(Value::as_str) == Some("linked_claim_topics")
}

fn topic_linked_attempts(
    graph: &KnowledgeGraph,
    topic_id: &str,
    by_id: &HashMap<&str, &Node>,
) -> HashMap<String, String> {
    let mut links = HashMap::new();
    for edge in &graph.edges {
        if !matches!(edge.kind.as_str(), "about" | "theme_link" | "theme_links") {
            continue;
        }
        let other = if edge.src == topic_id {
            Some(edge.dst.as_str())
        } else if edge.dst == topic_id {
            Some(edge.src.as_str())
        } else {
            None
        };
        let Some(attempt) = other
            .and_then(|id| by_id.get(id))
            .filter(|node| node.kind == "attempt")
        else {
            continue;
        };
        let basis = if is_inherited_attempt_topic_link(edge, Some(attempt)) {
            "inherited_from_claim_topic_membership"
        } else {
            "explicit_topic_edge"
        };
        // An explicit topic edge takes precedence if both forms exist.
        links
            .entry(attempt.id.clone())
            .and_modify(|old| {
                if basis == "explicit_topic_edge" {
                    *old = basis.to_owned();
                }
            })
            .or_insert_with(|| basis.to_owned());
    }
    links
}

fn claim_topic_role<'a>(node: &'a Node, topic_id: &str) -> Option<&'a str> {
    node.attrs
        .get("topic_roles")
        .and_then(|roles| roles.get(topic_id))
        .and_then(Value::as_str)
}

fn linked_structures_and_goals(
    graph: &KnowledgeGraph,
    root_id: &str,
    by_id: &HashMap<&str, &Node>,
) -> (HashSet<String>, HashSet<String>) {
    let mut structures = HashSet::new();
    let mut goals = HashSet::new();
    for edge in &graph.edges {
        let other = if edge.src == root_id {
            Some(edge.dst.as_str())
        } else if edge.dst == root_id {
            Some(edge.src.as_str())
        } else {
            None
        };
        let Some(other) = other else { continue };
        match by_id.get(other).map(|node| node.kind.as_str()) {
            Some("structure") => {
                structures.insert(other.to_owned());
            }
            Some("goal") => {
                goals.insert(other.to_owned());
            }
            _ => {}
        }
    }
    (structures, goals)
}

fn linked_structures_and_goals_typed(
    graph: &KnowledgeGraph,
    root_id: &str,
    by_id: &HashMap<&str, &Node>,
) -> (HashSet<String>, HashSet<String>) {
    let allowed = [
        "uses_structure",
        "about_structure",
        "implements",
        "checks",
        "addresses_open_goal",
        "leaves_open",
    ];
    let mut structures = HashSet::new();
    let mut goals = HashSet::new();
    for edge in graph
        .edges
        .iter()
        .filter(|edge| allowed.contains(&edge.kind.as_str()))
    {
        let other = if edge.src == root_id {
            Some(edge.dst.as_str())
        } else if edge.dst == root_id {
            Some(edge.src.as_str())
        } else {
            None
        };
        let Some(other) = other else { continue };
        match by_id.get(other).map(|node| node.kind.as_str()) {
            Some("structure") => {
                structures.insert(other.to_owned());
            }
            Some("goal") => {
                goals.insert(other.to_owned());
            }
            _ => {}
        }
    }
    (structures, goals)
}

fn normalize_claim_ref(value: &str, by_id: &HashMap<&str, &Node>) -> Option<String> {
    if by_id.get(value).is_some_and(|node| node.kind == "claim") {
        return Some(value.to_owned());
    }
    let candidate = format!("claim:{value}");
    by_id
        .get(candidate.as_str())
        .is_some_and(|node| node.kind == "claim")
        .then_some(candidate)
}

fn source_path_refs(
    graph: &KnowledgeGraph,
    ids: impl IntoIterator<Item = String>,
    limit: usize,
) -> Vec<Value> {
    let wanted: HashSet<String> = ids.into_iter().collect();
    let mut matching: Vec<&Source> = graph
        .sources
        .iter()
        .filter(|source| wanted.contains(&source.id))
        .collect();
    matching.sort_by(|a, b| {
        let priority = |source: &Source| {
            let path = source.path.to_lowercase();
            if path.ends_with("/proof.md") || path == "proof.md" || path.ends_with("/proof.txt") {
                0
            } else if path.ends_with("/contract_index.json") {
                1
            } else if path.ends_with("/checker.py") {
                2
            } else if path.ends_with("/readme.md") {
                3
            } else if path.ends_with("/results.json") {
                4
            } else {
                5
            }
        };
        priority(a)
            .cmp(&priority(b))
            .then_with(|| a.path.cmp(&b.path))
            .then_with(|| a.id.cmp(&b.id))
    });
    let refs: Vec<Value> = matching
        .into_iter()
        .take(limit)
        .map(|source| json!({"id":source.id,"path":source.path,"line_start":source.line_start,"line_end":source.line_end,"page":source.page}))
        .collect();
    refs
}

fn classified_attempt_findings(node: &Node) -> Value {
    let contract = node.attrs.get("original_contract").unwrap_or(&node.attrs);
    let mut buckets = HashMap::<&str, Vec<String>>::from([
        ("positive_results", Vec::new()),
        ("resolved_results", Vec::new()),
        ("conditional_results", Vec::new()),
        ("leaves_open", Vec::new()),
        ("scoped_exclusions", Vec::new()),
        ("not_claimed", Vec::new()),
        ("supporting_findings", Vec::new()),
        ("unclassified_findings", Vec::new()),
    ]);
    let mut totals = HashMap::<String, HashSet<String>>::new();
    let mut push = |bucket: &str, text: String| {
        let text = clip(text.trim(), 320);
        if !text.is_empty() {
            let unique = totals.entry(bucket.to_owned()).or_default();
            if !unique.insert(text.clone()) {
                return;
            }
            if let Some(values) = buckets.get_mut(bucket) {
                if values.len() < 8 {
                    values.push(text);
                }
            }
        }
    };
    for (key, default_bucket) in [
        ("positive_results", "positive_results"),
        ("resolved", "resolved_results"),
        ("findings", "unclassified_findings"),
        ("analytic_conditional_results", "conditional_results"),
        ("conditional_results", "conditional_results"),
        ("leaves_open", "leaves_open"),
        ("missing", "leaves_open"),
        ("not_derived", "leaves_open"),
        ("not_claimed", "not_claimed"),
        ("kills", "scoped_exclusions"),
        ("kills_free", "scoped_exclusions"),
        ("rules_out", "scoped_exclusions"),
        ("excludes", "scoped_exclusions"),
    ] {
        for container in [Some(contract), Some(&node.attrs)] {
            let Some(items) = container
                .and_then(|container| container.get(key))
                .and_then(Value::as_array)
            else {
                continue;
            };
            for item in items {
                let mut parsed_string = None;
                let record = if let Some(text) = item.as_str() {
                    if text.trim_start().starts_with('{') {
                        parsed_string = serde_json::from_str::<Value>(text).ok();
                    }
                    parsed_string.as_ref().unwrap_or(item)
                } else {
                    item
                };
                let kind = record
                    .get("kind")
                    .and_then(Value::as_str)
                    .map(str::to_ascii_lowercase);
                let is_verdict_summary = key == "findings"
                    && record
                        .get("text")
                        .or_else(|| record.get("statement"))
                        .and_then(Value::as_str)
                        .is_some_and(|text| {
                            text.to_ascii_lowercase().starts_with("overall verdict ")
                        });
                let bucket = if key == "not_claimed" {
                    // Explicit non-claims bound the record's scope. Keep them
                    // separate even when an importer attached a generic kind.
                    "not_claimed"
                } else {
                    match kind.as_deref() {
                        Some("open" | "missing" | "not_derived") => "leaves_open",
                        Some("refutes" | "refute" | "excludes" | "contradicts") => {
                            "scoped_exclusions"
                        }
                        Some("supports" | "support") if is_verdict_summary => "supporting_findings",
                        Some("supports" | "support" | "positive" | "established") => {
                            "positive_results"
                        }
                        Some("conditional") => "conditional_results",
                        Some(_) => "unclassified_findings",
                        None => default_bucket,
                    }
                };
                let text = record
                    .get("text")
                    .or_else(|| record.get("statement"))
                    .or_else(|| record.get("label"))
                    .or_else(|| record.get("result"))
                    .or_else(|| record.get("finding"))
                    .or_else(|| record.get("description"))
                    .map(|value| value_text(Some(value)))
                    .or_else(|| item.as_str().map(str::to_owned))
                    .unwrap_or_else(|| value_text(Some(record)));
                push(bucket, text);
            }
        }
    }
    let mut result = serde_json::Map::new();
    for (key, values) in buckets {
        result.insert(key.to_owned(), json!(values));
    }
    let mut totals_json = serde_json::Map::new();
    for (key, values) in totals {
        totals_json.insert(key, json!(values.len()));
    }
    result.insert("totals".to_owned(), Value::Object(totals_json));
    result.insert(
        "original_source_scope_correction".to_owned(),
        contract
            .get("original_source_scope_correction")
            .or_else(|| node.attrs.get("original_source_scope_correction"))
            .cloned()
            .unwrap_or(Value::Null),
    );
    Value::Object(result)
}

fn direct_typed_related_records(
    root: &Node,
    graph: &KnowledgeGraph,
    by_id: &HashMap<&str, &Node>,
) -> Vec<Value> {
    let mut links = HashSet::<(String, String)>::new();
    let allowed = ["uses", "extends", "refines", "references"];
    let contract = root.attrs.get("original_contract").unwrap_or(&root.attrs);
    for relation in contract
        .get("relations")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        let Some(kind) = relation
            .get("type")
            .and_then(Value::as_str)
            .filter(|kind| allowed.contains(kind))
        else {
            continue;
        };
        let Some(target) = relation.get("target").and_then(Value::as_str) else {
            continue;
        };
        let target_id = if by_id.contains_key(target) {
            target.to_owned()
        } else if by_id.contains_key(format!("contract:{target}").as_str()) {
            format!("contract:{target}")
        } else if by_id.contains_key(format!("attempt:{target}").as_str()) {
            format!("attempt:{target}")
        } else {
            continue;
        };
        if by_id
            .get(target_id.as_str())
            .is_some_and(|node| matches!(node.kind.as_str(), "attempt" | "contract"))
        {
            links.insert((kind.to_owned(), target_id));
        }
    }
    for edge in graph
        .edges
        .iter()
        .filter(|edge| edge.src == root.id || edge.dst == root.id)
    {
        if !allowed.contains(&edge.kind.as_str()) {
            continue;
        }
        let target = if edge.src == root.id {
            &edge.dst
        } else {
            &edge.src
        };
        if by_id
            .get(target.as_str())
            .is_some_and(|node| matches!(node.kind.as_str(), "attempt" | "contract"))
        {
            links.insert((edge.kind.clone(), target.clone()));
        }
    }
    let mut records = Vec::new();
    for (relation, target_id) in links {
        let Some(target) = by_id.get(target_id.as_str()).copied() else {
            continue;
        };
        let refs = node_source_ids(target);
        let findings = classified_attempt_findings(target);
        records.push(json!({
            "relation":relation,
            "id":target.id,
            "label":clip(&target.label,80),
            "verdict":attempt_verdict(target),
            "scope":clip(&attempt_scope(target),80),
            "positive_results":findings["positive_results"],
            "resolved_results":findings["resolved_results"],
            "conditional_results":findings["conditional_results"],
            "leaves_open":findings["leaves_open"],
            "scoped_exclusions":findings["scoped_exclusions"],
            "not_claimed":findings["not_claimed"],
            "supporting_findings":findings["supporting_findings"],
            "unclassified_findings":findings["unclassified_findings"],
            "original_source_scope_correction":findings["original_source_scope_correction"],
            "source_refs":source_path_refs(graph,refs.into_iter(),1)
        }));
    }
    records.sort_by(|a, b| {
        a["relation"]
            .as_str()
            .cmp(&b["relation"].as_str())
            .then_with(|| a["id"].as_str().cmp(&b["id"].as_str()))
    });
    records.truncate(8);
    records
}

fn goal_claim_ids(
    graph: &KnowledgeGraph,
    goal: &Node,
    by_id: &HashMap<&str, &Node>,
) -> HashSet<String> {
    let mut claims = HashSet::new();
    for edge in &graph.edges {
        if edge.kind != "addresses_open_goal" {
            continue;
        }
        let claim_id = if edge.src == goal.id {
            Some(edge.dst.as_str())
        } else if edge.dst == goal.id {
            Some(edge.src.as_str())
        } else {
            None
        };
        if let Some(claim) = claim_id
            .and_then(|claim_id| by_id.get(claim_id))
            .filter(|claim| claim.kind == "claim")
        {
            claims.insert(claim.id.clone());
        }
    }
    if let Some(claim_ids) = goal.attrs.get("claim_ids").and_then(Value::as_array) {
        for claim_id in claim_ids.iter().filter_map(Value::as_str) {
            let full_id = if claim_id.starts_with("claim:") {
                claim_id.to_owned()
            } else {
                format!("claim:{claim_id}")
            };
            if by_id
                .get(full_id.as_str())
                .is_some_and(|claim| claim.kind == "claim")
            {
                claims.insert(full_id);
            }
        }
    }
    claims
}

fn goal_required_steps(goal: &Node) -> Vec<String> {
    let mut steps: Vec<String> = goal
        .attrs
        .get("required_steps")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(Value::as_str)
        .filter(|step| !step.trim().is_empty())
        .map(str::to_owned)
        .collect();
    if steps.is_empty() {
        let fallback = goal_requirement(goal);
        if !fallback.is_empty() {
            steps.push(fallback);
        }
    }
    steps
}

fn topic_related_claim_score(claim: &Node, topic: &Node) -> usize {
    let mut topic_text = format!("{} {}", topic.label, topic.summary);
    for key in ["keywords", "topic_orientation"] {
        if let Some(value) = topic.attrs.get(key) {
            topic_text.push(' ');
            topic_text.push_str(&value_text(Some(value)));
        }
    }
    let mut claim_text = format!("{} {}", claim.label, claim.summary);
    for key in ["overview_de", "title_de", "interpretation_de"] {
        if let Some(value) = claim.attrs.pointer(&format!("/semantic/{key}")) {
            claim_text.push(' ');
            claim_text.push_str(&value_text(Some(value)));
        }
    }
    let topic_tokens: HashSet<String> = tokenize(&topic_text).into_iter().collect();
    let claim_tokens: HashSet<String> = tokenize(&claim_text).into_iter().collect();
    topic_tokens.intersection(&claim_tokens).count()
}

fn attempt_verdict(node: &Node) -> Option<String> {
    [
        node.attrs.get("verdict_raw"),
        node.attrs.get("verdict"),
        node.attrs.pointer("/original_contract/verdict"),
        node.attrs.pointer("/contract/verdict"),
        node.attrs.get("status"),
        node.attrs.get("state"),
        node.attrs.pointer("/assessment/category"),
    ]
    .into_iter()
    .flatten()
    .find_map(Value::as_str)
    .map(str::to_owned)
}

fn attempt_scope(node: &Node) -> String {
    [
        node.attrs.pointer("/assessment/scope"),
        node.attrs.get("scope"),
        node.attrs.pointer("/original_contract/scope"),
        node.attrs.pointer("/contract/scope"),
        node.attrs.pointer("/original_contract/checker_scope"),
    ]
    .into_iter()
    .flatten()
    .find_map(Value::as_str)
    .unwrap_or("")
    .to_owned()
}

fn attempt_question_or_scope(node: &Node) -> String {
    let scope = attempt_scope(node);
    if !scope.is_empty() {
        return scope;
    }
    [
        node.attrs.pointer("/original_contract/question"),
        node.attrs.pointer("/contract/question"),
        node.attrs.get("question"),
    ]
    .into_iter()
    .flatten()
    .find_map(Value::as_str)
    .unwrap_or("")
    .to_owned()
}

fn attempt_card(node: &Node, source_ids: &[String]) -> Value {
    let contract = node.attrs.get("original_contract").unwrap_or(&node.attrs);
    let compact_list = |keys: &[&str]| -> Vec<String> {
        let mut seen = HashSet::new();
        keys.iter()
            // Importers can retain an empty normalized field on the outer
            // attempt while the original contract carries the populated
            // field (or vice versa). Consider both and let a non-empty source
            // value win over an empty placeholder.
            .flat_map(|key| [contract.get(*key), node.attrs.get(*key)])
            .filter_map(|value| value.and_then(Value::as_array))
            .flatten()
            .filter_map(|value| {
                let item = clip(&value_text(Some(value)), 140);
                seen.insert(item.clone()).then_some(item)
            })
            .take(8)
            .collect()
    };
    let omitted_count = |keys: &[&str], emitted: usize| -> usize {
        let mut unique = HashSet::new();
        for key in keys {
            for container in [Some(contract), Some(&node.attrs)] {
                if let Some(items) = container
                    .and_then(|container| container.get(*key))
                    .and_then(Value::as_array)
                {
                    for item in items {
                        unique.insert(clip(&value_text(Some(item)), 140));
                    }
                }
            }
        }
        unique.len().saturating_sub(emitted)
    };
    let assumptions = compact_list(&["assumptions", "baseline_requirements", "required_baselines"]);
    let classified = classified_attempt_findings(node);
    let positive_results = classified["positive_results"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let resolved_results = classified["resolved_results"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let conditional_results = classified["conditional_results"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let rules_out = classified["scoped_exclusions"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let not_claimed = classified["not_claimed"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let missing_or_open = classified["leaves_open"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    // Keep the legacy `findings` alias for positive results while exposing
    // explicit outcome buckets separately. Consumers can deduplicate the
    // alias without losing compatibility with older dossier readers.
    let mut findings = classified["unclassified_findings"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    findings.extend(positive_results.iter().cloned());
    let supporting_findings = classified["supporting_findings"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    findings.extend(supporting_findings.iter().cloned());
    let missing_or_open_recorded = ["missing", "leaves_open", "not_derived"]
        .iter()
        .any(|key| contract.get(*key).is_some() || node.attrs.get(*key).is_some());
    let classified_omitted = |bucket: &str, emitted: usize| -> usize {
        classified
            .pointer(&format!("/totals/{bucket}"))
            .and_then(Value::as_u64)
            .unwrap_or(emitted as u64)
            .saturating_sub(emitted as u64) as usize
    };
    let assumptions_omitted = omitted_count(
        &["assumptions", "baseline_requirements", "required_baselines"],
        assumptions.len(),
    );
    let positive_omitted = classified_omitted("positive_results", positive_results.len());
    let resolved_omitted = classified_omitted("resolved_results", resolved_results.len());
    let conditional_omitted = classified_omitted("conditional_results", conditional_results.len());
    let rules_out_omitted = classified_omitted("scoped_exclusions", rules_out.len());
    let missing_omitted = classified_omitted("leaves_open", missing_or_open.len());
    let not_claimed_omitted = classified_omitted("not_claimed", not_claimed.len());
    json!({
        "id":node.id,
        "label":clip(&node.label,90),
        "verdict":attempt_verdict(node),
        "question_or_scope":clip(&attempt_question_or_scope(node),120),
        "recorded_conclusions":{"findings":findings,"positive_results":positive_results,"supporting_findings":supporting_findings,"resolved_results":resolved_results,"conditional_results":conditional_results,"assumptions":assumptions,"rules_out":rules_out,"missing_or_open":missing_or_open,"not_claimed":not_claimed,"omitted_counts":{"assumptions":assumptions_omitted,"positive_results":positive_omitted,"resolved_results":resolved_omitted,"conditional_results":conditional_omitted,"rules_out":rules_out_omitted,"missing_or_open":missing_omitted,"not_claimed":not_claimed_omitted},"missing_or_open_recorded":missing_or_open_recorded},
        "original_source_scope_correction":contract.get("original_source_scope_correction").or_else(||node.attrs.get("original_source_scope_correction")).map(|value|clip(&value_text(Some(value)),300)),
        "typed_relations":contract.get("relations").and_then(Value::as_array).into_iter().flatten().take(8).filter_map(|relation| {
            Some(json!({"type":relation.get("type")?,"target":relation.get("target")?}))
        }).collect::<Vec<_>>(),
        "source_ids":source_ids.iter().take(3).collect::<Vec<_>>(),
        "notice":"Stored attempt record; its verdict is not a proof of a linked claim."
    })
}

fn attempt_what_we_have(node: &Node, card: &Value) -> Vec<Value> {
    let conclusions = card
        .pointer("/recorded_conclusions")
        .unwrap_or(&Value::Null);
    let scope = attempt_scope(node);
    let source_ids: Vec<String> = node_source_ids(node).into_iter().take(4).collect();
    let mut items = Vec::new();
    for (field, category) in [
        ("positive_results", "recorded_positive_result"),
        ("resolved_results", "recorded_resolution"),
        ("conditional_results", "recorded_conditional_result"),
        ("findings", "recorded_finding"),
        ("supporting_findings", "recorded_supporting_finding"),
    ] {
        for statement in conclusions
            .get(field)
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
        {
            items.push(json!({
                "record_id":node.id,
                "category":category,
                "scope":if scope.is_empty(){Value::Null}else{json!(clip(&scope,240))},
                "statement":clip(statement,400),
                "source_ids":source_ids,
                "notice":"Recorded contract outcome within the stated scope; it does not upgrade the claim status or establish a broader proof."
            }));
            if items.len() >= 6 {
                return items;
            }
        }
    }
    items
}

fn remove_duplicate_findings(card: &mut Value) {
    let Some(conclusions) = card
        .get_mut("recorded_conclusions")
        .and_then(Value::as_object_mut)
    else {
        return;
    };
    let findings = conclusions
        .get("findings")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let positive = conclusions
        .get("positive_results")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    if findings == positive {
        conclusions.remove("findings");
    } else if let Some(items) = conclusions
        .get_mut("findings")
        .and_then(Value::as_array_mut)
    {
        items.retain(|item| !positive.contains(item));
        if items.is_empty() {
            conclusions.remove("findings");
        }
    }
}

fn attempt_rule_outs(node: &Node) -> Vec<String> {
    let contract = node.attrs.get("original_contract").unwrap_or(&node.attrs);
    let mut values = Vec::new();
    for key in [
        "rules_out",
        "rules_out_free",
        "kills_free",
        "kills",
        "excludes",
        "negative_results",
    ] {
        let value = contract.get(key).or_else(|| node.attrs.get(key));
        match value {
            Some(Value::Array(items)) => {
                values.extend(items.iter().map(|item| clip(&value_text(Some(item)), 240)))
            }
            Some(Value::String(text)) if !text.trim().is_empty() => values.push(clip(text, 240)),
            _ => {}
        }
    }
    values.retain(|text| !text.is_empty());
    values.dedup();
    values.truncate(4);
    values
}

fn attempt_topic_score(
    attempt: &Node,
    topic: &Node,
    topic_claims: &HashSet<String>,
    graph: &KnowledgeGraph,
) -> (usize, usize, usize) {
    let mut tested_claims = HashSet::new();
    for key in ["claims_tested"] {
        if let Some(values) = attempt.attrs.get(key).and_then(Value::as_array) {
            tested_claims.extend(values.iter().filter_map(Value::as_str).map(str::to_owned));
        }
        if let Some(values) = attempt
            .attrs
            .pointer(&format!("/original_contract/{key}"))
            .and_then(Value::as_array)
        {
            tested_claims.extend(values.iter().filter_map(Value::as_str).map(str::to_owned));
        }
    }
    for edge in graph.edges.iter().filter(|edge| edge.kind == "tests") {
        let claim = if edge.src == attempt.id {
            Some(edge.dst.as_str())
        } else if edge.dst == attempt.id {
            Some(edge.src.as_str())
        } else {
            None
        };
        if let Some(claim) = claim.filter(|claim| topic_claims.contains(*claim)) {
            tested_claims.insert(claim.to_owned());
        }
    }
    let claim_overlap = tested_claims.intersection(topic_claims).count();

    let mut topic_text = format!("{} {}", topic.label, topic.summary);
    if let Some(keywords) = topic.attrs.get("keywords").and_then(Value::as_array) {
        for keyword in keywords.iter().filter_map(Value::as_str) {
            topic_text.push(' ');
            topic_text.push_str(keyword);
        }
    }
    for key in ["what_we_have", "what_we_seek"] {
        if let Some(values) = topic.attrs.get(key).and_then(Value::as_array) {
            for value in values {
                topic_text.push(' ');
                topic_text.push_str(&value_text(Some(value)));
            }
        }
    }
    let mut attempt_text = format!("{} {}", attempt.id, attempt.label);
    for key in [
        "question",
        "summary",
        "summary_de",
        "contract",
        "verdict_raw",
    ] {
        if let Some(value) = attempt.attrs.get(key) {
            attempt_text.push(' ');
            attempt_text.push_str(&value_text(Some(value)));
        }
    }
    if let Some(contract) = attempt.attrs.get("original_contract") {
        for key in ["contract", "question", "summary", "verdict_raw"] {
            if let Some(value) = contract.get(key) {
                attempt_text.push(' ');
                attempt_text.push_str(&value_text(Some(value)));
            }
        }
    }
    if let Some(values) = attempt.attrs.get("claims_tested").and_then(Value::as_array) {
        attempt_text.push(' ');
        attempt_text.push_str(
            &values
                .iter()
                .filter_map(Value::as_str)
                .collect::<Vec<_>>()
                .join(" "),
        );
    }
    let stop_words = [
        "the", "and", "from", "with", "does", "this", "that", "into", "under", "there", "their",
        "which", "der", "die", "das", "und", "von", "mit", "aus", "des", "den", "dem", "ein",
        "eine", "einer", "eines", "wie", "ist", "sind", "wird", "werden", "für", "zur", "zum",
        "über", "hinaus", "nicht", "auch", "einem",
    ];
    let topic_tokens: HashSet<String> = tokenize(&topic_text)
        .into_iter()
        .filter(|token| !stop_words.contains(&token.as_str()))
        .map(topic_token_stem)
        .collect();
    let attempt_tokens: HashSet<String> = tokenize(&attempt_text)
        .into_iter()
        .map(topic_token_stem)
        .collect();
    let keyword_overlap = topic_tokens.intersection(&attempt_tokens).count();
    let verdict = attempt_verdict(attempt)
        .unwrap_or_default()
        .to_ascii_uppercase();
    let verdict_priority = if matches!(
        verdict.as_str(),
        "PARTIAL" | "FAIL" | "FAILED" | "NEGATIVE" | "INCONCLUSIVE"
    ) {
        2
    } else if !verdict.is_empty() {
        1
    } else {
        0
    };
    (claim_overlap, keyword_overlap, verdict_priority)
}

fn topic_token_stem(token: String) -> String {
    if token.starts_with("gravit") {
        "gravit".to_owned()
    } else if token.starts_with("raumzeit") {
        "raumzeit".to_owned()
    } else {
        token
    }
}

fn verdict_is_negative(verdict: &str) -> bool {
    matches!(
        verdict.trim().to_ascii_uppercase().as_str(),
        "FAIL" | "FAILED" | "NEGATIVE" | "FALSIFIED" | "REFUTED" | "INCONCLUSIVE_NEGATIVE"
    )
}

fn goal_requirement(node: &Node) -> String {
    ["requirement", "text", "description", "summary", "label"]
        .into_iter()
        .find_map(|key| node.attrs.get(key).and_then(Value::as_str))
        .unwrap_or(&node.label)
        .to_owned()
}

fn impact_matches_risk(value: &Value) -> bool {
    match value {
        Value::String(value) => matches!(
            value.as_str(),
            "conditional" | "circularity_risk" | "missing_derivation"
        ),
        Value::Array(values) => values.iter().any(impact_matches_risk),
        _ => false,
    }
}

fn set_intersection_strings(sets: &[HashSet<String>]) -> Vec<String> {
    let Some(first) = sets.first() else {
        return Vec::new();
    };
    let mut result: Vec<String> = first
        .iter()
        .filter(|value| sets.iter().skip(1).all(|set| set.contains(*value)))
        .cloned()
        .collect();
    result.sort();
    result
}

fn theory_patterns(theory: &Value) -> Vec<Value> {
    theory
        .get("patterns")
        .and_then(Value::as_array)
        .or_else(|| theory.get("items").and_then(Value::as_array))
        .or_else(|| theory.get("records").and_then(Value::as_array))
        .or_else(|| theory.get("nodes").and_then(Value::as_array))
        .cloned()
        .unwrap_or_default()
}

fn compact_theory_map(map: &Value, limit: usize) -> Value {
    let patterns = map
        .get("patterns")
        .and_then(Value::as_array)
        .or_else(|| map.get("items").and_then(Value::as_array))
        .or_else(|| map.get("records").and_then(Value::as_array))
        .or_else(|| map.get("nodes").and_then(Value::as_array));
    let pattern_total = patterns.map_or(0, Vec::len);
    let items: Vec<Value> = patterns.into_iter().flatten().take(limit).map(|entry| {
        let source_ids: Vec<String> = source_ref_ids(entry).into_iter().take(4).collect();
        let related_ids: Vec<String> = entry.get("related_ids").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str).take(8).map(str::to_owned).collect();
        let compact_texts = |key: &str| -> Vec<Value> {
            entry.get(key).and_then(Value::as_array).into_iter().flatten().take(2).map(|value| {
                if value.is_object() {
                    json!({"label_de":clip(&value_text(value.get("label_de").or_else(||value.get("label"))),120),"label_en":clip(&value_text(value.get("label_en")),120)})
                } else {
                    json!({"label_de":clip(&value_text(Some(value)),120)})
                }
            }).collect()
        };
        json!({
            "id":entry.get("id").cloned().unwrap_or(Value::Null),
            "type":entry.get("type").cloned().unwrap_or_else(||json!("pattern")),
            "kind":entry.get("kind").cloned().unwrap_or(Value::Null),
            "label":clip(value_text(entry.get("label_de").or_else(||entry.get("title_de")).or_else(||entry.get("name_de")).or_else(||entry.get("label")).or_else(||entry.get("title")).or_else(||entry.get("name")).or_else(||entry.get("label_en"))).as_str(),140),
            "summary":clip(value_text(entry.get("summary_de").or_else(||entry.get("description_de")).or_else(||entry.get("statement_de")).or_else(||entry.get("summary")).or_else(||entry.get("description")).or_else(||entry.get("statement")).or_else(||entry.get("summary_en"))).as_str(),240),
            "scope":clip(value_text(entry.get("scope")).as_str(),120),
            "source_ids":source_ids,
            "related_ids":related_ids,
            "conditions":compact_texts("conditions"),
            "counterevidence":compact_texts("counterevidence")
        })
    }).collect();
    // Keep the complete authored map structure in the compact MCP response, but
    // omit source excerpts and long descriptions so the map remains context-safe.
    let nodes: Vec<Value> = map
        .get("nodes")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .map(|node| {
            json!({
                "id": node.get("id").cloned().unwrap_or(Value::Null),
                "label": clip(value_text(node.get("label_de").or_else(|| node.get("label")).or_else(|| node.get("label_en"))).as_str(), 180),
                "state": node.get("state").cloned().unwrap_or(Value::Null),
                "stage": node.get("stage").cloned().unwrap_or(Value::Null),
                "qualification": clip(value_text(node.get("overview_de").or_else(|| node.get("overview_en"))).as_str(), 180)
            })
        })
        .collect();
    let edges: Vec<Value> = map
        .get("edges")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .map(|edge| {
            json!({
                "from": edge.get("from").cloned().unwrap_or(Value::Null),
                "to": edge.get("to").cloned().unwrap_or(Value::Null),
                "kind": edge.get("kind").cloned().unwrap_or(Value::Null)
            })
        })
        .collect();
    let pattern_index: Vec<Value> = patterns.into_iter().flatten().map(|pattern| {
        json!({"id":pattern["id"], "label":clip(value_text(
            pattern.get("label_de").or_else(|| pattern.get("label")).or_else(|| pattern.get("label_en"))
        ).as_str(), 96)})
    }).collect();
    let mut compact = json!({
        "available": true,
        "counts": {"nodes": nodes.len(), "edges": edges.len(), "patterns": pattern_total},
        "nodes": nodes,
        "edges": edges,
        "items": items,
        "pattern_index": pattern_index,
        "evidence_areas": crate::overview::areas("de", true),
        "truncated": pattern_total > limit,
        "omitted_patterns": 0,
        "notice": "Read-only theory map inventory; record presence is not proof validation. The evidence areas cover the whole theory; look up pattern_index IDs with get_dossier for details."
    });
    // Keep every node/edge and each retained pattern intact. Long source IDs can
    // make a fixed item count exceed the overview's context budget.
    while serde_json::to_vec(&compact).map_or(false, |s| s.len() >= 11900)
        && !compact["items"].as_array().unwrap().is_empty()
    {
        compact["items"].as_array_mut().unwrap().pop();
    }
    let retained = compact["items"].as_array().unwrap().len();
    compact["truncated"] = json!(pattern_total > retained);
    compact["omitted_patterns"] = json!(pattern_total - retained);
    compact
}

fn public_node(node: &Node) -> Value {
    let mut value = serde_json::to_value(node).unwrap_or(Value::Null);
    if let Some(attrs) = value.get_mut("attrs").and_then(Value::as_object_mut) {
        attrs.remove("embedding");
        attrs.remove("embedding_model");
    }
    value
}

fn compact_dossier_projection(dossier: &Value, max_chars: usize, neighbor_limit: usize) -> Value {
    let full_node = dossier.get("node").unwrap_or(&Value::Null);
    let attrs = full_node.get("attrs").unwrap_or(&Value::Null);
    let semantic = attrs.get("semantic").unwrap_or(&Value::Null);
    let original_contract = attrs
        .get("original_contract")
        .or_else(|| attrs.get("contract"))
        .unwrap_or(&Value::Null);
    let contract_fields = [
        "question",
        "verdict",
        "verdict_raw",
        "scope",
        "positive_results",
        "analytic_conditional_results",
        "findings",
        "assumptions",
        "baseline_requirements",
        "required_baselines",
        "known_baselines",
        "baseline_requirement",
        "kills",
        "kills_free",
        "leaves_open",
        "missing",
        "not_derived",
        "limitations",
    ];
    let mut merged_contract = original_contract.as_object().cloned().unwrap_or_default();
    for key in contract_fields {
        let Some(outer) = attrs.get(key) else {
            continue;
        };
        let nested_is_empty = merged_contract.get(key).is_none_or(|value| {
            value.is_null()
                || value.as_str().is_some_and(str::is_empty)
                || value.as_array().is_some_and(Vec::is_empty)
        });
        let outer_is_populated = !outer.is_null()
            && !outer.as_str().is_some_and(str::is_empty)
            && !outer.as_array().is_some_and(Vec::is_empty);
        if !merged_contract.contains_key(key) || (nested_is_empty && outer_is_populated) {
            merged_contract.insert(key.to_owned(), outer.clone());
        }
    }
    let contract = Value::Object(merged_contract);
    let evidence_text_limit = if max_chars <= 12_000 {
        180
    } else if max_chars <= 24_000 {
        320
    } else {
        600
    };
    let evidence_items_limit = if max_chars <= 12_000 {
        4
    } else if max_chars <= 24_000 {
        6
    } else {
        8
    };
    let mut omissions = serde_json::Map::new();
    let mut copy_fields = |source: &Value,
                           names: &[&str],
                           cap: usize,
                           target: &mut serde_json::Map<String, Value>,
                           prefix: &str,
                           bucket: &str| {
        let mut nested_fields = 0usize;
        let mut nested_characters = 0usize;
        if let Some(fields) = source.as_object() {
            let retained = names
                .iter()
                .filter(|name| fields.contains_key(**name))
                .count();
            let omitted = fields.len().saturating_sub(retained);
            if omitted > 0 {
                omissions.insert(format!("{bucket}_fields"), json!(omitted));
            }
        }
        for name in names {
            let Some(value) = source.get(*name) else {
                continue;
            };
            let output_key = format!("{prefix}{name}");
            if let Some(items) = value.as_array() {
                let selected: Vec<Value> = items
                    .iter()
                    .take(cap.min(evidence_items_limit))
                    .map(|item| compact_evidence_item(item, evidence_text_limit))
                    .collect();
                for item in items.iter().take(cap.min(evidence_items_limit)) {
                    let (fields, characters) =
                        compact_evidence_omission_stats(item, evidence_text_limit);
                    nested_fields += fields;
                    nested_characters += characters;
                }
                target.insert(output_key.clone(), Value::Array(selected));
                if items.len() > cap.min(evidence_items_limit) {
                    omissions.insert(
                        output_key,
                        json!(items.len() - cap.min(evidence_items_limit)),
                    );
                }
            } else if let Some(text) = value.as_str() {
                let clipped = clip(text, evidence_text_limit.max(240));
                if clipped.chars().count() < text.chars().count() {
                    omissions.insert(
                        format!("{output_key}_characters"),
                        json!(text.chars().count() - clipped.chars().count()),
                    );
                }
                target.insert(output_key, json!(clipped));
            } else {
                let (fields, characters) = compact_evidence_omission_stats(value, 600);
                nested_fields += fields;
                nested_characters += characters;
                target.insert(output_key, value.clone());
            }
        }
        if nested_fields > 0 {
            omissions.insert(format!("{bucket}_nested_fields"), json!(nested_fields));
        }
        if nested_characters > 0 {
            omissions.insert(
                format!("{bucket}_nested_characters"),
                json!(nested_characters),
            );
        }
    };

    let mut compact_attrs = serde_json::Map::new();
    for key in [
        "assessment",
        "enrichment_state",
        "cache_freshness",
        "provenance",
        "extraction_quality",
        "translations",
    ] {
        if let Some(value) = attrs.get(key) {
            compact_attrs.insert(
                key.to_owned(),
                if key == "cache_freshness" {
                    value.clone()
                } else {
                    compact_evidence_item(value, 500)
                },
            );
        }
    }
    compact_attrs.insert(
        "ledger_status".into(),
        attrs
            .pointer("/ledger/status")
            .cloned()
            .unwrap_or(Value::Null),
    );
    compact_attrs.insert(
        "source_refs".into(),
        attrs
            .get("source_refs")
            .cloned()
            .unwrap_or_else(|| json!([])),
    );

    let mut semantic_projection = serde_json::Map::new();
    copy_fields(
        semantic,
        &[
            "overview_de",
            "overview_en",
            "interpretation_de",
            "interpretation_en",
        ],
        1,
        &mut semantic_projection,
        "",
        "semantic",
    );
    copy_fields(
        semantic,
        &[
            "results",
            "premises",
            "derivations",
            "proof_gaps",
            "input_audit",
            "unknowns",
        ],
        evidence_items_limit,
        &mut semantic_projection,
        "",
        "semantic",
    );
    for name in ["script_analysis", "input_audit", "proof_gaps"] {
        if let Some(value) = attrs.pointer(&format!("/static_inputs/{name}")) {
            semantic_projection.insert(format!("static_{name}"), compact_evidence_item(value, 450));
        }
    }
    compact_attrs.insert("semantic".into(), Value::Object(semantic_projection));

    let mut contract_projection = serde_json::Map::new();
    copy_fields(
        &contract,
        &[
            "question",
            "verdict",
            "verdict_raw",
            "scope",
            "summary_de",
            "original_source_scope_correction",
            "analytic_results",
            "claims_tested",
            "gates",
            "no_t1_t8_closure",
            "positive_results",
            "analytic_conditional_results",
            "findings",
            "assumptions",
            "baseline_requirements",
            "required_baselines",
            "known_baselines",
            "baseline_requirement",
            "kills",
            "kills_free",
            "leaves_open",
            "missing",
            "not_derived",
            "limitations",
        ],
        evidence_items_limit,
        &mut contract_projection,
        "",
        "original_contract",
    );
    drop(copy_fields);
    for key in [
        "assessment",
        "provenance",
        "extraction_quality",
        "translations",
    ] {
        if let Some(value) = attrs.get(key) {
            let (fields, characters) = compact_evidence_omission_stats(value, 500);
            if fields > 0 {
                omissions.insert(format!("{key}_nested_fields"), json!(fields));
            }
            if characters > 0 {
                omissions.insert(format!("{key}_nested_characters"), json!(characters));
            }
        }
    }
    for name in ["script_analysis", "input_audit", "proof_gaps"] {
        if let Some(value) = attrs.pointer(&format!("/static_inputs/{name}")) {
            let (fields, characters) = compact_evidence_omission_stats(value, 450);
            if fields > 0 {
                omissions.insert(format!("static_{name}_nested_fields"), json!(fields));
            }
            if characters > 0 {
                omissions.insert(
                    format!("static_{name}_nested_characters"),
                    json!(characters),
                );
            }
        }
    }
    if !contract_projection.is_empty() {
        compact_attrs.insert(
            "original_contract".into(),
            Value::Object(contract_projection),
        );
    }
    if !omissions.is_empty() {
        compact_attrs.insert("omitted_counts".into(), Value::Object(omissions.clone()));
    }

    let label = full_node.get("label").and_then(Value::as_str).unwrap_or("");
    let summary = full_node
        .get("summary")
        .and_then(Value::as_str)
        .unwrap_or("");
    let mut compact_node = json!({
        "id":full_node.get("id"),
        "type":full_node.get("type"),
        "label":clip(label,240),
        "summary":clip(summary,evidence_text_limit.max(300)),
        "attrs":compact_attrs,
        "omissions":{"label_characters":label.chars().count().saturating_sub(clip(label,240).chars().count()),"summary_characters":summary.chars().count().saturating_sub(clip(summary,evidence_text_limit.max(300)).chars().count())}
    });
    let all_neighbors = dossier
        .get("related_nodes")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let related_nodes: Vec<Value> = all_neighbors.iter().take(neighbor_limit).map(|node| json!({
        "id":node.get("id"),"type":node.get("type"),"label":clip(node.get("label").and_then(Value::as_str).unwrap_or(""),120),"summary":clip(node.get("summary").and_then(Value::as_str).unwrap_or(""),180)
    })).collect();
    let edges: Vec<Value> = dossier
        .get("edges")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .take(neighbor_limit.saturating_mul(2))
        .map(|edge| json!({"src":edge.get("src"),"dst":edge.get("dst"),"type":edge.get("type")}))
        .collect();
    let all_sources = dossier
        .get("sources")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let node_id = full_node.get("id").and_then(Value::as_str).unwrap_or("");
    let mut source_id_set: HashSet<String> = source_ref_ids(attrs).into_iter().collect();
    for edge in dossier
        .get("edges")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        if !is_source_type(edge.get("type").and_then(Value::as_str).unwrap_or("")) {
            continue;
        }
        let src = edge.get("src").and_then(Value::as_str).unwrap_or("");
        let dst = edge.get("dst").and_then(Value::as_str).unwrap_or("");
        if src == node_id {
            source_id_set.insert(dst.to_owned());
        } else if dst == node_id {
            source_id_set.insert(src.to_owned());
        }
    }
    let mut source_ids: Vec<String> = source_id_set.into_iter().collect();
    source_ids.sort();
    let source_limit = if max_chars <= 12_000 { 8 } else { 20 };
    let retained_source_ids: Vec<String> = source_ids.iter().take(source_limit).cloned().collect();
    let node_source_ids: Vec<Value> = source_ids.iter().cloned().map(Value::String).collect();
    compact_node["attrs"]["source_refs"] = json!(retained_source_ids);
    let mut matching_sources: Vec<&Value> = all_sources
        .iter()
        .filter(|source| {
            node_source_ids
                .iter()
                .any(|id| id.as_str() == source.get("id").and_then(Value::as_str))
        })
        .collect();
    matching_sources.sort_by(|a, b| {
        let priority = |source: &Value| {
            let path = source.get("path").and_then(Value::as_str).unwrap_or("");
            if path.ends_with("/PROOF.md") || path == "PROOF.md" {
                0
            } else if path.ends_with("/contract_index.json") || path == "contract_index.json" {
                1
            } else if path.ends_with("/checker.py") || path == "checker.py" {
                2
            } else if path.ends_with("/README.md") || path == "README.md" {
                3
            } else {
                4
            }
        };
        priority(a).cmp(&priority(b)).then_with(|| {
            a.get("path")
                .and_then(Value::as_str)
                .cmp(&b.get("path").and_then(Value::as_str))
        })
    });
    let own_sources: Vec<Value> = matching_sources
        .into_iter()
        .take(source_limit)
        .map(|source| json!({"id":source.get("id"),"path":source.get("path"),"sha256":source.get("sha256"),"line_start":source.get("line_start"),"line_end":source.get("line_end"),"page":source.get("page")}))
        .collect();
    let count = |name: &str| {
        dossier
            .pointer(&format!("/counts/{name}"))
            .cloned()
            .unwrap_or(Value::Null)
    };
    let registered_own_sources = all_sources
        .iter()
        .filter(|source| {
            node_source_ids
                .iter()
                .any(|id| id.as_str() == source.get("id").and_then(Value::as_str))
        })
        .count();
    let result = json!({
        "node":compact_node,
        "related_nodes":related_nodes,
        "edges":edges,
        "sources":own_sources,
        "external_sources":dossier.get("external_sources").cloned().unwrap_or_else(||json!([])),
        "counts":{"related_nodes":count("related_nodes"),"edges":count("edges"),"sources":count("sources"),"own_sources_registered":registered_own_sources},
        "omissions":{"related_nodes":all_neighbors.len().saturating_sub(neighbor_limit),"edges":dossier.get("edges").and_then(Value::as_array).map_or(0,|items|items.len().saturating_sub(neighbor_limit.saturating_mul(2))),"own_sources":registered_own_sources.saturating_sub(own_sources.len()),"source_refs":source_ids.len().saturating_sub(source_limit),"source_refs_without_registered_path":node_source_ids.len().saturating_sub(registered_own_sources),"node_fields":omissions},
        "notice":"This MCP dossier is a bounded projection. The HTTP dossier retains full attributes; source IDs and paths are shown for targeted read_source calls. Omitted items are counted explicitly."
    });
    if json_char_count(&result) > max_chars {
        // The fixed projection stays useful under the minimum budget; shorten
        // only its display samples and report every removed collection item.
        let mut result = result;
        result["related_nodes"] = json!([]);
        result["edges"] = json!([]);
        result["omissions"]["related_nodes"] = json!(all_neighbors.len());
        result["omissions"]["edges"] = dossier
            .get("edges")
            .and_then(Value::as_array)
            .map_or(json!(0), |items| json!(items.len()));
        let retained_sources = result["sources"].as_array().map(Vec::len).unwrap_or(0);
        result["omissions"]["own_sources"] =
            json!(registered_own_sources.saturating_sub(retained_sources));
        let arrays = [
            ("/node/attrs/semantic/derivations", 0usize),
            ("/node/attrs/semantic/input_audit", 0),
            ("/node/attrs/semantic/proof_gaps", 0),
            ("/node/attrs/original_contract/analytic_results", 0),
            ("/node/attrs/original_contract/findings", 0),
            ("/node/attrs/original_contract/limitations", 0),
            ("/node/attrs/original_contract/assumptions", 1),
            ("/node/attrs/original_contract/positive_results", 1),
            ("/node/attrs/original_contract/leaves_open", 1),
            ("/node/attrs/semantic/results", 1),
            ("/node/attrs/semantic/premises", 1),
            ("/node/attrs/source_refs", 2),
            ("/sources", 1),
        ];
        let mut budget_items = serde_json::Map::new();
        while json_char_count(&result) > max_chars {
            let candidate = arrays
                .iter()
                .filter_map(|(path, minimum)| {
                    result
                        .pointer(path)
                        .and_then(Value::as_array)
                        .filter(|items| items.len() > *minimum)
                        .map(|items| (*path, json_char_count(items.last().unwrap())))
                })
                .max_by_key(|(_, size)| *size);
            let Some((path, _)) = candidate else { break };
            if let Some(items) = result.pointer_mut(path).and_then(Value::as_array_mut) {
                items.pop();
                let count = budget_items.get(path).and_then(Value::as_u64).unwrap_or(0) + 1;
                budget_items.insert(path.to_owned(), json!(count));
            }
        }
        if !budget_items.is_empty() {
            result["omissions"]["budget_items"] = Value::Object(budget_items);
        }
        let retained_sources = result["sources"].as_array().map_or(0, Vec::len);
        let retained_source_refs = result
            .pointer("/node/attrs/source_refs")
            .and_then(Value::as_array)
            .map_or(0, Vec::len);
        result["omissions"]["own_sources"] =
            json!(registered_own_sources.saturating_sub(retained_sources));
        result["omissions"]["source_refs"] =
            json!(source_ids.len().saturating_sub(retained_source_refs));
        for limit in [180usize, 120, 80, 50] {
            if json_char_count(&result) <= max_chars {
                break;
            }
            let removed = shrink_dossier_text(&mut result, limit);
            if removed > 0 {
                result["omissions"]["budget_text_characters"] = json!(
                    result["omissions"]["budget_text_characters"]
                        .as_u64()
                        .unwrap_or(0)
                        + removed as u64
                );
            }
        }
        result
    } else {
        result
    }
}

fn shrink_dossier_text(value: &mut Value, limit: usize) -> usize {
    match value {
        Value::String(text) => {
            let previous = text.chars().count();
            if previous > limit {
                *text = clip(text, limit);
            }
            previous.saturating_sub(text.chars().count())
        }
        Value::Array(items) => items
            .iter_mut()
            .map(|item| shrink_dossier_text(item, limit))
            .sum(),
        Value::Object(fields) => {
            let mut removed = 0usize;
            for (key, item) in fields.iter_mut() {
                if matches!(
                    key.as_str(),
                    "id" | "source_id"
                        | "path"
                        | "sha256"
                        | "type"
                        | "kind"
                        | "category"
                        | "authority"
                        | "status"
                        | "verdict"
                        | "origin"
                        | "enrichment_state"
                ) {
                    continue;
                }
                if let Value::String(text) = item {
                    let previous = text.chars().count();
                    if previous > limit {
                        *text = clip(text, limit);
                        removed += previous.saturating_sub(text.chars().count());
                    }
                } else {
                    removed += shrink_dossier_text(item, limit);
                }
            }
            removed
        }
        _ => 0,
    }
}

fn compact_evidence_item(value: &Value, max_text: usize) -> Value {
    match value {
        Value::String(text) => json!(clip(text, max_text)),
        Value::Array(items) => Value::Array(
            items
                .iter()
                .take(8)
                .map(|item| compact_evidence_item(item, max_text))
                .collect(),
        ),
        Value::Object(object) => {
            let mut result = serde_json::Map::new();
            for key in [
                "id",
                "label",
                "name",
                "category",
                "authority",
                "status",
                "notice",
                "text",
                "statement",
                "epistemic",
                "value",
                "interpretation",
                "scope",
                "limitations",
                "kind",
                "origin",
                "impact",
                "used_for",
                "missing_step",
                "why_checks_insufficient",
                "premise_ids",
                "upstream_claim_ids",
                "result_ids",
                "source_refs",
            ] {
                if let Some(item) = object.get(key) {
                    let item = if let Some(text) = item.as_str() {
                        json!(clip(text, max_text))
                    } else if item.is_array() {
                        let array_limit = if key == "source_refs" { 4 } else { 12 };
                        Value::Array(
                            item.as_array()
                                .unwrap()
                                .iter()
                                .take(array_limit)
                                .map(|entry| compact_evidence_item(entry, 120))
                                .collect(),
                        )
                    } else {
                        item.clone()
                    };
                    result.insert(key.to_owned(), item);
                }
            }
            Value::Object(result)
        }
        _ => value.clone(),
    }
}

fn compact_evidence_omission_stats(value: &Value, max_text: usize) -> (usize, usize) {
    const KEYS: &[&str] = &[
        "id",
        "label",
        "name",
        "category",
        "authority",
        "status",
        "notice",
        "text",
        "statement",
        "epistemic",
        "value",
        "interpretation",
        "scope",
        "limitations",
        "kind",
        "origin",
        "impact",
        "used_for",
        "missing_step",
        "why_checks_insufficient",
        "premise_ids",
        "upstream_claim_ids",
        "result_ids",
        "source_refs",
    ];
    match value {
        Value::String(text) => (0, text.chars().count().saturating_sub(max_text)),
        Value::Array(items) => {
            let mut fields = items.len().saturating_sub(12);
            let mut chars = 0usize;
            for item in items.iter().take(12) {
                let (nested_fields, nested_chars) = compact_evidence_omission_stats(item, max_text);
                fields += nested_fields;
                chars += nested_chars;
            }
            (fields, chars)
        }
        Value::Object(object) => {
            let mut fields = object
                .keys()
                .filter(|key| !KEYS.contains(&key.as_str()))
                .count();
            let mut chars = 0usize;
            for (key, item) in object
                .iter()
                .filter(|(key, _)| KEYS.contains(&key.as_str()))
            {
                let child_limit = if matches!(
                    key.as_str(),
                    "source_refs" | "premise_ids" | "result_ids" | "upstream_claim_ids"
                ) {
                    120
                } else {
                    max_text
                };
                let (nested_fields, nested_chars) =
                    compact_evidence_omission_stats(item, child_limit);
                fields += nested_fields;
                chars += nested_chars;
            }
            (fields, chars)
        }
        _ => (0, 0),
    }
}

fn compact_theme_value(value: &Value, depth: usize) -> Value {
    if depth >= 5 {
        return Value::Null;
    }
    match value {
        Value::String(text) => json!(clip(text, 280)),
        Value::Array(items) => Value::Array(
            items
                .iter()
                .take(12)
                .map(|item| compact_theme_value(item, depth + 1))
                .collect(),
        ),
        Value::Object(object) => {
            let mut result = serde_json::Map::new();
            for (key, item) in object.iter().take(20) {
                result.insert(key.clone(), compact_theme_value(item, depth + 1));
            }
            Value::Object(result)
        }
        _ => value.clone(),
    }
}

/// Theme pages are index cards. The dossier remains the full-record endpoint.
/// This projection keeps visible status, bilingual summary, and source links,
/// while avoiding repetition of the full ledger and AI annotation in each card.
fn compact_theme_node(node: &Node, topic_id: &str) -> Value {
    let mut attrs = serde_json::Map::new();
    let source_ids = node_source_ids(node);
    let assessment = node.attrs.get("assessment").and_then(Value::as_object);
    if let Some(assessment) = assessment {
        let mut compact = serde_json::Map::new();
        for key in [
            "category",
            "label",
            "label_en",
            "authority",
            "scope",
            "scope_en",
            "notice",
            "notice_en",
            "classification_method",
            "classification_uncertainty",
        ] {
            if let Some(value) = assessment.get(key) {
                let value = if matches!(key, "scope" | "scope_en") {
                    json!(clip(&value_text(Some(value)), 320))
                } else if matches!(key, "notice" | "notice_en") {
                    json!(clip(&value_text(Some(value)), 200))
                } else {
                    compact_theme_value(value, 0)
                };
                compact.insert(key.to_owned(), value);
            }
        }
        attrs.insert("assessment".to_owned(), Value::Object(compact));
    }
    if let Some(status) = node.attrs.pointer("/ledger/status") {
        attrs.insert("ledger_status".to_owned(), compact_theme_value(status, 0));
    }
    let display_fields: &[&str] = if node.kind == "claim" {
        &["enrichment_state"]
    } else {
        &[
            "enrichment_state",
            "label_en",
            "summary_en",
            "required_steps_en",
        ]
    };
    for key in display_fields {
        if let Some(value) = node.attrs.get(*key) {
            attrs.insert((*key).to_owned(), compact_theme_value(value, 0));
        }
    }
    if let Some(roles) = node.attrs.get("topic_roles") {
        if let Some(role) = roles.get(topic_id) {
            attrs.insert("topic_role".to_owned(), role.clone());
        }
    }
    let source_refs = node
        .attrs
        .get("source_refs")
        .or_else(|| node.attrs.pointer("/semantic/source_refs"))
        .or_else(|| node.attrs.pointer("/original_contract/source_refs"));
    if let Some(source_refs) = source_refs {
        let compact_refs = if let Some(items) = source_refs.as_array() {
            Value::Array(
                items
                    .iter()
                    .take(6)
                    .map(|item| compact_theme_value(item, 0))
                    .collect(),
            )
        } else {
            compact_theme_value(source_refs, 0)
        };
        attrs.insert("source_refs".to_owned(), compact_refs);
    } else if !source_ids.is_empty() {
        attrs.insert(
            "source_refs".to_owned(),
            Value::Array(
                source_ids
                    .iter()
                    .take(6)
                    .cloned()
                    .map(Value::String)
                    .collect(),
            ),
        );
    }
    if source_ids.len() > 6 {
        attrs.insert(
            "source_refs_omitted_count".to_owned(),
            json!(source_ids.len() - 6),
        );
    }
    if let Some(provenance) = node.attrs.get("provenance") {
        if let Some(reviewed) = provenance.get("reviewed") {
            attrs.insert("provenance".to_owned(), json!({"reviewed":reviewed}));
        }
    }
    let mut semantic = serde_json::Map::new();
    if let Some(source_semantic) = node.attrs.get("semantic").and_then(Value::as_object) {
        let semantic_keys: &[&str] = if node.kind == "claim" {
            &[
                "title_de",
                "title_en",
                "overview_de",
                "overview_en",
                "source_refs",
            ]
        } else {
            &[
                "title_de",
                "title_en",
                "overview_de",
                "overview_en",
                "interpretation_de",
                "interpretation_en",
                "statement_en",
                "source_refs",
            ]
        };
        for key in semantic_keys {
            if let Some(value) = source_semantic.get(*key) {
                semantic.insert((*key).to_owned(), compact_theme_value(value, 0));
            }
        }
    }
    if !semantic.is_empty() {
        attrs.insert("semantic".to_owned(), Value::Object(semantic));
    }
    if let Some(translations) = node.attrs.pointer("/translations/en") {
        let mut english = serde_json::Map::new();
        let have_title = node.attrs.pointer("/semantic/title_en").is_some()
            || node.attrs.get("label_en").is_some();
        let have_overview = node.attrs.pointer("/semantic/overview_en").is_some()
            || node.attrs.get("summary_en").is_some()
            || node.attrs.get("overview_en").is_some();
        for key in ["title", "overview"] {
            if (key == "title" && have_title) || (key == "overview" && have_overview) {
                continue;
            }
            if let Some(value) = translations.get(key) {
                english.insert(key.to_owned(), compact_theme_value(value, 0));
            }
        }
        if !english.is_empty() {
            attrs.insert("translations".to_owned(), json!({"en":english}));
        }
    }
    if let Some(orientation) = node
        .attrs
        .get("topic_orientation")
        .and_then(Value::as_object)
    {
        let mut compact = serde_json::Map::new();
        for key in [
            "role_de",
            "role_en",
            "have_de",
            "have_en",
            "seek_de",
            "seek_en",
            "connection_de",
            "connection_en",
            "next_bridge_de",
            "next_bridge_en",
        ] {
            if let Some(value) = orientation.get(key) {
                compact.insert(key.to_owned(), compact_theme_value(value, 0));
            }
        }
        attrs.insert("topic_orientation".to_owned(), Value::Object(compact));
    }
    if matches!(node.kind.as_str(), "attempt" | "contract") {
        for key in ["verdict", "verdict_raw", "status", "scope"] {
            if let Some(value) = node
                .attrs
                .get(key)
                .or_else(|| node.attrs.pointer(&format!("/original_contract/{key}")))
            {
                attrs.insert(key.to_owned(), compact_theme_value(value, 0));
            }
        }
        if let Some(contract) = node
            .attrs
            .get("original_contract")
            .and_then(Value::as_object)
        {
            let mut compact = serde_json::Map::new();
            for key in [
                "question",
                "scope",
                "verdict",
                "verdict_raw",
                "status",
                "findings",
                "positive_results",
                "conditional_results",
                "rules_out",
                "leaves_open",
                "not_derived",
                "limitations",
            ] {
                if let Some(value) = contract.get(key) {
                    compact.insert(key.to_owned(), compact_theme_value(value, 0));
                }
            }
            attrs.insert("original_contract".to_owned(), Value::Object(compact));
        }
    }
    json!({
        "id":node.id,
        "type":node.kind,
        "label":clip(&node.label, 180),
        "summary":clip(&node.summary, 700),
        "attrs":attrs
    })
}

/// A chain is an orientation aid, so retain the information needed to inspect
/// its recorded links while keeping large original ledgers out of every node.
/// Full attrs remain available through get_dossier.
fn chain_public_node(node: &Node) -> Value {
    const MAX_TEXT: usize = 1_200;
    const MAX_ARRAY: usize = 40;
    let compact_text = |value: &Value| -> Value {
        match value.as_str() {
            Some(text) if text.chars().count() > MAX_TEXT => Value::String(format!(
                "{}…",
                text.chars().take(MAX_TEXT).collect::<String>()
            )),
            Some(_) => value.clone(),
            None => value.clone(),
        }
    };
    let keep_keys: &[&str] = match node.kind.as_str() {
        "research_problem" | "research_finding" | "research_task" | "research_premise"
        | "research_strategy" => &[
            "status",
            "effective_status",
            "outcome",
            "scope_de",
            "scope_en",
            "source_state",
            "source_refs",
            "premise_ids",
            "depends_on",
            "resolution_ids",
            "provenance",
            "label_en",
            "summary_en",
        ],
        "claim" => &[
            "assessment",
            "ledger_status",
            "enrichment_state",
            "source_refs",
            "provenance",
            "extraction_quality",
            "translations",
        ],
        "assumption" => &[
            "text",
            "kind",
            "origin_claim_id",
            "source_refs",
            "provenance",
            "scope",
            "limitations",
        ],
        "derivation" => &[
            "argument",
            "scope",
            "premise_ids",
            "upstream_claim_ids",
            "result_ids",
            "source_refs",
            "provenance",
            "limitations",
        ],
        "result" => &[
            "statement",
            "epistemic",
            "value",
            "interpretation",
            "scope",
            "limitations",
            "source_refs",
            "provenance",
        ],
        "check" => &[
            "verdict",
            "status",
            "scope",
            "limitations",
            "result_ids",
            "source_refs",
            "provenance",
        ],
        "script" => &[
            "purpose_de",
            "inputs",
            "operations",
            "checks",
            "outputs",
            "limitations",
            "source_refs",
            "provenance",
        ],
        "attempt" => &[
            "verdict",
            "verdict_raw",
            "status",
            "assessment",
            "scope",
            "rules_out",
            "leaves_open",
            "not_derived",
            "limitations",
            "source_refs",
            "provenance",
        ],
        _ => &[
            "source_refs",
            "scope",
            "limitations",
            "provenance",
            "kind",
            "relation",
        ],
    };
    let mut attrs = serde_json::Map::new();
    if let Some(source_attrs) = node.attrs.as_object() {
        for key in keep_keys {
            if let Some(value) = source_attrs.get(*key) {
                let preserve_exactly = [
                    "text",
                    "argument",
                    "statement",
                    "scope",
                    "limitations",
                    "source_refs",
                    "premise_ids",
                    "upstream_claim_ids",
                    "result_ids",
                    "rules_out",
                    "leaves_open",
                    "not_derived",
                ]
                .contains(key);
                let value = if preserve_exactly {
                    value.clone()
                } else if let Some(items) = value.as_array() {
                    Value::Array(items.iter().take(MAX_ARRAY).map(&compact_text).collect())
                } else if let Some(obj) = value.as_object() {
                    // Assessment and provenance are structured, compact metadata.
                    let mut out = serde_json::Map::new();
                    for (nested_key, nested_value) in obj {
                        if [
                            "category",
                            "label",
                            "scope",
                            "authority",
                            "notice",
                            "reviewed",
                            "method",
                            "state",
                            "issues",
                            "source_refs",
                        ]
                        .contains(&nested_key.as_str())
                        {
                            let nested_value =
                                if nested_key == "scope" || nested_key == "source_refs" {
                                    nested_value.clone()
                                } else {
                                    compact_text(nested_value)
                                };
                            out.insert(nested_key.clone(), nested_value);
                        }
                    }
                    Value::Object(out)
                } else {
                    compact_text(value)
                };
                attrs.insert((*key).to_owned(), value);
            }
        }
    }
    json!({
        "id": node.id,
        "type": node.kind,
        "label": compact_text(&Value::String(node.label.clone())),
        "summary": compact_text(&Value::String(node.summary.clone())),
        "attrs": attrs
    })
}

fn theory_pattern_node(pattern: &Value) -> Result<Node> {
    let id = pattern
        .get("id")
        .and_then(Value::as_str)
        .context("theory pattern is missing its id")?;
    let source_ids: Vec<String> = pattern
        .get("source_refs")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(theory_pattern_source)
        .map(|source| source.id)
        .collect();
    Ok(Node {
        id: id.to_owned(),
        kind: "theory_pattern".to_owned(),
        label: value_text(pattern.get("label_de")),
        summary: value_text(pattern.get("summary_de")),
        attrs: json!({
            "kind":pattern.get("kind"),
            "label_en":pattern.get("label_en"),
            "summary_en":pattern.get("summary_en"),
            "reviewed":pattern.get("reviewed").and_then(Value::as_bool).unwrap_or(false),
            "epistemic_role":"authored_cross_source_pattern",
            "source_refs":source_ids,
            "related_ids":pattern.get("related_ids"),
            "conditions":pattern.get("conditions"),
            "counterevidence":pattern.get("counterevidence")
        }),
    })
}

fn theory_pattern_source(reference: &Value) -> Option<Source> {
    // External boundary references must remain explicit URLs, never local
    // read_source entries pointing at removed TFPT paths.
    if reference
        .get("source_location_state")
        .and_then(Value::as_str)
        == Some("moved_to_external_repository")
    {
        return None;
    }
    let id = reference
        .get("source_id")
        .and_then(Value::as_str)
        .or_else(|| reference.get("id").and_then(Value::as_str))
        .or_else(|| reference.as_str())?;
    Some(Source {
        id: id.to_owned(),
        path: reference.get("path").and_then(Value::as_str)?.to_owned(),
        sha256: reference.get("sha256").and_then(Value::as_str)?.to_owned(),
        line_start: reference
            .get("line_start")
            .and_then(Value::as_u64)
            .map(|line| line as usize),
        line_end: reference
            .get("line_end")
            .and_then(Value::as_u64)
            .map(|line| line as usize),
        line_count: reference
            .get("line_count")
            .and_then(Value::as_u64)
            .map(|lines| lines as usize),
        page: reference
            .get("page")
            .and_then(Value::as_u64)
            .map(|page| page as usize),
        text: reference
            .get("text")
            .or_else(|| reference.get("excerpt"))
            .and_then(Value::as_str)
            .map(str::to_owned),
    })
}
fn question_for_attempt<'a>(
    attempt: &Node,
    question_parent_by_id: &'a HashMap<&str, &Node>,
    nodes: &'a [Node],
) -> Option<&'a Node> {
    question_parent_by_id
        .iter()
        .find_map(|(question_id, parent)| {
            (parent.id == attempt.id)
                .then(|| nodes.iter().find(|node| node.id == *question_id))
                .flatten()
        })
}

fn compact_search_result(node: &Node, score: f32, matching_question: Option<&Node>) -> Value {
    let semantic = node.attrs.get("semantic");
    let counts = |key: &str| {
        semantic
            .and_then(|value| value.get(key))
            .and_then(Value::as_array)
            .map_or(0, Vec::len)
    };
    let summary = if node.summary.trim().is_empty() {
        value_text(
            node.attrs
                .pointer("/ledger/claim")
                .or_else(|| node.attrs.pointer("/original_contract/question"))
                .or_else(|| node.attrs.get("question")),
        )
    } else {
        node.summary.clone()
    };
    let translation = node.attrs.pointer("/translations/en");
    let label_en = node
        .attrs
        .get("label_en")
        .or_else(|| translation.and_then(|value| value.get("title")))
        .or_else(|| semantic.and_then(|value| value.get("title_en")));
    let summary_en = node
        .attrs
        .get("summary_en")
        .or_else(|| translation.and_then(|value| value.get("overview")))
        .or_else(|| semantic.and_then(|value| value.get("summary_en")))
        .or_else(|| node.attrs.pointer("/topic_orientation/role_en"));
    let overview_en = translation
        .and_then(|value| value.get("overview"))
        .or_else(|| semantic.and_then(|value| value.get("overview_en")))
        .or_else(|| node.attrs.pointer("/topic_orientation/role_en"));
    let orientation_en = node.attrs.pointer("/topic_orientation/role_en");
    let sources = node_source_ids(node);
    let source_total = sources.len();
    let mut result = json!({
        "id":node.id,
        "type":node.kind,
        "label":clip(&node.label,180),
        "label_en":label_en.and_then(Value::as_str).map(|text|clip(text,180)),
        "summary":clip(&summary,349),
        "summary_en":summary_en.and_then(Value::as_str).map(|text|clip(text,349)),
        "overview":clip(&value_text(semantic.and_then(|value| value.get("overview_de"))),240),
        "overview_en":overview_en.and_then(Value::as_str).map(|text|clip(text,240)),
        "orientation_en":orientation_en.and_then(Value::as_str).map(|text|clip(text,240)),
        "assessment":node.attrs.get("assessment").cloned().unwrap_or(Value::Null),
        "ledger_status":node.attrs.pointer("/ledger/status").cloned().unwrap_or(Value::Null),
        "enrichment_state":node.attrs.get("enrichment_state").cloned().unwrap_or(Value::Null),
        "source_ids":sources.into_iter().take(12).collect::<Vec<_>>(),
        "omissions":{"results":counts("results"),"proof_gaps":counts("proof_gaps"),"derivations":counts("derivations"),"input_audit":counts("input_audit"),"source_ids_truncated":source_total.saturating_sub(12),"detail":"Use get_dossier or audit_inputs for full recorded fields."},
        "publication_scope":node.attrs.get("publication_scope"),"content_scope":node.attrs.get("content_scope"),"summary_state":node.attrs.get("summary_state"),"source_current":node.attrs.get("source_current"),"source_freshness":node.attrs.get("source_freshness"),"context_current":node.attrs.get("context_current"),"reviewed":node.attrs.get("reviewed"),
        "score":score
    });
    if crate::papers::is_paper(node) {
        let current = crate::papers::summary_current(node, crate::papers::source_current(node));
        result["summary_current"] = json!(current);
        if !current {
            for field in ["summary", "summary_en", "overview", "overview_en"] {
                result[field] = Value::Null;
            }
        }
    }
    if let Some(question) = matching_question {
        result["contract_question_match"] = json!({
            "id":question.id,
            "label":clip(&question.label,180),
            "why_retained":"The related contract-question record matched this attempt's contract ID; it is shown here instead of as a duplicate search hit."
        });
    }
    result
}
fn absolute_path(path: PathBuf) -> Result<PathBuf> {
    if path.is_absolute() {
        Ok(path)
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}
fn tokenize(q: &str) -> Vec<String> {
    q.to_lowercase()
        .split(|c: char| !c.is_alphanumeric())
        .filter(|s| s.len() > 1)
        .map(str::to_owned)
        .collect()
}
fn count_impact(
    value: &str,
    conditional: &mut usize,
    circularity_risk: &mut usize,
    missing_derivation: &mut usize,
) {
    match value {
        "conditional" => *conditional += 1,
        "circularity_risk" => *circularity_risk += 1,
        "missing_derivation" => *missing_derivation += 1,
        _ => {}
    }
}
fn attempt_matches_status(node: &Node, requested: &str) -> bool {
    [
        node.attrs.get("status").and_then(Value::as_str),
        node.attrs.get("state").and_then(Value::as_str),
        node.attrs.get("verdict").and_then(Value::as_str),
        node.attrs.get("verdict_raw").and_then(Value::as_str),
        node.attrs
            .pointer("/assessment/status")
            .and_then(Value::as_str),
        node.attrs
            .pointer("/original_contract/verdict")
            .and_then(Value::as_str),
        node.attrs
            .pointer("/original_contract/verdict_raw")
            .and_then(Value::as_str),
    ]
    .into_iter()
    .flatten()
    .any(|value| value.eq_ignore_ascii_case(requested))
}
fn value_count(value: Option<&Value>) -> usize {
    match value {
        Some(Value::Array(items)) => items.len(),
        Some(Value::Object(items)) => items.len(),
        Some(Value::String(text)) if !text.trim().is_empty() => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(0) as usize,
        _ => 0,
    }
}
pub(crate) fn embedding(n: &Node) -> Option<Vec<f32>> {
    let v: Vec<f32> = n
        .attrs
        .get("embedding")?
        .as_array()?
        .iter()
        .filter_map(Value::as_f64)
        .map(|x| x as f32)
        .collect();
    if valid_embedding_vector(&v) {
        Some(v)
    } else {
        None
    }
}
pub(crate) fn cosine(a: &[f32], b: &[f32]) -> f32 {
    if a.len() != b.len() || a.is_empty() {
        return 0.0;
    }
    let (mut dot, mut aa, mut bb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        dot += x * y;
        aa += x * x;
        bb += y * y;
    }
    if aa == 0.0 || bb == 0.0 {
        0.0
    } else {
        dot / (aa.sqrt() * bb.sqrt())
    }
}
fn is_source_type(t: &str) -> bool {
    t.to_lowercase().contains("source")
        || t.to_lowercase().contains("cite")
        || t.to_lowercase().contains("evidence")
}
fn is_context_node(kind: &str) -> bool {
    matches!(kind, "topic" | "structure" | "goal")
}
fn bound_text(t: String) -> String {
    if t.len() > 200_000 {
        let mut end = 200_000;
        while !t.is_char_boundary(end) {
            end -= 1;
        }
        format!("{}\n[Excerpt truncated at 200 KB]", &t[..end])
    } else {
        t
    }
}
pub(crate) fn sensitive_path(p: &Path) -> bool {
    p.components().any(|c| {
        let s = c.as_os_str().to_string_lossy().to_lowercase();
        s.starts_with(".env")
            || s.contains("secret")
            || s.contains("credential")
            || s.ends_with(".pem")
            || s.ends_with(".key")
            || s == ".ssh"
    })
}
fn api_key() -> Result<String> {
    if let Ok(path) = std::env::var("OPENAI_API_KEY_FILE") {
        let contents = fs::read_to_string(path)?;
        return key_from_file_contents(&contents)
            .context("configured key file contains no usable API key");
    }
    if let Ok(key) = std::env::var("OPENAI_API_KEY") {
        if !key.trim().is_empty() {
            return Ok(key.trim().to_owned());
        }
    }
    if let Ok(text) = fs::read_to_string(std::env::current_dir()?.join(".env")) {
        if let Some(key) = dotenv_value(&text, "OPENAI_API_KEY") {
            return Ok(key);
        }
    }
    bail!("OpenAI API key is not configured")
}

fn key_from_file_contents(contents: &str) -> Option<String> {
    dotenv_value(contents, "OPENAI_API_KEY").or_else(|| {
        let raw = contents.trim();
        let looks_like_env_assignment = raw.split_once('=').is_some_and(|(name, _)| {
            !name.trim().is_empty()
                && name
                    .trim()
                    .chars()
                    .all(|ch| ch.is_ascii_uppercase() || ch.is_ascii_digit() || ch == '_')
        });
        (!raw.is_empty()
            && contents.lines().count() == 1
            && !raw.contains('\r')
            && !looks_like_env_assignment)
            .then(|| raw.to_owned())
    })
}

fn dotenv_value(contents: &str, key_name: &str) -> Option<String> {
    contents.lines().find_map(|line| {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return None;
        }
        let (name, value) = line.split_once('=')?;
        if name.trim() != key_name {
            return None;
        }
        let value = value.trim().trim_matches(['"', '\'']).trim();
        (!value.is_empty()).then(|| value.to_owned())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lexical_tokenization_is_case_insensitive() {
        assert_eq!(tokenize("Fixed POINT"), vec!["fixed", "point"]);
    }
    #[test]
    fn sensitive_paths_are_rejected() {
        assert!(sensitive_path(Path::new("/repo/.env")));
        assert!(sensitive_path(Path::new("/repo/private.pem")));
        assert!(!sensitive_path(Path::new("/repo/claims.csv")));
    }
    #[test]
    fn configured_key_file_accepts_dotenv_or_single_line_key_without_echoing_it() {
        assert_eq!(
            key_from_file_contents("# local config\nOTHER=value\nOPENAI_API_KEY=secret-value\n"),
            Some("secret-value".to_owned())
        );
        assert_eq!(
            key_from_file_contents("single-line-secret\n"),
            Some("single-line-secret".to_owned())
        );
        assert_eq!(key_from_file_contents("OTHER=value\n"), None);
    }
    #[test]
    fn public_nodes_omit_large_embedding_vectors() {
        let node: Node = serde_json::from_value(json!({"id":"claim:x","type":"claim","attrs":{"semantic":{"overview_de":"kept"},"embedding":vec![0.0;256],"embedding_model":"text-embedding-3-small"}})).unwrap();
        let public = public_node(&node);
        assert!(public.pointer("/attrs/embedding").is_none());
        assert!(public.pointer("/attrs/embedding_model").is_none());
        assert_eq!(
            public.pointer("/attrs/semantic/overview_de").unwrap(),
            "kept"
        );
    }
    #[test]
    fn attempt_card_uses_populated_outer_open_items_when_nested_contract_has_empty_placeholder() {
        let node: Node = serde_json::from_value(json!({
            "id":"contract:source-bridge",
            "type":"attempt",
            "label":"Source bridge",
            "attrs":{
                "original_contract":{"leaves_open":[],"missing":[]},
                "leaves_open":["The geometric horizon action is not source-derived."],
                "missing":[]
            }
        }))
        .unwrap();
        let card = attempt_card(&node, &[]);
        assert_eq!(
            card.pointer("/recorded_conclusions/missing_or_open/0"),
            Some(&json!(
                "The geometric horizon action is not source-derived."
            ))
        );
        assert_eq!(
            card.pointer("/recorded_conclusions/missing_or_open_recorded"),
            Some(&json!(true))
        );
    }
    #[test]
    fn attempt_card_marks_absent_open_fields_as_unrecorded() {
        let node: Node = serde_json::from_value(json!({
            "id":"contract:source-only",
            "type":"attempt",
            "label":"Source only",
            "attrs":{"original_contract":{"verdict":"PARTIAL"}}
        }))
        .unwrap();
        let card = attempt_card(&node, &[]);
        assert_eq!(
            card.pointer("/recorded_conclusions/missing_or_open"),
            Some(&json!([]))
        );
        assert_eq!(
            card.pointer("/recorded_conclusions/missing_or_open_recorded"),
            Some(&json!(false))
        );
    }
    #[test]
    fn attempt_card_retains_positive_and_conditional_results() {
        let node: Node = serde_json::from_value(json!({
            "id":"contract:positive-results",
            "type":"attempt",
            "label":"Scoped response",
            "attrs":{"original_contract":{
                "verdict":"PARTIAL",
                "positive_results":["The finite matrix has a one-dimensional kernel."],
                "analytic_conditional_results":["Under the declared seam condition, the index is stable."],
                "leaves_open":["The continuum limit remains unproved."]
            }}
        })).unwrap();
        let card = attempt_card(&node, &["source:result".to_owned()]);
        assert_eq!(
            card.pointer("/recorded_conclusions/findings/0"),
            Some(&json!("The finite matrix has a one-dimensional kernel."))
        );
        assert_eq!(
            card.pointer("/recorded_conclusions/positive_results/0"),
            Some(&json!("The finite matrix has a one-dimensional kernel."))
        );
        assert_eq!(
            card.pointer("/recorded_conclusions/findings")
                .and_then(Value::as_array)
                .map(Vec::len),
            Some(1)
        );
        assert_eq!(
            card.pointer("/recorded_conclusions/conditional_results/0"),
            Some(&json!(
                "Under the declared seam condition, the index is stable."
            ))
        );
    }
    #[tokio::test]
    async fn orientation_overlay_adds_english_structure_and_goal_fields_and_reloads() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-orientation-overlay-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let data = root.join("data/knowledge.json");
        fs::write(
            &data,
            serde_json::to_vec(&json!({"nodes":[
                {"id":"structure:p1-p2","type":"structure","label":"P1/P2"},
                {"id":"goal:T1","type":"goal","label":"Origin goal"}
            ],"edges":[],"sources":[]}))
            .unwrap(),
        )
        .unwrap();
        let orientation = root.join("orientation.json");
        fs::write(
            &orientation,
            serde_json::to_vec(&json!({
                "structures":[{"id":"structure:p1-p2","label_en":"Seam and carrier","summary_en":"The oriented seam and marked carrier."}],
                "goals":[{"id":"goal:T1","label_en":"Derive origin","summary_en":"Derive the selected physical origin.","required_steps_en":["Identify the source response."]}]
            })).unwrap(),
        ).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let graph = store.snapshot().await.unwrap();
        let structure = graph
            .nodes
            .iter()
            .find(|node| node.id == "structure:p1-p2")
            .unwrap();
        let goal = graph
            .nodes
            .iter()
            .find(|node| node.id == "goal:T1")
            .unwrap();
        assert_eq!(structure.attrs["label_en"], "Seam and carrier");
        assert_eq!(
            structure.attrs["summary_en"],
            "The oriented seam and marked carrier."
        );
        assert_eq!(
            goal.attrs["required_steps_en"][0],
            "Identify the source response."
        );

        fs::write(
            &orientation,
            serde_json::to_vec(&json!({
                "structures":[{"id":"structure:p1-p2","label_en":"Updated seam and carrier","summary_en":"The oriented seam and marked carrier."}],
                "goals":[{"id":"goal:T1","label_en":"Derive origin","summary_en":"Derive the selected physical origin.","required_steps_en":["Identify the source response."]}]
            })).unwrap(),
        ).unwrap();
        let reloaded = store.snapshot().await.unwrap();
        assert_eq!(
            reloaded
                .nodes
                .iter()
                .find(|node| node.id == "structure:p1-p2")
                .unwrap()
                .attrs["label_en"],
            "Updated seam and carrier"
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn graph_hydrates_only_exact_cached_embedding_matches_and_reports_coverage() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-vector-cache-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let label = "Cached claim label";
        let summary = "Exact summary used to derive vector cache key";
        let graph = json!({"meta":{"schema_version":1},"nodes":[
            {"id":"claim:cached","type":"claim","label":label,"summary":summary,"attrs":{}},
            {"id":"claim:unmatched","type":"claim","label":"Different label","summary":"Different summary","attrs":{}}
        ],"edges":[],"sources":[]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let before = store.status().await.unwrap();
        assert_eq!(before["meta"]["semantic_search"]["available"], false);
        let key = embedding_content_key(label, summary);
        let cache = json!({key:vec![0.25f32;256]});
        fs::write(
            root.join("data/embeddings-cache.json"),
            serde_json::to_vec(&cache).unwrap(),
        )
        .unwrap();
        let after = store.status().await.unwrap();
        assert_eq!(after["meta"]["semantic_search"]["available"], true);
        assert_eq!(after["meta"]["semantic_search"]["hydrated_nodes"], 1);
        assert_eq!(after["meta"]["semantic_search"]["vector_nodes"], 1);
        assert_eq!(after["meta"]["semantic_search"]["graph_nodes"], 2);
        assert_eq!(
            after["meta"]["semantic_search"]["model"],
            "text-embedding-3-small"
        );
        let snapshot = store.snapshot().await.unwrap();
        let cached_node = snapshot
            .nodes
            .iter()
            .find(|node| node.id == "claim:cached")
            .unwrap();
        assert_eq!(embedding(cached_node).unwrap(), vec![0.25f32; 256]);
        assert!(snapshot
            .nodes
            .iter()
            .find(|node| node.id == "claim:unmatched")
            .unwrap()
            .attrs
            .get("embedding")
            .is_none());
        fs::remove_dir_all(root).unwrap();
    }
    #[test]
    fn theory_map_overview_compacts_local_bilingual_patterns() {
        let map = json!({"nodes":[{"id":"map:origin","label_de":"Ursprung","state":"open","stage":"origin"}],"patterns":[{"id":"pattern:origin","label_de":"Ursprung","label_en":"Origin","summary_de":"Offene Herkunftsfrage","summary_en":"Open origin question","source_refs":["source:origin"]} ],"edges":[{"from":"map:origin","to":"map:source","kind":"requires"}]});
        let compact = compact_theory_map(&map, 20);
        assert_eq!(compact["items"][0]["label"], "Ursprung");
        assert_eq!(compact["items"][0]["summary"], "Offene Herkunftsfrage");
        assert_eq!(compact["items"][0]["source_ids"][0], "source:origin");
        assert_eq!(compact["nodes"][0]["state"], "open");
        assert_eq!(compact["edges"][0]["kind"], "requires");
        assert_eq!(compact["counts"]["nodes"], 1);
        assert_eq!(compact["counts"]["patterns"], 1);
        assert_eq!(compact["counts"]["edges"], 1);
    }
    #[test]
    fn real_theory_map_local_source_refs_are_current() {
        use sha2::{Digest, Sha256};
        let atlas = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let repo = atlas.parent().unwrap();
        let path = atlas.join("theory_map.json");
        if !path.exists() { return; }
        let map: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
        fn check(value: &Value, repo: &Path, checked: &mut usize) {
            match value {
                Value::Object(object) => {
                    if let (Some(path), Some(expected)) = (object.get("path").and_then(Value::as_str), object.get("sha256").and_then(Value::as_str)) {
                        if object.get("source_location_state").and_then(Value::as_str) == Some("moved_to_external_repository") {
                            assert!(object.get("source_url").and_then(Value::as_str).is_some_and(|url| url.starts_with("https://")));
                        } else {
                            let bytes = fs::read(repo.join(path)).unwrap_or_else(|error| panic!("{path}: {error}"));
                            assert_eq!(format!("{:x}", Sha256::digest(&bytes)), expected, "stale authored theory-map source: {path}");
                            let start = object.get("line_start").and_then(Value::as_u64).unwrap() as usize;
                            let end = object.get("line_end").and_then(Value::as_u64).unwrap() as usize;
                            assert!(start > 0 && end >= start && end <= String::from_utf8_lossy(&bytes).lines().count(), "invalid source extent: {path}");
                            *checked += 1;
                        }
                    }
                    for child in object.values() { check(child, repo, checked); }
                }
                Value::Array(array) => { for child in array { check(child, repo, checked); } }
                _ => {}
            }
        }
        let mut checked = 0;
        check(&map, repo, &mut checked);
        assert!(checked > 0);
    }

    #[test]
    fn real_theory_overview_keeps_complete_map_skeleton_without_excerpts() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("theory_map.json");
        if !path.exists() { return; }
        let map: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
        let compact = compact_theory_map(&map, 20);
        let node_count = map["nodes"].as_array().unwrap().len();
        let edge_count = map["edges"].as_array().unwrap().len();
        assert!(node_count > 0 && edge_count > 0);
        assert_eq!(compact["counts"]["nodes"], node_count);
        assert_eq!(compact["counts"]["edges"], edge_count);
        assert_eq!(compact["nodes"].as_array().unwrap().len(), node_count);
        assert_eq!(compact["edges"].as_array().unwrap().len(), edge_count);
        let pattern_count = map["patterns"].as_array().unwrap().len();
        let retained = compact["items"].as_array().unwrap().len();
        assert_eq!(compact["omitted_patterns"], pattern_count - retained);
        assert_eq!(compact["pattern_index"].as_array().unwrap().len(), pattern_count);
        for id in ["original-source", "gravity", "verification"] {
            assert!(compact["evidence_areas"].as_array().unwrap().iter().any(|area| area["id"] == id));
        }
        assert_eq!(compact["truncated"], pattern_count > retained);
        for (original, item) in map["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .zip(compact["nodes"].as_array().unwrap())
        {
            assert_eq!(original["id"], item["id"]);
            if let Some(qualification) = original.get("overview_de").and_then(Value::as_str) {
                assert_eq!(item["qualification"], clip(qualification, 180));
            }
        }
        let p1 = compact["nodes"].as_array().unwrap().iter().find(|n| n["id"] == "map:p1").unwrap();
        assert!(p1["qualification"].as_str().unwrap().contains("Randkern"));
        let encoded = serde_json::to_string(&compact).unwrap();
        assert!(
            encoded.len() < 12_000,
            "compact map is {} bytes",
            encoded.len()
        );
        for needle in ["P1", "P2", "E₈", "Universalraum", "Herkunftsherleitung"] {
            assert!(encoded.contains(needle), "missing {needle}");
        }
        assert!(compact["edges"].as_array().unwrap().iter().all(|edge| {
            edge.get("from").is_some_and(Value::is_string)
                && edge.get("to").is_some_and(Value::is_string)
        }));
        assert!(
            encoded.contains("open"),
            "open bridges should remain visible"
        );
        assert!(!encoded.contains("excerpt"));
    }
    #[test]
    fn malformed_edges_fail_validation() {
        let raw=br#"{"nodes":[{"id":"a","type":"claim"}],"edges":[{"id":"e","src":"a","dst":"missing","type":"supports"}]}"#;
        assert!(parse_graph(raw).is_err());
    }
    #[tokio::test]
    async fn graph_queries_preserve_refs_edges_and_source_ranges() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-service-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let source_text = "line one\nfixed point argument\nline three\n";
        fs::write(root.join("paper.md"), source_text).unwrap();
        let hash = format!("{:x}", Sha256::digest(source_text.as_bytes()));
        let graph = json!({
            "meta":{"schema_version":1,"generated_at":"now"},
            "nodes":[
                {"id":"claim:EM.FP.01","type":"claim","label":"Fixed point existence","summary":"Fixed point argument","attrs":{"source_refs":["source:paper"]}},
                {"id":"result:fixed-point","type":"result","label":"Fixed point","summary":"A fixed point exists","attrs":{}}
            ],
            "edges":[{"id":"edge:1","src":"claim:EM.FP.01","dst":"result:fixed-point","type":"yields","attrs":{"provenance":{"reviewed":false}}}],
            "sources":[{"id":"source:paper","path":"paper.md","sha256":hash,"line_start":2,"line_end":2,"page":3}]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let search = store
            .search(SearchRequest {
                query: "existence".into(),
                kind: None,
                limit: None,
            })
            .await
            .unwrap();
        assert_eq!(search["results"][0]["id"], "claim:EM.FP.01");
        let dossier = store.dossier("claim:EM.FP.01").await.unwrap();
        assert_eq!(dossier["related_nodes"][0]["id"], "result:fixed-point");
        assert_eq!(dossier["sources"][0]["id"], "source:paper");
        let trace = store.trace("claim:EM.FP.01", "outbound", 1).await.unwrap();
        assert_eq!(trace["edges"][0]["type"], "yields");
        let source = store.source("source:paper", None, None).await.unwrap();
        assert_eq!(source["text"], "fixed point argument");
        assert_eq!(source["page"], 3);
        assert_eq!(source["line_start"], 2);
        assert_eq!(source["line_end"], 2);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn search_collapses_contract_question_when_matching_attempt_is_also_found() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-search-contract-question-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let question = "Can the source normalization select the scalar Hamiltonian?";
        let graph = json!({
            "nodes":[
                {"id":"contract:source-normalization","type":"attempt","label":"Source normalization test","summary":"Recorded finite contract","attrs":{"name":"source-normalization","question":question,"original_contract":{"question":question},"verdict_raw":"PARTIAL"}},
                {"id":"contract_question:source-normalization","type":"contract_question","label":question,"summary":"Question copied from contract index","attrs":{"contract_id":"source-normalization","question":question}}
            ],
            "edges":[],"sources":[]
        });
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let result = store
            .search(SearchRequest {
                query: "source normalization scalar Hamiltonian".into(),
                kind: None,
                limit: Some(7),
            })
            .await
            .unwrap();
        let hits = result["results"].as_array().unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0]["id"], "contract:source-normalization");
        assert_eq!(
            hits[0]["contract_question_match"]["id"],
            "contract_question:source-normalization"
        );
        assert!(hits[0]["contract_question_match"]["why_retained"]
            .as_str()
            .unwrap()
            .contains("instead of as a duplicate"));
        assert_eq!(result["total"], 1);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn compact_dossier_keeps_useful_contract_evidence_and_counts_omissions() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-dossier-compact-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let huge = "repeated certification detail ".repeat(20_000);
        let mut nodes = vec![json!({
            "id":"attempt:charged-register","type":"attempt","label":"Charged source register","summary":"Bounded result",
            "attrs":{"source_refs":["source:register"],"verdict_raw":"PARTIAL","positive_results":["The finite register test passes."],"kills":["The uncharged model is not selected."],"leaves_open":["Derive the source normalization."],"findings":["Recorded outer finding."],"baseline_requirements":["Evaluate the native baseline before comparing shifted responses."],"original_contract":{"question":"Does the charged register select the source state?","verdict":"PARTIAL","positive_results":[],"kills":[],"leaves_open":[],"verbose_certificate":huge}}
        })];
        let mut edges = Vec::new();
        for i in 0..20 {
            let id = format!("claim:{i}");
            nodes.push(json!({"id":id,"type":"claim","label":format!("Neighbor {i}"),"summary":"neighbor context","attrs":{}}));
            edges.push(
                json!({"src":"attempt:charged-register","dst":id,"type":"documents","attrs":{}}),
            );
        }
        let graph = json!({"nodes":nodes,"edges":edges,"sources":[{"id":"source:register","path":"experiments/register.md","sha256":"abc","line_start":4,"line_end":10}]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let compact = store
            .dossier_compact("attempt:charged-register", 24_000, 3)
            .await
            .unwrap();
        assert!(json_char_count(&compact) < 24_000);
        assert_eq!(
            compact["node"]["attrs"]["original_contract"]["positive_results"][0],
            "The finite register test passes."
        );
        assert_eq!(
            compact["node"]["attrs"]["original_contract"]["kills"][0],
            "The uncharged model is not selected."
        );
        assert_eq!(
            compact["node"]["attrs"]["original_contract"]["leaves_open"][0],
            "Derive the source normalization."
        );
        assert_eq!(
            compact["node"]["attrs"]["original_contract"]["findings"][0],
            "Recorded outer finding."
        );
        assert_eq!(
            compact["node"]["attrs"]["original_contract"]["baseline_requirements"][0],
            "Evaluate the native baseline before comparing shifted responses."
        );
        assert_eq!(compact["sources"][0]["path"], "experiments/register.md");
        assert_eq!(compact["related_nodes"].as_array().unwrap().len(), 3);
        assert_eq!(compact["omissions"]["related_nodes"], 17);
        assert!(
            compact["node"]["attrs"]["omitted_counts"]["original_contract_fields"]
                .as_u64()
                .unwrap()
                > 0
        );
        let full = store.dossier("attempt:charged-register").await.unwrap();
        assert_eq!(
            full["node"]["attrs"]["original_contract"]["verbose_certificate"]
                .as_str()
                .unwrap()
                .len(),
            huge.len()
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn compact_dossier_meets_twelve_k_budget_with_premises_status_and_sources() {
        let long = "recorded premise evidence ".repeat(50);
        let source_ids: Vec<String> = (0..40).map(|index| format!("source:{index:02}")).collect();
        let sources: Vec<Value> = source_ids
            .iter()
            .map(|id| {
                json!({"id":id,"path":if id=="source:00"{"experiments/AX/PROOF.md"}else{"experiments/AX/source.md"},"sha256":"hash","line_start":1,"line_end":50,"text":long})
            })
            .collect();
        let evidence = || {
            (0..12)
                .map(|index| {
                    json!({"id":format!("evidence:{index}"),"kind":"declared_input","text":long,"statement":long,"scope":"finite algebra only","source_refs":source_ids})
                })
                .collect::<Vec<_>>()
        };
        let dossier = json!({
            "node":{"id":"claim:AX.P1.01","type":"claim","label":"AX.P1 declared identity","summary":long,"attrs":{"ledger":{"status":"Axiom"},"assessment":{"category":"recorded_established","scope":"finite arithmetic only","authority":"ledger"},"enrichment_state":"enriched","source_refs":source_ids,"semantic":{"overview_de":long,"overview_en":long,"premises":evidence(),"results":evidence(),"proof_gaps":evidence(),"input_audit":evidence()}}},
            "related_nodes":(0..20).map(|i|json!({"id":format!("result:{i}"),"type":"result","label":"related","summary":long})).collect::<Vec<_>>(),
            "edges":(0..40).map(|i|json!({"src":"claim:AX.P1.01","dst":format!("result:{i}"),"type":"has_result"})).collect::<Vec<_>>(),
            "sources":sources,
            "counts":{"related_nodes":20,"edges":40,"sources":40}
        });
        let compact = compact_dossier_projection(&dossier, 12_000, 1);
        assert!(
            json_char_count(&compact) <= 12_000,
            "chars={}",
            json_char_count(&compact)
        );
        assert_eq!(compact["node"]["attrs"]["ledger_status"], "Axiom");
        assert_eq!(
            compact["node"]["attrs"]["assessment"]["category"],
            "recorded_established"
        );
        assert!(!compact["node"]["attrs"]["semantic"]["premises"]
            .as_array()
            .unwrap()
            .is_empty());
        assert!(compact["sources"]
            .as_array()
            .unwrap()
            .iter()
            .any(|source| source["path"] == "experiments/AX/PROOF.md"));
        assert!(compact["omissions"]["source_refs"].as_u64().unwrap() > 0);
    }
    #[tokio::test]
    async fn whole_file_sources_support_bounded_windows_with_registered_extent() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-source-window-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let text = (1..=250)
            .map(|n| format!("line {n}"))
            .collect::<Vec<_>>()
            .join("\n");
        fs::write(root.join("corpus.tex"), &text).unwrap();
        let hash = format!("{:x}", Sha256::digest(text.as_bytes()));
        let graph = json!({"nodes":[],"edges":[],"sources":[
            {"id":"source:corpus","path":"corpus.tex","sha256":hash,"line_count":250},
            {"id":"source:bounded-ref","path":"corpus.tex","sha256":hash,"line_count":250,"line_start":85,"line_end":112},
            {"id":"source:excerpt-ref","path":"corpus.tex","sha256":hash,"line_count":250,"line_start":1,"line_end":80,"text":"registered excerpt"}
        ]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let default = store.source("source:corpus", None, None).await.unwrap();
        assert_eq!(default["line_start"], 1);
        assert_eq!(default["line_end"], 200);
        assert_eq!(default["text"].as_str().unwrap().lines().count(), 200);
        let selected = store
            .source("source:corpus", Some(85), Some(112))
            .await
            .unwrap();
        assert_eq!(selected["line_start"], 85);
        assert_eq!(selected["line_end"], 112);
        assert!(selected["text"].as_str().unwrap().starts_with("line 85\n"));
        assert!(selected["text"].as_str().unwrap().ends_with("line 112"));
        assert!(store
            .source("source:corpus", Some(85), Some(285))
            .await
            .is_err());
        assert!(store
            .source("source:corpus", Some(1), Some(201))
            .await
            .is_err());
        let bounded = store
            .source("source:bounded-ref", None, None)
            .await
            .unwrap();
        assert_eq!(bounded["line_start"], 85);
        assert_eq!(bounded["line_end"], 112);
        assert!(store
            .source("source:bounded-ref", Some(84), Some(100))
            .await
            .is_err());
        assert!(store
            .source("source:bounded-ref", Some(100), Some(113))
            .await
            .is_err());
        let excerpt_window = store
            .source("source:excerpt-ref", Some(1), Some(80))
            .await
            .unwrap();
        assert_eq!(excerpt_window["line_start"], 1);
        assert_eq!(excerpt_window["line_end"], 80);
        assert!(excerpt_window["text"]
            .as_str()
            .unwrap()
            .starts_with("line 1\n"));
        assert!(excerpt_window["text"]
            .as_str()
            .unwrap()
            .ends_with("line 80"));
        let excerpt_default = store
            .source("source:excerpt-ref", None, None)
            .await
            .unwrap();
        assert_eq!(excerpt_default["text"], excerpt_window["text"]);
        assert_eq!(excerpt_default["extent"]["kind"], "original_lines");
        assert_eq!(
            excerpt_default["extent"]["registered_extent_complete"],
            true
        );
        assert_eq!(default["extent"]["registered_extent_complete"], false);
        assert_eq!(default["extent"]["next_request"]["line_start"], 201);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn stored_pdf_preview_and_byte_truncation_have_explicit_actual_extent() {
        let root = std::env::temp_dir().join(format!(
            "atlas-source-extent-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let pdf = b"%PDF dummy fixture";
        let text = format!("{}\noperator relation follows here", "x".repeat(210_000));
        fs::write(root.join("source.pdf"), pdf).unwrap();
        fs::write(root.join("long.txt"), &text).unwrap();
        let graph = json!({"nodes":[],"edges":[],"sources":[
            {"id":"source:pdf","path":"source.pdf","sha256":format!("{:x}",Sha256::digest(pdf)),
                "line_start":1,"line_end":90,"text":"stored preview","page":2},
            {"id":"source:long","path":"long.txt","sha256":format!("{:x}",Sha256::digest(text.as_bytes())),
                "line_start":1,"line_end":2,"text":"unreliable short preview"}
        ]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let pdf = store.source("source:pdf", None, None).await.unwrap();
        assert_eq!(pdf["text"], "stored preview");
        assert_eq!(pdf["extent"]["kind"], "stored_pdf_excerpt");
        assert_eq!(pdf["extent"]["range_complete"], false);
        assert!(pdf["line_start"].is_null() && pdf["line_end"].is_null());
        assert!(store.source("source:pdf", Some(1), Some(2)).await.is_err());
        let long = store.source("source:long", None, None).await.unwrap();
        assert_eq!(long["extent"]["truncated"], true);
        assert_eq!(long["extent"]["range_complete"], false);
        assert_eq!(long["extent"]["last_line_complete"], false);
        assert_eq!(long["line_end"], 1);
        assert_eq!(long["extent"]["requested_line_end"], 2);
        fs::write(root.join("long.txt"), "changed").unwrap();
        assert!(store
            .source("source:long", None, None)
            .await
            .unwrap_err()
            .to_string()
            .contains("stale"));
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn live_plan_source_handle_follows_changes_without_graph_reload() {
        let root = std::env::temp_dir().join(format!(
            "atlas-plan-source-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let directory = root.join("docs/research/fundamental-plan");
        fs::create_dir_all(&directory).unwrap();
        fs::write(directory.join("progress.json"), "{\"revision\":1}").unwrap();
        let data = root.join("knowledge.json");
        fs::write(&data, b"{\"nodes\":[],\"edges\":[],\"sources\":[]}").unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let first = store
            .source("source:fundamental-plan:progress.json", None, None)
            .await
            .unwrap();
        fs::write(directory.join("progress.json"), "{\"revision\":2}").unwrap();
        let second = store
            .source("source:fundamental-plan:progress.json", None, None)
            .await
            .unwrap();
        assert_eq!(second["text"], "{\"revision\":2}");
        assert_ne!(first["sha256"], second["sha256"]);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn exact_claim_id_outranks_repeated_mentions_and_search_summary_is_bounded() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-exact-search-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let full_summary = "complete claim dossier summary ".repeat(20);
        let mut nodes = vec![
            json!({"id":"claim:EM.FP.01","type":"claim","label":"Fixed Point Claim","summary":full_summary,"attrs":{"enrichment_state":"enriched"}}),
        ];
        let repeated = "EM.FP.01 ".repeat(100);
        for i in 0..41 {
            let id = format!("result:QG.AMB.{i:02}");
            nodes.push(json!({"id":id,"type":"result","label":"Unrelated result","summary":repeated,"attrs":{}}));
        }
        let graph = json!({"nodes":nodes,"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let search = store
            .search(SearchRequest {
                query: "EM.FP.01".into(),
                kind: None,
                limit: Some(40),
            })
            .await
            .unwrap();
        assert_eq!(search["results"].as_array().unwrap().len(), 40);
        assert_eq!(search["results"][0]["id"], "claim:EM.FP.01");
        assert!(search["results"][0].get("attrs").is_none());
        assert!(
            search["results"][0]["summary"]
                .as_str()
                .unwrap()
                .chars()
                .count()
                <= 350
        );
        assert!(search["results"]
            .as_array()
            .unwrap()
            .iter()
            .any(|n| n["id"] == "claim:EM.FP.01"));
        assert_eq!(
            search["results"][0]["summary"]
                .as_str()
                .unwrap()
                .chars()
                .count(),
            350
        );
        let prefixed = store
            .search(SearchRequest {
                query: "claim:EM.FP.01".into(),
                kind: Some("claim".into()),
                limit: Some(10),
            })
            .await
            .unwrap();
        assert_eq!(prefixed["results"][0]["id"], "claim:EM.FP.01");
        let dossier = store.dossier("claim:EM.FP.01").await.unwrap();
        assert_eq!(dossier["node"]["summary"].as_str().unwrap(), full_summary);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn exact_claim_id_outranks_authored_pattern_hits() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-exact-pattern-search-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        fs::write(
            root.join("data/knowledge.json"),
            serde_json::to_vec(&json!({"nodes":[
                {"id":"claim:EM.FP.01","type":"claim","label":"Fixed point","summary":"Exact claim record","attrs":{}}
            ],"edges":[],"sources":[]})).unwrap(),
        )
        .unwrap();
        fs::write(
            root.join("theory_map.json"),
            serde_json::to_vec(&json!({"nodes":[],"edges":[],"patterns":[
                {"id":"pattern:readout-information-loss","label_de":"Readout mentions claim EM FP 01","summary_de":"Navigation pattern","source_refs":[]}
            ]})).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(root.join("data/knowledge.json"), root.clone()).unwrap();
        let result = store
            .search(SearchRequest {
                query: "claim:EM.FP.01".into(),
                kind: None,
                limit: Some(5),
            })
            .await
            .unwrap();
        assert_eq!(result["results"][0]["id"], "claim:EM.FP.01");
        assert!(result["results"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["id"] == "pattern:readout-information-loss"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn search_compacts_large_contracts_and_ranks_field_coverage_not_repetition() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-search-compact-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let repeated = "shared seed alpha flavor ".repeat(6000);
        let graph = json!({"nodes":[
            {"id":"structure:shared-seed","type":"structure","label":"Shared seed φ0/λ","summary":"Flavor and alpha readouts share a seed; this is a recorded correlation, not independent evidence.","attrs":{"keywords":["alpha","flavor","shared seed"],"label_en":"Shared seed phi0/lambda","summary_en":"Alpha and flavor readouts share one recorded seed.","translations":{"en":{"title":"Shared phi0 seed","overview":"The two readouts share a seed."}}}},
            {"id":"contract:noise","type":"attempt","label":"Unrelated contract","summary":"","attrs":{"original_contract":{"question":repeated,"verdict":"PARTIAL"}}}
        ],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let first = store
            .search(SearchRequest {
                query: "shared seed alpha flavor".into(),
                kind: None,
                limit: Some(5),
            })
            .await
            .unwrap();
        assert_eq!(first["results"][0]["id"], "structure:shared-seed");
        assert!(first["results"]
            .as_array()
            .unwrap()
            .iter()
            .all(|item| item.get("attrs").is_none()));
        assert_eq!(first["results"][0]["label_en"], "Shared seed phi0/lambda");
        assert_eq!(
            first["results"][0]["summary_en"],
            "Alpha and flavor readouts share one recorded seed."
        );
        assert_eq!(
            first["results"][0]["overview_en"],
            "The two readouts share a seed."
        );
        let english_only = store
            .search(SearchRequest {
                query: "readouts share recorded seed".into(),
                kind: Some("structure".into()),
                limit: Some(5),
            })
            .await
            .unwrap();
        assert_eq!(english_only["results"][0]["id"], "structure:shared-seed");
        assert!(serde_json::to_string(&first).unwrap().len() < 4_000);
        let second = store
            .search(SearchRequest {
                query: "shared seed alpha flavor".into(),
                kind: None,
                limit: Some(5),
            })
            .await
            .unwrap();
        assert_eq!(first, second);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn search_unprefixed_contract_name_ranks_canonical_record_first() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-search-canonical-contract-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let slug = "neutral-response-tomography-20261013";
        let graph = json!({"nodes":[
            {"id":format!("contract_project:{slug}"),"type":"contract_project","label":slug,"summary":"project container","attrs":{}},
            {"id":format!("experiment_record:{slug}"),"type":"experiment_record","label":slug,"summary":"experiment record","attrs":{}},
            {"id":format!("contract:{slug}"),"type":"contract","label":"Which native neutral-field responses repair the rank-one port?","summary":"Canonical scoped contract.","attrs":{"original_contract":{"name":slug,"verdict":"PARTIAL"}}}
        ],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let result = store
            .search(SearchRequest {
                query: slug.into(),
                kind: None,
                limit: Some(3),
            })
            .await
            .unwrap();
        assert_eq!(result["results"][0]["id"], format!("contract:{slug}"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn search_returns_authored_theory_patterns_with_source_refs_and_review_state() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-search-pattern-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        fs::create_dir_all(root.join("experiments")).unwrap();
        let source_text = "Normalized response and source map record.\n";
        fs::write(root.join("experiments/readout.md"), source_text).unwrap();
        let source_hash = format!("{:x}", Sha256::digest(source_text.as_bytes()));
        fs::write(
            root.join("data/knowledge.json"),
            serde_json::to_vec(&json!({"nodes":[
                {"id":"topic:gravity","type":"topic","label":"Gravity","summary":"Gravity topic","attrs":{}}
            ],"edges":[],"sources":[]})).unwrap(),
        )
        .unwrap();
        fs::write(
            root.join("theory_map.json"),
            serde_json::to_vec(&json!({
                "nodes":[],
                "edges":[],
                "patterns":[{
                    "id":"pattern:readout-information-loss",
                    "kind":"cross_source_inference",
                    "reviewed":false,
                    "label_de":"Arbeitsinferenz: Auslesungen verlieren Information",
                    "label_en":"Working inference: readouts lose information",
                    "summary_de":"Normierte Materieantworten bestimmen keinen rein geometrischen Zusatz.",
                    "summary_en":"Normalized matter responses do not determine a purely geometric addition.",
                    "source_refs":[{"source_id":"source:readout","path":"experiments/readout.md","sha256":source_hash,"line_start":1,"line_end":1,"excerpt":"Verified readout excerpt."}],
                    "related_ids":["topic:gravity"],
                    "conditions":[{"label_de":"Bedingt","label_en":"Conditional"}],
                    "counterevidence":[{"label_de":"Andere Antworten bleiben möglich","label_en":"Other responses remain possible"}]
                }]
            })).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(root.join("data/knowledge.json"), root.clone()).unwrap();
        let result = store
            .search(SearchRequest {
                query: "normalized matter response".into(),
                kind: None,
                limit: Some(5),
            })
            .await
            .unwrap();
        assert_eq!(
            result["results"][0]["id"],
            "pattern:readout-information-loss"
        );
        assert_eq!(result["results"][0]["type"], "theory_pattern");
        assert_eq!(result["results"][0]["kind"], "cross_source_inference");
        assert_eq!(result["results"][0]["reviewed"], false);
        assert_eq!(result["results"][0]["source_ids"][0], "source:readout");
        assert!(result["results"][0]["summary_en"]
            .as_str()
            .unwrap()
            .contains("matter responses"));
        let typed = store
            .search(SearchRequest {
                query: "normalized matter response".into(),
                kind: Some("theory_pattern".into()),
                limit: Some(5),
            })
            .await
            .unwrap();
        assert_eq!(typed["count"], 1);
        assert_eq!(
            typed["results"][0]["id"],
            "pattern:readout-information-loss"
        );
        let dossier = store
            .dossier("pattern:readout-information-loss")
            .await
            .unwrap();
        assert_eq!(dossier["node"]["type"], "theory_pattern");
        assert_eq!(dossier["node"]["attrs"]["reviewed"], false);
        assert_eq!(dossier["related_nodes"][0]["id"], "topic:gravity");
        assert_eq!(dossier["edges"][0]["type"], "about_topic");
        assert_eq!(dossier["sources"][0]["id"], "source:readout");
        let source = store.source("source:readout", None, None).await.unwrap();
        assert_eq!(source["text"], source_text.trim_end());
        assert_eq!(source["extent"]["kind"], "original_lines");
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn research_context_is_bounded_deterministic_and_groups_only_explicit_shared_identity() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-context-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({
            "meta":{"schema_version":1,"graph_integrity":{"unresolved_semantic_links":2,"issues":[{"src":"claim:a","dst":"result:missing","type":"has_result","reason":"missing_endpoint"}],"notice":"Some links unresolved."}},
            "nodes":[
                {"id":"topic:test","type":"topic","label":"Topic question","summary":"Recorded topic","attrs":{"topic_orientation":{"role_de":"role","next_bridge_de":"derive dynamics"}}},
                {"id":"claim:a","type":"claim","label":"First claim","summary":"claim scope","attrs":{"assessment":{"category":"mixed","scope":"conditional statement"},"enrichment_state":"enriched","semantic":{"topic_ids":["topic:test"],"seed_ids":["seed:common"],"objects":[{"id":"object:local-a","name":"X","identity_scope":"claim_local"},{"name":"same-symbol","identity_scope":"unspecified"}],"input_audit":[{"name":"scale","origin":"free_parameter","impact":"conditional","source_refs":["source:a"]}],"proof_gaps":[{"missing_step":"dynamics"}],"results":[{"id":"result:a","statement":"A scoped result","scope":"narrow","source_refs":["source:a"]}]}}},
                {"id":"claim:b","type":"claim","label":"Second claim","summary":"another claim","attrs":{"assessment":{"category":"open"},"enrichment_state":"pending","semantic":{"topic_ids":["topic:test"],"seed_ids":["seed:common"],"objects":[{"id":"object:local-b","name":"X","identity_scope":"claim_local"},{"name":"same-symbol","identity_scope":"unspecified"}]}}},
                {"id":"attempt:a","type":"attempt","label":"Recorded finite test","summary":"historical","attrs":{"verdict_raw":"PARTIAL","original_contract":{"rules_out":["Only within finite scope"]},"source_refs":["source:a"]}},
                {"id":"structure:a","type":"structure","label":"Shared structure","summary":"Explicit link","attrs":{}}
            ],
            "edges":[
                {"id":"t1","src":"topic:test","dst":"claim:a","type":"about","attrs":{}},
                {"id":"t2","src":"topic:test","dst":"claim:b","type":"about","attrs":{}},
                {"id":"s1","src":"claim:a","dst":"structure:a","type":"uses_structure","attrs":{}},
                {"id":"s2","src":"claim:b","dst":"structure:a","type":"uses_structure","attrs":{}},
                {"id":"a1","src":"attempt:a","dst":"claim:a","type":"tests","attrs":{}}
            ],
            "sources":[{"id":"source:a","path":"missing.txt","sha256":""}]
        });
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let context = store
            .research_context("topic:test", 12_000, 0)
            .await
            .unwrap();
        let encoded = serde_json::to_string(&context).unwrap();
        assert!(encoded.len() <= 12_000);
        assert_eq!(context["coverage"]["claims_total"], 2);
        assert_eq!(context["scope"], Value::Null);
        assert_eq!(context["scope_recorded"], false);
        assert_eq!(
            context["recorded_shared_inputs"]["seeds"][0]["seed_id"],
            "seed:common"
        );
        assert!(context["recorded_shared_inputs"]["objects"]
            .as_array()
            .unwrap()
            .is_empty());
        assert!(context["recorded_shared_inputs"]["identity_notice"]
            .as_str()
            .unwrap()
            .contains("labels, symbols"));
        assert!(context["recorded_connections"]
            .as_array()
            .unwrap()
            .iter()
            .any(|edge| edge["edge_type"] == "uses_structure"));
        assert_eq!(
            context["recorded_shared_inputs"]["structures"][0]["structure_id"],
            "structure:a"
        );
        assert_eq!(context["prior_attempts"]["items"][0]["id"], "attempt:a");
        assert_eq!(context["scoped_negatives"][0]["verdict"], "PARTIAL");
        assert_eq!(
            context["relevant_graph_integrity_issues"][0]["dst"],
            "result:missing"
        );
        assert_eq!(
            context,
            store
                .research_context("topic:test", 12_000, 0)
                .await
                .unwrap()
        );
        let small = store.research_context("topic:test", 1000, 0).await.unwrap();
        assert!(serde_json::to_string(&small).unwrap().chars().count() <= 1000);
        assert_eq!(small["id"], "topic:test");
        assert!(small.get("coverage").is_some());
        assert!(small.get("claim_status").is_some());
        assert!(store.research_context("topic:test", 0, 0).await.is_err());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn research_context_budget_fallback_keeps_status_premises_and_source_paths() {
        let long = "recorded detail ".repeat(120);
        let attempt = json!({
            "id":"attempt:large",
            "label":"Large prior attempt",
            "verdict":"PARTIAL",
            "question_or_scope":"Bounded conditional calculation",
            "source_refs":[{"id":"source:attempt","path":"experiments/attempt/PROOF.md","line_start":1,"line_end":40}],
            "recorded_conclusions":{
                "assumptions":[long.clone(),long.clone(),long.clone(),long.clone()],
                "positive_results":[long.clone(),long.clone(),long.clone(),long.clone()],
                "rules_out":[long.clone(),long.clone(),long.clone(),long.clone()],
                "missing_or_open":[long.clone(),long.clone(),long.clone(),long.clone()],
                "omitted_counts":{"assumptions":0,"positive_results":0,"rules_out":0,"missing_or_open":0}
            }
        });
        let packet = json!({
            "id":"claim:AX.P1.01",
            "kind":"claim",
            "question":"Declared finite claim",
            "claim_status":{"category":"recorded_established","scope":"finite algebra only"},
            "scope":"Finite arithmetic reduction; no physical derivation.",
            "scope_recorded":true,
            "coverage":{"claims_total":1,"claims_in_page":1,"recorded_attempt_total":3,"premise_records_total":1,"premise_records_in_packet":1,"premise_records_omitted":0,"graph_integrity_issue_total":0},
            "prior_attempts":{"items":[attempt.clone(),attempt.clone(),attempt.clone()],"total":3,"offset":0,"next_offset":3,"truncated":false},
            "premise_records":[{"claim_id":"claim:AX.P1.01","id":"premise:anchor","kind":"declared_input","text":"The anchor is an explicit premise, not a derived value.","source_ids":["source:ax"]}],
            "assumption_records":[{"claim_id":"claim:AX.P1.01","name":"c3","origin":"hardcoded","impact":"conditional","used_for":"The finite reduction."}],
            "what_we_have":[{"claim_id":"claim:AX.P1.01","category":"recorded_established","scope":"finite algebra only","summary":"The declared arithmetic reduction holds in scope."}],
            "next_discriminators":[{"claim_id":"claim:AX.P1.01","requirement":"Derive c3 from the original principles."}],
            "source_refs":[{"id":"source:ax","path":"proof/AX-P1.md","line_start":10,"line_end":18}],
            "recorded_connections":[{"src":"claim:AX.P1.01","dst":"result:one","edge_type":"has_result"}]
        });
        let compact = fit_research_context_budget(packet, 8_000);
        assert!(json_char_count(&compact) <= 8_000);
        assert_eq!(compact["claim_status"]["category"], "recorded_established");
        assert_eq!(compact["scope_recorded"], true);
        assert_eq!(
            compact["premise_records"][0]["text"],
            "The anchor is an explicit premise, not a derived value."
        );
        assert_eq!(compact["source_refs"][0]["path"], "proof/AX-P1.md");
        assert_eq!(
            compact["what_we_have"][0]["summary"],
            "The declared arithmetic reduction holds in scope."
        );
        assert!(compact["truncated"].as_bool().unwrap());
    }

    #[test]
    fn attempt_context_budget_keeps_root_results_and_counts_every_cut() {
        let positive_results: Vec<Value> = (0..6)
            .map(|i| {
                json!(format!(
                    "Recorded positive result {i}: the scoped calculation yields value {i}."
                ))
            })
            .collect();
        let assumptions: Vec<Value> = (0..4)
            .map(|i| json!(format!("Explicit baseline assumption {i}.")))
            .collect();
        let oversized_noise = "historical negative record ".repeat(10_000);
        let packet = json!({
            "id":"contract:native-source-reconstruction-20261014",
            "kind":"contract",
            "question":"How far does the native source reconstruction establish the scoped response?",
            "claim_status":null,
            "scope":"Only the registered finite source reconstruction; no broader physical derivation.",
            "scope_recorded":true,
            "coverage":{"claims_total":0,"recorded_attempt_total":1,"premise_records_total":0,"graph_integrity_issue_total":0},
            "attempt_record":{
                "id":"contract:native-source-reconstruction-20261014",
                "label":"Native source reconstruction",
                "verdict":"PARTIAL",
                "question_or_scope":"Registered finite source reconstruction.",
                "source_refs":[{"id":"source:proof","path":"experiments/native/PROOF.md","line_start":1,"line_end":90}],
                "recorded_conclusions":{
                    "positive_results":positive_results,
                    "resolved_results":[],
                    "conditional_results":[],
                    "assumptions":assumptions,
                    "rules_out":[],
                    "missing_or_open":[],
                    "not_claimed":[],
                    "omitted_counts":{"positive_results":0,"assumptions":0},
                    "missing_or_open_recorded":false
                }
            },
            "what_we_have":positive_results.iter().map(|statement|json!({"statement":statement,"scope":"registered finite source reconstruction"})).collect::<Vec<_>>(),
            "prior_attempts":{"items":[{"id":"attempt:old","verdict":"PARTIAL"}],"total":1,"offset":0},
            "scoped_negatives":[{"statement":oversized_noise}],
            "recorded_connections":[{"src":"attempt:old","dst":"claim:old","edge_type":"tests"}],
            "source_refs":[{"id":"source:proof","path":"experiments/native/PROOF.md","line_start":1,"line_end":90}],
            "source_ids":["source:proof"]
        });

        for budget in [6_500, 12_000] {
            let compact = fit_research_context_budget(packet.clone(), budget);
            assert!(json_char_count(&compact) <= budget);
            let card = &compact["attempt_record"];
            assert_eq!(card["id"], "contract:native-source-reconstruction-20261014");
            assert_eq!(card["verdict"], "PARTIAL");
            assert_eq!(compact["scope_recorded"], true);
            let conclusions = &card["recorded_conclusions"];
            assert_eq!(conclusions["positive_results"].as_array().unwrap().len(), 6);
            assert_eq!(conclusions["assumptions"].as_array().unwrap().len(), 4);
            assert_eq!(conclusions["omitted_counts"]["positive_results"], 0);
            assert_eq!(conclusions["omitted_counts"]["assumptions"], 0);
            assert!(compact.get("scoped_negatives").is_none());
            assert!(card["source_refs"][0]["path"]
                .as_str()
                .unwrap()
                .ends_with("PROOF.md"));
        }

        let mut expanded = packet.clone();
        expanded["attempt_record"]["recorded_conclusions"]["positive_results"] =
            json!((0..8).map(|i| format!("Positive {i}")).collect::<Vec<_>>());
        expanded["attempt_record"]["recorded_conclusions"]["assumptions"] = json!((0..8)
            .map(|i| format!("Assumption {i}"))
            .collect::<Vec<_>>());
        expanded["attempt_record"]["recorded_conclusions"]["omitted_counts"]["positive_results"] =
            json!(2);
        expanded["attempt_record"]["recorded_conclusions"]["omitted_counts"]["assumptions"] =
            json!(1);
        let compact = fit_research_context_budget(expanded, 6_500);
        let conclusions = &compact["attempt_record"]["recorded_conclusions"];
        assert_eq!(conclusions["positive_results"].as_array().unwrap().len(), 6);
        assert_eq!(conclusions["assumptions"].as_array().unwrap().len(), 6);
        assert_eq!(conclusions["omitted_counts"]["positive_results"], 4);
        assert_eq!(conclusions["omitted_counts"]["assumptions"], 3);

        let brief = fit_char_budget(packet, 6_500);
        assert!(json_char_count(&brief) <= 6_500);
        assert_eq!(brief["attempt_record"]["verdict"], "PARTIAL");
        assert_eq!(
            brief["attempt_record"]["recorded_conclusions"]["positive_results"]
                .as_array()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn claim_brief_budget_fallback_keeps_root_overview_and_does_not_claim_attempts_missing() {
        let overview = "Declared finite arithmetic identity; this is an input to the compiler and does not derive a physical source. ".repeat(14);
        let oversized = "historical record ".repeat(250);
        let brief = json!({
            "id":"claim:AX.P1.01",
            "kind":"claim",
            "question":"Declared finite claim",
            "claim_status":{"category":"declared_axiom","scope":"finite arithmetic only"},
            "scope":"Finite arithmetic reduction; no physical derivation.",
            "enrichment_state":"error",
            "what_we_have":[{"claim_id":"claim:AX.P1.01","category":"declared_axiom","scope":"finite arithmetic only","statement":overview}],
            "counts":{"claims":1,"attempts":19,"recorded_proof_gaps":null,"claims_without_enriched_gap_audit":1},
            "proof_gap_coverage":{"recorded_gaps_from_enriched_claims":null,"claims_without_enriched_gap_audit":1,"complete":false,"scope":"claim_annotations_in_scope"},
            "scoped_results":[{"statement":oversized}],
            "typed_related_records":[{"positive_results":[oversized.clone()]},{"positive_results":[oversized.clone()]}],
            "attempts":{"items":[],"total":19,"offset":0,"limit":4,"next_offset":0,"truncated":true},
            "source_ids":["source:ax"],
            "next_testable_requirements":[{"requirement":"Derive this input from the original principles."}]
        });
        let compact = fit_char_budget(brief, 8_000);
        assert!(json_char_count(&compact) <= 8_000);
        assert_eq!(compact["claim_status"]["category"], "declared_axiom");
        assert_eq!(compact["proof_gap_coverage"]["complete"], false);
        assert!(compact["what_we_have"][0]["statement"]
            .as_str()
            .unwrap()
            .contains("does not derive a physical source"));
        assert!(compact.get("attempts_page").is_none());
    }

    #[tokio::test]
    async fn research_context_keeps_topic_linked_attempts_when_claim_is_related() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-context-topic-attempt-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({
            "nodes":[
                {"id":"topic:universalraum","type":"topic","label":"Universalraum","summary":"Question about a shared universal space"},
                {"id":"claim:primary","type":"claim","label":"Primary claim","summary":"Primary","attrs":{"topic_roles":{"topic:universalraum":"primary"},"assessment":{"category":"open","scope":"only this claim"}}},
                {"id":"claim:related","type":"claim","label":"Related claim","summary":"Related","attrs":{"topic_roles":{"topic:universalraum":"related"},"assessment":{"category":"open","scope":"related only"}}},
                {"id":"claim:related-universalraum","type":"claim","label":"Universalraum shared operator bridge","summary":"A recorded bridge question","attrs":{"topic_roles":{"topic:universalraum":"related"},"assessment":{"category":"mixed","scope":"bridge only"}}},
                {"id":"attempt:universalraum","type":"attempt","label":"Universalraum contract","summary":"Historical bounded attempt","attrs":{"verdict_raw":"PARTIAL","original_contract":{"leaves_open":["derive the map"]}}}
            ],
            "edges":[
                {"src":"topic:universalraum","dst":"claim:primary","type":"about","attrs":{}},
                {"src":"topic:universalraum","dst":"claim:related","type":"about","attrs":{}},
                {"src":"topic:universalraum","dst":"claim:related-universalraum","type":"about","attrs":{}},
                {"src":"attempt:universalraum","dst":"topic:universalraum","type":"about","attrs":{"method":"linked_claim_topics"}},
                {"src":"attempt:universalraum","dst":"claim:related","type":"tests","attrs":{}}
            ],
            "sources":[]
        });
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let context = store
            .research_context("topic:universalraum", 12_000, 0)
            .await
            .unwrap();

        assert_eq!(context["coverage"]["claims_total"], 1);
        assert_eq!(context["coverage"]["related_claims_total"], 2);
        assert_eq!(context["related_claims"]["total"], 2);
        assert_eq!(
            context["related_claims"]["items"][0]["id"],
            "claim:related-universalraum"
        );
        assert_eq!(
            context["related_claims"]["items"][0]["topic_role"],
            "related"
        );
        assert_eq!(context["coverage"]["recorded_attempt_total"], 1);
        assert_eq!(
            context["prior_attempts"]["items"][0]["id"],
            "attempt:universalraum"
        );
        assert_eq!(
            context["prior_attempts"]["items"][0]["topic_link_basis"],
            "inherited_from_claim_topic_membership"
        );
        assert_eq!(context["prior_attempts"]["items"][0]["verdict"], "PARTIAL");
        assert!(context["prior_attempts"]["notice"]
            .as_str()
            .unwrap()
            .contains("inherited"));
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn research_context_attempt_pagination_uses_full_stable_topic_membership() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-context-pages-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let mut nodes =
            vec![json!({"id":"topic:test","type":"topic","label":"Topic","summary":"Question"})];
        let mut edges = Vec::new();
        for i in 0..45 {
            let id = format!("claim:{i:02}");
            nodes.push(json!({"id":id,"type":"claim","label":format!("Claim {i}"),"summary":"Member","attrs":{"semantic":{"topic_ids":["topic:test"]}}}));
            edges.push(json!({"src":id,"dst":"topic:test","type":"about","attrs":{}}));
        }
        for (group, start, end) in [("A", 0, 5), ("B", 40, 45)] {
            for i in 0..3 {
                let id = format!("attempt:{group}{i}");
                let claims: Vec<String> = (start..end).map(|n| format!("claim:{n:02}")).collect();
                nodes.push(json!({"id":id,"type":"attempt","label":format!("Run {group}{i}"),"attrs":{"claims_tested":claims,"verdict":"PARTIAL"}}));
                for claim in start..end {
                    edges.push(json!({"src":id,"dst":format!("claim:{claim:02}"),"type":"tests","attrs":{}}));
                }
            }
        }
        let data = root.join("knowledge.json");
        fs::write(
            &data,
            serde_json::to_vec(&json!({"nodes":nodes,"edges":edges,"sources":[]})).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let first = store
            .research_context("topic:test", 12_000, 0)
            .await
            .unwrap();
        let second = store
            .research_context("topic:test", 12_000, 3)
            .await
            .unwrap();
        let first_ids: Vec<&str> = first["prior_attempts"]["items"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x["id"].as_str().unwrap())
            .collect();
        let second_ids: Vec<&str> = second["prior_attempts"]["items"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x["id"].as_str().unwrap())
            .collect();
        assert_eq!(first["prior_attempts"]["total"], 6);
        assert_eq!(second["prior_attempts"]["total"], 6);
        assert!(
            first_ids.iter().all(|id| !second_ids.contains(id)),
            "pages must not repeat attempts: {first_ids:?} then {second_ids:?}"
        );
        assert_eq!(first_ids, vec!["attempt:A0", "attempt:A1", "attempt:A2"]);
        assert_eq!(second_ids, vec!["attempt:B0", "attempt:B1", "attempt:B2"]);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn attempt_research_context_and_brief_are_exact_and_keep_original_evidence() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-attempt-research-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let source_text = "The declared baseline must be evaluated before the comparison.\n";
        fs::write(root.join("baseline.md"), source_text).unwrap();
        let source_hash = format!("{:x}", Sha256::digest(source_text.as_bytes()));
        let graph = json!({
            "nodes":[
                {"id":"topic:shared","type":"topic","label":"Broad topic","summary":"A much wider question"},
                {"id":"claim:TARGET","type":"claim","label":"Target claim","summary":"Target summary","attrs":{"semantic":{"topic_ids":["topic:shared"],"results":[{"id":"r1","statement":"Bounded result","source_refs":["source:baseline"]}]}}},
                {"id":"claim:UNRELATED","type":"claim","label":"Unrelated topic member","summary":"Must not fan out","attrs":{"semantic":{"topic_ids":["topic:shared"]}}},
                {"id":"attempt:native-sector-calibration-20261013","type":"attempt","label":"Native sector calibration","summary":"Calibration record","attrs":{"verdict_raw":"PARTIAL","question":"Does calibration select the sector?","claims_tested":["TARGET"],"source_refs":["source:baseline"],"original_contract":{"positive_results":["A bounded calibration was recorded."],"analytic_conditional_results":["Only under the stated baseline."],"assumptions":["The reference baseline is fixed first."],"kills":["The uncalibrated alternative."],"leaves_open":["Whether the baseline follows from the original principles."],"original_source_scope_correction":"An algebraic scalar shift is not thereby admissible in the original profile-free source functional.","claims_tested":["TARGET"],"source_refs":["source:baseline"],"relations":[{"type":"refines","target":"source-observation-identifiability-test"}]}}},
                {"id":"contract:source-observation-identifiability-test","type":"attempt","label":"Prior source identifiability record","attrs":{"verdict_raw":"PARTIAL","scope":"A registered source witness only.","source_refs":["source:baseline"],"original_contract":{"findings":[{"kind":"supports","text":"A source witness is registered."},{"kind":"open","text":"The source functional remains to be derived."},{"kind":"refutes","text":"The specific scalar identification is not established."}]}}},
                {"id":"attempt:other","type":"attempt","label":"Unrelated attempt","attrs":{"claims_tested":["UNRELATED"]}},
                {"id":"attempt:standalone","type":"attempt","label":"Standalone historical record","attrs":{"verdict_raw":"PARTIAL","positive_results":["Recorded finding without a claim membership."]}},
                {"id":"structure:baseline","type":"structure","label":"Registered baseline","summary":"Explicit baseline"}
            ],
            "edges":[
                {"src":"attempt:native-sector-calibration-20261013","dst":"claim:TARGET","type":"tests","attrs":{}},
                {"src":"attempt:native-sector-calibration-20261013","dst":"structure:baseline","type":"uses_structure","attrs":{}},
                {"src":"claim:TARGET","dst":"topic:shared","type":"about","attrs":{}},
                {"src":"claim:UNRELATED","dst":"topic:shared","type":"about","attrs":{}},
                {"src":"attempt:other","dst":"claim:UNRELATED","type":"tests","attrs":{}}
            ],
            "sources":[{"id":"source:baseline","path":"baseline.md","sha256":source_hash,"line_start":1,"line_end":1,"text":"The declared baseline must be evaluated before the comparison."}]
        });
        fs::write(
            root.join("data/knowledge.json"),
            serde_json::to_vec(&graph).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(root.join("data/knowledge.json"), root.clone()).unwrap();
        let id = "attempt:native-sector-calibration-20261013";
        let context = store.research_context(id, 12_000, 0).await.unwrap();
        assert_eq!(context["coverage"]["claims_total"], 1);
        assert_eq!(context["coverage"]["recorded_attempt_total"], 1);
        assert_eq!(context["prior_attempts"]["items"][0]["id"], id);
        assert_eq!(
            context["prior_attempts"]["items"][0]["topic_link_basis"],
            "root_attempt_id"
        );
        assert_eq!(
            context["prior_attempts"]["items"][0]["recorded_conclusions"]["findings"],
            Value::Null
        );
        assert_eq!(context["prior_attempts"]["items"][0]["original_source_scope_correction"], "An algebraic scalar shift is not thereby admissible in the original profile-free source functional.");
        assert_eq!(
            context["prior_attempts"]["items"][0]["recorded_conclusions"]["positive_results"][0],
            "A bounded calibration was recorded."
        );
        assert_eq!(
            context["prior_attempts"]["items"][0]["recorded_conclusions"]["missing_or_open"][0],
            "Whether the baseline follows from the original principles."
        );
        assert_eq!(context["source_refs"][0]["path"], "baseline.md");
        assert_eq!(context["typed_related_records"][0]["relation"], "refines");
        assert_eq!(
            context["typed_related_records"][0]["id"],
            "contract:source-observation-identifiability-test"
        );
        assert_eq!(
            context["typed_related_records"][0]["positive_results"][0],
            "A source witness is registered."
        );
        assert_eq!(
            context["typed_related_records"][0]["leaves_open"][0],
            "The source functional remains to be derived."
        );
        assert_eq!(
            context["typed_related_records"][0]["scoped_exclusions"][0],
            "The specific scalar identification is not established."
        );
        assert!(!context["recorded_connections"]
            .as_array()
            .unwrap()
            .iter()
            .any(|edge| edge["dst"] == "topic:shared" || edge["src"] == "topic:shared"));
        let brief = store.research_brief(id, 12_000, 0).await.unwrap();
        assert_eq!(brief["counts"]["claims"], 1);
        assert_eq!(brief["counts"]["attempts"], 1);
        assert_eq!(brief["proof_gap_coverage"]["complete"], Value::Null);
        assert_eq!(
            brief["proof_gap_coverage"]["scope"],
            "not_a_full_claim_audit"
        );
        assert_eq!(
            brief["attempt_record"]["recorded_conclusions"]["assumptions"][0],
            "The reference baseline is fixed first."
        );
        assert_eq!(brief["attempt_record"]["original_source_scope_correction"], "An algebraic scalar shift is not thereby admissible in the original profile-free source functional.");
        assert_eq!(
            brief["attempt_record"]["source_refs"][0]["path"],
            "baseline.md"
        );
        assert_eq!(
            brief["typed_related_records"][0]["id"],
            "contract:source-observation-identifiability-test"
        );
        let standalone = store
            .research_brief("attempt:standalone", 12_000, 0)
            .await
            .unwrap();
        assert_eq!(standalone["proof_gap_coverage"]["complete"], Value::Null);
        assert_eq!(
            standalone["proof_gap_coverage"]["scope"],
            "not_applicable_no_claim_membership"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn contract_context_surfaces_resolved_results_before_large_exclusion_lists() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-resolved-contracts-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(SystemTime::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let proof = "Registered contract evidence.\n";
        fs::write(root.join("PROOF.md"), proof).unwrap();
        let proof_hash = format!("{:x}", Sha256::digest(proof.as_bytes()));
        let resolved_cp = "Within the stated marked orbit-scalar one-step class, complete positivity uniquely selects the recorded sector weights.";
        let resolved_holonomy = "B_Y is the second holonomy derivative along y of the same function whose zeroth order is 12T+R.";
        let kills: Vec<String> = (0..106)
            .map(|i| format!("Scoped exclusion {i:03}"))
            .collect();
        let mut not_claimed =
            vec![json!({"kind":"open","text":"Not claimed is not an open-item classification."})];
        not_claimed.extend((0..12).map(|i| json!(format!("Not claimed boundary {i:02}"))));
        let graph = json!({
            "nodes":[
                {"id":"contract:marked-seam-cp-lift-20260921","type":"attempt","label":"Marked seam CP lift","attrs":{"verdict_raw":"PARTIAL","source_refs":["source:marked"],"original_contract":{"scope":"Only within the stated marked orbit-scalar one-step class.","resolved":[resolved_cp,"The original words provide a positive operator-system readout.","The selected channel is not a nontrivial product-preserving charged-field time.","No norm-continuous embedding exists in this class."],"findings":["{\"kind\":\"supports\",\"text\":\"An explicitly conditional finite CP connection is recorded.\"}"],"kills":kills,"source_refs":["source:marked"]}}},
                {"id":"contract:seam-holonomy-generating-20261001","type":"attempt","label":"Seam holonomy generating function","attrs":{"verdict_raw":"PARTIAL","question":"Does the holonomy response determine the selection?","source_refs":["source:holonomy"],"original_contract":{"question":"Does the holonomy response determine the selection?","resolved":[resolved_holonomy,"S=48 is selected conditional on port consistency.","The family deck action lifts linearly only as mu8, with g^4=(-1)^F."],"not_claimed":not_claimed,"source_refs":["source:holonomy"]}}}
            ],
            "edges":[],
            "sources":[
                {"id":"source:marked","path":"PROOF.md","sha256":proof_hash,"line_count":1},
                {"id":"source:holonomy","path":"PROOF.md","sha256":proof_hash,"line_count":1}
            ]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();

        for (id, expected, scope_recorded) in [
            ("contract:marked-seam-cp-lift-20260921", resolved_cp, true),
            (
                "contract:seam-holonomy-generating-20261001",
                resolved_holonomy,
                false,
            ),
        ] {
            let context = store.research_context(id, 12_000, 0).await.unwrap();
            assert!(!context["what_we_have"].as_array().unwrap().is_empty());
            assert!(context["what_we_have"]
                .as_array()
                .unwrap()
                .iter()
                .any(|item| item["statement"].as_str().unwrap().contains(expected)));
            assert_eq!(context["attempt_record"]["verdict"], "PARTIAL");
            assert_eq!(context["scope_recorded"], scope_recorded);
            if !scope_recorded {
                assert!(context["scope"].is_null());
            }
            assert!(context["source_refs"]
                .as_array()
                .unwrap()
                .iter()
                .any(|source| source["path"] == "PROOF.md"));
            if id == "contract:seam-holonomy-generating-20261001" {
                let conclusions = &context["attempt_record"]["recorded_conclusions"];
                assert_eq!(conclusions["not_claimed"].as_array().unwrap().len(), 8);
                assert_eq!(conclusions["omitted_counts"]["not_claimed"], 5);
                assert!(conclusions["not_claimed"][1]
                    .as_str()
                    .unwrap()
                    .contains("Not claimed boundary 00"));
                assert!(conclusions["not_claimed"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|item| item
                        .as_str()
                        .unwrap()
                        .contains("Not claimed is not an open-item classification")));
                assert!(!conclusions["missing_or_open"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|item| item
                        .as_str()
                        .unwrap()
                        .contains("Not claimed is not an open-item classification")));
                assert!(!conclusions["rules_out"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|item| item
                        .as_str()
                        .unwrap()
                        .contains("Not claimed is not an open-item classification")));
                let broader = store.research_context(id, 18_000, 0).await.unwrap();
                assert!(json_char_count(&broader) <= 18_000);
                assert!(broader["what_we_have"].as_array().unwrap().len() >= 3);
                assert_eq!(broader["attempt_record"]["verdict"], "PARTIAL");
            }
        }

        let cp = store
            .research_brief("contract:marked-seam-cp-lift-20260921", 12_000, 0)
            .await
            .unwrap();
        let conclusions = &cp["attempt_record"]["recorded_conclusions"];
        assert_eq!(conclusions["resolved_results"][0], resolved_cp);
        assert!(conclusions["positive_results"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item == "An explicitly conditional finite CP connection is recorded."));
        assert_eq!(conclusions["rules_out"].as_array().unwrap().len(), 8);
        assert!(conclusions["omitted_counts"]["rules_out"].as_u64().unwrap() >= 98);
        assert!(conclusions
            .get("findings")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .all(|item| !item.as_str().unwrap_or("").starts_with('{')));
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn research_brief_and_context_support_exact_goal_scope_without_topic_noise() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-goal-research-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({
            "nodes":[
                {"id":"goal:universalraum","type":"goal","label":"Universalraum","summary":"Unify the specified structures","attrs":{"status":"open","claim_ids":["claim:linked"],"required_steps":["derive the shared operator map","identify its source action"]}},
                {"id":"topic:universalraum","type":"topic","label":"Universalraum","summary":"Broad question"},
                {"id":"claim:linked","type":"claim","label":"Goal-linked scoped result","summary":"A relevant bounded result","attrs":{"assessment":{"category":"conditional","scope":"only the specified representation"},"semantic":{"topic_ids":["topic:universalraum"]}}},
                {"id":"claim:edge-linked","type":"claim","label":"Explicitly addressed claim","summary":"A second exact goal member","attrs":{"assessment":{"category":"open","scope":"open derivation"}}},
                {"id":"claim:topic-noise","type":"claim","label":"Unrelated Universalraum question","summary":"Should not be pulled through topic membership","attrs":{"assessment":{"category":"open"},"semantic":{"topic_ids":["topic:universalraum"]}}},
                {"id":"attempt:goal-open","type":"attempt","label":"Goal-specific attempt","attrs":{"verdict_raw":"PARTIAL","original_contract":{"leaves_open":["derive the shared operator map"]}}},
                {"id":"attempt:tested","type":"attempt","label":"Test of exact goal claim","attrs":{"verdict_raw":"PASS"}},
                {"id":"attempt:noise","type":"attempt","label":"Unrelated topic attempt","attrs":{"verdict_raw":"FAIL"}}
            ],
            "edges":[
                {"src":"goal:universalraum","dst":"claim:edge-linked","type":"addresses_open_goal","attrs":{}},
                {"src":"topic:universalraum","dst":"claim:topic-noise","type":"about","attrs":{}},
                {"src":"goal:universalraum","dst":"attempt:goal-open","type":"leaves_open","attrs":{}},
                {"src":"attempt:tested","dst":"claim:linked","type":"tests","attrs":{}},
                {"src":"attempt:noise","dst":"claim:topic-noise","type":"tests","attrs":{}}
            ],
            "sources":[]
        });
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();

        let brief = store
            .research_brief("goal:universalraum", 12_000, 0)
            .await
            .unwrap();
        assert_eq!(brief["kind"], "goal");
        assert_eq!(brief["goal_status"], "open");
        assert_eq!(brief["counts"]["claims"], 2);
        assert_eq!(brief["counts"]["attempts"], 2);
        assert!(brief["next_testable_requirements"]
            .as_array()
            .unwrap()
            .iter()
            .any(|step| {
                step["requirement"] == "derive the shared operator map"
                    && step["origin"] == "explicit goal required step"
            }));

        let context = store
            .research_context("goal:universalraum", 12_000, 0)
            .await
            .unwrap();
        assert_eq!(context["goal_status"], "open");
        assert_eq!(context["coverage"]["claims_total"], 2);
        assert_eq!(context["coverage"]["recorded_attempt_total"], 2);
        let attempt_ids: Vec<&str> = context["prior_attempts"]["items"]
            .as_array()
            .unwrap()
            .iter()
            .map(|item| item["id"].as_str().unwrap())
            .collect();
        assert!(attempt_ids.contains(&"attempt:goal-open"));
        assert!(attempt_ids.contains(&"attempt:tested"));
        assert!(!attempt_ids.contains(&"attempt:noise"));
        assert!(context["next_discriminators"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["requirement"] == "derive the shared operator map"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn compact_attempt_brief_keeps_primary_conclusions_within_ten_k_budget() {
        let long = "A source-grounded finding with detailed scope and conditions. ".repeat(35);
        let related: Vec<Value> = (0..5)
            .map(|i| {
                json!({
                    "id":format!("contract:related-{i}"),
                    "relation":"refines",
                    "label":format!("Related contract {i}"),
                    "verdict":"PARTIAL",
                    "scope":long,
                    "positive_results":[long,long],
                    "conditional_results":[],
                    "leaves_open":[long,long],
                    "scoped_exclusions":[],
                    "unclassified_findings":[],
                    "original_source_scope_correction":long,
                    "source_refs":[{"id":"source:one","path":"proof/one.md"},{"id":"source:two","path":"proof/two.md"}]
                })
            })
            .collect();
        let mut brief = json!({
            "id":"contract:native-sector-calibration-test",
            "kind":"contract",
            "question":"Native sector calibration",
            "verdict":"PARTIAL",
            "attempt_record":{"recorded_conclusions":{"assumptions":["Baseline fixed.","D^k preparation fixed."],"positive_results":["Relative calibration only."],"missing_or_open":["Absolute normalization remains open."]},"original_source_scope_correction":"Relative invisibility does not prove admissibility in the source functional.","source_refs":[{"id":"source:own","path":"native/PROOF.txt"}]},
            "typed_related_records":related,
            "recorded_connections":[],
            "linked_claims":{"items":[]},
            "source_refs":[],
            "source_ids":[],
            "counts":{"typed_related_records":5}
        });
        compact_attempt_research_brief(&mut brief);
        let compact = fit_char_budget(brief, 10_000);
        assert!(json_char_count(&compact) <= 10_000);
        assert_eq!(compact["verdict"], "PARTIAL");
        assert_eq!(
            compact["attempt_record"]["recorded_conclusions"]["positive_results"][0],
            "Relative calibration only."
        );
        assert_eq!(
            compact["typed_related_records"].as_array().unwrap().len(),
            5
        );
        assert_eq!(
            compact["typed_related_records"][0]["positive_results"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            compact["typed_related_records"][0]["omitted_counts"]["positive_results"],
            1
        );
        assert!(compact["typed_related_records"][0]["source_refs"][0]["path"] == "proof/one.md");
    }

    #[tokio::test]
    async fn research_brief_is_deterministic_bounded_and_keeps_evidence_buckets_distinct() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-brief-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let mut nodes = vec![
            json!({"id":"topic:test","type":"topic","label":"Model origin","summary":"Where does the model come from?"}),
        ];
        let mut edges = Vec::new();
        for i in 0..500 {
            let id = format!("claim:{i:03}");
            let category = if i % 2 == 0 { "mixed" } else { "open" };
            nodes.push(json!({
                "id":id,"type":"claim","label":format!("Claim {i}"),"summary":"Short recorded summary",
                "attrs":{
                    "assessment":{"category":category,"scope":"Only the explicitly stated narrow sector"},
                    "enrichment_state":"enriched",
                    "semantic":{
                        "overview_de":"A bounded claim overview.",
                        "results":[{"id":format!("result:{i}"),"statement":"Recorded statement.","scope":"Narrow scope","source_refs":["source:paper"]}],
                        "proof_gaps":[{"missing_step":"derive the missing map","why_checks_insufficient":"finite check only","source_refs":["source:paper"]}],
                        "input_audit":[{"name":"k","origin":"free_parameter","value":"chosen","used_for":"sets normalization","impact":"conditional","source_refs":["source:paper"]}],
                        "topic_ids":["topic:test"]
                    }
                }
            }));
            edges.push(json!({"id":format!("about:{i}"),"src":"topic:test","dst":id,"type":"about","attrs":{}}));
        }
        nodes.push(json!({"id":"script:test","type":"script","label":"Static analysis","attrs":{"static_inputs":{"syntactic_input_candidates":[{"name":"scale","line":22}]}}}));
        edges.push(json!({"id":"script-link","src":"claim:000","dst":"script:test","type":"checks","attrs":{}}));
        for i in 0..15 {
            let attempt_id = format!("attempt:{i:03}");
            nodes.push(json!({"id":attempt_id,"type":"attempt","label":format!("Finite check {i}"),"attrs":{"verdict_raw":"PASS","question":format!("Finite scope {}", "long detail ".repeat(18)),"source_refs":["source:paper"]}}));
            edges.push(json!({"id":format!("attempt-link:{i}"),"src":format!("attempt:{i:03}"),"dst":format!("claim:{i:03}"),"type":"tests","attrs":{}}));
        }
        let graph = json!({"meta":{"schema_version":1},"nodes":nodes,"edges":edges,"sources":[{"id":"source:paper","path":"paper.tex","sha256":"","line_start":1,"line_end":2}]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let brief = store.research_brief("topic:test", 12_000, 0).await.unwrap();
        let encoded = serde_json::to_string(&brief).unwrap();
        assert!(
            encoded.chars().count() <= 12_000,
            "{} chars",
            encoded.chars().count()
        );
        assert_eq!(brief["counts"]["claims"], 500);
        assert_eq!(brief["counts"]["claims_by_assessment"]["mixed"], 250);
        assert_eq!(brief["counts"]["claims_by_assessment"]["open"], 250);
        assert_eq!(brief["attempts"]["total"], 15);
        let first_page = brief["attempts"]["items"].as_array().unwrap();
        assert!(brief["scoped_results"].as_array().unwrap().len() >= 3);
        assert!(
            first_page.len() >= 3,
            "budget compaction should preserve several attempt records"
        );
        assert!(first_page.iter().any(|item| item["verdict"] == "PASS"));
        assert_eq!(brief["attempts"]["next_offset"], first_page.len());
        assert_eq!(brief["attempts"]["truncated"], first_page.len() < 15);
        assert_eq!(brief["counts"]["syntactic_input_candidates"], 1);
        assert_eq!(brief["what_we_seek"]["recorded_proof_gaps"], 500);
        assert!(brief["next_testable_requirements"]
            .as_array()
            .unwrap()
            .iter()
            .all(|item| item["origin"] != "syntactic_input_candidate"));
        assert!(brief["hidden_input_risks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["origin"] == "free_parameter"));
        assert!(!encoded.contains("input_audit"));
        assert!(!encoded.contains("embedding"));
        assert_eq!(
            store.research_brief("topic:test", 12_000, 0).await.unwrap(),
            brief
        );
        let next_offset = brief["attempts"]["next_offset"].as_u64().unwrap() as usize;
        let page = store
            .research_brief("topic:test", 12_000, next_offset)
            .await
            .unwrap();
        assert_eq!(page["attempts"]["offset"], next_offset);
        if let (Some(first), Some(next)) = (
            first_page.last(),
            page["attempts"]["items"]
                .as_array()
                .and_then(|items| items.first()),
        ) {
            assert_ne!(first["id"], next["id"]);
        }
        assert!(store.research_brief("topic:test", 30_001, 0).await.is_err());
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn stale_annotation_cannot_report_complete_gap_coverage() {
        let root = std::env::temp_dir().join(format!(
            "semantic-stale-brief-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({"nodes":[{"id":"claim:test","type":"claim","label":"Test","attrs":{"enrichment_state":"enriched","cache_freshness":{"state":"source_stale","sources_current":false},"semantic":{"overview_de":"Known bounded result","proof_gaps":[{"label":"Old missing source","missing_step":"Source absent from old packet"}]}}}],"edges":[],"sources":[]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let brief = store.research_brief("claim:test", 12_000, 0).await.unwrap();
        assert_eq!(brief["proof_gap_coverage"]["complete"], false);
        assert_eq!(
            brief["proof_gap_coverage"]["claims_without_enriched_gap_audit"],
            0
        );
        assert_eq!(
            brief["proof_gap_coverage"]["claims_without_current_annotation_sources"],
            1
        );
        assert_eq!(
            brief["next_testable_requirements"][0]["annotation_freshness"]["state"],
            "source_stale"
        );
        let small = store.research_brief("claim:test", 1000, 0).await.unwrap();
        assert!(json_char_count(&small) <= 1000);
        assert_eq!(small["proof_gap_coverage"]["complete"], false);
        assert_eq!(
            small["proof_gap_coverage"]["claims_without_current_annotation_sources"],
            1
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn moved_prime_boundary_is_an_explicit_external_source() {
        let reference = json!({"source_id":"source:old","path":"removed.tex","sha256":"old","source_url":"https://github.com/example/prime/blob/commit/paper.tex","reference_role":"external_scope_boundary","source_location_state":"moved_to_external_repository"});
        assert!(theory_pattern_source(&reference).is_none());
        let dossier = json!({"node":{"id":"pattern:test","type":"theory_pattern","attrs":{}},"related_nodes":[],"edges":[],"sources":[],"external_sources":[reference]});
        let compact = compact_dossier_projection(&dossier, 4_000, 0);
        assert_eq!(
            compact["external_sources"][0]["reference_role"],
            "external_scope_boundary"
        );
        assert_eq!(
            compact["external_sources"][0]["source_url"],
            dossier["external_sources"][0]["source_url"]
        );
        assert!(json_char_count(&compact) <= 4_000);
    }

    #[test]
    fn compact_dossier_retains_annotation_source_freshness() {
        let dossier = json!({"node":{"id":"claim:test","type":"claim","attrs":{"enrichment_state":"enriched","cache_freshness":{"state":"source_stale","sources_current":false,"stale":1},"semantic":{"unknowns":["Old packet did not contain this proof"]}}},"related_nodes":[],"edges":[],"sources":[]});
        let compact = compact_dossier_projection(&dossier, 4_000, 0);
        assert_eq!(
            compact["node"]["attrs"]["cache_freshness"],
            dossier["node"]["attrs"]["cache_freshness"]
        );
        assert!(json_char_count(&compact) <= 4_000);
    }

    #[tokio::test]
    async fn minimum_size_claim_brief_keeps_unavailable_gap_coverage_explicit() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-brief-min-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({"meta":{"schema_version":1},"nodes":[{
            "id":"claim:AX.P1.01","type":"claim","label":"P1 claim",
            "summary":"A recorded claim.","attrs":{"enrichment_state":"error","ledger":{"status":"O"},
            "assessment":{"category":"open","scope":"Recorded claim scope"},"semantic":{"proof_gaps":[]}}
        }],"edges":[],"sources":[]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let brief = store
            .research_brief("claim:AX.P1.01", 1000, 0)
            .await
            .unwrap();
        assert!(json_char_count(&brief) <= 1000);
        assert_eq!(brief["enrichment_state"], "error");
        assert_eq!(brief["counts"]["recorded_proof_gaps"], Value::Null);
        assert_eq!(brief["proof_gap_coverage"]["complete"], false);
        assert_eq!(
            brief["proof_gap_coverage"]["claims_without_enriched_gap_audit"],
            1
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn compare_aspects_reports_shared_explicit_links_without_concluding_proof() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-compare-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({"nodes":[
            {"id":"claim:a","type":"claim","label":"A","attrs":{"assessment":{"category":"mixed"},"semantic":{"derivations":[{"upstream_claim_ids":["claim:seed"]}],"seed_ids":["seed:x"],"input_audit":[{"name":"c","origin":"free_parameter"}]}}},
            {"id":"claim:b","type":"claim","label":"B","attrs":{"assessment":{"category":"open"},"semantic":{"derivations":[{"upstream_claim_ids":["claim:seed"]}],"seed_ids":["seed:x"],"input_audit":[{"name":"d","origin":"hardcoded"}]}}},
            {"id":"topic:alpha","type":"topic","label":"Alpha","attrs":{}},
            {"id":"topic:flavor","type":"topic","label":"Flavor","attrs":{}},
            {"id":"claim:alpha-member","type":"claim","label":"Alpha member","attrs":{"assessment":{"category":"mixed"},"semantic":{"topic_ids":["topic:alpha"],"seed_ids":["seed:phi0"],"derivations":[{"upstream_claim_ids":["claim:shared"]}]}}},
            {"id":"claim:flavor-member","type":"claim","label":"Flavor member","attrs":{"assessment":{"category":"conditional"},"semantic":{"topic_ids":["topic:flavor"],"seed_ids":["seed:phi0"],"derivations":[{"upstream_claim_ids":["claim:shared"]}]}}},
            {"id":"structure:shared-seed","type":"structure","label":"Shared seed φ0/λ","attrs":{}}
        ],"edges":[
            {"id":"about-structure-alpha","src":"structure:shared-seed","dst":"topic:alpha","type":"about","attrs":{}},
            {"id":"about-structure-flavor","src":"structure:shared-seed","dst":"topic:flavor","type":"about","attrs":{}},
            {"id":"shares-seed-alpha","src":"topic:alpha","dst":"structure:shared-seed","type":"shares_seed","attrs":{}},
            {"id":"shares-seed-flavor","src":"topic:flavor","dst":"structure:shared-seed","type":"shares_seed","attrs":{}},
            {"id":"about-alpha-claim","src":"topic:alpha","dst":"claim:alpha-member","type":"about","attrs":{}},
            {"id":"about-flavor-claim","src":"topic:flavor","dst":"claim:flavor-member","type":"about","attrs":{}}
        ],"sources":[]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let compared = store
            .compare_aspects(&["claim:a".into(), "claim:b".into()], 12_000)
            .await
            .unwrap();
        assert_eq!(
            compared["shared_explicit"]["upstream_claim_ids"][0],
            "claim:seed"
        );
        assert_eq!(compared["shared_explicit"]["seeds"][0], "seed:x");
        assert!(compared["notice"]
            .as_str()
            .unwrap()
            .contains("do not establish"));
        assert_eq!(
            compared["items"][0]["assumptions"][0]["origin"],
            "free_parameter"
        );
        assert_eq!(compared["items"][1]["ledger_status"], Value::Null);
        assert_eq!(
            compared,
            store
                .compare_aspects(&["claim:a".into(), "claim:b".into()], 12_000)
                .await
                .unwrap(),
            "HashSet-backed fields must serialize deterministically"
        );
        let topic_comparison = store
            .compare_aspects(&["topic:alpha".into(), "topic:flavor".into()], 12_000)
            .await
            .unwrap();
        assert!(topic_comparison["shared_explicit"]["recorded_graph_nodes"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["id"] == "structure:shared-seed"));
        assert!(topic_comparison["shared_explicit"]["recorded_graph_links"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["edge_type"] == "shares_seed"
                && item["node_id"] == "structure:shared-seed"));
        assert_eq!(topic_comparison["shared_explicit"]["seeds"][0], "seed:phi0");
        assert_eq!(
            topic_comparison["shared_explicit"]["upstream_claim_ids"][0],
            "claim:shared"
        );
        assert_eq!(
            topic_comparison["items"][0]["topic_contents"]["claim_count"],
            1
        );
        let compact_comparison = store
            .compare_aspects(&["topic:alpha".into(), "topic:flavor".into()], 1000)
            .await
            .unwrap();
        assert!(json_char_count(&compact_comparison) <= 1000);
        assert_eq!(compact_comparison["ids"].as_array().unwrap().len(), 2);
        assert_eq!(compact_comparison["items"].as_array().unwrap().len(), 2);
        assert!(
            compact_comparison["shared_explicit"]["recorded_graph_links"]
                .as_array()
                .unwrap()
                .iter()
                .any(|item| item["edge_type"] == "shares_seed")
        );
        assert!(
            compact_comparison["shared_explicit"]["recorded_graph_nodes"]
                .as_array()
                .unwrap()
                .iter()
                .any(|item| item["id"] == "structure:shared-seed")
        );
        assert!(compact_comparison["status_comparison"].as_str().is_some());
        assert!(compact_comparison["items"]
            .as_array()
            .unwrap()
            .iter()
            .all(|item| {
                item.as_object().is_some_and(|fields| {
                    fields.contains_key("id") && fields.values().all(|field| !field.is_null())
                })
            }));
        assert!(compact_comparison["items"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["topic_contents"]["assessment_counts"].is_object()));
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn compare_aspects_keeps_recorded_contract_assumptions_scope_and_correction() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-compare-contracts-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({
            "nodes":[
                {"id":"contract:native","type":"attempt","label":"Native sector calibration","summary":"A bounded calibration record.","attrs":{"verdict_raw":"PARTIAL","original_contract":{"verdict_raw":"PARTIAL","scope":"Only the declared native sector and baseline subtraction.","assumptions":["The native H baseline is fixed before subtraction.","D^k preparation is held fixed."],"positive_results":["The shifted response is algebraically invisible in this setup."],"leaves_open":["Admissibility in the original profile-free functional remains open."],"original_source_scope_correction":"Algebraic invisibility does not establish admissibility in the original functional.","source_refs":["source:native"]}}},
                {"id":"contract:neutral","type":"attempt","label":"Neutral port identifiability","summary":"A separate bounded identifiability record.","attrs":{"verdict_raw":"PARTIAL","original_contract":{"verdict_raw":"PARTIAL","scope":"Rank-one neutral port with fixed output map.","assumptions":["The output map is held fixed.","The reference channel is calibrated."],"positive_results":["The same output admits a nonzero D5 coefficient."],"leaves_open":["A carrier-resolved observable is still required."],"source_refs":["source:neutral"]}}},
                {"id":"contract:shared-observation","type":"contract","label":"Shared observation contract","summary":"A common explicit record.","attrs":{}}
            ],
            "edges":[
                {"src":"contract:native","dst":"contract:shared-observation","type":"references"},
                {"src":"contract:neutral","dst":"contract:shared-observation","type":"references"}
            ],
            "sources":[
                {"id":"source:native","path":"native/PROOF.md","sha256":"a"},
                {"id":"source:neutral","path":"neutral/PROOF.md","sha256":"b"}
            ]
        });
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let compared = store
            .compare_aspects(
                &["contract:native".into(), "contract:neutral".into()],
                12_000,
            )
            .await
            .unwrap();
        let native = &compared["items"][0];
        let neutral = &compared["items"][1];
        assert_eq!(native["verdict"], "PARTIAL");
        assert_eq!(
            native["scope"],
            "Only the declared native sector and baseline subtraction."
        );
        assert_eq!(native["assumptions_total"], 2);
        assert_eq!(native["assumptions"].as_array().unwrap().len(), 2);
        assert_eq!(
            native["original_source_scope_correction"],
            "Algebraic invisibility does not establish admissibility in the original functional."
        );
        assert_eq!(
            native["recorded_conclusions"]["positive_results"][0],
            "The shifted response is algebraically invisible in this setup."
        );
        assert_eq!(neutral["verdict"], "PARTIAL");
        assert_eq!(neutral["assumptions_total"], 2);
        assert_eq!(neutral["assumptions"].as_array().unwrap().len(), 2);
        assert_eq!(
            neutral["recorded_conclusions"]["missing_or_open"][0],
            "A carrier-resolved observable is still required."
        );
        assert!(compared["shared_explicit"]["recorded_graph_nodes"]
            .as_array()
            .unwrap()
            .iter()
            .any(|node| node["id"] == "contract:shared-observation"));
        assert!(json_char_count(&compared) <= 12_000);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn research_packets_prioritize_load_bearing_inputs_over_benign_values() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-input-priority-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let benign = ["NY", "ZERO_CUT", "SX", "SY", "SZ", "MASS"]
            .into_iter()
            .map(|name| json!({"name":name,"origin":"upstream_result","impact":"benign","value":"recorded output"}))
            .collect::<Vec<_>>();
        let mut inputs = benign;
        inputs.extend([
            json!({"name":"eta_external","origin":"free_parameter","impact":"conditional","value":"chosen scale","used_for":"sets source normalization","source_refs":["source:eta"]}),
            json!({"name":"unused_c3v","origin":"unknown","impact":"missing_derivation","value":"no source use recorded","used_for":"script argument is unused"}),
            json!({"name":"self_input","origin":"free_parameter","impact":"circularity_risk","value":"fed from current output"}),
        ]);
        let graph = json!({"nodes":[
            {"id":"claim:inputs","type":"claim","label":"Input-sensitive result","summary":"Scoped output","attrs":{"assessment":{"category":"conditional","scope":"test only"},"enrichment_state":"enriched","semantic":{"input_audit":inputs}}},
            {"id":"claim:other","type":"claim","label":"Other claim","summary":"Other","attrs":{"assessment":{"category":"open"},"semantic":{}}}
        ],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let context = store
            .research_context("claim:inputs", 12_000, 0)
            .await
            .unwrap();
        assert_eq!(context["assumption_records"][0]["name"], "self_input");
        assert_eq!(context["assumption_records"][1]["name"], "unused_c3v");
        assert_eq!(context["assumption_records"][2]["name"], "eta_external");
        assert_eq!(context["coverage"]["assumption_records_total"], 3);
        let comparison = store
            .compare_aspects(&["claim:inputs".into(), "claim:other".into()], 8_000)
            .await
            .unwrap();
        assert_eq!(
            comparison["items"][0]["assumptions"][0]["name"],
            "self_input"
        );
        assert_eq!(
            comparison["items"][0]["assumptions"][1]["name"],
            "unused_c3v"
        );
        assert_eq!(
            comparison["items"][0]["assumptions"][2]["name"],
            "eta_external"
        );
        assert_eq!(comparison["items"][0]["assumptions_omitted"], 0);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn compact_compare_budget_keeps_scope_inputs_and_derivation_per_item() {
        let verbose = "x".repeat(12_000);
        let items: Vec<Value> = ["claim:BULK.PRIME_COMPLETION.01", "claim:GRAV.NONCIRCULAR.01", "contract:source-observation-identifiability-20261013"]
            .iter()
            .map(|id| json!({
                "id":id,
                "type":if id.starts_with("contract:") {"attempt"} else {"claim"},
                "label":format!("Evidence for {id}"),
                "scope":"Finite indexed core only; no continuum conclusion.",
                "summary":verbose,
                "assessment":{"category":"conditional","scope":"Finite indexed core only"},
                "ledger_status":"Partial; additional conditions remain open",
                "enrichment_state":"enriched",
                "upstream_claim_ids":["claim:upstream-1"],
                "seeds":["object:shared-input"],
                "assumptions":[{"name":"normalization","origin":"free_parameter","value":"chosen scale","impact":"changes source response","source_ids":["source:1"]}],
                "derivations":[{"id":"derivation:1","label":"Conditional response bound","premise_ids":["claim:premise-1"],"upstream_claim_ids":["claim:upstream-1"],"result_ids":["result:1"],"argument":"Only under the recorded finite-index assumption.","scope":"Finite indexed core"}],
                "source_ids":["source:1"],
                "recorded_graph_links":[],
                "topic_contents":null
            }))
            .collect();
        let value = json!({
            "ids":["claim:BULK.PRIME_COMPLETION.01","claim:GRAV.NONCIRCULAR.01","contract:source-observation-identifiability-20261013"],
            "items":items,
            "shared_explicit":{"upstream_claim_ids":[],"seeds":[],"recorded_graph_nodes":[],"recorded_graph_links":[]},
            "shared_source_ids":[],
            "status_comparison":"Recorded labels only.",
            "divergent_assumptions":"Compare each item's recorded input list.",
            "notice":"Comparison aid only.",
            "max_chars":10_000
        });
        let compact = fit_compare_budget(value, 10_000);
        assert!(json_char_count(&compact) <= 10_000);
        assert_eq!(compact["items"].as_array().unwrap().len(), 3);
        for item in compact["items"].as_array().unwrap() {
            assert!(item["scope"]
                .as_str()
                .is_some_and(|scope| scope.contains("Finite indexed core")));
            assert!(item["summary"]
                .as_str()
                .is_some_and(|summary| !summary.is_empty()));
            assert_eq!(item["assumptions"][0]["origin"], "free_parameter");
            assert_eq!(item["derivations"][0]["premise_ids"][0], "claim:premise-1");
            assert_eq!(
                item["derivations"][0]["argument"],
                "Only under the recorded finite-index assumption."
            );
        }
    }

    #[tokio::test]
    async fn gravity_brief_prioritizes_specific_seam_bridge_over_generic_contracts() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-gravity-rank-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({
            "nodes":[
                {"id":"topic:gravity","type":"topic","label":"Gravitation und Raumzeit","summary":"Dynamische Materie, lokale Gravitation und Geometrie bleiben offen.","attrs":{"keywords":["Gravitation","Einstein","ADM","Spin-2","Metrik","Raumzeit"]}},
                {"id":"claim:GRAV.ENTROPY.01","type":"claim","label":"Entropy normalization","attrs":{"assessment":{"category":"mixed"}}},
                {"id":"claim:GRAV.LOCAL.01","type":"claim","label":"Local source response","attrs":{"assessment":{"category":"open"}}},
                {"id":"claim:GENERIC.01","type":"claim","label":"Generic compiler relation","attrs":{"assessment":{"category":"recorded_established"}}},
                {"id":"contract:compiler-bond-response-20260919","type":"attempt","label":"Compiler bond response","attrs":{"verdict_raw":"PARTIAL","question":"Generic compiler bond observable contract","claims_tested":["claim:GENERIC.01"]}},
                {"id":"contract:seam-gravity-source-bridge-20261013","type":"attempt","label":"Does the seam gravity source bridge uniquely select a geometric horizon action?","attrs":{"verdict_raw":"PARTIAL","question":"Does finite modular covariance establish a geometric horizon action, entropy normalization, and local source area?","claims_tested":["claim:GRAV.ENTROPY.01","claim:GRAV.LOCAL.01"]}}
            ],
            "edges":[
                {"id":"about:g1","src":"topic:gravity","dst":"claim:GRAV.ENTROPY.01","type":"about","attrs":{}},
                {"id":"about:g2","src":"topic:gravity","dst":"claim:GRAV.LOCAL.01","type":"about","attrs":{}},
                {"id":"about:g3","src":"topic:gravity","dst":"claim:GENERIC.01","type":"about","attrs":{}},
                {"id":"test:generic","src":"contract:compiler-bond-response-20260919","dst":"claim:GENERIC.01","type":"tests","attrs":{}},
                {"id":"test:gravity1","src":"contract:seam-gravity-source-bridge-20261013","dst":"claim:GRAV.ENTROPY.01","type":"tests","attrs":{}},
                {"id":"test:gravity2","src":"contract:seam-gravity-source-bridge-20261013","dst":"claim:GRAV.LOCAL.01","type":"tests","attrs":{}}
            ],
            "sources":[]
        });
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let brief = store
            .research_brief("topic:gravity", 12_000, 0)
            .await
            .unwrap();
        assert_eq!(
            brief["attempts"]["items"][0]["id"],
            "contract:seam-gravity-source-bridge-20261013"
        );
        assert!(brief["attempts"]["items"][0]["question_or_scope"]
            .as_str()
            .unwrap()
            .contains("local source area"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn inline_pdf_excerpt_checks_repo_path_and_hash_first() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-pdf-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let pdf_bytes = b"%PDF-1.4\0binary source bytes";
        fs::write(root.join("paper.pdf"), pdf_bytes).unwrap();
        fs::write(root.join(".env"), b"private value").unwrap();
        let env_bytes = b"private value";
        let graph = json!({"nodes":[],"edges":[],"sources":[
            {"id":"source:pdf","path":"paper.pdf","sha256":format!("{:x}",Sha256::digest(pdf_bytes)),"page":3,"text":"supplied extracted quote"},
            {"id":"source:env","path":".env","sha256":format!("{:x}",Sha256::digest(env_bytes)),"text":"must not be served"},
            {"id":"source:escape","path":"../outside.txt","sha256":"ignored","text":"must not be served"}
        ]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let source = store.source("source:pdf", None, None).await.unwrap();
        assert_eq!(source["text"], "supplied extracted quote");
        assert_eq!(source["page"], 3);
        assert!(store
            .source("source:env", None, None)
            .await
            .unwrap_err()
            .to_string()
            .contains("unsafe"));
        assert!(store
            .source("source:escape", None, None)
            .await
            .unwrap_err()
            .to_string()
            .contains("unsafe"));
        fs::write(root.join("paper.pdf"), b"changed binary file").unwrap();
        assert!(store
            .source("source:pdf", None, None)
            .await
            .unwrap_err()
            .to_string()
            .contains("hash changed"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn trace_and_dossier_bound_large_fanout() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-fanout-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let mut nodes = vec![json!({"id":"topic:root","type":"topic","label":"Root","attrs":{}})];
        let mut edges = Vec::new();
        for i in 0..220 {
            let id = format!("result:{i}");
            nodes.push(json!({"id":id,"type":"result","label":id,"attrs":{}}));
            edges.push(json!({"id":format!("edge:{i}"),"src":"topic:root","dst":id,"type":"related","attrs":{}}));
        }
        let graph = json!({"nodes":nodes,"edges":edges,"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let trace = store.trace("topic:root", "outbound", 1).await.unwrap();
        assert_eq!(trace["nodes"].as_array().unwrap().len(), 150);
        assert_eq!(trace["truncated"], true);
        let dossier = store.dossier("topic:root").await.unwrap();
        assert_eq!(dossier["related_nodes"].as_array().unwrap().len(), 200);
        assert_eq!(dossier["edges"].as_array().unwrap().len(), 200);
        assert_eq!(dossier["truncated"], true);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn claim_trace_keeps_derivation_and_does_not_expand_topic_hub() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-topic-trace-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let mut nodes = vec![
            json!({"id":"claim:root","type":"claim","label":"Root claim","attrs":{}}),
            json!({"id":"topic:shared","type":"topic","label":"Shared topic","attrs":{}}),
            json!({"id":"result:r1","type":"result","label":"Result","attrs":{}}),
            json!({"id":"derivation:d1","type":"derivation","label":"Derivation","attrs":{}}),
            json!({"id":"assumption:a1","type":"assumption","label":"Premise","attrs":{}}),
        ];
        let mut edges = vec![
            json!({"id":"about-root","src":"claim:root","dst":"topic:shared","type":"about","attrs":{}}),
            json!({"id":"has-result","src":"claim:root","dst":"result:r1","type":"has_result","attrs":{}}),
            json!({"id":"has-derivation","src":"claim:root","dst":"derivation:d1","type":"has_derivation","attrs":{}}),
            json!({"id":"required-together","src":"assumption:a1","dst":"derivation:d1","type":"required_together","attrs":{}}),
            json!({"id":"derives","src":"derivation:d1","dst":"result:r1","type":"derives","attrs":{}}),
        ];
        for i in 0..220 {
            let id = format!("claim:unrelated-{i}");
            nodes.push(json!({"id":id,"type":"claim","label":id,"attrs":{}}));
            edges.push(json!({"id":format!("about-{i}"),"src":id,"dst":"topic:shared","type":"about","attrs":{}}));
        }
        let graph = json!({"nodes":nodes,"edges":edges,"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let trace = store.trace("claim:root", "both", 3).await.unwrap();
        let ids: HashSet<_> = trace["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|n| n["id"].as_str())
            .collect();
        let edge_types: HashSet<_> = trace["edges"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|e| e["type"].as_str())
            .collect();
        assert!(ids.contains("topic:shared"));
        assert!(ids.contains("derivation:d1"));
        assert!(ids.contains("assumption:a1"));
        assert!(!ids.contains("claim:unrelated-0"));
        assert!(edge_types.contains("required_together"));
        assert!(edge_types.contains("derives"));
        assert_eq!(trace["truncated"], false);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn claim_context_keeps_scoped_later_formal_link_ahead_of_seed_correlations() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-context-formal-link-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let overview = "Das Ledger führt den E8-Glue-Claim im Status Lattice. Der ältere v1-Skriptpfad berechnet die Cartan-Determinanten det(A3)=det(D5)=4 und setzt die Diskriminantenformwerte ein; diese Rechnungen allein leiten weder die abstrakten Z4-Gruppen noch eine konkrete Glue-Einbettung her. Ein späteres, unabhängiges FORM.GLUE.01-Resultat klassifiziert die isotropen Glue-Untergruppen im explizit deklarierten Modell Z4×Z4 mit qZ=(5x²+3y²) mod 8. Dieses separate formale Resultat beweist keine physische Auswahl von Carrier, Naht oder v1-Basis.";
        let mut nodes = vec![
            json!({"id":"claim:E8.GLU.01","type":"claim","label":"E8 glue lattice","summary":overview,"attrs":{"assessment":{"category":"recorded_established","scope":"finite lattice model"},"semantic":{"overview_de":overview,"results":[]}}}),
            json!({"id":"claim:FORM.GLUE.01","type":"claim","label":"Formal isotropic glue classification","summary":"Classification in the declared discriminant-form model.","attrs":{"assessment":{"category":"recorded_established","scope":"declared Z4×Z4 model"}}}),
        ];
        let mut edges = vec![json!({
            "id":"edge:formal-later",
            "src":"claim:E8.GLU.01",
            "dst":"claim:FORM.GLUE.01",
            "type":"mentions",
            "attrs":{"reason":"Later independent formal result: classifies isotropic glue subgroups in the declared carrier model; not a proof imported into the older script or a physical-selection result.","conditions":["Keep the formal qZ model explicit; det=4 alone is not Z4.","Do not identify the formal model with the v1 basis or physical selection."],"source_refs":["source:formal"]}
        })];
        for i in 0..20 {
            let id = format!("claim:seed-neighbor-{i:02}");
            nodes.push(json!({"id":id,"type":"claim","label":"Shared seed neighbor","summary":"A broad structural correlation.","attrs":{}}));
            edges.push(json!({"src":"claim:E8.GLU.01","dst":id,"type":"shares_seed","attrs":{}}));
        }
        let graph = json!({"nodes":nodes,"edges":edges,"sources":[{"id":"source:formal","path":"formal/FORM.GLUE.lean","sha256":"abc"}]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let context = store
            .research_context("claim:E8.GLU.01", 14_000, 0)
            .await
            .unwrap();
        let connection = context["recorded_connections"]
            .as_array()
            .unwrap()
            .iter()
            .find(|edge| edge["dst"] == "claim:FORM.GLUE.01")
            .unwrap();
        assert_eq!(connection["edge_type"], "mentions");
        assert!(connection["reason"]
            .as_str()
            .unwrap()
            .contains("Later independent formal result"));
        assert_eq!(connection["source_ids"][0], "source:formal");
        assert_eq!(
            context["what_we_have"][0]["statement"].as_str().unwrap(),
            overview
        );
        assert!(context["what_we_have"][0]["statement"]
            .as_str()
            .unwrap()
            .contains("qZ=(5x²+3y²)"));
        assert!(json_char_count(&context) <= 14_000);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn audit_counts_only_recorded_enriched_claim_annotations() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-audit-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let graph = json!({
            "meta":{"errors":[
                {"claim_id":"a","error":"old extraction failure superseded by enrichment"},
                {"claim_id":"b","error":"current extraction failure"},
                {"claim_id":"missing","error":"unlinked stale record"}
            ]},
            "nodes":[
                {"id":"script:a","type":"script","label":"Script A","attrs":{"static_inputs":{"syntactic_input_candidates":[{"text":"threshold = 7"}]}}},
                {"id":"result:a","type":"result","label":"Result A","attrs":{}},
                {"id":"claim:a","type":"claim","label":"Claim A","attrs":{"enrichment_state":"enriched","semantic":{"input_audit":[{"origin":"hardcoded","impact":"conditional"},{"origin":"upstream_result","impact":"benign"},{"origin":"unknown","impact":["circularity_risk","missing_derivation"]}],"proof_gaps":[{"kind":"missing_derivation"}]}}},
                {"id":"claim:b","type":"claim","label":"Claim B","attrs":{"enrichment_state":"error"}}
            ],
            "edges":[{"id":"e1","src":"script:a","dst":"result:a","type":"produces","attrs":{}},{"id":"e2","src":"result:a","dst":"claim:a","type":"used_by","attrs":{}}],
            "sources":[]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let status = store.status().await.unwrap();
        assert_eq!(status["audit"]["enriched_claims"], 1);
        assert_eq!(
            status["audit"]["input_audit_entries"]["by_origin"]["hardcoded"],
            1
        );
        assert_eq!(
            status["audit"]["input_audit_entries"]["by_origin"]["upstream_result"],
            1
        );
        assert_eq!(
            status["audit"]["input_audit_entries"]["by_origin"]["unknown"],
            1
        );
        assert_eq!(status["audit"]["by_impact"]["conditional"], 1);
        assert_eq!(status["audit"]["by_impact"]["circularity_risk"], 1);
        assert_eq!(status["audit"]["by_impact"]["missing_derivation"], 1);
        assert_eq!(status["audit"]["proof_gaps"], 1);
        assert_eq!(status["audit"]["input_audit_entries"]["total"], 3);
        assert!(status["meta"].get("errors").is_none());
        assert_eq!(status["meta"]["error_summary"]["recorded_total"], 3);
        assert_eq!(
            status["meta"]["error_summary"]["superseded_by_current_enrichment"],
            1
        );
        assert_eq!(
            status["meta"]["error_summary"]["unresolved_current_claims"],
            1
        );
        assert_eq!(status["meta"]["error_summary"]["unlinked_records"], 1);
        let script_audit = store.audit_inputs("script:a").await.unwrap();
        assert_eq!(
            script_audit["syntactic_input_candidates"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(script_audit["related_claims"][0]["id"], "claim:a");
        assert_eq!(
            script_audit["related_claims"][0]["proof_gaps"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        assert!(script_audit["interpretation_note"]
            .as_str()
            .unwrap()
            .contains("not proof gaps"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn overlay_preserves_main_claim_and_theme_keeps_mixed_open_and_attempt_separate() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-overlay-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let main = json!({"meta":{"schema_version":1},"nodes":[
            {"id":"topic:alpha","type":"topic","label":"Alpha","summary":"Topic text stays","attrs":{"source_refs":["source:primary"]}},
            {"id":"claim:est","type":"claim","label":"Established","summary":"Main summary","attrs":{"ledger":{"full":"kept"},"semantic":{"topic_ids":["topic:alpha"],"overview_de":"Main semantic"}}},
            {"id":"claim:mixed","type":"claim","label":"Mixed","summary":"Mixed claim","attrs":{"semantic":{"topic_ids":["topic:alpha"]}}},
            {"id":"claim:open","type":"claim","label":"Open","summary":"Open claim","attrs":{"assessment":{"category":"open","label":"Unresolved","scope":"full claim","authority":"ledger","notice":"open"}}}
        ],"edges":[{"id":"about-open","src":"topic:alpha","dst":"claim:open","type":"about"}],"sources":[{"id":"source:primary","path":"paper.md","sha256":"abc","text":"main source"}]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&main).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let overlay = json!({"nodes":[
            {"id":"claim:est","type":"claim","label":"overlay label must not replace","summary":"overlay summary must not replace","attrs":{"ledger":{"wrong":true},"semantic":{"wrong":true},"assessment":{"category":"recorded_established","label":"Scoped identity","scope":"T1-T8 identity","authority":"ledger","notice":"not physical closure"}}},
            {"id":"claim:mixed","type":"claim","attrs":{"assessment":{"category":"mixed","label":"Numeric E; physical O","scope":"numeric result only","authority":"ledger","notice":"physical derivation open"}}},
            {"id":"attempt:pass","type":"attempt","label":"Historical pass","attrs":{"status":"PASS"}},
            {"id":"structure:formal","type":"structure","label":"Formal structures","attrs":{}},
            {"id":"goal:physical","type":"goal","label":"Physical derivation","attrs":{}}
        ],"edges":[
            {"id":"theme-est","src":"topic:alpha","dst":"claim:est","type":"theme_links"},
            {"id":"theme-mixed","src":"topic:alpha","dst":"claim:mixed","type":"about"},
            {"id":"test-pass","src":"attempt:pass","dst":"claim:open","type":"tests"},
            {"id":"theme-structure","src":"topic:alpha","dst":"structure:formal","type":"about"},
            {"id":"theme-goal","src":"topic:alpha","dst":"goal:physical","type":"about"},
            {"id":"references-only","src":"claim:est","dst":"attempt:pass","type":"references"}
        ],"sources":[{"id":"source:context","path":"paper.md","sha256":"abc","text":"overlay source"}]});
        fs::write(
            root.join("context.json"),
            serde_json::to_vec(&overlay).unwrap(),
        )
        .unwrap();
        let theme = store.theme("topic:alpha", 100, 0, false).await.unwrap();
        assert_eq!(theme["solved"].as_array().unwrap().len(), 1);
        assert_eq!(theme["mixed"].as_array().unwrap().len(), 1);
        assert_eq!(theme["open"].as_array().unwrap().len(), 1);
        assert_eq!(theme["summary"]["recorded_established"], 1);
        assert_eq!(theme["summary"]["mixed"], 1);
        assert_eq!(theme["summary"]["open"], 1);
        assert_eq!(
            theme["pagination"]["recorded_established"]["truncated"],
            false
        );
        assert_eq!(theme["attempts"].as_array().unwrap().len(), 1);
        assert_eq!(theme["structures"][0]["id"], "structure:formal");
        assert_eq!(theme["goals"][0]["id"], "goal:physical");
        let dossier = store.dossier("claim:est").await.unwrap();
        assert_eq!(dossier["node"]["label"], "Established");
        assert_eq!(dossier["node"]["summary"], "Main summary");
        assert_eq!(dossier["node"]["attrs"]["ledger"]["full"], "kept");
        assert_eq!(
            dossier["node"]["attrs"]["semantic"]["overview_de"],
            "Main semantic"
        );
        assert_eq!(
            dossier["node"]["attrs"]["assessment"]["scope"],
            "T1-T8 identity"
        );
        let attempts = store
            .list_attempts("topic:alpha", None, 100, 0)
            .await
            .unwrap();
        assert_eq!(attempts["total"], 1);
        let status = store.status().await.unwrap();
        assert_eq!(status["assessment_counts"]["recorded_established"], 1);
        assert_eq!(status["assessment_counts"]["mixed"], 1);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn solution_chain_keeps_references_separate_from_reasoning_edges() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-chain-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let graph = json!({"nodes":[{"id":"claim:c","type":"claim","label":"Claim","attrs":{}},{"id":"assumption:a","type":"assumption","label":"Assumption","attrs":{}},{"id":"attempt:historical","type":"attempt","label":"Old pass","attrs":{"status":"PASS"}}],"edges":[{"id":"requires","src":"claim:c","dst":"assumption:a","type":"requires"},{"id":"ref","src":"claim:c","dst":"attempt:historical","type":"references"}],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let chain = store.solution_chain("claim:c", 4).await.unwrap();
        assert_eq!(chain["edges"].as_array().unwrap().len(), 1);
        assert_eq!(chain["edges"][0]["type"], "requires");
        assert_eq!(chain["references"].as_array().unwrap().len(), 1);
        let node_ids: HashSet<_> = chain["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|n| n["id"].as_str())
            .collect();
        assert!(chain["edges"]
            .as_array()
            .unwrap()
            .iter()
            .all(|e| node_ids.contains(e["src"].as_str().unwrap())
                && node_ids.contains(e["dst"].as_str().unwrap())));
        assert_eq!(chain["root"]["id"], "claim:c");
        assert_eq!(chain["traversal_complete"], true);
        assert_eq!(chain["all_premises_complete"], Value::Null);
        assert!(chain["interpretation"]["premise_completeness"]
            .as_str()
            .unwrap()
            .contains("Not assessed"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn solution_chain_projects_nodes_without_losing_premise_evidence() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-chain-compact-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let bulky = "original ledger text ".repeat(15_000);
        let premise_text = "the recorded normalized boundary condition ".repeat(40);
        let premise_sources: Vec<String> = (0..45).map(|i| format!("source:premise-{i}")).collect();
        let graph = json!({
            "nodes":[
                {"id":"claim:c","type":"claim","label":"Claim","summary":"Claim summary","attrs":{"ledger":{"full":bulky},"semantic":{"overview_de":bulky},"assessment":{"category":"conditional","scope":"within finite model"},"source_refs":["source:claim"]}},
                {"id":"derivation:d","type":"derivation","label":"Derivation","summary":"Derivation summary","attrs":{"argument":"The second premise fixes the boundary term.","scope":"within finite model","premise_ids":["assumption:a"],"source_refs":["source:derivation"],"limitations":["does not establish physical dynamics"]}},
                {"id":"assumption:a","type":"assumption","label":"Boundary condition","summary":"Explicit premise","attrs":{"text":premise_text.clone(),"kind":"boundary_condition","source_refs":premise_sources.clone()}}
            ],
            "edges":[
                {"id":"e1","src":"claim:c","dst":"derivation:d","type":"has_derivation"},
                {"id":"e2","src":"derivation:d","dst":"assumption:a","type":"requires"}
            ],
            "sources":[]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let chain = store.solution_chain("claim:c", 4).await.unwrap();
        let rendered = serde_json::to_string(&chain).unwrap();
        assert!(
            rendered.len() < 12_000,
            "chain response grew to {} bytes",
            rendered.len()
        );
        assert!(!rendered.contains("original ledger text"));
        assert_eq!(
            chain["root"]["attrs"]["assessment"]["scope"],
            "within finite model"
        );
        let derivation = chain["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|n| n["id"] == "derivation:d")
            .unwrap();
        assert_eq!(
            derivation["attrs"]["argument"],
            "The second premise fixes the boundary term."
        );
        assert_eq!(derivation["attrs"]["premise_ids"][0], "assumption:a");
        assert_eq!(
            derivation["attrs"]["limitations"][0],
            "does not establish physical dynamics"
        );
        let assumption = chain["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|n| n["id"] == "assumption:a")
            .unwrap();
        assert_eq!(assumption["attrs"]["text"], premise_text);
        assert_eq!(
            assumption["attrs"]["source_refs"].as_array().unwrap().len(),
            45
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn solution_chain_marks_incomplete_derivation_when_node_cap_cuts_premises() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-chain-cap-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let mut nodes = vec![
            json!({"id":"derivation:d","type":"derivation","label":"Shared derivation","attrs":{}}),
        ];
        let mut edges = Vec::new();
        for i in 0..100 {
            let id = format!("assumption:{i}");
            nodes.push(json!({"id":id,"type":"assumption","label":id,"attrs":{}}));
            edges.push(json!({"id":format!("e{i}"),"src":id,"dst":"derivation:d","type":"required_together"}));
        }
        let data = root.join("knowledge.json");
        fs::write(
            &data,
            serde_json::to_vec(&json!({"nodes":nodes,"edges":edges,"sources":[]})).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let chain = store.solution_chain("derivation:d", 2).await.unwrap();
        assert_eq!(chain["truncated"], true);
        assert_eq!(chain["all_premises_complete"], false);
        assert!(chain["reason"]
            .as_str()
            .unwrap()
            .contains("omitted premises"));
        assert!(chain["truncation_reasons"]
            .as_array()
            .unwrap()
            .iter()
            .any(|x| x == "node_limit"));
        assert!(chain["affected_derivations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|x| x["id"] == "derivation:d" && x["all_premises_complete"] == false));
        let nodes: HashSet<_> = chain["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|n| n["id"].as_str())
            .collect();
        assert!(chain["edges"]
            .as_array()
            .unwrap()
            .iter()
            .all(|e| nodes.contains(e["src"].as_str().unwrap())
                && nodes.contains(e["dst"].as_str().unwrap())));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn attempts_filter_contract_verdict_fields_for_pass_and_partial() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-verdict-filter-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let graph = json!({"nodes":[{"id":"claim:c","type":"claim","attrs":{}},{"id":"attempt:pass","type":"attempt","attrs":{"verdict_raw":"PASS"}},{"id":"attempt:partial","type":"attempt","attrs":{"original_contract":{"verdict":"PARTIAL"}}}],"edges":[{"id":"p","src":"attempt:pass","dst":"claim:c","type":"tests"},{"id":"q","src":"claim:c","dst":"attempt:partial","type":"tests"}],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let pass = store
            .list_attempts("claim:c", Some("PASS"), 100, 0)
            .await
            .unwrap();
        let partial = store
            .list_attempts("claim:c", Some("PARTIAL"), 100, 0)
            .await
            .unwrap();
        assert_eq!(pass["total"], 1);
        assert_eq!(pass["attempts"][0]["id"], "attempt:pass");
        assert_eq!(partial["total"], 1);
        assert_eq!(partial["attempts"][0]["id"], "attempt:partial");
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn claim_attempt_listing_includes_compact_contracts_and_orders_narrow_records_first() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-claim-contract-attempts-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let oversized = "large unrelated contract metadata ".repeat(5_000);
        let graph = json!({
            "nodes":[
                {"id":"claim:c","type":"claim","label":"Target claim","attrs":{}},
                {"id":"claim:other","type":"claim","label":"Other claim","attrs":{}},
                {"id":"contract:narrow","type":"contract","label":"Narrow source-specific test","summary":"One-claim scoped result.","attrs":{"assessment":{"category":"conditional","scope":"Only the target response."},"original_contract":{"verdict_raw":"PARTIAL","scope":"Only the target response.","positive_results":["A bounded response reconstruction."],"leaves_open":["Physical source selection remains open."],"source_refs":["source:narrow"]}}},
                {"id":"contract:wide","type":"contract","label":"Broad omnibus record","summary":"Tests more claims.","attrs":{"original_contract":{"verdict_raw":"PARTIAL","scope":"A multi-claim survey.","positive_results":["Many items summarized."],"metadata":oversized,"source_refs":["source:wide"]}}}
            ],
            "edges":[
                {"src":"contract:narrow","dst":"claim:c","type":"tests"},
                {"src":"contract:wide","dst":"claim:c","type":"tests"},
                {"src":"contract:wide","dst":"claim:other","type":"tests"}
            ],
            "sources":[{"id":"source:narrow","path":"narrow/PROOF.txt","sha256":"a"},{"id":"source:wide","path":"wide/PROOF.txt","sha256":"b"}]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let listed = store.list_attempts("claim:c", None, 10, 0).await.unwrap();
        assert_eq!(listed["total"], 2);
        assert_eq!(listed["attempts"][0]["id"], "contract:narrow");
        assert_eq!(listed["attempts"][1]["id"], "contract:wide");
        assert_eq!(
            listed["attempts"][0]["attrs"]["original_contract"]["positive_results"][0],
            "A bounded response reconstruction."
        );
        assert_eq!(
            listed["attempts"][0]["attrs"]["original_contract"]["leaves_open"][0],
            "Physical source selection remains open."
        );
        assert_eq!(
            listed["attempts"][0]["attrs"]["source_refs"][0],
            "source:narrow"
        );
        assert!(listed["claim_link_basis"]
            .as_str()
            .unwrap()
            .contains("not inferred chronological predecessors"));
        assert!(json_char_count(&listed) < 5_000);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn theme_categories_page_with_offset_and_report_all_totals() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-theme-pages-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let mut nodes = vec![json!({"id":"topic:t","type":"topic","attrs":{}})];
        let mut edges = Vec::new();
        for i in 0..5 {
            let id = format!("claim:m{i}");
            nodes.push(json!({"id":id,"type":"claim","attrs":{"assessment":{"category":"mixed"}}}));
            edges.push(json!({"id":format!("e{i}"),"src":"topic:t","dst":id,"type":"about"}));
        }
        let data = root.join("knowledge.json");
        fs::write(
            &data,
            serde_json::to_vec(&json!({"nodes":nodes,"edges":edges,"sources":[]})).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let page = store.theme("topic:t", 2, 2, false).await.unwrap();
        assert_eq!(page["mixed"].as_array().unwrap().len(), 2);
        assert_eq!(page["totals"]["mixed"], 5);
        assert_eq!(page["pagination"]["mixed"]["offset"], 2);
        assert_eq!(page["pagination"]["mixed"]["truncated"], true);
        assert_eq!(page["pagination"]["unknown"]["truncated"], false);
        fs::remove_dir_all(root).unwrap();
    }

    #[tokio::test]
    async fn topic_overview_compacts_cards_and_preserves_display_and_source_fields() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-theme-compact-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let bulky = "full original ledger and extraction details ".repeat(8_000);
        let graph = json!({
            "nodes":[
                {"id":"topic:gravity","type":"topic","label":"Gravity","summary":"Topic summary","attrs":{"topic_orientation":{"role_en":"Scope in English","role_de":"Geltungsbereich"}}},
                {"id":"claim:gravity","type":"claim","label":"Gravity result","summary":"Short card summary","attrs":{"ledger":{"status":"Conditional","claim":bulky.clone()},"assessment":{"category":"conditional","authority":"ledger","scope":"finite model","notice":"Recorded scope","source_refs":["source:ledger"]},"enrichment_state":"enriched","semantic":{"overview_de":"Korter Überblick","overview_en":"Short overview","interpretation_de":"Deutung","source_refs":[{"source_id":"source:gravity","line_start":12,"line_end":18}]},"translations":{"en":{"title":"Gravity title","overview":"Gravity overview","interpretation":"English interpretation"}},"source_refs":[{"source_id":"source:gravity","line_start":12,"line_end":18}]}},
                {"id":"attempt:gravity","type":"attempt","label":"Partial attempt","summary":"Attempt card","attrs":{"original_contract":{"question":bulky.clone(),"verdict":"PARTIAL","scope":"finite model","rules_out":["No universal closure"],"leaves_open":["Continuum limit"],"limitations":["finite example only"]},"source_refs":["source:attempt"]}}
            ],
            "edges":[
                {"id":"about-claim","src":"topic:gravity","dst":"claim:gravity","type":"about"},
                {"id":"about-attempt","src":"topic:gravity","dst":"attempt:gravity","type":"about"}
            ],
            "sources":[]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let theme = store.theme("topic:gravity", 100, 0, false).await.unwrap();
        let rendered = serde_json::to_string(&theme).unwrap();
        assert!(
            rendered.len() < 15_000,
            "theme payload grew to {} bytes",
            rendered.len()
        );
        assert_eq!(
            theme["conditional"][0]["attrs"]["assessment"]["scope"],
            "finite model"
        );
        assert_eq!(
            theme["conditional"][0]["attrs"]["semantic"]["overview_en"],
            "Short overview"
        );
        assert_eq!(
            theme["conditional"][0]["attrs"]["source_refs"][0]["line_start"],
            12
        );
        assert_eq!(
            theme["attempts"][0]["attrs"]["original_contract"]["leaves_open"][0],
            "Continuum limit"
        );
        assert!(!rendered.contains(&bulky));
        assert_eq!(
            theme["topic"]["attrs"]["topic_orientation"]["role_en"],
            "Scope in English"
        );
        let dossier = store.dossier("claim:gravity").await.unwrap();
        assert!(
            dossier["node"]["attrs"]["ledger"]["claim"]
                .as_str()
                .unwrap()
                .len()
                > 100_000
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn topic_overview_keeps_linked_attempts_when_related_claims_are_hidden() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-topic-attempts-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let graph = json!({
            "nodes":[
                {"id":"topic:universalraum","type":"topic","label":"Universalraum","attrs":{}},
                {"id":"claim:primary","type":"claim","label":"Primary","attrs":{"topic_roles":{"topic:universalraum":"primary"}}},
                {"id":"claim:related","type":"claim","label":"Related","attrs":{"topic_roles":{"topic:universalraum":"related"}}},
                {"id":"attempt:direct","type":"attempt","label":"Direct topic attempt","summary":"Direct","attrs":{"original_contract":{"verdict":"PARTIAL"}}},
                {"id":"attempt:inherited","type":"attempt","label":"Inherited topic attempt","summary":"Inherited","attrs":{"original_contract":{"verdict":"PASS"}}},
                {"id":"attempt:related-only","type":"attempt","label":"Only tested by related claim","attrs":{"original_contract":{"verdict":"FAIL"}}}
            ],
            "edges":[
                {"id":"direct-link","src":"topic:universalraum","dst":"attempt:direct","type":"about"},
                {"id":"inherited-link","src":"attempt:inherited","dst":"topic:universalraum","type":"about","attrs":{"method":"linked_claim_topics"}},
                {"id":"tests-hidden-claim","src":"attempt:related-only","dst":"claim:related","type":"tests"}
            ],
            "sources":[]
        });
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let theme = store
            .theme("topic:universalraum", 100, 0, false)
            .await
            .unwrap();
        assert_eq!(theme["summary"]["attempts"], 2);
        let basis: HashMap<String, String> = theme["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|item| {
                (
                    item["id"].as_str().unwrap().to_owned(),
                    item["topic_link_basis"].as_str().unwrap().to_owned(),
                )
            })
            .collect();
        assert_eq!(basis["attempt:direct"], "explicit_topic_edge");
        assert_eq!(
            basis["attempt:inherited"],
            "inherited_from_claim_topic_membership"
        );
        assert!(!basis.contains_key("attempt:related-only"));
        let attempts = store
            .list_attempts("topic:universalraum", None, 100, 0)
            .await
            .unwrap();
        assert_eq!(attempts["total"], 2);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn related_topic_claims_are_opt_in_and_do_not_drag_in_their_attempts() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-topic-roles-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(root.join("data")).unwrap();
        let graph = json!({"nodes":[
            {"id":"topic:e8","type":"topic","label":"E8","attrs":{}},
            {"id":"claim:e8-primary","type":"claim","label":"Primary E8 result","attrs":{"topic_roles":{"topic:e8":"primary"},"assessment":{"category":"recorded_established"}}},
            {"id":"claim:prime-related","type":"claim","label":"Related prime factor","attrs":{"topic_roles":{"topic:e8":"related"},"assessment":{"category":"open"}}},
            {"id":"claim:legacy-unclassified","type":"claim","label":"Legacy unclassified claim","attrs":{"assessment":{"category":"mixed"}}},
            {"id":"attempt:primary","type":"attempt","label":"Primary attempt","attrs":{"verdict_raw":"PASS"}},
            {"id":"attempt:related","type":"attempt","label":"Related attempt","attrs":{"verdict_raw":"PARTIAL"}},
                {"id":"attempt:direct-related","type":"attempt","label":"Directly linked review","attrs":{"verdict_raw":"PARTIAL"}},
                {"id":"attempt:inherited-about","type":"attempt","label":"Inherited topic link","attrs":{"verdict_raw":"PARTIAL"}}
        ],"edges":[
            {"id":"about:p","src":"topic:e8","dst":"claim:e8-primary","type":"about","attrs":{}},
            {"id":"about:r","src":"topic:e8","dst":"claim:prime-related","type":"about","attrs":{}},
            {"id":"about:l","src":"topic:e8","dst":"claim:legacy-unclassified","type":"about","attrs":{}},
            {"id":"tests:p","src":"attempt:primary","dst":"claim:e8-primary","type":"tests","attrs":{}},
            {"id":"tests:r","src":"attempt:related","dst":"claim:prime-related","type":"tests","attrs":{}},
                {"id":"about:direct","src":"topic:e8","dst":"attempt:direct-related","type":"about","attrs":{}},
                {"id":"about:inherited","src":"topic:e8","dst":"attempt:inherited-about","type":"about","attrs":{"method":"linked_claim_topics"}}
        ],"sources":[]});
        let data = root.join("data/knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let primary = store.theme("topic:e8", 100, 0, false).await.unwrap();
        assert_eq!(primary["related_claims_total"], 1);
        assert!(primary["membership_notice"]
            .as_str()
            .unwrap()
            .contains("no topic role"));
        assert_eq!(primary["summary"]["recorded_established"], 1);
        assert_eq!(primary["summary"]["open"], 0);
        assert_eq!(primary["summary"]["mixed"], 1);
        let primary_attempts: Vec<_> = primary["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x["id"].as_str().unwrap())
            .collect();
        assert!(primary_attempts.contains(&"attempt:primary"));
        assert!(primary_attempts.contains(&"attempt:direct-related"));
        assert!(!primary_attempts.contains(&"attempt:related"));
        assert!(primary_attempts.contains(&"attempt:inherited-about"));
        let expanded = store.theme("topic:e8", 100, 0, true).await.unwrap();
        assert_eq!(expanded["summary"]["open"], 1);
        assert!(expanded["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .any(|x| x["id"] == "attempt:related"));
        assert!(expanded["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .any(|x| x["id"] == "attempt:inherited-about"));
        let brief = store.research_brief("topic:e8", 12_000, 0).await.unwrap();
        assert_eq!(brief["counts"]["claims"], 2);
        assert_eq!(brief["counts"]["related_claims_total"], 1);
        assert_eq!(brief["counts"]["attempts"], 2);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn context_overlay_is_reloaded_when_added_after_startup() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-reload-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let graph = json!({"nodes":[{"id":"topic:t","type":"topic","label":"Topic","attrs":{}}],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&graph).unwrap()).unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        assert_eq!(store.status().await.unwrap()["counts"]["nodes"], 1);
        let overlay = json!({"nodes":[{"id":"claim:new","type":"claim","label":"New","attrs":{"assessment":{"category":"open"}}}],"edges":[{"id":"about","src":"topic:t","dst":"claim:new","type":"about"}],"sources":[]});
        fs::write(
            root.join("context.json"),
            serde_json::to_vec(&overlay).unwrap(),
        )
        .unwrap();
        let status = store.status().await.unwrap();
        assert_eq!(status["counts"]["nodes"], 2);
        assert_eq!(status["assessment_counts"]["open"], 1);
        assert_eq!(
            store.theme("topic:t", 100, 0, false).await.unwrap()["open"][0]["id"],
            "claim:new"
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn corpus_overlay_adds_original_sections_and_skips_unknown_edges_safely() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-corpus-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let main = json!({"nodes":[
            {"id":"claim:C.01","type":"claim","label":"Main claim","attrs":{"ledger":{"status":"[O]"},"semantic":{"title_de":"Authoritative"}}}
        ],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&main).unwrap()).unwrap();
        let context = json!({"nodes":[
            {"id":"attempt:A.01","type":"attempt","label":"Context attempt","attrs":{"assessment":{"category":"historical"}}}
        ],"edges":[],"sources":[]});
        fs::write(
            root.join("context.json"),
            serde_json::to_vec(&context).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let corpus = json!({"nodes":[
            {"id":"claim:C.01","type":"claim","label":"Corpus replacement","attrs":{"ledger":{"status":"[E]"},"semantic":{"title_de":"Wrong"},"corpus_hint":"retained as extra"}},
            {"id":"attempt:A.01","type":"experiment_record","label":"Corpus attempt replacement","attrs":{"assessment":{"category":"solved"},"corpus_hint":"extra"}},
            {"id":"paper_section:S.01","type":"paper_section","label":"Original source section","summary":"A section cited by an experiment record.","attrs":{"path":"paper.tex","line_start":12,"line_end":20}}
        ],"edges":[
            {"id":"known","src":"paper_section:S.01","dst":"claim:C.01","type":"explicit_claim_mention"},
            {"id":"unknown","src":"paper_section:S.01","dst":"claim:missing","type":"explicit_claim_mention"}
        ],"sources":[{"id":"source:corpus","path":"paper.tex","sha256":"abc","line_start":12,"line_end":20}]});
        fs::write(
            root.join("corpus.json"),
            serde_json::to_vec(&corpus).unwrap(),
        )
        .unwrap();
        let graph = store.snapshot().await.unwrap();
        let claim = graph
            .nodes
            .iter()
            .find(|node| node.id == "claim:C.01")
            .unwrap();
        let attempt = graph
            .nodes
            .iter()
            .find(|node| node.id == "attempt:A.01")
            .unwrap();
        let section = graph
            .nodes
            .iter()
            .find(|node| node.id == "paper_section:S.01")
            .unwrap();
        assert_eq!(claim.attrs["ledger"]["status"], "[O]");
        assert_eq!(claim.attrs["semantic"]["title_de"], "Authoritative");
        assert_eq!(claim.attrs["corpus_hint"], "retained as extra");
        assert_eq!(attempt.kind, "attempt");
        assert_eq!(attempt.label, "Context attempt");
        assert_eq!(attempt.attrs["assessment"]["category"], "historical");
        assert_eq!(attempt.attrs["corpus_hint"], "extra");
        assert_eq!(section.label, "Original source section");
        assert!(graph.edges.iter().any(|edge| edge.id == "known"));
        assert!(!graph.edges.iter().any(|edge| edge.id == "unknown"));
        assert_eq!(graph.meta["corpus_merge"]["nodes_added"], 1);
        assert_eq!(graph.meta["corpus_merge"]["nodes_merged"], 2);
        assert_eq!(
            graph.meta["corpus_merge"]["edges_skipped_unknown_endpoint"],
            1
        );
        assert_eq!(graph.sources.len(), 1);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn translations_are_fingerprint_checked_and_reloaded_without_touching_semantics() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-translations-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let semantic = json!({
            "title_de":"Titel",
            "overview_de":"Überblick",
            "interpretation_de":"Deutung",
            "results":[{"id":"R1","statement":"unchanged"}]
        });
        let data = root.join("knowledge.json");
        fs::write(
            &data,
            serde_json::to_vec(&json!({"nodes":[{"id":"claim:t","type":"claim","label":"Claim","summary":"Original","attrs":{"semantic":semantic,"ledger":{"status":"[O]"}}}],"edges":[],"sources":[]})).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let fingerprint = {
            let mut hasher = Sha256::new();
            hasher.update("Titel\nÜberblick\nDeutung".as_bytes());
            format!("{:x}", hasher.finalize())
        };
        let stale = json!({"nodes":[{"id":"claim:t","type":"claim","attrs":{"translations":{"source_semantic_sha256":"stale","en":{"title":"Wrong"}}}}],"edges":[],"sources":[]});
        fs::write(
            root.join("translations.json"),
            serde_json::to_vec(&stale).unwrap(),
        )
        .unwrap();
        let dossier = store.dossier("claim:t").await.unwrap();
        assert!(dossier["node"]["attrs"].get("translations").is_none());
        let fresh = json!({"nodes":[{"id":"claim:t","type":"claim","attrs":{"translations":{"source_semantic_sha256":fingerprint,"en":{"title":"Title"},"provenance":{"reviewed":false}}}}],"edges":[],"sources":[]});
        fs::write(
            root.join("translations.json"),
            serde_json::to_vec(&fresh).unwrap(),
        )
        .unwrap();
        let dossier = store.dossier("claim:t").await.unwrap();
        assert_eq!(
            dossier["node"]["attrs"]["translations"]["en"]["title"],
            "Title"
        );
        assert_eq!(
            dossier["node"]["attrs"]["semantic"]["results"][0]["statement"],
            "unchanged"
        );
        assert_eq!(dossier["node"]["attrs"]["ledger"]["status"], "[O]");
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn extraction_quality_overlay_is_metadata_and_self_edges_are_not_proof_edges() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-quality-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let semantic = json!({
            "title_de":"Original",
            "overview_de":"",
            "interpretation_de":""
        });
        let main = json!({"nodes":[
            {"id":"claim:X.01","type":"claim","label":"Claim","attrs":{"enrichment_state":"enriched","ledger":{"status":"[O]"},"semantic":semantic}},
            {"id":"claim:Y.01","type":"claim","label":"Independent predecessor","attrs":{}},
            {"id":"derivation:X.01:D1","type":"derivation","attrs":{"upstream_claim_ids":["X.01","Y.01"]}},
            {"id":"derivation:opaque","type":"derivation","attrs":{"upstream_claim_ids":["Y.01"]}}
        ],"edges":[
            {"id":"self","src":"claim:X.01","dst":"derivation:X.01:D1","type":"uses_claim_under_scope"},
            {"id":"upstream","src":"claim:Y.01","dst":"derivation:X.01:D1","type":"uses_claim_under_scope"},
            {"id":"explicit-owner","src":"claim:X.01","dst":"derivation:opaque","type":"has_derivation"},
            {"id":"opaque-upstream","src":"claim:Y.01","dst":"derivation:opaque","type":"uses_claim_under_scope"},
            {"id":"opaque-self","src":"claim:X.01","dst":"derivation:opaque","type":"uses_claim_under_scope"}
        ],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&main).unwrap()).unwrap();
        let annotation_fingerprint = {
            let mut hasher = Sha256::new();
            hasher.update(serde_json::to_vec(&semantic).unwrap());
            format!("{:x}", hasher.finalize())
        };
        let summary_fingerprint = {
            let mut hasher = Sha256::new();
            hasher.update("Original\n\n".as_bytes());
            format!("{:x}", hasher.finalize())
        };
        let quality = json!({"nodes":[{"id":"claim:X.01","type":"claim","attrs":{"extraction_quality":{"state":"needs_review","annotation_sha256":annotation_fingerprint,"source_semantic_sha256":"stale-but-ignored","issues":[{"code":"self_upstream"}]}}}],"edges":[],"sources":[{"id":"source:quality-review","path":"experiments/review.md","sha256":"abc","line_start":1,"line_end":2}]});
        fs::write(
            root.join("quality.json"),
            serde_json::to_vec(&quality).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let graph = store.snapshot().await.unwrap();
        let claim = graph
            .nodes
            .iter()
            .find(|node| node.id == "claim:X.01")
            .unwrap();
        assert_eq!(claim.attrs["semantic"]["title_de"], "Original");
        assert_eq!(claim.attrs["ledger"]["status"], "[O]");
        assert_eq!(claim.attrs["extraction_quality"]["state"], "needs_review");
        assert!(graph
            .sources
            .iter()
            .any(|source| source.id == "source:quality-review"));
        assert_eq!(graph.edges[0].kind, "self_reference");
        assert_eq!(
            graph.edges[0].attrs["epistemic_role"],
            "non_evidentiary_self_reference"
        );
        assert_eq!(graph.edges[1].kind, "uses_claim_under_scope");
        assert_eq!(graph.edges[3].kind, "uses_claim_under_scope");
        assert_eq!(graph.edges[4].kind, "self_reference");
        assert!(graph.edges[1].attrs.get("epistemic_role").is_none());
        let stale_quality = json!({"nodes":[{"id":"claim:X.01","type":"claim","attrs":{"extraction_quality":{"state":"needs_review","annotation_sha256":"stale","source_semantic_sha256":summary_fingerprint,"issues":[{"code":"old_issue"}]}}}],"edges":[],"sources":[]});
        fs::write(
            root.join("quality.json"),
            serde_json::to_vec(&stale_quality).unwrap(),
        )
        .unwrap();
        let graph = store.snapshot().await.unwrap();
        let claim = graph
            .nodes
            .iter()
            .find(|node| node.id == "claim:X.01")
            .unwrap();
        assert_eq!(
            claim.attrs["extraction_quality"]["state"],
            "outdated_review"
        );
        assert!(claim.attrs["extraction_quality"].get("issues").is_none());
        assert_eq!(
            claim.attrs["extraction_quality"]["stale_issues"][0]["code"],
            "old_issue"
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn enriched_priority_overlay_promotes_pending_semantics_but_never_replaces_enriched() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-priority-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let main = json!({"nodes":[
            {"id":"claim:pending","type":"claim","label":"Pending","attrs":{"enrichment_state":"pending","ledger":{"original":true},"semantic":{"old":true},"provenance":{"reviewed":false}}},
            {"id":"claim:done","type":"claim","label":"Done","attrs":{"enrichment_state":"enriched","ledger":{"original":true},"semantic":{"kept":true},"provenance":{"reviewed":true}}}
        ],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&main).unwrap()).unwrap();
        let overlay = json!({"nodes":[
            {"id":"claim:pending","type":"claim","attrs":{"priority_annotation":{"enrichment_state":"enriched","semantic":{"promoted":true},"provenance":{"reviewed":false},"source_refs":["source:new"]},"ledger":{"replace":true}}},
            {"id":"claim:done","type":"claim","attrs":{"priority_annotation":{"enrichment_state":"enriched","semantic":{"downgrade":true},"provenance":{"reviewed":false}}}}
        ],"edges":[],"sources":[]});
        fs::write(
            root.join("context.json"),
            serde_json::to_vec(&overlay).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let pending = store.dossier("claim:pending").await.unwrap()["node"].clone();
        let done = store.dossier("claim:done").await.unwrap()["node"].clone();
        assert_eq!(pending["attrs"]["semantic"]["promoted"], true);
        assert_eq!(pending["attrs"]["ledger"]["original"], true);
        assert_eq!(pending["attrs"]["enrichment_state"], "enriched");
        assert_eq!(pending["attrs"]["source_refs"][0], "source:new");
        assert!(pending["attrs"].get("priority_annotation").is_none());
        assert_eq!(done["attrs"]["semantic"]["kept"], true);
        assert_eq!(done["attrs"]["provenance"]["reviewed"], true);
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn priority_replaces_healthy_semantics_only_with_an_exact_declared_fingerprint() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-priority-replace-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let old = json!({
            "title_de":"Old title",
            "overview_de":"Old overview",
            "interpretation_de":"Old interpretation",
            "results":[{"id":"old","statement":"old structure"}]
        });
        let replacement = json!({
            "title_de":"Reviewed title",
            "overview_de":"Reviewed overview",
            "interpretation_de":"Reviewed interpretation",
            "results":[{"id":"new","statement":"corrected structure"}]
        });
        let matching_hash = semantic_translation_fingerprint(&Node {
            attrs: json!({"semantic":old}),
            ..Node::default()
        })
        .unwrap();
        let main = json!({"nodes":[
            {"id":"claim:replace","type":"claim","label":"Replace","attrs":{"enrichment_state":"enriched","ledger":{"status":"[O]"},"semantic":old,"provenance":{"reviewed":true}}},
            {"id":"claim:retain","type":"claim","label":"Retain","attrs":{"enrichment_state":"enriched","ledger":{"status":"[E]"},"semantic":old,"provenance":{"reviewed":true}}}
        ],"edges":[],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&main).unwrap()).unwrap();
        let overlay = json!({"nodes":[
            {"id":"claim:replace","type":"claim","attrs":{"priority_annotation":{"enrichment_state":"enriched","replaces_semantic_sha256":matching_hash,"semantic":replacement,"provenance":{"reviewed":false},"source_refs":["source:priority"],"source_coverage":{"status":"current"}}}},
            {"id":"claim:retain","type":"claim","attrs":{"priority_annotation":{"enrichment_state":"enriched","replaces_semantic_sha256":"wrong-fingerprint","semantic":{"title_de":"must not replace"},"provenance":{"reviewed":false}}}},
            {"id":"derivation:replace:D1","type":"derivation","label":"Current local derivation","attrs":{"priority_owner":"claim:replace"}},
            {"id":"derivation:retain:D1","type":"derivation","label":"Stale local derivation","attrs":{"priority_owner":"claim:retain"}}
        ],"edges":[
            {"id":"replacement-edge","src":"claim:replace","dst":"derivation:replace:D1","type":"has_derivation","attrs":{"priority_owner":"claim:replace"}},
            {"id":"stale-edge","src":"claim:retain","dst":"derivation:retain:D1","type":"has_derivation","attrs":{"priority_owner":"claim:retain"}}
        ],"sources":[]});
        fs::write(
            root.join("context.json"),
            serde_json::to_vec(&overlay).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let graph = store.snapshot().await.unwrap();
        let replaced = graph
            .nodes
            .iter()
            .find(|node| node.id == "claim:replace")
            .unwrap();
        let retained = graph
            .nodes
            .iter()
            .find(|node| node.id == "claim:retain")
            .unwrap();
        assert_eq!(
            replaced.attrs["semantic"]["results"][0]["statement"],
            "corrected structure"
        );
        assert_eq!(replaced.attrs["ledger"]["status"], "[O]");
        assert_eq!(replaced.attrs["provenance"]["reviewed"], false);
        assert_eq!(
            retained.attrs["semantic"]["results"][0]["statement"],
            "old structure"
        );
        assert_eq!(retained.attrs["provenance"]["reviewed"], true);
        assert!(graph
            .nodes
            .iter()
            .any(|node| node.id == "derivation:replace:D1"));
        assert!(!graph
            .nodes
            .iter()
            .any(|node| node.id == "derivation:retain:D1"));
        assert!(graph.edges.iter().any(|edge| edge.id == "replacement-edge"));
        assert!(!graph.edges.iter().any(|edge| edge.id == "stale-edge"));
        fs::remove_dir_all(root).unwrap();
    }
    #[tokio::test]
    async fn priority_local_subgraph_is_imported_only_for_promoted_claims() {
        let root = std::env::temp_dir().join(format!(
            "semantic-explorer-priority-owner-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let main = json!({"nodes":[
            {"id":"claim:healthy","type":"claim","label":"Healthy","attrs":{"enrichment_state":"enriched","ledger":{"original":true},"semantic":{"main":true}}},
            {"id":"claim:pending","type":"claim","label":"Pending","attrs":{"enrichment_state":"pending","ledger":{"original":true},"semantic":{"old":true}}},
            {"id":"derivation:shared","type":"derivation","label":"Main D1","summary":"Authoritative D1","attrs":{"owner_claim":"claim:healthy"}}
        ],"edges":[],"sources":[]});
        let overlay = json!({"nodes":[
            {"id":"claim:healthy","type":"claim","attrs":{"priority_annotation":{"enrichment_state":"enriched","semantic":{"priority":true}}}},
            {"id":"claim:pending","type":"claim","attrs":{"priority_annotation":{"enrichment_state":"enriched","semantic":{"promoted":true}}}},
            {"id":"derivation:shared","type":"derivation","label":"Wrong priority D1","summary":"Must not replace main D1","attrs":{"priority_owner":"claim:healthy"}},
            {"id":"derivation:priority","type":"derivation","label":"Priority D1 for pending","attrs":{"priority_owner":"claim:pending"}}
        ],"edges":[
            {"id":"edge:healthy","src":"claim:healthy","dst":"derivation:shared","type":"has_derivation","attrs":{"priority_owner":"claim:healthy"}},
            {"id":"edge:pending","src":"claim:pending","dst":"derivation:priority","type":"has_derivation","attrs":{"priority_owner":"claim:pending"}}
        ],"sources":[]});
        let data = root.join("knowledge.json");
        fs::write(&data, serde_json::to_vec(&main).unwrap()).unwrap();
        fs::write(
            root.join("context.json"),
            serde_json::to_vec(&overlay).unwrap(),
        )
        .unwrap();
        let store = GraphStore::new(data, root.clone()).unwrap();
        let graph = store.snapshot().await.unwrap();
        let d1 = graph
            .nodes
            .iter()
            .find(|n| n.id == "derivation:shared")
            .unwrap();
        let d2 = graph
            .nodes
            .iter()
            .find(|n| n.id == "derivation:priority")
            .unwrap();
        assert_eq!(d1.label, "Main D1");
        assert_eq!(d1.summary, "Authoritative D1");
        assert_eq!(d2.label, "Priority D1 for pending");
        assert!(!graph.edges.iter().any(|e| e.id == "edge:healthy"));
        assert!(graph.edges.iter().any(|e| e.id == "edge:pending"));
        assert_eq!(
            graph
                .nodes
                .iter()
                .find(|n| n.id == "claim:healthy")
                .unwrap()
                .attrs["semantic"]["main"],
            true
        );
        assert_eq!(
            graph
                .nodes
                .iter()
                .find(|n| n.id == "claim:pending")
                .unwrap()
                .attrs["semantic"]["promoted"],
            true
        );
        fs::remove_dir_all(root).unwrap();
    }
}
