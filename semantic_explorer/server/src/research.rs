//! Durable, source-pinned research decisions. This layer never promotes ledger proofs.
use crate::service::{Edge, KnowledgeGraph, Node, Source};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{HashMap, HashSet},
    fs,
    path::{Component, Path, PathBuf},
};

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Text {
    pub de: String,
    pub en: String,
}
impl Text {
    fn get(&self, lang: &str) -> &str {
        if lang == "en" {
            &self.en
        } else {
            &self.de
        }
    }
    fn valid(&self) -> bool {
        !self.de.trim().is_empty() && !self.en.trim().is_empty()
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Premise {
    pub id: String,
    pub title: Text,
    pub statement: Text,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Citation {
    pub path: String,
    pub sha256: String,
    pub line_start: usize,
    pub line_end: usize,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Finding {
    pub id: String,
    pub title: Text,
    pub statement: Text,
    pub status: String,
    pub premise_ids: Vec<String>,
    pub scope: Text,
    pub sources: Vec<Citation>,
    pub attempt_ids: Vec<String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Task {
    #[serde(default = "construction_kind")]
    pub kind: String,
    pub id: String,
    pub title: Text,
    pub status: String,
    pub depends_on: Vec<String>,
    pub resolution_ids: Vec<String>,
    pub required_premise_ids: Vec<String>,
    pub success_criterion: Text,
    pub next_action: Text,
}
fn construction_kind() -> String {
    "construction".into()
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Strategy {
    pub id: String,
    pub title: Text,
    pub aliases: Vec<String>,
    pub outcome: String,
    pub premise_ids: Vec<String>,
    pub finding_ids: Vec<String>,
    pub explanation: Text,
    pub changes_required: Text,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Connection {
    pub id: String,
    pub target_problem_id: String,
    pub kind: String,
    pub premise_ids: Vec<String>,
    pub finding_ids: Vec<String>,
    pub explanation: Text,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Problem {
    pub id: String,
    pub title: Text,
    pub question: Text,
    pub topics: Vec<String>,
    pub goals: Vec<String>,
    pub accepted_premise_ids: Vec<String>,
    pub findings: Vec<Finding>,
    pub tasks: Vec<Task>,
    pub strategies: Vec<Strategy>,
    pub reviewed_attempt_ids: Vec<String>,
    #[serde(default)]
    pub connections: Vec<Connection>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Program {
    pub schema_version: u32,
    pub revision: String,
    pub premises: Vec<Premise>,
    pub problems: Vec<Problem>,
}

pub fn load(path: &Path) -> Result<Program> {
    if !path.exists() {
        return Ok(Program {
            schema_version: 1,
            revision: "empty".into(),
            premises: vec![],
            problems: vec![],
        });
    }
    let mut p: Program =
        serde_json::from_slice(&fs::read(path)?).context("invalid research program")?;
    validate(&p)?;
    p.revision.clear();
    p.revision = digest(&serde_json::to_vec(&p)?)[..24].to_owned();
    Ok(p)
}
fn require(condition: bool, message: impl AsRef<str>) -> Result<()> {
    if !condition {
        bail!("research validation: {}", message.as_ref())
    }
    Ok(())
}
pub fn validate(p: &Program) -> Result<()> {
    require(p.schema_version == 1, "unsupported schema_version")?;
    require(!p.revision.is_empty(), "revision required")?;
    let mut ids = HashSet::new();
    let mut register = |id: &str, prefix: &str| -> Result<()> {
        require(
            id.starts_with(prefix)
                && id.len() > prefix.len()
                && !id.chars().any(char::is_whitespace),
            format!("invalid ID {id}"),
        )?;
        require(ids.insert(id.to_owned()), format!("duplicate ID {id}"))
    };
    let premises: HashSet<_> = p.premises.iter().map(|x| &x.id).collect();
    for x in &p.premises {
        register(&x.id, "premise:")?;
        require(
            x.title.valid() && x.statement.valid(),
            format!("bilingual premise required: {}", x.id),
        )?;
    }
    for q in &p.problems {
        register(&q.id, "problem:")?;
        require(
            q.title.valid() && q.question.valid(),
            format!("bilingual question required: {}", q.id),
        )?;
        let fs: HashMap<_, _> = q.findings.iter().map(|f| (f.id.as_str(), f)).collect();
        let ts: HashMap<_, _> = q.tasks.iter().map(|t| (t.id.as_str(), t)).collect();
        for id in &q.accepted_premise_ids {
            require(premises.contains(id), format!("unknown premise {id}"))?;
        }
        for f in &q.findings {
            register(&f.id, "finding:")?;
            require(
                f.title.valid() && f.statement.valid() && f.scope.valid(),
                format!("bilingual finding/scope required: {}", f.id),
            )?;
            require(
                ["established_scoped", "conditional", "excluded_scoped"]
                    .contains(&f.status.as_str()),
                "invalid finding status",
            )?;
            require(
                !f.sources.is_empty(),
                format!("finding without sources: {}", f.id),
            )?;
            for id in &f.premise_ids {
                require(premises.contains(id), format!("unknown premise {id}"))?;
            }
            for s in &f.sources {
                require(
                    !s.path.is_empty()
                        && Path::new(&s.path)
                            .components()
                            .all(|c| matches!(c, Component::Normal(_))),
                    "source must be repository-relative without traversal",
                )?;
                require(
                    s.sha256.len() == 64 && s.sha256.bytes().all(|c| c.is_ascii_hexdigit()),
                    "source needs full SHA256",
                )?;
                require(
                    s.line_start > 0
                        && s.line_end >= s.line_start
                        && s.line_end - s.line_start < 200,
                    "source needs inclusive range of at most 200 lines",
                )?;
            }
        }
        for t in &q.tasks {
            register(&t.id, "task:")?;
            require(
                t.title.valid() && t.success_criterion.valid() && t.next_action.valid(),
                "bilingual task required",
            )?;
            require(
                ["open", "resolved", "excluded"].contains(&t.status.as_str()),
                "invalid task status",
            )?;
            require(
                ["construction", "assessment"].contains(&t.kind.as_str()),
                "invalid task kind",
            )?;
            require(
                t.status == "open" || !t.resolution_ids.is_empty(),
                "closed task needs supporting findings",
            )?;
            for id in &t.required_premise_ids {
                require(premises.contains(id), format!("unknown premise {id}"))?;
            }
            for id in &t.depends_on {
                require(
                    ts.contains_key(id.as_str()),
                    format!("unknown prerequisite task {id}"),
                )?;
            }
            for id in &t.resolution_ids {
                require(
                    fs.contains_key(id.as_str()),
                    format!("unknown resolution {id}"),
                )?;
            }
            for id in &t.resolution_ids {
                for premise in &fs[id.as_str()].premise_ids {
                    require(
                        t.required_premise_ids.contains(premise)
                            || q.accepted_premise_ids.contains(premise),
                        format!("task {} omits resolution premise {premise}", t.id),
                    )?;
                }
            }
            if t.status == "excluded" {
                require(
                    t.resolution_ids
                        .iter()
                        .any(|id| fs[id.as_str()].status == "excluded_scoped"),
                    "excluded task needs a scoped exclusion",
                )?;
            }
            if t.status == "resolved" && t.kind == "construction" {
                require(
                    t.resolution_ids
                        .iter()
                        .all(|id| fs[id.as_str()].status != "excluded_scoped"),
                    "an exclusion cannot resolve a construction task",
                )?;
            }
        }
        fn visit<'a>(
            id: &'a str,
            ts: &HashMap<&'a str, &'a Task>,
            active: &mut HashSet<&'a str>,
            done: &mut HashSet<&'a str>,
        ) -> Result<()> {
            if done.contains(id) {
                return Ok(());
            }
            require(active.insert(id), format!("cyclic task dependency: {id}"))?;
            for d in &ts[id].depends_on {
                visit(d, ts, active, done)?;
            }
            active.remove(id);
            done.insert(id);
            Ok(())
        }
        let mut done = HashSet::new();
        for t in &q.tasks {
            visit(&t.id, &ts, &mut HashSet::new(), &mut done)?;
        }
        for c in &q.connections {
            register(&c.id, "connection:")?;
            require(c.explanation.valid(), "bilingual connection required")?;
            require(
                p.problems
                    .iter()
                    .any(|other| other.id == c.target_problem_id && other.id != q.id),
                "unknown or self research connection target",
            )?;
            require(
                [
                    "requires_result",
                    "shares_unresolved_input",
                    "constrains_under",
                ]
                .contains(&c.kind.as_str()),
                "invalid connection kind",
            )?;
            require(
                !c.finding_ids.is_empty(),
                "connection needs supporting source-reviewed findings",
            )?;
            for id in &c.premise_ids {
                require(premises.contains(id), format!("unknown premise {id}"))?;
            }
            for id in &c.finding_ids {
                require(
                    fs.contains_key(id.as_str()),
                    format!("unknown connection finding {id}"),
                )?;
                for premise in &fs[id.as_str()].premise_ids {
                    require(
                        c.premise_ids.contains(premise) || q.accepted_premise_ids.contains(premise),
                        format!("connection {} omits finding premise {premise}", c.id),
                    )?;
                }
            }
        }
        for s in &q.strategies {
            register(&s.id, "strategy:")?;
            require(
                s.title.valid() && s.explanation.valid() && s.changes_required.valid(),
                "bilingual strategy required",
            )?;
            require(
                ["reuse_result", "excluded_scoped", "insufficient_for_goal"]
                    .contains(&s.outcome.as_str()),
                "invalid strategy outcome",
            )?;
            require(
                !s.finding_ids.is_empty(),
                "strategy decision needs findings",
            )?;
            for id in &s.finding_ids {
                require(
                    fs.contains_key(id.as_str()),
                    format!("unknown strategy finding {id}"),
                )?;
            }
            for id in &s.premise_ids {
                require(premises.contains(id), format!("unknown premise {id}"))?;
            }
            // Never allow a strategy to silently drop a premise of the result it reuses.
            for id in &s.finding_ids {
                for premise in &fs[id.as_str()].premise_ids {
                    require(
                        s.premise_ids.contains(premise) || q.accepted_premise_ids.contains(premise),
                        format!("strategy {} omits finding premise {premise}", s.id),
                    )?;
                }
            }
            if s.outcome == "excluded_scoped" {
                require(
                    s.finding_ids
                        .iter()
                        .any(|id| fs[id.as_str()].status == "excluded_scoped"),
                    "excluded strategy needs an exclusion",
                )?;
            }
        }
    }
    Ok(())
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn source_id(s: &Citation) -> String {
    format!(
        "research-source:{}",
        &digest(format!("{}:{}:{}:{}", s.path, s.sha256, s.line_start, s.line_end).as_bytes())
            [..24]
    )
}
fn source_current(
    s: &Citation,
    root: &Path,
    cache: &mut HashMap<String, Option<(String, usize)>>,
) -> bool {
    let entry = cache.entry(s.path.clone()).or_insert_with(|| {
        let canonical_root = fs::canonicalize(root).ok()?;
        let path = fs::canonicalize(canonical_root.join(&s.path)).ok()?;
        if !path.starts_with(&canonical_root) {
            return None;
        }
        let bytes = fs::read(path).ok()?;
        let lines = std::str::from_utf8(&bytes).ok()?.lines().count();
        Some((digest(&bytes), lines))
    });
    entry
        .as_ref()
        .is_some_and(|(h, n)| h == &s.sha256 && *n >= s.line_end)
}
fn freshness(p: &Program, root: &Path) -> HashMap<String, bool> {
    let mut cache = HashMap::new();
    p.problems
        .iter()
        .flat_map(|q| &q.findings)
        .map(|f| {
            (
                f.id.clone(),
                f.sources
                    .iter()
                    .all(|s| source_current(s, root, &mut cache)),
            )
        })
        .collect()
}
fn text_attrs(title: &Text, statement: &Text) -> Value {
    json!({"label_de":title.de,"label_en":title.en,"summary_de":statement.de,"summary_en":statement.en,"semantic":{"overview_de":statement.de,"overview_en":statement.en}})
}
fn node(id: &str, kind: &str, title: &Text, statement: &Text, mut attrs: Value) -> Node {
    for (k, v) in text_attrs(title, statement).as_object().unwrap() {
        attrs[k] = v.clone();
    }
    attrs["provenance"] = json!("source_reviewed_research_record_not_ledger_proof");
    Node {
        id: id.into(),
        kind: kind.into(),
        label: title.de.clone(),
        summary: statement.de.clone(),
        attrs,
    }
}
fn edge(g: &mut KnowledgeGraph, src: &str, dst: &str, kind: &str) {
    g.edges.push(Edge {
        id: format!(
            "research-edge:{}",
            &digest(format!("{src}:{kind}:{dst}").as_bytes())[..24]
        ),
        src: src.into(),
        dst: dst.into(),
        kind: kind.into(),
        attrs: json!({"origin":"research_program","proof_inference":false}),
    });
}
const PLAN_DIRECTORY: &str = "docs/research/fundamental-plan";

// The execution plan is live operational state, separate from historical
// scientific decisions. Never turn an old open task into a new work directive.
pub(crate) fn current_plan_sources(root: &Path) -> Vec<Source> {
    let Ok(canonical_root) = fs::canonicalize(root) else {
        return Vec::new();
    };
    [
        "progress.json",
        "PLAN.md",
        "RESULTS.md",
        "SEARCH-METHODS.md",
    ]
    .into_iter()
    .filter_map(|name| {
        let relative = format!("{PLAN_DIRECTORY}/{name}");
        let path = fs::canonicalize(root.join(&relative)).ok()?;
        if !path.starts_with(&canonical_root) || crate::service::sensitive_path(&path) {
            return None;
        }
        if fs::metadata(&path).ok()?.len() > 8 * 1024 * 1024 {
            return None;
        }
        let bytes = fs::read(path).ok()?;
        let contents = std::str::from_utf8(&bytes).ok()?;
        Some(Source {
            id: format!("source:fundamental-plan:{name}"),
            path: relative,
            sha256: digest(&bytes),
            line_count: Some(contents.lines().count()),
            ..Default::default()
        })
    })
    .collect()
}

fn current_plan(root: &Path, atlas_revision: &str, max_chars: usize) -> Value {
    let sources = current_plan_sources(root);
    let Some(source) = sources.iter().find(|s| s.path.ends_with("/progress.json")) else {
        return json!({"state":"unavailable","path":format!("{PLAN_DIRECTORY}/progress.json"),
            "execution_priority":"unknown","notice":"No readable current plan; historical task readiness is dependency status only."});
    };
    let loaded = fs::read(root.join(&source.path))
        .ok()
        .filter(|bytes| digest(bytes) == source.sha256)
        .and_then(|bytes| serde_json::from_slice::<Value>(&bytes).ok());
    let Some(plan) = loaded.filter(|plan| {
        plan.is_object()
            && plan["tasks"].is_array()
            && plan["next_task"].is_string()
            && plan["revision"].is_number()
    }) else {
        return json!({"state":"invalid","source_id":source.id,"execution_priority":"unknown",
            "notice":"Current plan could not be validated; do not infer execution priority from historical readiness."});
    };
    let next = plan["tasks"]
        .as_array()
        .unwrap()
        .iter()
        .find(|task| task["id"] == plan["next_task"]);
    let recent: Vec<_> = plan["result_ids"]
        .as_array()
        .into_iter()
        .flatten()
        .rev()
        .take(3)
        .cloned()
        .collect();
    let mut out = json!({"state":"current","revision":plan["revision"],"updated_at":plan["updated_at"],
        "source_id":source.id,"sha256":source.sha256,"next_task":plan["next_task"],
        "task_status":next.map(|task| &task["status"]),"task_title":next.map(|task| &task["title"]),
        "continuation":plan["continuation"],"recent_result_ids":recent,
        "atlas_revision_at_update":plan["current_atlas_revision"],
        "atlas_revision_matches":plan["current_atlas_revision"].as_str()==Some(atlas_revision),
        "sources":sources.iter().map(|s|json!({"id":s.id,"path":s.path,"sha256":s.sha256,"line_count":s.line_count})).collect::<Vec<_>>(),
        "content_complete":true,
        "notice":"Live execution priority. Historical task status and scoped scientific verdicts are preserved independently."});
    if serde_json::to_string(&out).unwrap().chars().count() > max_chars {
        out.as_object_mut().unwrap().remove("sources");
        out.as_object_mut().unwrap().remove("task_title");
        out.as_object_mut().unwrap().remove("notice");
        out.as_object_mut().unwrap().remove("sha256");
        out.as_object_mut()
            .unwrap()
            .remove("atlas_revision_at_update");
        out["content_complete"] = json!(false);
        out["next_request"] = json!({"tool":"read_source","arguments":{"source_id":source.id}});
        let continuation = plan["continuation"].as_str().unwrap_or_default();
        out["continuation"] = Value::Null;
        out["continuation_complete"] = json!(false);
        let overhead = serde_json::to_string(&out).unwrap().chars().count() + 16;
        let available = max_chars.saturating_sub(overhead);
        out["continuation"] = json!(continuation.chars().take(available).collect::<String>());
        out["continuation_complete"] = json!(continuation.chars().count() <= available);
        loop {
            let excess = serde_json::to_string(&out)
                .unwrap()
                .chars()
                .count()
                .saturating_sub(max_chars);
            if excess == 0 {
                break;
            }
            let text = out["continuation"].as_str().unwrap_or_default();
            if text.is_empty() {
                break;
            }
            out["continuation"] = json!(text
                .chars()
                .take(text.chars().count().saturating_sub(excess))
                .collect::<String>());
            out["continuation_complete"] = json!(false);
        }
    }
    out
}

fn execution_priority(plan: &Value) -> &'static str {
    if plan["state"] == "current" {
        "historical_consult_current_plan"
    } else {
        "unassigned"
    }
}

// Result identifiers are anchors, not ordinary vocabulary. A finding that owns
// R123 must precede a large historical problem that merely mentions its words.
fn result_identifiers(query: &str) -> HashSet<String> {
    query
        .to_lowercase()
        .split(|c: char| !c.is_alphanumeric())
        .filter(|token| {
            token.starts_with('r')
                && token.len() > 1
                && token[1..].chars().all(|c| c.is_ascii_digit())
        })
        .map(str::to_owned)
        .collect()
}

fn owns_result(id: &str, title: &str, requested: &HashSet<String>) -> bool {
    let headline = format!("{id} {title}").to_lowercase();
    !requested.is_empty()
        && headline
            .split(|c: char| !c.is_alphanumeric())
            .any(|token| requested.contains(token))
}

fn asks_current_priority(query: &str) -> bool {
    if query.contains(':') && !query.contains(' ') {
        return false;
    }
    let lower = query.to_lowercase();
    let words: Vec<_> = lower.split(|c: char| !c.is_alphanumeric()).collect();
    words.iter().any(|word| {
        word.starts_with("priorit")
            || word.starts_with("priorität")
            || word.starts_with("nächste")
            || *word == "next"
            || word.starts_with("fortsetz")
    }) || (words
        .iter()
        .any(|word| word.starts_with("aktuell") || *word == "current" || *word == "latest")
        && words.iter().any(|word| {
            word.starts_with("fundament")
                || word.starts_with("forschung")
                || ["research", "theory", "theorie", "stand"].contains(word)
        }))
}

pub fn augment(g: &mut KnowledgeGraph, p: &Program, root: &Path) {
    let current = freshness(p, root);
    let plan = current_plan(root, &p.revision, 6000);
    for source in current_plan_sources(root) {
        if !g.sources.iter().any(|existing| existing.id == source.id) {
            g.sources.push(source);
        }
    }
    let existing: HashSet<_> = g.nodes.iter().map(|n| n.id.clone()).collect();
    let mut pending_total = HashSet::new();
    let stale = current.values().filter(|v| !**v).count();
    for premise in &p.premises {
        g.nodes.push(node(
            &premise.id,
            "research_premise",
            &premise.title,
            &premise.statement,
            json!({"accepted_globally":false}),
        ));
    }
    for q in &p.problems {
        let pending = pending_attempts(q, g);
        pending_total.extend(pending.iter().map(|n| n.id.clone()));
        g.nodes.push(node(&q.id,"research_problem",&q.title,&q.question,json!({"revision":p.revision,"status":"research_question","pending_attempts":pending.len()})));
        for topic in q.topics.iter().chain(&q.goals) {
            if existing.contains(topic) {
                edge(g, &q.id, topic, "addresses");
            }
        }
        for premise in &q.accepted_premise_ids {
            edge(g, &q.id, premise, "accepts_premise");
        }
        for f in &q.findings {
            let refs: Vec<_> = f.sources.iter().map(source_id).collect();
            g.nodes.push(node(&f.id,"research_finding",&f.title,&f.statement,json!({"status":f.status,"scope_de":f.scope.de,"scope_en":f.scope.en,"source_state":if current[&f.id]{"current"}else{"stale"},"source_refs":refs,"premise_ids":f.premise_ids})));
            edge(g, &q.id, &f.id, "has_finding");
            for premise in &f.premise_ids {
                edge(g, &f.id, premise, "requires_premise");
            }
            for attempt in &f.attempt_ids {
                if existing.contains(attempt) {
                    edge(g, &f.id, attempt, "evidenced_by");
                }
            }
            for s in &f.sources {
                let id = source_id(s);
                if !g.sources.iter().any(|x| x.id == id) {
                    g.sources.push(Source {
                        id,
                        path: s.path.clone(),
                        sha256: s.sha256.clone(),
                        line_start: Some(s.line_start),
                        line_end: Some(s.line_end),
                        ..Default::default()
                    });
                }
            }
        }
        for t in &q.tasks {
            g.nodes.push(node(&t.id,"research_task",&t.title,&t.success_criterion,json!({"status":t.status,"effective_status":task_status(t,q,&current),"execution_priority":execution_priority(&plan),"task_kind":t.kind,"required_premise_ids":t.required_premise_ids,"next_action_de":t.next_action.de,"next_action_en":t.next_action.en,"depends_on":t.depends_on,"resolution_ids":t.resolution_ids})));
            edge(g, &q.id, &t.id, "has_task");
            for d in &t.depends_on {
                edge(g, &t.id, d, "requires_task");
            }
            for f in &t.resolution_ids {
                edge(g, f, &t.id, "resolves_in_scope");
            }
            for premise in &t.required_premise_ids {
                edge(g, &t.id, premise, "requires_premise");
            }
        }
        for c in &q.connections {
            g.edges.push(Edge{id:c.id.clone(),src:q.id.clone(),dst:c.target_problem_id.clone(),kind:c.kind.clone(),attrs:json!({"origin":"research_program","explanation_de":c.explanation.de,"explanation_en":c.explanation.en,"premise_ids":c.premise_ids,"finding_ids":c.finding_ids,"proof_inference":false})});
        }
        for s in &q.strategies {
            g.nodes.push(node(&s.id,"research_strategy",&s.title,&s.explanation,json!({"outcome":s.outcome,"premise_ids":s.premise_ids,"finding_ids":s.finding_ids})));
            edge(g, &q.id, &s.id, "records_strategy");
            for f in &s.finding_ids {
                edge(g, &s.id, f, "evaluated_by");
            }
        }
    }
    g.meta["research_program"] = json!({"revision":p.revision,"reviewed_problems":p.problems.len(),"reviewed_findings":current.len(),"stale_findings":stale,"pending_attempts":pending_total.len(),"coverage":"explicit source review only; pending records are not new scientific open questions","entry_tool":"research_state","plan_tool":"assess_research_plan","update_tool":"record_research_update"});
}
fn task_status<'a>(t: &'a Task, q: &'a Problem, current: &HashMap<String, bool>) -> &'a str {
    if t.resolution_ids
        .iter()
        .any(|id| !current.get(id).copied().unwrap_or(false))
    {
        return "needs_review";
    }
    if t.depends_on.iter().any(|id| {
        q.tasks
            .iter()
            .find(|x| &x.id == id)
            .is_none_or(|x| task_status(x, q, current) != "resolved")
    }) {
        return if t.status == "open" {
            "blocked"
        } else {
            "needs_review"
        };
    }
    match t.status.as_str() {
        "open" => "ready",
        other => other,
    }
}
fn pending_attempts<'a>(q: &Problem, g: &'a KnowledgeGraph) -> Vec<&'a Node> {
    let roots: HashSet<_> = q
        .topics
        .iter()
        .chain(&q.goals)
        .map(String::as_str)
        .collect();
    let linked: HashSet<_> = g
        .edges
        .iter()
        .filter_map(|e| {
            if roots.contains(e.dst.as_str()) {
                Some(e.src.as_str())
            } else if roots.contains(e.src.as_str()) {
                Some(e.dst.as_str())
            } else {
                None
            }
        })
        .collect();
    let reviewed: HashSet<_> = q
        .reviewed_attempt_ids
        .iter()
        .chain(q.findings.iter().flat_map(|f| &f.attempt_ids))
        .map(String::as_str)
        .collect();
    let mut out: Vec<_> = g
        .nodes
        .iter()
        .filter(|n| {
            n.kind == "attempt"
                && linked.contains(n.id.as_str())
                && !reviewed.contains(n.id.as_str())
        })
        .collect();
    out.sort_by(|a, b| b.id.cmp(&a.id));
    out
}
fn tokens(s: &str) -> HashSet<String> {
    let lower = s.to_lowercase().replace('α', " alpha ");
    let mut out: HashSet<String> = lower
        .split(|c: char| !c.is_alphanumeric())
        .filter(|w| {
            (w.chars().count() > 2
                || (w.chars().count() >= 2 && w.chars().any(|c| c.is_ascii_digit())))
                && ![
                    "the",
                    "and",
                    "der",
                    "die",
                    "das",
                    "was",
                    "wie",
                    "von",
                    "aus",
                    "für",
                    "mit",
                    "und",
                    "ist",
                    "schon",
                    "what",
                    "how",
                    "which",
                    "welche",
                    "welcher",
                    "welchen",
                    "wieso",
                    "warum",
                    "that",
                    "this",
                    "dass",
                    "eine",
                    "einer",
                    "einen",
                    "einem",
                    "den",
                    "dem",
                    "des",
                    "werden",
                    "wird",
                    "sind",
                    "sich",
                    "unter",
                    "oder",
                    "auch",
                    "noch",
                    "bereits",
                    "darüber",
                    "bekannt",
                    "muss",
                    "nächster",
                    "forschungsschritt",
                    "tatsächlich",
                    "gezeigt",
                    "herzuleiten",
                    "über",
                    "about",
                    "known",
                    "need",
                    "needs",
                    "next",
                    "from",
                    "given",
                    "does",
                    "have",
                    "been",
                ]
                .contains(w)
        })
        .map(|w| {
            if w == "alpha" || w.starts_with("feinstruktur") || w.starts_with("fine_structure") {
                "alpha".into()
            } else if w.starts_with("coefficient") || w.starts_with("koeffiz") {
                "coefficient".into()
            } else if w.starts_with("factor") || w.starts_with("faktor") || w.starts_with("korrekturfaktor") {
                "factor".into()
            } else if w.starts_with("gedächtn") || w.starts_with("memory") || w.starts_with("memories") {
                "memory".into()
            } else if w.starts_with("zukunft")
                || w.starts_with("zukünft")
                || w.starts_with("future")
            {
                "future".into()
            } else if w.starts_with("zustand")
                || w.starts_with("zuständ")
                || w == "state"
                || w == "states"
            {
                "state".into()
            } else if w.starts_with("zeit")
                || w == "time"
                || w == "dynamics"
                || w.starts_with("dynamik")
            {
                "time".into()
            } else if w.contains("projektor") || w.contains("projector") {
                "projector".into()
            } else if w.starts_with("markier") || w.starts_with("marking") {
                "marking".into()
            } else if w.starts_with("famil") {
                "family".into()
            } else if w.starts_with("naht") || w == "seam" {
                "seam".into()
            } else if w.starts_with("urspr") || w == "origin" || w == "original" {
                "origin".into()
            } else if w.starts_with("quelle") || w == "source" || w == "sources" {
                "source".into()
            } else if w.starts_with("positiv") {
                "positive".into()
            } else if w.starts_with("monodrom") {
                "monodromy".into()
            } else if ["matrix", "matrizen", "matrices"].contains(&w) {
                "matrix".into()
            } else if ["identität", "identitaet", "identity"].contains(&w) {
                "identity".into()
            } else {
                w.to_owned()
            }
        })
        .collect();
    let compact = lower.replace(' ', "");
    if compact.contains("2+2") || compact.contains("2plus2") {
        out.insert("two_plus_two".into());
    }
    if compact.contains("3+1") {
        out.insert("three_plus_one".into());
    }
    out
}
fn match_tokens(q: &Problem) -> (HashSet<String>, HashSet<String>) {
    let title = tokens(&format!(
        "{} {} {} {}",
        q.title.de, q.title.en, q.question.de, q.question.en
    ));
    let mut body = title.clone();
    for f in &q.findings {
        body.extend(tokens(&format!(
            "{} {} {} {}",
            f.title.de, f.title.en, f.statement.de, f.statement.en
        )));
    }
    for s in &q.strategies {
        body.extend(tokens(&format!(
            "{} {} {}",
            s.title.de,
            s.title.en,
            s.aliases.join(" ")
        )));
    }
    (title, body)
}
fn matches_problem(q: &Problem, query: &str) -> usize {
    if query.is_empty() {
        return 1;
    }
    if q.id == query
        || q.accepted_premise_ids
            .iter()
            .chain(q.findings.iter().flat_map(|f| &f.premise_ids))
            .any(|id| id == query)
        || q.topics.iter().chain(&q.goals).any(|s| s == query)
        || q.findings
            .iter()
            .any(|f| f.id == query || f.attempt_ids.iter().any(|id| id == query))
        || q.tasks.iter().any(|t| t.id == query)
        || q.strategies.iter().any(|s| s.id == query)
    {
        return 10000;
    }
    if query.contains(':') && !query.contains(' ') {
        return 0;
    }
    let (title, body) = match_tokens(q);
    let needle = tokens(query);
    needle.intersection(&body).count() + 3 * needle.intersection(&title).count()
}
fn finding_view(f: &Finding, lang: &str, current: &HashMap<String, bool>) -> Value {
    json!({"id":f.id,"title":f.title.get(lang),"statement":f.statement.get(lang),"status":f.status,"scope":f.scope.get(lang),"premise_ids":f.premise_ids,"source_state":if current[&f.id]{"current"}else{"stale"},"sources":f.sources.iter().map(|s|json!({"source_id":source_id(s),"path":s.path,"line_start":s.line_start,"line_end":s.line_end,"sha256":s.sha256})).collect::<Vec<_>>(),"attempt_ids":f.attempt_ids})
}
fn premise_view(p: &Premise, lang: &str) -> Value {
    json!({"id":p.id,"title":p.title.get(lang),"statement":p.statement.get(lang)})
}
fn problem_premise_ids(q: &Problem) -> HashSet<&String> {
    q.accepted_premise_ids
        .iter()
        .chain(q.findings.iter().flat_map(|f| &f.premise_ids))
        .chain(q.tasks.iter().flat_map(|t| &t.required_premise_ids))
        .chain(q.strategies.iter().flat_map(|s| &s.premise_ids))
        .chain(q.connections.iter().flat_map(|c| &c.premise_ids))
        .collect()
}

fn problem_view(
    q: &Problem,
    p: &Program,
    g: &KnowledgeGraph,
    lang: &str,
    current: &HashMap<String, bool>,
) -> Value {
    let used = problem_premise_ids(q);
    let pending = pending_attempts(q, g);
    json!({"id":q.id,"title":q.title.get(lang),"question":q.question.get(lang),"topics":q.topics,"goals":q.goals,
        "accepted_premises":p.premises.iter().filter(|x|q.accepted_premise_ids.contains(&x.id)).map(|x|premise_view(x,lang)).collect::<Vec<_>>(),
        "premises":p.premises.iter().filter(|x|used.contains(&x.id)).map(|x|premise_view(x,lang)).collect::<Vec<_>>(),
        "findings":q.findings.iter().map(|f|finding_view(f,lang,current)).collect::<Vec<_>>(),
        "tasks":q.tasks.iter().map(|t|json!({"id":t.id,"title":t.title.get(lang),"kind":t.kind,"status":t.status,"effective_status":task_status(t,q,current),"depends_on":t.depends_on,"resolution_ids":t.resolution_ids,"required_premise_ids":t.required_premise_ids,"success_criterion":t.success_criterion.get(lang),"next_action":t.next_action.get(lang)})).collect::<Vec<_>>(),
        "connections":p.problems.iter().flat_map(|origin|origin.connections.iter().filter(|c|origin.id==q.id || c.target_problem_id==q.id).map(move |c|(origin,c))).map(|(origin,c)|{let target=if origin.id==q.id {&c.target_problem_id}else{&origin.id};json!({"id":c.id,"source_problem_id":origin.id,"target_problem_id":c.target_problem_id,"other_problem_id":target,"other_title":p.problems.iter().find(|other|&other.id==target).map(|other|other.title.get(lang)),"direction":if origin.id==q.id {"outgoing"}else{"incoming"},"target_title":p.problems.iter().find(|other|&other.id==target).map(|other|other.title.get(lang)),"kind":c.kind,"premise_ids":c.premise_ids,"finding_ids":c.finding_ids,"explanation":c.explanation.get(lang),"source_state":if c.finding_ids.iter().all(|id|current.get(id).copied().unwrap_or(false)){"current"}else{"stale"},"proof_inference":false})}).collect::<Vec<_>>(),
        "strategies":q.strategies.iter().map(|s|json!({"id":s.id,"title":s.title.get(lang),"outcome":s.outcome,"premise_ids":s.premise_ids,"finding_ids":s.finding_ids,"explanation":s.explanation.get(lang),"changes_required":s.changes_required.get(lang)})).collect::<Vec<_>>(),
        "pending_attempts":{"total":pending.len(),"items":pending.iter().take(3).map(|n|json!({"id":n.id,"title":n.attrs.get(if lang=="en"{"label_en"}else{"label_de"}).and_then(Value::as_str).unwrap_or(&n.label)})).collect::<Vec<_>>()}})
}

const DETAIL_FIELDS: [&str; 4] = ["findings", "tasks", "strategies", "connections"];

fn detail_tokens(s: &str) -> HashSet<String> {
    let mut out = tokens(s);
    out.extend(
        s.to_lowercase()
            .split(|c: char| !c.is_alphanumeric())
            .filter(|word| word.chars().count() >= 2 && word.chars().any(|c| c.is_ascii_digit()))
            .map(str::to_owned),
    );
    out
}

fn detail_records(full: &Value, query: &str) -> Result<Vec<(&'static str, usize)>> {
    let needle = detail_tokens(query);
    let requested_results = result_identifiers(query);
    let mut documents = Vec::new();
    for (field_order, field) in DETAIL_FIELDS.iter().enumerate() {
        for (index, item) in full[*field]
            .as_array()
            .context("problem detail field must be an array")?
            .iter()
            .enumerate()
        {
            let id = item["id"].as_str().unwrap_or_default();
            let headline = detail_tokens(&format!(
                "{} {}",
                id,
                item["title"].as_str().unwrap_or_default()
            ));
            let mut content = headline.clone();
            for key in [
                "statement",
                "scope",
                "explanation",
                "changes_required",
                "success_criterion",
                "next_action",
            ] {
                content.extend(detail_tokens(&serde_json::to_string(&item[key])?));
            }
            let metadata = detail_tokens(&serde_json::to_string(item)?);
            documents.push((field_order, index, *field, id, headline, content, metadata));
        }
    }
    // Sources and dependency IDs are useful fallbacks, but must not outrank
    // the result itself. Rare content terms distinguish a precise finding
    // inside a large problem; metadata-only terms receive no rarity bonus.
    let frequencies: HashMap<_, _> = needle
        .iter()
        .map(|token| {
            (
                token,
                documents
                    .iter()
                    .filter(|doc| doc.5.contains(token))
                    .count()
                    .max(1),
            )
        })
        .collect();
    let mut records = Vec::new();
    for (field_order, index, field, id, headline, content, metadata) in &documents {
        let score = if *id == query {
            usize::MAX
        } else if needle.is_empty() || query == full["id"].as_str().unwrap_or_default() {
            0
        } else {
            needle
                .iter()
                .map(|token| {
                    let rarity = (documents.len() / frequencies[token]).max(1);
                    if headline.contains(token) {
                        20 * rarity
                    } else if content.contains(token) {
                        5 * rarity
                    } else if metadata.contains(token) {
                        1
                    } else {
                        0
                    }
                })
                .sum()
        };
        let coverage = if needle.is_empty() || query == full["id"].as_str().unwrap_or_default() {
            0
        } else {
            needle.intersection(content).count()
        };
        let anchored = if *id == query {
            2
        } else if owns_result(
            id,
            full[*field][*index]["title"].as_str().unwrap_or_default(),
            &requested_results,
        ) {
            1
        } else {
            0
        };
        records.push((anchored, coverage, score, *field_order, *index, *field));
    }
    records.sort_by(|a, b| {
        b.0.cmp(&a.0)
            .then(b.1.cmp(&a.1))
            .then(b.2.cmp(&a.2))
            .then(a.3.cmp(&b.3))
            .then(a.4.cmp(&b.4))
    });
    Ok(records
        .into_iter()
        .map(|(_, _, _, _, index, field)| (field, index))
        .collect())
}
fn referenced_ids(item: &Value, keys: &[&str]) -> Vec<String> {
    keys.iter()
        .flat_map(|key| {
            item[*key]
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(Value::as_str)
                .map(str::to_owned)
        })
        .collect()
}

fn record_handle(full: &Value, id: &str, included: &HashSet<String>) -> Value {
    for (kind, field) in [
        ("finding", "findings"),
        ("task", "tasks"),
        ("strategy", "strategies"),
        ("connection", "connections"),
    ] {
        if let Some(item) = full[field]
            .as_array()
            .into_iter()
            .flatten()
            .find(|item| item["id"] == id)
        {
            return json!({
                "id": id,
                "kind": kind,
                "title": item.get("title"),
                "status": item.get("effective_status").or_else(||item.get("status")),
                "source_state": item.get("source_state"),
                "included": included.contains(id),
            });
        }
    }
    json!({"id":id,"kind":"external","included":false})
}

fn partial_problem_view(
    full: &Value,
    p: &Program,
    query: &str,
    lang: &str,
    detail_offset: usize,
    take: usize,
) -> Result<Value> {
    let records = detail_records(full, query)?;
    let start = detail_offset.min(records.len());
    let end = (start + take).min(records.len());
    let selected = &records[start..end];
    let mut out = json!({
        "id": full["id"],
        "title": full["title"],
        "question": full["question"],
        "topics": full["topics"],
        "goals": full["goals"],
        "execution_priority": full["execution_priority"],
        "packet_state": "partial",
        "accepted_premise_ids": [],
        "accepted_premises": [],
        "premises": [],
        "premise_handles": [],
        "dependency_handles": [],
        "source_handles": [],
        "active_task_handles": [],
        "findings": [],
        "tasks": [],
        "strategies": [],
        "connections": [],
        "pending_attempts": full["pending_attempts"],
        "notice": if lang=="en" {
            "Partial problem page. Every returned record retains its complete scoped text, premise IDs and source handles; premise entries are compact ID/title handles. Unreturned records may contain additional conclusions or prerequisites; this page is not a proof of the whole problem."
        } else {
            "Partielle Problemseite. Jeder ausgegebene Datensatz behält seinen vollständigen Geltungsbereich, seine Premissen-IDs und Quellenverweise; Premisseneinträge sind kompakte ID- und Titelverweise. Nicht ausgegebene Datensätze können weitere Schlussfolgerungen oder Voraussetzungen enthalten; diese Seite ist kein Beweis des Gesamtproblems."
        },
    });
    let mut included = HashSet::new();
    let mut premise_ids = HashSet::new();
    let mut dependency_ids = HashSet::new();
    let mut source_handles = Vec::new();
    let mut returned = HashMap::<&str, usize>::new();
    for (field, index) in selected {
        let item = full[*field][*index].clone();
        if let Some(id) = item["id"].as_str() {
            included.insert(id.to_owned());
        }
        premise_ids.extend(referenced_ids(
            &item,
            &["premise_ids", "required_premise_ids"],
        ));
        dependency_ids.extend(referenced_ids(
            &item,
            &["depends_on", "resolution_ids", "finding_ids"],
        ));
        if *field == "findings" {
            for source in item["sources"].as_array().into_iter().flatten() {
                let mut handle = source.clone();
                handle["finding_id"] = item["id"].clone();
                source_handles.push(handle);
            }
        }
        out[*field].as_array_mut().unwrap().push(item);
        *returned.entry(*field).or_default() += 1;
    }
    // A returned finding may be the resolution of a task not present on this page.
    for task in full["tasks"].as_array().into_iter().flatten() {
        if referenced_ids(task, &["resolution_ids"])
            .iter()
            .any(|id| included.contains(id))
        {
            if let Some(id) = task["id"].as_str() {
                dependency_ids.insert(id.to_owned());
            }
        }
    }
    let premise_handles: Vec<_> = p
        .premises
        .iter()
        .filter(|premise| premise_ids.contains(&premise.id))
        .map(|premise| json!({"id":premise.id,"title":premise.title.get(lang)}))
        .collect();
    let all_accepted: HashSet<_> = full["accepted_premises"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(|premise| premise["id"].as_str())
        .collect();
    let mut accepted: Vec<_> = premise_ids
        .iter()
        .filter(|id| all_accepted.contains(id.as_str()))
        .cloned()
        .collect();
    accepted.sort();
    let accepted_returned = accepted.len();
    let accepted_complete = accepted_returned == all_accepted.len();
    out["accepted_premise_ids"] = json!(accepted);
    out["accepted_premises"] = json!(premise_handles
        .iter()
        .filter(|premise| premise["id"]
            .as_str()
            .is_some_and(|id| all_accepted.contains(id)))
        .cloned()
        .collect::<Vec<_>>());
    out["premises"] = json!(premise_handles);
    let mut premise_ids_sorted: Vec<_> = premise_ids.iter().collect();
    premise_ids_sorted.sort();
    let mut dependency_ids_sorted: Vec<_> = dependency_ids.iter().collect();
    dependency_ids_sorted.sort();
    out["premise_handles"] = json!(premise_ids_sorted
        .into_iter()
        .map(|id| json!({"id":id,"included":true}))
        .collect::<Vec<_>>());
    out["dependency_handles"] = json!(dependency_ids_sorted
        .into_iter()
        .map(|id| record_handle(full, id, &included))
        .collect::<Vec<_>>());
    out["source_handles"] = json!(source_handles);
    let active_tasks: Vec<_> = records
        .iter()
        .filter(|(field, index)| *field == "tasks" && full["tasks"][*index]["status"] == "open")
        .map(|(_, index)| &full["tasks"][*index])
        .collect();
    out["active_task_handles"] = json!(active_tasks
        .iter()
        .take(3)
        .map(|task|json!({
            "id":task["id"],
            "title":task["title"],
            "status":task["status"],
            "effective_status":task["effective_status"],
            "execution_priority":task["execution_priority"],
            "required_premise_ids":task["required_premise_ids"],
            "included":task["id"].as_str().is_some_and(|id|included.contains(id)),
            "next_query":task["id"],
            "next_request":{"query":task["id"],"lang":lang,"max_chars":12_000,"offset":0,"detail_offset":0},
        }))
        .collect::<Vec<_>>());
    let total = records.len();
    out["partial_coverage"] = json!({
        "complete": end == total,
        "selection": "query_relevance_then_stable_record_order",
        "detail_offset": start,
        "detail_returned": end-start,
        "detail_total": total,
        "next_detail_offset": if end<total {json!(end)} else {Value::Null},
        "returned": {
            "findings": returned.get("findings").copied().unwrap_or(0),
            "tasks": returned.get("tasks").copied().unwrap_or(0),
            "strategies": returned.get("strategies").copied().unwrap_or(0),
            "connections": returned.get("connections").copied().unwrap_or(0),
        },
        "totals": {
            "findings": full["findings"].as_array().map_or(0, Vec::len),
            "tasks": full["tasks"].as_array().map_or(0, Vec::len),
            "strategies": full["strategies"].as_array().map_or(0, Vec::len),
            "connections": full["connections"].as_array().map_or(0, Vec::len),
        },
        "accepted_premise_ids_complete": accepted_complete,
        "accepted_premise_ids_returned": accepted_returned,
        "accepted_premise_ids_total": all_accepted.len(),
        "premise_ids_complete_for_returned_records": true,
        "premise_text_state": "handles_only",
        "returned_record_content_complete": true,
        "unreturned_record_content_complete": false,
        "active_task_handles_returned": active_tasks.len().min(3),
        "active_task_handles_total": active_tasks.len(),
        "active_task_handles_omitted": active_tasks.len().saturating_sub(3),
        "active_task_handle_content_complete": false,
    });
    Ok(out)
}

fn overview_state(mut value: Value, graph: &KnowledgeGraph, cap: usize, offset: usize, lang: &str) -> Result<Value> {
    value["view"] = json!("theory_overview");
    value["notice"] = json!(if lang == "en" {
        "Scoped results and their evidence are separate from the current work plan. Open a problem by ID for details."
    } else {
        "Ergebnisse mit Geltungsbereich und Quellen stehen getrennt vom Arbeitsplan. Einzelthemen über ihre genaue Problem-ID öffnen."
    });
    value.as_object_mut().unwrap().remove("topic_coverage");
    value["coverage"].as_object_mut().unwrap().remove("pending_attempts_scope");
    if let Some(total) = graph.meta.pointer("/coverage/scripts_total").and_then(Value::as_u64) {
        value["coverage"]["indexed_verification_scripts"] = json!(total);
    }
    let plan = &value["current_plan"];
    value["current_plan"] = json!({"next_task":plan["next_task"], "source":plan["next_request"]});
    let mut index = value["problem_index"].take();
    let index = index.as_array_mut().unwrap();
    let total = index.len();
    value["problem_index"] = json!([]);
    value["evidence_areas"] = crate::overview::areas(lang, false);
    value["budget"]["evidence_areas_compact"] = json!(false);
    if serde_json::to_string(&value)?.chars().count() + 400 > cap {
        value["evidence_areas"] = crate::overview::areas(lang, true);
        value["budget"]["evidence_areas_compact"] = json!(true);
    }
    value["budget"]["index_total"] = json!(total);
    let mut next = offset.min(total);
    for item in index.iter().skip(offset) {
        let cards = value["problem_index"].as_array_mut().unwrap();
        cards.push(json!({"id":item["id"], "title":item["title"]}));
        if serde_json::to_string(&value)?.chars().count() + 350 > cap {
            value["problem_index"].as_array_mut().unwrap().pop();
            break;
        }
        next += 1;
    }
    value["budget"]["next_offset"] = json!(next);
    value["budget"]["truncated"] = json!(next < total || value["budget"]["evidence_areas_compact"] == true);
    if next < total {
        value["budget"]["continuation"] = json!({"query":value["query"], "lang":lang,
            "max_chars":if next == offset { 12_000 } else { cap },
            "offset":next, "expected_revision":value["revision"]});
    }
    for _ in 0..3 {
        value["budget"]["serialized_chars"] = json!(serde_json::to_string(&value)?.chars().count());
    }
    require(serde_json::to_string(&value)?.chars().count() <= cap, "overview exceeds max_chars")?;
    Ok(value)
}

pub fn state(
    p: &Program,
    g: &KnowledgeGraph,
    root: &Path,
    query: &str,
    lang: &str,
    max_chars: usize,
    offset: usize,
    detail_offset: usize,
) -> Result<Value> {
    require(["de", "en"].contains(&lang), "lang must be de or en")?;
    require(
        (2000..=60000).contains(&max_chars),
        "max_chars must be between 2000 and 60000",
    )?;
    let requested_results = result_identifiers(query);
    let overview = crate::overview::requested(query) && requested_results.is_empty() && detail_offset == 0;
    let current = freshness(p, root);
    let needle = tokens(query);
    let plan = current_plan(
        root,
        &p.revision,
        if overview || (query.contains(':') && !query.contains(' ')) {
            650
        } else {
            (max_chars / 3).min(6000)
        },
    );
    let documents: Vec<_> = p.problems.iter().map(match_tokens).collect();
    let priority_results: Vec<_> = if !overview && requested_results.is_empty() && asks_current_priority(query) {
        plan["recent_result_ids"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .map(result_identifiers)
            .collect()
    } else {
        Vec::new()
    };
    let mut matches: Vec<_> = p
        .problems
        .iter()
        .zip(&documents)
        .map(|(q, (title, _body))| {
            let exact = matches_problem(q, query);
            let score = if overview {
                1
            } else if exact == 10000 {
                30_000
            } else if q.findings.iter().any(|finding| {
                owns_result(
                    &finding.id,
                    &format!("{} {}", finding.title.de, finding.title.en),
                    &requested_results,
                )
            }) {
                20_000
            } else if let Some(position) = priority_results.iter().position(|ids| {
                q.findings.iter().any(|finding| {
                    owns_result(
                        &finding.id,
                        &format!("{} {}", finding.title.de, finding.title.en),
                        ids,
                    )
                })
            }) {
                15_000usize.saturating_sub(position)
            } else if query.is_empty() {
                1 + 100
                    * p.problems
                        .iter()
                        .flat_map(|x| &x.connections)
                        .filter(|c| c.target_problem_id == q.id)
                        .count()
            } else if query.contains(':') && !query.contains(' ') {
                0
            } else {
                // A large historical problem must not combine unrelated findings
                // into a better match than a focused answer to the actual question.
                let score_passage = |headline: &HashSet<String>, content: &HashSet<String>| {
                    needle.iter().filter(|t| content.contains(*t)).map(|t| {
                        let frequency = documents.iter()
                            .filter(|(_, words)| words.contains(t)).count().max(1);
                        (if headline.contains(t) { 400 } else { 100 }) / frequency
                    }).sum::<usize>()
                };
                let problem_title = tokens(&format!("{} {}", q.title.de, q.title.en));
                let heading_score = score_passage(&problem_title, title) * 4;
                let finding_score = q.findings.iter().map(|finding| {
                    let headline = tokens(&format!("{} {}", finding.title.de, finding.title.en));
                    let mut content = headline.clone();
                    content.extend(tokens(&format!("{} {}", finding.statement.de, finding.statement.en)));
                    score_passage(&headline, &content)
                }).max().unwrap_or(0);
                let strategy_score = q.strategies.iter().map(|strategy| {
                    let headline = tokens(&format!("{} {} {}", strategy.title.de, strategy.title.en,
                        strategy.aliases.join(" ")));
                    score_passage(&headline, &headline)
                }).max().unwrap_or(0);
                let topic_score = q.topics.iter().map(|topic| {
                    let words = tokens(topic);
                    needle.intersection(&words).count()
                }).sum::<usize>() * 1000;
                topic_score + heading_score + 4 * finding_score.max(strategy_score)
            };
            (score, q)
        })
        .filter(|(n, _)| *n > 0)
        .collect();
    matches.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.id.cmp(&b.1.id)));
    let summary_cards:Vec<_>=matches.iter().map(|(_,q)|json!({"id":q.id,"title":q.title.get(lang),"findings":q.findings.len(),"resolved_tasks":q.tasks.iter().filter(|t|task_status(t,q,&current)=="resolved").count(),"execution_priority":execution_priority(&plan),"dependency_ready_tasks":q.tasks.iter().filter(|t|task_status(t,q,&current)=="ready").map(|t|json!({"id":t.id,"title":t.title.get(lang),"execution_priority":execution_priority(&plan)})).collect::<Vec<_>>()})).collect();
    let pending: HashSet<_> = matches
        .iter()
        .flat_map(|(_, q)| pending_attempts(q, g))
        .map(|n| &n.id)
        .collect();
    let mut result = json!({"revision":p.revision,"query":query,"lang":lang,"coverage":{"reviewed_problems":p.problems.len(),"reviewed_findings":current.len(),"stale_findings":current.values().filter(|v|!**v).count(),"matched_problems":matches.len(),"pending_attempts":pending.len(),"corpus_complete":false},"problem_index":summary_cards,"problems":[],"budget":{"max_chars":max_chars,"offset":offset,"detail_offset":detail_offset,"next_offset":offset,"truncated":false},"notice":if lang=="en"{"Source-reviewed decisions in stated scopes. New or unreviewed attempts are a curation backlog, not new scientific open questions. Graph membership and source hashes do not prove conclusions. Empty matches mean missing reviewed coverage, not absence of results."}else{"Quellengeprüfte Einordnungen mit Geltungsbereich. Noch nicht eingearbeitete Versuche sind ein Prüfbestand, keine neuen wissenschaftlichen offenen Fragen. Graphzuordnung und Quellenhash beweisen keine Schlussfolgerung. Keine Treffer bedeutet fehlende geprüfte Abdeckung, nicht fehlende Ergebnisse."}});
    result["current_plan"] = plan;
    result["topic_coverage"]=json!(g.nodes.iter().filter(|n|n.kind=="topic").map(|n|{
        let count=p.problems.iter().filter(|q|q.topics.contains(&n.id)).count();
        json!({"id":n.id,"title":n.attrs.get(if lang=="en"{"label_en"}else{"label_de"}).and_then(Value::as_str).unwrap_or(&n.label),"reviewed_problems":count,"reviewed":count>0})
    }).collect::<Vec<_>>());
    if !query.is_empty() {
        result.as_object_mut().unwrap().remove("topic_coverage");
    }
    result["coverage"]["pending_attempts_scope"]=json!("union of attempts directly linked to topics/goals of matched reviewed problems; not scientific open questions");
    if overview {
        return overview_state(result, g, max_chars, offset, lang);
    }
    if matches.is_empty() {
        let root_node = g.nodes.iter().find(|n| n.id == query);
        let linked: HashSet<_> = g
            .edges
            .iter()
            .filter_map(|e| {
                if e.src == query {
                    Some(e.dst.as_str())
                } else if e.dst == query {
                    Some(e.src.as_str())
                } else {
                    None
                }
            })
            .collect();
        let attempts = g
            .nodes
            .iter()
            .filter(|n| n.kind == "attempt" && linked.contains(n.id.as_str()))
            .count();
        let claims = g
            .nodes
            .iter()
            .filter(|n| n.kind == "claim" && linked.contains(n.id.as_str()))
            .count();
        result["corpus_fallback"] = json!({"reviewed_decision_state":"not_compiled","id":root_node.map(|n|&n.id),"directly_linked_attempts":if root_node.is_some(){json!(attempts)}else{Value::Null},"directly_linked_claims":if root_node.is_some(){json!(claims)}else{Value::Null},"next_tool":if root_node.is_some(){"research_context"}else{"search_knowledge"},"arguments":if root_node.is_some(){json!({"id":query,"max_chars":12000})}else{json!({"query":query,"limit":5})},"web_route":root_node.map(|n|format!("#topic/{}",n.id))});
        if root_node.is_some() {
            result["coverage"]["pending_attempts"] = json!(attempts);
            result["coverage"]["pending_attempts_scope"] =
                json!("directly linked attempts for this unreviewed query node");
        }
    }
    // Prefer complete problem packets. If the first match is intrinsically too large,
    // return a record-complete inner page with explicit coverage and continuation.
    let mut partial = None;
    for (_, q) in matches.iter().skip(offset) {
        let mut view = problem_view(q, p, g, lang, &current);
        view["execution_priority"] = json!(execution_priority(&result["current_plan"]));
        for task in view["tasks"].as_array_mut().unwrap() {
            task["execution_priority"] = json!(execution_priority(&result["current_plan"]));
        }
        let first = result["problems"].as_array().unwrap().is_empty();
        if !first || detail_offset == 0 {
            result["problems"]
                .as_array_mut()
                .unwrap()
                .push(view.clone());
        }
        let mut over_budget = serde_json::to_string(&result)?.chars().count() + 200 > max_chars;
        if over_budget && first {
            // The first complete decision has priority over redundant navigation.
            result.as_object_mut().unwrap().remove("topic_coverage");
            let compact_index: Vec<_> = result["problem_index"].as_array().unwrap().iter()
                .map(|x|json!({"id":x["id"],"title":x["title"],"findings":x["findings"],"resolved_tasks":x["resolved_tasks"]})).collect();
            result["problem_index"] = json!(compact_index);
            result["budget"]["index_compacted"] = json!(true);
            over_budget = serde_json::to_string(&result)?.chars().count() + 200 > max_chars;
        }
        if !first && over_budget {
            result["problems"].as_array_mut().unwrap().pop();
            break;
        }
        if first && detail_offset == 0 && !over_budget {
            continue;
        }
        if first && detail_offset == 0 {
            result["problems"].as_array_mut().unwrap().pop();
        }
        if result["problems"].as_array().unwrap().is_empty() {
            if result["problem_index"].as_array().unwrap().len() > 1 {
                result["problem_index"].as_array_mut().unwrap().truncate(1);
                result["budget"]["index_scoped_to_partial_problem"] = json!(true);
            }
            let full_chars = serde_json::to_string(&view)?.chars().count();
            let records = detail_records(&view, query)?;
            let start = detail_offset.min(records.len());
            let mut returned = 0;
            let mut best = partial_problem_view(&view, p, query, lang, start, 0)?;
            // Count continuation metadata while packing, rather than evicting
            // the sole problem-index card after the page has already been filled.
            result["budget"]["partial_problem"] = json!(true);
            result["budget"]["full_problem_chars"] = json!(full_chars);
            result["budget"]["next_detail_offset"] = json!(start);
            result["budget"]["next_query"] = json!(q.id);
            result["budget"]["continuation"] = json!({"query":query,"lang":lang,
                "max_chars":max_chars,"offset":offset,"detail_offset":start,
                "expected_revision":p.revision});
            result["problems"]
                .as_array_mut()
                .unwrap()
                .push(best.clone());
            while start + returned < records.len() {
                let candidate = partial_problem_view(&view, p, query, lang, start, returned + 1)?;
                result["problems"].as_array_mut().unwrap()[0] = candidate.clone();
                if serde_json::to_string(&result)?.chars().count() + 200 > max_chars {
                    result["problems"].as_array_mut().unwrap()[0] = best;
                    break;
                }
                best = candidate;
                returned += 1;
            }
            let end = start + returned;
            let detail_complete = end == records.len();
            result["budget"]["partial_problem"] = json!(true);
            result["budget"]["full_problem_chars"] = json!(full_chars);
            result["budget"]["next_detail_offset"] = if detail_complete {
                Value::Null
            } else {
                json!(end)
            };
            result["budget"]["next_query"] = json!(q.id);
            result["budget"]["continuation"] = if detail_complete {
                json!({"query":query,"lang":lang,"max_chars":max_chars,"offset":offset+1,"detail_offset":0,"expected_revision":p.revision})
            } else {
                json!({"query":query,"lang":lang,"max_chars":max_chars,"offset":offset,"detail_offset":end,"expected_revision":p.revision})
            };
            if !detail_complete {
                let one = partial_problem_view(&view, p, query, lang, end, 1)?;
                let mut required = result.clone();
                required["problems"].as_array_mut().unwrap()[0] = one;
                required["budget"]
                    .as_object_mut()
                    .unwrap()
                    .remove("required_chars_for_next_detail");
                let needed = serde_json::to_string(&required)?.chars().count() + 200;
                if needed > max_chars {
                    result["budget"]["required_chars_for_next_detail"] = json!(needed);
                    result["budget"]["continuation"]["max_chars"] = json!(needed.min(60_000));
                }
            }
            partial = Some(detail_complete);
            break;
        }
        break;
    }
    let count = result["problems"].as_array().unwrap().len();
    if let Some(detail_complete) = partial {
        result["budget"]["next_offset"] = json!(if detail_complete { offset + 1 } else { offset });
        result["budget"]["truncated"] = json!(!detail_complete || offset + 1 < matches.len());
    } else {
        result["budget"]["next_offset"] = json!(offset + count);
        result["budget"]["truncated"] = json!(offset + count < matches.len());
    }
    if count == 0 && offset < matches.len() && partial.is_none() {
        let needed = serde_json::to_string(&problem_view(matches[offset].1, p, g, lang, &current))?
            .chars()
            .count()
            + serde_json::to_string(&result)?.chars().count()
            + 300;
        result["budget"]["required_chars_for_next_problem"] = json!(needed);
        result["budget"]["next_query"] = json!(matches[offset].1.id);
    }
    // The index is navigation only; paginate it explicitly if a very small budget cannot hold it.
    while serde_json::to_string(&result)?.chars().count() + 100 > max_chars
        && !result["problem_index"].as_array().unwrap().is_empty()
    {
        result["problem_index"].as_array_mut().unwrap().pop();
    }
    result["budget"]["index_omitted"] =
        json!(matches.len() - result["problem_index"].as_array().unwrap().len());
    let topic_total = result
        .get("topic_coverage")
        .and_then(Value::as_array)
        .map_or(0, Vec::len);
    while serde_json::to_string(&result)?.chars().count() + 200 > max_chars
        && result
            .get("topic_coverage")
            .and_then(Value::as_array)
            .is_some_and(|v| !v.is_empty())
    {
        result["topic_coverage"].as_array_mut().unwrap().pop();
    }
    result["budget"]["topic_index_omitted"] = json!(
        topic_total
            - result
                .get("topic_coverage")
                .and_then(Value::as_array)
                .map_or(0, Vec::len)
    );
    for _ in 0..3 {
        result["budget"]["serialized_chars"] =
            json!(serde_json::to_string(&result)?.chars().count());
    }
    require(
        serde_json::to_string(&result)?.chars().count() <= max_chars,
        "internal response budget exceeded",
    )?;
    Ok(result)
}
pub fn assess(
    p: &Program,
    root: &Path,
    problem_id: &str,
    strategy: &str,
    supplied: &[String],
    lang: &str,
) -> Result<Value> {
    require(["de", "en"].contains(&lang), "lang must be de or en")?;
    let q = p
        .problems
        .iter()
        .find(|q| q.id == problem_id)
        .context("unknown research problem")?;
    let Some(s) = q.strategies.iter().find(|s| s.id == strategy) else {
        let needle = tokens(strategy);
        let candidates: Vec<_> = q
            .strategies
            .iter()
            .filter(|s| !tokens(&serde_json::to_string(s).unwrap_or_default()).is_disjoint(&needle))
            .map(|s| json!({"id":s.id,"title":s.title.get(lang),"premise_ids":s.premise_ids}))
            .collect();
        return Ok(
            json!({"decision":"unclassified","candidates":candidates,"explanation":if lang=="en"{"Text similarity is not applicability. Select an exact strategy ID and state its premises; no match is not evidence of novelty."}else{"Textähnlichkeit belegt keine Anwendbarkeit. Eine konkrete Strategie-ID und deren Voraussetzungen angeben; kein Treffer ist kein Neuheitsbeleg."}}),
        );
    };
    let known: HashSet<_> = p.premises.iter().map(|p| p.id.as_str()).collect();
    for id in supplied {
        require(
            known.contains(id.as_str()),
            format!(
                "unknown premise {id}; record its definition before evaluating a changed premise"
            ),
        )?;
    }
    let have: HashSet<_> = supplied.iter().chain(&q.accepted_premise_ids).collect();
    let missing: Vec<_> = s
        .premise_ids
        .iter()
        .filter(|id| !have.contains(id))
        .collect();
    let extra: Vec<_> = supplied
        .iter()
        .filter(|id| !s.premise_ids.contains(id) && !q.accepted_premise_ids.contains(id))
        .collect();
    let current = freshness(p, root);
    let stale = s
        .finding_ids
        .iter()
        .any(|id| !current.get(id).copied().unwrap_or(false));
    let decision = if stale {
        "needs_review"
    } else if !missing.is_empty() {
        "scope_not_established"
    } else {
        &s.outcome
    };
    Ok(
        json!({"revision":p.revision,"problem_id":q.id,"strategy_id":s.id,"decision":decision,"explanation":s.explanation.get(lang),"missing_premise_ids":missing,"additional_premise_ids":extra,"premises":p.premises.iter().filter(|x|s.premise_ids.contains(&x.id)).map(|x|premise_view(x,lang)).collect::<Vec<_>>(),"changes_required":s.changes_required.get(lang),"prior_findings":q.findings.iter().filter(|f|s.finding_ids.contains(&f.id)).map(|f|finding_view(f,lang,&current)).collect::<Vec<_>>(),"notice":if lang=="en"{"Applicability uses explicit premise IDs only. Premise truth and compatibility are not inferred. An exclusion applies only to the recorded strategy and scope, never to the whole theory."}else{"Anwendbarkeit verwendet nur explizite Voraussetzung-IDs. Wahrheit und Verträglichkeit der Voraussetzungen werden nicht erschlossen. Ein Ausschluss gilt nur für die dokumentierte Strategie und ihren Geltungsbereich, nie für die gesamte Theorie."}}),
    )
}

pub fn editable_record(p: &Program, id: &str, max_chars: usize) -> Result<Value> {
    require(
        (2000..=60000).contains(&max_chars),
        "max_chars must be between 2000 and 60000",
    )?;
    let q = p
        .problems
        .iter()
        .find(|q| q.id == id)
        .context("editable=true requires an exact problem ID")?;
    let used = problem_premise_ids(q);
    let premises: Vec<_> = p.premises.iter().filter(|x| used.contains(&x.id)).collect();
    let v = json!({"revision":p.revision,"problem":q,"premises":premises,"update_tool":"record_research_update","instructions":"Review original sources, edit this complete problem, supply expected_revision and a reason. Add new immutable premise definitions if needed. Findings and tasks cannot silently disappear. Source hashes and line ranges must match current files. This changes research navigation, not the original proof ledger."});
    let size = serde_json::to_string(&v)?.chars().count();
    require(
        size <= max_chars,
        format!("editable record needs {size} characters; increase max_chars"),
    )?;
    Ok(v)
}

/// A process-independent compare-and-swap update. Snapshots retain prior decisions for audit.
pub fn update(
    path: &Path,
    root: &Path,
    g: &KnowledgeGraph,
    expected: &str,
    problem: Value,
    new_premises: Vec<Premise>,
    reason: &str,
) -> Result<Value> {
    require(!reason.trim().is_empty(), "reason required")?;
    let parent = path.parent().context("research directory missing")?;
    fs::create_dir_all(parent)?;
    let lock = parent.join(".update.lock");
    let handle = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&lock)
        .context("research update in progress; retry after the writer completes")?;
    struct Guard(PathBuf);
    impl Drop for Guard {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.0);
        }
    }
    let _guard = Guard(lock);
    drop(handle);
    let mut p = load(path)?;
    require(
        p.revision == expected,
        format!("revision conflict: current revision is {}", p.revision),
    )?;
    let q: Problem = serde_json::from_value(problem)
        .context("problem must be a complete bilingual research record")?;
    let old = serde_json::to_vec_pretty(&p)?;
    // Reviewed findings are append-only. Correct a record explicitly by supplying a replacement
    // source/statement with the same stable ID; omission cannot erase an earlier exclusion.
    if let Some(prior) = p.problems.iter().find(|x| x.id == q.id) {
        for f in &prior.findings {
            require(
                q.findings.iter().any(|x| x.id == f.id),
                format!("cannot silently remove finding {}", f.id),
            )?;
        }
        for t in &prior.tasks {
            require(
                q.tasks.iter().any(|x| x.id == t.id),
                format!("cannot silently remove task {}", t.id),
            )?;
        }
        for s in &prior.strategies {
            require(
                q.strategies.iter().any(|x| x.id == s.id),
                format!("cannot silently remove strategy {}", s.id),
            )?;
        }
        for c in &prior.connections {
            require(
                q.connections.iter().any(|x| x.id == c.id),
                format!("cannot silently remove connection {}", c.id),
            )?;
        }
    }
    for premise in new_premises {
        if let Some(prior) = p.premises.iter().find(|x| x.id == premise.id) {
            require(
                serde_json::to_value(prior)? == serde_json::to_value(&premise)?,
                "premise identity is immutable; use a new ID for changed assumptions",
            )?;
        } else {
            p.premises.push(premise);
        }
    }
    let ids: HashSet<_> = g.nodes.iter().map(|n| n.id.as_str()).collect();
    for id in q
        .topics
        .iter()
        .chain(&q.goals)
        .chain(&q.reviewed_attempt_ids)
        .chain(q.findings.iter().flat_map(|f| &f.attempt_ids))
    {
        require(
            ids.contains(id.as_str()),
            format!("unregistered graph reference {id}"),
        )?;
    }
    let id = q.id.clone();
    if let Some(i) = p.problems.iter().position(|x| x.id == q.id) {
        p.problems[i] = q;
    } else {
        p.problems.push(q);
    }
    validate(&p)?;
    let mut cache = HashMap::new();
    for f in &p.problems.iter().find(|q| q.id == id).unwrap().findings {
        for s in &f.sources {
            require(
                source_current(s, root, &mut cache),
                format!("stale or missing source {}", s.path),
            )?;
        }
    }
    p.revision.clear();
    p.revision = digest(&serde_json::to_vec(&p)?)[..24].to_owned();
    let history = parent.join("history");
    fs::create_dir_all(&history)?;
    // Only hash-derived names are used; the caller cannot choose a file path.
    let old_hash = digest(&old);
    fs::write(history.join(format!("{old_hash}.json")), old)?;
    let tmp = parent.join(format!(".{}.tmp", p.revision));
    fs::write(&tmp, serde_json::to_vec_pretty(&p)?)?;
    fs::write(
        history.join(format!("{}.change.json", p.revision)),
        serde_json::to_vec_pretty(
            &json!({"revision":p.revision,"previous":expected,"problem_id":id,"reason":reason,"activation":"prepared; the canonical program or a later retained snapshot must contain this revision to establish activation"}),
        )?,
    )?;
    fs::rename(&tmp, path)?;
    Ok(
        json!({"revision":p.revision,"problem_id":id,"persisted":true,"ledger_changed":false,"notice":"Stored source-reviewed record; this operation does not establish mathematical validity or promote a ledger claim."}),
    )
}

#[cfg(test)]
#[path = "research_tests.rs"]
mod tests;
