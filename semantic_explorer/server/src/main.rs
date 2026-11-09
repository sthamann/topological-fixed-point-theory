use anyhow::{Context, Result};
use axum::{
    extract::{Query, State},
    http::StatusCode,
    response::{IntoResponse, Response},
    routing::get,
    Json, Router,
};
use clap::{Parser, Subcommand};
use serde::Deserialize;
use serde_json::{json, Value};
use std::{path::PathBuf, sync::Arc};
use tower_http::services::{ServeDir, ServeFile};

mod overview;
mod papers;
mod research;
mod service;
use service::{GraphStore, SearchRequest};

#[derive(Parser)]
#[command(
    name = "semantic-explorer",
    version,
    about = "TFPT semantic knowledge explorer"
)]
struct Args {
    #[arg(long, default_value = "../data/knowledge.json")]
    data: PathBuf,
    #[arg(long, default_value = "../..")]
    repo: PathBuf,
    #[arg(long, default_value = "../web")]
    web: PathBuf,
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    Serve {
        #[arg(long, default_value = "127.0.0.1")]
        host: String,
        #[arg(long, default_value_t = 8731)]
        port: u16,
    },
    Mcp,
    /// Write the current merged Atlas graph as JSON to stdout.
    Export,
}

#[derive(Clone)]
struct AppState {
    store: Arc<GraphStore>,
}

#[derive(Deserialize)]
struct SearchQuery {
    q: Option<String>,
    #[serde(rename = "type")]
    kind: Option<String>,
    limit: Option<usize>,
}
#[derive(Deserialize)]
struct IdQuery {
    id: String,
}
#[derive(Deserialize)]
struct SourceQuery {
    id: String,
    line_start: Option<usize>,
    line_end: Option<usize>,
}
#[derive(Deserialize)]
struct TraceQuery {
    id: String,
    depth: Option<usize>,
}
#[derive(Deserialize)]
struct ThemeQuery {
    id: String,
    limit: Option<usize>,
    offset: Option<usize>,
    include_related: Option<bool>,
}
#[derive(Deserialize)]
struct ChainQuery {
    id: String,
    depth: Option<usize>,
}
#[derive(Deserialize)]
struct AttemptsQuery {
    id: String,
    status: Option<String>,
    limit: Option<usize>,
    offset: Option<usize>,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse();
    let base = std::env::current_dir()?;
    let resolve = |p: PathBuf| if p.is_absolute() { p } else { base.join(p) };
    let store = Arc::new(GraphStore::new(resolve(args.data), resolve(args.repo))?);
    match args.command {
        Command::Serve { host, port } => {
            serve(AppState { store }, resolve(args.web), host, port).await
        }
        Command::Mcp => mcp::run(store).await,
        Command::Export => {
            let graph = store.export_graph().await?;
            println!("{}", serde_json::to_string(&graph)?);
            Ok(())
        }
    }
}

async fn serve(state: AppState, web_root: PathBuf, host: String, port: u16) -> Result<()> {
    let index = web_root.join("index.html");
    let app = app_router(state, web_root, index);
    let addr: std::net::SocketAddr = format!("{host}:{port}")
        .parse()
        .context("invalid listen address")?;
    let listener = tokio::net::TcpListener::bind(addr).await?;
    axum::serve(listener, app).await?;
    Ok(())
}

fn app_router(state: AppState, web_root: PathBuf, index: PathBuf) -> Router {
    Router::new()
        .route("/api/papers", get(api_list_papers))
        .route("/api/papers/outline", get(api_paper_outline))
        .route("/api/papers/section", get(api_paper_section))
        .route("/api/papers/search", get(api_search_papers))
        .route("/api/papers/evidence", get(api_paper_evidence))
        .route("/api/papers/coverage", get(api_paper_coverage))
        .route("/api/status", get(api_status))
        .route("/api/research", get(api_research))
        .route("/api/research-plan", get(api_research_plan))
        .route("/api/topics", get(api_topics))
        .route("/api/search", get(api_search).post(api_search_post))
        .route("/api/node", get(api_node))
        .route("/api/audit", get(api_audit))
        .route("/api/trace", get(api_trace))
        .route("/api/source", get(api_source))
        .route("/api/graph", get(api_graph))
        .route("/api/theme", get(api_theme))
        .route("/api/chain", get(api_chain))
        .route("/api/attempts", get(api_attempts))
        .route("/api/brief", get(api_brief))
        .route("/api/context", get(api_context))
        .route("/api/theory-map", get(api_theory_map))
        .fallback_service(ServeDir::new(&web_root).not_found_service(ServeFile::new(index)))
        .with_state(state)
}

async fn api_list_papers(
    State(s): State<AppState>,
    Query(q): Query<papers::PaperArgs>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .paper_tool("list_papers", &serde_json::to_value(q)?)
            .await?,
    ))
}
async fn api_paper_outline(
    State(s): State<AppState>,
    Query(q): Query<papers::PaperArgs>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .paper_tool("paper_outline", &serde_json::to_value(q)?)
            .await?,
    ))
}
async fn api_paper_section(
    State(s): State<AppState>,
    Query(q): Query<papers::PaperArgs>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .paper_tool("paper_section", &serde_json::to_value(q)?)
            .await?,
    ))
}
async fn api_search_papers(
    State(s): State<AppState>,
    Query(q): Query<papers::PaperArgs>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .paper_tool("search_papers", &serde_json::to_value(q)?)
            .await?,
    ))
}
async fn api_paper_evidence(
    State(s): State<AppState>,
    Query(q): Query<papers::PaperArgs>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .paper_tool("paper_evidence", &serde_json::to_value(q)?)
            .await?,
    ))
}
async fn api_paper_coverage(
    State(s): State<AppState>,
    Query(q): Query<papers::PaperArgs>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .paper_tool("paper_coverage", &serde_json::to_value(q)?)
            .await?,
    ))
}

#[derive(Deserialize)]
struct ResearchQuery {
    query: Option<String>,
    lang: Option<String>,
    max_chars: Option<usize>,
    offset: Option<usize>,
    detail_offset: Option<usize>,
    expected_revision: Option<String>,
}
#[derive(Deserialize)]
struct PlanQuery {
    problem_id: String,
    strategy_id: String,
    premise_ids: Option<String>,
    lang: Option<String>,
}
async fn api_research(
    State(s): State<AppState>,
    Query(q): Query<ResearchQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .research_state(
                q.query.as_deref().unwrap_or(""),
                q.lang.as_deref().unwrap_or("de"),
                q.max_chars.unwrap_or(24000),
                q.offset.unwrap_or(0),
                q.detail_offset.unwrap_or(0),
                q.expected_revision.as_deref(),
            )
            .await?,
    ))
}
async fn api_research_plan(
    State(s): State<AppState>,
    Query(q): Query<PlanQuery>,
) -> Result<Json<Value>, ApiError> {
    let premises = q
        .premise_ids
        .unwrap_or_default()
        .split(',')
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect::<Vec<_>>();
    Ok(Json(
        s.store
            .assess_research_plan(
                &q.problem_id,
                &q.strategy_id,
                &premises,
                q.lang.as_deref().unwrap_or("de"),
            )
            .await?,
    ))
}
async fn api_status(State(s): State<AppState>) -> Result<Json<Value>, ApiError> {
    Ok(Json(s.store.status().await?))
}
async fn api_topics(State(s): State<AppState>) -> Result<Json<Value>, ApiError> {
    Ok(Json(s.store.topics().await?))
}
async fn api_search(
    State(s): State<AppState>,
    Query(q): Query<SearchQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .search(SearchRequest {
                query: q.q.unwrap_or_default(),
                kind: q.kind,
                limit: q.limit,
            })
            .await?,
    ))
}
async fn api_search_post(
    State(s): State<AppState>,
    Json(q): Json<SearchRequest>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(s.store.search(q).await?))
}
async fn api_node(
    State(s): State<AppState>,
    Query(q): Query<IdQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(s.store.dossier(&q.id).await?))
}
async fn api_audit(
    State(s): State<AppState>,
    Query(q): Query<IdQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(s.store.audit_inputs(&q.id).await?))
}
async fn api_trace(
    State(s): State<AppState>,
    Query(q): Query<TraceQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store.trace(&q.id, "both", q.depth.unwrap_or(2)).await?,
    ))
}
async fn api_source(
    State(s): State<AppState>,
    Query(q): Query<SourceQuery>,
) -> Result<Response, ApiError> {
    Ok(Json(s.store.source(&q.id, q.line_start, q.line_end).await?).into_response())
}
async fn api_graph(
    State(s): State<AppState>,
    Query(q): Query<GraphQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .graph_slice(q.id.as_deref(), q.limit.unwrap_or(500))
            .await?,
    ))
}
async fn api_theme(
    State(s): State<AppState>,
    Query(q): Query<ThemeQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .theme(
                &q.id,
                q.limit.unwrap_or(100),
                q.offset.unwrap_or(0),
                q.include_related.unwrap_or(false),
            )
            .await?,
    ))
}
async fn api_chain(
    State(s): State<AppState>,
    Query(q): Query<ChainQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store.solution_chain(&q.id, q.depth.unwrap_or(4)).await?,
    ))
}
async fn api_attempts(
    State(s): State<AppState>,
    Query(q): Query<AttemptsQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .list_attempts(
                &q.id,
                q.status.as_deref(),
                q.limit.unwrap_or(100),
                q.offset.unwrap_or(0),
            )
            .await?,
    ))
}
async fn api_brief(
    State(s): State<AppState>,
    Query(q): Query<BriefQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .research_brief(&q.id, q.max_chars.unwrap_or(12_000), q.offset.unwrap_or(0))
            .await?,
    ))
}
async fn api_context(
    State(s): State<AppState>,
    Query(q): Query<BriefQuery>,
) -> Result<Json<Value>, ApiError> {
    Ok(Json(
        s.store
            .research_context(&q.id, q.max_chars.unwrap_or(12_000), q.offset.unwrap_or(0))
            .await?,
    ))
}
async fn api_theory_map(State(s): State<AppState>) -> Result<Json<Value>, ApiError> {
    Ok(Json(s.store.theory_map().await?))
}
#[derive(Deserialize)]
struct BriefQuery {
    id: String,
    max_chars: Option<usize>,
    offset: Option<usize>,
}
#[derive(Deserialize)]
struct GraphQuery {
    id: Option<String>,
    limit: Option<usize>,
}

struct ApiError(anyhow::Error);
impl<E: Into<anyhow::Error>> From<E> for ApiError {
    fn from(e: E) -> Self {
        Self(e.into())
    }
}
impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        let msg = self.0.to_string();
        let status = if msg.contains("not found") {
            StatusCode::NOT_FOUND
        } else if msg.contains("limit")
            || msg.contains("invalid")
            || msg.contains("unsafe")
            || msg.contains("depth must")
            || msg.contains("max_chars")
        {
            StatusCode::BAD_REQUEST
        } else {
            StatusCode::INTERNAL_SERVER_ERROR
        };
        (status, Json(json!({"error": msg}))).into_response()
    }
}

mod mcp {
    use super::*;
    use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};

    const VERSIONS: &[&str] = &["2024-11-05", "2025-03-26", "2025-06-18"];

    pub async fn run(store: Arc<GraphStore>) -> Result<()> {
        let stdin = tokio::io::stdin();
        let mut lines = BufReader::new(stdin).lines();
        let mut stdout = tokio::io::stdout();
        while let Some(line) = lines.next_line().await? {
            if line.trim().is_empty() {
                continue;
            }
            let response = match serde_json::from_str::<Value>(&line) {
                Ok(req) => handle(&store, req).await,
                Err(e) => Some(
                    json!({"jsonrpc":"2.0","id":null,"error":{"code":-32700,"message":format!("Parse error: {e}")}}),
                ),
            };
            if let Some(response) = response {
                let mut bytes = serde_json::to_vec(&response)?;
                bytes.push(b'\n');
                stdout.write_all(&bytes).await?;
                stdout.flush().await?;
            }
        }
        Ok(())
    }

    async fn handle(store: &GraphStore, req: Value) -> Option<Value> {
        let id = req.get("id").cloned();
        let method = req.get("method").and_then(Value::as_str).unwrap_or("");
        let params = req.get("params").cloned().unwrap_or_else(|| json!({}));
        let Some(id) = id else {
            return None;
        };
        let response = match method {
            "initialize" => {
                let requested = params
                    .get("protocolVersion")
                    .and_then(Value::as_str)
                    .unwrap_or(VERSIONS[0]);
                let selected = if VERSIONS.contains(&requested) {
                    requested
                } else {
                    VERSIONS[0]
                };
                Ok(
                    json!({"protocolVersion":selected,"capabilities":{"tools":{"listChanged":false}},"serverInfo":{"name":"tfpt-semantic-explorer","version":env!("CARGO_PKG_VERSION")},"instructions":"Use this as the sole current TFPT/Universalraum theory atlas. For overall TFPT status start with the scope and evidence areas in theory_overview, separate current work from established results, and then use research_state(query in natural language, lang, max_chars=12000). Read its full problem packet when packet_state is absent. For packet_state=partial, use budget.continuation with detail_offset and treat only each returned scoped record as complete; the page is not proof of the whole problem. Use assess_research_plan with an exact strategy ID and explicit premise IDs before repeating work. Persist source-reviewed changes through record_research_update with expected_revision. Missing reviewed coverage is not scientific openness. For papers use list_papers/search_papers, paper_outline, paper_section, paper_evidence and read_source; paper_coverage separates source completeness from unreviewed AI summaries. For other unreviewed corpus material use research_context and read_source. Preserve scope, conjunctive premises, hidden inputs and negative results. Graph paths and AI annotations are navigation, not proofs. Consult prior attempts before proposing another run; name the changed premise. Source or coverage gaps remain explicit. Original ledger and cited sources are authoritative."}),
                )
            }
            "ping" => Ok(json!({})),
            "tools/list" => Ok(json!({"tools":tool_specs()})),
            "tools/call" => call_tool(store, &params).await,
            _ => {
                return Some(
                    json!({"jsonrpc":"2.0","id":id,"error":{"code":-32601,"message":"Method not found"}}),
                )
            }
        };
        Some(match response {
            Ok(result) => json!({"jsonrpc":"2.0","id":id,"result":result}),
            Err(error) => {
                json!({"jsonrpc":"2.0","id":id,"error":{"code":-32602,"message":error.to_string()}})
            }
        })
    }

    fn tool_specs() -> Vec<Value> {
        let mut specs = vec![
            tool("research_state", "Read durable source-reviewed question state using natural language or topic/goal/problem ID. Separates resolved subresults, conditional findings, scoped exclusions, task prerequisites and unreviewed attempts. Whole-theory queries in German or English return evidence areas and a paged topic index before any research detail; use an exact problem ID for details. Oversized problems return explicitly partial, record-complete pages with exact premise, dependency and source handles; follow every argument in budget.continuation. A changed expected_revision rejects the continuation instead of mixing pages. A partial page is never proof of the whole problem. editable=true still returns only the complete authoring record. problem_index lists matches, not proofs. Original sources remain authoritative.", json!({"type":"object","properties":{"query":{"type":"string","default":""},"editable":{"type":"boolean","default":false,"description":"Return a complete editable authoring record for an exact problem ID; use max_chars=60000. detail_offset does not apply."},"lang":{"type":"string","enum":["de","en"],"default":"de"},"max_chars":{"type":"integer","default":12000,"minimum":2000,"maximum":60000},"offset":{"type":"integer","default":0,"minimum":0},"detail_offset":{"type":"integer","default":0,"minimum":0,"description":"Inner record offset for the oversized problem at offset. Use the exact value in budget.continuation."},"expected_revision":{"type":"string","description":"Revision from budget.continuation; rejects pagination across an Atlas update."}}})),
            tool("assess_research_plan", "Check a recorded strategy against explicitly supplied premise IDs. Exact strategy IDs permit reuse/exclusion only in the recorded scope; free text returns candidates and never an exclusion. Missing premises or changed sources block decisive reuse. No match is not proof of novelty; premise truth/compatibility is not inferred.", json!({"type":"object","properties":{"problem_id":{"type":"string"},"strategy_id":{"type":"string"},"premise_ids":{"type":"array","items":{"type":"string"},"default":[]},"lang":{"type":"string","enum":["de","en"],"default":"de"}},"required":["problem_id","strategy_id"]})),
            tool("record_research_update", "Persist one complete reviewed problem record in the canonical Atlas, with immutable premise identities, source SHA256/line ranges, dependency validation, optimistic revision check, and retained history. Read semantic_explorer/research/program.json for the authoring schema; research_state supplies the current revision. Never use this to promote ledger proof status. Findings/tasks cannot silently disappear. The caller must review original sources; passing validation alone does not prove the result.", json!({"type":"object","properties":{"expected_revision":{"type":"string"},"problem":{"type":"object"},"premises":{"type":"array","items":{"type":"object"},"default":[]},"reason":{"type":"string"}},"required":["expected_revision","problem","reason"]})),
            tool("search_knowledge", "Search graph records with bilingual lexical matching and optional configured embeddings; directly matching authored theory-map patterns may appear as theory_pattern results and preserve their review state and source references.", json!({"type":"object","properties":{"query":{"type":"string"},"type":{"type":"string","description":"Optional node type filter; use theory_pattern to search authored theory-map patterns only."},"limit":{"type":"integer","maximum":50}},"required":["query"]})),
            tool("get_dossier", "Return a compact node dossier with its key results, premises, assumptions, positive findings, exclusions, explicit open steps, own registered source paths, neighbor counts, and omission counts. Use read_source for cited passages. Set expanded=true for the complete connected dossier when its serialized size fits max_chars. The HTTP dossier remains full.", json!({"type":"object","properties":{"id":{"type":"string"},"max_chars":{"type":"integer","default":24000,"minimum":1000,"maximum":500000},"neighbor_limit":{"type":"integer","default":8,"minimum":1,"maximum":50},"expanded":{"type":"boolean","default":false}},"required":["id"]})),
            tool("audit_inputs", "Inspect recorded input origins, impacts and proof gaps for a claim, or syntactic input candidates plus bounded related claim audits for a script.", json!({"type":"object","properties":{"id":{"type":"string"}},"required":["id"]})),
            tool("trace_reasoning", "Traverse explicit graph edges. Edge labels are reported as stored and are not upgraded into proof claims.", json!({"type":"object","properties":{"id":{"type":"string"},"direction":{"type":"string","enum":["inbound","outbound","both"]},"depth":{"type":"integer","maximum":5}},"required":["id"]})),
            tool("list_topics", "List topic nodes.", json!({"type":"object","properties":{}})),
            tool("graph_status", "Summarize graph metadata and node/edge counts.", json!({"type":"object","properties":{}})),
            tool("read_source", "Read only registered source records after repository-path and hash validation. Optional inclusive line ranges stay within the registered source extent and are limited to 200 lines; otherwise a bounded default window is returned.", json!({"type":"object","properties":{"source_id":{"type":"string"},"line_start":{"type":"integer","minimum":1},"line_end":{"type":"integer","minimum":1}},"required":["source_id"]})),
            tool("topic_overview", "Show explicit theme membership and its recorded assessment categories; established status remains scoped and historical attempts do not prove the claim. Related-role claims are excluded by default and can be included explicitly.", json!({"type":"object","properties":{"id":{"type":"string"},"limit":{"type":"integer","maximum":100},"offset":{"type":"integer","minimum":0},"include_related":{"type":"boolean","default":false}},"required":["id"]})),
            tool("solution_chain", "Return a bounded graph of explicit reasoning edges. References are separate. `traversal_complete` reports only whether this traversal hit a bound; `all_premises_complete` is null unless a bound omitted stored edges, because traversal cannot establish premise completeness or validity.", json!({"type":"object","properties":{"id":{"type":"string"},"depth":{"type":"integer","maximum":4}},"required":["id"]})),
            tool("list_attempts", "List explicitly linked historical attempts or records for a topic or claim; attempt status does not establish the linked claim.", json!({"type":"object","properties":{"id":{"type":"string"},"status":{"type":"string"},"limit":{"type":"integer","maximum":100},"offset":{"type":"integer","minimum":0}},"required":["id"]})),
            tool("research_brief", "Create a compact deterministic evidence pack for a topic, claim, explicit goal, or exact attempt/contract ID. Attempt/contract roots include only their directly tested claims and their own record; they do not fan out through topics or related attempts. Includes recorded annotations, scoped assessments, gaps, inputs, compact attempt conclusions, and registered source paths. Goal IDs use only explicit claim membership and goal-linked attempts. This is not proof and makes no LLM calls. max_chars bounds the inner JSON; the transport envelope adds framing.", json!({"type":"object","properties":{"id":{"type":"string"},"max_chars":{"type":"integer","default":12000,"minimum":1000,"maximum":30000},"offset":{"type":"integer","minimum":0}},"required":["id"]})),
            tool("research_context", "Prepare a compact deterministic handoff for a topic, claim, explicit goal, or exact attempt/contract ID. Attempt/contract roots include only directly tested claims and their own record plus typed direct relations, never topic fanout. Includes status/scope, shared IDs, recorded connections, prior attempts, scoped exclusions, input records, source IDs and registered paths. Goal membership is explicit; related topic claims remain separately sampled. Equal labels are not identity and stored links are not proof. max_chars bounds the inner JSON, with offset for pages.", json!({"type":"object","properties":{"id":{"type":"string"},"max_chars":{"type":"integer","default":12000,"minimum":1000,"maximum":30000},"offset":{"type":"integer","minimum":0}},"required":["id"]})),
            tool("compare_aspects", "Compare two to six graph aspects by explicit upstream links, seeds, status, recorded assumptions and source IDs. Shared features are not a logical conclusion. max_chars bounds the compact inner JSON content; the JSON-RPC transport envelope adds its ordinary framing.", json!({"type":"object","properties":{"ids":{"type":"array","items":{"type":"string"},"minItems":2,"maxItems":6},"max_chars":{"type":"integer","default":12000,"minimum":1000,"maximum":30000}},"required":["ids"]})),
            tool("theory_overview", "Return a compact readonly summary of the locally saved theory map, if available.", json!({"type":"object","properties":{}})),
        ];
        specs.extend(crate::papers::tool_specs());
        specs
    }
    fn tool(name: &str, description: &str, input_schema: Value) -> Value {
        json!({"name":name,"description":description,"inputSchema":input_schema})
    }

    async fn call_tool(store: &GraphStore, params: &Value) -> Result<Value> {
        let name = params
            .get("name")
            .and_then(Value::as_str)
            .context("tools/call requires name")?;
        let a = params
            .get("arguments")
            .cloned()
            .unwrap_or_else(|| json!({}));
        let result = match name {
            "research_state" if a.get("editable").and_then(Value::as_bool).unwrap_or(false) => {
                store
                    .research_record(
                        arg_str(&a, "query")?,
                        a.get("max_chars").and_then(Value::as_u64).unwrap_or(12000) as usize,
                    )
                    .await?
            }
            "research_state" => {
                store
                    .research_state(
                        a.get("query").and_then(Value::as_str).unwrap_or(""),
                        a.get("lang").and_then(Value::as_str).unwrap_or("de"),
                        a.get("max_chars").and_then(Value::as_u64).unwrap_or(12000) as usize,
                        a.get("offset").and_then(Value::as_u64).unwrap_or(0) as usize,
                        a.get("detail_offset").and_then(Value::as_u64).unwrap_or(0) as usize,
                        a.get("expected_revision").and_then(Value::as_str),
                    )
                    .await?
            }
            "assess_research_plan" => {
                store
                    .assess_research_plan(
                        arg_str(&a, "problem_id")?,
                        arg_str(&a, "strategy_id")?,
                        &serde_json::from_value::<Vec<String>>(
                            a.get("premise_ids").cloned().unwrap_or(json!([])),
                        )?,
                        a.get("lang").and_then(Value::as_str).unwrap_or("de"),
                    )
                    .await?
            }
            "record_research_update" => {
                store
                    .record_research_update(
                        arg_str(&a, "expected_revision")?,
                        a.get("problem").cloned().context("problem required")?,
                        serde_json::from_value(a.get("premises").cloned().unwrap_or(json!([])))?,
                        arg_str(&a, "reason")?,
                    )
                    .await?
            }
            "list_papers" | "paper_outline" | "paper_section" | "search_papers"
            | "paper_evidence" | "paper_coverage" => store.paper_tool(name, &a).await?,
            "search_knowledge" => {
                store
                    .search(SearchRequest {
                        query: arg_str(&a, "query")?.to_owned(),
                        kind: a.get("type").and_then(Value::as_str).map(str::to_owned),
                        limit: a.get("limit").and_then(Value::as_u64).map(|n| n as usize),
                    })
                    .await?
            }
            "get_dossier" => {
                let id = arg_str(&a, "id")?;
                let max_chars =
                    a.get("max_chars").and_then(Value::as_u64).unwrap_or(24_000) as usize;
                if !(1_000..=500_000).contains(&max_chars) {
                    anyhow::bail!("max_chars must be between 1000 and 500000");
                }
                if a.get("expanded").and_then(Value::as_bool).unwrap_or(false) {
                    let full = store.dossier(id).await?;
                    let chars = serde_json::to_string(&full)?.chars().count();
                    if chars > max_chars {
                        anyhow::bail!("full dossier is {chars} characters, above max_chars={max_chars}; use the compact dossier or HTTP dossier");
                    }
                    full
                } else {
                    store
                        .dossier_compact(
                            id,
                            max_chars,
                            a.get("neighbor_limit").and_then(Value::as_u64).unwrap_or(8) as usize,
                        )
                        .await?
                }
            }
            "audit_inputs" => store.audit_inputs(arg_str(&a, "id")?).await?,
            "trace_reasoning" => {
                store
                    .trace(
                        arg_str(&a, "id")?,
                        a.get("direction").and_then(Value::as_str).unwrap_or("both"),
                        a.get("depth").and_then(Value::as_u64).unwrap_or(2) as usize,
                    )
                    .await?
            }
            "list_topics" => store.topics().await?,
            "graph_status" => {
                let mut status = store.status().await?;
                status["service_identity"] = json!({
                    "name":"tfpt-semantic-explorer",
                    "version":env!("CARGO_PKG_VERSION"),
                    "tool_catalog_count":tool_specs().len(),
                    "notice":"The catalog count helps distinguish server generations; use tools/list for the exact tool names."
                });
                status
            }
            "read_source" => {
                store
                    .source(
                        arg_str(&a, "source_id")?,
                        a.get("line_start")
                            .and_then(Value::as_u64)
                            .map(|n| n as usize),
                        a.get("line_end")
                            .and_then(Value::as_u64)
                            .map(|n| n as usize),
                    )
                    .await?
            }
            "topic_overview" => {
                store
                    .theme(
                        arg_str(&a, "id")?,
                        a.get("limit").and_then(Value::as_u64).unwrap_or(100) as usize,
                        a.get("offset").and_then(Value::as_u64).unwrap_or(0) as usize,
                        a.get("include_related")
                            .and_then(Value::as_bool)
                            .unwrap_or(false),
                    )
                    .await?
            }
            "solution_chain" => {
                store
                    .solution_chain(
                        arg_str(&a, "id")?,
                        a.get("depth").and_then(Value::as_u64).unwrap_or(4) as usize,
                    )
                    .await?
            }
            "list_attempts" => {
                store
                    .list_attempts(
                        arg_str(&a, "id")?,
                        a.get("status").and_then(Value::as_str),
                        a.get("limit").and_then(Value::as_u64).unwrap_or(100) as usize,
                        a.get("offset").and_then(Value::as_u64).unwrap_or(0) as usize,
                    )
                    .await?
            }
            "research_brief" => {
                store
                    .research_brief(
                        arg_str(&a, "id")?,
                        a.get("max_chars").and_then(Value::as_u64).unwrap_or(12_000) as usize,
                        a.get("offset").and_then(Value::as_u64).unwrap_or(0) as usize,
                    )
                    .await?
            }
            "research_context" => {
                store
                    .research_context(
                        arg_str(&a, "id")?,
                        a.get("max_chars").and_then(Value::as_u64).unwrap_or(12_000) as usize,
                        a.get("offset").and_then(Value::as_u64).unwrap_or(0) as usize,
                    )
                    .await?
            }
            "compare_aspects" => {
                let ids = a
                    .get("ids")
                    .and_then(Value::as_array)
                    .context("compare_aspects requires ids array")?
                    .iter()
                    .map(|id| {
                        id.as_str()
                            .map(str::to_owned)
                            .context("ids must be strings")
                    })
                    .collect::<Result<Vec<_>>>()?;
                store
                    .compare_aspects(
                        &ids,
                        a.get("max_chars").and_then(Value::as_u64).unwrap_or(12_000) as usize,
                    )
                    .await?
            }
            "theory_overview" => store.theory_overview().await?,
            _ => anyhow::bail!("unknown tool: {name}"),
        };
        let compact_text = matches!(
            name,
            "list_papers"
                | "paper_outline"
                | "paper_section"
                | "search_papers"
                | "paper_evidence"
                | "paper_coverage"
                | "research_state"
                | "assess_research_plan"
                | "record_research_update"
                | "get_dossier"
                | "search_knowledge"
                | "graph_status"
                | "research_brief"
                | "research_context"
                | "compare_aspects"
                | "theory_overview"
                | "solution_chain"
                | "topic_overview"
                | "list_attempts"
                | "read_source"
        );
        let text = if compact_text {
            serde_json::to_string(&result)?
        } else {
            serde_json::to_string_pretty(&result)?
        };
        let mut tool_result = json!({"content":[{"type":"text","text":text}],"isError":false});
        if !compact_text {
            tool_result["structuredContent"] = result;
        }
        Ok(tool_result)
    }
    fn arg_str<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
        v.get(key)
            .and_then(Value::as_str)
            .filter(|s| !s.is_empty())
            .with_context(|| format!("missing string argument {key}"))
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use axum::body::{to_bytes, Body};
        use http::Request;
        use sha2::{Digest, Sha256};
        use tower::ServiceExt;
        #[test]
        fn supported_protocol_versions_are_explicit() {
            assert!(VERSIONS.contains(&"2024-11-05"));
            assert!(VERSIONS.contains(&"2025-03-26"));
            assert!(VERSIONS.contains(&"2025-06-18"));
        }
        #[test]
        fn tool_list_has_expected_tools() {
            let specs = tool_specs();
            let names: Vec<_> = specs.iter().filter_map(|x| x["name"].as_str()).collect();
            assert_eq!(names.len(), 23);
            for name in crate::papers::TOOLS {
                assert!(names.contains(name));
            }
            assert!(names.contains(&"research_state"));
            assert!(names.contains(&"assess_research_plan"));
            assert!(names.contains(&"record_research_update"));
            assert!(names.contains(&"trace_reasoning"));
            assert!(names.contains(&"get_dossier"));
            assert!(names.contains(&"audit_inputs"));
            assert!(names.contains(&"topic_overview"));
            assert!(names.contains(&"solution_chain"));
            assert!(names.contains(&"list_attempts"));
            assert!(names.contains(&"research_brief"));
            assert!(names.contains(&"research_context"));
            assert!(names.contains(&"compare_aspects"));
            assert!(names.contains(&"theory_overview"));
        }
        #[tokio::test]
        async fn paper_tools_are_dispatched_through_mcp_and_http() {
            let root = std::env::temp_dir().join(format!(
                "tfpt-paper-api-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            std::fs::create_dir_all(&root).unwrap();
            std::fs::write(root.join("source.tex"), "a registered source\n").unwrap();
            let sha = format!("{:x}", Sha256::digest(b"a registered source\n"));
            let data = root.join("knowledge.json");
            std::fs::write(&data, r#"{"meta":{},"nodes":[],"edges":[],"sources":[]}"#).unwrap();
            std::fs::write(root.join("papers.json"),serde_json::to_vec(&json!({"meta":{"coverage":{"expected_registry_coverage":{"sections":1}}},"nodes":[{"id":"document:test","type":"paper_document","label":"Paper","attrs":{"scope":"public","source_refs":["paper-source:test"],"ordered_section_ids":["section:test"]}},{"id":"section:test","type":"paper_section","label":"A native clock","summary":"An unreviewed summary","attrs":{"scope":"public","source_refs":["paper-source:test"],"source_text":"native clock in complete body","summary_state":"pending"}}],"edges":[],"sources":[{"id":"paper-source:test","path":"source.tex","sha256":sha,"line_start":1,"line_end":1}]})).unwrap()).unwrap();
            let store = Arc::new(GraphStore::new(data, root.clone()).unwrap());
            for tool in crate::papers::TOOLS {
                let reply=handle(&store,json!({"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":tool,"arguments":{"query":"clock","document_id":"document:test","section_id":"section:test","id":"section:test","max_chars":4000}}})).await.unwrap();
                assert!(reply.get("error").is_none(), "{tool}: {reply}");
                let text = reply["result"]["content"][0]["text"].as_str().unwrap();
                assert!(text.chars().count() <= 4000);
                let value: Value = serde_json::from_str(text).unwrap();
                assert_eq!(value["budget"]["actual_chars"], text.chars().count());
            }
            let app = app_router(AppState { store }, root.clone(), root.join("index.html"));
            let response = app
                .oneshot(
                    Request::builder()
                        .uri("/api/papers/section?section_id=section:test&max_chars=4000")
                        .body(Body::empty())
                        .unwrap(),
                )
                .await
                .unwrap();
            assert_eq!(response.status(), StatusCode::OK);
            let body = to_bytes(response.into_body(), 10000).await.unwrap();
            let value: Value = serde_json::from_slice(&body).unwrap();
            assert_eq!(value["section"]["id"], "section:test");
            std::fs::remove_dir_all(root).unwrap();
        }
        #[tokio::test]
        async fn mcp_initialize_and_search_return_protocol_results() {
            let root = std::env::temp_dir().join(format!(
                "semantic-explorer-mcp-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            std::fs::create_dir_all(&root).unwrap();
            let source_text = (1..=10)
                .map(|n| format!("source line {n}"))
                .collect::<Vec<_>>()
                .join("\n");
            std::fs::write(root.join("source.txt"), &source_text).unwrap();
            let source_hash = format!("{:x}", Sha256::digest(source_text.as_bytes()));
            let data = root.join("knowledge.json");
            std::fs::write(&data, serde_json::to_vec(&json!({"meta":{"schema_version":1},"nodes":[{"id":"claim:EM.FP.01","type":"claim","label":"Fixed point existence","summary":"The fixed point exists","attrs":{"ledger":{"claim":"fixed point"},"assessment":{"category":"open","scope":"test scope","authority":"ledger"}}},{"id":"claim:SECOND","type":"claim","label":"Second related claim","summary":"A second claim","attrs":{}},{"id":"topic:test","type":"topic","label":"Test topic","summary":"A small test topic","attrs":{}}],"edges":[{"id":"topic-claim","src":"topic:test","dst":"claim:EM.FP.01","type":"about"}],"sources":[{"id":"source:paper","path":"source.txt","sha256":source_hash,"line_count":10}]})).unwrap()).unwrap();
            let store = GraphStore::new(data, root.clone()).unwrap();
            let init=handle(&store,json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}})).await.unwrap();
            assert_eq!(init["result"]["protocolVersion"], "2025-06-18");
            let search=handle(&store,json!({"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"search_knowledge","arguments":{"query":"fixed point"}}})).await.unwrap();
            assert_eq!(search["result"]["isError"], false);
            assert!(search["result"].get("structuredContent").is_none());
            let search_text = search["result"]["content"][0]["text"].as_str().unwrap();
            let parsed_search: Value = serde_json::from_str(search_text).unwrap();
            assert_eq!(parsed_search["results"][0]["id"], "claim:EM.FP.01");
            assert!(parsed_search["results"][0].get("attrs").is_none());
            let dossier = handle(&store, json!({"jsonrpc":"2.0","id":11,"method":"tools/call","params":{"name":"get_dossier","arguments":{"id":"claim:EM.FP.01"}}})).await.unwrap();
            assert_eq!(dossier["result"]["isError"], false);
            assert!(dossier["result"].get("structuredContent").is_none());
            let parsed_dossier: Value =
                serde_json::from_str(dossier["result"]["content"][0]["text"].as_str().unwrap())
                    .unwrap();
            assert_eq!(parsed_dossier["node"]["id"], "claim:EM.FP.01");
            assert!(parsed_dossier.get("omissions").is_some());
            let status = handle(&store, json!({"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"graph_status","arguments":{}}})).await.unwrap();
            assert!(status["result"].get("structuredContent").is_none());
            let status_text = status["result"]["content"][0]["text"].as_str().unwrap();
            let parsed_status: Value = serde_json::from_str(status_text).unwrap();
            assert_eq!(parsed_status["counts"]["nodes"], 3);
            let chain = handle(&store, json!({"jsonrpc":"2.0","id":9,"method":"tools/call","params":{"name":"solution_chain","arguments":{"id":"claim:EM.FP.01","depth":2}}})).await.unwrap();
            assert!(chain["result"].get("structuredContent").is_none());
            let chain_text = chain["result"]["content"][0]["text"].as_str().unwrap();
            assert_eq!(
                serde_json::from_str::<Value>(chain_text).unwrap()["root"]["id"],
                "claim:EM.FP.01"
            );
            let theme = handle(&store, json!({"jsonrpc":"2.0","id":10,"method":"tools/call","params":{"name":"topic_overview","arguments":{"id":"topic:test","limit":100}}})).await.unwrap();
            assert!(theme["result"].get("structuredContent").is_none());
            let theme_text = theme["result"]["content"][0]["text"].as_str().unwrap();
            let parsed_theme: Value = serde_json::from_str(theme_text).unwrap();
            assert_eq!(
                parsed_theme["open"][0]["attrs"]["assessment"]["scope"],
                "test scope"
            );
            assert!(parsed_theme["open"][0]["attrs"].get("ledger").is_none());
            let source = handle(&store, json!({"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"read_source","arguments":{"source_id":"source:paper","line_start":4,"line_end":6}}})).await.unwrap();
            assert!(source["result"].get("structuredContent").is_none());
            let source_text = source["result"]["content"][0]["text"].as_str().unwrap();
            let source_value: Value = serde_json::from_str(source_text).unwrap();
            assert_eq!(source_value["line_start"], 4);
            assert_eq!(source_value["line_end"], 6);
            assert_eq!(
                source_value["text"],
                "source line 4\nsource line 5\nsource line 6"
            );
            // A consumer must be able to forward the returned continuation
            // unchanged to the real MCP dispatcher.
            let source_next=handle(&store,json!({"jsonrpc":"2.0","id":12,"method":"tools/call",
                "params":{"name":"read_source","arguments":source_value["extent"]["next_request"]}})).await.unwrap();
            assert!(source_next.get("error").is_none(), "{source_next}");
            assert_eq!(source_next["result"]["isError"], false);
            let source_next: Value = serde_json::from_str(
                source_next["result"]["content"][0]["text"]
                    .as_str()
                    .unwrap(),
            )
            .unwrap();
            assert_eq!(source_next["line_start"], 7);
            assert_eq!(source_next["line_end"], 10);
            assert_eq!(
                source_next["text"],
                "source line 7\nsource line 8\nsource line 9\nsource line 10"
            );
            assert!(source_next["extent"]["next_request"].is_null());

            let plan_directory = root.join("docs/research/fundamental-plan");
            std::fs::create_dir_all(&plan_directory).unwrap();
            std::fs::write(plan_directory.join("progress.json"),serde_json::to_vec(&json!({
                "revision":1,"next_task":"F1","tasks":[],"continuation":"original source ".repeat(1000)
            })).unwrap()).unwrap();
            let plan_state=handle(&store,json!({"jsonrpc":"2.0","id":13,"method":"tools/call",
                "params":{"name":"research_state","arguments":{"query":"current priority","lang":"en","max_chars":12000}}})).await.unwrap();
            assert!(plan_state.get("error").is_none(), "{plan_state}");
            let plan_state: Value =
                serde_json::from_str(plan_state["result"]["content"][0]["text"].as_str().unwrap())
                    .unwrap();
            let plan_next = plan_state["current_plan"]["next_request"].clone();
            assert_eq!(plan_next["tool"], "read_source");
            let plan_source = handle(
                &store,
                json!({"jsonrpc":"2.0","id":14,"method":"tools/call",
                "params":{"name":plan_next["tool"],"arguments":plan_next["arguments"]}}),
            )
            .await
            .unwrap();
            assert!(plan_source.get("error").is_none(), "{plan_source}");
            assert_eq!(plan_source["result"]["isError"], false);
            let plan_source: Value = serde_json::from_str(
                plan_source["result"]["content"][0]["text"]
                    .as_str()
                    .unwrap(),
            )
            .unwrap();
            assert_eq!(plan_source["id"], "source:fundamental-plan:progress.json");
            assert!(plan_source["text"]
                .as_str()
                .unwrap()
                .contains("original source"));
            let brief=handle(&store,json!({"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"research_brief","arguments":{"id":"claim:EM.FP.01","max_chars":12000}}})).await.unwrap();
            assert_eq!(brief["result"]["isError"], false);
            assert!(brief["result"].get("structuredContent").is_none());
            let brief_text = brief["result"]["content"][0]["text"].as_str().unwrap();
            assert!(brief_text.chars().count() <= 12000);
            let parsed_brief: Value = serde_json::from_str(brief_text).unwrap();
            assert_eq!(parsed_brief["id"], "claim:EM.FP.01");
            assert!(parsed_brief["notice"]
                .as_str()
                .unwrap()
                .contains("Deterministic"));
            let context=handle(&store,json!({"jsonrpc":"2.0","id":8,"method":"tools/call","params":{"name":"research_context","arguments":{"id":"claim:EM.FP.01","max_chars":12000}}})).await.unwrap();
            assert_eq!(context["result"]["isError"], false);
            assert!(context["result"].get("structuredContent").is_none());
            let context_text = context["result"]["content"][0]["text"].as_str().unwrap();
            assert!(context_text.chars().count() <= 12000);
            assert_eq!(
                serde_json::from_str::<Value>(context_text).unwrap()["id"],
                "claim:EM.FP.01"
            );
            let complete_response = serde_json::to_string(&brief).unwrap();
            assert!(
                complete_response.len() < brief_text.len() * 2,
                "complete MCP response duplicated too much content"
            );
            let comparison=handle(&store,json!({"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"compare_aspects","arguments":{"ids":["claim:EM.FP.01","claim:SECOND"],"max_chars":12000}}})).await.unwrap();
            assert!(comparison["result"].get("structuredContent").is_none());
            let compare_text = comparison["result"]["content"][0]["text"].as_str().unwrap();
            let parsed_compare: Value = serde_json::from_str(compare_text).unwrap();
            assert_eq!(parsed_compare["items"].as_array().unwrap().len(), 2);
            let theory=handle(&store,json!({"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"theory_overview","arguments":{}}})).await.unwrap();
            assert!(theory["result"].get("structuredContent").is_none());
            let theory_text = theory["result"]["content"][0]["text"].as_str().unwrap();
            let parsed_theory: Value = serde_json::from_str(theory_text).unwrap();
            assert_eq!(parsed_theory["available"], false);
            std::fs::remove_dir_all(root).unwrap();
        }
        #[tokio::test]
        async fn theory_map_http_endpoint_returns_full_authored_map() {
            let root = std::env::temp_dir().join(format!(
                "semantic-explorer-theory-map-http-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            let explorer = root.join("semantic_explorer");
            let data_dir = explorer.join("data");
            let web = explorer.join("web");
            std::fs::create_dir_all(&data_dir).unwrap();
            std::fs::create_dir_all(&web).unwrap();
            let nodes: Vec<_> = (0..28)
                .map(|i| {
                    json!({"id":format!("map:n{i}"),"label_de":format!("Knoten {i}"),"label_en":format!("Node {i}"),"source_refs":[{"source_id":"source:ledger","path":"verification/status_ledger.csv","line_start":10+i,"line_end":12+i,"excerpt":"Original excerpt"}]})
                })
                .collect();
            let edges: Vec<_> = (0..29)
                .map(|i| {
                    json!({"from":format!("map:n{}",i%28),"to":format!("map:n{}",(i+1)%28),"kind":"conditional","label_de":"bedingt","label_en":"conditional","source_refs":[{"source_id":"source:ledger","path":"verification/status_ledger.csv","line_start":20+i,"line_end":20+i,"excerpt":"Edge quote"}]})
                })
                .collect();
            let theory_map = json!({"title_de":"Gesamtkarte","title_en":"Theory map","caption_de":"Quelle bleibt sichtbar","caption_en":"Sources stay visible","nodes":nodes,"edges":edges,"lanes":[{"id":"lane:one","label_de":"Linie","label_en":"Lane"}],"stages":[{"id":"stage:one","label_de":"Stufe","label_en":"Stage"}],"states":[{"id":"open","label_de":"Offen","label_en":"Open"}],"patterns":[{"id":"pattern:one","label_de":"Muster","label_en":"Pattern","summary_de":"Bedingte Brücke","summary_en":"Conditional bridge","source_refs":[{"source_id":"source:ledger","path":"verification/status_ledger.csv","line_start":4,"line_end":5,"excerpt":"Pattern source quote"}]}]});
            std::fs::write(
                explorer.join("theory_map.json"),
                serde_json::to_vec(&theory_map).unwrap(),
            )
            .unwrap();
            let data = data_dir.join("knowledge.json");
            std::fs::write(
                &data,
                r#"{"meta":{"schema_version":1},"nodes":[],"edges":[],"sources":[]}"#,
            )
            .unwrap();
            let store = Arc::new(GraphStore::new(data, root.clone()).unwrap());
            let app = app_router(AppState { store }, web.clone(), web.join("index.html"));
            let response = app
                .oneshot(
                    Request::builder()
                        .uri("/api/theory-map")
                        .body(Body::empty())
                        .unwrap(),
                )
                .await
                .unwrap();
            assert_eq!(response.status(), StatusCode::OK);
            let payload: Value =
                serde_json::from_slice(&to_bytes(response.into_body(), 2_000_000).await.unwrap())
                    .unwrap();
            assert_eq!(payload["nodes"].as_array().unwrap().len(), 28);
            assert_eq!(payload["edges"].as_array().unwrap().len(), 29);
            assert_eq!(payload["nodes"][0]["label_en"], "Node 0");
            assert_eq!(
                payload["nodes"][0]["source_refs"][0]["path"],
                "verification/status_ledger.csv"
            );
            assert_eq!(payload["nodes"][0]["source_refs"][0]["line_start"], 10);
            assert_eq!(
                payload["patterns"][0]["source_refs"][0]["excerpt"],
                "Pattern source quote"
            );
            std::fs::remove_dir_all(root).unwrap();
        }

        #[tokio::test]
        async fn invalid_brief_budget_and_chain_depth_are_http_400() {
            let root = std::env::temp_dir().join(format!(
                "semantic-explorer-http-bounds-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            std::fs::create_dir_all(&root).unwrap();
            let data = root.join("knowledge.json");
            std::fs::write(
                &data,
                r#"{"meta":{"schema_version":1},"nodes":[],"edges":[],"sources":[]}"#,
            )
            .unwrap();
            let store = Arc::new(GraphStore::new(data, root.clone()).unwrap());
            let app = app_router(AppState { store }, root.clone(), root.join("index.html"));
            for uri in [
                "/api/brief?id=topic%3Aalpha&max_chars=0",
                "/api/context?id=topic%3Aalpha&max_chars=0",
                "/api/chain?id=claim%3AAX.P1.01&depth=999",
            ] {
                let response = app
                    .clone()
                    .oneshot(Request::builder().uri(uri).body(Body::empty()).unwrap())
                    .await
                    .unwrap();
                assert_eq!(response.status(), StatusCode::BAD_REQUEST, "{uri}");
            }
            std::fs::remove_dir_all(root).unwrap();
        }
    }
}
