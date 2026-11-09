use super::*;
use serde_json::json;
use std::{
    fs,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempDir(PathBuf);

impl TempDir {
    fn new(label: &str) -> Self {
        let nonce = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "tfpt-research-{label}-{}-{nanos}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&path).expect("create temporary research directory");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn installed_program() -> Option<Program> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../research/program.json");
    path.exists().then(|| serde_json::from_slice(&fs::read(path).unwrap()).unwrap())
}

fn text(value: &str) -> Text {
    Text {
        de: value.to_owned(),
        en: value.to_owned(),
    }
}

fn premise(id: &str, statement: &str) -> Premise {
    Premise {
        id: id.to_owned(),
        title: text(statement),
        statement: text(statement),
    }
}

fn citation(root: &Path, relative: &str, contents: &str) -> Citation {
    let path = root.join(relative);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("create source parent");
    }
    fs::write(&path, contents).expect("write pinned source");
    Citation {
        path: relative.to_owned(),
        sha256: digest(contents.as_bytes()),
        line_start: 1,
        line_end: contents.lines().count().max(1),
    }
}

fn finding(id: &str, status: &str, premise_ids: &[&str], source: Citation) -> Finding {
    Finding {
        id: id.to_owned(),
        title: text(id),
        statement: text("Reviewed result"),
        status: status.to_owned(),
        premise_ids: premise_ids.iter().map(|id| (*id).to_owned()).collect(),
        scope: text("Only under the explicitly named premises"),
        sources: vec![source],
        attempt_ids: vec![],
    }
}

fn task(id: &str, status: &str, resolutions: &[&str], required_premises: &[&str]) -> Task {
    Task {
        kind: "construction".to_owned(),
        id: id.to_owned(),
        title: text(id),
        status: status.to_owned(),
        depends_on: vec![],
        resolution_ids: resolutions.iter().map(|id| (*id).to_owned()).collect(),
        required_premise_ids: required_premises
            .iter()
            .map(|id| (*id).to_owned())
            .collect(),
        success_criterion: text("The scoped decision is supported"),
        next_action: text("Continue with the next unresolved dependency"),
    }
}

#[test]
fn negative_finding_can_resolve_assessment_but_not_construction() {
    let temp = TempDir::new("task-kind");
    let mut program = fixture(temp.path());
    let task = &mut program.problems[0].tasks[0];
    task.status = "resolved".to_owned();

    let error = validate(&program)
        .expect_err("a scoped exclusion cannot resolve a construction task")
        .to_string();
    assert!(error.contains("exclusion cannot resolve a construction task"));

    program.problems[0].tasks[0].kind = "assessment".to_owned();
    validate(&program).expect("the same exclusion may decide an assessment task");
}

fn fixture(root: &Path) -> Program {
    let pin = citation(root, "sources/result.txt", "reviewed source\nsecond line\n");
    let parallel = premise("premise:parallel", "The two pairs are parallel");
    let common = premise("premise:common", "The common source convention is retained");
    let excluded = finding(
        "finding:parallel-excluded",
        "excluded_scoped",
        &["premise:parallel"],
        pin.clone(),
    );
    let established = finding(
        "finding:known-scoped",
        "established_scoped",
        &["premise:common"],
        pin,
    );
    let problem = Problem {
        id: "problem:operator-map".to_owned(),
        title: text("Operator map"),
        question: text("Which operator map is supported?"),
        topics: vec!["claim:topic".to_owned()],
        goals: vec!["claim:goal".to_owned()],
        accepted_premise_ids: vec![],
        findings: vec![excluded, established],
        tasks: vec![
            task(
                "task:parallel-route",
                "excluded",
                &["finding:parallel-excluded"],
                &["premise:parallel"],
            ),
            task(
                "task:known-route",
                "resolved",
                &["finding:known-scoped"],
                &["premise:common"],
            ),
        ],
        strategies: vec![
            Strategy {
                id: "strategy:parallel-two-plus-two".to_owned(),
                title: text("Parallel two plus two"),
                aliases: vec!["parallel split".to_owned()],
                outcome: "excluded_scoped".to_owned(),
                premise_ids: vec!["premise:parallel".to_owned()],
                finding_ids: vec!["finding:parallel-excluded".to_owned()],
                explanation: text("Excluded only when the pairs are parallel"),
                changes_required: text("Supply a different, explicit premise"),
            },
            Strategy {
                id: "strategy:reuse-known".to_owned(),
                title: text("Reuse the known scoped result"),
                aliases: vec![],
                outcome: "reuse_result".to_owned(),
                premise_ids: vec!["premise:common".to_owned()],
                finding_ids: vec!["finding:known-scoped".to_owned()],
                explanation: text("Reuse only with the common-source premise"),
                changes_required: text("State that premise explicitly"),
            },
        ],
        reviewed_attempt_ids: vec![],
        connections: vec![],
    };
    Program {
        schema_version: 1,
        revision: "fixture-revision".to_owned(),
        premises: vec![parallel, common],
        problems: vec![problem],
    }
}

fn graph() -> KnowledgeGraph {
    KnowledgeGraph {
        nodes: vec![
            Node {
                id: "claim:topic".to_owned(),
                kind: "claim".to_owned(),
                ..Default::default()
            },
            Node {
                id: "claim:goal".to_owned(),
                kind: "claim".to_owned(),
                ..Default::default()
            },
        ],
        ..Default::default()
    }
}

fn add_connected_problem(program: &mut Program) {
    program.problems.push(Problem {
        id: "problem:target".to_owned(),
        title: text("Target problem"),
        question: text("What does the source problem constrain?"),
        topics: vec!["claim:topic".to_owned()],
        goals: vec!["claim:goal".to_owned()],
        accepted_premise_ids: vec![],
        findings: vec![],
        tasks: vec![],
        strategies: vec![],
        reviewed_attempt_ids: vec![],
        connections: vec![],
    });
    program.problems[0].connections.push(Connection {
        id: "connection:operator-target".to_owned(),
        target_problem_id: "problem:target".to_owned(),
        kind: "constrains_under".to_owned(),
        premise_ids: vec!["premise:parallel".to_owned()],
        finding_ids: vec!["finding:parallel-excluded".to_owned()],
        explanation: text("The exclusion constrains the target only under parallelity"),
    });
}

#[test]
fn editable_record_scopes_premises_without_truncating_the_problem() {
    let temp = TempDir::new("editable-scope");
    let mut program = fixture(temp.path());
    add_connected_problem(&mut program);
    for (id, statement) in [
        ("premise:accepted-only", "Accepted only"),
        ("premise:task-only", "Task only"),
        ("premise:strategy-only", "Strategy only"),
        ("premise:connection-only", "Connection only"),
    ] {
        program.premises.push(premise(id, statement));
    }
    program
        .premises
        .push(premise("premise:unrelated", &"x".repeat(60000)));
    let q = &mut program.problems[0];
    q.accepted_premise_ids.push("premise:accepted-only".into());
    q.tasks[0]
        .required_premise_ids
        .push("premise:task-only".into());
    q.strategies[0]
        .premise_ids
        .push("premise:strategy-only".into());
    q.connections[0]
        .premise_ids
        .push("premise:connection-only".into());
    let record = editable_record(&program, "problem:operator-map", 60000).unwrap();
    assert_eq!(
        record["problem"],
        serde_json::to_value(&program.problems[0]).unwrap()
    );
    let ids: HashSet<_> = record["premises"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x["id"].as_str().unwrap())
        .collect();
    assert_eq!(
        ids,
        HashSet::from([
            "premise:parallel",
            "premise:common",
            "premise:accepted-only",
            "premise:task-only",
            "premise:strategy-only",
            "premise:connection-only",
        ])
    );
    let current = freshness(&program, temp.path());
    let view = problem_view(&program.problems[0], &program, &graph(), "en", &current);
    assert!(view["premises"]
        .as_array()
        .unwrap()
        .iter()
        .any(|x| x["id"] == "premise:connection-only"));
    program.problems[0].question = text(&"x".repeat(60000));
    assert!(editable_record(&program, "problem:operator-map", 60000).is_err());
}

#[test]
fn scoped_decisions_require_exact_strategy_and_explicit_premises() {
    let temp = TempDir::new("premises");
    let program = fixture(temp.path());
    validate(&program).expect("fixture must validate");

    let free_text = assess(
        &program,
        temp.path(),
        "problem:operator-map",
        "parallel two plus two",
        &[],
        "en",
    )
    .expect("free-text assessment");
    assert_eq!(free_text["decision"], "unclassified");
    assert!(free_text["candidates"]
        .as_array()
        .expect("candidate list")
        .iter()
        .any(|candidate| candidate["id"] == "strategy:parallel-two-plus-two"));

    let without_parallel = assess(
        &program,
        temp.path(),
        "problem:operator-map",
        "strategy:parallel-two-plus-two",
        &[],
        "en",
    )
    .expect("exact assessment without premise");
    assert_eq!(without_parallel["decision"], "scope_not_established");
    assert_eq!(
        without_parallel["missing_premise_ids"],
        json!(["premise:parallel"])
    );

    let with_parallel = assess(
        &program,
        temp.path(),
        "problem:operator-map",
        "strategy:parallel-two-plus-two",
        &["premise:parallel".to_owned()],
        "en",
    )
    .expect("exact assessment with premise");
    assert_eq!(with_parallel["decision"], "excluded_scoped");

    let reuse_without_scope = assess(
        &program,
        temp.path(),
        "problem:operator-map",
        "strategy:reuse-known",
        &[],
        "en",
    )
    .expect("known result without premise");
    assert_eq!(reuse_without_scope["decision"], "scope_not_established");
}

#[test]
fn resolved_task_cannot_drop_a_resolution_finding_premise() {
    let temp = TempDir::new("task-premise");
    let mut program = fixture(temp.path());
    let task = program.problems[0]
        .tasks
        .iter_mut()
        .find(|task| task.id == "task:known-route")
        .expect("known task");
    task.required_premise_ids.clear();

    let error = validate(&program).expect_err("resolved task must retain finding scope");
    assert!(
        error
            .to_string()
            .contains("task task:known-route omits resolution premise premise:common"),
        "unexpected scope-loss diagnostic: {error:#}"
    );
}

#[test]
fn stale_source_revokes_resolved_and_reuse_decisions() {
    let temp = TempDir::new("stale");
    let program = fixture(temp.path());
    validate(&program).expect("fixture must validate");

    let current = state(
        &program,
        &graph(),
        temp.path(),
        "problem:operator-map",
        "en",
        30_000,
        0,
        0,
    )
    .expect("current state");
    let known = current["problems"][0]["tasks"]
        .as_array()
        .expect("tasks")
        .iter()
        .find(|task| task["id"] == "task:known-route")
        .expect("known task");
    assert_eq!(known["effective_status"], "resolved");

    fs::write(temp.path().join("sources/result.txt"), "changed source\n").expect("make pin stale");
    let stale = state(
        &program,
        &graph(),
        temp.path(),
        "problem:operator-map",
        "en",
        30_000,
        0,
        0,
    )
    .expect("stale state");
    let known = stale["problems"][0]["tasks"]
        .as_array()
        .expect("tasks")
        .iter()
        .find(|task| task["id"] == "task:known-route")
        .expect("known task");
    assert_eq!(known["effective_status"], "needs_review");

    let decision = assess(
        &program,
        temp.path(),
        "problem:operator-map",
        "strategy:reuse-known",
        &["premise:common".to_owned()],
        "en",
    )
    .expect("stale assessment");
    assert_eq!(decision["decision"], "needs_review");
}

#[test]
fn dangling_and_cyclic_task_dependencies_are_rejected() {
    let temp = TempDir::new("dependencies");
    let program = fixture(temp.path());

    let mut dangling = program.clone();
    dangling.problems[0].tasks[0].depends_on = vec!["task:missing".to_owned()];
    let error = validate(&dangling).expect_err("dangling dependency must fail");
    assert!(error
        .to_string()
        .contains("unknown prerequisite task task:missing"));

    let mut cycle = program;
    cycle.problems[0].tasks[0].depends_on = vec!["task:known-route".to_owned()];
    cycle.problems[0].tasks[1].depends_on = vec!["task:parallel-route".to_owned()];
    let error = validate(&cycle).expect_err("dependency cycle must fail");
    assert!(error.to_string().contains("cyclic task dependency"));
}

#[test]
fn connections_require_existing_targets_and_preserve_finding_scope() {
    let temp = TempDir::new("connections");
    let mut program = fixture(temp.path());
    add_connected_problem(&mut program);
    validate(&program).expect("scoped connection must validate");

    let mut missing_scope = program.clone();
    missing_scope.problems[0].connections[0].premise_ids.clear();
    let error = validate(&missing_scope).expect_err("connection cannot drop finding premise");
    assert!(error
        .to_string()
        .contains("connection connection:operator-target omits finding premise premise:parallel"));

    let mut dangling = program;
    dangling.problems[0].connections[0].target_problem_id = "problem:missing".to_owned();
    let error = validate(&dangling).expect_err("connection target must exist");
    assert!(error
        .to_string()
        .contains("unknown or self research connection target"));
}

#[test]
fn oversized_problem_returns_scoped_partial_page_and_continuation() {
    let temp = TempDir::new("budget");
    let mut program = fixture(temp.path());
    let expected_scope = program.problems[0].findings[0].scope.en.clone();
    let template = program.problems[0].findings[0].clone();
    for index in 0..20 {
        let mut extra = template.clone();
        extra.id = format!("finding:budget-extra-{index}");
        extra.title = text(&format!("Additional reviewed result {index}"));
        program.problems[0].findings.push(extra);
    }
    validate(&program).expect("long fixture must validate");

    let compact = state(
        &program,
        &graph(),
        temp.path(),
        "problem:operator-map",
        "en",
        12_000,
        0,
        0,
    )
    .expect("compact state");
    assert_eq!(compact["problems"].as_array().expect("problems").len(), 1);
    let page = &compact["problems"][0];
    assert_eq!(page["packet_state"], "partial");
    assert_eq!(page["partial_coverage"]["complete"], false);
    assert_eq!(
        page["partial_coverage"]["returned_record_content_complete"],
        true
    );
    assert!(page["findings"].as_array().unwrap().iter().any(|finding| {
        finding["scope"] == expected_scope
            && finding["premise_ids"] == json!(["premise:parallel"])
            && !finding["sources"].as_array().unwrap().is_empty()
    }));
    assert!(page["premises"]
        .as_array()
        .unwrap()
        .iter()
        .any(|premise| { premise["id"] == "premise:parallel" }));
    assert!(!page["dependency_handles"].as_array().unwrap().is_empty());
    assert_eq!(compact["budget"]["truncated"], true);
    assert_eq!(compact["budget"]["next_query"], "problem:operator-map");
    let next = compact["budget"]["next_detail_offset"].as_u64().unwrap() as usize;
    assert!(next > 0);
    let continuation = state(
        &program,
        &graph(),
        temp.path(),
        "problem:operator-map",
        "en",
        12_000,
        0,
        next,
    )
    .expect("continued compact state");
    assert_eq!(continuation["problems"][0]["packet_state"], "partial");
    assert_eq!(
        continuation["problems"][0]["partial_coverage"]["detail_offset"],
        next
    );

    let complete = state(
        &program,
        &graph(),
        temp.path(),
        "problem:operator-map",
        "en",
        60_000,
        0,
        0,
    )
    .expect("complete state");
    assert_eq!(complete["problems"].as_array().expect("problems").len(), 1);
    assert_eq!(
        complete["problems"][0]["findings"][0]["scope"],
        expected_scope
    );
    assert_eq!(
        complete["problems"][0]["findings"][0]["premise_ids"],
        json!(["premise:parallel"])
    );
}

#[test]
fn update_persists_and_fresh_load_observes_new_revision() {
    let temp = TempDir::new("persistence");
    let root = temp.path().join("repository");
    fs::create_dir_all(&root).expect("create repository root");
    let program = fixture(&root);
    let path = temp.path().join("state/program.json");
    fs::create_dir_all(path.parent().expect("state parent")).expect("create state directory");
    let initial = Program {
        schema_version: 1,
        revision: "empty-revision".to_owned(),
        premises: vec![],
        problems: vec![],
    };
    fs::write(&path, serde_json::to_vec_pretty(&initial).unwrap()).expect("write initial state");
    let initial_revision = load(&path).expect("load content-derived revision").revision;

    let result = update(
        &path,
        &root,
        &graph(),
        &initial_revision,
        serde_json::to_value(&program.problems[0]).unwrap(),
        program.premises.clone(),
        "Record the reviewed scoped decisions",
    )
    .expect("persist update");
    let revision = result["revision"].as_str().expect("new revision");
    assert_ne!(revision, "empty-revision");

    let loaded = load(&path).expect("fresh process-style load");
    assert_eq!(loaded.revision, revision);
    assert_eq!(loaded.problems.len(), 1);
    assert_eq!(loaded.problems[0].findings.len(), 2);

    let conflict = update(
        &path,
        &root,
        &graph(),
        &initial_revision,
        serde_json::to_value(&program.problems[0]).unwrap(),
        vec![],
        "Stale writer",
    )
    .expect_err("stale compare-and-swap must fail");
    assert!(conflict.to_string().contains("revision conflict"));
}

#[test]
fn update_cannot_silently_remove_strategy_or_connection_decisions() {
    let temp = TempDir::new("decision-retention");
    let root = temp.path().join("repository");
    fs::create_dir_all(&root).expect("create repository root");
    let mut program = fixture(&root);
    add_connected_problem(&mut program);
    let path = temp.path().join("state/program.json");
    fs::create_dir_all(path.parent().expect("state parent")).expect("create state directory");
    fs::write(&path, serde_json::to_vec_pretty(&program).unwrap()).expect("write program");
    let revision = load(&path).expect("load revision").revision;

    let mut without_strategy = program.problems[0].clone();
    without_strategy.strategies.remove(0);
    let error = update(
        &path,
        &root,
        &graph(),
        &revision,
        serde_json::to_value(without_strategy).unwrap(),
        vec![],
        "Accidental strategy omission",
    )
    .expect_err("strategy decision must be retained");
    assert!(error
        .to_string()
        .contains("cannot silently remove strategy strategy:parallel-two-plus-two"));

    let mut without_connection = program.problems[0].clone();
    without_connection.connections.clear();
    let error = update(
        &path,
        &root,
        &graph(),
        &revision,
        serde_json::to_value(without_connection).unwrap(),
        vec![],
        "Accidental connection omission",
    )
    .expect_err("connection decision must be retained");
    assert!(error
        .to_string()
        .contains("cannot silently remove connection connection:operator-target"));
}

#[test]
fn content_edit_with_unchanged_embedded_revision_invalidates_cas() {
    let temp = TempDir::new("content-cas");
    let root = temp.path().join("repository");
    fs::create_dir_all(&root).expect("create repository root");
    let mut program = fixture(&root);
    let path = temp.path().join("state/program.json");
    fs::create_dir_all(path.parent().expect("state parent")).expect("create state directory");
    fs::write(&path, serde_json::to_vec_pretty(&program).unwrap()).expect("write program");
    let prior_revision = load(&path).expect("load original revision").revision;

    program.problems[0].title = text("Manually edited title");
    fs::write(&path, serde_json::to_vec_pretty(&program).unwrap()).expect("manual edit");
    let edited_revision = load(&path).expect("load edited revision").revision;
    assert_ne!(edited_revision, prior_revision);

    let error = update(
        &path,
        &root,
        &graph(),
        &prior_revision,
        serde_json::to_value(&program.problems[0]).unwrap(),
        vec![],
        "Stale writer after manual edit",
    )
    .expect_err("content-derived CAS must detect manual edit");
    assert!(error.to_string().contains("revision conflict"));
}

#[test]
fn update_validates_traversal_before_attempting_source_freshness() {
    let temp = TempDir::new("traversal");
    let root = temp.path().join("repository");
    fs::create_dir_all(&root).expect("create repository root");
    let outside = temp.path().join("outside.txt");
    fs::write(&outside, "outside\n").expect("write outside sentinel");

    let path = temp.path().join("state/program.json");
    fs::create_dir_all(path.parent().expect("state parent")).expect("create state directory");
    let initial = Program {
        schema_version: 1,
        revision: "empty-revision".to_owned(),
        premises: vec![],
        problems: vec![],
    };
    fs::write(&path, serde_json::to_vec_pretty(&initial).unwrap()).expect("write initial state");
    let initial_revision = load(&path).expect("load content-derived revision").revision;

    let bad_finding = finding(
        "finding:bad-path",
        "established_scoped",
        &[],
        Citation {
            path: "../outside.txt".to_owned(),
            sha256: digest(b"outside\n"),
            line_start: 1,
            line_end: 1,
        },
    );
    let bad_problem = Problem {
        id: "problem:bad-path".to_owned(),
        title: text("Bad path"),
        question: text("Must not escape the repository"),
        topics: vec!["claim:topic".to_owned()],
        goals: vec!["claim:goal".to_owned()],
        accepted_premise_ids: vec![],
        findings: vec![bad_finding],
        tasks: vec![],
        strategies: vec![],
        reviewed_attempt_ids: vec![],
        connections: vec![],
    };

    let error = update(
        &path,
        &root,
        &graph(),
        &initial_revision,
        serde_json::to_value(bad_problem).unwrap(),
        vec![],
        "Traversal regression",
    )
    .expect_err("traversal must be rejected");
    assert!(
        error
            .to_string()
            .contains("source must be repository-relative without traversal"),
        "schema/path validation must run before source freshness: {error:#}"
    );
}

#[test]
fn actual_program_handles_natural_questions_and_preserves_connection_endpoints() {
    let Some(p) = installed_program() else { return; };
    validate(&p).expect("authored research records satisfy scope and dependency contracts");
    let temp = TempDir::new("actual-query");
    let mut g = KnowledgeGraph::default();
    for i in 0..14 {
        g.nodes.push(Node {
            id: format!("topic:fixture-{i}"),
            kind: "topic".into(),
            label: format!("Long topic title for navigation {i}"),
            ..Default::default()
        });
    }
    for i in 0..8 {
        let id = format!("attempt:pending-{i}");
        g.nodes.push(Node {
            id: id.clone(),
            kind: "attempt".into(),
            label: format!(
                "Long historical attempt title {}",
                "scope and original-source review ".repeat(8)
            ),
            ..Default::default()
        });
        g.edges.push(Edge {
            src: id,
            dst: "topic:universalraum".into(),
            kind: "theme_links".into(),
            ..Default::default()
        });
    }
    let answer=state(&p,&g,temp.path(),"Was ist bereits darüber bekannt, Gedächtnis und Zustand aus zukünftigem Verhalten herzuleiten, und was muss als nächster Forschungsschritt tatsächlich noch gezeigt werden?","de",12000,0,0).unwrap();
    assert_eq!(
        answer["problem_index"][0]["id"],
        "problem:future-quotient-and-gns-reconstruction"
    );
    assert!(serde_json::to_string(&answer).unwrap().chars().count() <= 12000);
    let complete = if answer["problems"].as_array().unwrap().is_empty() {
        assert_eq!(
            answer["budget"]["next_query"],
            "problem:future-quotient-and-gns-reconstruction"
        );
        let needed = answer["budget"]["required_chars_for_next_problem"]
            .as_u64()
            .unwrap() as usize;
        assert!(needed > 12000 && needed <= 60000);
        state(
            &p,
            &g,
            temp.path(),
            answer["budget"]["next_query"].as_str().unwrap(),
            "de",
            needed,
            0,
            0,
        )
        .unwrap()
    } else {
        answer
    };
    assert_eq!(
        complete["problems"][0]["id"],
        "problem:future-quotient-and-gns-reconstruction"
    );
    let answer = state(
        &p,
        &g,
        temp.path(),
        "problem:common-source-state-event-selection",
        "en",
        60000,
        0,
        0,
    )
    .unwrap();
    assert!(!answer["problems"].as_array().unwrap().is_empty());
    assert!(serde_json::to_string(&answer).unwrap().chars().count() <= 60000);
    if answer["problems"][0]["packet_state"] == "partial" {
        assert_eq!(answer["budget"]["partial_problem"], true);
        assert_eq!(answer["problems"][0]["partial_coverage"]["complete"], false);
        assert!(
            answer["problems"][0]["partial_coverage"]["detail_returned"]
                .as_u64()
                .unwrap()
                > 0
        );
        return;
    }
    let connections = answer["problems"][0]["connections"].as_array().unwrap();
    assert!(!connections.is_empty());
    for c in connections {
        assert_ne!(c["source_problem_id"], c["target_problem_id"]);
        if c["direction"] == "incoming" {
            assert_eq!(
                c["target_problem_id"],
                "problem:common-source-state-event-selection"
            );
            assert_eq!(c["other_problem_id"], c["source_problem_id"]);
        } else {
            assert_eq!(c["direction"], "outgoing");
            assert_eq!(
                c["source_problem_id"],
                "problem:common-source-state-event-selection"
            );
            assert_eq!(c["other_problem_id"], c["target_problem_id"]);
        }
    }
}

#[test]
fn actual_operator_problem_remains_editable_with_its_own_premises() {
    let Some(p) = installed_program() else { return; };
    match editable_record(&p, "problem:automated-operator-source-synthesis", 60000) {
        Ok(record) => {
            assert!(!record["problem"]["findings"].as_array().unwrap().is_empty());
            assert!(record["premises"].as_array().unwrap().len() < p.premises.len());
            assert!(serde_json::to_string(&record).unwrap().chars().count() <= 60000);
        }
        Err(error) => {
            // Oversized records remain intact; the API rejects rather than clipping them.
            assert!(error.to_string().contains("editable record needs"));
            let temp = TempDir::new("actual-operator-budget");
            let state = state(
                &p,
                &KnowledgeGraph::default(),
                temp.path(),
                "problem:automated-operator-source-synthesis",
                "de",
                60000,
                0,
                0,
            )
            .unwrap();
            assert_eq!(state["budget"]["truncated"], true);
            assert_eq!(state["budget"]["partial_problem"], true);
            assert_eq!(state["problems"][0]["packet_state"], "partial");
            assert!(state["problem_index"]
                .as_array()
                .unwrap()
                .iter()
                .any(|r| r["id"] == "problem:automated-operator-source-synthesis"));
        }
    }
}

#[test]
fn actual_large_clock_problem_returns_relevant_scoped_records_at_normal_budget() {
    let Some(p) = installed_program() else { return; };
    let temp = TempDir::new("actual-large-clock-page");
    let query = "P1 lokaler Operator W8 ursprünglicher Stress Quelle bereits gelöst";
    let twelve_k = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        query,
        "de",
        12_000,
        0,
        0,
    )
    .expect("twelve-k canonical state");
    assert_eq!(twelve_k["problems"][0]["id"], "problem:atlas-clock-defect-and-operational-reconstruction");
    assert_eq!(twelve_k["problems"][0]["packet_state"], "partial");
    assert!(
        twelve_k["problems"][0]["partial_coverage"]["detail_returned"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert!(twelve_k["problems"][0]["findings"]
        .as_array()
        .unwrap()
        .iter()
        .any(|finding| finding["id"]
            == "finding:full-original-c1-stress-cycle-reconstructs-neutral-w8-response"));
    assert!(serde_json::to_string(&twelve_k).unwrap().chars().count() <= 12_000);
    let active = twelve_k["problems"][0]["active_task_handles"]
        .as_array()
        .unwrap();
    let main_task = active
        .iter()
        .find(|task| {
            task["id"] == "task:derive-original-seam-selection-of-native-observable-process"
        })
        .expect("the active seam-selection goal must survive record pagination");
    assert_eq!(main_task["status"], "open");
    assert_eq!(
        main_task["required_premise_ids"],
        json!(["premise:declared-p1-p2"])
    );
    assert_eq!(main_task["next_query"], main_task["id"]);
    let twelve_next = &twelve_k["budget"]["continuation"];
    let twelve_continued = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        query,
        "de",
        twelve_next["max_chars"].as_u64().unwrap() as usize,
        twelve_next["offset"].as_u64().unwrap() as usize,
        twelve_next["detail_offset"].as_u64().unwrap() as usize,
    )
    .expect("twelve-k continuation raises the budget when the next record needs it");
    assert!(
        twelve_continued["problems"][0]["partial_coverage"]["detail_returned"]
            .as_u64()
            .unwrap()
            > 0
    );

    let first = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        query,
        "de",
        24_000,
        0,
        0,
    )
    .expect("bounded canonical state");
    assert_eq!(
        first["problems"][0]["id"],
        "problem:atlas-clock-defect-and-operational-reconstruction"
    );
    assert_eq!(first["problems"][0]["packet_state"], "partial");
    assert!(first["problems"][0]["findings"]
        .as_array()
        .unwrap()
        .iter()
        .any(|finding| finding["id"]
            == "finding:full-original-c1-stress-cycle-reconstructs-neutral-w8-response"));
    assert!(first["problems"][0]["active_task_handles"]
        .as_array()
        .unwrap()
        .iter()
        .any(
            |task| task["id"] == "task:derive-original-seam-selection-of-native-observable-process"
        ));
    for finding in first["problems"][0]["findings"].as_array().unwrap() {
        assert!(finding["scope"].is_string());
        assert!(finding["premise_ids"].is_array());
        assert!(!finding["sources"].as_array().unwrap().is_empty());
    }
    assert_eq!(
        first["problems"][0]["partial_coverage"]["returned_record_content_complete"],
        true
    );
    assert_eq!(first["budget"]["partial_problem"], true);
    assert!(serde_json::to_string(&first).unwrap().chars().count() <= 24_000);

    let next = first["budget"]["next_detail_offset"]
        .as_u64()
        .expect("large problem needs a continuation") as usize;
    let second = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        query,
        "de",
        24_000,
        0,
        next,
    )
    .expect("continued canonical state");
    assert_eq!(
        second["problems"][0]["partial_coverage"]["detail_offset"],
        next
    );
    let first_ids: HashSet<_> = DETAIL_FIELDS
        .iter()
        .flat_map(|field| first["problems"][0][*field].as_array().unwrap())
        .filter_map(|record| record["id"].as_str())
        .collect();
    let second_ids: HashSet<_> = DETAIL_FIELDS
        .iter()
        .flat_map(|field| second["problems"][0][*field].as_array().unwrap())
        .filter_map(|record| record["id"].as_str())
        .collect();
    assert!(first_ids.is_disjoint(&second_ids));
}

#[test]
fn uncovered_topic_points_to_existing_corpus_and_counts_attempts() {
    let temp = TempDir::new("corpus-fallback");
    let p = fixture(temp.path());
    let mut g = KnowledgeGraph::default();
    g.nodes.push(Node {
        id: "topic:cosmology".into(),
        kind: "topic".into(),
        ..Default::default()
    });
    g.nodes.push(Node {
        id: "attempt:existing".into(),
        kind: "attempt".into(),
        ..Default::default()
    });
    g.edges.push(Edge {
        src: "attempt:existing".into(),
        dst: "topic:cosmology".into(),
        kind: "in_topic".into(),
        ..Default::default()
    });
    let v = state(&p, &g, temp.path(), "topic:cosmology", "en", 4000, 0, 0).unwrap();
    assert_eq!(v["coverage"]["matched_problems"], 0);
    assert_eq!(v["coverage"]["pending_attempts"], 1);
    assert_eq!(v["corpus_fallback"]["next_tool"], "research_context");
    assert_eq!(v["corpus_fallback"]["arguments"]["id"], "topic:cosmology");
}

#[test]
fn research_detail_ranking_prefers_result_over_source_metadata() {
    let full = json!({
        "id":"problem:large",
        "findings":[
            {"id":"finding:source-rich-unrelated","title":"Neutral response",
             "statement":"Another reviewed result",
             "sources":[{"path":"D4/monodromy/vierte/Matrix/Identität/C3xC3/rigid-family.txt"}]},
            {"id":"finding:archived-d4-repair","title":"Archivierte D4 Punctureformel",
             "statement":"Die vierte Monodromie ist falsch; vorhandene Darstellung repariert sie"}
        ],
        "tasks":[],"strategies":[],"connections":[]
    });
    let records = detail_records(
        &full,
        "D4 monodromy vierte Matrix Identität C3xC3 rigid-family",
    )
    .unwrap();
    assert_eq!(records[0], ("findings", 1));
}

#[test]
fn research_detail_ranking_preserves_exact_record_id() {
    let full = json!({
        "id":"problem:large",
        "findings":[
            {"id":"finding:other","title":"finding archived d4 repair",
             "statement":"finding archived d4 repair"},
            {"id":"finding:archived-d4-repair","title":"Correct reviewed result"}
        ],
        "tasks":[],"strategies":[],"connections":[]
    });
    assert_eq!(
        detail_records(&full, "finding:archived-d4-repair").unwrap()[0],
        ("findings", 1)
    );
}

#[test]
fn research_detail_ranking_recovers_actual_d4_solution_deterministically() {
    let Some(p) = installed_program() else { return; };
    let temp = TempDir::new("actual-d4-ranking");
    let query = "D4 monodromy vierte Matrix Identität C3xC3 rigid-family";
    let g = KnowledgeGraph::default();
    let first = state(&p, &g, temp.path(), query, "de", 24_000, 0, 0).unwrap();
    let second = state(&p, &g, temp.path(), query, "de", 24_000, 0, 0).unwrap();
    assert_eq!(first, second);
    assert_eq!(
        first["problems"][0]["findings"][0]["id"],
        "finding:archived-d4-puncture-formula-fails-and-existing-v117-repairs-existence"
    );
    assert!(
        first["problems"][0]["findings"][0]["sources"]
            .as_array()
            .unwrap()
            .len()
            >= 2
    );
    assert!(serde_json::to_string(&first).unwrap().chars().count() <= 24_000);
}

#[test]
fn live_plan_precedes_historical_readiness_and_follows_progress_updates() {
    let temp = TempDir::new("live-plan");
    let mut p = fixture(temp.path());
    p.problems[0].title = text("aktuelle fundamentale Priorität");
    p.problems[0]
        .tasks
        .push(task("task:old-open", "open", &[], &[]));
    let directory = temp.path().join(PLAN_DIRECTORY);
    fs::create_dir_all(&directory).unwrap();
    let mut plan = json!({"revision":101,"updated_at":"today","next_task":"F1",
        "tasks":[{"id":"F1","title":"Current source","status":"in_progress"}],
        "result_ids":["R999"],"continuation":"Use the original sewing map.",
        "current_atlas_revision":p.revision});
    fs::write(
        directory.join("progress.json"),
        serde_json::to_vec(&plan).unwrap(),
    )
    .unwrap();
    let g = KnowledgeGraph::default();
    let first = state(
        &p,
        &g,
        temp.path(),
        "aktuelle fundamentale Priorität",
        "de",
        12_000,
        0,
        0,
    )
    .unwrap();
    assert_eq!(first["current_plan"]["next_task"], "F1");
    assert_eq!(
        first["current_plan"]["continuation"],
        "Use the original sewing map."
    );
    assert_eq!(first["current_plan"]["recent_result_ids"][0], "R999");
    assert_eq!(
        first["problem_index"][0]["execution_priority"],
        "historical_consult_current_plan"
    );
    assert!(first["problem_index"][0].get("ready_tasks").is_none());
    let historical = first["problems"][0]["tasks"]
        .as_array()
        .unwrap()
        .iter()
        .find(|task| task["id"] == "task:old-open")
        .unwrap();
    assert_eq!(historical["status"], "open");
    assert_eq!(historical["effective_status"], "ready");
    assert_eq!(
        historical["execution_priority"],
        "historical_consult_current_plan"
    );
    plan["revision"] = json!(102);
    plan["continuation"] = json!("Follow the corrected original field map.");
    plan["result_ids"] = json!(["R999", "R1000"]);
    fs::write(
        directory.join("progress.json"),
        serde_json::to_vec(&plan).unwrap(),
    )
    .unwrap();
    let second = state(
        &p,
        &g,
        temp.path(),
        "aktuelle fundamentale Priorität",
        "de",
        12_000,
        0,
        0,
    )
    .unwrap();
    assert_eq!(second["current_plan"]["revision"], 102);
    assert_eq!(second["current_plan"]["recent_result_ids"][0], "R1000");
    assert_eq!(
        second["current_plan"]["continuation"],
        "Follow the corrected original field map."
    );
    assert!(serde_json::to_string(&second).unwrap().chars().count() <= 12_000);
    let exact = state(
        &p,
        &g,
        temp.path(),
        "finding:known-scoped",
        "en",
        12_000,
        0,
        0,
    )
    .unwrap();
    assert_eq!(exact["current_plan"]["revision"], 102);
    assert!(
        serde_json::to_string(&exact["current_plan"])
            .unwrap()
            .chars()
            .count()
            <= 650
    );
}

#[test]
fn result_identifier_owns_problem_and_detail_even_in_natural_query() {
    let temp = TempDir::new("result-anchor");
    let mut p = fixture(temp.path());
    let mut current = p.problems[0].clone();
    current.id = "problem:recent-result".to_owned();
    current.title = text("Latest result");
    current.findings[0].id = "finding:r1234-original-result".to_owned();
    current.findings[0].title = text("R1234 original field map");
    p.problems[0].title = text("Was sagt R1234 zur aktuellen fundamentalen Quelle");
    p.problems[0].findings[0].statement = text(&"R1234 aktuelle fundamentale Quelle ".repeat(20));
    current.findings[1].statement =
        text(&"Was sagt R1234 zur aktuellen fundamentalen Quelle ".repeat(20));
    p.problems.push(current);
    let answer = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        "Was sagt R1234 zur aktuellen fundamentalen Quelle?",
        "de",
        12_000,
        0,
        0,
    )
    .unwrap();
    assert_eq!(answer["problems"][0]["id"], "problem:recent-result");
    let view = problem_view(
        &p.problems[1],
        &p,
        &KnowledgeGraph::default(),
        "de",
        &freshness(&p, temp.path()),
    );
    assert_eq!(
        detail_records(&view, "Was sagt R1234 zur aktuellen fundamentalen Quelle?").unwrap()[0],
        ("findings", 0)
    );
}

#[test]
fn compact_current_plan_never_hides_omitted_continuation() {
    let temp = TempDir::new("bounded-plan");
    let directory = temp.path().join(PLAN_DIRECTORY);
    fs::create_dir_all(&directory).unwrap();
    fs::write(
        directory.join("progress.json"),
        serde_json::to_vec(&json!({"revision":7,
        "next_task":"F1","tasks":[],"continuation":"ganze Quelle ".repeat(3000)}))
        .unwrap(),
    )
    .unwrap();
    let plan = current_plan(temp.path(), "revision", 650);
    assert_eq!(plan["content_complete"], false);
    assert_eq!(plan["continuation_complete"], false);
    assert_eq!(
        plan["next_request"]["arguments"]["source_id"],
        "source:fundamental-plan:progress.json"
    );
    assert!(serde_json::to_string(&plan).unwrap().chars().count() <= 650);
}

#[test]
fn explicit_priority_question_uses_live_plan_result_owner_without_overriding_exact_ids() {
    let temp = TempDir::new("priority-owner");
    let mut p = fixture(temp.path());
    p.problems[0].title =
        text("Was ist aktuelle fundamentale Forschungspriorität und nächster Schritt");
    let mut recent = p.problems[0].clone();
    recent.id = "problem:current-original-source".to_owned();
    recent.title = text("Original source synthesis");
    recent.findings[0].id = "finding:r1234-original-result".to_owned();
    recent.findings[0].title = text("R1234 typed source result");
    p.problems.push(recent);
    let directory = temp.path().join(PLAN_DIRECTORY);
    fs::create_dir_all(&directory).unwrap();
    fs::write(
        directory.join("progress.json"),
        serde_json::to_vec(&json!({
        "revision":1,"next_task":"F1","tasks":[],"continuation":"Continue the original source.",
        "result_ids":["R1234"]}))
        .unwrap(),
    )
    .unwrap();
    let answer = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        "Was ist aktuelle fundamentale Forschungspriorität und nächster Schritt?",
        "de",
        12_000,
        0,
        0,
    )
    .unwrap();
    assert_eq!(
        answer["problems"][0]["id"],
        "problem:current-original-source"
    );
    assert_eq!(
        answer["current_plan"]["continuation"],
        "Continue the original source."
    );
    let exact = state(
        &p,
        &KnowledgeGraph::default(),
        temp.path(),
        "problem:operator-map",
        "de",
        12_000,
        0,
        0,
    )
    .unwrap();
    assert_eq!(exact["problems"][0]["id"], "problem:operator-map");
    assert_eq!(exact["coverage"]["matched_problems"], 1);
}


#[test]
fn research_problem_ranking_returns_alpha_provenance_for_natural_question() {
    let Some(p) = installed_program() else { return; };
    let temp = TempDir::new("actual-alpha-provenance");
    for query in [
        "Welche Faktoren der Alpha Herleitung sind schon hergeleitet und was fehlt noch?",
        "alpha coefficient primitive 5/4 double cover e^-2alpha 48 c3",
        "Woher kommen die Korrekturfaktoren der Feinstrukturkonstante?",
    ] {
        let v = state(&p, &KnowledgeGraph::default(), temp.path(), query, "de", 24_000, 0, 0).unwrap();
        assert_eq!(v["problems"][0]["id"], "problem:alpha-fixed-input-root", "{query}");
        assert!(v["problems"][0]["findings"].as_array().unwrap().iter().any(|f|
            f["id"] == "finding:alpha-archive-term-provenance-and-master-functional-scope"), "{query}");
    }
}


#[test]
fn overall_status_preserves_evidence_and_pages_the_index() {
    let temp = TempDir::new("portfolio");
    let mut program = Program {
        revision: "overview-fixture".into(),
        problems: (0..34).map(|i| serde_json::from_value(json!({
            "id": format!("problem:tfpt-{i:02}"),
            "title": {"de": format!("TFPT Forschungsfrage {i}"), "en":format!("TFPT research question {i}")},
            "question": {"de": "TFPT", "en": "TFPT"},
            "topics": [], "goals": [], "accepted_premise_ids": [], "findings": [],
            "tasks": [], "strategies": [], "reviewed_attempt_ids": []
        })).unwrap()).collect(),
        ..serde_json::from_value(json!({"schema_version":1, "revision":"base", "premises":[], "problems":[]})).unwrap()
    };
    let graph = KnowledgeGraph { meta: json!({"coverage":{"scripts_total":715}}), ..Default::default() };
    for (query, lang) in [
        ("TFPT Gesamtstatus ursprüngliche Quelle Herleitungen Verifikationen Gravitation Raumzeit Quantenfeldtheorie Alpha Flavor Clocks", "de"),
        ("What is the overall status of TFPT?", "en"),
        ("TFPT: wo stehen wir insgesamt?", "de"),
        ("", "en"),
    ] {
        let answer = state(&program, &graph, temp.path(), query, lang, 12_000, 0, 0).unwrap();
        assert_eq!(answer["view"], "theory_overview");
        assert_eq!(answer["coverage"]["indexed_verification_scripts"], 715);
        assert_eq!(answer["problems"], json!([]));
        assert_eq!(answer["problem_index"].as_array().unwrap().len(), 34);
        for id in ["original-source", "matter-alpha", "gravity", "verification", "completion"] {
            let area = answer["evidence_areas"].as_array().unwrap().iter().find(|a| a["id"] == id).unwrap();
            assert!(area["result"].is_string() && area["scope"].is_string() && area["evidence"].is_string());
        }
    }
    let small = state(&program, &graph, temp.path(), "", "de", 2_000, 0, 0).unwrap();
    assert_eq!(small["budget"]["index_total"], 34);
    assert_eq!(small["evidence_areas"].as_array().unwrap().len(), 7);
    let next = &small["budget"]["continuation"];
    let continued = state(&program, &graph, temp.path(), "", "de",
        next["max_chars"].as_u64().unwrap() as usize, next["offset"].as_u64().unwrap() as usize, 0).unwrap();
    assert!(continued["budget"]["next_offset"].as_u64().unwrap()
        > small["budget"]["next_offset"].as_u64().unwrap());
    program.problems.truncate(1);
    let detail = state(&program, &graph, temp.path(), "problem:tfpt-00", "de", 12_000, 0, 0).unwrap();
    assert_eq!(detail["problems"][0]["id"], "problem:tfpt-00");
    assert!(detail.get("evidence_areas").is_none());
}
