//! Compact retrieval of authored relations and scoped research records.
use crate::research::{Finding, Problem, Program};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

fn array(value: &Value) -> impl Iterator<Item = &Value> {
    value.as_array().into_iter().flatten()
}

fn strings(value: &Value) -> impl Iterator<Item = &str> {
    array(value).filter_map(Value::as_str)
}

fn localized<'a>(value: &'a Value, de: &str, en: &str, lang: &str) -> &'a Value {
    if lang == "en" && !value[en].is_null() {
        &value[en]
    } else {
        &value[de]
    }
}

fn normalize(value: &str) -> String {
    let value = value
        .to_lowercase()
        .replace("2⁄3", " zwei drittel ")
        .replace("2/3", " zwei drittel ")
        .replace('α', " alpha ")
        .replace('₀', "0")
        .replace('₁', "1")
        .replace('₂', "2")
        .replace('₃', "3")
        .replace('₄', "4")
        .replace('₅', "5")
        .replace('₆', "6")
        .replace('₇', "7")
        .replace('₈', "8")
        .replace('₉', "9")
        .replace("μ", "mu");
    value
        .split(|c: char| !c.is_alphanumeric())
        .filter(|w| !w.is_empty())
        .collect::<Vec<_>>()
        .join(" ")
}

struct Concepts(Vec<(String, Vec<String>)>);
impl Concepts {
    fn new(map: &Value) -> Self {
        Self(
            map["query_concepts"]
                .as_object()
                .into_iter()
                .flatten()
                .map(|(id, aliases)| (id.clone(), strings(aliases).map(normalize).collect()))
                .collect(),
        )
    }
    fn identify(&self, text: &str) -> BTreeSet<String> {
        let text = normalize(text);
        let text = format!(" {text} ");
        self.0
            .iter()
            .filter(|(_, aliases)| {
                aliases.iter().any(|alias| {
                    if alias.is_empty() {
                        return false;
                    }
                    if alias.chars().count() <= 3 || alias.chars().any(|c| c.is_ascii_digit()) {
                        text.contains(&format!(" {alias} "))
                    } else {
                        text.contains(&format!(" {alias}"))
                    }
                })
            })
            .map(|(id, _)| id.clone())
            .collect()
    }
}

pub fn moved_source(reference: &Value) -> Option<Value> {
    if reference["source_location_state"] != "moved_to_external_repository" {
        return None;
    }
    let url = reference["source_url"].as_str()?;
    if !url.starts_with("https://") {
        return None;
    }
    Some(json!({"id":reference["source_id"], "source_url":url,
        "repository":reference["external_repository"],"revision":reference["external_revision"],
        "path":reference["external_path"],"sha256":reference["external_sha256"],
        "line_start":reference["line_start"],"line_end":reference["line_end"],
        "source_location_state":"moved_to_external_repository",
        "reference_role":reference["reference_role"]}))
}

pub fn map_records(map: &Value) -> impl Iterator<Item = &Value> {
    array(&map["patterns"]).chain(array(&map["nodes"]))
}

pub fn external_sources(value: &Value) -> Vec<Value> {
    array(&value["source_refs"])
        .filter_map(moved_source)
        .chain(array(&value["external_sources"]).cloned())
        .collect()
}

pub fn source_ids(value: &Value) -> Vec<&str> {
    let mut seen = HashSet::new();
    array(&value["source_refs"])
        .filter(|source| source["source_location_state"] != "moved_to_external_repository")
        .filter_map(|source| source["source_id"].as_str().or_else(|| source.as_str()))
        .filter(|id| seen.insert(*id))
        .take(12)
        .collect()
}

fn pattern_card(pattern: &Value, lang: &str) -> Value {
    let guidance = |key: &str| {
        array(&pattern[key])
            .map(|item| {
                if item.is_string() {
                    item.clone()
                } else {
                    localized(item, "label_de", "label_en", lang).clone()
                }
            })
            .collect::<Vec<_>>()
    };
    json!({
        "id":pattern["id"],"relation_type":pattern["relation_type"],"reviewed_as_recorded":pattern["reviewed"],
        "title":localized(pattern,"label_de","label_en",lang),
        "statement":localized(pattern,"summary_de","summary_en",lang),
        "conditions":guidance("conditions"),"counterevidence":guidance("counterevidence"),
        "readings":array(&pattern["readings"]).map(|reading|json!({"concept":reading["concept"],
            "scope":reading["scope"],"statement":localized(reading,"label_de","label_en",lang),
            "source_id":reading["source_id"]})).collect::<Vec<_>>(),
        "related_ids":pattern["related_ids"],"source_ids":source_ids(pattern),
        "external_sources":external_sources(pattern)
    })
}

fn ranked_patterns<'a>(map: &'a Value, wanted: &BTreeSet<String>) -> Vec<(usize, &'a Value)> {
    let mut result: Vec<_> = array(&map["patterns"])
        .filter_map(|pattern| {
            let hits = strings(&pattern["concepts"])
                .filter(|term| wanted.contains(*term))
                .count();
            (hits > 0).then_some((hits, pattern))
        })
        .collect();
    result.sort_by(|(na, a), (nb, b)| {
        nb.cmp(na)
            .then_with(|| {
                array(&a["concepts"])
                    .count()
                    .cmp(&array(&b["concepts"]).count())
            })
            .then_with(|| a["id"].as_str().cmp(&b["id"].as_str()))
    });
    result
}

pub fn pattern_search(map: &Value, query: &str, limit: usize) -> Option<Vec<Value>> {
    let wanted = Concepts::new(map).identify(query);
    if query.trim().is_empty() {
        return Some(
            array(&map["patterns"])
                .take(limit)
                .map(|pattern| {
                    let mut card = pattern_card(pattern, "de");
                    card["label"] = card["title"].clone();
                    card["summary"] = card["statement"].clone();
                    card["type"] = json!("theory_pattern");
                    card
                })
                .collect(),
        );
    }
    if wanted.is_empty() {
        return None;
    }
    let ranks = ranked_patterns(map, &wanted);
    let best = ranks.first().map_or(0, |(hits, _)| *hits);
    Some(
        ranks
            .into_iter()
            .filter(|(hits, _)| *hits >= best.min(2))
            .take(limit)
            .map(|(hits, pattern)| {
                let mut card = pattern_card(pattern, "de");
                card["label"] = card["title"].clone();
                card["summary"] = card["statement"].clone();
                card["type"] = json!("theory_pattern");
                card["score"] = json!(hits as f64 / wanted.len() as f64);
                card
            })
            .collect(),
    )
}

pub fn comparison(map: &Value, ids: &[String]) -> Vec<Value> {
    array(&map["patterns"])
        .filter_map(|pattern| {
            let participants: Vec<_> = ids
                .iter()
                .filter(|id| {
                    pattern["id"].as_str() == Some(id.as_str())
                        || strings(&pattern["related_ids"]).any(|related| related == id.as_str())
                })
                .collect();
            if participants.len() < 2 {
                return None;
            }
            let mut card = pattern_card(pattern, "de");
            card["compared_ids"] = json!(participants);
            Some(card)
        })
        .take(8)
        .collect()
}

pub fn inventory(p: &Program, current: &HashMap<String, bool>) -> Value {
    let mut total = BTreeMap::<&str, usize>::new();
    let mut by_topic = BTreeMap::<&str, BTreeMap<&str, usize>>::new();
    let mut tasks = BTreeMap::<&str, usize>::new();
    for problem in &p.problems {
        for f in &problem.findings {
            let state = if current.get(&f.id) == Some(&true) {
                f.status.as_str()
            } else {
                "stale_source"
            };
            *total.entry(state).or_default() += 1;
            for topic in &problem.topics {
                *by_topic.entry(topic).or_default().entry(state).or_default() += 1;
            }
        }
        for task in &problem.tasks {
            *tasks
                .entry(crate::research::task_status(task, problem, current))
                .or_default() += 1;
        }
    }
    json!({"findings_by_status":total,"by_recorded_topic":by_topic,
        "tasks_by_effective_status":tasks,"problems":p.problems.len(),
        "notice":"Finding totals count each record once. Topic groups overlap through their recorded problem assignments. Stale sources and open tasks are counted separately."})
}

fn vocab(text: &str) -> HashSet<String> {
    normalize(text)
        .split_whitespace()
        .filter(|word| {
            word.chars().count() >= 4
                && !matches!(
                    *word,
                    "welche"
                        | "welcher"
                        | "welchen"
                        | "tatsächlich"
                        | "haben"
                        | "kann"
                        | "aus"
                        | "derselben"
                        | "bereits"
                        | "einem"
                        | "einen"
                        | "einer"
                        | "gleichen"
                        | "dieselbe"
                        | "damit"
                        | "what"
                        | "which"
                        | "with"
                        | "same"
                        | "does"
                        | "that"
                        | "this"
                        | "from"
                        | "have"
                        | "under"
                        | "zwischen"
                        | "gemeinsam"
                        | "gemeinsame"
                        | "fragen"
                        | "ergebnisse"
                        | "befunde"
                        | "muster"
                        | "ausschlüsse"
                        | "positiven"
                        | "negativen"
                )
        })
        .map(str::to_owned)
        .collect()
}

fn fit_text(value: &str, cap: usize) -> String {
    let mut cut: String = value.chars().take(cap).collect();
    if cut.len() < value.len() {
        cut.push('…');
    }
    cut
}

struct Match<'a> {
    finding: &'a Finding,
    problem: &'a Problem,
    score: usize,
}

fn finding_card(row: &Match<'_>, lang: &str) -> Value {
    let f = row.finding;
    json!({"id":f.id,"problem_id":row.problem.id,"title":f.title.get(lang),
        "excerpt":fit_text(f.statement.get(lang),420),"scope":fit_text(f.scope.get(lang),210),
        "status":f.status,"premise_ids":f.premise_ids,"topics":row.problem.topics,
        "sources":f.sources.iter().take(3).map(|source|json!({"id":crate::research::source_id(source),
            "line_start":source.line_start,"line_end":source.line_end})).collect::<Vec<_>>(),
        "source_count":f.sources.len(),"detail_query":f.id})
}

fn shared_source_groups(matches: &[Match<'_>]) -> Vec<Value> {
    // An exact source span shared by distinct problems is a candidate to inspect.
    // Keep each finding's premise set. Shared provenance is not a shared conclusion.
    let mut groups = BTreeMap::<String, Vec<&Match<'_>>>::new();
    for row in matches {
        for citation in &row.finding.sources {
            groups
                .entry(crate::research::source_id(citation))
                .or_default()
                .push(row);
        }
    }
    let mut ranked: Vec<_> = groups
        .into_iter()
        .filter(|(_, rows)| {
            rows.iter()
                .map(|row| row.problem.id.as_str())
                .collect::<BTreeSet<_>>()
                .len()
                >= 2
        })
        .collect();
    ranked.sort_by(|(id_a, a), (id_b, b)| b.len().cmp(&a.len()).then(id_a.cmp(id_b)));
    ranked
        .into_iter()
        .take(3)
        .map(|(id, rows)| {
            json!({
                "relation_type":"same_source_span_candidate","source_id":id,
                "matching_findings":rows.len(),
                "findings":rows.iter().take(4).map(|row|json!({"id":row.finding.id,
                    "problem_id":row.problem.id,"status":row.finding.status,
                    "premise_ids":row.finding.premise_ids})).collect::<Vec<_>>()
            })
        })
        .collect()
}

fn room(packet: &Value, max_chars: usize) -> bool {
    serde_json::to_string(packet).is_ok_and(|s| s.chars().count() + 100 <= max_chars)
}

pub fn question(
    p: &Program,
    root: &std::path::Path,
    map: &Value,
    query: &str,
    lang: &str,
    max_chars: usize,
) -> Option<Value> {
    let concepts = Concepts::new(map);
    let wanted = concepts.identify(query);
    let explicit_patterns = [
        "muster",
        "patterns",
        "gegenbefunde",
        "ausschlüsse",
        "counterexamples",
        "positive and negative",
        "positive und negative",
    ]
    .iter()
    .any(|term| query.to_lowercase().contains(term));
    if wanted.len() < 2 && !explicit_patterns {
        return None;
    }
    let patterns = ranked_patterns(map, &wanted);
    let best = patterns.first().map_or(0, |(hits, _)| *hits);
    if best < 2 && !explicit_patterns {
        return None;
    }
    let needle = vocab(query);
    let current = crate::research::freshness(p, root);
    let mut matches = Vec::new();
    let mut statuses = BTreeMap::<&str, usize>::new();
    for problem in &p.problems {
        for f in &problem.findings {
            let title = format!("{} {}", f.title.de, f.title.en);
            let body = format!("{title} {} {}", f.statement.de, f.statement.en);
            let title_tags = concepts.identify(&title);
            let body_tags = concepts.identify(&body);
            let title_vocab = vocab(&title);
            let title_words = needle.intersection(&title_vocab).count();
            let body_words = needle.intersection(&vocab(&body)).count();
            let tag_count = wanted.intersection(&body_tags).count();
            if tag_count == 0 && title_words < 2
                || tag_count == 1 && body_words < 2 && wanted.len() > 2
            {
                continue;
            }
            if current.get(&f.id) != Some(&true) {
                *statuses.entry("stale_source").or_default() += 1;
                continue;
            }
            *statuses.entry(f.status.as_str()).or_default() += 1;
            matches.push(Match {
                finding: f,
                problem,
                score: tag_count * 16
                    + wanted.intersection(&title_tags).count() * 9
                    + title_words * 6
                    + body_words,
            });
        }
    }
    matches.sort_by(|a, b| b.score.cmp(&a.score).then(a.finding.id.cmp(&b.finding.id)));
    let mut packet = json!({"revision":p.revision,"view":"cross_topic_evidence","query":query,"lang":lang,
    "concepts":wanted,"matched_findings_by_status":statuses,"matched_finding_total":matches.len(),
    "patterns":[],"findings":[],"recorded_connections":[],"shared_source_candidates":[],"problem_index":[],
    "budget":{"max_chars":max_chars},"notice":if lang=="en" {
        "Navigation across recorded scopes. Common labels and keywords do not establish a connection. Use exact IDs for full premises and sources."
    } else {
        "Querverweise zwischen dokumentierten Geltungsbereichen. Übereinstimmende Wörter begründen keine Verbindung. Vollständige Befunde über ihre Kennung öffnen."
    }});
    if explicit_patterns && wanted.is_empty() {
        let counts = inventory(p, &current);
        packet["evidence_inventory"] = json!({"findings_by_status":counts["findings_by_status"],
            "tasks_by_effective_status":counts["tasks_by_effective_status"],
            "topic_breakdown_tool":"theory_overview",
            "catalog":{"tool":"search_knowledge","type":"theory_pattern","query":""}});
    }
    for (_, pattern) in patterns
        .iter()
        .filter(|(hits, _)| *hits >= best.min(2))
        .take(4)
    {
        packet["patterns"]
            .as_array_mut()?
            .push(pattern_card(pattern, lang));
        // Reserve space for both a positive/conditional and a negative record.
        if serde_json::to_string(&packet).ok()?.chars().count() + 1900 > max_chars {
            packet["patterns"].as_array_mut()?.pop();
            break;
        }
    }
    let mut selected = BTreeSet::new();
    let mut problems = BTreeSet::new();
    for status in ["established_scoped", "excluded_scoped", "conditional"] {
        if let Some(row) = matches.iter().find(|row| row.finding.status == status) {
            packet["findings"]
                .as_array_mut()?
                .push(finding_card(row, lang));
            if !room(&packet, max_chars) {
                packet["findings"].as_array_mut()?.pop();
            } else {
                selected.insert(row.finding.id.as_str());
                problems.insert(row.problem.id.as_str());
            }
        }
    }
    for row in &matches {
        if selected.len() >= 7 {
            break;
        }
        if selected.contains(row.finding.id.as_str()) || problems.contains(row.problem.id.as_str())
        {
            continue;
        }
        packet["findings"]
            .as_array_mut()?
            .push(finding_card(row, lang));
        if !room(&packet, max_chars) {
            packet["findings"].as_array_mut()?.pop();
            break;
        }
        selected.insert(&row.finding.id);
        problems.insert(&row.problem.id);
    }
    for problem in &p.problems {
        for connection in &problem.connections {
            if packet["recorded_connections"].as_array()?.len() >= 4 {
                break;
            }
            if !connection
                .finding_ids
                .iter()
                .any(|id| selected.contains(id.as_str()))
            {
                continue;
            }
            packet["recorded_connections"]
                .as_array_mut()?
                .push(json!({"id":connection.id,
                "from":problem.id,"to":connection.target_problem_id,"type":connection.kind,
                "finding_ids":connection.finding_ids,"premise_ids":connection.premise_ids}));
            if !room(&packet, max_chars) {
                packet["recorded_connections"].as_array_mut()?.pop();
                break;
            }
        }
    }
    for row in matches
        .iter()
        .filter(|row| selected.contains(row.finding.id.as_str()))
    {
        if array(&packet["problem_index"]).any(|item| item["id"] == row.problem.id) {
            continue;
        }
        packet["problem_index"]
            .as_array_mut()?
            .push(json!({"id":row.problem.id,"title":row.problem.title.get(lang)}));
        if !room(&packet, max_chars) {
            packet["problem_index"].as_array_mut()?.pop();
            break;
        }
    }
    for group in shared_source_groups(&matches) {
        packet["shared_source_candidates"]
            .as_array_mut()?
            .push(group);
        if !room(&packet, max_chars) {
            packet["shared_source_candidates"].as_array_mut()?.pop();
            break;
        }
    }
    let remaining = matches.len().saturating_sub(selected.len());
    packet["budget"]["omitted_matching_findings"] = json!(remaining);
    packet["budget"]["serialized_chars"] =
        json!(serde_json::to_string(&packet).ok()?.chars().count());
    (serde_json::to_string(&packet).ok()?.chars().count() <= max_chars).then_some(packet)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compares_authored_relation_without_inventing_one_for_unrelated_ids() {
        let map = json!({"patterns":[{
            "id":"pattern:test", "label_de":"Gemeinsamer Wert", "summary_de":"Zwei Observablen",
            "relation_type":"shared_value_distinct_observables", "related_ids":["map:left","claim:RIGHT"],
            "counterevidence":[{"label_de":"Verschiedene Operatoren"}]
        }]});
        let result = comparison(&map, &["map:left".into(), "claim:RIGHT".into()]);
        assert_eq!(
            result[0]["relation_type"],
            "shared_value_distinct_observables"
        );
        assert_eq!(result[0]["counterevidence"][0], "Verschiedene Operatoren");
        assert!(comparison(&map, &["map:left".into(), "claim:OTHER".into()]).is_empty());
    }

    #[test]
    fn recognizes_spoken_fraction_and_subscript_operator_names() {
        let map = json!({"query_concepts":{"fraction":["2/3"],"clock":["c₆"],"alpha":["α"]}});
        let tags = Concepts::new(&map);
        assert_eq!(
            tags.identify("zwei Drittel an C6 und α"),
            BTreeSet::from(["fraction".into(), "clock".into(), "alpha".into()])
        );
        assert!(tags.identify("AFC61").is_empty());
    }
}
