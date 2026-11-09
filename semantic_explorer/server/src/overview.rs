use serde_json::{json, Value};

const AREAS: &str = include_str!("../../../website/lib/atlas-status.json");

pub fn areas(lang: &str, compact: bool) -> Value {
    let entries: Vec<Value> = serde_json::from_str(AREAS).expect("checked-in TFPT status areas");
    json!(entries.iter().map(|area| {
        let mut entry = json!({"id":area["id"], "title":area["title"][lang], "evidence":area["evidence"]});
        if !compact {
            entry["result"] = area["result"][lang].clone();
            entry["scope"] = area["scope"][lang].clone();
            entry["public_href"] = area["href"].clone();
        }
        entry
    }).collect::<Vec<_>>())
}

pub fn requested(query: &str) -> bool {
    let lower = query.trim().to_lowercase();
    if lower.is_empty() {
        return true;
    }
    let words: Vec<_> = lower.split(|c: char| !c.is_alphanumeric()).collect();
    words.iter().any(|word| matches!(*word, "gesamtstatus" | "gesamtstand" | "gesamtüberblick" | "gesamtbild"))
        || ["overall status", "overall picture", "big picture", "theory overview", "state of tfpt", "status of tfpt",
            "stand von tfpt", "tfpt stand", "tfpt status", "tfpt überblick", "tfpt overview",
            "wo stehen wir insgesamt"].iter().any(|phrase| lower.contains(phrase))
        || (words.contains(&"tfpt")
            && words.iter().any(|word| matches!(*word, "insgesamt" | "overall" | "überblick" | "overview")))
}
