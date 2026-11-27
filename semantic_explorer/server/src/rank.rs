//! Query analysis and relevance weights for the Atlas search.
//!
//! The search used to test raw substrings and to add a fixed number per token.
//! Function words therefore matched inside identifiers ("der" in "derivation"),
//! every token counted the same however common it was, and a German question
//! could not reach an English record.  This module supplies the four missing
//! pieces: stop words, whole-word matching with light stems, a bilingual
//! glossary of the theory's vocabulary, and rarity weights.
//!
//! Nothing here judges content.  The weights order navigation hits; the ledger
//! and the original sources own every status.

use std::sync::OnceLock;

/// One query concept: the word as typed plus the forms that count as a match.
#[derive(Debug, Clone)]
pub struct Term {
    pub token: String,
    pub variants: Vec<String>,
    /// Forms of the head or tail of a German compound ("Twistorbrücke" names
    /// a twistor); they count for less than the word itself.
    pub parts: Vec<String>,
}

const STOPWORDS: &[&str] = &[
    // German
    "der", "die", "das", "dem", "den", "des", "ein", "eine", "einer", "einem", "einen", "eines",
    "und", "oder", "aber", "als", "auch", "nicht", "ist", "sind", "war", "wird", "werden", "wurde",
    "hat", "haben", "es", "sie", "er", "wir", "ich", "dass", "ob", "wenn", "dann", "zu", "zur",
    "zum", "von", "vom", "im", "in", "an", "am", "auf", "aus", "bei", "mit", "nach", "ueber",
    "unter", "fuer", "durch", "wie", "was", "wo", "welche", "welcher", "welches", "gibt", "sich",
    "noch", "schon", "nur", "man", "so", "da", "steht", "stehen", "stand", "fest", "ohne", "bereits",
    "bekannt", "kann", "koennen", "muss", "soll", "hier", "dort", "dazu", "dabei", "alle", "alles",
    "mehr", "sehr", "dies", "diese", "dieser", "dieses", "jede", "jeder", "kein", "keine", "sein",
    "seine", "ihr", "ihre", "zwischen", "gegen", "bis", "um", "vor", "wegen", "weil", "damit", "denn",
    "doch", "nun", "also", "etwa", "bzw", "sowie", "genau", "eigentlich", "wurden", "worden", "habe",
    "welchen", "welchem", "wieso", "warum", "woher", "wohin", "kommen", "kommt", "liegt", "geht",
    // English
    "the", "a", "of", "and", "or", "but", "is", "are", "was", "were", "be", "been", "by", "for",
    "from", "has", "have", "had", "into", "it", "its", "on", "that", "this", "these", "those",
    "to", "which", "with", "does", "do", "did", "can", "about", "there", "their", "than", "then",
    "if", "not", "what", "how", "where", "why", "when", "we", "at", "as", "any", "all", "some", "will",
    "would", "should", "could", "may", "might", "must", "our", "your", "they", "them", "you", "also",
    "only", "just", "more", "most", "such", "other", "being", "between", "after", "before", "because",
    "while", "during", "without", "already", "known", "come", "comes",
];

/// Groups of forms that name the same concept in the theory's two working
/// languages.  Entries are folded on first use, so they may be written with
/// umlauts.  A phrase entry matches as a phrase.
const GLOSSARY: &[&[&str]] = &[
    &["doppelbrechung", "birefringenz", "birefringence", "polarisationsdrehung", "polarization rotation", "polarisation rotation"],
    &["herleitung", "ableitung", "derivation", "hergeleitet", "derived"],
    &["vorhersage", "prediction", "prognose", "predicted"],
    &["feinstruktur", "feinstrukturkonstante", "fine structure", "fine-structure", "alpha"],
    &["naht", "seam"],
    &["rand", "boundary"],
    &["randkern", "boundary kernel"],
    &["uhr", "clock", "uhren", "clocks"],
    &["gitter", "lattice"],
    &["wurzel", "root", "wurzeln", "roots"],
    &["spiegelung", "reflection", "reflexion"],
    &["reflexionspositivität", "reflexionspositiv", "reflection positivity", "reflection positive", "reflection-positive", "rp"],
    &["zustand", "state"],
    &["zeit", "time"],
    &["schwerkraft", "gravitation", "gravity", "gravitativ", "gravitational"],
    &["raumzeit", "spacetime", "space-time"],
    &["krümmung", "curvature"],
    &["fläche", "area"],
    &["horizont", "horizon"],
    &["dunkle energie", "dark energy", "kosmologische konstante", "cosmological constant"],
    &["dunkle materie", "dark matter"],
    &["baryon", "baryonen", "baryons", "baryonic"],
    &["materie", "matter"],
    &["masse", "mass", "massen", "masses"],
    &["mischung", "mixing", "mischungswinkel", "mixing angle"],
    &["familie", "family", "generation", "familien", "families", "generationen", "generations"],
    &["eichgruppe", "gauge group", "eichsymmetrie", "gauge symmetry"],
    &["eichfeld", "gauge field", "eichung", "gauge"],
    &["standardmodell", "standard model"],
    &["chiralität", "chirality", "chiral", "händigkeit", "handedness"],
    &["anomalie", "anomaly", "anomalien", "anomalies"],
    &["doppelüberlagerung", "double cover", "doppelte überlagerung", "zweiblättrig", "überlagerung", "covering"],
    &["träger", "carrier"],
    &["marke", "mark", "marken", "marks", "markierung", "marking"],
    &["viermarken", "four marks", "fourmark", "vier marken", "four-mark"],
    &["windung", "winding", "windungszahl", "winding number"],
    &["kern", "kernel"],
    &["saat", "seed", "keim"],
    &["ausschluss", "exclusion", "ausgeschlossen", "excluded", "no-go", "nogo"],
    &["offen", "open"],
    &["lücke", "gap"],
    &["bedingt", "conditional"],
    &["beweis", "proof", "bewiesen", "proven", "proved"],
    &["satz", "theorem"],
    &["messung", "measurement", "messwert"],
    &["gemessen", "measured"],
    &["daten", "data"],
    &["kopplung", "coupling", "kopplungskonstante", "coupling constant"],
    &["skala", "scale"],
    &["vereinheitlichung", "unification"],
    &["protonzerfall", "protonenzerfall", "proton decay"],
    &["skalaron", "scalaron"],
    &["kosmologie", "cosmology", "kosmologisch", "kosmisch", "cosmological", "cosmic"],
    &["hintergrundstrahlung", "mikrowellenhintergrund", "cosmic microwave background", "cmb"],
    &["verschränkung", "entanglement"],
    &["entropie", "entropy"],
    &["wirkung", "action"],
    &["bewegungsgleichung", "equation of motion", "feldgleichung", "field equation"],
    &["einsteingleichung", "einstein equation", "einstein equations", "einsteingleichungen", "einstein-gleichung", "einstein eq", "einstein tensor", "einstein-hilbert"],
    &["gleichung", "equation", "gleichungen", "equations"],
    &["zerlegung", "decomposition"],
    &["darstellung", "representation", "darstellungen", "representations"],
    &["spinor", "spinoren", "spinors"],
    &["halbspinor", "half-spinor", "half spinor"],
    &["signatur", "signature"],
    &["twistor", "twistoren", "twistors"],
    &["twistorraum", "twistor space"],
    &["brücke", "bridge"],
    &["spiegel", "mirror"],
    &["spiegelfermion", "spiegelfermionen", "mirror fermion", "mirror fermions", "nomirror", "no mirror"],
    &["fermion", "fermionen", "fermions"],
    &["boson", "bosonen", "bosons"],
    &["kreis", "circle"],
    &["gerade", "line"],
    &["quelle", "source"],
    &["himmelskugel", "himmelssphäre", "celestial sphere", "celestial"],
    &["helizität", "helicity"],
    &["elektromagnetismus", "electromagnetism", "elektromagnetisch", "electromagnetic"],
    &["starke wechselwirkung", "strong interaction", "qcd"],
    &["farbe", "color", "colour"],
    &["schwache wechselwirkung", "weak interaction"],
    &["schwach", "weak"],
    &["wechselwirkung", "interaction"],
    &["lepton", "leptonen", "leptons"],
    &["cp-verletzung", "cp violation", "cp-phase", "cp phase"],
    &["fehlerkorrektur", "error correction", "qec"],
    &["quantisierung", "quantization", "quantisation"],
    &["kontinuum", "continuum"],
    &["kontinuumslimes", "continuum limit"],
    &["lokal", "local"],
    &["lokalität", "locality"],
    &["netz", "net"],
    &["kausal", "causal", "kausalität", "causality"],
    &["kausalstruktur", "causal structure"],
    &["lichtkegel", "light cone", "lightcone"],
    &["temperatur", "temperature"],
    &["thermisch", "thermal"],
    &["modularfluss", "modular flow", "tomita"],
    &["vakuum", "vacuum"],
    &["energie", "energy"],
    &["impuls", "momentum"],
    &["drehung", "rotation"],
    &["vierteldrehung", "quarter turn", "quarter rotation"],
    &["zentrum", "centre", "center"],
    &["klebung", "verklebung", "glue", "gluing", "sewing"],
    &["ursprung", "origin", "herkunft"],
    &["universalraum", "universal space"],
    &["prozessgesetz", "process law"],
    &["auswahl", "selection"],
    &["selektor", "selector", "wähler"],
    &["eindeutigkeit", "uniqueness", "eindeutig", "unique"],
    &["existenz", "existence"],
    &["fixpunkt", "fixed point", "fixed-point"],
    &["verhältnis", "ratio"],
    &["winkel", "angle"],
    &["determinante", "determinant"],
    &["spur", "trace"],
    &["eigenwert", "eigenvalue"],
    &["geschlecht", "genus"],
    &["kurve", "curve"],
    &["kugel", "sphäre", "sphere"],
    &["streumatrix", "s-matrix", "s matrix"],
    &["streuung", "scattering"],
    &["einschleifen", "one-loop", "one loop"],
    &["renormierung", "renormalization", "renormalisation"],
    &["schwelle", "threshold"],
    &["dichte", "density"],
    &["hubble-konstante", "hubble constant", "hubble"],
    &["spektralindex", "spectral index"],
    &["tensor-zu-skalar", "tensor-to-scalar", "tensor to scalar"],
    &["gravitationswelle", "gravitationswellen", "gravitational wave", "gravitational waves"],
    &["schwarzes loch", "schwarze löcher", "black hole", "black holes"],
    &["killtest", "kill-test", "kill test", "falsifikation", "falsification"],
    &["eingefroren", "frozen"],
    &["prüfung", "check", "test"],
    &["widerspruch", "contradiction"],
    &["spannung", "tension"],
    &["realform", "reelle form", "real form"],
    &["orientierung", "orientation", "orientiert", "oriented"],
    &["zeitrichtung", "arrow of time", "zeitpfeil"],
    &["gleichgewicht", "equilibrium"],
];

fn glossary() -> &'static Vec<Vec<String>> {
    static FOLDED: OnceLock<Vec<Vec<String>>> = OnceLock::new();
    FOLDED.get_or_init(|| {
        GLOSSARY
            .iter()
            .map(|group| group.iter().map(|entry| fold(entry)).collect())
            .collect()
    })
}

/// Function words of both working languages, in any spelling.
pub fn is_stopword(word: &str) -> bool {
    STOPWORDS.contains(&fold(word).as_str())
}

/// Lower case, transliterated umlauts, spelled-out Greek letters and plain
/// digits for subscripts, so that "φ₀", "phi0" and "Änderung"/"aenderung" meet.
pub fn fold(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 8);
    for c in text.chars() {
        match c {
            'ä' | 'Ä' => out.push_str("ae"),
            'ö' | 'Ö' => out.push_str("oe"),
            'ü' | 'Ü' => out.push_str("ue"),
            'ß' => out.push_str("ss"),
            'α' => out.push_str("alpha"),
            'β' => out.push_str("beta"),
            'γ' | 'Γ' => out.push_str("gamma"),
            'δ' | 'Δ' => out.push_str("delta"),
            'θ' | 'Θ' => out.push_str("theta"),
            'λ' | 'Λ' => out.push_str("lambda"),
            'μ' | 'µ' => out.push_str("mu"),
            'π' => out.push_str("pi"),
            'ρ' => out.push_str("rho"),
            'σ' | 'Σ' => out.push_str("sigma"),
            'τ' => out.push_str("tau"),
            'φ' | 'ϕ' | 'Φ' => out.push_str("phi"),
            'ω' | 'Ω' => out.push_str("omega"),
            'ℤ' => out.push('z'),
            '₀'..='₉' => out.push(char::from(b'0' + (c as u32 - '₀' as u32) as u8)),
            _ => out.extend(c.to_lowercase()),
        }
    }
    out
}

/// Strip one common German or English ending from a long word.
fn stem(word: &str) -> String {
    let length = word.chars().count();
    if length < 6 || !word.chars().all(char::is_alphabetic) {
        return word.to_owned();
    }
    for suffix in ["en", "er", "es", "em", "e", "n", "s"] {
        if let Some(rest) = word.strip_suffix(suffix) {
            if rest.chars().count() >= 5 {
                return rest.to_owned();
            }
        }
    }
    word.to_owned()
}

fn same_word(entry: &str, word: &str, stem: &str) -> bool {
    if entry == word || entry == stem {
        return true;
    }
    let near = |long: &str, short: &str| {
        short.chars().count() >= 5
            && long.starts_with(short)
            && long.chars().count() - short.chars().count() <= 3
    };
    near(entry, stem) || near(word, entry)
}

/// Split a query into concepts.  Function words are dropped unless nothing
/// else is left; glossary phrases are recognised before single words.
pub fn query_terms(query: &str) -> Vec<Term> {
    let folded = fold(query);
    let words: Vec<&str> = folded
        .split(|c: char| !c.is_alphanumeric())
        .filter(|word| !word.is_empty())
        .collect();
    let mut used = vec![false; words.len()];
    let mut found: Vec<(usize, Term)> = Vec::new();
    for group in glossary() {
        for entry in group.iter().filter(|entry| entry.contains([' ', '-'])) {
            let parts: Vec<&str> = entry
                .split(|c: char| !c.is_alphanumeric())
                .filter(|part| !part.is_empty())
                .collect();
            if parts.len() < 2 || parts.len() > words.len() {
                continue;
            }
            for start in 0..=words.len() - parts.len() {
                if (0..parts.len()).all(|k| !used[start + k] && words[start + k] == parts[k]) {
                    (0..parts.len()).for_each(|k| used[start + k] = true);
                    let mut variants = vec![entry.clone()];
                    variants.extend(group.iter().filter(|other| *other != entry).cloned());
                    found.push((start, Term { token: entry.clone(), variants, parts: Vec::new() }));
                }
            }
        }
    }
    let single = |word: &str| -> Term {
        let stem = stem(word);
        let mut variants = vec![word.to_owned()];
        if stem != word {
            variants.push(stem.clone());
        }
        for group in glossary() {
            if group
                .iter()
                .any(|entry| !entry.contains(' ') && same_word(entry, word, &stem))
            {
                for entry in group {
                    if !variants.contains(entry) {
                        variants.push(entry.clone());
                    }
                }
            }
        }
        // A long word that the glossary does not know may be a compound.  Only
        // its longest known head and its longest known tail count: "…gleichung"
        // is an equation and must not be read as "…eichung" (gauge).
        let mut parts = Vec::new();
        if word.chars().count() >= 9 && variants.len() <= 2 {
            let known = |entry: &&String| !entry.contains([' ', '-']) && entry.chars().count() >= 4 && entry.as_str() != word;
            let longest = |head: bool| {
                glossary()
                    .iter()
                    .flat_map(|group| group.iter().filter(known).map(move |entry| (entry, group)))
                    .filter(|(entry, _)| if head { word.starts_with(entry.as_str()) } else { word.ends_with(entry.as_str()) })
                    .max_by_key(|(entry, _)| entry.chars().count())
                    .map(|(_, group)| group)
            };
            for group in [longest(true), longest(false)].into_iter().flatten() {
                for entry in group {
                    if !variants.contains(entry) && !parts.contains(entry) {
                        parts.push(entry.clone());
                    }
                }
            }
        }
        Term { token: word.to_owned(), variants, parts }
    };
    for (index, word) in words.iter().enumerate() {
        if used[index] || word.chars().count() < 2 || STOPWORDS.contains(word) {
            continue;
        }
        found.push((index, single(word)));
    }
    if found.is_empty() {
        for (index, word) in words.iter().enumerate() {
            if word.chars().count() >= 2 {
                found.push((index, single(word)));
            }
        }
    }
    found.sort_by_key(|(index, _)| *index);
    let mut terms: Vec<Term> = Vec::new();
    for (_, term) in found {
        if !terms.iter().any(|known| known.token == term.token) {
            terms.push(term);
        }
    }
    terms
}

/// How well one form occurs in a folded text: 1 for a whole word, less for a
/// word beginning or a part of a compound, 0 for a chance substring.
fn occurrence(text: &str, form: &str) -> f32 {
    if form.is_empty() {
        return 0.0;
    }
    let length = form.chars().count();
    let letters_only = !form.ends_with(|c: char| c.is_ascii_digit());
    let mut best = 0.0f32;
    for (start, _) in text.match_indices(form) {
        let before = text[..start].chars().next_back();
        let after = text[start + form.len()..].chars().next();
        let at_start = before.is_none_or(|c| !c.is_alphanumeric());
        let at_end = after.is_none_or(|c| !c.is_alphanumeric());
        let strength = if at_start && at_end {
            1.0
        } else if at_start && letters_only && after.is_some_and(|c| c.is_ascii_digit()) {
            1.0
        } else if at_start && length >= 5 {
            0.85
        } else if length >= 7 {
            0.6
        } else {
            0.0
        };
        if strength > best {
            best = strength;
            if best >= 1.0 {
                break;
            }
        }
    }
    best
}

/// Match strength of a concept in a folded text.  The typed word counts in
/// full, its stem and its glossary partners slightly less.
pub fn term_hit(text: &str, term: &Term) -> f32 {
    let mut best = 0.0f32;
    for (index, form) in term.variants.iter().enumerate() {
        let factor = if index == 0 { 1.0 } else { 0.9 };
        let value = occurrence(text, form) * factor;
        if value > best {
            best = value;
        }
        if best >= 1.0 {
            break;
        }
    }
    if best < 0.6 {
        for form in &term.parts {
            let value = occurrence(text, form) * 0.6;
            if value > best {
                best = value;
            }
        }
    }
    best
}

/// Best weighted field for a concept.  `fields` must be ordered by falling
/// weight; the scan stops when no later field can win.
pub fn best_field(term: &Term, fields: &[(f32, &str)]) -> f32 {
    let mut best = 0.0f32;
    for (weight, text) in fields {
        if best >= *weight {
            break;
        }
        if text.is_empty() {
            continue;
        }
        let value = weight * term_hit(text, term);
        if value > best {
            best = value;
        }
    }
    best
}

/// Rarity weight per concept, 1 for the rarest one that occurs at all.  A
/// concept that occurs nowhere gets 0 and cannot dilute the others.
pub fn weights(frequency: &[usize], population: usize) -> Vec<f32> {
    let rarity: Vec<f32> = frequency
        .iter()
        .map(|&count| {
            if count == 0 {
                0.0
            } else {
                ((population as f32 + 1.0) / (count as f32 + 1.0)).ln().max(0.05)
            }
        })
        .collect();
    let top = rarity.iter().copied().fold(0.0f32, f32::max);
    if top <= 0.0 {
        return vec![1.0; frequency.len()];
    }
    rarity
        .iter()
        .map(|&value| if value <= 0.0 { 0.0 } else { (value / top).max(0.12) })
        .collect()
}

/// Lexical score of one record.  A record that matches every concept in its
/// identifier reaches 9 per concept; missing the rare concepts costs more than
/// missing the common ones, and partial coverage is damped.
pub fn combine(hits: &[f32], weights: &[f32]) -> f32 {
    let total: f32 = weights.iter().sum();
    if total <= 0.0 || hits.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0f32;
    let mut covered = 0.0f32;
    for (hit, weight) in hits.iter().zip(weights) {
        sum += hit * weight;
        if *hit > 0.0 {
            covered += weight;
        }
    }
    let scale = hits.len() as f32 / total;
    sum * scale * (0.35 + 0.65 * covered / total)
}

/// Navigation rank of a record type.  Ledger claims and reviewed research
/// records lead; sub-records of a claim and archived paper passages follow.
/// This orders hits and says nothing about truth.
pub fn authority(kind: &str, content_scope: Option<&str>, paper: bool) -> f32 {
    if paper {
        return match content_scope {
            Some("archive" | "historical") => 0.7,
            _ => 1.0,
        };
    }
    match kind {
        "claim" => 1.5,
        "research_finding" | "theory_pattern" | "theory_map_aspect" => 1.35,
        "research_problem" | "research_task" | "topic" | "goal" | "structure" | "theory_aspect" => 1.2,
        "attempt" | "contract" | "experiment" => 1.1,
        "assumption" | "result" | "check" | "object" | "derivation" => 0.85,
        "script" | "lean_theorem" | "lean_declaration" | "experiment_source" | "tex_source"
        | "paper_source" | "experiment_script" | "experiment_record" => 0.9,
        _ => 1.0,
    }
}

/// Record types whose summary is their own statement and not a description.
pub fn statement_kind(kind: &str) -> bool {
    matches!(
        kind,
        "claim" | "research_finding" | "research_problem" | "research_task" | "attempt" | "contract"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokens(query: &str) -> Vec<String> {
        query_terms(query).into_iter().map(|term| term.token).collect()
    }

    #[test]
    fn function_words_are_dropped_and_kept_only_when_alone() {
        assert_eq!(
            tokens("Herleitung der kosmischen Doppelbrechung"),
            ["herleitung", "kosmischen", "doppelbrechung"]
        );
        assert_eq!(tokens("der die das"), ["der", "die", "das"]);
    }

    #[test]
    fn a_german_question_reaches_english_records() {
        let terms = query_terms("Herleitung der kosmischen Doppelbrechung");
        let text = fold("Cosmic birefringence: derivation of the rotation angle");
        assert!(terms.iter().all(|term| term_hit(&text, term) >= 0.85));
        let phrase = query_terms("fine structure constant");
        assert_eq!(phrase[0].token, "fine structure");
        assert!(term_hit(&fold("Feinstrukturkonstante und α-Fixpunkt"), &phrase[0]) > 0.8);
    }

    #[test]
    fn a_function_word_no_longer_matches_inside_an_identifier() {
        let der = Term { token: "der".into(), variants: vec!["der".into()], parts: Vec::new() };
        assert_eq!(term_hit("derivation:celest.wp5e.beta.01:d1", &der), 0.0);
        let beta = &query_terms("beta")[0];
        assert_eq!(term_hit("woit.beta3.duality.01", beta), 1.0);
        assert_eq!(term_hit("alphabetagamma", beta), 0.0);
    }

    #[test]
    fn a_compound_reaches_records_on_its_parts() {
        let bridge = &query_terms("Twistorbrücke")[0];
        let hit = term_hit(&fold("clock and Woit structure on twistor space force signature (2,2)"), bridge);
        assert!(hit > 0.5 && hit < 0.9);
        let question = query_terms("Was steht zur Chiralität ohne Spiegelfermionen fest?");
        assert_eq!(question.iter().map(|term| term.token.as_str()).collect::<Vec<_>>(), ["chiralitaet", "spiegelfermionen"]);
        let text = fold("CHIRAL4D.NOMIRROR.01: one chiral generation, no mirror fermions");
        assert!(question.iter().all(|term| term_hit(&text, term) >= 0.85));
    }

    #[test]
    fn an_equation_is_not_a_gauge() {
        let field = &query_terms("Wellengleichungen")[0];
        assert!(field.parts.iter().any(|part| part == "equation"));
        assert!(!field.parts.iter().any(|part| part == "gauge"));
        let einstein = &query_terms("Einsteingleichung")[0];
        assert!(einstein.parts.is_empty());
        assert!(term_hit(&fold("the full covariant Einstein equation, parameter-free"), einstein) >= 0.9);
        assert_eq!(term_hit(&fold("gauge.detline.fixpoint.01"), einstein), 0.0);
    }

    #[test]
    fn folding_joins_symbols_and_spellings() {
        assert_eq!(fold("φ₀ und E₈, Änderung"), "phi0 und e8, aenderung");
        let term = &query_terms("phi0")[0];
        assert_eq!(term_hit(&fold("Gemeinsamer Seed φ₀"), term), 1.0);
    }

    #[test]
    fn rare_concepts_outweigh_common_ones() {
        let weight = weights(&[150, 3000, 0], 65_000);
        assert_eq!(weight[0], 1.0);
        assert!(weight[1] < 0.6 && weight[1] > 0.12);
        assert_eq!(weight[2], 0.0);
        let rare_only = combine(&[3.0, 0.0, 0.0], &weight);
        let common_only = combine(&[0.0, 9.0, 0.0], &weight);
        assert!(rare_only > common_only * 0.5);
        let both = combine(&[3.0, 3.0, 0.0], &weight);
        assert!(both > rare_only && both > common_only);
    }

    #[test]
    fn claims_lead_archived_paper_passages() {
        assert!(authority("claim", None, false) > authority("paper_section", Some("active"), true));
        assert!(authority("paper_section", Some("active"), true) > authority("paper_section", Some("archive"), true));
        assert!(authority("derivation", None, false) < 1.0);
    }
}
