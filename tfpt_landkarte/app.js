/* TFPT Landkarte · Oberfläche. Ohne Abhängigkeiten, ohne Build. */
(function () {
  "use strict";
  const T = window.TFPT;
  const A = window.ATLAS || { problems: [], revision: "", generated: "" };
  const main = document.getElementById("main");
  const drawer = document.getElementById("drawer");
  const byId = Object.fromEntries(T.nodes.map((n) => [n.id, n]));
  const edges = [];
  const nbr = {};
  T.nodes.forEach((n) => (nbr[n.id] = []));
  T.nodes.forEach((n) => (n.l || []).forEach(([to, type]) => {
    if (!byId[to]) return;
    edges.push({ a: n.id, b: to, type });
    nbr[n.id].push({ id: to, type, out: true });
    nbr[to].push({ id: n.id, type, out: false });
  }));
  const state = { num: null, tag: null, kinds: new Set(Object.keys(T.kinds)), dynSort: null, dynSt: null, strat: "", stratQ: "", atlasQ: "", mode: "vor" };

  /* ── kleine Helfer ── */
  function h(tag, a, ...c) {
    const e = document.createElement(tag);
    for (const k in a || {}) {
      const v = a[k];
      if (v == null || v === false) continue;
      if (k === "class") e.className = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v === true ? "" : v);
    }
    for (const x of c.flat(Infinity)) if (x != null && x !== false) e.append(x.nodeType ? x : document.createTextNode(String(x)));
    return e;
  }
  function s(tag, a, ...c) {
    const e = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (const k in a || {}) if (a[k] != null) e.setAttribute(k, a[k]);
    for (const x of c.flat()) if (x != null) e.append(x.nodeType ? x : document.createTextNode(String(x)));
    return e;
  }
  const pill = (st) => h("span", { class: "pill st-" + st, title: T.status[st][1] }, T.status[st][0]);
  const kindTag = (k) => h("span", { class: "kind k-" + k }, T.kinds[k]);
  const go = (view) => (location.hash = "#/" + view);
  function card(n) {
    return h("button", { class: "card " + n.st, onclick: () => open(n.id) },
      h("div", { class: "row" }, kindTag(n.k), pill(n.st)), h("h3", {}, n.t), h("p", {}, n.s));
  }
  function legendStatus() {
    return h("div", { class: "legend" }, Object.keys(T.status).map((k) => h("span", {}, pill(k), " ", T.status[k][1])));
  }

  /* ── Detailansicht ── */
  function open(id) {
    const n = byId[id];
    if (!n) return;
    drawer.hidden = false;
    drawer.replaceChildren(...[
      h("button", { class: "x", onclick: close, "aria-label": "Schließen" }, "Schließen"),
      kindTag(n.k), h("h2", {}, n.t), pill(n.st),
      h("p", { class: "short" }, n.s),
      (n.b || []).map((p) => h("p", {}, p)),
      n.req && h("div", {}, h("h4", {}, "Was dieser Ansatz liefert"),
        T.req.map(([k, name]) => h("span", { class: "chip", title: cellWord(n.req[k]) }, cellSym(n.req[k]), " ", name))),
      n.k === "hyp" && hypBody(n),
      ["same", "diff", "open"].some((k) => n[k] && n[k].length) && triple(n),
      n.n && n.n.length ? h("div", {}, h("h4", {}, "Zahlen · antippen zeigt, wo sie noch auftauchen"),
        n.n.map((x) => h("button", { class: "chip", onclick: () => { state.num = x; state.tag = null; close(); go("ideen"); render(); } }, x))) : null,
      n.g && n.g.length ? h("div", {}, h("h4", {}, "Stichworte"),
        n.g.map((x) => h("button", { class: "chip", onclick: () => { state.tag = x; state.num = null; close(); go("ideen"); render(); } }, x))) : null,
      nbr[id].length ? h("div", {}, h("h4", {}, "Verbunden mit"),
        nbr[id].map((r) => h("button", { class: "hop", onclick: () => open(r.id) },
          h("small", {}, (r.out ? "→ " : "← ") + shortType(r.type)), h("span", {}, h("b", {}, byId[r.id].t), " ", kindTag(byId[r.id].k))))) : null,
      n.src && n.src.length ? h("div", {}, h("h4", {}, "Quellen im Repository"), n.src.map((x) => h("div", { class: "src" }, x))) : null,
      h("h4", {}, "Weiterdenken"),
      h("div", { class: "bar" },
        h("button", { class: "btn", onclick: () => { addNote("Zu „" + n.t + "“: "); close(); go("ideen"); } }, "Notiz dazu"),
        h("button", { class: "btn", onclick: () => { state.atlasQ = (n.g && n.g[0]) || n.t.split(" ")[0]; close(); go("atlas"); render(); } }, "Im Atlas nachschlagen"))].flat(Infinity).filter((x) => x != null && x !== false));
    drawer.scrollTop = 0;
  }
  function close() { drawer.hidden = true; }
  const shortType = (t) => ({ leitet: "führt zu", bedingt: "bedingt", seed: "speist", offen: "fehlt für", rueck: "rückwärts", gleich: "gleiche Zahl", karte: "verbunden", dockt: "dockt an", testet: "prüft", grenze: "begrenzt" }[t] || t);
  const cellSym = (c) => ({ j: "●", t: "◐", n: "✕", o: "·" }[c] || "·");
  const cellWord = (c) => ({ j: "liefert im Geltungsbereich", t: "teilweise oder bedingt", n: "scheitert dort", o: "nicht geprüft" }[c] || "nicht geprüft");
  function triple(n) {
    const col = (k, title) => h("div", { class: k }, h("h4", {}, title), (n[k] || []).length ? h("ul", {}, n[k].map((x) => h("li", {}, x))) : h("p", { class: "scope" }, "nichts vermerkt"));
    return h("div", { class: "tri" }, col("same", "Stimmt exakt überein"), col("diff", "Unterscheidet sich"), col("open", "Offen"));
  }

  /* ── 1 Überblick ── */
  function vStart() {
    const cnt = (st) => T.nodes.filter((n) => n.st === st).length;
    const F = A.problems.reduce((x, p) => x + p.f.length, 0), S = A.problems.reduce((x, p) => x + p.s.length, 0);
    const tile = (view, title, text) => h("a", { class: "card", href: "#/" + view }, h("h3", {}, title), h("p", {}, text));
    return [
      h("h1", {}, "Zwei Eingaben, eine Kette, eine offene Stelle"),
      h("p", { class: "lead" }, "TFPT setzt zwei Dinge voraus: eine Naht mit fester Normierung (P1) und einen Träger mit fünf Plätzen (P2). Daraus folgt Schritt für Schritt das Gerüst des Standardmodells samt Zahlenwerten. Was fehlt, ist die Regel, die das Ganze in Bewegung setzt."),
      h("div", { class: "steps" },
        [["1 · Eingaben", "Naht mit c₃ = 1/(8π), Träger mit 3 + 2 Plätzen"], ["2 · Geometrie", "Kugel mit vier Marken: Rang 3 und Dimension 5"],
         ["3 · Abschluss", "D₅ und A₃ verkleben zu E₈, die Uhren greifen ineinander"], ["4 · Auslesungen", "Ladungen, drei Familien, α, Mischungswinkel, Kosmologie"],
         ["5 · Offen", "Zeit, Ort, 3+1 Dimensionen und die gemeinsame Quelle"]].map(([a, b]) => h("div", {}, h("b", {}, a), b))),
      h("div", { class: "stats" },
        [[T.nodes.length, "Karten in dieser Landkarte"], [cnt("belegt"), "belegt im Geltungsbereich"], [cnt("bedingt"), "bedingt"],
         [cnt("aus"), "ausgeschlossen"], [cnt("offen") + cnt("neu"), "offen oder nicht probiert"], [F, "Befunde im Atlas"], [S, "geprüfte Strategien im Atlas"]]
          .map(([n, t]) => h("div", {}, h("b", {}, n), h("span", {}, t)))),
      h("h2", {}, "Wohin zuerst?"),
      h("div", { class: "grid" },
        tile("herleitung", "Herleitung", "Die ganze Kette von P1 und P2 bis zu Standardmodell, Vorhersagen und Rückbestimmung."),
        tile("dynamik", "Wo steckt die Dynamik?", "Alle Kandidaten nebeneinander: was jeder liefert, wo er scheitert, was nie geprüft wurde."),
        tile("muster", "Muster", "Qubits, Fehlerkorrektur, Primzahlen als Ereignisprotokoll, Ikosaeder und wiederkehrende Zahlen."),
        tile("ursprung", "Ursprung von P1 und P2", "Alle Ansätze zur Herkunft der beiden Eingaben und zum Universalraum."),
        tile("verbindungen", "Origin, OPH, andere Theorien", "Was exakt übereinstimmt, was sich unterscheidet, was offen ist."),
        tile("ideen", "Ideenraum", "Alles als Netz. Zahlen verfolgen, Brücken würfeln, Einfälle notieren."),
        tile("hypothesen", "Hypothesen", "Zwanzig Hypothesen in zwei Runden, alle nachgerechnet. Jede mit Test, Abbruchkriterium und Ergebnis."),
        tile("atlas", "Atlas", "Der vollständige Forschungsstand: jede Frage, jeder Befund, jede versuchte Strategie.")),
      h("h2", {}, "Wo es hängt"),
      h("div", { class: "grid" }, ["quelle", "gap-r113", "org-diagnose", "gap-lepton"].map((id) => card(byId[id]))),
      h("h2", {}, "So liest man die Farben"), legendStatus(),
      h("p", { class: "note" }, "Diese Seite ordnet ein und vereinfacht. Die Einstufungen sind meine Lesart des Korpus. Wo es darauf ankommt, zählt die genannte Quelle."),
    ];
  }

  /* ── 2 Herleitung ── */
  const COLS = ["Eingaben", "Geometrie", "Bausteine", "Abschluss", "Auslesungen", "Physik und Lücken"];
  function vHerleitung() {
    const W = 196, NW = 170, NH = 58, RH = 66, X0 = 12, Y0 = 34;
    const ch = T.nodes.filter((n) => n.k === "chain");
    const pos = (n) => ({ x: X0 + n.col * W, y: Y0 + n.row * RH });
    const H = Y0 + (Math.max(...ch.map((n) => n.row)) + 1) * RH + 10;
    const svg = s("svg", { class: "flow", viewBox: `0 0 ${X0 * 2 + 6 * W - (W - NW)} ${H}`, role: "img", "aria-label": "Herleitungskette" });
    COLS.forEach((c, i) => svg.append(s("text", { class: "col-h", x: X0 + i * W, y: 16 }, c)));
    const fwd = new Set(["leitet", "bedingt", "seed"]);
    const E = edges.filter((e) => byId[e.a].k === "chain" && byId[e.b].k === "chain");
    const eEls = E.map((e) => {
      const a = pos(byId[e.a]), b = pos(byId[e.b]);
      let d;
      if (byId[e.b].col > byId[e.a].col) { const x1 = a.x + NW, y1 = a.y + NH / 2, x2 = b.x, y2 = b.y + NH / 2, m = (x2 - x1) / 2; d = `M${x1},${y1} C${x1 + m},${y1} ${x2 - m},${y2} ${x2},${y2}`; }
      else if (byId[e.b].col === byId[e.a].col) { const x = a.x + NW, y1 = a.y + NH / 2, y2 = b.y + NH / 2; d = `M${x},${y1} C${x + 22},${y1} ${x + 22},${y2} ${x},${y2}`; }
      else { const x1 = a.x, y1 = a.y + NH - 6, x2 = b.x + NW, y2 = b.y + NH - 6, dip = 26; d = `M${x1},${y1} C${x1 - 40},${y1 + dip} ${x2 + 40},${y2 + dip} ${x2},${y2}`; }
      const el = s("path", { class: "fedge t-" + e.type, d }, s("title", {}, byId[e.a].t + " " + T.linkTypes[e.type] + " " + byId[e.b].t));
      svg.append(el);
      return el;
    });
    const nEls = {};
    ch.forEach((n) => {
      const p = pos(n);
      const fo = s("foreignObject", { class: "fnode st-" + n.st, x: p.x, y: p.y, width: NW, height: NH, tabindex: 0 });
      fo.append(h("div", {}, h("b", {}, n.t), h("small", {}, T.status[n.st][0])));
      fo.addEventListener("click", () => open(n.id));
      fo.addEventListener("keydown", (ev) => ev.key === "Enter" && open(n.id));
      fo.addEventListener("mouseenter", () => trace(n.id));
      fo.addEventListener("mouseleave", paint);
      svg.append(fo);
      nEls[n.id] = fo;
    });
    function reach(id, dir) {
      const seen = new Set([id]), q = [id];
      while (q.length) { const c = q.pop(); E.forEach((e) => { if (!fwd.has(e.type)) return; const [f, t] = dir ? [e.a, e.b] : [e.b, e.a]; if (f === c && !seen.has(t)) { seen.add(t); q.push(t); } }); }
      return seen;
    }
    function trace(id) {
      const keep = new Set([...reach(id, true), ...reach(id, false)]);
      ch.forEach((n) => nEls[n.id].classList.toggle("dim", !keep.has(n.id)));
      E.forEach((e, i) => { const on = keep.has(e.a) && keep.has(e.b) && fwd.has(e.type); eEls[i].classList.toggle("dim", !on); eEls[i].classList.toggle("hot", on); });
    }
    function paint() {
      const m = state.mode;
      ch.forEach((n) => nEls[n.id].classList.toggle("dim", m === "luecke" && !["offen", "gesetzt"].includes(n.st)));
      E.forEach((e, i) => {
        const on = m === "vor" ? fwd.has(e.type) : m === "rueck" ? e.type === "rueck" : e.type === "offen";
        eEls[i].classList.toggle("hot", m !== "vor" && on); eEls[i].classList.toggle("dim", !on);
      });
    }
    const modes = [["vor", "Vorwärts lesen"], ["rueck", "Rückwärts bestimmen"], ["luecke", "Wo es fehlt"]];
    const hint = { vor: "Von links nach rechts: Jede Karte folgt aus denen links von ihr. Mit der Maus über eine Karte fahren zeigt alles, wovon sie abhängt und was aus ihr folgt.",
      rueck: "Der Abschluss schaut zurück: Verlangt man, dass am Ende alles zusammenpasst, bleiben nur g = 5 und die Vier übrig. Das zeigt, dass die Eingaben stimmig sind. Es erklärt nicht, woher sie kommen.",
      luecke: "Gestrichelte Karten sind offen, graue sind gesetzt. Die roten Linien zeigen, was für den Ursprung fehlt: eine Quelle, die P1 und P2 selbst erzeugt." };
    const bar = h("div", { class: "bar" }, h("span", { class: "lbl" }, "Lesart"), modes.map(([k, t]) => h("button", { class: "btn" + (state.mode === k ? " on" : ""), onclick: () => { state.mode = k; render(); } }, t)));
    paint();
    return [
      h("h1", {}, "Die Herleitung"),
      h("p", { class: "lead" }, "Zwei Eingaben links, die Physik rechts. Jede Karte lässt sich antippen."),
      bar, h("p", { class: "note" }, hint[state.mode]), h("div", { class: "scroll" }, svg),
      h("div", { class: "legend" }, [["leitet", "folgt"], ["bedingt", "folgt unter Annahmen"], ["seed", "speist als Zahl"], ["rueck", "bestimmt rückwärts"], ["offen", "fehlt"]]
        .map(([k, t]) => h("span", {}, h("i", { style: `border-color: var(${{ leitet: "--ink", bedingt: "--s-bedingt", seed: "--k-pat", rueck: "--k-ur", offen: "--s-offen" }[k]}); border-top-style: ${k === "leitet" || k === "seed" ? "solid" : "dashed"}` }), t))),
      legendStatus(),
      h("h2", {}, "Was die Kette ausgibt"),
      h("p", {}, "Eingefrorene Werte vom ", (A.frozen && A.frozen.freeze_date) || "9. Juni 2026", ". Der Balken zeigt den Abstand zur Messung in Standardabweichungen. Innerhalb des hellen Bandes liegt ein Wert bei höchstens 2 σ."),
      h("div", { class: "scroll" }, outputsTable()),
      h("p", { class: "note warn" }, "Fast alle Zeilen lesen denselben Seed φ₀. Ihre Treffer hängen zusammen und zählen nicht als unabhängige Bestätigungen. Die Messwerte stammen aus dem Lauf vom 18. November 2026 und wurden dort überwiegend aus dem Gedächtnis zitiert."),
    ];
  }
  function outputsTable() {
    const dot = (o) => {
      if (o.pull == null) return h("span", { class: "scope" }, o.pct ? o.pct + " zur Polmasse" : "noch kein Messvergleich");
      const w = 150, x = (v) => w / 2 + Math.max(-3, Math.min(3, v)) * (w / 6.4);
      return h("span", { style: "white-space:nowrap" }, s("svg", { width: w, height: 18, viewBox: `0 0 ${w} 18`, role: "img", "aria-label": o.pull + " Sigma" },
        s("rect", { x: x(-2), y: 3, width: x(2) - x(-2), height: 12, fill: "#eef2f5" }), s("rect", { x: x(-1), y: 3, width: x(1) - x(-1), height: 12, fill: "#dfe6ec" }),
        s("line", { x1: w / 2, y1: 1, x2: w / 2, y2: 17, stroke: "#b3c1cf" }),
        s("circle", { cx: x(o.pull), cy: 9, r: 4.5, fill: Math.abs(o.pull) > 2 ? "#a8402f" : "#20384f" })),
        " " + (o.pull > 0 ? "+" : "−") + Math.abs(o.pull).toFixed(2).replace(".", ",") + " σ");
    };
    return h("table", {}, h("thead", {}, h("tr", {}, ["Größe", "Formel", "TFPT", "Messung", "Abstand", "Ebene"].map((x) => h("th", {}, x)))),
      h("tbody", {}, T.outputs.map((o) => h("tr", {}, h("td", {}, h("b", {}, o.t), o.note && h("div", { class: "scope" }, o.note)),
        h("td", { class: "num" }, o.f), h("td", { class: "num" }, o.v), h("td", { class: "num" }, o.m), h("td", {}, dot(o)), h("td", {}, o.layer)))));
  }

  /* ── 3 Dynamik ── */
  function vDynamik() {
    let rows = T.nodes.filter((n) => n.k === "dyn");
    const stCount = {}; rows.forEach((n) => (stCount[n.st] = (stCount[n.st] || 0) + 1));
    const all = rows;
    if (state.dynSt) rows = rows.filter((n) => n.st === state.dynSt);
    const rank = { j: 0, t: 1, n: 2, o: 3 };
    const body = [];
    if (state.dynSort) rows = rows.slice().sort((a, b) => rank[a.req[state.dynSort]] - rank[b.req[state.dynSort]]);
    let grp = null;
    rows.forEach((n) => {
      if (!state.dynSort && n.grp !== grp) { grp = n.grp; body.push(h("tr", { class: "grp" }, h("td", { colspan: 9 }, grp))); }
      body.push(h("tr", { class: "hit", onclick: () => open(n.id), tabindex: 0, onkeydown: (e) => e.key === "Enter" && open(n.id) },
        h("td", {}, h("b", {}, n.t), h("small", {}, n.s)), h("td", {}, pill(n.st)),
        T.req.map(([k]) => h("td", { class: "cell c-" + n.req[k], title: cellWord(n.req[k]) }, cellSym(n.req[k])))));
    });
    body.push(h("tr", { class: "grp" }, h("td", { colspan: 2 }, "Liefert voll · von " + all.length), T.req.map(([k]) => h("td", { class: "cell" }, all.filter((n) => n.req[k] === "j").length))));
    const outs = { reuse_result: "wiederverwendbar", excluded_scoped: "ausgeschlossen", insufficient_for_goal: "reicht nicht fürs Ziel" };
    const strat = [];
    A.problems.forEach((p) => p.s.forEach((x) => strat.push({ p, x })));
    const q = state.stratQ.toLowerCase();
    const shown = strat.filter(({ p, x }) => (!state.strat || x.o === state.strat) && (!q || (x.t + " " + x.x + " " + p.title).toLowerCase().includes(q)));
    return [
      h("h1", {}, "Wo steckt die Dynamik?"),
      h("p", { class: "lead" }, "Die Kette beschreibt, was ist. Gesucht ist, was sie in Bewegung setzt. Hier stehen alle Kandidaten nebeneinander, gemessen an sieben Dingen, die eine echte Dynamik liefern muss."),
      h("div", { class: "bar" }, h("span", { class: "lbl" }, "Zeigen"),
        h("button", { class: "btn" + (!state.dynSt ? " on" : ""), onclick: () => { state.dynSt = null; render(); } }, "alle " + all.length),
        Object.keys(stCount).map((k) => h("button", { class: "btn" + (state.dynSt === k ? " on" : ""), onclick: () => { state.dynSt = k; render(); } }, T.status[k][0] + " " + stCount[k])),
        state.dynSort && h("button", { class: "btn", onclick: () => { state.dynSort = null; render(); } }, "Sortierung aufheben")),
      h("div", { class: "scroll" }, h("table", { class: "mtx" },
        h("thead", {}, h("tr", {}, h("th", {}, "Ansatz"), h("th", {}, "Status"),
          T.req.map(([k, name, tip]) => h("th", { class: "cell", title: tip }, h("button", { onclick: () => { state.dynSort = k; render(); }, title: tip + " Antippen sortiert danach." }, name))))),
        h("tbody", {}, body))),
      h("div", { class: "legend" }, [["j", "liefert im Geltungsbereich"], ["t", "teilweise oder bedingt"], ["n", "scheitert dort"], ["o", "nicht geprüft"]]
        .map(([c, t]) => h("span", {}, h("b", { class: "c-" + c }, cellSym(c)), " ", t))),
      h("p", { class: "note warn" }, "Die Punkte in der Tabelle sind meine Zusammenfassung der Verträge, keine Ledgereinträge. Sie zeigen auf einen Blick: Kein Ansatz liefert Ort und 3+1 Dimensionen voll, und fast keiner kommt ohne Zusatzwahl aus. Dort liegt die Lücke."),
      h("h2", {}, "Was noch niemand ausgeführt hat"),
      h("div", { class: "grid" }, T.nodes.filter((n) => n.k === "gap").map(card)),
      h("h2", {}, "Alle " + strat.length + " geprüften Strategien aus dem Atlas"),
      h("p", {}, "Vollständige Liste dessen, was versucht wurde, mit Ausgang und dem, was sich ändern müsste, damit ein neuer Versuch lohnt."),
      h("div", { class: "bar" },
        h("button", { class: "btn" + (!state.strat ? " on" : ""), onclick: () => { state.strat = ""; render(); } }, "alle"),
        Object.keys(outs).map((k) => h("button", { class: "btn" + (state.strat === k ? " on" : ""), onclick: () => { state.strat = k; render(); } }, outs[k] + " " + strat.filter((z) => z.x.o === k).length)),
        h("input", { class: "btn", type: "search", placeholder: "Strategien durchsuchen", value: state.stratQ, "aria-label": "Strategien durchsuchen",
          onchange: (e) => { state.stratQ = e.target.value; render(); } })),
      h("p", { class: "scope" }, shown.length + " Treffer" + (shown.length > 80 ? ", die ersten 80 werden gezeigt. Suche eingrenzen für mehr." : "")),
      shown.slice(0, 80).map(({ p, x }) => h("div", { class: "rec" }, h("b", {}, h("span", { class: "tag a-" + x.o }, outs[x.o] || x.o), x.t),
        h("p", {}, x.x), x.c && h("p", { class: "scope" }, "Müsste sich ändern: " + x.c), h("p", { class: "scope" }, "Frage: " + p.title))),
    ];
  }

  /* ── 4 Muster ── */
  function groups(list, key) { const g = {}; list.forEach((n) => (g[n[key]] = g[n[key]] || []).push(n)); return g; }
  function vMuster() {
    const g = groups(T.nodes.filter((n) => n.k === "pat"), "grp");
    const pats = (A.patterns || []);
    return [
      h("h1", {}, "Muster und Verbindungen"),
      h("p", { class: "lead" }, "Strukturen, die an mehreren Stellen auftauchen. Manche sind bewiesene Gleichheiten, manche nur dieselbe Zahl an zwei Orten, manche sind Warnungen."),
      h("p", { class: "note warn" }, "Faustregel: Gleiche Zahl, gleiches Spektrum oder gleiche Gruppenordnung liefern noch keine Abbildung zwischen den Dingen. Ein Muster wird erst zum Ergebnis, wenn jemand die Abbildung hinschreibt."),
      Object.keys(g).map((k) => [h("h2", {}, k), h("div", { class: "grid" }, g[k].map(card))]),
      pats.length ? [h("h2", {}, "Musterregister des Atlas"), h("p", {}, "Die " + pats.length + " im Atlas geführten Muster mit ihren Bedingungen und Gegenbelegen."),
        pats.map((p) => h("details", { class: "prob" }, h("summary", {}, h("b", {}, p.label), h("span", {}, p.kind)),
          h("div", { class: "in" }, h("p", {}, p.summary),
            p.conditions.length ? h("p", { class: "scope" }, "Bedingungen: " + p.conditions.map(txt).join(" · ")) : null,
            p.counter.length ? h("p", { class: "scope" }, "Gegenbelege: " + p.counter.map(txt).join(" · ")) : null)))] : null,
    ];
  }
  const txt = (x) => (typeof x === "string" ? x : x.de || x.label_de || x.text || x.id || JSON.stringify(x));

  /* ── 5 Ursprung ── */
  function vUrsprung() {
    const org = T.nodes.filter((n) => n.k === "org");
    const lanes = [["p1", "Woher kommt P1?", "Naht, Normierung c₃, Windung"], ["both", "Beide zugleich", "Ein Objekt, zwei Lesungen"],
      ["p2", "Woher kommt P2?", "Fünf Plätze, 3 + 2"], ["sel", "Was wählt aus?", "Marken, Modul, Suchverfahren"]];
    return [
      h("h1", {}, "Ursprung von P1 und P2"),
      h("p", { class: "lead" }, "Die Kette beginnt mit zwei Setzungen. Hier stehen alle Versuche, sie aus etwas Tieferem zu gewinnen."),
      h("div", { class: "grid" }, card(byId.p1), card(byId.p2), card(byId["org-diagnose"])),
      h("h2", {}, "Die Ansätze"),
      h("div", { class: "lanes" }, lanes.map(([k, t, sub]) => h("div", { class: "lane" }, h("h3", {}, t), h("small", {}, sub), org.filter((n) => n.lane === k).map(card)))),
      h("h2", {}, "Universalraum"),
      h("p", {}, "Der Universalraum ist der Versuch, den größten Raum zu benennen, in dem TFPT ausgeführt wird oder aus dem es herausfällt."),
      h("div", { class: "grid" }, T.nodes.filter((n) => n.k === "ur").map(card)),
      h("h2", {}, "Dynamikkandidaten, die am Ursprung ansetzen"),
      h("div", { class: "grid" }, ["dyn-holo", "dyn-reader", "dyn-openword", "dyn-selfdual", "dyn-lorentz", "dyn-reset"].map((id) => card(byId[id]))),
      legendStatus(),
    ];
  }

  /* ── 6 Verbindungen ── */
  function vVerbindungen() {
    const g = groups(T.nodes.filter((n) => n.k === "ext"), "sub");
    const intro = { "Origin Idee": "Die physische Erzählung zum Compiler: Was bedeutet die Naht in der Welt?",
      "OPH": "Ein fremdes Programm mit anderem Startpunkt. Der Vergleich zeigt, welche Teile von TFPT unabhängig wiedergefunden werden.",
      "Andere Theorien": "Wo TFPT bekannte Ideen berührt, und wo es sich bewusst abgrenzt." };
    return [
      h("h1", {}, "Verbindungen nach außen"),
      h("p", { class: "lead" }, "Jede Verbindung in drei Spalten: was exakt gleich ist, was sich unterscheidet, was offen bleibt."),
      Object.keys(g).map((k) => [h("h2", {}, k), h("p", {}, intro[k]),
        g[k].map((n) => h("div", { class: "ext" }, h("header", {}, h("button", { onclick: () => open(n.id) }, h("h3", {}, n.t)), pill(n.st)), h("p", {}, n.s), triple(n)))]),
    ];
  }

  /* ── 7 Ideenraum ── */
  let layout = null;
  function computeLayout() {
    const W = 1000, H = 720, ns = T.nodes, kinds = Object.keys(T.kinds);
    const P = ns.map((n, i) => { const a = (kinds.indexOf(n.k) / kinds.length) * 6.283 + (i % 7) * 0.11, r = 200 + ((i * 37) % 120); return { x: W / 2 + r * Math.cos(a), y: H / 2 + r * Math.sin(a), vx: 0, vy: 0 }; });
    const idx = Object.fromEntries(ns.map((n, i) => [n.id, i]));
    for (let it = 0; it < 380; it++) {
      const cool = 1 - it / 380;
      for (let i = 0; i < P.length; i++) for (let j = i + 1; j < P.length; j++) {
        let dx = P[i].x - P[j].x, dy = P[i].y - P[j].y, d2 = dx * dx + dy * dy + 0.01, f = 9000 / d2; const d = Math.sqrt(d2);
        dx /= d; dy /= d; P[i].vx += dx * f; P[i].vy += dy * f; P[j].vx -= dx * f; P[j].vy -= dy * f;
      }
      edges.forEach((e) => { const a = P[idx[e.a]], b = P[idx[e.b]]; const dx = b.x - a.x, dy = b.y - a.y, d = Math.sqrt(dx * dx + dy * dy) + 0.01, f = (d - 120) * 0.03;
        a.vx += (dx / d) * f; a.vy += (dy / d) * f; b.vx -= (dx / d) * f; b.vy -= (dy / d) * f; });
      P.forEach((p) => { p.vx += (W / 2 - p.x) * 0.005; p.vy += (H / 2 - p.y) * 0.008; const sp = Math.min(1, 14 * cool / (Math.hypot(p.vx, p.vy) + 0.01));
        p.x += p.vx * sp; p.y += p.vy * sp; p.vx = 0; p.vy = 0; });
    }
    const xs = P.map((p) => p.x), ys = P.map((p) => p.y), x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
    P.forEach((p) => { p.x = 40 + ((p.x - x0) / (x1 - x0)) * (W - 150); p.y = 24 + ((p.y - y0) / (y1 - y0)) * (H - 48); });
    return { W, H, P: Object.fromEntries(ns.map((n, i) => [n.id, P[i]])) };
  }
  const notes = { key: "tfpt-landkarte-notizen",
    all() { try { return JSON.parse(localStorage.getItem(this.key) || "[]"); } catch (e) { return []; } },
    save(a) { try { localStorage.setItem(this.key, JSON.stringify(a)); } catch (e) { /* ohne Speicher weiter */ } } };
  let draft = "";
  function addNote(text) { draft = text; }
  function vIdeen() {
    layout = layout || computeLayout();
    const { W, H, P } = layout;
    const match = (n) => (state.num ? (n.n || []).includes(state.num) : state.tag ? (n.g || []).includes(state.tag) : false);
    const filt = state.num || state.tag;
    const svg = s("svg", { class: "graph", viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Netz aller Karten" });
    const vis = (n) => state.kinds.has(n.k);
    const eEls = edges.filter((e) => vis(byId[e.a]) && vis(byId[e.b])).map((e) => {
      const el = s("line", { class: "gedge t-" + e.type, x1: P[e.a].x, y1: P[e.a].y, x2: P[e.b].x, y2: P[e.b].y }); svg.append(el); return { e, el };
    });
    const nEls = {};
    T.nodes.filter(vis).forEach((n) => {
      const deg = nbr[n.id].length, g = s("g", { class: "gnode", transform: `translate(${P[n.id].x},${P[n.id].y})`, tabindex: 0 },
        s("circle", { r: 5 + Math.min(9, deg * 0.9), fill: `var(--k-${n.k})`, opacity: n.st === "aus" ? 0.45 : 1 }),
        s("text", { x: 10 + Math.min(9, deg * 0.9), y: 4 }, n.t.length > 34 ? n.t.slice(0, 32) + "…" : n.t), s("title", {}, n.t + " · " + T.status[n.st][0]));
      g.addEventListener("click", () => open(n.id));
      g.addEventListener("keydown", (ev) => ev.key === "Enter" && open(n.id));
      g.addEventListener("mouseenter", () => focus(n.id)); g.addEventListener("mouseleave", () => focus(null));
      svg.append(g); nEls[n.id] = g;
    });
    function focus(id) {
      const keep = id ? new Set([id, ...nbr[id].map((r) => r.id)]) : null;
      Object.keys(nEls).forEach((k) => {
        const on = keep ? keep.has(k) : filt ? match(byId[k]) : true;
        nEls[k].classList.toggle("dim", !on);
        nEls[k].querySelector("text").style.display = keep ? (keep.has(k) ? "" : "none") : filt ? (on ? "" : "none") : nbr[k].length >= 8 ? "" : "none";
      });
      eEls.forEach(({ e, el }) => { const on = id ? e.a === id || e.b === id : false; el.classList.toggle("hot", on); el.classList.toggle("dim", id ? !on : !!filt); });
    }
    focus(null);
    const count = (key) => { const c = {}; T.nodes.forEach((n) => (n[key] || []).forEach((x) => (c[x] = (c[x] || 0) + 1))); return Object.entries(c).filter((x) => x[1] >= 2).sort((a, b) => b[1] - a[1]); };
    const chips = (list, cur, set) => list.map(([x, c]) => h("button", { class: "chip" + (cur === x ? " on" : ""), onclick: () => { set(cur === x ? null : x); render(); } }, x, h("sup", {}, c)));
    const hits = filt ? T.nodes.filter(match) : [];
    const N = notes.all();
    const ta = h("textarea", { placeholder: "Einfall, Frage oder Brücke festhalten", "aria-label": "Notiz" }); ta.value = draft;
    ta.addEventListener("input", () => (draft = ta.value));
    return [
      h("h1", {}, "Ideenraum"),
      h("p", { class: "lead" }, "Alles auf einer Fläche. Eine Zahl oder ein Stichwort antippen zeigt, wo es sonst noch vorkommt. Daraus entstehen Fragen."),
      h("div", { class: "bar" }, h("span", { class: "lbl" }, "Schichten"), Object.keys(T.kinds).map((k) => h("button", { class: "btn" + (state.kinds.has(k) ? " on" : ""),
        onclick: () => { state.kinds.has(k) ? state.kinds.delete(k) : state.kinds.add(k); render(); } }, h("span", { class: "k-" + k, style: state.kinds.has(k) ? "color:#fff" : "" }, "● "), T.kinds[k]))),
      h("div", { class: "ideen" },
        h("div", {}, svg,
          filt ? h("div", { class: "panel", style: "margin-top:12px" }, h("h3", {}, (state.num ? "Die Zahl " : "Stichwort ") + filt + " an " + hits.length + " Stellen"),
            h("p", { class: "note warn" }, "Gleiche Zahl heißt noch nicht gleiche Sache. Die spannende Frage ist jeweils: Gibt es eine Abbildung, die das eine ins andere überführt?"),
            h("div", { class: "grid" }, hits.map(card))) : null),
        h("div", { class: "side" },
          h("div", { class: "panel" }, h("h3", {}, "Zahlen, die wiederkehren"), chips(count("n"), state.num, (x) => { state.num = x; state.tag = null; })),
          h("div", { class: "panel" }, h("h3", {}, "Stichworte"), chips(count("g"), state.tag, (x) => { state.tag = x; state.num = null; })),
          h("div", { class: "panel", id: "bridge" }, bridgePanel()),
          h("div", { class: "panel" }, h("h3", {}, "Notizen"), ta,
            h("div", { class: "bar" },
              h("button", { class: "btn go", onclick: () => { if (!draft.trim()) return; const a = notes.all(); a.unshift({ ts: new Date().toISOString().slice(0, 16).replace("T", " "), text: draft.trim() }); notes.save(a); draft = ""; render(); } }, "Festhalten"),
              N.length ? h("button", { class: "btn", onclick: exportNotes }, "Als Markdown sichern") : null),
            h("ul", { class: "plain notes" }, N.map((n, i) => h("li", {}, h("small", {}, n.ts + " "), h("button", { title: "Notiz löschen", onclick: () => { const a = notes.all(); a.splice(i, 1); notes.save(a); render(); } }, "entfernen"), h("div", {}, n.text)))),
            h("p", { class: "scope" }, "Notizen bleiben in diesem Browser gespeichert.")))),
    ];
  }
  let bridge = null;
  function rollBridge() {
    const linked = (a, b) => nbr[a.id].some((r) => r.id === b.id);
    const ns = T.nodes.filter((n) => n.st !== "aus" || Math.random() < 0.3);
    for (let i = 0; i < 400; i++) {
      const a = ns[Math.floor(Math.random() * ns.length)], b = ns[Math.floor(Math.random() * ns.length)];
      if (a === b || a.k === b.k || linked(a, b)) continue;
      const sn = (a.n || []).filter((x) => (b.n || []).includes(x)), sg = (a.g || []).filter((x) => (b.g || []).includes(x));
      if (sn.length + sg.length) { bridge = { a, b, sn, sg }; return; }
    }
  }
  function bridgePanel() {
    const out = [h("h3", {}, "Brücke würfeln"), h("p", { class: "scope" }, "Zwei Karten, die etwas teilen und noch nicht verbunden sind.")];
    if (bridge) {
      const { a, b, sn, sg } = bridge, shared = [...sn, ...sg].join(", ");
      out.push(h("div", { class: "bridge" }, card(a), h("div", { class: "mid" }, "↔"), card(b)),
        h("p", {}, "Gemeinsam: ", h("b", {}, shared)),
        h("ul", {}, ["Welche Abbildung müsste es geben, damit beide dasselbe Objekt meinen?", "Tragen beide dieselbe Gruppenwirkung, denselben Zustand, dieselbe Normierung?",
          "Welche Rechnung würde die Brücke widerlegen?"].map((x) => h("li", {}, x))),
        h("button", { class: "btn", onclick: () => { addNote("Brücke: " + a.t + " ↔ " + b.t + " (gemeinsam: " + shared + "). "); render(); } }, "Als Notiz übernehmen"));
    }
    out.push(h("div", { class: "bar" }, h("button", { class: "btn go", onclick: () => { rollBridge(); document.getElementById("bridge").replaceChildren(...bridgePanel()); } }, bridge ? "Neu würfeln" : "Würfeln")));
    return out;
  }
  function exportNotes() {
    const md = "# TFPT Landkarte · Notizen\n\n" + notes.all().map((n) => "## " + n.ts + "\n\n" + n.text + "\n").join("\n");
    const a = h("a", { href: URL.createObjectURL(new Blob([md], { type: "text/markdown" })), download: "tfpt-notizen.md" });
    document.body.append(a); a.click(); a.remove();
  }

  /* ── 8 Atlas ── */
  const AL = { established_scoped: "belegt im Bereich", conditional: "bedingt", excluded_scoped: "ausgeschlossen im Bereich", stale_source: "Quelle veraltet",
    reuse_result: "wiederverwendbar", insufficient_for_goal: "reicht nicht fürs Ziel", resolved: "erledigt", open: "offen", ready: "bereit", blocked: "blockiert" };
  function vAtlas() {
    const q = state.atlasQ.toLowerCase();
    const hit = (p) => !q || (p.title + " " + p.q + " " + p.topics.join(" ")).toLowerCase().includes(q) || p.f.some((f) => (f.t + " " + f.x).toLowerCase().includes(q)) || p.s.some((x) => (x.t + " " + x.x).toLowerCase().includes(q));
    const list = A.problems.filter(hit);
    const rec = (title, st, text, extra) => h("div", { class: "rec" }, h("b", {}, h("span", { class: "tag a-" + st }, AL[st] || st), title), text && h("p", {}, text), extra);
    const inner = (p) => {
      const sub = (arr, fn) => (q ? arr.filter((x) => (x.t + " " + (x.x || x.n || "")).toLowerCase().includes(q)) : arr).concat().slice(0, 400).map(fn);
      const fq = sub(p.f, (f) => rec(f.t, f.s, f.x, [f.g && h("p", { class: "scope" }, "Geltungsbereich: " + f.g), f.p.length ? h("p", { class: "scope src" }, f.p.join("  ·  ")) : null]));
      const sq = sub(p.s, (x) => rec(x.t, x.o, x.x, x.c && h("p", { class: "scope" }, "Müsste sich ändern: " + x.c)));
      const kq = sub(p.k, (k) => rec(k.t, k.s, k.n));
      return h("div", { class: "in" }, h("p", {}, p.q), h("p", { class: "scope src" }, p.id),
        sq.length ? [h("h3", {}, "Versuchte Strategien"), sq] : null, kq.length ? [h("h3", {}, "Aufgaben"), kq] : null, fq.length ? [h("h3", {}, "Befunde"), fq] : null,
        q && !fq.length && !sq.length && !kq.length ? h("p", { class: "scope" }, "Treffer nur in Titel oder Frage.") : null);
    };
    return [
      h("h1", {}, "Atlas"),
      h("p", { class: "lead" }, "Der vollständige Stand aus dem Forschungsprogramm: " + A.problems.length + " Fragen mit allen Befunden, Strategien und Aufgaben. Unverändert übernommen."),
      h("div", { class: "bar" }, h("input", { class: "btn", style: "min-width:min(420px,100%)", type: "search", placeholder: "Atlas durchsuchen, zum Beispiel KMS, Higgs, Windung", value: state.atlasQ, "aria-label": "Atlas durchsuchen",
        onchange: (e) => { state.atlasQ = e.target.value; render(); } }), h("span", { class: "lbl" }, list.length + " Fragen")),
      A.problems.length ? null : h("p", { class: "note warn" }, "Keine Atlasdaten geladen. Mit python3 -B tfpt_landkarte/build_data.py erzeugen."),
      list.map((p) => { const d = h("details", { class: "prob" }, h("summary", {}, h("b", {}, p.title), h("span", {}, p.f.length + " Befunde · " + p.s.length + " Strategien · " + p.k.filter((k) => k.s !== "resolved").length + " offen")));
        d.addEventListener("toggle", () => { if (d.open && d.children.length < 2) d.append(inner(p)); }); return d; }),
    ];
  }

  /* ── Hypothesen ── */
  function hypBody(n) {
    const names = Object.fromEntries(T.req.map(([k, name]) => [k, name]));
    return [n.res && h("div", { class: "verdict st-" + n.st }, h("h4", {}, "Ergebnis der Prüfung · " + n.chk + " Rechnungen"), h("p", {}, h("b", {}, n.res)),
        n.found.length ? [h("h4", {}, "Trägt"), h("ul", {}, n.found.map((x) => h("li", {}, x)))] : null,
        n.falls.length ? [h("h4", {}, "Fällt"), h("ul", {}, n.falls.map((x) => h("li", {}, x)))] : null,
        n.opn ? [h("h4", {}, "Offen"), h("p", {}, n.opn)] : null),
      n.pre && [h("h4", {}, "Schon gerechnet"), h("p", {}, n.pre)],
      h("h4", {}, "Woraus sie folgt"), h("ul", {}, n.why.map((x) => h("li", {}, x))),
      h("h4", {}, "Würde liefern"), n.gives.map((k) => h("span", { class: "chip" }, names[k])),
      h("h4", {}, "Entscheidender Test"), h("p", {}, n.test), h("h4", {}, "Woran sie scheitert"), h("p", {}, n.kill)];
  }
  function hypCard(n) {
    return h("div", { class: "ext" },
      h("header", {}, h("button", { onclick: () => open(n.id) }, h("h3", {}, n.nr + " · " + n.t)), n.gives.map((k) => h("span", { class: "chip" }, T.req.find((r) => r[0] === k)[1]))),
      h("p", {}, n.s),
      n.res && h("p", { class: "verdictline" }, h("span", { class: "st st-" + n.st }, "● " + T.status[n.st][0]), " ", h("b", {}, n.res), " · " + n.chk),
      n.pre && h("p", { class: "verdictline" }, h("span", { class: "st st-neu" }, "Schon gerechnet: "), n.pre),
      h("div", { class: "tri" },
        h("div", { class: "same" }, h("h4", {}, "Woraus sie folgt"), h("ul", {}, n.why.map((x) => h("li", {}, x)))),
        h("div", { class: "open" }, h("h4", {}, "Entscheidender Test"), h("p", {}, n.test)),
        h("div", { class: "diff" }, h("h4", {}, "Woran sie scheitert"), h("p", {}, n.kill))));
  }
  function vHypothesen() {
    const hy = T.nodes.filter((n) => n.k === "hyp").sort((a, b) => a.nr - b.nr);
    const r2 = hy.filter((n) => n.runde === 2), r1 = hy.filter((n) => n.runde !== 2);
    return [
      h("h1", {}, "Hypothesen"),
      h("p", { class: "lead" }, "Zwei Runden mit je zehn Hypothesen, beide nachgerechnet. Keine löst die offene Frage. Unter jeder steht, was trägt, was fällt und was offen bleibt."),
      h("h2", {}, "Runde 2 · zehn Hypothesen, geprüft"),
      h("p", { class: "note warn" }, "Am 19. November 2026 nachgerechnet (111 Rechnungen, Vertrag HYPCHECKB20261119). Vier fallen, eine steht schon im Original, drei tragen im Kern ohne Bindung an das Original, eine trägt halb, eine ist ohne fremde Daten nicht entscheidbar. Nichts davon steht im Ledger oder im Atlas."),
      r2.map(hypCard),
      h("h2", {}, "Parallele Prüfung · andere Sitzung"),
      h("p", { class: "note" }, "Zehn andere Hypothesen, in einer anderen Sitzung geprüft. Die Unterlagen liegen nicht im Repository. " + (T.parallelNote || "")),
      h("div", { class: "ext" }, h("ul", { class: "par" }, T.parallel.map(([t, r, to]) => h("li", {}, h("b", {}, t + ". "), r, to.length ? h("span", { class: "soft" }, " → Hypothese " + to.join(", ")) : null)))),
      h("h2", {}, "Runde 1 · zehn Hypothesen, geprüft"),
      h("p", { class: "note warn" }, "Alle zehn sind am 19. November 2026 nachgerechnet (110 Rechnungen, Vertrag HYPCHECK20261119). Keine löst die offene Frage. Fünf fallen, bei einer fällt die Formel und eine exakte ersetzt sie, vier tragen in engerer Fassung. Fünf Befunde stehen im Atlas, nichts im Ledger."),
      r1.map(hypCard),
    ];
  }

  /* ── Navigation, Suche ── */
  const VIEWS = [["start", "Überblick", vStart], ["herleitung", "Herleitung", vHerleitung], ["dynamik", "Dynamik", vDynamik], ["muster", "Muster", vMuster],
    ["ursprung", "Ursprung P1 · P2", vUrsprung], ["verbindungen", "Origin · OPH · andere", vVerbindungen], ["hypothesen", "Hypothesen", vHypothesen], ["ideen", "Ideenraum", vIdeen], ["atlas", "Atlas", vAtlas]];
  function render() {
    const cur = (location.hash.replace("#/", "") || "start");
    const v = VIEWS.find((x) => x[0] === cur) || VIEWS[0];
    document.getElementById("nav").replaceChildren(...VIEWS.map(([k, t]) => h("a", { href: "#/" + k, class: k === v[0] ? "on" : "" }, t)));
    main.replaceChildren(...[v[2]()].flat(Infinity).filter((x) => x != null && x !== false));
    document.title = v[1] + " · TFPT Landkarte";
  }
  window.addEventListener("hashchange", () => { close(); render(); window.scrollTo(0, 0); });
  document.addEventListener("keydown", (e) => e.key === "Escape" && (close(), (hitsEl.hidden = true)));
  const qEl = document.getElementById("q"), hitsEl = document.getElementById("hits");
  qEl.addEventListener("input", () => {
    const q = qEl.value.trim().toLowerCase();
    if (q.length < 2) { hitsEl.hidden = true; return; }
    const ns = T.nodes.filter((n) => (n.t + " " + n.s + " " + (n.g || []).join(" ") + " " + (n.n || []).join(" ") + " " + (n.b || []).join(" ")).toLowerCase().includes(q)).slice(0, 10);
    const ps = A.problems.filter((p) => (p.title + " " + p.q).toLowerCase().includes(q)).slice(0, 5);
    hitsEl.replaceChildren(
      ...ns.map((n) => h("button", { onclick: () => { hitsEl.hidden = true; open(n.id); } }, n.t, h("small", {}, T.kinds[n.k] + " · " + T.status[n.st][0]))),
      ...ps.map((p) => h("button", { onclick: () => { hitsEl.hidden = true; state.atlasQ = qEl.value.trim(); go("atlas"); render(); } }, p.title, h("small", {}, "Atlas · " + p.f.length + " Befunde"))),
      h("button", { onclick: () => { hitsEl.hidden = true; state.atlasQ = qEl.value.trim(); go("atlas"); render(); } }, "Im ganzen Atlas nach „" + qEl.value.trim() + "“ suchen"));
    hitsEl.hidden = false;
  });
  document.addEventListener("click", (e) => { if (!e.target.closest(".search")) hitsEl.hidden = true; });
  document.getElementById("stand").textContent = A.generated ? " Atlasstand " + A.generated + ", Revision " + String(A.revision).slice(0, 8) + "." : "";
  render();
})();
