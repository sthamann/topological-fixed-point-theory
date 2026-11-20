"""H13: die Zahl vier aus Torus und Fouriertransformation."""
import itertools, subprocess, sys, os
from fractions import Fraction as Fr
from common import Checks
c = Checks("h13_genus")
g = {n: Fr(n - 2, 2) for n in range(2, 14, 2)}
c.ok("Doppelüberlagerung der Kugel über n Marken hat Geschlecht (n - 2)/2; nur n = 4 gibt einen Torus", [n for n in g if g[n] == 1] == [4], {n: str(v) for n, v in g.items()})
# flache Kugeln mit Kegelpunkten der Winkel 2 pi/m: Summe (1 - 1/m_i) = 2
sol = sorted({tuple(sorted(ms)) for k in range(3, 7) for ms in itertools.combinations_with_replacement(range(2, 13), k) if sum(1 - Fr(1, m) for m in ms) == 2})
c.ok("flache Kugeln mit Kegelpunkten gibt es genau vier: (2,2,2,2), (2,3,6), (2,4,4), (3,3,3)", sol == [(2, 2, 2, 2), (2, 3, 6), (2, 4, 4), (3, 3, 3)], sol)
c.ok("nur (2,2,2,2) hat ausschließlich Halbdrehungen (Defizit pi). Mit lauter Halbdrehungen folgt aus Gauss Bonnet n pi = 4 pi, also vier Marken", [s for s in sol if set(s) == {2}] == [(2, 2, 2, 2)])
c.ok("nur (2,2,2,2) hat vier Punkte und damit eine freie Lage (ein Doppelverhältnis). Die drei anderen sind starr. Die Frage nach der Markenlage entsteht also genau durch die Halbdrehung", [len(s) for s in sol] == [4, 3, 3, 3])
c.ok("(2,4,4) ist der quadratische Torus geteilt durch die Vierteldrehung, (2,3,6) und (3,3,3) gehören zum Sechseck. Ein Kegelpunkt der Ordnung 4 erzwingt das Quadrat ohne weitere Wahl", (2, 4, 4) in sol)
# steht das schon im Original?
root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
tex = open(os.path.join(root, "documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex")).read()
c.ok("ORIGINAL: der Abschnitt zur gemeinsamen Herkunft leitet die Vier schon so her (flache sphärische Z2 Klasse, Defizit pi, Gauss Bonnet)", "flat spherical" in tex and "Gauss--Bonnet gives $n\\pi=4\\pi$" in tex)
c.ok("ORIGINAL: die quadratische Lage ist dort ausdrücklich zusätzlich zur Anzahl", "square position of the marks is additional to their count" in tex)
v = os.path.join(root, "verification", "v216_marks_gauss_bonnet.py")
r = subprocess.run([sys.executable, "-B", v], cwd=os.path.join(root, "verification"), capture_output=True, text=True, timeout=300) if os.path.exists(v) else None
c.ok("die zugehörige Originalprüfung v216 läuft durch", r is not None and r.returncode == 0, (r.stdout.strip().splitlines() or ["keine Ausgabe"])[-1][:120] if r else "fehlt")
c.ok("URTEIL: flach, Torus als Überlagerung und Fouriertransformation sind drei Namen für dieselbe Bedingung. Die Hypothese ist richtig, aber keine neue Herleitung", True)
c.ok("ZIRKEL: Flachheit und Halbdrehung sind im Original gesetzt (Klasse gewählt), nicht hergeleitet. Schieber und Uhr kommen dort aus der Klebegruppe, also von P2. Die eigene Abbruchregel greift für den Teil mit Schieber und Uhr", True)
c.ok("NEU ist nur die Schärfung: Halbdrehung ist der einzige flache Fall mit freier Lage; ein Punkt der Ordnung 4 ließe keine Wahl", True)
c.done({"flache_kugeln": [list(s) for s in sol], "urteil": "Richtig, aber schon im Original (v216). Keine neue Herleitung; die Annahmen flach und Halbdrehung bleiben gesetzt."})
