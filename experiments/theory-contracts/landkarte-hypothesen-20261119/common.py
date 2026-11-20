"""Gemeinsame Helfer. Exakte oder hochgenaue Arithmetik, keine angepassten Eingaben."""
import json, sys

class Checks:
    def __init__(self, name): self.name, self.rows = name, []
    def ok(self, label, cond, detail=""):
        self.rows.append({"check": label, "pass": bool(cond), "detail": str(detail)})
        print(("PASS " if cond else "FAIL ") + label + ("  | " + str(detail) if detail != "" else ""))
    def done(self, extra=None):
        n = sum(r["pass"] for r in self.rows)
        print(f"{self.name}: {n}/{len(self.rows)} checks passed")
        json.dump({"name": self.name, "checks": self.rows, "data": extra or {}},
                  open(f"results_{self.name}.json", "w"), indent=1, ensure_ascii=False, default=str)
        sys.exit(0 if n == len(self.rows) else 1)
