"""A closed, recursive population of native fields and Kohonen memories.

Rows are cells. Every row retains all 35 complex symmetric quartic components.
K is the original signed native Gram matrix, not an energy operator.
The feedback law below is a chosen experiment, not a derived TFPT law.
"""
from dataclasses import asdict, dataclass
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def unit_rows(x):
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    # A zero response stays zero; no hidden reseeding or injected information.
    return x / np.where(norms > 1e-14, norms, 1)


@dataclass
class Config:
    seed: int = 7
    coupling: float = 0.5
    feedback: float = 0.25
    learning: float = 0.15
    neighborhood: float = 0.25
    initial_noise: float = 0.05
    kernel: str = 'signed'

    def validate(self):
        for key in ['coupling', 'feedback', 'learning', 'neighborhood', 'initial_noise']:
            value = getattr(self, key)
            if not np.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{key} must be in [0,1]')
        if self.kernel not in ['signed', 'absolute', 'none']:
            raise ValueError('Unknown kernel')


class Machine:
    def __init__(self, config=None):
        self.config = config or Config()
        self.config.validate()
        with np.load(HERE / 'native_operators.npz') as data:
            self.native = {k: data[k].copy() for k in data.files}
        k = self.native['kernel']
        if self.config.kernel == 'absolute':
            k = np.abs(k)
        if self.config.kernel == 'none':
            k = np.eye(60)
        self.coupler = k / np.linalg.eigvalsh(k)[-1]
        self.neighbors = np.eye(60) + self.config.neighborhood * (
            (self.native['fidelity'] > 0.5 - 1e-12) & (self.native['fidelity'] < 1 - 1e-12))
        self.rng = np.random.default_rng(self.config.seed)
        f = self.native['features']
        noise = self.rng.normal(size=f.shape) + 1j * self.rng.normal(size=f.shape)
        self.x = unit_rows(f + self.config.initial_noise * unit_rows(noise))
        # Memories start as a shuffled copy, without five-code targets or labels.
        self.memory = self.x[self.rng.permutation(60)].copy()
        self.step_number = 0
        self.history = []
        self.events = []
        self.last_winners = np.arange(60)
        self.last_quantization_error = 0.0
        self.last_memory_change = 0.0
        self.last_delta = 0.0
        self.measure()

    def advance(self, steps=1, measure_every=10):
        c = self.config
        for _ in range(steps):
            old = self.x
            # Collective interaction with inertia. No five-code projection here.
            y = unit_rows((1 - c.coupling) * old + c.coupling * (self.coupler @ old))
            # Euclidean BMU in complex space = ordinary SOM in 70 real dimensions.
            distances = (np.sum(abs(y)**2, axis=1)[:, None]
                         + np.sum(abs(self.memory)**2, axis=1)[None, :]
                         - 2 * (y.conj() @ self.memory.T).real)
            winners = np.argmin(distances, axis=1)
            weights = self.neighbors[winners]
            mass = weights.sum(axis=0)
            target = (weights.T @ y) / np.where(mass[:, None] > 0, mass[:, None], 1)
            active = mass > 0
            updated = self.memory.copy()
            updated[active] += c.learning * (target[active] - self.memory[active])
            self.last_memory_change = float(np.linalg.norm(updated - self.memory) / np.sqrt(60))
            self.memory = updated
            # This is the recursion: the learned output becomes the next input.
            self.x = unit_rows((1 - c.feedback) * y + c.feedback * self.memory[winners])
            self.last_delta = float(np.linalg.norm(self.x - old) / np.sqrt(60))
            self.last_quantization_error = float(np.sqrt(np.maximum(0, distances[np.arange(60), winners])).mean())
            self.last_winners = winners
            self.step_number += 1
            if not np.isfinite(self.x).all() or not np.isfinite(self.memory).all():
                raise RuntimeError('Nonfinite state')
            if self.step_number % measure_every == 0:
                self.measure()
        if self.history[-1]['step'] != self.step_number:
            self.measure()
        return self.history[-1]

    def measure(self):
        eigenvalues = np.linalg.eigvalsh(self.x.conj().T @ self.x).clip(0)
        p = eigenvalues / max(eigenvalues.sum(), 1e-30)
        nonzero = p[p > 1e-15]
        rank = float(np.exp(-np.sum(nonzero * np.log(nonzero))))
        code_weight = float(np.sum(abs(self.x @ self.native['code'])**2) / max(np.sum(abs(self.x)**2), 1e-30))
        gram = abs(self.x.conj() @ self.x.T)**2
        alignment = float((gram.sum() - np.trace(gram)) / (60 * 59))
        probe = self.x[0] @ self.native['code']
        values = dict(step=self.step_number, effective_rank=rank, code_weight=code_weight,
                      alignment=alignment, change=self.last_delta,
                      memory_change=self.last_memory_change,
                      quantization_error=self.last_quantization_error,
                      winners=int(len(np.unique(self.last_winners))),
                      zero_cells=int(np.sum(np.linalg.norm(self.x, axis=1) < 1e-14)),
                      probe_a=float(probe[0].real), probe_b=float(probe[1].real))
        self.history.append(values)
        return values

    def perturb(self, strength=0.1):
        if not 0 <= strength <= 1:
            raise ValueError('Perturbation must be in [0,1]')
        noise = self.rng.normal(size=self.x.shape) + 1j * self.rng.normal(size=self.x.shape)
        self.x = unit_rows(self.x + strength * unit_rows(noise))
        self.events.append(dict(step=self.step_number, action='perturb', strength=strength))
        self.measure()

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        meta = dict(config=asdict(self.config), step=self.step_number, history=self.history,
                    events=self.events, rng=self.rng.bit_generator.state,
                    last_quantization_error=self.last_quantization_error,
                    last_memory_change=self.last_memory_change, last_delta=self.last_delta)
        with path.open('wb') as f:
            np.savez_compressed(f, x=self.x, memory=self.memory, winners=self.last_winners,
                                metadata=np.array(json.dumps(meta)))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            meta = json.loads(str(data['metadata']))
            result = cls(Config(**meta['config']))
            result.x = data['x'].copy()
            result.memory = data['memory'].copy()
            result.last_winners = data['winners'].copy()
        result.step_number = meta['step']
        result.history = meta['history']
        result.events = meta['events']
        result.rng.bit_generator.state = meta['rng']
        for key in ['last_quantization_error', 'last_memory_change', 'last_delta']:
            setattr(result, key, meta[key])
        return result

    def public_state(self):
        coefficients = self.x @ self.native['code']
        return dict(config=asdict(self.config), metrics=self.history[-1], history=self.history[-2000:],
                    code_components=abs(coefficients).tolist(),
                    gram=(abs(self.x.conj() @ self.x.T)**2).tolist(),
                    winners=self.last_winners.tolist(), events=self.events)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--feedback', type=float, default=0.25)
    parser.add_argument('--learning', type=float, default=0.15)
    parser.add_argument('--kernel', choices=['signed', 'absolute', 'none'], default='signed')
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--save', type=Path, default=HERE / 'runs' / 'latest.npz')
    args = parser.parse_args()
    m = Machine.load(args.resume) if args.resume else Machine(Config(
        seed=args.seed, feedback=args.feedback, learning=args.learning, kernel=args.kernel))
    print(json.dumps(m.advance(args.steps), indent=2))
    m.save(args.save)
    args.save.with_suffix('.json').write_text(json.dumps(m.public_state(), indent=2) + '\n')


if __name__ == '__main__':
    main()
