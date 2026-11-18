# Six origin routes and physical reconstruction

OriginLab evaluates matrix sewing and adjoints, exact Petz recovery, finite
moment closure, regional algebras and the joint cone, finite chiral index and
abelian anomaly/inflow, and independent scattering constraints. It tests all
57 nontrivial subsets of these routes for each of 32 generated control
variations. A conjunction of necessary checks does not establish a physical
model. For Markov plus moment closure, it also replays the actual word
columns, generator products and marked forms on the recovered state.

The reconstruction lane finds minimal-support Pauli separators, computes
state-dependent Gram rank and reconstructs the null quotient with its metric.
It tests operator descent, adjoints and one and two step transfer forms.
Exact all-word comparison either returns a distinguishing word or closes
the joint reachable subspace. For cyclic star representations it also
returns and verifies a metric-unitary intertwiner. This is a claim about the
declared common interface. It says nothing about unobserved sectors.

The generated complex lane exhausts finite labelled complexes at declared
bounds. The rational transfer lane exhausts the configured coefficient grid.
Use a manifest to add exact inputs beyond the built-in families.

```sh
python -B research/originlab/autonomy.py --out /path/to/campaign \
  --cases research/originlab/reconstruction-cases.json --workers 6 \
  --max-resolution 5 --lean-toolchain leanprover/lean4:v4.29.1 --watch-seconds 90
```

Custom jobs support `custom_complex` with simplices,
`compare_reconstructions` with marked coherent documents, and
`route_candidate` with the six-route input format in routes.py.
The supplied examples include a three-sphere triangulation and two
different coordinate descriptions of the same finite quantum process.

Original R113 bundles may supply `origin_route_inputs` with identity
`domain_id:native_inputs_sha256`. The declared moments must equal the
original Gram and word labels. A separate `physical_reconstruction` record
accepts `identity_index`, `forms`, optional `step`, `twice`, `adjoints`,
`density_model`, and `regional`. All original forms use the source's own
word basis and Gaussian rational pairs. Regional data use four translation
forms ordered P0, Px, Py, Pz, named algebra bases and transport forms.
The finite regional test verifies closure, whole algebra transport and the
joint cone. Linking those transports to the continuous translation group
remains an explicit proof obligation.

Unbounded scenario coverage is not a decidable acceptance criterion.
An exact verdict covers the stated class and input. Open work includes
the original weighted P1 packet, physical continuum reconstruction,
global anomalies, a physical mirror gap, the complete scattering problem,
and equivalence of full interacting theories.

References: [Walck and Lyons](https://arxiv.org/abs/0808.0859),
[Hayden, Jozsa, Petz and Winter](https://arxiv.org/abs/quant-ph/0304007),
[Schützenberger](https://www-igm.univ-mlv.fr/~berstel/SiteSchutzenberger/Travaux/A/1961-4DefinitionFamilyAutInfCtl.pdf).
