# FOURTH PASS: CK algebra isomorphism and generic III1 null model

**20 November 2026. Negative discrimination test, not a proof of F1.**

The two primitive 2x2 binary adjacency matrices Fib=[[0,1],[1,1]] and Full=[[1,1],[1,1]] both satisfy det(I−Aᵀ)=−1; the exact Smith normal form is diag(1,1). Thus their Cuntz-Krieger C* algebras have K0=K1=0. By the established unital simple Cuntz–Krieger classification, both are **abstractly isomorphic to the same O2**, even though their graph shifts, canonical gauge actions and KMS states are different. The bare algebra loses the information that distinguishes the two clocks. Any original process must retain a distinguished automorphism flow, state and local inclusions; a GNS isomorphism is not automatically a physical identification of flows.

More importantly, **type III1 itself is generic**, not a fingerprint of Fibonacci. For the Cuntz algebra O2 with quasi-free evolution and normalized branch probabilities (p,1−p), the known quasi-free KMS factor theorem gives type III1 when log(p)/log(1−p) is irrational. The standard explicit example p=1/3,1-p=2/3 has type III1 independently of the Fibonacci rule. If p=a/b is a reduced rational in (0,1), and positive m,n satisfied p^m=(1-p)^n, then pairwise gcd(a,b)=gcd(b-a,b)=gcd(a,b-a)=1 forces a=b−a=1, hence p=1/2. **Every other rational p is irrational-log**. Among all reduced p with denominator at most20, 126/127 yield the nonperiodic type III1 class. The remaining p=1/2 is the standard type III_{1/2} case.

Proof and exact 127-point baseline: `cuntz_oq_no_go.py`.

**Consequences for TFPT:** two KMS clocks are not an automatic physical time choice, and algebra type III1 is shared with many unrelated processes. The relevant discriminator must instead be a unique original **local inclusion net**, a dynamical automorphism, a faithful source state, positive-energy translations and gravitational metric response. This branch does not reconstruct those objects. Status remains **Stufe 1 mathematical candidate; F1, P1, P2 physical, gravity open**.

References:
- Cuntz-Krieger K theory and unital classification: https://publications.mfo.de/bitstream/handle/mfo/3289/OWR_2012_20.pdf
- Explicit O2 quasi-free state p=1/3,2/3 yielding III1: https://citeseerx.ist.psu.edu/document?doi=84ff7b6315e755814f216eaf6c23d8252d3ea160&repid=rep1&type=pdf
