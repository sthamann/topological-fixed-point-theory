"""Exact all-grade weight counts for the original affine E8 root return.

Fix a root gamma with norm squared two. On the line of full Cartan charges
Z gamma, e=E_gamma,-1, f=e* and h=H_gamma,0-I form an sl(2) triple.
The conserved grade is d=L0-H_gamma,0/2. At charge k gamma the h weight
is 2k-1, the lattice energy is k^2, and the oscillator degree is d+k-k^2.
Each d sector is finite and the positive form makes the triple unitary.
The weight counts therefore determine its finite irreducible multiplicities.
On the central pair h=-1,+1, the spin j+1/2 ladder has amplitude j+1.
"""
from collections import Counter
from math import isqrt

import sympy as sp


def coloured_partitions(limit):
    counts = [1, *([0]*limit)]
    for part in range(1,limit+1):
        for _ in range(8):
            for degree in range(part,limit+1):
                counts[degree] += counts[degree-part]
    return counts


def grade_return(d, counts):
    def oscillator_count(degree):
        return counts[degree] if degree>=0 else 0

    last = (isqrt(4*d+1)-1)//2
    multiplicities = {
        j:oscillator_count(d-j*(j+1))-oscillator_count(d-(j+1)*(j+2))
        for j in range(last+1)
    }
    assert all(n>=0 for n in multiplicities.values())
    for k in range(-last-1,last+3):
        weight = 2*k-1
        from_irreducibles = sum(m for j,m in multiplicities.items() if abs(weight)<=2*j+1)
        assert from_irreducibles==oscillator_count(d+k-k*k)
    spectrum = Counter({(j+1)**2:n for j,n in multiplicities.items() if n})
    assert sum(spectrum.values())==counts[d]
    return dict(spectrum)


def first_excited_return():
    # [E_gamma,-1,H_a,-1] = -<gamma,a> E_gamma,-2 and
    # E_gamma,-2 Omega = gamma(-1)|gamma>. In the oscillator basis
    # the return is the root reflection on the Cartan polarization.
    identity = sp.eye(8)
    gamma = identity[:,0]+identity[:,5]
    reflection = identity-gamma*gamma.T
    assert reflection.T*reflection==identity
    assert reflection*gamma==-gamma
    first,second = identity[:,1],identity[:,2]
    assert reflection*first==first and reflection*second==second and first.dot(second)==0
    return reflection.rank()


if __name__=="__main__":
    counts = coloured_partitions(48)
    spectra = [grade_return(d,counts) for d in range(len(counts))]
    assert spectra[0]=={1:1}
    assert spectra[1]=={1:8}
    assert spectra[2]=={1:43,4:1}
    assert spectra[6][9]==1
    for d,spectrum in enumerate(spectra[:7]):
        print(f"conserved grade {d}: dimension {counts[d]}, F E spectrum {spectrum}")
    print(f"Grade-one neutral to grade-two charged map: root reflection, rank {first_excited_return()}.")
    print("49 exact grade controls passed. Original affine vacuum scope; R113 binding remains open.")
