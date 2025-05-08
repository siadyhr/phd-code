var('mu x y t q r')
assume(q>0)
q = function('q')(t, x, y)
r = function('r')(t, x, y)
"""Sol with coordinates (t, x, y),
multiplication
	(s, x1, y1)(t, x2, y2)
	=
	(
		s+t,
		x1 + e^s x2,
		y1 + e^(-s) y2
	)

Cosymplectic structure is
	(mu dt, dx dy),
and critical metrics
	mu^2 dt^2 + e^(-2t) dx^2 + e^(2t) dy^2
"""

import sys
sys.path.append("../SAGE-lib/")
import Riemanngeometry

myvars = [t, x, y]
G = matrix([
	[mu**2, 0, 0],
	[0, (1+r**2)/q, r],
	[0, r, q]
])
G = diagonal_matrix([mu**2, exp(-2*t), exp(2*t)])

def print_Christoffel_symbols():
	for i in range(len(myvars)):
		for j in range(len(myvars)):
			for k in range(len(myvars)):
				if k<j:
					continue
				print(
						"\Gamma^{%s}_{%s, %s} =" % (myvars[i], myvars[j], myvars[k]),
						Riemanngeometry.Christoffelsymbol(G, myvars, i, j, k)
					)

def print_Ricci_tensor():
	for i in range(len(myvars)):
		for j in range(len(myvars)):
			if j<i:
				continue
			print(
				"Ric_{%s, %s} =" % (myvars[i], myvars[j]),
				Riemanngeometry.Riccitensor(G, myvars, i, j)
			)

def print_sectional_curvature():
	for i in range(len(myvars)):
		for j in range(len(myvars)):
			if j <= i:
				continue
			print(
				"K(d/d%s, d/d%s) =" % (myvars[i], myvars[j]),
				Riemanngeometry.sectional_curvature(G, myvars, i, j).simplify_full()
			)

print_Christoffel_symbols()
print()
print_Ricci_tensor()
print()
print_sectional_curvature()
