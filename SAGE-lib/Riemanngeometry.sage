def Christoffelsymbol(G, myvars, i, j, k):
    return sum(
            (1/2) * (1/G)[i,m] * (
                diff(G[m,j], myvars[k])
                +
                diff(G[m,k], myvars[j])
                -
                diff(G[j,j], myvars[m])
            )
            for m in range(len(myvars))
        )

def Riemanntensor(G, myvars, i, j, k, l):
	"""Calculate R_(ijk)^l. That is, the d/dx^l
	coefficient of
		g(R(d/dx^i, d/dx^j) d/dx^k
	"""
	return (
		diff(Christoffelsymbol(G, myvars, i,l,j), myvars[k])
		-
		diff(Christoffelsymbol(G, myvars, i,k,j), myvars[l])
		+
		sum(
			Christoffelsymbol(G, myvars, i, k, m) * Christoffelsymbol(G, myvars, m, l, j)
			+
			Christoffelsymbol(G, myvars, i, l, m) * Christoffelsymbol(G, myvars, m, k, j)
			for m in range(len(myvars))
		)
	)


def Riccitensor(G, myvars, i, j):
    return sum(
            Riemanntensor(G, myvars, m, i, m, j)
            for m in range(len(myvars))
            )

def sectional_curvature(G, myvars, i, j):
	"""Calculate sectional curvature
		K(d/dx^i, d/dx^j)

	Sectional curvature in the (u, v) plane
	is given as
		g(R(u, v)v, u)/(|u|^2|v|^2 - g(u, v)^2)
	`Riemanntensor` calculates R_ijk^l. We get
	the inner product
		g(R(e_i, e_j)e_j, e_i)
	as the sum
		R_ijj^a g_ai
	"""
	return (
		sum(
			Riemanntensor(G, myvars, i, j, j, a) * G[a, i]
			for a in range(len(myvars))
		)/(
			G[i,i]*G[j,j] - G[i, j]**2
		)
	)
