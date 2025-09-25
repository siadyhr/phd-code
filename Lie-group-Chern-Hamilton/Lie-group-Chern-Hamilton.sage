generators = ", ".join("e%s" % i for i in range(1, 4))
structure_coefficients = {
	('e1', 'e2') : {'e2' : 1},
	('e1', 'e3') : {'e3' : -1},
}

# Create Lie algebra over QQ with generators given labels
# and dictionary of structure coefficients
lie_algebra = LieAlgebra(QQ, generators, structure_coefficients)

# This is a list of the generators
generators = list(lie_algebra.basis())
print(lie_algebra)

# This is a list of the names of the generators
print(lie_algebra.basis().keys())
print(generators)

# Can we solve R = ad_R(X)?
# Equivalently, what is (ad_R)^{-1}(R)?

def ad_e(i):
    """Return a metrix representing
    ad(generators[i])"""
	out = matrix(len(generators), len(generators))
	for j in range(len(generators)):
		for k in range(len(generators)):
#			print(lie_algebra([generators[i], generators[j]]))
            # The (expression)[key] gives the coefficient to the
            # Lie algebra element key
			out[j,k] = (lie_algebra([generators[i], generators[j]])[
				lie_algebra.basis().keys()[k]
			])
	return out)

# Try creating the matrix representing ad(R). Can we calculate
# if its range is consistent with what it should be?
var('R1 R2 R3')
R_coefficients = [R1, R2, R3]
adR = sum(ad_e(i) * R_coefficients[i] for i in range(3))
print(adR)
print(adR.column_space())
#adR.column_space().coordinate_vector([R1, R2, R3])
adR.column_space().coordinate_vector([0, 1, 0])
#print(help(adR.column_space()))
