

from sage.all import *


def gen_small(s, n):
	"""
	s+1 entries of 1s and s entries of -1s
	"""
	deg = n
	coeff_vector = deg*[0]
	coeff_vector[deg-1] = 1
	coeff_vector[0] = 1
	index_set = set({0,deg-1})
	for i in range(s-2):
	# add 1's
		while True:
			index1 = ZZ.random_element(1,deg-1)
			if not index1 in index_set:
				coeff_vector[index1] = 1
				index_set = index_set.union({index1})
				break
	# add -1's
	for i in range(s):
		while True:
			index2 = ZZ.random_element(1,deg-1)
			if not index2 in index_set:
				coeff_vector[index2] = -1
				index_set = index_set.union({index2})
				break
	return coeff_vector



def print_ntru(q, h, variable_x, filename):
	n = len(list(h))
	f = open(filename, 'w')
	f.write(str(q)+'\n')
	#f.write('[')
	HMat = [0]*n
	for i in range(n):
		hvector = list(h* variable_x**i)
		HMat[i] = hvector
		f.write( str(hvector).replace(',','') +'\n')
		#f.write( str(hvector)+'\n')
	#f.write(']')
	f.close()

	return HMat

def all_rotations(g, variable_x, q):
	n = len(list(g))
	rotations = [0]*(2*n)
	i = 0
	while i < n:
		rotations[2*i] = list(g*variable_x**i)
		rotations[2*i+1] = [-rotations[2*i][j] for j in range(len(rotations[2*i]))]
		for j in range(len(rotations[2*i])):
			if rotations[2*i][j] == q-1:
				rotations[2*i][j] = -1
			if rotations[2*i+1][j] == q-1:
				rotations[2*i+1][j] = -1
		i +=1
	return rotations

def gen_ntru_challenge(n):

	K = CyclotomicField(2*n)

	P = Primes()
	q = next_prime(25*n)


	F = GF(q)
	Fx = PolynomialRing(F, 'x')
	Fx_qou = Fx.quotient(K.polynomial(), 'x')
	variable_x = Fx_qou.gen()

	sparsity = ceil(n/3.)
	f_poly = (gen_small(sparsity, n))
	g_poly = (gen_small(sparsity, n))
	h = Fx_qou(f_poly)/Fx_qou(g_poly)

	rotations = all_rotations(Fx_qou(g_poly),variable_x,q)

	#print('g*h', Fx_qou(g_poly)*h)


	filename = 'ntru_n_'+str(n)+'.txt'
	Hmat = print_ntru(q, h, variable_x, filename)
	Hmat = matrix(ZZ,[hrow for hrow in Hmat])


	qvec = vector(ZZ,g_poly)*Hmat - vector(f_poly)
	assert(len(qvec) == n)
	#print("qvec:", qvec)
	qvec_red = [0]*int(n)
	for i in range(n):
		assert qvec_red[i] % q == 0
		qvec_red[i]  = -qvec[i] / q
	#print("qvec_red:", qvec_red)
	B = matrix(ZZ, 2*n, 2*n)

	for i in range(n):
		B[i,i] = 1
		for j in range(n):
			B[i,n+j] = Hmat[i, j]
		B[i+n, i+n] = q
	#print("B:")
	#print(B)
	f_check = vector(list(g_poly) + list(qvec))*B
	f_check = vector(ZZ, [f_check[i] for i in range(2*n)])
	#print(g_poly, f_poly)
	#print(f_check[:n1])


	"""
	B = B.LLL()
	b0 = B[0][:n1]
	print('b0:', b0, norm(b0))


	for i in range(len(rotations)):
		if vector(b0) == vector(rotations[i]):
			print(i, rotations[i])
			break
	"""
	filename = 'ntru_n_'+str(n)+'_solution.txt'
	f = open(filename, 'w')
	f.write(str(list(f_poly))+'\n')
	f.write(str(list(g_poly)))
	f.close()

	return h, q



if __name__ == '__main__':
	n = 64
	gen_ntru_challenge(64)
