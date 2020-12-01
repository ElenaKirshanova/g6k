from math import sqrt, pi, e, log, floor, exp, ceil


def get_delta(beta):
	small = (( 2, 1.02190),  # noqa
			 ( 5, 1.01862),  # noqa
			 (10, 1.01616),
			 (15, 1.01485),
			 (20, 1.01420),
			 (25, 1.01342),
			 (28, 1.01331),
			 (40, 1.01295))

	beta = float(beta)
	if beta <= 2:
		return (1.0219)
	elif beta < 40:
		for i in range(1, len(small)):
			if small[i][0] > beta:
				return (small[i-1][1])
	elif beta== 40:
		return (small[-1][1])
	else:
		return ( beta/(2*pi*e) * (pi*beta)**(1./beta))**(1./(2*(beta-1.)))

def getGSA(q, nqs, nbasis, beta):
	assert beta!=0

	delta = get_delta(beta)
	d = nqs+nbasis
	log_q = log(q)

	zone1 = nqs*[log_q]
	slope = - 2 * log(get_delta(beta))
	zone2_length = int(floor(log_q/-slope))
	zone2 = [log_q + i * slope for i in range(1, zone2_length+1)]
	zone3 = nbasis*[0]
	GSA = zone1+zone2+zone3
	lattice_vol =nqs*log_q
	current_vol = sum([GSA[i] for i in range(d)])

	#correct the volume since now current_vol > lattice_vol
	ind = 0
	while current_vol>lattice_vol:
		current_vol -= GSA[ind]
		current_vol += GSA[ind+d]
		ind += 1


	assert ind<=zone2_length
	GSA = GSA[ind:ind+d]

	i_index = max(0, nqs-ind)
	j_index = min(zone2_length, d - i_index)


	#error we make by overshooting the while loop
	err = lattice_vol - current_vol
	for i in range(i_index, i_index+j_index):
		GSA[i] += err / j_index

	current_vol = sum(GSA)
	assert abs(current_vol/lattice_vol - 1) < 1e-6

	return GSA, i_index, i_index+j_index, current_vol


def sievig0292(dim):
	return dim*log(sqrt(3./2))/log(2.)

def sievig0349(dim):
	return dim*0.349 + 16


def find_beta(n, q, nsamples, svp_alg = sievig0349):
	rt_min = 10000
	nsamples_opt = nsamples
	beta_opt = 10000
	GSA_opt = []
	#"""
	for b in range(40, n+nsamples, 2):
		if svp_alg(b) > rt_min:
			break
		#for m in range(max(b-n, int(nsamples/2.)), nsamples, 5):
		m = nsamples
		GSA, i, j, vol  = getGSA(q, m, n, b)
		if exp(GSA[m+n-b]) > sqrt(2.0/3.0 * b)+1 and svp_alg(b) < rt_min:
			rt_min = svp_alg(b)
			nsamples_opt = m
			beta_opt = b
			GSA_opt = GSA
			vol_opt = vol
	"""
	for nsample in range( int(nsamples/3.), nsamples, 5):
		dim = n+nsample
		for beta in range(20,dim,2):
			if svp_alg(beta)>=rt_min:
				break
			GSA, i, j  = getGSA(q, nsample, n, beta)
			if exp(GSA[dim-beta]) > sqrt(2.0/3.0 * beta) and svp_alg(beta)<rt_min: #2/3 is the error sparsity (assume the error is uniform from {+-1, 0}
				beta_opt = beta
				nsamples_opt = nsample
				rt_min = svp_alg(beta)
				GSA_opt = GSA
				#print('find_beta', beta_opt, nsamples_opt, rt_min)
				break
	"""
	return beta_opt, nsamples_opt, rt_min, GSA_opt, vol
