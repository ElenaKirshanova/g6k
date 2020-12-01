#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import absolute_import
from math import e, lgamma, log, pi

from fpylll import BKZ as fplll_bkz, GSO, IntegerMatrix, LLL
from fpylll.tools.bkz_simulator import simulate
from fpylll.util import gaussian_heuristic

from g6k.algorithms.bkz import default_dim4free_fun
from g6k.utils.util import load_lwe_challenge
from six.moves import range

from g6k.utils.ntru_gsa import find_beta

from math import sqrt, pi, e, log, floor, exp, ceil
import scipy.special

def log_gh_svp(d, delta_bkz, svp_dim, n, q):
	"""
	Calculates the log of the Gaussian heuristic of the context in which
	SVP will be ran to try and discover the projected embedded error.

	The volume component of the Gaussian heuristic (in particular the lengths
	of the appropriate Gram--Schmidt vectors) is estimated using the GSA
	[Schnorr03] with the multiplicative factor = delta_bkz ** -2.

	NB, here we use the exact volume of an n dimensional sphere to calculate
	the ``ball_part`` rather than the usual approximation in the Gaussian
	heuristic.

	:param d: the dimension of the embedding lattice = n + m + 1
	:param delta_bkz: the root Hermite factor given by the BKZ reduction
	:param svp_dim: the dimension of the SVP call in context [d-svp_dim:d]
	:param n: the dimension of the LWE secret
	:param q: the modulus of the LWE instance

	"""
	d = float(d)
	svp_dim = float(svp_dim)
	ball_part = ((1./svp_dim)*lgamma((svp_dim/2.)+1))-(.5*log(pi))
	vol_part = ((1./d)*(d-n-1)*log(q))+((svp_dim-d)*log(delta_bkz))
	return ball_part + vol_part

def multinom(n, c):
	assert sum(c) == n, 'bad input to multinom!'
	res = 1
	n_ = n
	for i in range(len(c)):
		res*=scipy.special.binom(n_, c[i])
		n_ = n_ - c[i]
	return res

def BabaiRT(n):
	return n**3

def plain_hybrid_compleixty(paramset, verbose = False):

	q = paramset['q']
	n = paramset['n']
	w = paramset['w']


	best_g = 0
	best_nsamples = 0
	best_beta, best_nsamples, best_prep, best_GSA, best_vol = find_beta(n, q, n)
	best_rt = 2**best_prep

	for g in range(3, int(n/2)):
		beta, nsamples, prep_rt, GSA, vol = find_beta(n-g, q, n) #g determines beta
		w_scaled = (float(w*g) / n) # assume the weight is uniformly distributed over s
		S = multinom(g, [ceil(w_scaled/3.), ceil(w_scaled/3.), g - 2.*ceil(w_scaled/3.)]) # number of CVP batches
		#print('g:', g, beta, w_scaled, S)
		rt_CVP = S*BabaiRT(n)
		rt_log = max(prep_rt, log(rt_CVP, 2)) # not precise
		#print(prep_rt, log(rt_CVP, 2))
		rt = 2**(prep_rt) + rt_CVP
		if rt < best_rt:
			best_g = g
			best_rt = rt
			best_rt_log = rt_log
			best_beta = beta
			best_nsamples = nsamples
			best_GSA = GSA
			best_vol = vol
			if verbose:
				print('rt_log:', best_rt_log, 'beta:', beta,'g:', g)
	return best_beta, best_g, log(best_rt,2), best_nsamples, best_GSA, best_vol



def ntru_plain_hybrid_basis(A, g, q, nsamples):
	"""
		Construct ntru lattice basis
	"""
	n = A.ncols
	ell = n - g
	print(n, ell, nsamples)
	B = IntegerMatrix((nsamples+ell), (ell+nsamples))
	Bg = IntegerMatrix(g, n)


	for i in range(ell):
		B[i,i] = 1
		for j in range(nsamples):
			B[i,j+n] = A[i, j]
	for i in range(nsamples):
		B[i+ell, i+ell] = q

	for i in range(g):
		for j in range(n):
			Bg[i,j] = A[i+ell, j]

	B = LLL.reduction(B)
	return B, Bg

def sim_params(n, alpha):
	A, c, q = load_lwe_challenge(n, alpha)
	stddev = alpha*q
	winning_params = []
	for m in range(60, min(2*n+1, A.nrows+1)):
		B = primal_lattice_basis(A, c, q, m=m)
		M = GSO.Mat(B)
		M.update_gso()
		beta_bound = min(m+1, 110+default_dim4free_fun(110)+1)
		svp_bound = min(m+1, 151)
		rs = [M.get_r(i, i) for i in range(M.B.nrows)]
		for beta in range(40, beta_bound):
			rs, _ = simulate(rs, fplll_bkz.EasyParam(beta, max_loops=1))
			for svp_dim in range(40, svp_bound):
				gh = gaussian_heuristic(rs[M.B.nrows-svp_dim:])
				if svp_dim*(stddev**2) < gh:
					winning_params.append([beta, svp_dim, m+1])
					break
	min_param = find_min_complexity(winning_params)
	return min_param
