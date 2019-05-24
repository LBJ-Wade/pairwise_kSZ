import numpy as np
from scipy.spatial import distance as dist
from scipy import sparse
from scipy.spatial import cKDTree as KDTree
import numba, argparse, gzip
from tqdm import tqdm
from utils import RaDec2XYZ
from IPython import embed

# Several routines I've tried to calculate 

#=========================================================================================================
@numba.jit()
def CalculatePairwiseCDist(ra, dec, com_dists, field, sep_min=0, sep_max=300, nbins=20, flip_sign=False):
	assert(len(ra) == len(com_dists))

	nclusts = len(ra)
	delta_sep = (sep_max-sep_min)/nbins
	bins = np.arange(sep_min,sep_max+delta_sep,delta_sep)

	vec_unit = RaDec2XYZ(ra,dec) # Get unit vectors pointing to the clusters
	vec_dist = (vec_unit.T * com_dists).T # Mpc

	costheta = 1 - dist.cdist(vec_unit, vec_unit, 'cosine') # matrix with cos(theta_{ij}) as entries
	com_sep  = dist.cdist(vec_dist, vec_dist, 'euclidean')  # matrix with comoving separation between pairs of clusters
	sep_max  = min(sep_max,bins[-1])

	pairwise    = np.zeros(len(bins)-1)
	denominator = np.zeros(len(bins)-1)
	numerator   = np.zeros(len(bins)-1)

	if flip_sign: # NULL TEST !!!!
		for i in range(1, nclusts):
			for j in range(i):
				if (com_sep[i,j] > sep_min) and (com_sep[i,j] < sep_max):
					# this_bin = np.digitize(com_sep[i,j], bins) - 1 # find the right comoving separation bin
					this_bin = np.int(np.floor((com_sep[i,j]-sep_min)/delta_sep))
					f_ij = field[i] + field[j] # Note the plus instead of minus
					c_ij = (com_dists[i]-com_dists[j]) * (1+costheta[i,j]) / (2*com_sep[i,j])

					numerator[this_bin]   += f_ij * c_ij
					denominator[this_bin] += c_ij**2
	else:
		for i in range(1, nclusts):
			for j in range(i):
				if (com_sep[i,j] > sep_min) and (com_sep[i,j] < sep_max):
					# this_bin = np.digitize(com_sep[i,j], bins) - 1 # find the right comoving separation bin
					this_bin = np.int(np.floor((com_sep[i,j]-sep_min)/delta_sep))
					f_ij = field[i] - field[j]
					c_ij = (com_dists[i]-com_dists[j]) * (1+costheta[i,j]) / (2*com_sep[i,j])
					# print T_ij, c_ij

					numerator[this_bin]   += f_ij * c_ij
					denominator[this_bin] += c_ij**2

	pairwise = numerator / denominator

	return -1 * pairwise


def GetBootstrap(ra, dec, com_dists, field, sep_min=0, sep_max=300, nbins=20, nboot=25, flip_sign=False):
	nclusts = len(ra)
	delta_sep = (sep_max-sep_min)/nbins
	bins = np.arange(sep_min,sep_max+delta_sep,delta_sep)

	pw_bs = np.zeros((nboot,len(bins)-1))

	# Bootstrap estimation of covariance matrix
	print("...calculate covariance through bootstrap...")
	for i in tqdm(range(nboot)):
		idx = np.random.choice(np.arange(len(ra)), size=len(ra))
		ra_tmp = ra[idx]
		dec_tmp = dec[idx]
		com_dists_tmp = com_dists[idx]
		field_tmp = field[idx]
		pw_bs[i] = CalculatePairwiseCDist(ra_tmp, dec_tmp, com_dists_tmp, field_tmp, sep_min=sep_min, sep_max=sep_max, flip_sign=flip_sign)

	print("==>done...")
	return pw_bs

#=========================================================================================================
# This stuff below is not currently used !!!!!!!!!!!
#=========================================================================================================

@numba.jit(parallel=True)
def CalculatePairwiseFlender(ra, dec, com_dists, field, sep_min=0, sep_max=300, nbins=20, flip_sign=False, nboot=20):
	assert(len(ra) == len(com_dists))

	nclusts = len(ra)
	delta_sep = (sep_max-sep_min)/nbins
	bins = np.arange(sep_min,sep_max+delta_sep,delta_sep)

	vec_unit = RaDec2XYZ(ra,dec) # Get unit vectors pointing to the clusters
	vec_dist = (vec_unit.T * com_dists).T # Mpc

	bootstrap_weight = np.zeros((nboot,nclusts))

	pairwise    = np.zeros((nboot,nbins))
	denominator = np.zeros((nboot,nbins))
	numerator   = np.zeros((nboot,nbins))
	npairs      = np.zeros((nboot,nbins))

	for b in range(nboot):
		for i in range(nclusts):
			r = np.random.randint(nclusts)
			bootstrap_weight[b][r] += 1

	if flip_sign: # NULL TEST !!!!
		for i in range(1, nclusts):
			npartners = 0
			for j in range(i):
				costheta = np.dot(vec_unit[i],vec_unit[j])
				com_sep  = np.sqrt(com_dists[i]**2 + com_dists[j]**2 -2*com_dists[i]*com_dists[j]*costheta)#; // comoving 

				if (com_sep >= sep_min) and (com_sep < sep_max):
					this_bin = np.int(np.floor((com_sep-sep_min)/delta_sep))
					f_ij = field[i] + field[j] # Note the plus instead of minus
					c_ij = (com_dists[i]-com_dists[j]) * (1+costheta) / (2*com_sep)

					for b in range(nboot):
						numerator[b][this_bin]   += f_ij * c_ij * bootstrap_weight[b][i] * bootstrap_weight[b][j] 
						denominator[b][this_bin] += c_ij **2 * bootstrap_weight[b][i] * bootstrap_weight[b][j] 
						npairs[b][this_bin]      += bootstrap_weight[b][i] * bootstrap_weight[b][j]
			print("percent complete: %.2f"%(1.*i/(nclusts)*100)) 
	else:
		for i in range(1, nclusts):
			npartners = 0
			for j in range(i):
				costheta = np.dot(vec_unit[i],vec_unit[j])
				com_sep  = np.sqrt(com_dists[i]**2 + com_dists[j]**2 -2*com_dists[i]*com_dists[j]*costheta)#; // comoving 

				if (com_sep >= sep_min) and (com_sep < sep_max):
					this_bin = np.int(np.floor((com_sep-sep_min)/delta_sep))
					f_ij = field[i] - field[j] # Note the plus instead of minus
					c_ij = (com_dists[i]-com_dists[j]) * (1+costheta) / (2*com_sep)

					for b in range(nboot):
						numerator[b][this_bin]   += f_ij * c_ij * bootstrap_weight[b][i] * bootstrap_weight[b][j] 
						denominator[b][this_bin] += c_ij **2 * bootstrap_weight[b][i] * bootstrap_weight[b][j] 
						npairs[b][this_bin]      += bootstrap_weight[b][i] * bootstrap_weight[b][j] 
			print("percent complete: %.2f"%(1.*i/(nclusts)*100))


	pw_bar = np.zeros(nbins)
	mu_nupairs = np.zeros(nbins)

	for k in range(nbins):
		for b in range(nboot):
			pairwise[b][k] = numerator[b][k]/denominator[b][k]
			pw_bar[k] += pairwise[b][k]
			mu_nupairs[k] += npairs[b][k]
		mu_nupairs[k] /= nboot
		pw_bar[k]     /= nboot


	# pw_cov = np.zeros((nbins,nbins))

	# for k in range(nbins):
	# 	for l in range(k):
	# 		for b in range(nboot):
	# 			pw_cov[k+1,l+1] += (pairwise[b][k] - pw_bar[k]) * (pairwise[b][l] - pw_bar[l])
	# 	pw_cov[k+1,l+1] /= nboot
	# 	pw_cov[l+1,k+1] = pw_cov[k+1,l+1]

	return pw_bar, pairwise

def CalculatePairwiseKDTree(ra, dec, com_dists, field, sep_min=0, sep_max=300, nbins=20, flip_sign=False):
	assert(len(ra) == len(com_dists))

	nclusts = len(ra)
	delta_sep = (sep_max-sep_min)/nbins
	bins = np.arange(sep_min,sep_max+delta_sep,delta_sep)

	vec_unit = RaDec2XYZ(ra,dec) # Get unit vectors pointing to the clusters
	vec_dist = (vec_unit.T * com_dists).T # Mpc

	tree = KDTree(vec_dist)
	D    = tree.sparse_distance_matrix(tree, sep_max, p=2.0, output_type='ndarray')
	DU   = D[D['i'] < D['j']]

	pw = PairwiseKDTreeCore(DU, vec_unit, com_dists, field, bins, flip_sign)

	return pw

@numba.jit()
def PairwiseKDTreeCore(DU, vec_unit, com_dists, field, bins, flip_sign):
	pairwise    = np.zeros(len(bins)-1)
	denominator = np.zeros(len(bins)-1)
	numerator   = np.zeros(len(bins)-1)

	if flip_sign:
		# Loop over separation bins
		for ir in range(nbins-1):
			DU_tmp = DU[(DU['v'] > bins[ir]) & (DU['v'] < bins[ir+1])]

			# Loop over pairs in a given r-bin
			for idx in range(len(DU_tmp)):
				costheta = np.dot(vec_unit[idx]['i'],vec_unit[idx]['j'])
				com_sep  = DU_tmp[idx]['v']

				f_ij = field[DU_tmp[idx]['i']] + field[DU_tmp[idx]['j']]
				c_ij = (com_dists[DU_tmp[idx]['i']]-com_dists[DU_tmp[idx]['j']]) * (1+costheta) / (2*com_sep)

				numerator[ir]   += f_ij * c_ij
				denominator[ir] += c_ij**2
	else:
		# Loop over separation bins
		for ir in range(nbins-1):
			DU_tmp = DU[(DU['v'] > bins[ir]) & (DU['v'] < bins[ir+1])]

			# Loop over pairs in a given r-bin
			for idx in range(len(DU_tmp)):
				costheta = np.dot(vec_unit[idx]['i'],vec_unit[idx]['j'])
				com_sep  = DU_tmp[idx]['v']

				f_ij = field[DU_tmp[idx]['i']] - field[DU_tmp[idx]['j']]
				c_ij = (com_dists[DU_tmp[idx]['i']]-com_dists[DU_tmp[idx]['j']]) * (1+costheta) / (2*com_sep)

				numerator[ir]   += f_ij * c_ij
				denominator[ir] += c_ij**2


	pairwise = numerator / denominator

	return -1 * pairwise
