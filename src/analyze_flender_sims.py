import numpy as np
import pandas as pd
from scipy.spatial import distance as dist
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
import numba, argparse, gzip
from tqdm import tqdm
from IPython import embed
import pickle, sys

from pairwise import CalculatePairwiseCDist as CalculatePairwise
from pairwise import GetBootstrap
from utils import weight_func

def main(args):

	# Read in data/settings/...
	print("...reading catalog...")
	cat = pd.read_pickle(args.fname_cat)
	cat = cat[(cat.Z >= args.z_min) & (cat.Z <= args.z_max)]
	cat = cat[(cat.M200 >= args.M_min) & (cat.M200 <= args.M_max)]
	if args.N_clusts > 0:
		if args.N_clusts < len(cat):
			cat = cat.sample(args.N_clusts)
	args.N_clusts = len(cat)

	ra    = cat.RA.values
	dec   = cat.DEC.values
	zs    = cat.Z.values
	v_los = cat.VLOS.values
	T_kSZ = cat.TKSZ.values
	print("==> done...")
	print("...there are %d good clusters between %.2f < z < %.2f and %.2f < log10(M200) < %.2f !"%(len(cat), args.z_min, args.z_max, np.log10(args.M_min), np.log10(args.M_max)))

	# Comoving separation bins 
	delta_sep = (args.r_max-args.r_min)/args.nbins
	bins      = np.arange(args.r_min, args.r_max+delta_sep, delta_sep)
	r         = 0.5*(bins[1:]+bins[:-1])

	# Set up cosmo class
	cosmo = FlatLambdaCDM(H0=args.H0, Om0=args.Om0)

	# embed()
	# sys.exit()

	# Apply photo-z errors
	if args.sigma_photoz > 0.:
		print("...applying photo-z errors with scatter sigma_z = %.3f"%args.sigma_photoz)
		zerr = np.random.normal(loc=0., scale=args.sigma_photoz*(1+zs))
		zs += zerr
		print("==> done...")

	# Compute comoving distance to each cluster
	print("...calculating comoving distances...")
	com_dists = cosmo.comoving_distance(zs).value # Mpc
	print("==> done...")

	if args.apply_debiasing == 1:
		print("...applying z-dependent correction...")
		T_kSZ_new = np.zeros_like(T_kSZ)
		for i in range(len(T_kSZ)):
			T_kSZ_new[i] = T_kSZ[i] - np.sum( weight_func(zs[i], zs, args.sigma_z) * T_kSZ ) / np.sum( weight_func(zs[i], zs, args.sigma_z) ) 

		T_kSZ = T_kSZ_new
		del T_kSZ_new
		print("==> done...")


	# ======= VELOCITY ================================================
	if args.do_velocity == 1:
		# Calculate the pairwise statistic
		print("...calculate Pairwise Estimator [vel]...")
		pw_vel = CalculatePairwise(ra, dec, com_dists, v_los, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=False)
		pw_vel_sims = GetBootstrap(ra, dec, com_dists, v_los, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False)
		pw_vel_cov  = np.cov(pw_vel_sims.T)
		pw_vel_corr = np.corrcoef(pw_vel_sims.T)
		print("==> done...")

		print("...calculate Pairwise Estimator [vel](flip-sign)...")
		pw_flip_vel = CalculatePairwise(ra, dec, com_dists, v_los, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=True)
		pw_flip_vel_sims = GetBootstrap(ra, dec, com_dists, v_los, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False)
		pw_flip_vel_cov  = np.cov(pw_flip_vel_sims.T)
		pw_flip_vel_corr = np.corrcoef(pw_flip_vel_sims.T)
		print("==> done...")

		# Reshuffle the (corrected) temper: this is a NULL TEST!
		if args.reshuffle_temperatures:
			print("...reshuffling temperatures")
			v_los_reshuffle = np.random.permutation(v_los)
			print("==> done...")
			print("...calculate Pairwise Estimator (T-reshuffle)...")
			pw_resh_vlos_vel = CalculatePairwise(ra, dec, com_dists, v_los_reshuffle, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=False)
			pw_resh_vlos_vel_sims = GetBootstrap(ra, dec, com_dists, v_los_reshuffle, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False)
			pw_resh_vlos_vel_cov  = np.cov(pw_resh_vlos_vel_sims.T)
			pw_resh_vlos_vel_corr = np.corrcoef(pw_resh_vlos_vel_sims.T)
			print("==> done...")
		else:
			pw_resh_vlos_vel = None

		# Reshuffle the cluster redshifts: this is a NULL TEST!
		if args.reshuffle_redshift:
			print("...reshuffling redshifts")
			# zs_reshuffle = np.random.permutation(zs)
			com_dists_reshuffle = np.random.permutation(com_dists)
			print("==> done...")
			print("...calculate Pairwise Estimator (z-reshuffle)...")
			pw_resh_zs_vel = CalculatePairwise(ra, dec, com_dists_reshuffle, v_los, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=False)
			pw_resh_zs_vel_sims = GetBootstrap(ra, dec, com_dists_reshuffle, v_los, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False)
			pw_resh_zs_vel_cov  = np.cov(pw_resh_zs_vel_sims.T)
			pw_resh_zs_vel_corr = np.corrcoef(pw_resh_zs_vel_sims.T)

			print("==> done...")
		else:
			pw_resh_zs_vel = None
	else:
		pw_vel                = None
		pw_vel_sims           = None
		pw_vel_cov            = None
		pw_vel_corr           = None
		pw_flip_vel           = None
		pw_flip_vel_sims      = None
		pw_flip_vel_cov       = None
		pw_flip_vel_corr      = None
		pw_resh_vlos_vel      = None
		pw_resh_vlos_vel_sims = None
		pw_resh_vlos_vel_cov  = None
		pw_resh_vlos_vel_corr = None
		pw_resh_zs_vel        = None
		pw_resh_zs_vel_sims   = None
		pw_resh_zs_vel_cov    = None
		pw_resh_zs_vel_corr   = None


	# ======= TEMPERATURE ================================================
	# Calculate the pairwise statistic
	if args.do_temperature:
		print("...calculate Pairwise Estimator [T]...")
		pw_T = CalculatePairwise(ra, dec, com_dists, T_kSZ, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=False)
		pw_T_sims = GetBootstrap(ra, dec, com_dists, T_kSZ, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False)
		pw_T_cov  = np.cov(pw_T_sims.T)
		pw_T_corr = np.corrcoef(pw_T_sims.T)
		print("==> done...")

		print("...calculate Pairwise Estimator [T](flip-sign)...")
		pw_flip_T = CalculatePairwise(ra, dec, com_dists, T_kSZ, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=True)
		pw_flip_T_sims = GetBootstrap(ra, dec, com_dists, T_kSZ, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=True)
		pw_flip_T_cov  = np.cov(pw_flip_T_sims.T)
		pw_flip_T_corr = np.corrcoef(pw_flip_T_sims.T)
		print("==> done...")

		# Reshuffle the (corrected) temper: this is a NULL TEST!
		if args.reshuffle_temperatures:
			print("...reshuffling temperatures")
			T_kSZ_reshuffle = np.random.permutation(T_kSZ)
			print("==> done...")
			print("...calculate Pairwise Estimator (T-reshuffle)...")
			pw_resh_T_kSZ_T = CalculatePairwise(ra, dec, com_dists, T_kSZ_reshuffle, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=False)
			pw_resh_T_kSZ_T_sims = GetBootstrap(ra, dec, com_dists, T_kSZ_reshuffle, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False) 
			pw_resh_T_kSZ_T_cov  = np.cov(pw_resh_T_kSZ_T_sims.T)
			pw_resh_T_kSZ_T_corr = np.corrcoef(pw_resh_T_kSZ_T_sims.T)
			print("==> done...")
		else:
			pw_resh_T_kSZ_T = None

		# Reshuffle the cluster redshifts: this is a NULL TEST!
		if args.reshuffle_redshift:
			print("...reshuffling redshifts")
			# zs_reshuffle = np.random.permutation(zs)
			com_dists_reshuffle = np.random.permutation(com_dists)
			print("==> done...")
			print("...calculate Pairwise Estimator (z-reshuffle)...")
			pw_resh_zs_T = CalculatePairwise(ra, dec, com_dists_reshuffle, T_kSZ, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, flip_sign=False)
			pw_resh_zs_T_sims = GetBootstrap(ra, dec, com_dists_reshuffle, T_kSZ, sep_min=args.r_min, sep_max=args.r_max, nbins=args.nbins, nboot=args.N_bootstraps, flip_sign=False) 
			pw_resh_zs_T_cov  = np.cov(pw_resh_zs_T_sims.T)
			pw_resh_zs_T_corr = np.corrcoef(pw_resh_zs_T_sims.T)
			print("==> done...")
		else:
			pw_resh_zs_T = None
	else:
		pw_T                 = None
		pw_T_sims            = None
		pw_T_cov             = None
		pw_T_corr            = None
		pw_flip_T            = None
		pw_flip_T_sims       = None
		pw_flip_T_cov        = None
		pw_flip_T_corr       = None
		pw_resh_T_kSZ_T      = None
		pw_resh_T_kSZ_T_sims = None
		pw_resh_T_kSZ_T_cov  = None
		pw_resh_T_kSZ_T_corr = None
		pw_resh_zs_T         = None
		pw_resh_zs_T_sims    = None
		pw_resh_zs_T_cov     = None
		pw_resh_zs_T_corr    = None


	embed()


	data={	'pw_T': pw_T,
			'pw_T_cov': pw_T_cov,
			'pw_T_corr': pw_T_corr,
			'pw_flip_T': pw_flip_T,
			'pw_flip_T_cov': pw_flip_T_cov,
			'pw_flip_T_corr': pw_flip_T_corr,
			'pw_resh_zs_T':pw_resh_zs_T,
			'pw_resh_zs_T_cov':pw_resh_zs_T_cov,
			'pw_resh_zs_T_corr':pw_resh_zs_T_corr,
			'pw_resh_T_kSZ_T':pw_resh_T_kSZ_T,
			'pw_resh_T_kSZ_T_cov':pw_resh_T_kSZ_T_cov,
			'pw_resh_T_kSZ_T_corr':pw_resh_T_kSZ_T_corr,
			'pw_vel': pw_vel,
			'pw_vel_cov': pw_vel_cov,
			'pw_vel_corr': pw_vel_corr,
			'pw_flip_vel': pw_flip_vel,
			'pw_flip_vel_cov': pw_flip_vel_cov,
			'pw_flip_vel_corr': pw_flip_vel_corr,
			'pw_resh_vlos_vel'     :pw_resh_vlos_vel,
			'pw_resh_vlos_vel_cov' :pw_resh_vlos_vel_cov ,
			'pw_resh_vlos_vel_corr':pw_resh_vlos_vel_corr,
			'pw_resh_zs_vel'       :pw_resh_zs_vel,
			'pw_resh_zs_vel_cov'   :pw_resh_zs_vel_cov,
			'pw_resh_zs_vel_corr'  :pw_resh_zs_vel_corr,
			'bins':bins,
			'r':r,}
			# 'z_min':args.z_min,
			# 'z_max':args.z_max,
			# 'H0':args.H0,
			# 'Om0':args.Om0,
			# }

	for key, val in vars(args).iteritems():
		data[key] = val


	print("...dump results to file...")
	pickle.dump(data, gzip.open(args.fname_output,'wb'))
	print("==> done...")

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-fname_cat', dest='fname_cat', action='store', help='Name of the pickle file containing the catalog', type=str, default='/Users/fbianchini/Research/kSZ/src/T_kSZ_only_Flender_fwhm1.2arcmin_SPTSZ_footprint_theta2.0arcmin_AP.pkl.gz')
	parser.add_argument('-r_min', dest='r_min', action='store', help='Minimum comoving separation [Mpc]', type=float, default=0.)
	parser.add_argument('-r_max', dest='r_max', action='store', help='Maximum comoving separation [Mpc]', type=float, default=300.)
	parser.add_argument('-nbins', dest='nbins', action='store', help='Number of comoving separation bins',type=int, default=20)
	parser.add_argument('-z_min', dest='z_min', action='store', help='Minimum redshift', type=float, default=0.1)
	parser.add_argument('-z_max', dest='z_max', action='store', help='Maximum redshift', type=float, default=0.8)
	parser.add_argument('-M_min', dest='M_min', action='store', help='Minimum mass [M_sun]', type=float, default=1e14)
	parser.add_argument('-M_max', dest='M_max', action='store', help='Maximum mass [M_sun]', type=float, default=3e14)
	parser.add_argument('-do_velocity', dest='do_velocity', action='store', help='Estimate mean pairwise VELOCITY?', type=int, default=1)
	parser.add_argument('-do_temperature', dest='do_temperature', action='store', help='Estimate mean pairwise kSZ temperature?', type=int, default=1)
	parser.add_argument('-apply_debiasing', dest='apply_debiasing', action='store', help='', type=int, default=1)
	parser.add_argument('-sigma_z', dest='sigma_z', action='store', help='', type=float, default=0.01)
	parser.add_argument('-sigma_photoz', dest='sigma_photoz', action='store', help='', type=float, default=0.0)
	parser.add_argument('-Om0', dest='Om0', action='store', help='Maximum redshift', type=float, default=0.264)
	parser.add_argument('-H0', dest='H0', action='store', help='Hubble constant [km/s/Mpc]', type=float, default=71.)
	parser.add_argument('-reshuffle_redshift', dest='reshuffle_redshift', action='store', help='Reshuffle redshift (for null test)', type=bool, default=True)
	parser.add_argument('-reshuffle_temperatures', dest='reshuffle_temperatures', action='store', help='Reshuffle temperatures (for null test)', type=bool, default=True)
	parser.add_argument('-N_bootstraps', dest='N_bootstraps', action='store', help='Number of bootstraps realization', type=int, default=25)
	parser.add_argument('-N_clusts', dest='N_clusts', action='store', help='Number of clusters', type=int, default=0)
	parser.add_argument('-fname_output', dest='fname_output', action='store', help='Name of the pickle file containing the pairwise estimator', type=str, default='results_pw.pkl.gz')

	args = parser.parse_args()
	main(args)
