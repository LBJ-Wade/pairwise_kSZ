import numpy as np
import sys, gzip, os, pickle, argparse
from scipy import ndimage
import pandas as pd
from tqdm import tqdm

from IPython import embed

from utils import GetAP, GetAvg, GaussSmooth, MakeTempMap, MakeNoiseMap

"""
example: python make_T_kSZ_table_flender_sims.py -tSZ True -cmb True -noise 5. -fname_output T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_theta2.0arcmin_AP.pkl.gz
"""

RA_SPTSZ_min  = 0
RA_SPTSZ_max  = 165
DEC_SPTSZ_min = -65
DEC_SPTSZ_max = -40

RA_SPTpol_min  = 30
RA_SPTpol_max  = 100
DEC_SPTpol_min = -65
DEC_SPTpol_max = -50

RA_SPTpol_deep_min  = 45
RA_SPTpol_deep_max  = 65
DEC_SPTpol_deep_min = -60
DEC_SPTpol_deep_max = -50

RA_SPTSZ_Soergel_min  = 0
RA_SPTSZ_Soergel_max  = 110
DEC_SPTSZ_Soergel_min = -60
DEC_SPTSZ_Soergel_max = -40


def main(args):
	# Read in the cutouts file
	print("...reading cutouts file %s..."%args.fname_cutouts)
	cutouts = pickle.load(gzip.open(args.fname_cutouts,'rb'))
	print("==> done...")

	# Selecting cutouts in main footprint 
	good_cutouts = np.where(cutouts[args.footprint]==True)[0]

	for key in cutouts.iterkeys():
		cutouts[key] = cutouts[key][good_cutouts]

	# Make sure footprint2 is contained in footprint 1
	if args.footprint2 is not None:
		idx_footprint2 = np.where(cutouts[args.footprint2] == True)[0]
		assert( np.all(cutouts[args.footprint][idx_footprint2] == True) ) 
	else:
		args.footprint2 = args.footprint
		args.noise2     = args.noise


	# Extracting kSZ/tSZ cutouts
	cuts_kSZ  = cutouts['cutouts_kSZ']
	if args.tSZ:
		print("...including tSZ...")
		cuts_tSZ  = cutouts['cutouts_tSZ']

	print("...found %d good cutouts over footprint %s..."%(len(cuts_kSZ),args.footprint))

	np.random.seed(args.seed)

	if args.cmb:
		print("...including primary CMB...")
		l, cl = np.loadtxt('../data/CMB_spectra.dat',unpack=1)
		cltt = np.nan_to_num(cl/l/(l+1)*2*np.pi*(2.725e6)**2)

	if args.noise > 0:
		print("...including instrumental CMB noise at %.2f microK-arcmin..."%args.noise)
		if args.noise2 > 0 and args.footprint2 != args.footprint:
			print("...including instrumental CMB noise for second footprint %s at %.2f microK-arcmin..."%(args.footprint2,args.noise2))

	# Where we dump kSZ estimate
	T_kSZ = []

	# embed()
	# sys.exit()

	print("...starting AP calculation...")
	for idx in tqdm(range(len(cuts_kSZ))):
		if cutouts[args.footprint][idx] == 1:
			cut_tmp = GaussSmooth(cuts_kSZ[idx], args.beam, args.reso) #in \muK
			
			# Add tSZ emission (@150 GHz)?
			if args.tSZ:
				cut_tmp += GaussSmooth(cuts_tSZ[idx], args.beam, args.reso) #in \muK

			# Add primary lensed CMB?
			if args.cmb:
				cut_tmp += GaussSmooth(MakeTempMap(cltt, cut_tmp.shape[0], args.reso), args.beam, args.reso)#in \muK

			# Add instrumental noise?
			if args.noise > 0:
				if cutouts[args.footprint2][idx] == 1:
					cut_tmp += MakeNoiseMap(cut_tmp.shape[0], args.reso, args.noise2, args.atmospheric_noise_level, args.one_over_f_noise_level)
				else:
					cut_tmp += MakeNoiseMap(cut_tmp.shape[0], args.reso, args.noise, args.atmospheric_noise_level, args.one_over_f_noise_level)

			# Performing photometry
			if args.AP_type == 'AP':
				T_kSZ.append( GetAP(cut_tmp, args.theta_arcmin, reso=args.reso) )
			elif args.AP_type == 'AVG':
				T_kSZ.append( GetAvg(cut_tmp, args.theta_arcmin, reso=args.reso) )

			del cut_tmp

	print("==> done...")

	try:
		T_kSZ = np.asarray(T_kSZ)[:,0]
	except:
		pass

	# embed()

	dataset = {"RA":   cutouts['RA'],
			   "DEC":  cutouts['DEC'],
			   "Z":    cutouts['Z'],
			   "M200": cutouts['M200'] ,
			   "VLOS": cutouts['VLOS'] ,
			   "TKSZ": T_kSZ
			   }

	for key, val in vars(args).iteritems():
		dataset[key] = val

	# WE NEED TO CONVERT LAMBA -> MASS -> VIRAL RADIUS TO MAKE AN ADAPTIVE AP !!!

	# Dumping the dataset to a pandas dataframe and to file
	print("...dump zusammen to file %s"%args.fname_output)
	df = pd.DataFrame(data=dataset)
	df.to_pickle(args.fname_output, compression='gzip')
	print("...done...")

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-fname_cutouts', dest='fname_cutouts', action='store', help='Name of the pickle file containing the cutouts', type=str, default='../data/kSZ_tSZ_SPTSZ_footprint_zsm1.pkl.gz')
	parser.add_argument('-AP_type', dest='AP_type', action='store', help='Type of aperture photometry', type=str, default='AP')
	parser.add_argument('-theta_arcmin', dest='theta_arcmin', action='store', help='Radius of the aperture photometry in arcmin ', type=float, default=2.)
	parser.add_argument('-reso', dest='reso', action='store', help='Pixel resolution ', type=float, default=.25)
	parser.add_argument('-beam', dest='beam', action='store', help='FWHM in arcmin ', type=float, default=1.2)
	parser.add_argument('-atmospheric_noise_level', dest='atmospheric_noise_level', action='store', help='atmospheric_noise_level ', type=float, default=0.)
	parser.add_argument('-one_over_f_noise_level', dest='one_over_f_noise_level', action='store', help='one_over_f_noise_level ', type=float, default=0.)
	parser.add_argument('-cmb', dest='cmb', action='store', help='Include primary CMB', type=bool, default=False)
	parser.add_argument('-tSZ', dest='tSZ', action='store', help='Include thermal Sunyaev-Zeldovich cluster emission', type=bool, default=False)
	parser.add_argument('-noise', dest='noise', action='store', help='Include instrumental noise', type=float, default=0.)
	parser.add_argument('-noise2', dest='noise2', action='store', help='Include instrumental noise (if the second footprint is observed with another experiment)', type=float, default=0.)
	parser.add_argument('-seed', dest='seed', action='store', type=int, default=123)
	parser.add_argument('-footprint', dest='footprint', help='Larger footprint', action='store', type=str, default='SPT-SZ')
	parser.add_argument('-footprint2', dest='footprint2',help='Smaller footprint (inside footprint)', action='store', type=str, default=None)
	parser.add_argument('-fname_output', dest='fname_output', action='store', help='Name of the pickle file containing the cutouts', type=str, default='../data/kSZ_tables/results_T_kSZ.pkl.gz')

	args = parser.parse_args()
	main(args)
