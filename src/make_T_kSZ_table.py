import numpy as np
import sys, gzip, os, pickle, argparse
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
import pandas as pd
from tqdm import tqdm

def GetAP(cut, theta, reso=0.5):
	theta_R_pix = theta/reso
	positions = [(cut.shape[1]/2., cut.shape[0]/2.)]
	apertures = CircularAperture(positions, r=theta_R_pix)
	annulus_apertures = CircularAnnulus(positions, r_in=theta_R_pix, r_out=theta_R_pix*np.sqrt(2))
	apers = [apertures, annulus_apertures]

	phot_table = aperture_photometry(cut, apers)
	bkg_mean   = phot_table['aperture_sum_1'] / annulus_apertures.area()
	bkg_sum    = bkg_mean * apertures.area()
	final_sum  = phot_table['aperture_sum_0'] - bkg_sum
	phot_table['residual_aperture_sum'] = final_sum

	return phot_table['residual_aperture_sum'][0]/apertures.area()

def GetAvg(cut, theta, reso=0.5):
	theta_R_pix = theta/reso
	apertures = CircularAperture([(cut.shape[1]/2., cut.shape[0]/2.)], r=theta_R_pix)
	phot_table = aperture_photometry(cut, apertures)
	bkg_mean = phot_table['aperture_sum']
	return bkg_mean/apertures.area()

def main(args):
	# Read in the cutouts file
	print("...reading cutouts file %s..."%args.fname_cutouts)
	cuts = pickle.load(gzip.open(args.fname_cutouts,'rb'))
	print("...done...\n")

	reso     = cuts['resol_arcmins']
	keys     = cuts['cutouts'].keys()
	keys_arr = np.asarray(keys)
	cuts     = cuts['cutouts']

	# Indeces of "good" clusters
	print("...calculating number of good cutouts w/ weight > %.2f..."%args.weight_cut)
	good_clusts = np.where( keys_arr[:,4] > args.weight_cut )[0]
	print("...found %d good cutout...\n"%len(good_clusts))

	# Where we dump kSZ estimate
	T_kSZ = []

	print("...starting AP calculation...")
	for idx in tqdm(good_clusts):
		cut_tmp = cuts[keys[idx]][0][80:120,80:120] * 1e6 #in \muK
		T_kSZ.append( GetAP(cut_tmp, args.theta_arcmin, reso=reso) )

	print("...done...\n")

	T_kSZ = np.asarray(T_kSZ)

	dataset = {"RA":     keys_arr[good_clusts,0],
			   "DEC":    keys_arr[good_clusts,1],
			   "ZPH":    keys_arr[good_clusts,2],
			   "LAMBDA": keys_arr[good_clusts,3],
			   "WEIGHT": keys_arr[good_clusts,4],
			   "TKSZ":   T_kSZ
			   }

	# WE NEED TO CONVERT LAMBA -> MASS -> VIRAL RADIUS TO MAKE AN ADAPTIVE AP !!!

	# Dumping the dataset to a pandas dataframe and to file
	print("...dump zusammen to file %s"%args.fname_output)
	df = pd.DataFrame(data=dataset)
	df.to_pickle(args.fname_output, compression='gzip')
	print("...done...")

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-fname_cutouts', dest='fname_cutouts', action='store', help='Name of the pickle file containing the cutouts', type=str, default='/Users/fbianchini/Research/kSZ/data/map_150.pkl.gz_cutouts_150ghz_no_downsampling_y3_v6.4.22_lgt5_redmapper_clusters_full')
	parser.add_argument('-weight_cut', dest='weight_cut', action='store', help='Cutout weight threshold: discard cutout below this value', type=float, default=0.45)
	parser.add_argument('-theta_arcmin', dest='theta_arcmin', action='store', help='Radius of the aperture photometry in arcmin ', type=float, default=2.)
	parser.add_argument('-fname_output', dest='fname_output', action='store', help='Name of the pickle file containing the cutouts', type=str, default='results_T_kSZ.pkl.gz')

	args = parser.parse_args()
	main(args)
