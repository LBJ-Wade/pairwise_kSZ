import numpy as np
import sys, gzip, os, pickle, argparse
from scipy import ndimage
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from tqdm import tqdm

from IPython import embed

arcmin2rad = np.pi / 180. / 60. 
rad2arcmin = 1./arcmin2rad

def GetAP(cut, theta, reso=0.5):
	"""
	Performs aperture photometry.

	Calculates the average value of the cutout within a circle of radius theta minus the average in the outer ring [theta, sqrt(2)*theta]

	Params
	------
	- cut:  the cutout (numpy array)
	- theta: radius (float, in arcmin)
	- reso: pixel resolution (float, in arcmin/pix)
	"""
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
	""" 
	Calculates the average value of the cutout within a circle of radius theta.

	Params
	------
	- cut:  the cutout (numpy array)
	- theta: radius (float, in arcmin)
	- reso: pixel resolution (float, in arcmin/pix)
	"""
	theta_R_pix = theta/reso
	apertures = CircularAperture([(cut.shape[1]/2., cut.shape[0]/2.)], r=theta_R_pix)
	phot_table = aperture_photometry(cut, apertures)
	bkg_mean = phot_table['aperture_sum']
	return bkg_mean/apertures.area()

def GaussSmooth(map, fwhm, reso, order=0):
	"""
	Smooth the map with a Gaussian beam specified by its FWHM (in arcmin).
	- fwhm: 
	- reso: pixel resolution (in arcmin)
	"""
	# reso_  = reso * 180.*60./np.pi # in arcmin
	sigma  = fwhm / np.sqrt(8*np.log(2)) / reso
	# print("\t smoothing map with sigma = %4f" %sigma)
	return ndimage.gaussian_filter(map, sigma=sigma, order=order)

class Pix(object):
	""" 
	Class that contains the pixelization scheme (considers rectangular pixels).
	This is at the core of both the Maps and FFT classes.

	Params
	------
	- nx:  # of pixels in x dir
	- dx: pixel resolution in *arcmin*
	"""
	def __init__(self, nx, dx, ny=None, dy=None):
		if ny is None:
			ny = nx
		if dy is None:
			dy = dx

		# Converting pixel resolution from arcmin -> rad
		dx *= arcmin2rad 
		dy *= arcmin2rad

		self.nx    = nx 
		self.dx    = dx 
		self.ny    = ny 
		self.dy    = dy
		self.npix  = self.nx * self.ny
		self.shape = (self.ny, self.nx)

		# Area (in radians) and sky fraction of the patch
		self.area = (self.nx * self.dx) * (self.ny * self.dy)
		self.fsky = self.area / 4. / np.pi

	def GetLxLy(self, shift=True):
		""" 
		Returns two grids with the (lx, ly) pair associated with each Fourier mode in the map. 
		If shift=True (default), \ell = 0 is centered in the grid
		~ Note: already multiplied by 2\pi 
		"""
		if shift:
			return np.meshgrid( np.fft.fftshift(np.fft.fftfreq(self.nx, self.dx))*2.*np.pi, np.fft.fftshift(np.fft.fftfreq(self.ny, self.dy))*2.*np.pi )
		else:
			return np.meshgrid( np.fft.fftfreq(self.nx, self.dx)*2.*np.pi, np.fft.fftfreq(self.ny, self.dy)*2.*np.pi )

	def GetLxLyAngle(self, shift=True):
		""" 
		Returns a grid with the angle between L and the lx axis. 
		If shift=True (default), \ell = 0 is centered in the grid
		~ Note: already multiplied by 2\pi 
		"""
		lx, ly = self.GetLxLy(shift=shift)
		return 2*np.arctan2(lx, -ly)

	def GetL(self, shift=True):
		""" 
		Returns a grid with the wavenumber l = \sqrt(lx**2 + ly**2) for each Fourier mode in the map. 
		If shift=True (default), \ell = 0 is centered in the grid
		"""
		lx, ly = self.GetLxLy(shift=shift)
		return np.sqrt(lx**2 + ly**2)


	def Compatible(self, other):	
		return ( (self.nx == other.nx) and
			     (self.ny == other.ny) and
				 (self.dx == other.dx) and
				 (self.dy == other.dy) )

def MakeTempMap(cltt, nx, dx, ny=None, dy=None, buff=1):
	"""
	Generate flat-sky temperature map

	Params
	------
	- cltt:  array with the CMB temperature power spectrum
	- nx: number of pixels along the x-axis
	- dx: pixel resolution in *arcmin*
	- if buff > 1, then extract the cutout from a larger box (to avoid periodic boundaries)
	"""
	ls = np.arange(len(cltt))

	Nx = int(nx * buff)
	if ny is None:
		ny = nx
	Ny = int(ny * buff)

	# Assuming all pixels equal !
	pix     = Pix(Nx, dx, ny=Ny, dy=dy) 
	npix    = pix.npix
	shape   = pix.shape
	lx, ly  = pix.GetLxLy(shift=False)
	L       = pix.GetL(shift=False)
	idx     = np.where(L > ls.max())
	xx      = np.interp(L, ls, cltt)
	xx      =  xx / pix.area * (pix.nx*pix.ny)**2
	xx[idx] = 0.
	A       = np.sqrt(xx)
	u       = np.random.normal(size=npix).reshape(shape)+1.0j*np.random.normal(size=npix).reshape(shape)
	kMapX   = A*u
	mapX    = (np.fft.ifft2(kMapX)).real
	mapX    = mapX[int((buff-1)/2)*ny:int((buff+1)/2)*ny,int((buff-1)/2)*nx:int((buff+1)/2)*nx]

	return mapX

def MakeNoiseMap(
    N, pix_size,
    white_noise_level, atmospheric_noise_level,
    one_over_f_noise_level, seed=457349875):
    """
    (From Julien Peloton's BlobFinder)
    Makes a realization of instrument noise, atmosphere and 1/f noise level set at 1 degrees
    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * white_noise_level: float, level of white noise [uk.arcmin]
        * atmospheric_noise_level: float, level of atmospheric noise [uk.arcmin]
        * one_over_f_noise_level: float, level of 1/f noise [uk.arcmin]
    """
    ## make a white noise map
    # state_initial = np.random.RandomState(seed)
    # white_noise = state_initial.normal(
    #     0, 1, (N, N)) * white_noise_level / pix_size
    white_noise = np.random.normal(
        0, 1, (N, N)) * white_noise_level / pix_size

    ## make an atmospheric noise map
    atmospheric_noise = 0.
    if (atmospheric_noise_level != 0):
        ones = np.ones(N)
        inds  = (np.arange(N) + .5 - N / 2.)
        X = np.outer(ones, inds)
        Y = np.transpose(X)

        ## angles relative to 1 degrees
        R = np.sqrt(X**2. + Y**2.) * pix_size /60.

        ## 0.01 is a regularization factor
        mag_k = 2 * np.pi/(R + .01)
        atmospheric_noise = np.fft.fft2(np.random.normal(0, 1, (N, N)))
        atmospheric_noise  = np.fft.ifft2(
            atmospheric_noise * np.fft.fftshift(
                mag_k**(5/3.))) * atmospheric_noise_level / pix_size

    ## make a 1/f map, along a single direction to illustrate striping
    oneoverf_noise = 0.
    if (one_over_f_noise_level != 0):
        ones = np.ones(N)
        inds  = (np.arange(N) + .5 - N/2.)

        ## angles relative to 1 degrees
        X = np.outer(ones, inds) * pix_size /60.
        ## 0.01 is a regularization factor
        kx = 2 * np.pi / (X + .01)
        oneoverf_noise = np.fft.fft2(np.random.normal(0, 1, (N, N)))
        oneoverf_noise = np.fft.ifft2(
            oneoverf_noise * np.fft.fftshift(kx)) * one_over_f_noise_level / pix_size

    ## return the noise map
    noise_map = np.real(white_noise + atmospheric_noise + oneoverf_noise)
    return noise_map

def RaDec2XYZ(ra,dec):
	"""
	From (ra,dec) -> unit vector on the sphere
	"""
	rar  = np.radians(ra)
	decr = np.radians(dec)

	x = np.cos(rar) * np.cos(decr)
	y = np.sin(rar) * np.cos(decr)
	z = np.sin(decr)

	vec = np.array([x,y,z]).T

	return vec
	
def weight_func(z_i, z_j, sigma_z):
	'''
	See Eq. 
	'''
	return np.exp(-0.5*(z_i-z_j)**2/sigma_z**2)