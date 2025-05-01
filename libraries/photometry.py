""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.centroids import centroid_sources
import numpy as np
import pandas as pd
import warnings
from astropy.wcs import WCS
from astropy.time import Time
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')


class Photometry:

    @staticmethod
    def single_frame_aperture_photometry(star_list, img_name, fin_name):
        """ This function will find the subtraction stars to use for the differencing, they will be the same stars for
        every frame. This will help in detrending later.

        :parameter star_list - The data frame with the list of stars to use for the subtraction
        :parameter img_name - The file to extract flux from
        """

        # get the image for photometry
        img, header = fits.getdata(img_name, header=True)

        # get the various important header information
        # jd = 2460579.154861 + (int(img_name.split('_')[-2])-168) * (300. / (24. * 60. * 60.))
        time = Time(header['DATE'], format='isot', scale='utc')
        jd = time.jd
        exp_time = header['EXPTIME']

        # get the stellar positions from the master frame
        positions = np.transpose((star_list['x'], star_list['y']))

        aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
        aperture_annulus = CircularAnnulus(positions,
                                           r_in=Configuration.ANNULI_INNER,
                                           r_out=Configuration.ANNULI_OUTER)
        apers = [aperture, aperture_annulus]

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, apers, method='exact')

        # extract the sky background for each annuli based on either a global or local subtraction
        # sky = phot_table['aperture_sum_1'] / aperture_annulus.area
        sky = np.median(img)

        # subtract the sky background to get the stellar flux and square root of total flux to get the photometric error
        img_flux = np.array(phot_table['aperture_sum_0'] - (sky * aperture.area))

        # calculate the expected photometric error
        star_error = np.sqrt(np.abs(phot_table['aperture_sum_0']))

        # combine sky and signal error in quadrature
        img_flux_er = np.array(np.sqrt(star_error ** 2))

        # combine the fluxes
        flux = img_flux.astype(float) + star_list['master_flux'].to_numpy().astype(float)
        flux_er = np.sqrt(img_flux_er.astype(float) ** 2 + star_list['master_flux_er'].to_numpy().astype(float) ** 2)

        # convert to magnitude
        mag = 25. - 2.5 * np.log10(flux)
        mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

        # get the zeropoint
        m_mag = star_list['master_mag'].to_numpy()

        dmag = mag[~np.isnan(mag)] - m_mag[~np.isnan(mag)]

        f_mags = np.arange(np.floor(mag[~np.isnan(mag)].min()), np.floor(mag[~np.isnan(mag)].max())) + 0.5
        n_mags = np.zeros(len(f_mags))
        for mag_idx, m_lw in enumerate(np.arange(np.floor(mag[~np.isnan(mag)].min()),
                                                 np.floor(mag[~np.isnan(mag)].max()))):
            dmag_bin = dmag[np.argwhere((mag[~np.isnan(mag)] > m_lw) & (mag[~np.isnan(mag)] < m_lw + 1))]
            dmag_mn, dmag_md, dmag_sg = sigma_clipped_stats(np.array(dmag_bin, dtype=float), sigma=2.5)
            n_mags[mag_idx] = dmag_md

        # replace nans with -9.999999
        mag = np.where(np.isnan(mag), -9.999999, mag)
        off = np.interp(mag, f_mags, n_mags)

        # generate the final flux file
        flux_file = star_list.copy().reset_index(drop=True)
        flux_file['flux'] = flux
        flux_file['flux_er'] = flux_er
        flux_file['mag'] = mag
        flux_file['mag_er'] = mag_er
        flux_file['sky'] = sky
        flux_file['jd'] = jd
        flux_file['zpt'] = off
        flux_file['cln'] = np.where(np.isnan(mag), -9.999999, mag - off)

        flux_file.to_csv(fin_name, header=True, index=False)
        return

    @staticmethod
    def combine_flux_files(star_list):
        """ This function combines all of the flux files in a given directory into a single data frame.

        :parameter star_list- The master frame star list

        :return - Nothing is being returned, but the raw files are output to disk
        """

        # get the flux files to read in
        files, dates = Utils.get_all_files_per_field(Configuration.DIFFERENCED_DIRECTORY,
                                                     Configuration.FIELD,
                                                     '.flux')
        nfiles = len(files)

        num_rrows = len(star_list)

        # make the holders for the light curves
        jd = np.zeros(nfiles)
        mag = np.zeros((num_rrows, nfiles))
        er = np.zeros((num_rrows, nfiles))
        cln = np.zeros((num_rrows, nfiles))
        zpt = np.zeros((num_rrows, nfiles))
        for idy, file in enumerate(files):

            # read in the data frame with the flux information
            img_flux = pd.read_csv(file, header=0)

            # set the data to the numpy array
            jd[idy] = img_flux.loc[0, 'jd']
            mag[:, idy] = img_flux['mag'].to_numpy()
            er[:, idy] = img_flux['mag_er'].to_numpy()
            cln[:, idy] = img_flux['cln'].to_numpy()
            zpt[:, idy] = img_flux['zpt'].to_numpy()

        # write out the light curve data
        Photometry.write_light_curves(num_rrows, jd, mag, er, cln, zpt)

        return

    @staticmethod
    def write_light_curves(nstars, jd, mag, er, cln, zpt):
        """ This function will write the ETSI columns to a text file for later

        :return - Nothing is returned, but the light curve files are written
        """

        # initialize the light curve data frame
        lc = pd.DataFrame(columns=['jd', 'mag', 'er', 'cln', 'zpt'])

        Utils.log("Starting light curve writing...", "info")

        for idx in range(0, nstars):
            if idx >= 1000:
                star_id = str(idx)
            elif (idx < 1000) & (idx >= 100):
                star_id = '0' + str(idx)
            elif (idx < 100) & (idx >= 10):
                star_id = '00' + str(idx)
            else:
                star_id = '000' + str(idx)

            # add the time, magnitude and error to the data frame
            lc['jd'] = np.around(jd, decimals=6)
            lc['mag'] = np.around(mag[idx, :], decimals=6)
            lc['er'] = np.around(er[idx, :], decimals=6)
            lc['cln'] = np.around(cln[idx, :], decimals=6)
            lc['zpt'] = np.around(zpt[idx, :], decimals=6)
            
            # check or create field specific light curve directory
            Utils.make_field_specific_lightcurve_directory(Configuration.FIELD)

            # write the new file
            lc[['jd', 'cln', 'er', 'mag', 'zpt']].to_csv(Configuration.LIGHTCURVE_DIRECTORY +
                                                         Configuration.FIELD + "/" +
                                                         star_id + ".lc", sep=" ", index=False, na_rep='9.999999')

        return
