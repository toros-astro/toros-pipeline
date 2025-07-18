""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from libraries.photometry import Photometry
import os
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.centroids import centroid_sources
import numpy as np
import pandas as pd
import warnings
from astropy.io import fits
from astropy.wcs import WCS
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')

## Missing os.mkdir for the field of interest.
class Lightcurves:

    @staticmethod
    def generate_flux_files(star_list):
        """ This function will generate light curves for all of the stars in a given star list.

        :parameter star_list - The star list generated byt he master frame

        :return Nothing is returned, but light curve files are generated for all stars
        """

        # get the image list to difference
        files, dates = Utils.get_all_files_per_field(Configuration.DIFFERENCED_DIRECTORY,
                                                     Configuration.FIELD,
                                                     Configuration.FILE_EXTENSION)

        if Configuration.TRANSIENT_LC == 'Y':
            trans = pd.Series(index=star_list.columns.to_list())
            master, master_head = fits.getdata(
                Configuration.MASTER_DIRECTORY + Configuration.FIELD + "_master.fits",
                header=True)
            w = WCS(master_head)
            ra = Configuration.TRANSIENT_RA
            dec = Configuration.TRANSIENT_DEC

            # convert to x, y
            x, y = w.all_world2pix(ra, dec, 0)

            # add the x/y to the star data frame
            trans['x'] = x
            trans['y'] = y

            # centroid on the master frame, but make sure it doesn't return an invalid result
            xcen, ycen = centroid_sources(master, trans.x.item(), trans.y.item(), box_size=5)
            if np.isnan(xcen) | np.isnan(ycen):
                trans['xcen'] = x
                trans['ycen'] = y
            else:
                trans['xcen'] = xcen[0]
                trans['ycen'] = ycen[0]

            # fill in the remainder of the data
            trans['star_id'] = Configuration.TRANSIENT_NAME
            trans['toros_field_id'] = Configuration.FIELD
            trans['source_id'] = '-999'
            trans['ra'] = Configuration.TRANSIENT_RA
            trans['dec'] = Configuration.TRANSIENT_DEC

            positions = np.transpose((trans['x'], trans['y']))

            aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
            aperture_annulus = CircularAnnulus(positions,
                                               r_in=Configuration.ANNULI_INNER,
                                               r_out=Configuration.ANNULI_OUTER)
            apers = [aperture, aperture_annulus]
            phot_table = aperture_photometry(master, apers, method='exact')
            sky = np.median(master)
            flux = np.array(phot_table['aperture_sum_0'] - (sky * aperture.area))
            trans['master_flux'] = np.around(flux, decimals=6)
            star_error = np.sqrt(np.abs(phot_table['aperture_sum_0']))
            flux_er = np.array(np.sqrt(star_error ** 2))
            trans['master_flux_er'] = np.around(flux_er, decimals=6)

            # convert to magnitude
            trans['master_mag'] = np.around(25. - 2.5 * np.log10(flux), decimals=6)
            trans['master_mag_er'] = np.around((np.log(10.) / 2.5) * (flux_er / flux), decimals=6)

            star_list.loc[len(star_list)] = trans

        # begin the algorithm to produce the photometry
        for idx, file in enumerate(files):
            Utils.log(f"File name: {file}", "info")
            fin_nme = file.split('.fits')[0] + '.flux'
            Utils.log(f"fin_nme: {file} \n os.path.isfile(): {os.path.isfile(fin_nme)}", "info")
            if os.path.isfile(fin_nme) == True:
                Utils.log("Flux file " + fin_nme + " found. Skipping...", "info")

            # check to see if the differenced file already exists
            if os.path.isfile(fin_nme) == False:

                Utils.log("Working to extract flux from " + file + ".", "info")
                Photometry.single_frame_aperture_photometry(star_list,
                                                            file,
                                                            fin_nme)

        Utils.log("Differencing complete for " + Configuration.FIELD + ".", "info")

    @staticmethod
    def mk_raw_lightcurves(star_list):
        """ This function will create the individual raw light curve files for each star in the specific star list

        :parameter flux_dir - A data frame with the master frame flux data

        :return nothing is returned, but each light curve is output
        """

        # combine the flux from the flux files, and write the raw light curves
        Photometry.combine_flux_files(star_list)

        return
