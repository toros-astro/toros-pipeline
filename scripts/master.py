from config import Configuration
from libraries.utils import Utils
import numpy as np
import os
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astropy.io import fits
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils.aperture import aperture_photometry
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')
from scripts.difference import BigDiff

class Master:

    @staticmethod
    def pull_master():
        """ This script will generate the master file and photometry file for the data reduction.

        :return - The master frame is returned and the star list is printed
        """

        master = Master.mk_master(combine_type='median')

        star_list = Master.master_phot(master)

        kernel_stars = BigDiff.find_subtraction_stars(star_list)

        return master, star_list, kernel_stars

    @staticmethod
    def master_phot(master):
        """ This program will generate the star lists for the master frame and provide a photometry file."""

        if os.path.isfile(Configuration.MASTER_DIRECTORY + 'star_list.txt') == 0:

            # get the image statistics
            master_mn, master_mdn, master_std = sigma_clipped_stats(master, sigma=3)
            daofind = DAOStarFinder(fwhm=7., threshold=10. * master_std)
            sources = daofind(master - master_mdn)

            for col in sources.colnames:
                sources[col].info.format = '%.6g'

            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

            aperture = CircularAperture(positions, r=Configuration.CIRC_APER_SIZE)
            aperture_annulus = CircularAnnulus(positions,
                                               r_in=Configuration.CIRC_ANNULI_INNER,
                                               r_out=Configuration.CIRC_ANNULI_OUTER)
            apers = [aperture, aperture_annulus]

            # run the photometry to get the data table
            phot_table = aperture_photometry(master, apers, method='exact')

            # extract the sky background for each annuli based on either a global or local subtraction
            sky = phot_table['aperture_sum_1'] / aperture_annulus.area

            # subtract the sky background to get the stellar flux and square root of total flux to get the error
            flux = np.array(phot_table['aperture_sum_0'] - sky * aperture.area) / Configuration.EXP_TIME

            # calculate the expected photometric error
            star_error = np.sqrt((phot_table['aperture_sum_0'] - sky * aperture.area) * Configuration.GAIN)
            sky_error = np.sqrt(aperture.area * sky * Configuration.GAIN)

            # combine sky and signal error in quadrature
            flux_er = np.array(np.sqrt(star_error ** 2 + sky_error ** 2))

            # convert to magnitude
            mag = 25. - 2.5 * np.log10(flux)
            mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

            # initialize the light curve data frame
            star_list = pd.DataFrame(columns={'x', 'y', 'master_flux', 'master_flux_er',
                                              'master_mag', 'master_mag_er', 'sky'})
            star_list['x'] = sources['xcentroid']
            star_list['y'] = sources['ycentroid']
            star_list['master_mag'] = mag
            star_list['master_mag_er'] = mag_er
            star_list['master_flux'] = flux
            star_list['master_flux_er'] = flux_er
            star_list['sky'] = sky
            star_list = star_list.reset_index()
            star_list = star_list.rename(columns={'index': 'star_id'})
            star_list = star_list[star_list['master_mag_er'] < 0.5].copy().reset_index(drop=True)

            star_list.to_csv(Configuration.MASTER_DIRECTORY + 'star_list.txt', sep=' ', index=False)

        else:
            star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + 'star_list.txt', delimiter=' ', header=0)

        return star_list

    @staticmethod
    def mk_master(combine_type='median'):
        """ This function will make the master frame that will be used for the differencing

        :return - The master frame is returned and written to the master directory
        """

        file_name = 'master' + Configuration.FILE_EXTENSION

        if os.path.isfile(Configuration.MASTER_DIRECTORY + file_name) == 0:

            # get the image list
            image_list = Utils.get_file_list(Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE,
                                             Configuration.FILE_EXTENSION)

            # determine the number of loops we need to move through for each image
            nfiles = len(image_list)

            if combine_type == 'mean':
                # update the log
                Utils.log("Generating the master frame from multiple files using a mean combination. There are "
                          + str(nfiles) + " images to combine.", "info")

                for kk in range(0, nfiles):
                    # read in the images for the master frame
                    img_tmp = fits.getdata(Configuration.CLEAN_DIRECTORY +
                                           Configuration.REF_DATE + '/' + image_list[kk])

                    # initialize if its the first file, otherwise....
                    if kk == 0:
                        master = img_tmp
                    else:
                        master = master + img_tmp

                # generate the mean master file
                master = master / nfiles

                # pull the header information from the first file
                master_header = fits.getheader(Configuration.CLEAN_DIRECTORY +
                                               Configuration.REF_DATE + '/' + image_list[0])

                master_header['COMB'] = 'mean'
                master_header['NUM_COMB'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             master, master_header, overwrite=True)

            elif combine_type == 'median':

                # determine the number of loops we need to move through for each image
                nbulk = 100

                # get the integer and remainder for the combination
                full_bulk = nfiles // nbulk
                part_bulk = nfiles % nbulk

                if part_bulk > 0:
                    hold_bulk = full_bulk + 1
                else:
                    hold_bulk = full_bulk

                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXIS_Y, Configuration.AXIS_X))

                # update the log
                Utils.log("Generating a master frame from multiple files in bulks of " + str(nbulk) +
                          " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                          str(hold_bulk) + " mini-files to median combine.", "info")

                for kk in range(0, hold_bulk):

                    # loop through the images in sets of nbulk
                    if kk < full_bulk:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(nbulk, Configuration.AXIS_Y, Configuration.AXIS_X))

                        # generate the max index
                        mx_index = nbulk
                    else:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(part_bulk, Configuration.AXIS_Y, Configuration.AXIS_X))

                        # generate the max index
                        mx_index = part_bulk

                    # make the starting index
                    loop_start = kk * nbulk
                    idx_cnt = 0

                    Utils.log("Making mini file " + str(kk) + ".", "info")

                    # now loop through the images
                    for jj in range(loop_start, mx_index + loop_start):
                        # read in the image directly into the block_hold

                        master_tmp = fits.getdata(
                            Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE + '/' + image_list[kk])

                        block_hold[idx_cnt] = master_tmp

                        # increase the iteration
                        idx_cnt += 1

                    # median the data into a single file
                    hold_data[kk] = np.median(block_hold, axis=0)

                # median the mini-images into one large image
                master = np.median(hold_data, axis=0)

                # pull the header information from the first file of the set
                master_header = fits.getheader(Configuration.CLEAN_DIRECTORY +
                                               Configuration.REF_DATE + '/' + image_list[0])

                master_header['MAST_COMB'] = 'median'
                master_header['NUM_MAST'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             master, master_header, overwrite=True)
        else:
            master = fits.getdata(Configuration.MASTER_DIRECTORY + file_name, 0)

        return master
