from config import Configuration
from libraries.utils import Utils
import numpy as np
import os
from astropy.io import fits
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.centroids import centroid_sources
import pandas as pd
import warnings
from FITS_tools.hcongrid import hcongrid
from astroquery.mast import Catalogs
from astropy.wcs import WCS
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')


class Master:

    @staticmethod
    def pull_master():
        """ This script will generate the master file and photometry file for the data reduction.

        :return - The master frame is returned and the star list is printed
        """

        master, master_header = Master.mk_master(combine_type='median')

        star_list = Master.master_phot(master, master_header)

        # kernel_stars = BigDiff.find_subtraction_stars(star_list)

        return master, star_list  # kernel_stars

    @staticmethod
    def master_phot(master, master_header):
        """ This program will generate the star lists for the master frame and provide a photometry file."""

        if os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt') == 0:

            # create the string useful for query_region
            field = str(Configuration.RA) + " " + str(Configuration.DEC)

            # select the columns we want to import into the data table
            columns = ["toros_field_id", "source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag",
                       "phot_rp_mean_mag", "teff_val", "parallax", "parallax_error", "pmra",
                       "pmra_error", "pmdec", "pmdec_error"]

            # run the query
            Utils.log('Querying MAST for all stars within the toros field: ' + str(Configuration.FIELD), 'info')
            catalog_data = Catalogs.query_region(field,
                                                 radius=Configuration.SEARCH_DIST/1.5,
                                                 catalog="Gaia").to_pandas()
            Utils.log('Query finished. ' + str(len(catalog_data)) + ' stars found.', 'info')

            # add the toros field to the catalog data
            catalog_data['toros_field_id'] = Configuration.FIELD

            # pull out the necessary columns
            star_list = catalog_data[columns]

            # get the header file and convert to x/y pixel positions
            w = WCS(master_header)
            ra = star_list.ra.to_numpy()
            dec = star_list.dec.to_numpy()

            # convert to x, y
            x, y = w.all_world2pix(ra, dec, 0)

            # add the x/y to the star data frame
            star_list['x'] = x
            star_list['y'] = y
            star_list = star_list[(star_list.x >= 10) & (star_list.x < (Configuration.AXS_X - 10)) &
                                  (star_list.y >= 10) & (star_list.y < (Configuration.AXS_X - 10))].copy().reset_index(drop=True)

            star_list['xcen'], star_list['ycen'] = centroid_sources(master,
                                                                    star_list.x.to_numpy(),
                                                                    star_list.y.to_numpy(),
                                                                    box_size=5)
            bd_idxs = np.where(np.isnan(star_list.xcen) | np.isnan(star_list.ycen))
            if len(bd_idxs[0]) > 0:
                for bd_idx in bd_idxs[0]:
                    star_list.loc[bd_idx, 'xcen'] = star_list.loc[bd_idx, 'x']
                    star_list.loc[bd_idx, 'ycen'] = star_list.loc[bd_idx, 'y']

            # centroid the positions
            positions = star_list[['xcen', 'ycen']].copy().reset_index(drop=True)  # positions = (x, y)

            # run aperture photometry
            # set up the star aperture and sky annuli
            aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
            annulus_aperture = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER,
                                               r_out=Configuration.ANNULI_OUTER)
            apers = [aperture, annulus_aperture]

            # run the photometry to get the data table
            phot_table = aperture_photometry(master, apers, method='exact')

            # extract the sky background for each annuli based on either a global or local subtraction
            sky = phot_table['aperture_sum_1'] / annulus_aperture.area

            # subtract the sky background to get the stellar flux and square root of total flux to get the error
            flux = np.array(phot_table['aperture_sum_0'] - sky * aperture.area)

            # calculate the expected photometric error
            flux_er = np.sqrt((phot_table['aperture_sum_0']))

            # convert to magnitude
            mag = 25. - 2.5 * np.log10(flux)
            mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

            # initialize the light curve data frame
            star_list['master_mag'] = mag
            star_list['master_mag_er'] = mag_er
            star_list['master_flux'] = flux
            star_list['master_flux_er'] = flux_er
            star_list['sky'] = sky
            star_list = star_list.reset_index()
            star_list = star_list.rename(columns={'index': 'star_id'})

            star_list.to_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', sep=' ', index=False)

        else:
            star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', delimiter=' ', header=0)

        return star_list

    @staticmethod
    def mk_master(combine_type='median'):
        """ This function will make the master frame that will be used for the differencing

        :return - The master frame is returned and written to the master directory
        """

        file_name = Configuration.FIELD +'_master' + Configuration.FILE_EXTENSION

        if os.path.isfile(Configuration.MASTER_DIRECTORY + file_name) == 0:

            # get the image list
            image_list, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY,
                                                              Configuration.FIELD,
                                                              Configuration.FILE_EXTENSION)
            # determine the number of loops we need to move through for each image
            nfiles = len(image_list)
            if combine_type == 'mean':
                # update the log
                Utils.log("Generating the master frame from multiple files using a mean combination. There are "
                          + str(nfiles) + " images to combine.", "info")

                for kk in range(0, nfiles):
                    # read in the images for the master frame
                    img_tmp = fits.getdata(image_list[kk])

                    # initialize if it's the first file, otherwise....
                    if kk == 0:
                        master = img_tmp
                    else:
                        master = master + img_tmp

                # generate the mean master file
                master = master / nfiles

                # pull the header information from the first file
                master_header = fits.getheader(image_list[0])

                master_header['COMB'] = 'mean'
                master_header['NUM_COMB'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             master, master_header, overwrite=True)

            elif combine_type == 'median':

                # determine the number of loops we need to move through for each image
                nbulk = 20

                # get the integer and remainder for the combination
                full_bulk = nfiles // nbulk
                part_bulk = nfiles % nbulk

                if part_bulk > 0:
                    hold_bulk = full_bulk + 1
                else:
                    hold_bulk = full_bulk

                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_Y, Configuration.AXS_X))

                # update the log
                Utils.log("Generating a master frame from multiple files in bulks of " + str(nbulk) +
                          " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                          str(hold_bulk) + " mini-files to median combine.", "info")

                for kk in range(0, hold_bulk):

                    # loop through the images in sets of nbulk
                    if kk < full_bulk:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(nbulk, Configuration.AXS_Y, Configuration.AXS_X))

                        # generate the max index
                        mx_index = nbulk
                    else:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS_Y, Configuration.AXS_X))

                        # generate the max index
                        mx_index = part_bulk

                    # make the starting index
                    loop_start = kk * nbulk
                    idx_cnt = 0

                    Utils.log("Making mini file " + str(kk) + ".", "info")

                    # now loop through the images
                    for jj in range(loop_start, mx_index + loop_start):
                        # read in the image directly into the block_hold

                        master_tmp, master_head = fits.getdata(image_list[jj], header=True)

                        if (kk == 0) & (jj == 0):
                            block_hold[idx_cnt] = master_tmp
                            master_hold = master_tmp
                            master_header = master_head
                        else:
                            # block_hold[idx_cnt], footprint = aa.register(master_hold, master_tmp)
                            tmp = hcongrid(master_tmp, master_head, master_header)
                            block_hold[idx_cnt] = tmp

                        # increase the iteration
                        idx_cnt += 1

                    # median the data into a single file
                    hold_data[kk] = np.median(block_hold, axis=0)

                # median the mini-images into one large image
                master = np.median(hold_data, axis=0)

                master_header['MAST_COMB'] = 'median'
                master_header['NUM_MAST'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             master, master_header, overwrite=True)
        else:
            master, master_header = fits.getdata(Configuration.MASTER_DIRECTORY + file_name, header=True)

        return master, master_header
