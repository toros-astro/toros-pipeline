""" This script will prepare images for subtraction and then difference the images from the master frame. """
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from config import Configuration
import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils.aperture import aperture_photometry
from astropy.stats import sigma_clipped_stats


class BigDiff:

    @staticmethod
    def difference_images(star_list, kernel_stars_og):
        """ This function will generate the master frame and generates position files, or if a single frame is chosen,
        then only the position file are generated.

        :parameter star_list - A data frame with the aperture photometry from the master image
        :parameter kernel_stars - The list of kernel stars used in the subtraction

        :return - Nothing is returned, however, the images are differenced
        """

        # get the image list to difference
        files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY + Configuration.DATE + '/',
                                    Configuration.FILE_EXTENSION)
        nfiles = len(files)

        # read in the master frame information
        master = fits.getdata(Configuration.MASTER_DIRECTORY + "master.fits")
        master_header = fits.getheader(Configuration.MASTER_DIRECTORY + "master.fits")

        # prepare the oisdifference.c file for differencing
        BigDiff.prep_ois(master, master_header, len(kernel_stars_og))

        # begin with the algorithm to difference the images
        for ii in range(0, nfiles):

            fin_nme = Preprocessing.mk_nme(files[ii], 'Y', 'N', 'N', 'N', 'N')

            if os.path.isfile(Configuration.DIFFERENCED_DIRECTORY + Configuration.DATE + '/' + fin_nme) == 1:
                Utils.log("File " + files[ii] + " found. Skipping...", "info")

            # check to see if the differenced file already exists
            if os.path.isfile(Configuration.DIFFERENCED_DIRECTORY + Configuration.DATE + '/' + fin_nme) == 0:
                Utils.log("Working to difference file " + files[ii] + ".", "info")
                BigDiff.diff_img(star_list,
                                 Configuration.CLEAN_DIRECTORY + Configuration.DATE + '/' + files[ii],
                                 fin_nme)

        Utils.log("Differencing complete for " + Configuration.STAR + ".", "info")

        return

    @staticmethod
    def diff_img(star_list, file, out_name):
        """ This function will check for and determine the reference stars. It will then difference the image.

        :parameter file - The file name to difference.
        :parameter out_name - The final file name.

        :return - Nothing is returned but the image is differenced
        """

        # read in the image
        org_img, org_header = fits.getdata(file, header=True)
        img_sky_mean, img_sky_median, img_sky_std = sigma_clipped_stats(org_img, sigma=3.0)

        # write the new image file
        img_sbkg = org_img - img_sky_median
        fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.fits',
                     img_sbkg, org_header, overwrite=True)

        # get the kernel stars for the subtraction
        kernel_stars = BigDiff.find_subtraction_stars_img(org_img, star_list)

        BigDiff.ois_difference(out_name, org_header)

        return

    @staticmethod
    def ois_difference(out_name, header):
        """ This function will run the c code oisdifference

        :parameter out_name - The file name for the difference file
        :parameter header - The header of the image

        :return Nothing is returned, the image is differenced
        """
        Utils.log('Now starting image subtraction.', 'info')
        Utils.log('The kernel size is: ' + str(Configuration.KRNL * 2 + 1) + 'x' + str(Configuration.KRNL * 2 + 1) +
                  '; the stamp size is: ' + str(Configuration.STMP * 2 + 1) + 'x' + str(Configuration.STMP * 2 + 1) +
                  '; the polynomial is: ' + str(Configuration.ORDR) + ' order; and ' + str(Configuration.NRSTARS) +
                  ' stars were used in the subtraction.', 'info')

        # change to the directory
        os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)

        # run the c code
        os.system('./a.out')

        # update the header file
        dimg, diff_header = fits.getdata('dimg.fits', header=True)
        header['diffed'] = 'Y'

        # update the image with the new file header
        fits.writeto('dimg.fits', dimg[0:Configuration.AXIS_Y, 0:Configuration.AXIS_X], header, overwrite=True)

        # move the differenced file to the difference directory
        os.system('mv dimg.fits ' + Configuration.DIFFERENCED_DATE_DIRECTORY + '/' + out_name)

        # change back to the working directory
        os.chdir(Configuration.WORKING_DIRECTORY)

        # get the photometry from the differenced image
        Utils.log('Image subtraction complete.', 'info')
        return

    @staticmethod
    def prep_ois(master, master_header, n_kstars):
        """ This function will prepare the files necessary for the ois difference.

        :parameter master - The master image for differencing
        :parameter master_header - The header file for the master image

        :return - Nothing is returned but the necessary text files are written,
                    and the code is compiled for differencing
        """
        # compile the oisdifference.c code
        os.system('cp oisdifference.c ' + Configuration.CODE_DIFFERENCE_DIRECTORY)
        os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)
        os.system('gcc oisdifference.c -lcfitsio -lm')
        os.chdir(Configuration.WORKING_DIRECTORY)

        # prepare the master frame
        # clean up the master frame and write to directory
        master_sky_mean, master_sky_median, master_sky_std = sigma_clipped_stats(master, sigma=3.0)

        # subtract the background from the master frame
        master_sbkg = master - master_sky_median

        # write the new master file
        fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.fits', master_sbkg, master_header, overwrite=True)

        # prepare the text files
        # write the parameter file now that we have the stars
        #Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'parms.txt', 'w', "%1d %1d %1d %4d\n" %
        #                 (Configuration.STMP, Configuration.KRNL, Configuration.ORDR, n_kstars))
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.txt', 'w', "ref.fits")
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.txt', 'w', "img.fits")

        return

    @staticmethod
    def find_subtraction_stars(star_list):
        """ This function will find the subtraction stars to use for the differencing, they will be the same stars for
        every frame. This will help in detrending later.

        :parameter star_list - The data frame with the list of stars to use for the subtraction
        """

        Utils.log('Finding stars for kernel from the star list.', 'info')

        # now clip the stars to begin the subtraction
        diff_list = star_list[(star_list['x'] > Configuration.AXS_LIMIT) &
                              (star_list['x'] < Configuration.AXIS_X - Configuration.AXS_LIMIT) &
                              (star_list['y'] > Configuration.AXS_LIMIT) &
                              (star_list['y'] < Configuration.AXIS_Y - Configuration.AXS_LIMIT)].copy().reset_index(drop=True)

        if len(diff_list) > Configuration.NRSTARS:
            diff_list = diff_list.sample(n=Configuration.NRSTARS)
        else:
            Utils.log("There are not enough stars on the frame to use " + str(Configuration.NRSTARS) +
                      " in the subtraction. Using all available stars.", "info")

        # add 1 for indexing in C vs indexing in python
        diff_list['x'] = np.around(diff_list['x'] + 1, decimals=0)
        diff_list['y'] = np.around(diff_list['y'] + 1, decimals=0)

        # export the differencing stars
        diff_list[['x', 'y']].astype(int).to_csv(Configuration.CODE_DIFFERENCE_DIRECTORY +
                                                 'refstars.txt', index=0, header=0, sep=" ")
        diff_list[['star_id', 'x', 'y']].to_csv(Configuration.MASTER_DIRECTORY + "kernel_stars.txt",
                                                index=0, header=0, sep=" ")
        return diff_list

    @staticmethod
    def find_subtraction_stars_img(img, star_list):
        """ This function will find the subtraction stars to use for the differencing, they will be the same stars for
        every frame. This will help in detrending later.

        :parameter star_list - The data frame with the list of stars to use for the subtraction
        """

        Utils.log('Finding stars for kernel from the star list.', 'info')

        positions = np.transpose((star_list['x'], star_list['y']))

        aperture = CircularAperture(positions, r=Configuration.CIRC_APER_SIZE)
        aperture_annulus = CircularAnnulus(positions,
                                           r_in=Configuration.CIRC_ANNULI_INNER,
                                           r_out=Configuration.CIRC_ANNULI_OUTER)
        apers = [aperture, aperture_annulus]

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, apers, method='exact')

        # extract the sky background for each annuli based on either a global or local subtraction
        sky = phot_table['aperture_sum_1'] / aperture_annulus.area

        # subtract the sky background to get the stellar flux and square root of total flux to get the photometric error
        flux = np.array(phot_table['aperture_sum_0'] - (sky * aperture.area)) / Configuration.EXP_TIME

        # calculate the expected photometric error
        star_error = np.sqrt((phot_table['aperture_sum_0'] - (sky * aperture.area)) * Configuration.GAIN)
        sky_error = np.sqrt(aperture.area * sky * Configuration.GAIN)

        # combine sky and signal error in quadrature
        flux_er = np.array(np.sqrt(star_error ** 2 + sky_error ** 2))

        # convert to magnitude
        mag = 25 - 2.5 * np.log10(flux)
        mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

        diff_list = star_list.copy().reset_index(drop=True)
        diff_list['min_dist'] = diff_list.apply(lambda x: np.sort(np.sqrt((x['x'] - diff_list['x']) ** 2 +
                                                                          x['y'] - diff_list['y']) ** 2)[1],
                                                axis=1)
        dist_cut = 2 * Configuration.STMP + 1
        diff_list['dmag'] = np.abs(diff_list['master_mag'].to_numpy() - mag)
        mn, md, sg = sigma_clipped_stats(diff_list.dmag, sigma=2)
        mag_cut = md + sg

        # now clip the stars to begin the subtraction
        diff_list = star_list[(star_list['x'] > Configuration.AXS_LIMIT) &
                              (star_list['x'] < Configuration.AXIS_X - Configuration.AXS_LIMIT) &
                              (star_list['y'] > Configuration.AXS_LIMIT) &
                              (star_list['y'] < Configuration.AXIS_Y - Configuration.AXS_LIMIT) &
                              (diff_list['min_dist'] > dist_cut) &
                              (diff_list['dmag'] < mag_cut)].copy().reset_index(drop=True)

        if len(diff_list) > Configuration.NRSTARS:
            diff_list = diff_list.sample(n=Configuration.NRSTARS)
        else:
            Utils.log("There are not enough stars on the frame to use " + str(Configuration.NRSTARS) +
                      " in the subtraction. Using all available stars.", "info")

        # add 1 for indexing in C vs indexing in python
        diff_list['x'] = np.around(diff_list['x'] + 1, decimals=0)
        diff_list['y'] = np.around(diff_list['y'] + 1, decimals=0)

        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'parms.txt', 'w', "%1d %1d %1d %4d\n" %
                        (Configuration.STMP, Configuration.KRNL, Configuration.ORDR, len(diff_list)))

        # export the differencing stars
        diff_list[['x', 'y']].astype(int).to_csv(Configuration.CODE_DIFFERENCE_DIRECTORY +
                                                 'refstars.txt', index=0, header=0, sep=" ")

        return diff_list
