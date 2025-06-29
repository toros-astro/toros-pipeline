import pandas as pd
from config import Configuration
from libraries.utils import Utils
import numpy as np
import os
import astropy
import astropy.stats
from astropy.nddata.utils import Cutout2D
from astropy.time import Time
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import twirl
from photutils.detection import DAOStarFinder
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from scipy.interpolate import griddata
import astroalign as aa
import matplotlib
matplotlib.use('TkAgg')


class Preprocessing:

    @staticmethod
    def mk_bias(bias_overwrite, combine_type='mean'):
        """ This function will make the master bias frame using the provided image list and desired method.
        :parameter bias_overwrite - Y/N if you want to force an overwrite for the bias frame
        :parameter combine_type - Right now the combination type is mean, but it can be updated for whatever you
                                method is desired.

        :return - The bias frame is returned and written to the calibration directory
        """

        # set the master bias file name
        file_name = 'bias' + Configuration.FILE_EXTENSION

        # check to see if the master bias frame exists, and if it does save some time by skipping
        if (os.path.isfile(Configuration.CALIBRATION_DIRECTORY + file_name) == 0) | (bias_overwrite == 'Y'):

            if combine_type == 'mean':
                # check if the temporary files exist
                chk_bias_list = Utils.get_file_list(Configuration.BIAS_DIRECTORY, Configuration.FILE_EXTENSION)

                # get the image list
                images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'bias_list.csv', sep=',')
                image_list = images.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.File, axis=1).to_list()

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)

                if len(chk_bias_list) ==0:

                    # update the log
                    Utils.log("Generating a master bias frame from multiple files using a mean combination. There are "
                        + str(nfiles) + " images to combine.", "info")

                     # set the temporary number file count
                     tmp_num = 0

                    # pull the header information from the first file of the set
                    bias_header = fits.getheader(image_list[0])

                    for kk in range(0, nfiles):
                        # read in the bias frame
                        bias_tmp = fits.getdata(image_list[kk]).astype('float')

                        # initialize if it is the first file, otherwise....
                        try:
                            bias_image = bias_image + bias_tmp
                        except:
                            bias_image = bias_tmp

                        # update the log to keep an eye on progress
                        if (kk % 20) == 0:
                            Utils.log(str(kk+1) + " files read in. " + 
                                      str(nfiles - kk -1) + "files remain.", "info")

                # write out the temporary files to save memory when possible
                if ((kk % 100 ==0) & (kk > 0)) | (kk == nfiles - 1 ):
                    # update the log
                    Utils.log("Writing temporary file " + str(tmp_num) + ". " +
                              str(nfiles - kk - 1) + " files remain.", "info")

                    # write out the temporary file
                    fits.writeto(Configuration.BIAS_DIRECTORY + str(tmp_num) + "_tmp_bias.fits",
                                 bias_image, bias_header, overwrite=True)
                    tmp_num = tmp_num + 1
                    del bias_image # save memory
                del bias_tmp # save memory
            else:
                Utils.log("Temporary bias files found. Skipping generating files. "
                          "Delete them if you need to!", "info")
            # get the temporary file list
            tmp_bias_list = Utils.get_file_list(Configuration.BIAS_DIRECTORY, Configuration.FILE_EXTENSION)

            # combine all of the temporary files
            for kk in range (0, len(tmp_bias_list)):
                # read in the bias frame
                bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + tmp_bias_list[kk]).astype('float')

                # initialize if it is the first file, otherwise....
                try:
                    bias_image_fin = bias_image_fin + bias_tmp
                except:
                    bias_image_fin = bias_tmp
            del bias_tmp

            # generate the mean bias file
            bias_image_mean = bias_image _fin / nfiles
            del bias_image_fin

                # update the header with relevant information
                bias_header['BIAS_COMB'] = 'mean'
                bias_header['NUM_BIAS'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.CALIBRATION_DIRECTORY + file_name,
                             bias_image_mean, bias_header, overwrite=True)
            else:
                Utils.log("Specific bias-combination method is not available, try again.", "info")
        else:
            # the bias frame already exists, so go ahead and read that in
            bias_image_mean = fits.getdata(Configuration.CALIBRATION_DIRECTORY + file_name, 0)

        return bias_image_mean

    @staticmethod
    def mk_dark(time_scale, overwrite_dark, combine_type='median'):
        """ This function will make the master dark frame using the provided image list.
        :time_scale - the length of the exposure for the dark frames
        :parameter overwrite_dark - Y/N if you want to force the current file to be overwritten
        :parameter combine_type - Either median or mean depending on how you want to combine the files

        :return - The dark frame is returned and written to the calibration directory
        """

        # make the file name
        file_name = str(time_scale) + 's_dark_' + Configuration.FILE_EXTENSION

        if (os.path.isfile(Configuration.CALIBRATION_DIRECTORY + file_name) == 0) | (overwrite_dark == 'Y'):

            if combine_type == 'mean':
                # check if the temporary files exist
                check_dark_list = Utils.get_file_list(Configuration.DARK_DIRECTORY,
                                                      str(time_scale) + Configuration.FILE_EXTENSION)

                # get the image list
                images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'dark_list.csv', sep=',')
                image_list = images.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.Files, axis=1).to_list()

                nfiles = len(image_list)
                if len(chk_dark_list) == 0:

                    # pull in the bias frame
                    bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY, combine_type='mean')

                    # update the log
                    Utils.log("Generating a master dark frame for time scale " + str(time_scale) +
                    "s from multiple files using a mean combination. There are "
                          + str(nfiles) + " images to combine.", "info")

                    # set the temporary number file count
                    tmp_num = 0

                    # pull the header information from the first file of the set
                    dark_header = fits.getheader(image_list[0])

                    for kk in range(0, nfiles):
                        # read in the dark frame
                        dark_tmp = fits.getdata(image_list[kk]).astype('float')

                        # initialize if it is the first file, otherwise....
                        try:
                            dark_image = (dark_image - bias) + dark_tmp
                        except:
                            dark_image = dark_tmp - bias

                        # update the log to keep an eye on progress
                        if (kk % 20) == 0:
                            Utils.log(str(kk+1) + " files read in. " +
                                      str(nfiles - kk - 1) + " files remain", "info")

                        # write out the temporary files to save memory when possible
                        if ((kk % 100 == 0) & (kk > 0)) | (kk == nfiles -1):
                            #update the log
                            Utils.log("Writing temporary file " + str(tmp_num) + ". " + 
                                      str(nfiles - kk - 1) + " files remain", "info")

                            # write out the temporary file
                            fits.writeto(Configuration.DARK_DIRECTORY + str(tmp_num) + "_"
                                         + str(time_scale) + "s" + "_tmp_dark.fits",
                                         dark_image, dark_header, overwrite=True)
                            del dark_image # save memory
                        del dark_tmp # save memory
                    else:
                        Utils.log("Temporary dark fiels found. Skipping generating files. "
                                  "Delete them if you need to!", "info")
                        # get the temporary file list
                        tmp_dark_list = Utils.get_file_list(Configuration.DARK_DIRECTORY, Configuration.FILE_EXTENSION)

                        # combine all of the temporary files
                        for kk in range(0, len(tmp_dark_list)):
                            # read in the bias frame
                            dark_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + tmp_dark_list[kk]).astype('float')

                            # initialize if it is the first file, otherwise....
                            try:
                                dark_image_fin = dark_image_fin + dark_tmp
                            except:
                                dark_image_fin = dark_tmp
                        del dark_tmp

                        # generate the mean dark file
                        dark_image_mean = dark_image_fin / nfiles
                        del dark_image_fin

                        # update the header with relevant information
                dark_header['DARK_COMB'] = 'mean'
                dark_header['NUM_DARK'] = nfiles
                dark_header['EXPTIME'] = time_scale

                # write the image out to the master directory
                fits.writeto(Configuration.CALIBRATION_DIRECTORY + file_name,
                             dark_image_mean, dark_header, overwrite=True)
            else:
                Utils.log("Specific bias-combination method is not available, try again.", "info")
        else:
            # the bias frame already exists, so go ahead and read that in
            dark_image_mean = fits.getdata(Configuration.CALIBRATION_DIRECTORY + file_name, 0)

        return dark_image_mean

    @staticmethod
    def mk_flat(flat_exp=5, dark_exp=300):
        """ This function will make the master flat frame using the provided image list.
        :parameter flat_exp - The exposure time for the flat frames
        :parameter dark_exp - The exposure time for the dark frames
        :parameter combine_type - Mean/Median depending on how you want to combine the flat frame

        :return - The flat field for the given date is returned and written to the calibration directory
        """

        if os.path.isfile(Configuration.CALIBRATION_DIRECTORY + "flat.fits") == 0:

            # read in the bias file
            bias = Preprocessing.mk_bias(bias_overwrite='N', combine_type='mean')
            dark = Preprocessing.mk_dark(overwrite_dark='N', combine_type='mean')

            # get the image list
            images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'flat_list.csv', sep=',')
            image_list = images.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.Files, axis=1).to_list()

            # determine the number of loops we need to move through for each image
            nfiles = len(image_list)
            nbulk = 30

            # get the integer and remainder for the combination
            full_bulk = nfiles // nbulk
            part_bulk = nfiles % nbulk

            if part_bulk > 0:
                hold_bulk = full_bulk + 1
            else:
                hold_bulk = full_bulk

            # here is the 'holder'
            hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_Y, Configuration.AXS_X))

            #update the log
            Utils.log("Generating a master flat field from multiple files in bulks of "  + str(nbulk) + 
                      " images. There are " + str(nfiles) + " images to combine, which means there should be " + 
                      str(hold_bulk) + " mini-files to combine." , "info")

            for kk in range(0, hold_bulk):

                # loop through the images in sets of nbulk
                if kk < full_bulk:
                    #generate the image holder
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

                Utils.log("Making mini-flat field frame number " + str(kk) + ".", "info")

                # now loop through the images
                for jj in range(loop_start, mx_index + loop_start):

                    # read in the flat file
                    flat_tmp, flat_head = fits.getdata(image_list[jj], header=True)

                    # get the scale factor for the dark frame
                    dark_scale = dark_exp / flat_exp

                    # remove the bias frame from the temporary flat field and scale based on exposure time
                    block_hold[idx_cnt] = ((flat_tmp - bias) -(dark / dark_scale))

                    # increase the iteration
                    idx_cnt += 1

                # median the data into a single file
                hold_data[kk] = np.median(block_hold, axis=0)

            # median the mini-images into one large image
            flat_image = np.median(hold_data, axis=0)
            nflat_image = flat_image / np.median(flat_image[5950:5950 + 3900, 3200:3200 + 900])

            # pull the header information from the first file of the set
            flat_header = fits.getheader(image_list[0])
            flat_header['comb_typ'] = 'median'
            flat_header['median_val'] = np.median(flat_image[5950:5950 + 3900, 3200:3200 + 900])
            flat_header['norm_pix'] = 'median'
            flat_header['num_comb'] = len(image_list)
            flat_header['mean_pix'] = np.mean(nflat_image)
            flat_header['std_pix'] = np.std(nflat_image)
            flat_header['max_pix'] = np.max(nflat_image)
            flat_header['min_pix'] = np.min(nflat_image)
            flat_header['mean'] = np.mean(flat_image)
            flat_header['std'] = np.std(flat_image)
            flat_header['max'] = np.max(flat_image)
            flat_header['min'] = np.min(flat_image)

            # write the image out to the master directory
            fits.writeto(Configuration.CALIBRATION_DIRECTORY + "flat.fits",
                         nflat_image, flat_header, overwrite=True)

        else:
            nflat_image = fits.getdata(Configuration.CALIBRATION_DIRECTORY + "flat.fits")

        return nflat_image

    @staticmethod
    def correct_header(img, header):
        """ This function will plate solve the image, add the time stamp, and exposure time to the header if need be.

        :parameter img - This is the image you want to plate solve
        :parameter header - the header of the image you want to correct

        return img, header - The corrected image will be sent back
        """

        # get the approximate center of the image
        center = SkyCoord(Configuration.RA, Configuration.DEC, unit=['deg', 'deg'])
        pixel = Configuration.PIXEL_SIZE * u.arcsec
        fov = np.max(np.shape(img)) * pixel.to(u.deg)

        # now query the gaia region
        all_stars = twirl.gaia_radecs(center, 1.25 * fov)

        # keep the isolated stars
        all_stars = twirl.geometry.sparsify(all_stars, 0.01)[0:30] #Original: 0.01 degree = 0.6 arcmin 

        # get the stars in the image
        xy = twirl.find_peaks(img[500:,500:])[0:30] + 500

        if len(xy) > 30:
            # now compute the new wcs
            wcs = twirl.compute_wcs(xy[0:30], all_stars[0:30], tolerance=10)

            # add the WCS to the header
            h = wcs.to_header()

            for idx, v in enumerate(h):
                header[v] = (h[idx], h.comments[idx])

            # add additional information such as exposure time and time of exposure
            try:
                header['DATE']
            except:
                header['EXPTIME'] = Configuration.EXP_TIME
                header['DATE'] = Time.now().iso
        else:
            Utils.log("Bad image!", "info")
            header['BADIMAGE'] = 'Y'

        return img, header

    @staticmethod
    def sky_subtract(img, header, sky_write='N'):
        """  The function has been updated to include the photutils background subtraction routine. 

        :parameter img - The image to be cleaned
        :parameter header - The header object to be updated
        :parameter sky_write - Y/N if you want to write the residual background for de-bugging

        :return img_sub, header - The cleaned image and updated header file
        """
        # set up the background clipping functions
        sigma_clip = SigmaClip(sigma=3)
        bkg_estimator = MedianBackground()

        # find the most likely position of the cluster
        daofind = DAOStarFinder(fwhm=3.0, threshold=50)
        sources = daofind(img[3000:, 3000:])

        # find the cluster
        x_cen = (np.sum(sources[sources['flux'] > 0]['xcentroid'] * sources[sources['flux'] > 0]['flux']) /
                 np.sum(sources[sources['flux'] > 0]['flux'])) + 3000
        y_cen = (np.sum(sources[sources['flux'] > 0]['ycentroid'] * sources[sources['flux'] > 0]['flux']) /
                 np.sum(sources[sources['flux'] > 0]['flux'])) + 3000

        # make the masked image
        mask_img = np.zeros((Configuration.AXS_Y, Configuration.AXS_X))

        mask_img[int(y_cen - 1000):int(y_cen + 1000), int(x_cen - 1000):int(x_cen + 1000)] = 1  # TUC-47
        if ((x_cen - 5200) - 200 > 0) & ((x_cen - 5200) > 0):
            mask_img[int((y_cen - 500) - 200):int((y_cen - 500) + 200),
            int((x_cen - 5200) - 200):int((x_cen - 5200) + 200)] = 1 # Other GC

        # do the 2D background estimation
        bkg = Background2D(img, (Configuration.PIX, Configuration.PIX), filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask_img)
        sky = bkg.background

        # subtract the sky gradient and add back the median background
        img_sub = img - sky
        # fin_img = img_sub + np.quantile(sky, 0.25)

        # update the header
        header['sky_medn'] = bkg.background_median
        header['sky_sig'] = bkg.background_rms_median
        header['sky'] = np.quantile(sky, 0.25)
        header['sky_sub'] = 'yes'

        # if desired, write out the sky background to the working directory
        if sky_write == 'Y':
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'sky_background' + Configuration.FILE_EXTENSION, sky, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img' + Configuration.FILE_EXTENSION, img, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img_sub' + Configuration.FILE_EXTENSION, img_sub, header=header, overwrite=True)

        return fin_img, header

    @staticmethod
    def sky_subtract_old(img, header, sky_write='N'):
        """  The function breaks the image into bxs x bxs sections to save memeory, then cleans each section
        the sections are then recombined and smoothed to remove the transitions. The residual image is then subtractedAdd commentMore actions
        from the image and the header is updated appropriately.

        :parameter img - The image to be cleaned
        :parameter header - The header object to be updated
        :parameter sky_write - Y/N if you want to write the residual background for de-bugging

        :return img_sub, header - The cleaned image and updated header file
        """

        # use the sampling space to make the appropriate size vectors
        lop = 2 * Configuration.PIX

        # size holder for later
        sze = int((Configuration.AXS_X / Configuration.PIX) * (Configuration.AXS_Y / Configuration.PIX) +
                  (Configuration.AXS_X / Configuration.PIX) + (Configuration.AXS_Y / Configuration.PIX) + 1)

        # calculate the sky statistics
        sky_mean, sky_median, sky_sig = astropy.stats.sigma_clipped_stats(img, sigma=2.5)

        # create holder arrays for good and bad pixels
        x = np.zeros(shape=sze)
        y = np.zeros(shape=sze)
        v = np.zeros(shape=sze)
        s = np.zeros(shape=sze)
        nd = int(0)

        # begin the sampling of the "local" sky value
        for jj in range(0, Configuration.AXS_X + Configuration.PIX, Configuration.PIX):
            for kk in range(0, Configuration.AXS_Y + Configuration.PIX, Configuration.PIX):
                il = np.amax([jj - lop, 0])
                ih = np.amin([jj + lop, Configuration.AXS_X - 1])
                jl = np.amax([kk - lop, 0])
                jh = np.amin([kk + lop, Configuration.AXS_Y - 1])
                c = img[jl:jh, il:ih]

                # select the median value with clipping
                lsky_mean, lsky, ssky = astropy.stats.sigma_clipped_stats(c, sigma=2.5)

                x[nd] = np.amin([jj, Configuration.AXS_X - 1])  # determine the pixel to input
                y[nd] = np.amin([kk, Configuration.AXS_Y - 1])  # determine the pixel to input
                v[nd] = lsky  # median sky
                s[nd] = ssky  # sigma sky
                nd = nd + 1

        # now we want to remove any possible values which have bad sky values
        rj = np.argwhere(v <= 0)  # stuff to remove
        kp = np.argwhere(v > 0)  # stuff to keep

        if len(rj) > 0:

            # keep only the good points
            xgood = x[kp]
            ygood = y[kp]
            vgood = v[kp]

            for jj in range(0, len(rj[0])):
                # select the bad point
                xbad = x[rj[jj]]
                ybad = y[rj[jj]]

                # use the distance formula to get the closest points
                rd = np.sqrt((xgood - xbad) ** 2. + (ygood - ybad) ** 2.)

                # sort the radii
                pp = sorted(range(len(rd)), key=lambda k: rd[k])

                # use the closest 10 points to get a median
                vnear = vgood[pp[0:9]]
                ave = np.median(vnear)

                # insert the good value into the array
                v[rj[jj]] = ave

        # now we want to remove any possible values which have bad sigmas
        rj = np.argwhere(s >= 2 * sky_sig)
        kp = np.argwhere(s < 2 * sky_sig)

        if len(rj) > 0:
            # keep only the good points
            xgood = np.array(x[kp])
            ygood = np.array(y[kp])
            vgood = np.array(v[kp])

            for jj in range(0, len(rj)):
                # select the bad point
                xbad = x[rj[jj]]
                ybad = y[rj[jj]]

                # use the distance formula to get the closest points
                rd = np.sqrt((xgood - xbad) ** 2. + (ygood - ybad) ** 2.)

                # sort the radii
                pp = sorted(range(len(rd)), key=lambda k: rd[k])

                # use the closest 10 points to get a median
                vnear = vgood[pp[0:9]]
                ave = np.median(vnear)
                if np.isfinite(ave) == 0:
                    ave = np.median(v[np.isfinite(v)])

                # insert the good value into the array
                v[rj[jj]] = ave

        # set up a meshgrid to interpolate to
        xi = np.linspace(0, Configuration.AXS_X - 1, Configuration.AXS_X)
        yi = np.linspace(0, Configuration.AXS_Y - 1, Configuration.AXS_Y)
        xx, yy = np.meshgrid(xi, yi)

        # remove any nan of inf values
        if len(v[~np.isfinite(v)]) > 0:
            v[~np.isfinite(v)] = np.median(v[np.isfinite(v)])

        # now we interpolate to the rest of the image with a cubic interpolation
        res = griddata((x, y), v, (xx, yy), method='cubic')

        # subtract the sky gradient and add back the median background
        img_sub = img - res
        fin_img = img_sub + np.quantile(res, 0.25)

        # update the header
        header['sky_medn'] = sky_median
        header['sky_sig'] = sky_sig
        header['sky'] = np.quantile(res, 0.25)
        header['sky_sub'] = 'yes'

        # if desired, write out the sky background to the working directory
        if sky_write == 'Y':
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'sky_background' + Configuration.FILE_EXTENSION, res, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img' + Configuration.FILE_EXTENSION, img, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img_sub' + Configuration.FILE_EXTENSION, fin_img, header=header, overwrite=True)

        return fin_img, header

    @staticmethod
    def align_img(img, header, ref_path):
        """ This function will align to the ref_path image
        :parameter img - The image to flatten
        :parameter header - The image header file
        :parameter ref_path - The path to the reference iamge

        :return align_img, header - The updated image and header """

        # read in the reference image for alignment
        ref, ref_header = fits.getdata(ref_path, header=True)

        ref_files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE + '/',
                                        Configuration.FILE_EXTENSION)
        ref1, ref1_header = fits.getdata(Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE + '/' + ref_files[0],
                                         header=True)
        # align the image
        # align_img, footprint = aa.register(img, ref)
        transf, (source_list, target_list) = aa.find_transform(ref, ref1)
        align_img = aa.apply_transform(transf, img, ref1)
        return align_img[0], header

    @staticmethod
    def bias_subtract(img, header):
        """ This function will subtract a bias frame
        :parameter img - The image to de-bias / de-dark
        :parameter header - The image header file

        :return bias_sub, header - The updated image and header """

        # read in the bias frame and subtract
        bias = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'bias.fits')

        # subtract the bias from the image
        bias_sub = img - bias

        # update the header
        header['BIAS_SUBT'] = 'Y'

        return bias_sub, header

    @staticmethod
    def dark_subtract(img, header):
        """ This function will subtract a dark frame
        :parameter img - The image to DARK SUBSTRACT
        :parameter header - The image header file

        :return dark_sub, header - The updated image and header """

        # read in the bias frame and subtract
        dark = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'dark.fits')

        # subtract the bias from the image
        dark_sub = img - dark

        # update the header
        header['DARK_SUBT'] = 'Y'

        return dark_sub, header

    @staticmethod
    def clip_image(img, header):
        """ This function will clip the image and remove any overscan regions. This function is written for TOROS
        specifically, and you will need to update it for any given CCD.

        :parameter - image - The image to clip
        :parameter - header - the header of the image

        :return image_clip, header - The clipped image and the new header
        """

        # make the clipped image
        image_clip = np.zeros((Configuration.AXS_Y, Configuration.AXS_X))

        # what are the sizes of the over scan?
        ovs_x_size = 180
        ovs_y_size = 40

        # what are the sizes of each chip (not including the over scan)?
        chip_x_size = 1320
        chip_y_size = 5280

        # what is the full size of the chip (including over scan)
        full_chip_x = chip_x_size + ovs_x_size
        full_chip_y = chip_y_size + ovs_y_size

        # move through x and y
        idx = 0
        for x in range(0, Configuration.AXS_X_RW, full_chip_x):

            idy = 0
            for y in range(0, Configuration.AXS_Y_RW, full_chip_y):

                # put the clipped image into the holder image
                image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y:y + chip_y_size, x:x + chip_x_size]

                # increase the size of the yclip
                idy = idy + chip_y_size

            # increase the size of the xclip
            idx = idx + chip_x_size

        # update the header
        header['OVERSCAN'] = 'removed'
        header['X_CLIP'] = ovs_x_size
        header['Y_CLIP'] = ovs_y_size
        header['NAXIS1'] = image_clip.shape[1]
        header['NAXIS2'] = image_clip.shape[0]

        return image_clip, header

    @staticmethod
    def flat_divide(img, header):
        """ This function will divide a flat field.
        :parameter img - The image to flatten
        :parameter header - The image header file

        :return flat_div, header - The updated image and header """
        # read in the flat frame
        flat = fits.getdata(Configuration.CALIBRATION_DIRECTORY + "flat.fits")

        # subtract the bias from the image
        flat_div = img / flat

        # update the header
        header['FLATTEN'] = 'Y'

        return flat_div, header

    @staticmethod
    def mk_nme(file, difference_image='N', image_clip='N', sky_subtract='N', bias_subtract='N',
               flat_divide='N', dark_subtract="N", plate_solve='N'):
        """ This function will create the appropriate name for the file based on while steps are taken.
        :parameter file - The string with the file name
        :parameter image_clip - Y/N if the image was clipped
        :parameter bias_subtract - Y/N if a bias is subtracted
        :parameter flat_divide - Y/N if a flat field is divided
        :parameter dark_subtract - Y/N if the image was dark subtracted
        :parameter sky_subtract - Y/N if sky subtraction was taken
        :parameter difference_image - Y/N if image subtraction occurred
        :parameter plate_solve - Y/N if the plate solving occurred

        :return file_name - A string with the new file name
        """

        # if everything is N then the file name is the original filename
        file_name = file

        # update the file name with a 'd' if at the differencing step
        if difference_image == 'Y':
            file = file_name.replace("/clean/", "/diff/")
            nme_hld = file.split('.fits')
            file_name = nme_hld[0]  + 'ad' + Configuration.FILE_EXTENSION

        # otherwise...
        if difference_image == 'N':
            # first replace the "raw" directory with the "clean" directory
            file_hld = file_name.split('/')
            file = Configuration.CLEAN_DIRECTORY + file_hld[-2] + "/" + Configuration.FIELD + "/" + file_hld[-1]
            nme_hld = file.split('.fits')

            # update the name to be appropriate for what was done to the file
            # nothing occurs
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'N') \
                    and (dark_subtract == 'N') and (image_clip == 'N'):
                file_name =  nme_hld[0]  + Configuration.FILE_EXTENSION
            # bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_b' + Configuration.FILE_EXTENSION
            # flat
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_f' + Configuration.FILE_EXTENSION
            # sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_s' + Configuration.FILE_EXTENSION
            # dark subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_k' + Configuration.FILE_EXTENSION
            # image clipping only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_c' + Configuration.FILE_EXTENSION
            
            # bias and clip only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bc' + Configuration.FILE_EXTENSION
            # bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bf' + Configuration.FILE_EXTENSION
            # bias and sky_subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bs' + Configuration.FILE_EXTENSION
            # bias and dark only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bk' + Configuration.FILE_EXTENSION
            # bias and dark and sky subtract and image clip only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bkcs' + Configuration.FILE_EXTENSION
            # bias and flat and sky subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bfs' + Configuration.FILE_EXTENSION
            # bias and flat and dark
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bkf' + Configuration.FILE_EXTENSION
            # bias and flat and image clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bcf' + Configuration.FILE_EXTENSION

            # bias and flat and sky and dark
            if (bias_subtract == 'Y') and (flat_divide == 'Y')and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bkfs' + Configuration.FILE_EXTENSION
            # bias and flat and sky and clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y')and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bcfs' + Configuration.FILE_EXTENSION
            # bias and flat and sky and clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y')and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bkcf' + Configuration.FILE_EXTENSION
            # bias and flat and dark and sky and clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bkcfs' + Configuration.FILE_EXTENSION

            # bias and flat and dark and sky and clip and plate_solve
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='Y') and (plate_solve == 'Y'):
                file_name =  nme_hld[0]  + '_bkcfsp' + Configuration.FILE_EXTENSION
                
        return file_name

    @staticmethod
    def remove_overscan(img, header):
        """ This function will remove the over scan region from the images
        :argument img - The np.array with the image
        :argument header - The header object

        :return cut_img, header - Return the new header and cut image"""

        # pull out the wcs in the header
        w = WCS(header)

        # cut the image, but maintain the appropriate header information and positioning - X & Y need to be reversed
        cut = Cutout2D(img, (Configuration.X_CENT, Configuration.Y_CENT),
                       (Configuration.AXIS_Y, Configuration.AXIS_X), wcs=w)

        # create a new image based on the cut data
        cut_img = cut.data

        # update the header
        header['CLIPPED'] = 'Y'

        return cut_img, header
