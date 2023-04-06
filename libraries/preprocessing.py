from config import Configuration
from libraries.utils import Utils
import numpy as np
import os
from FITS_tools.hcongrid import hcongrid
import astropy
import astropy.stats
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits
import scipy
import scipy.ndimage
from scipy.interpolate import Rbf
from scipy.interpolate import griddata
from scipy.interpolate import interpn
import astroalign as aa
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class Preprocessing:

    @staticmethod
    def mk_bias(image_directory, dark="N", combine_type='median'):
        """ This function will make the master bias frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination
        :parameter dark - if the file being generated is a dark, update the name
        :parameter combine_type - Either median or mean depending on how you want to combine the files

        :return - The bias frame is returned and written to the master directory
        """

        if dark == 'Y':
            file_name = 'dark' + Configuration.FILE_EXTENSION
        else:
            file_name = 'bias' + Configuration.FILE_EXTENSION

        if os.path.isfile(Configuration.MASTER_DIRECTORY + file_name) == 0:

            if combine_type == 'median':
                # get the image list
                image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)
                nbulk = 20

                # get the integer and remainder for the combination
                full_bulk = nfiles // nbulk
                part_bulk = nfiles % nbulk
                if part_bulk > 0:
                    hold_bulk = full_bulk + 1
                else:
                    hold_bulk = full_bulk

                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_X, Configuration.AXS_Y))

                # update the log
                Utils.log("Generating a master bias/dark frame from multiple files in bulks of " + str(nbulk) +
                          " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                          str(hold_bulk) + " mini-files to combine.", "info")

                for kk in range(0, hold_bulk):

                    # loop through the images in sets of nbulk
                    if kk < full_bulk:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(nbulk, Configuration.AXS_X, Configuration.AXS_Y))

                        # generate the max index
                        mx_index = nbulk
                    else:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS_X, Configuration.AXS_Y))

                        # generate the max index
                        mx_index = part_bulk

                    # make the starting index
                    loop_start = kk * nbulk
                    idx_cnt = 0

                    Utils.log("Making mini-bias/dark file " + str(kk) + ".", "info")

                    # now loop through the images
                    for jj in range(loop_start, mx_index + loop_start):
                        # read in the image directly into the block_hold
                        if dark == 'Y':
                            bias_tmp = fits.getdata(Configuration.DARKS_DIRECTORY + image_list[jj])
                        else:
                            bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + image_list[jj])

                        if np.ndim(bias_tmp) > 2:
                            bias_tmp = bias_tmp[0]

                        block_hold[idx_cnt] = bias_tmp

                        # increase the iteration
                        idx_cnt += 1

                    # median the data into a single file
                    hold_data[kk] = np.median(block_hold, axis=0)

                # median the mini-images into one large image
                bias_image = np.median(hold_data, axis=0)

                # pull the header information from the first file of the set
                if dark == 'Y':
                    bias_header = fits.getheader(Configuration.DARKS_DIRECTORY + image_list[0])
                else:
                    bias_header = fits.getheader(Configuration.BIAS_DIRECTORY + image_list[0])
                bias_header['BIAS_COMB'] = 'median'
                bias_header['NUM_BIAS'] = nfiles
                bias_header["BIASSUB"] = 'Y'

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             bias_image, bias_header, overwrite=True)
            else:
                # get the image list
                image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)

                # update the log
                Utils.log("Generating a master bias/dark frame from multiple files using a mean combination. There are "
                          + str(nfiles) + " images to combine.", "info")

                for kk in range(0, nfiles):
                    # read in the bias frame
                    bias_tmp = fits.getdata(Configuration.DARKS_DIRECTORY + image_list[kk])

                    # initialize if its the first file, otherwise....
                    if kk == 0:
                        bias_image = bias_tmp
                    else:
                        bias_image = bias_image + bias_tmp

                # generate the mean bias file
                bias_image = bias_image / nfiles

                # pull the header information from the first file of the set
                if dark == 'Y':
                    bias_header = fits.getheader(Configuration.DARKS_DIRECTORY + image_list[0])
                else:
                    bias_header = fits.getheader(Configuration.BIAS_DIRECTORY + image_list[0])

                bias_header['BIAS_COMB'] = 'mean'
                bias_header['NUM_BIAS'] = nfiles
                bias_header["BIASSUB"] = 'Y'

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             bias_image, bias_header, overwrite=True)
        else:
            bias_image = fits.getdata(Configuration.MASTER_DIRECTORY + file_name, 0)

        return bias_image

    @staticmethod
    def mk_flat(image_directory):
        """ This function will make the master flat frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination

        :return - The bias frame is returned and written to the master directory
        """

        if os.path.isfile(Configuration.MASTER_DIRECTORY + "_flat" + ".fits") == 0:

            # read in the bias file
            bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY)

            # get the image list
            image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

            # here is the 'holder'
            hold_data = np.zeros((Configuration.AXS_X, Configuration.AXS_Y))

            Utils.log("Making the flat frame using a mean.", "info")

            # now loop through the images
            for idx, image in enumerate(image_list):
                # read in the flat file
                flat_tmp = fits.getdata(image_directory + "\\" + image, 0)
                if np.ndim(flat_tmp) > 2:
                    flat_tmp = flat_tmp[0]

                # read in the image directly into the block_hold
                hold_data += flat_tmp - bias

            # median the mini-images into one large image
            flat_image = hold_data / float(len(image_list))
            nflat_image = flat_image / np.median(flat_image)

            # pull the header information from the first file of the set
            flat_header = fits.getheader(image_directory + "\\" + image_list[0])
            flat_header['comb_typ'] = 'mean'
            flat_header['norm_pix'] = 'mean'
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
            fits.writeto(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat" + Configuration.FILE_EXTENSION,
                         nflat_image, flat_header, overwrite=True)
        else:
            nflat_image = fits.getdata(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat" + Configuration.FILE_EXTENSION, 0)

        return nflat_image

    @staticmethod
    def sky_subtract(img, header, sky_write='N'):
        """  The function breaks the image into bxs x bxs sections to save memeory, then cleans each section
        the sections are then recombined and smoothed to remove the transitions. The residual image is then subtracted
        from the image and the header is updated appropriately.

        :parameter img - The image to be cleaned
        :parameter header - The header object to be updated
        :parameter sky_write - Y/N if you want to write the residual background for de-bugging

        :return img_sub, header - The cleaned image and updated header file
        """

        # use the sampling space to make the appropriate size vectors
        lop = 2 * Configuration.PIX

        # size holder for later
        sze = int((Configuration.AXIS_X / Configuration.PIX) * (Configuration.AXIS_Y / Configuration.PIX) +
                  (Configuration.AXIS_X / Configuration.PIX) + (Configuration.AXIS_Y / Configuration.PIX) + 1)

        # generate an empty matrix to hold the background
        res = np.zeros(shape=(Configuration.AXIS_Y, Configuration.AXIS_X))

        # calculate the sky statistics
        sky_mean, sky_median, sky_sig = astropy.stats.sigma_clipped_stats(img, sigma=2.5)

        # create holder arrays for good and bad pixels
        x = np.zeros(shape=sze)
        y = np.zeros(shape=sze)
        v = np.zeros(shape=sze)
        s = np.zeros(shape=sze)
        nd = int(0)

        # begin the sampling of the "local" sky value
        for jj in range(0, Configuration.AXIS_X + Configuration.PIX, Configuration.PIX):
            for kk in range(0, Configuration.AXIS_Y + Configuration.PIX, Configuration.PIX):
                il = np.amax([jj - lop, 0])
                ih = np.amin([jj + lop, Configuration.AXIS_X - 1])
                jl = np.amax([kk - lop, 0])
                jh = np.amin([kk + lop, Configuration.AXIS_Y - 1])
                c = img[jl:jh, il:ih]

                # select the median value with clipping
                lsky_mean, lsky, ssky = astropy.stats.sigma_clipped_stats(c, sigma=2.5)

                x[nd] = np.amin([jj, Configuration.AXIS_X - 1])  # determine the pixel to input
                y[nd] = np.amin([kk, Configuration.AXIS_Y - 1])  # determine the pixel to input
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
        xi = np.linspace(0, Configuration.AXIS_X - 1, Configuration.AXIS_X)
        yi = np.linspace(0, Configuration.AXIS_Y - 1, Configuration.AXIS_Y)
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
    def bias_subtract(img, header, dark):
        """ This function will subtract a bias frame
        :parameter img - The image to de-bias / de-dark
        :parameter header - The image header file
        :parameter dark - If Y then the dark frame will be subtracted rather than a bias

        :return bias_sub, header - The updated image and header """

        if dark == 'Y':
            file_name = "dark" + Configuration.FILE_EXTENSION
        else:
            file_name = 'bias' + Configuration.FILE_EXTENSION

        # read in the bias frame and subtract
        bias = fits.getdata(Configuration.MASTER_DIRECTORY + file_name)

        # subtract the bias from the image
        bias_sub = img - bias

        # update the header
        if dark == 'Y':
            header['dark_sub'] = 'Y'
        else:
            header['bias_sub'] = 'Y'

        return bias_sub, header

    @staticmethod
    def flat_divide(img, header):
        """ This function will divide a flat field.
        :parameter img - The image to flatten
        :parameter header - The image header file

        :return flat_div, header - The updated image and header """
        # read in the flat frame
        flat = fits.getdata(Configuration.MASTER_DIRECTORY + 'flat' + Configuration.FILE_EXTENSION)

        # subtract the bias from the image
        flat_div = img / flat

        # update the header
        header['flat_div'] = 'Y'

        return flat_div, header

    @staticmethod
    def mk_nme(file, difference_image='N', sky_subtract='N', bias_subtract='N',
               flat_divide='N', alignment='N', dark_subtract="N"):
        """ This function will create the appropriate name for the file based on while steps are taken.
        :argument file - The string with the file name
        :argument difference_image - Y/N if image subtraction occured
        :argument sky_subtract - Y/N if sky subtraction was taken
        :argument bias_subtract - Y/N if a bias is subtracted
        :argument flat_divide - Y/N if a flat field is divided
        :argument alignment - Y/N if the image is aligned
        :argument dark_subtract - Y/N if the image was dark subtracted

        :return file_name - A string with the new file name
        """
        # make a new name for the file based on which actions are taken
        nme_hld = file.split('.')

        # if everything is N then the file name is the original filename
        file_name = file

        # update the file name with a 'd' if at the differencing step
        if difference_image == 'Y':
            file_name = nme_hld[0] + 'd' + Configuration.FILE_EXTENSION

        # otherwise...
        if difference_image == 'N':
            # update the name to be appropriate for what was done to the file
            # nothing occurs
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N') \
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + Configuration.FILE_EXTENSION
            # bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_b' + Configuration.FILE_EXTENSION
            # flat
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_f' + Configuration.FILE_EXTENSION
            # align
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_a' + Configuration.FILE_EXTENSION
            # sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_s' + Configuration.FILE_EXTENSION
            # dark subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_k' + Configuration.FILE_EXTENSION

            # bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bf' + Configuration.FILE_EXTENSION
            # bias and align only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_ba' + Configuration.FILE_EXTENSION
            # bias and sky_subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bs' + Configuration.FILE_EXTENSION
            # bias and flat and align only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bfa' + Configuration.FILE_EXTENSION

            # flat and align only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_fa' + Configuration.FILE_EXTENSION
            # flat and sky_subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_fs' + Configuration.FILE_EXTENSION
            # flat and align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_fsa' + Configuration.FILE_EXTENSION

            # align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_sa' + Configuration.FILE_EXTENSION

            # flat and sky subtract and bias only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bfs' + Configuration.FILE_EXTENSION

            # all steps taken
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkfsa' + Configuration.FILE_EXTENSION

            # dark and bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bk' + Configuration.FILE_EXTENSION
            # dark and flat only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N') \
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_kf' + Configuration.FILE_EXTENSION
            # dark and align only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_ka' + Configuration.FILE_EXTENSION
            # dark and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_ks' + Configuration.FILE_EXTENSION
            # dark and bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkf' + Configuration.FILE_EXTENSION
            # dark and bias and align only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bka' + Configuration.FILE_EXTENSION
            # dark and bias and sky subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bks' + Configuration.FILE_EXTENSION
            # dark and bias and flat and align only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkfa' + Configuration.FILE_EXTENSION
            # dark and bias and flat and sky subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkfs' + Configuration.FILE_EXTENSION
            # dark and flat and align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_kfas' + Configuration.FILE_EXTENSION

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
