""" This script will clean the FFI images, it has been significantly reduced to only
show specific steps. The majority of the functional library can be found in the libraries
directory in utils, fits, photometry, image, and preprocessing libraries."""
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from config import Configuration
import os
from astropy.io import fits
import time
import numpy as np


class Clean:

    @staticmethod
    def clean_images(sky_subtract="N", bias_subtract="N", flat_divide='N', alignment='N', dark_subtract="N"):
        """ This is the main function script to clean multiple images, alternatively clean_img can be used to clean
        a single image.
        :parameter sky_subtract - Y/N if you want to subtract the sky from the images (default = Y)
        :parameter bias_subtract - Y/N if you want to subtract the bias from the images (default = N)
        :parameter flat_divide - Y/N if you want to flatten the images (default = N)
        :parameter alignment - Y/N if you want to align the images (default = N)
        :parameter dark_subtract -Y/N if you want to subtract a dark frame from the science images (default = N)

        return no value is returned, the values images from in_path are cleaned and deposited in out_path
        """
        st = time.time()  # clock started

        # get the file list
        Utils.log("Getting file list...", "info")
        files = Utils.get_file_list(Configuration.RAW_DIRECTORY + Configuration.DATE + '/',
                                    Configuration.RAW_FILE_EXTENSION)

        # break if there are no files
        if len(files) == 0:
            Utils.log("No .fits files found in " + Configuration.RAW_DIRECTORY + Configuration.DATE + ". Breaking...",
                      "debug")
            return()
        
        Utils.log("Starting to clean " + str(len(files)) + " images.", "info")
        for idx, file in enumerate(files):

            # make a new name for the file based on which actions are taken
            file_name = Preprocessing.mk_nme(file, 'N', sky_subtract, bias_subtract,
                                             flat_divide, alignment, dark_subtract)

            # only create the files that don't exist
            if os.path.isfile(Configuration.CLEAN_DIRECTORY + Configuration.DATE + '/' + file_name) == 1:
                Utils.log("Image " + Configuration.CLEAN_DIRECTORY + file_name +
                          " already exists. Skipping for now...", "info")

            # if the image does not exist then clean
            if os.path.isfile(Configuration.CLEAN_DIRECTORY + Configuration.DATE + '/' + file_name) == 0:

                if (alignment == 'Y') & (idx == 0) & (Configuration.DATE == Configuration.REF_DATE):
                    # clean the image
                    clean_img, header, bd_flag = Clean.clean_img(file,
                                                                 Configuration.RAW_DIRECTORY + Configuration.REF_DATE,
                                                                 sky_subtract, bias_subtract, flat_divide,
                                                                 'N', dark_subtract)
                else:
                    # clean the image
                    ref_files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE + '/',
                                                    Configuration.FILE_EXTENSION)
                    # clean_img, header, bd_flag = Clean.clean_img(file,
                    #                                             Configuration.CLEAN_DIRECTORY +
                    #                                             Configuration.REF_DATE + '/' + ref_files[0],
                    #                                             sky_subtract, bias_subtract, flat_divide,
                    #                                             'Y', dark_subtract)
                    rr_img, header, bd_flag = Clean.clean_img(Configuration.RAW_DIRECTORY +
                                                              Configuration.DATES[1] + '/' + file,
                                                              Configuration.RAW_DIRECTORY +
                                                              Configuration.DATES[1] + '/' + file,
                                                              sky_subtract, bias_subtract, flat_divide,
                                                              'N', dark_subtract)
                    fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'fake.fits', rr_img, header, overwrite=True)
                    clean_img, header, bd_flag = Clean.clean_img(Configuration.RAW_DIRECTORY + Configuration.DATE +
                                                                 '/' + file,
                                                                 Configuration.ANALYSIS_DIRECTORY + 'fake.fits',
                                                                 sky_subtract, bias_subtract, flat_divide,
                                                                 'Y', dark_subtract)

                # write out the file
                if bd_flag == 0:
                    fits.writeto(Configuration.CLEAN_DIRECTORY + Configuration.DATE + '/' + file_name,
                                 clean_img, header, overwrite=True)

                    # print an update to the cleaning process
                    Utils.log("Cleaned image written as " +
                              Configuration.CLEAN_DIRECTORY + Configuration.DATE + '/' + file_name + ".", "info")
                else:
                    Utils.log(file_name + " is a bad image. Not written.", "info")

            Utils.log(str(len(files) - idx - 1) + " images remain to be cleaned.",  "info")

        fn = time.time()  # clock stopped
        Utils.log("Imaging cleaning complete in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

    @staticmethod
    def clean_img(file, ref_path, sky_subtract="N", bias_subtract='N',
                  flat_divide='N', alignment='N', dark_subtract="N"):
        """ This function is the primary script to clean the image, various other functions found in this class
        can be found in the various libraries imported.

        :parameter  file - The file name of the image you would like to clean
        :parameter ref_path - The path to the reference frame
        :parameter sky_subtract - Y/N if you want to subtract the sky from the image (default = Y)
        :parameter bias_subtract - Y/N if you want to remove a bias frame (default = N)
        :parameter flat_divide - Y/N if you want to flatten the image (default = N)
        :parameter alignment - Y/N if you want to align the image (default = N)
        :parameter dark_subtract -Y/N if you want to subtract the dark frame (default = N)
        """

        Utils.log("Now cleaning " + file + ".", "info")

        # read in the image
        img, header = fits.getdata(file, 1, header=True)
        if np.ndim(img) > 2:
            img = img[0]

        # bias subtract if necessary
        if bias_subtract == 'Y':
            st = time.time()
            bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY, dark='N', combine_type='mean')
            img, header = Preprocessing.bias_subtract(img, header, dark='N')
            fn = time.time()
            Utils.log("Image bias corrected in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping bias correction....", "info")

        # dark subtract if necessary
        if dark_subtract == 'Y':
            st = time.time()
            dark = Preprocessing.mk_bias(Configuration.DARKS_DIRECTORY, dark='Y', combine_type='mean')
            img, header = Preprocessing.bias_subtract(img, header, dark='Y')
            fn = time.time()
            Utils.log("Image dark corrected in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping dark correction....", "info")

        # flat divide if necessary
        if flat_divide == 'Y':
            st = time.time()
            flat = Preprocessing.mk_flat(Configuration.FLATS_DIRECTORY)
            img, header = Preprocessing.flat_divide(img, header)
            fn = time.time()
            Utils.log("Image flattened in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping image flattening....", "info")

        if Configuration.CUT_IMAGE == 'Y':
            Utils.log("Clipping the image from " + str(Configuration.AXS_X) + "x" + str(Configuration.AXS_Y) + " to " +
                      str(Configuration.AXIS_X) + "x" + str(Configuration.AXIS_Y), "info")
            img, header = Preprocessing.remove_overscan(img, header)
        else:
            header['CLIPPED'] = 'N'
            Utils.log("Skipping image clipping...", "info")

        # sky subtract if necessary
        if sky_subtract == 'Y':
            st = time.time()
            # the background sample size is set to pix x pix pixels, and bxs x bxs sub images
            # this should not be hard coded...update for later

            Utils.log("A background box of " + str(Configuration.PIX) + " x " + str(Configuration.PIX) +
                      " will be used for background subtraction.", "info")

            img, header = Preprocessing.sky_subtract(img, header, Configuration.WRITE_SKY)
            fn = time.time()
            Utils.log("Sky subtracted in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

        if sky_subtract == 'N':
            Utils.log("Skipping sky subtraction...", "info")

        # align the image if necessary
        if alignment == 'Y':
            st = time.time()
            img, header = Preprocessing.align_img(img, header, ref_path)
            fn = time.time()
            Utils.log("Image aligned in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

        if alignment == 'N':
            Utils.log("Skipping image alignment....", "info")

        Utils.log("Cleaning finished.", "info")
        bd_flag = 0

        return img, header, bd_flag
