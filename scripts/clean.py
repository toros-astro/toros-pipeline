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
import subprocess
import shutil

class Clean:

    @staticmethod
    def clean_images(image_clip='N', bias_subtract="N", dark_subtract="N",
                     flat_divide='N', sky_subtract="N", plate_solve='N'):
        """ This is the main function script to clean multiple images, alternatively clean_img can be used to clean
        a single image.
        :parameter image_clip - Y/N if you want to clip any excess from the images (default = N)
        :parameter bias_subtract - Y/N if you want to subtract the bias from the images (default = N)
        :parameter dark_subtract -Y/N if you want to subtract a dark frame from the science images (default = N)
        :parameter flat_divide - Y/N if you want to flatten the images (default = N)
        :parameter sky_subtract - Y/N if you want to subtract the sky from the images (default = N)
        :parameter plate_solve -Y/N if you want to plate_solve the image (default = N)

        return no value is returned, the values images from in_path are cleaned and deposited in out_path
        """
        st = time.time()  # clock started

        # get the file list for all dates the FIELD was observed
        Utils.log("Getting file list...", "info")
        files, date_dirs = Utils.get_all_files_per_field(Configuration.RAW_DIRECTORY,
                                                         Configuration.FIELD,
                                                         Configuration.RAW_FILE_EXTENSION)

        # make the output directories (clean and diff)
        output_dirs = []

        for dte in date_dirs:
            output_dirs.append(Configuration.DATA_DIRECTORY + "clean/" + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + "clean/" + dte + "/" + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + "diff/" + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + "diff/" + dte + "/" + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + "review/" + dte) 
            output_dirs.append(Configuration.DATA_DIRECTORY + "review/" + dte + "/" + Configuration.FIELD)

        Utils.create_directories(output_dirs)
        # break if there are no files
        if len(files) == 0:
            Utils.log("No .fits files found for " + Configuration.FIELD + "!" +  ". Breaking...",
                      "debug")
            return()

        Utils.log("Starting to clean " + str(len(files)) + " images.", "info")
        for idx, file in enumerate(files):

            # make a new name for the file based on which actions are taken
            file_name = Preprocessing.mk_nme(file,
                                             difference_image='N',
                                             image_clip = image_clip,
                                             bias_subtract = bias_subtract,
                                             dark_subtract = dark_subtract,
                                             flat_divide  = flat_divide,
                                             sky_subtract = sky_subtract,
                                             plate_solve = plate_solve)

            # only create the files that don't exist
            # skip if in review pile
            review_file_name = file_name.replace("/clean/", "/review/")
            if os.path.isfile(file_name) == 1:
                Utils.log("Image " + file_name +
                          " already exists. Skipping for now...", "info")
            
            if os.path.isfile(review_file_name):
                Utils.log("Image " + file_name +
                          " is marked for review. Skipping for now...", "info")
            # if the image does not exist then clean
            if os.path.isfile(file_name) == 0 and os.path.isfile(review_file_name) == 0:

                # clean the image
                clean_img, header, bd_flag = Clean.clean_img(file, image_clip,
                                                             bias_subtract, dark_subtract, flat_divide,
                                                             sky_subtract, plate_solve)

                # write out the file
                if bd_flag == 0:
                    #fits.writeto(file_name,
                    #             clean_img, header, overwrite=True)

                    # print an update to the cleaning process
                    Utils.log("Cleaned image written as " + file_name + ".", "info")
                else:
                    Utils.log(file_name + " is a bad image. Not written.", "info")

            Utils.log(str(len(files) - idx - 1) + " images remain to be cleaned.",  "info")

        fn = time.time()  # clock stopped
        Utils.log("Imaging cleaning complete in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

    @staticmethod
    def clean_img(file, image_clip='Y', bias_subtract='N',
                  dark_subtract="N", flat_divide='N', sky_subtract="N", plate_solve="N"):
        """ This function is the primary script to clean the image, various other functions found in this class
        can be found in the various libraries imported.

        :parameter  file - The file name of the image you would like to clean
        :parameter image_clip - Y/N if you want to clip the image (default = Y)
        :parameter sky_subtract - Y/N if you want to subtract the sky from the image (default = Y)
        :parameter bias_subtract - Y/N if you want to remove a bias frame (default = N)
        :parameter flat_divide - Y/N if you want to flatten the image (default = N)
        :parameter dark_subtract -Y/N if you want to subtract the dark frame (default = N)
        :parameter plate_solve -Y/N if you want to plate solve the image (default = N)
        """

        Utils.log("Now cleaning " + file + ".", "info")

        # read in the image
        img, header = fits.getdata(file, header=True)

        # get cleaned file path
        clean_file_name = Preprocessing.mk_nme(file, difference_image='N',
                                         image_clip=image_clip,
                                         sky_subtract=sky_subtract,
                                         bias_subtract=bias_subtract,
                                         flat_divide=flat_divide,
                                         dark_subtract=dark_subtract,
                                         plate_solve=plate_solve)

        # bias subtract if necessary
        if bias_subtract == 'Y':
            st = time.time()
            bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY, bias_overwrite='N', combine_type='mean')
            img, header = Preprocessing.bias_subtract(img, header)
            fn = time.time()
            Utils.log("Image bias corrected in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping bias correction....", "info")

        # dark subtract if necessary
        if dark_subtract == 'Y':
            st = time.time()
            dark = Preprocessing.mk_dark(Configuration.DARKS_DIRECTORY, overwrite_dark='N', combine_type='mean')
            img, header = Preprocessing.dark_subtract(img, header)
            fn = time.time()
            Utils.log("Image dark corrected in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping dark correction....", "info")

        if image_clip == 'Y':
            st = time.time()
            img, header = Preprocessing.clip_image(img, header)
            fn = time.time()
            Utils.log("Image over scan removed in " + str(np.around(fn - st, decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping over scan removal...", "info")

        # flat divide if necessary
        if flat_divide == 'Y':
            st = time.time()
            flat = Preprocessing.mk_flat(Configuration.FLATS_DIRECTORY, combine_type='median')
            img, header = Preprocessing.flat_divide(img, header)
            fn = time.time()
            Utils.log("Image flattened in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping image flattening....", "info")

        # sky subtract if necessary
        if sky_subtract == 'Y':
            st = time.time()
            # the background sample size is set to pix x pix pixels, and bxs x bxs sub images
            # this should not be hard coded...update for later

            Utils.log("A background box of " + str(Configuration.PIX) + " x " + str(Configuration.PIX) +
                      " will be used for background subtraction.", "info")
            fits.writeto(Configuration.DATA_DIRECTORY+"debug/debug_beforeskysub.fits", img, header, overwrite=True)
            img, header = Preprocessing.sky_subtract(img, header, Configuration.WRITE_SKY)
            fn = time.time()
            Utils.log("Sky subtracted in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

        else:
            Utils.log("Skipping sky subtraction...", "info")

        if plate_solve == 'Y':
            st = time.time()
            Utils.log("Now plate solving and correcting the header.", "info")
            fits.writeto(Configuration.DATA_DIRECTORY+"debug/debug_platesolve.fits", img, header, overwrite=True)
            img, header = Preprocessing.correct_header(img, header)
            fits.writeto(clean_file_name+".temp" , img, header, overwrite=True)

            # set output directory for astrometry
            output_dir = os.path.dirname(clean_file_name)

            astrometry_command = "solve-field --scale-units arcsecperpix --scale-low 0.48 --scale-high 0.5 --no-plots --temp-axy --index-xyls none --match none --rdls none --solved none --corr none --dir " + output_dir + " --new-fits " + clean_file_name + " " + clean_file_name + ".temp"

            subprocess.run(astrometry_command, shell=True)

            # Clean up astrometry output, remove .temp and .wcs files
            astrometry_extensions = [ ".temp", ".wcs"]
            
            # First check that astrometry generated plate solved file (.fits)
            if os.path.exists(clean_file_name):
                # if file exists then remove extra files
                for extension in astrometry_extensions:
                    file_to_remove = clean_file_name + extension
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                        Utils.log(f"File {file_to_remove} deleted.", "info")
            else:
                # if file does not exist then move the temp file out of clean for review
                temp_file_name = clean_file_name.replace("/clean/","/review/")
                shutil.move(clean_file_name+".temp", temp_file_name)
                Utils.log(f"Plate Solve Failed. Moving file to review directory for human inspection.", "info")
            
            fn = time.time()
            Utils.log("Plate solve step completed in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        Utils.log("Cleaning finished.", "info")

        bd_flag = 0

        return img, header, bd_flag
