""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'tolar'
    RAW_FILE_EXTENSION = '.fits.fz'
    FILE_EXTENSION = '.fits'

    # update for different data products
    STAR = 'NGC2243'
    RA = 97.3950000
    DEC = -31.2819444
    DATES = ['20220107', '20220108', 'injection']
    REF_DATE = DATES[0]
    DATE = DATES[2]
    PHOTOMETRY = 'PSF'
    APERTURE_SHAPE = 'circ'

    # steps to skip
    CLEAN_SKIP = 'N'
    WRITE_SKY = 'N'
    CALIBRATE_SKIP = 'N'
    MASTER_SKIP = 'N'
    DIFFERENCE_SKIP = 'N'
    PHOTOMETRY_SKIP = 'N'
    LIGHTCURVE_SKIP = 'N'
    COLOR_SKIP = 'N'

    # get image information
    AXS_X = 3352  # 1676
    AXS_Y = 2532  # 1266
    AXS = 2048
    GAIN = 0.380  # in e-/ADU
    FOV = 5  # size of the image in arc minutes
    SEARCH_DIST = FOV / 60.0
    EXP_TIME = 60.

    # update the differencing information, primarily the number of stars to use, and the kernel size
    KRNL = 2  # kernel size 2 * KNRL + 1
    STMP = 11  # stamp size ot use 2 * STMP + 1
    ORDR = 0  # order of the kernel to use, 0 is stationary, 1 or 2 is spatially varying
    NRSTARS = 1000  # number of stars used to solve for kernel
    BRIGHT_STARS = 20000  # the top stars to search for in kernel stars
    KERNEL_LIMIT = 0.5  # the maximum allowable offset in zeropoint in magnitudes
    AXS_LIMIT = 50  # the number of pixel close to the edge of the frame to use
    RMS_LOW_LIMIT = 0.005  # the lower limit on precision to use for the kernel stars
    RMS_UP_LIMIT = 0.02  # the upper limit on precision to use for the kernel stars

    # update sky subtraction specific information
    PIX_BOX = 128
    PIX = 64

    # update the image axes to work for a PIX_BOX setting
    X_CUT = int(AXS_X % PIX_BOX)  # work to make it divisible by a PIX x PIX box
    Y_CUT = int(AXS_Y % PIX_BOX)
    X_CENT = int(AXS_X / 2)  # get the center with respect to the old image
    Y_CENT = int(AXS_Y / 2)
    AXIS_X = int(AXS_X - X_CUT)  # get the size of the new image
    AXIS_Y = int(AXS_Y - Y_CUT)

    # if the image is divisible by the box, then don't worry
    if (Y_CUT != 0) | (X_CUT != 0):
        CUT_IMAGE = 'Y'
    else:
        CUT_IMAGE = 'N'

    # a photometry configuration
    FWHM = 15.  # fwhm of the image
    THRESHOLD = 5.  # the threshold for a source above the background

    # aperture information
    CIRC_APER_SIZE = 10  # circular aperture

    # aperture annulus for the sky background automatically determined from the main aperture
    CIRC_ANNULI_INNER = CIRC_APER_SIZE + 2
    CIRC_ANNULI_OUTER = CIRC_APER_SIZE + 4

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = "/home/ryan.oelkers/Development/toros/"
    ALERTS_DIRECTORY = WORKING_DIRECTORY + 'alerts/'
    ANALYSIS_DIRECTORY = WORKING_DIRECTORY + 'analysis/'
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs/'
    QUERIES_DIRECTORY = WORKING_DIRECTORY + 'queries/'
    CODE_DIFFERENCE_DIRECTORY = WORKING_DIRECTORY + 'difference/'

    # input paths for data etc
    DATA_DIRECTORY = WORKING_DIRECTORY + "data/"
    DARKS_DIRECTORY = DATA_DIRECTORY + "darks/"
    BIAS_DIRECTORY = DATA_DIRECTORY + "bias/"
    FLATS_DIRECTORY = DATA_DIRECTORY + "flats/"
    RAW_DIRECTORY = DATA_DIRECTORY + "raw/"
    CLEAN_DIRECTORY = DATA_DIRECTORY + "clean/"
    MASTER_DIRECTORY = DATA_DIRECTORY + "master/"
    CENTROID_DIRECTORY = MASTER_DIRECTORY + "centroids/"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc/" + DATE + "/"
    DIFFERENCED_DIRECTORY = DATA_DIRECTORY + "diff/"
    DIFFERENCED_DATE_DIRECTORY = DATA_DIRECTORY + "diff/" + DATE + "/"
    CLEAN_DATE_DIRECTORY = CLEAN_DIRECTORY + DATE + "/"

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, MASTER_DIRECTORY, LIGHTCURVE_DIRECTORY,
                   CENTROID_DIRECTORY, RAW_DIRECTORY, FLATS_DIRECTORY, BIAS_DIRECTORY, DARKS_DIRECTORY,
                   DIFFERENCED_DIRECTORY, DIFFERENCED_DATE_DIRECTORY, CLEAN_DATE_DIRECTORY]

    # BROKER CONFIGURATION SPECIFICS
    LISTEN_NED_WAIT = 1

    # observing conditions
    SEEING = 0.93  # assumes 2 pix FWHM

    # sky brightness at TOLAR
    SKY = [22.1, 21.1, 20.1, 18.7]

    # bandpass values in nm
    BP = [147, 141, 147, 147]

    # cwl
    CWL = [473.5, 638.5, 775.5, 922.5]

    # telescope information
    TOROS_MIRROR_D = 0.610  # m
    TOROS_MIRROR_R = TOROS_MIRROR_D / 2  # nm

    # CCD information
    READOUT_NOISE = 5.0  # electrons / pixel
    CCD_QE = 0.85
    FILTER_QE = 0.9
    TELESCOPE_SEC_QE = 0.96
    TELECSCOPE_PRI_QE = 0.96
    VIGNETTING = 0.756
    ATMOSPHERE_QE = [0.8, 0.9, 0.9, 0.9]

    # total throughput needs to be multiplied by atmospheric quantum efficiency
    TOTAL_THROUGHPUT = CCD_QE * FILTER_QE * TELESCOPE_SEC_QE * TELECSCOPE_PRI_QE * VIGNETTING

    # telescope information
    PIXEL_SIZE = 0.468  # arcsec per pixel
    NUM_PIXELS = 10560  # pixels per side
    TOROS_DEC_LIMIT = 33.84  # declination limit of the telescope in degrees

    # force toros field generation?
    FEILD_GENERATION = 'N'
    FIELD_SIZE = 1.19

    TOROS_LONGITUDE = -67.32833333
    TOROS_LATITUDE = -24.62055556
    TOROS_ELEVATION = 4650
    MOON_DISTANCE = 60

    EXPOSURE_TIME = 300
    EXPOSURE_TIME_DAY = EXPOSURE_TIME / 60. / 60. / 24.
    NUM_EXPOSURES = 1
    READ_TIME = 90
    OVERHEAD = 30
    TOTAL_EXPOSURE = (EXPOSURE_TIME + READ_TIME) * NUM_EXPOSURES + OVERHEAD