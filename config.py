""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'mcast'
    RAW_FILE_EXTENSION = '.fits'
    FILE_EXTENSION = '.fits'

    # update for different data products
    fields = [ 'FIELD_0e.001', 'FIELD_2b.022', 'FIELD_30.000']
    FIELD = fields[0] #'FIELD_2b.022' #'FIELD_42.074' #'FIELD_2b.022' #'FIELD_36.007' # 'FIELD_28.01c'  # 'FIELD_2f.0d2'  # # FIELD_30.000
    RA = 6.207 #19.1811417 #52.759 #148.198 #10.18868 # 46.636  # 307.317 # # 0.723
    DEC = -72.292 #-59.9845556  #-35.611 #-6.519 #-21.69714 # -39.405  # -30.551 #  # -29.286
    DATE = '2024-09-30'
    PHOTOMETRY = 'PSF'
    APERTURE_SHAPE = 'circ'

    # is there a transient you want the light curve for?
    TRANSIENT_LC = 'N'
    TRANSIENT_NAME = 'SN2024aaur'#'AT2024xsq' # 'AT2024xhm'  # 'AT2024ykx' #
    TRANSIENT_RA = 148.614629 #10.1958012 # 46.9968525  # 307.422951333 #
    TRANSIENT_DEC = -6.952958 #-21.9258309  # -38.93038375  # -30.269894 #

    # steps to skip
    CLEAN_SKIP = 'N'
    WRITE_SKY = 'Y'
    CALIBRATE_SKIP = 'N'
    MASTER_SKIP = 'N'
    DIFFERENCE_SKIP = 'N'
    PHOTOMETRY_SKIP = 'N'
    LIGHTCURVE_SKIP = 'N'
    COLOR_SKIP = 'Y'

    # telescope information
    PIXEL_SIZE = 0.4959  #  0.47  # arcsec per pixel
    NUM_PIXELS = 10560  # pixels per side
    TOROS_DEC_LIMIT = 26.66  # declination limit of the telescope in degrees
    FOV = (PIXEL_SIZE * NUM_PIXELS) / 3600.
    SEARCH_DIST = FOV
    EXP_TIME = 300
    GAIN = 0.380  # in e-/ADU

    # get image information
    AXS_X_RW = 12000
    AXS_Y_RW = 10600
    AXS_X = 10560
    AXS_Y = 10560
    AXS = 10560

    # update the differencing information, primarily the number of stars to use, and the kernel size
    KRNL = 2  # kernel size 2 * KNRL + 1
    STMP = 15  # stamp size ot use 2 * STMP + 1
    ORDR = 0  # order of the kernel to use, 0 is stationary, 1 or 2 is spatially varying
    NRSTARS = 3000  # number of stars used to solve for kernel
    BRIGHT_STARS = 20000  # the top stars to search for in kernel stars
    KERNEL_LIMIT = 0.5  # the maximum allowable offset in zeropoint in magnitudes
    AXS_LIMIT = 100  # the number of pixel close to the edge of the frame to use
    RMS_LOW_LIMIT = 0.005  # the lower limit on precision to use for the kernel stars
    RMS_UP_LIMIT = 0.02  # the upper limit on precision to use for the kernel stars

    # update sky subtraction specific information
    PIX = 220

    # a photometry configuration
    FWHM = 15.  # fwhm of the image
    THRESHOLD = 5.  # the threshold for a source above the background

    # aperture information
    APER_SIZE = 16  # circular aperture

    # aperture annulus for the sky background automatically determined from the main aperture
    ANNULI_INNER = APER_SIZE + 2
    ANNULI_OUTER = APER_SIZE + 4

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = "/Volumes/datadrive/working_directory/"
    ALERTS_DIRECTORY = WORKING_DIRECTORY + 'alerts/'
    ANALYSIS_DIRECTORY = WORKING_DIRECTORY + 'analysis/'
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs/'
    QUERIES_DIRECTORY = WORKING_DIRECTORY + 'queries/'
    CODE_DIFFERENCE_DIRECTORY = WORKING_DIRECTORY + 'difference/'

    # input paths for data etc
    DATA_DIRECTORY = "/Volumes/datadrive/coldharbor/"
    DARKS_DIRECTORY = DATA_DIRECTORY + "darks/"
    BIAS_DIRECTORY = DATA_DIRECTORY + "bias/"
    FLATS_DIRECTORY = DATA_DIRECTORY + "flats/"
    RAW_DIRECTORY = DATA_DIRECTORY + "raw/"
    CLEAN_DIRECTORY = DATA_DIRECTORY + "clean/"
    REVIEW_DIRECTORY = DATA_DIRECTORY + "review/"
    MASTER_DIRECTORY = DATA_DIRECTORY + "master/"
    CALIBRATION_DIRECTORY = DATA_DIRECTORY + "calibration/"
    CENTROID_DIRECTORY = MASTER_DIRECTORY + "centroids/"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc/"
    DIFFERENCED_DIRECTORY = DATA_DIRECTORY + "diff/"
    CLEAN_DATE_DIRECTORY = CLEAN_DIRECTORY + DATE + "/"
    REVIEW_DATE_DIRECTORY = REVIEW_DIRECTORY + DATE + "/"

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY, CALIBRATION_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, MASTER_DIRECTORY, LIGHTCURVE_DIRECTORY,
                   CENTROID_DIRECTORY, RAW_DIRECTORY, FLATS_DIRECTORY, BIAS_DIRECTORY, DARKS_DIRECTORY,
                   DIFFERENCED_DIRECTORY, CLEAN_DATE_DIRECTORY, REVIEW_DATE_DIRECTORY, CODE_DIFFERENCE_DIRECTORY]

    # BROKER CONFIGURATION SPECIFICS
    LISTEN_NED_WAIT = 1

    # observing conditions
    SEEING = 0.93  # assumes 2 pix FWHM

    # sky brightness at TOLAR in SDSS griz
    SKY = [22.1, 21.1, 20.1, 18.7]

    # SDSS griz bandpass values in nm (width of the filter)
    BP = [147, 141, 147, 147]

    # SDSS griz central wavelength
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

    # force toros field generation?
    FIELD_GENERATION = 'N'

    # The Felix Aguilar Observatory is more Southern than Tolar
    TOROS_LONGITUDE = -69.3265  # -67.32833333
    TOROS_LATITUDE = -31.8023  # -24.62055556
    TOROS_ELEVATION = 2420
    UTC = -3
    MOON_DISTANCE = 60

    EXPOSURE_TIME = 300
    EXPOSURE_TIME_DAY = EXPOSURE_TIME / 60. / 60. / 24.
    NUM_EXPOSURES = 1
    READ_TIME = 90
    OVERHEAD = 30
    TOTAL_EXPOSURE = (EXPOSURE_TIME + READ_TIME) * NUM_EXPOSURES + OVERHEAD
