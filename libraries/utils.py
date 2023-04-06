""" This class is used for basic functions not specific to this code base such as:
            logging, writing to files, making directories, and generating file lists"""
from config import Configuration
import logging
import os


class Utils:

    @staticmethod
    def log(statement, level):
        """ The logger function to log all activity from the program to both the screen, and a file.

        :argument statement: A string which shows want needs to be logged
        :argument level:- The type of statement: info, debug, etc

        :return - Nothing is returned, but the log is updated and possibly the log is printed to the screen
        """

        # create the logger
        logging.basicConfig(format="%(asctime)s - %(levelname)s: %(message)s",
                            filename=Configuration.LOG_DIRECTORY + "toros.log",
                            filemode='a')
        logger = logging.getLogger()

        if not getattr(logger, 'handler_set', None):
            logger.setLevel(logging.DEBUG)

            # create console handler and set level to debug
            ch = logging.StreamHandler()
            ch.setLevel(logging.DEBUG)

            # create the formatter
            formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")

            # add formatter to ch
            ch.setFormatter(formatter)

            # add ch to logger
            logger.addHandler(ch)

            # 'set' Handler
            logger.handler_set = True

        if level == 'info':
            logger.info(statement)
        if level == 'debug':
            logger.debug(statement)
        if level == 'warning':
            logger.warning(statement)
        if level == 'error':
            logger.error(statement)
        if level == 'critical':
            logger.critical(statement)

    @staticmethod
    def write_txt(path, typ, line):
        """ This function will either write or append a text file, primarily used for results or logging.

        :parameter path: The location for the file
        :parameter typ: The write type either 'w', or 'a'
        :parameter line: The line you want to write to the file

        :return: Nothing is returned, the file is simply written or appended.

        """

        # open the file for writing
        file = open(path, typ)

        # write the line to the file
        file.write(line)

        # close the file
        file.close()

    @staticmethod
    def get_file_list(path, file_ext):
        """ This function will return the files in a given directory without their extension.
        :parameter path - A string with the path the file.
        :parameter file_ext - The file type to make a list, generally something like *.fits.

        :return files_no_ext - A list of files without the extension."""

        # get the files in the given path directory with the given file extension
        file_list = [f for f in os.listdir(path) if f.endswith(file_ext)]

        # sort based on the number of the image, first taken image will be first
        file_list.sort(key=len)

        return file_list

    @staticmethod
    def check_alert_exists(path, file_ext):
        """ This function will check to see if an alert file exists, and if so, it will return the alert, delete the
        file, and then flag the user."""

        # set initial exist flag
        exst = 0
        # check for alert files

    @staticmethod
    def create_directories(directory_list):
        """ This function will check for each directory in directory list, and create it if it doesn't already exist

        :parameter directory_list - The list of directories to create

        :return - Nothing is returned but directories are created if necessary
        """

        for path in directory_list:
            # if the path does not exist then create it!
            if os.path.exists(path) is False:
                os.mkdir(path)
                Utils.log(path + ' created.', 'info')
