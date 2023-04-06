""" This class is used to house all of the alert type of scripts: listen for alerts, send out alerts, play alerts..."""
import requests
import pandas as pd
import smtplib
from time import sleep
from config import Configuration
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from libraries.utils import Utils
from playsound import playsound
import logging
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


class Alerts:

    @staticmethod
    def check_for_gcn(alert_directory):
        """ This function will check to see if an alert file is in the ALERTS directory. If so, it will send an alert
        text and email, and send the fields to the broker to prioritize.

        :parameter alert_directory - The directory where the various LIGO, GRB, etc alerts live
        """

        # pull out the alert files
        files = Utils.get_file_list(alert_directory, '_alert.txt')

        if len(files) > 0:
            Utils.log('Alert found, querying NED for list of galaxies.', 'info')
            # read in the event information deposited by the listener
            event_name = pd.read_csv(alert_directory + files[0])

            # query NED
            galaxies = Alerts.query_ned(event_name)
        else:
            # no alerts found so ignore
            Utils.log('No alerts found. Ignoring request.', "info")

        return files

    @staticmethod
    def play_alert_sound(play_count=10):
        """ This function will play an alert sound on the host computer when an event is detected. By default it will
        play the sound 10 times, but for some alerts it will play continuously until a user presses a key.

        :parameter play_count - By default the program will play an alert sound 10 times, but it can be any number of
        times, or set to 'inf'. If set to 'inf', a user must press a key to stop the alert sound.

        :return Nothing is returned, but a sound is played for the desired length of time.
        """
        # start the playcount iteration check at 0
        sound_iter = 0

        # run through the playcount
        while sound_iter < play_count:
            playsound('broker/mixkit-critical-alarm-1004.wav')
            sound_iter = sound_iter + 1
            Utils.log("ALERT! A new OOT event has been detected, please confirm broker is prioritizing correctly.",
                      "warning")
            try:
                Utils.log("Confirm you are here with Crtl+C", "warning")
                sleep(10)
            except KeyboardInterrupt:
                Utils.log("Thanks for confirming you are here, please check broker is prioritizing correctly.",
                          "warning")
                break
            else:
                Utils.log("No observer present, alerting again.", "warning")
        if sound_iter == play_count:
            Utils.log("No observer present. Alerting team to contact observer directly.", "warning")
            Alerts.alert_toros_team("observer")

        Utils.log("Resuming normal broker operations.", "info")
        return

    @staticmethod
    def query_ned(ligo_name, voe_number='latest', file_format='csv'):
        """ This function will query the NED API for the best ligo galaxies to observe in the search radius.

        :parameter ligo_name - The name of the ligo event, ex GW170817
        :parameter voe_number - The number for hte VOEvent serial number, default = latest
        :parameter file_format - The file type to download, default = csv

        :return ned_result - A pandas data frame containing the field information is returned
        """

        # set up the URL for the download
        url = 'https://ned.ipac.caltech.edu/uri/NED::GWFglist/' + file_format + '/' + ligo_name + '/' + voe_number

        # attempt to pull down the data, and if broken, then skip
        try:
            # pull and write the request
            r = requests.get(url, allow_redirects=True)
            open(Configuration.ANALYSIS_DIRECTORY + ligo_name + '_NED_result.' + file_format, 'wb').write(r.content)

            # read in the csv file
            ned_result = pd.read_csv(Configuration.ANALYSIS_DIRECTORY + ligo_name + '_NED_result.' + file_format)

            # pull out the galaxies
            return ned_result
        except:

            Utils.log("Query to NED produced no result.", "info")

            return

    @staticmethod
    def alert_toros_team(alert_type='test', event_name='test'):
        """ This function send text messages and emails to the TOROS team when an alert is detected.

        :parameter - alert_type - The example alert type for the text to send (example GW for gravitational wave), the
        default is to simply send a test alert.
        :parameter - event_name - The name of the event for the header of the message

        :return Nothing is returned, however a message is sent to devices
        """

        # TOROS alerts email requirements
        email = "toros.alerts@yahoo.com"
        pas = "ncnrlqwdthofhoch"
        smtp = "smtp.mail.yahoo.com"
        port = 465

        # list of phone numbers and emails to send the alert message to
        sms_gateway = ['2672619337@vtext.com', 'ryan.oelkers@tamu.edu']
        # Moises(phone), Richard (phone), Ryan (phone), Ryan (email), Lucas (phone)
        # ['9564342399@tmomail.net', '2674218771@tmomail.net', '2672619337@vtext.com',
        # 'ryan.oelkers@tamu.edu', '9797390438@vtext.com']

        # login to the email server using the appropriate credentials
        server = smtplib.SMTP_SSL(smtp, port)
        server.ehlo()
        server.login(email, pas)

        # set up the message based on the type of alert
        msg = MIMEMultipart()
        msg['From'] = email
        msg['To'] = ", ".join(sms_gateway)

        # gravitational wave
        if alert_type == 'LVC Preliminary':
            msg['Subject'] = "New Preliminary LIGO-GW ALERT: " + event_name + "\n"
            body = "Ole! A new GW-event has been detected by LIGO. Check TOROS queue and results.\n"
        elif alert_type == 'LVC Initial':
            msg['Subject'] = "New Initial LIGO-GW ALERT: " + event_name + "\n"
            body = "Ole! A new GW-event has been detected by LIGO. Check TOROS queue and results.\n"
        # gamma-ray-burst
        elif alert_type == 'GRB':
            msg['Subject'] = "New Gamma-ray Burst ALERT: " + event_name + "\n"
            body = "Ole! A new GRB-event has been detected. Check TOROS queue and results.\n"
        # supernovae
        elif alert_type == 'SNE':
            msg['Subject'] = "New Supernovae ALERT: " + event_name + "\n"
            body = "Ole! A new Supernova has been detected. Check TOROS queue and results.\n"
        # missing observer
        elif alert_type == 'OBSERVER':
            msg['Subject'] = "Observer not present\n"
            body = "A new TOROS alert has been detected, but the observer may not be present.\n"
        # retraction
        elif alert_type == 'LVC Retraction':
            msg['Subject'] = "Event: " + event_name + " retracted\n"
            body = "The previous alert has been retracted.\n"
        # no database connection
        elif alert_type == 'DATABASE':
            msg['Subject'] = "Database connection failed\n"
            body = "TOROS computer <" + Configuration.MACHINE + "> cannot connect to the database.\n"
        # test the system
        else:
            msg['Subject'] = "TOROS Alert Test for event: " + event_name + "\n"
            body = "This is a test of the TOROS alert message system.\n"

        # convert message to text and send
        msg.attach(MIMEText(body, 'plain'))
        sms = msg.as_string()
        server.sendmail(email, sms_gateway, sms)

        # clean up
        server.quit()

        # log that the message was sent
        Utils.log("The following alert message was sent to the TOROS team: " + body, "info")

        return
