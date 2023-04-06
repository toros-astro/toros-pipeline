import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from libraries.utils import Utils
import logging
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


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
    if alert_type == 'GW':
        msg['Subject'] = "New LIGO-GW ALERT: " + event_name + "\n"
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
    elif alert_type == 'observer':
        msg['Subject'] = "Observer not present\n"
        body = "A new TOROS alert has been detected, but the observer may not be present.\n"
    # retraction
    elif alert_type == 'RETRACTION':
        msg['Subject'] = "Event: " + event_name + " retracted\n"
        body = "The previous alert has been retracted.\n"
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