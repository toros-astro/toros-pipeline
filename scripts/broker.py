from libraries.utils import Utils
from libraries.priority import Priority
from libraries.alerts import Alerts
from config import Configuration
from libraries.dbaccess import DBaccess
import time

# the broker will run at all times, there should always be at least 10 fields in the observation queue, if the number of
# fields in the observing queue drops below 10, then the broker will update the list with the appropriate amount
while True:

    # pull in the main toros fields, try the database, and if no connection, pull from local file
    try:
        sql_query = DBaccess.get_query(Configuration.QUERIES_DIRECTORY + 'field_list.sql')
        survey_fields = DBaccess.query_torosdb(sql_query, Configuration.MACHINE)
    except ConnectionError:
        Alerts.alert_toros_team("DATABASE")
        Utils.log("No connection to TOROSdb! Pulling from legacy file.", "info")
        survey_fields = Priority.toros_field_generator(Configuration.FIELD_SIZE)

    Utils.log("TOROS broker is selecting fields for observation...", "info")

    # get the best field for observation
    obs_field = Priority.toros_field_selector(survey_fields)

    # send fields to observing table in TOROSdb
    # pull in the main toros fields, try the database, and if no connection, pull from local file
    try:
        sql_query = DBaccess.get_update_obs_query(Configuration.QUERIES_DIRECTORY + 'upd_obs_queue.sql', obs_field)
        survey_fields = DBaccess.update_torosdb(sql_query, Configuration.MACHINE)
    except ConnectionError:
        Alerts.alert_toros_team("DATABASE")
        Utils.log("No connection to TOROSdb! Field observations cannot be updated!", "info")

    # send the field to observe to the observing table
    Utils.log("Field to observe now is: " + obs_field['toros_field_id'], "info")
    Utils.log("Pausing for 5 minutes to allow for observations.", "info")
    time.sleep(300)
