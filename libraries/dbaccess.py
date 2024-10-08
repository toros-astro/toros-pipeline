""" This class is used for basic functions not specific to this code base such as:
 logging, file writing, and testing."""
from libraries.utils import Utils
import psycopg2
import pandas as pd
import os
from sqlalchemy import create_engine


class DBaccess:

    @staticmethod
    def get_query(path):
        """ This function will read in the appropriate SQL query verbatim. Assumes no update to query!!

        :parameter path - The location of the SQL file

        :return A string with the SQL
        """

        # read in teh appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        return sql_cmd

    @staticmethod
    def get_ins_obs_query(path, toros_field):
        """ This function will generate the sql to update the observing queue table in TOROSdb.

        :parameter path - The location of the SQL file
        :parameter toros_field - A data frame with the line to update the observing table

        :return A string with the SQL
        """
        # read in the appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        sql_cmd = sql_cmd.replace('%(toros_field_id)s', toros_field.loc[0, 'toros_field_id'])
        sql_cmd = sql_cmd.replace('%(ra)s', str(toros_field.loc[0, 'ra']))
        sql_cmd = sql_cmd.replace('%(dec)s', str(toros_field.loc[0, 'dec']))
        sql_cmd = sql_cmd.replace('%(program)s', toros_field.loc[0, 'program'])
        sql_cmd = sql_cmd.replace('%(exposure_time)s', str(toros_field.loc[0, 'exposure_time']))
        sql_cmd = sql_cmd.replace('%(cadence)s', str(toros_field.loc[0, 'cadence']))
        sql_cmd = sql_cmd.replace('%(ephemeris)s', str(toros_field.loc[0, 'ephemeris']))
        sql_cmd = sql_cmd.replace('%(period)s', str(toros_field.loc[0, 'period']))
        sql_cmd = sql_cmd.replace('%(obs_date)s', str(toros_field.loc[0, 'obs_date']))

        return sql_cmd

    @staticmethod
    def insert_into_gaia_table(star_list, db_table_num):
        """ This function will generate the sql to update the observing queue table in TOROSdb.

        :parameter star_list - The gaia data queried from MAST to update on TOROSdb
        :parameter db_table_num - The table number for the Gaia insert

        :return - no data is returned, but the GAIA rows are sent to the TOROSdb
        """

        # set up the connection object based on which computer is being used
        conn = psycopg2.connect(host="127.0.0.1", port=5432, database="torosdb", user="vaquero", password="b44rnd0r")

        # set up the cursor object
        cur = conn.cursor()

        # do block inserts of 1000 rows
        for idx in range(0, len(star_list), 1000):

            # make sure we don't go out of bounds
            try:
                block_list = star_list[idx:idx+1000]
            except:
                block_list = star_list[idx::]

            # set up the argument string based on the TOROS query
            args_str = ','.join(cur.mogrify("(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", x).decode('utf-8')
                                for x in block_list.iterrows())

            # execute and commit the insert command
            cur.execute("INSERT INTO torosdb.toros_star_list_" + db_table_num +
                        " (toros_field_id, source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, "
                        "teff_val, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error) VALUES " +
                        args_str)
            conn.commit()

        # shut it down
        cur.close()
        conn.close()

        return

    @staticmethod
    def get_update_obs_query(path, toros_field):
        """ This function will generate the sql to update the observing queue table in TOROSdb.

        :parameter path - The location of the SQL file
        :parameter toros_field - A data frame with the line to update the observing table

        :return A string with the SQL
        """
        # read in the appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        sql_cmd = sql_cmd.replace('%(toros_field_id)s', toros_field['toros_field_id'])
        sql_cmd = sql_cmd.replace('%(ra)s', str(toros_field['ra']))
        sql_cmd = sql_cmd.replace('%(dec)s', str(toros_field['dec']))
        sql_cmd = sql_cmd.replace('%(program)s', toros_field['program'])
        sql_cmd = sql_cmd.replace('%(exposure_time)s', str(toros_field['exposure_time']))
        sql_cmd = sql_cmd.replace('%(cadence)s', str(toros_field['cadence']))
        sql_cmd = sql_cmd.replace('%(ephemeris)s', str(toros_field['ephemeris']))
        sql_cmd = sql_cmd.replace('%(period)s', str(toros_field['period']))

        return sql_cmd

    @staticmethod
    def get_update_field_query(path, toros_field):
        """ This function will generate the sql to update the observing queue table in TOROSdb.

        :parameter path - The location of the SQL file
        :parameter toros_field - A data frame with the line to update the observing table

        :return A string with the SQL
        """
        # read in the appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        sql_cmd = sql_cmd.replace('%(toros_field_id)s', toros_field.loc[0, 'toros_field_id'])

        return sql_cmd

    @staticmethod
    def query_torosdb(sql_cmd, machine):
        """ This function will query the TOROS database and return a data frame based on the query used.
        This function will confirm a dumped .csv file does not exist prior to accessing the database.

        :parameter sql_cmd - The sql query to use
        :parameter machine - The current computer being used for the program

        :return df - a data frame with the query results
        """
        # set up the connection object based on which computer is being used
        if machine == 'tolar':
            url = 'postgresql+psycopg2://vaquero:b44rnd0r@127.0.0.1:5432/torosdb'

        # set up connection
        db = create_engine(url)
        conn = db.connect()

        # generate the data frame with the queried results
        df = pd.read_sql_query(sql_cmd, conn)

        # shut it down
        conn.close()

        return df

    @staticmethod
    def update_torosdb(sql_cmd, machine):
        """ This function will update data in the TOROS database. No data is returned.

        :parameter sql_cmd - The sql query to use
        :parameter machine - The current computer being used for the program

        :return df - nothing is returned, but the SQL statement is executed
        """
        # set up the connection object based on which computer is being used
        if machine == 'tolar':
            conn = psycopg2.connect(host="127.0.0.1",
                                    port=5432,
                                    database="torosdb",
                                    user="vaquero",
                                    password="b44rnd0r")

        # set up the cursor object
        cur = conn.cursor()

        # execute and commit the insert command
        cur.execute(sql_cmd)
        conn.commit()

        # shut it down
        cur.close()
        conn.close()

        return

    @staticmethod
    def insert_torosdb(sql_cmd, machine):
        """ This function will insert data into the TOROS database. No data is returned.

        :parameter sql_cmd - The sql query to use
        :parameter machine - The current computer being used for the program

        :return df - nothing is returned, but the SQL statement is executed
        """
        # set up the connection object based on which computer is being used
        if machine == 'tolar':
            conn = psycopg2.connect(host="127.0.0.1",
                                    port=5432,
                                    database="torosdb",
                                    user="vaquero",
                                    password="b44rnd0r")

        # set up the cursor object
        cur = conn.cursor()

        # execute and commit the insert command
        cur.execute(sql_cmd)
        conn.commit()

        # shut it down
        cur.close()
        conn.close()

        return

    @staticmethod
    def query_n_dump_torosdb(sql_cmd, out_path, file_name, machine):
        """ This function will query the TOROS database, return a data frame, and dump the result. This function
        will confirm a dumped .csv file does not exist prior to accessing the database.

        :parameter sql_cmd - The sql query to use
        :parameter out_path - The output for the file
        :parameter file_name - The desired filename
        :parameter machine - The current computer being used for the program

        :return df - a data frame with the query results
        """
        # set up the connection object based on whether you are on tessdev
        if machine == 'tolar':
            conn = psycopg2.connect(host="127.0.0.1",
                                    port=5432,
                                    database="torosdb",
                                    user="vaquero",
                                    password="b@rnd00r")

        # set up the cursor object
        cur = conn.cursor()

        if os.path.isfile(out_path + file_name) == 1:
            Utils.log("Legacy file found, no reason to query TOROSdb", "info")
            # read in from a file
            df = pd.read_csv(out_path + file_name, index_col=0)
            Utils.log("CSV read complete.", "info")

        if os.path.isfile(out_path + file_name) == 0:
            Utils.log("No legacy file found. Querying TOROSdb...", "info")
            # generate the data frame with the queried results
            df = pd.read_sql_query(sql_cmd, conn)

            Utils.log("Query complete. Dumping to .csv file " + out_path + file_name, "info")
            # dump the file to csv
            df.to_csv(out_path + file_name)

            Utils.log("Dump complete.", "info")

        # shut it down
        cur.close()
        conn.close()

        return df

