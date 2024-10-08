from config import Configuration
from libraries.priority import Priority
from libraries.dbaccess import DBaccess
from astroquery.mast import Catalogs
from libraries.utils import Utils

# get the toros field information
toros_fields = Priority.toros_field_generator(1.19)

for idx, row in toros_fields.iterrows():

    if (idx >= 100) & (idx < 1000):
        # create the string useful for query_region
        field = str(row.ra) + " " + str(row.dec)

        # select the columns we want to import into the data table
        columns = ["toros_field_id", "source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
                   "teff_val", "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error"]

        # run the query
        Utils.log('Querying MAST for all stars within the toros field: ' + str(row.toros_field_id), 'info')
        catalog_data = Catalogs.query_region(field, radius=Configuration.SEARCH_DIST/1.5, catalog="Gaia")
        Utils.log('Query finished. ' + str(len(catalog_data)) + ' stars found.', 'info')

        # add the toros field to the catalog data
        catalog_data['toros_field_id'] = row.toros_field_id

        # pull out the necessary columns
        star_list = catalog_data[columns]

        # generate the table to insert the data into
        table_num = int(idx / 1000) + 1
        if table_num < 10:
            table_nme = '00' + str(table_num)
        elif (table_num >= 10) & (table_num < 100):
            table_nme = '0' + str(table_num)
        else:
            table_nme = str(table_num)

        Utils.log("Starting to insert targets for field: " + str(row.toros_field_id) + ".", "info")
        # generate SQL command and send to database
        DBaccess.insert_into_gaia_table(star_list, table_nme)

        Utils.log("Targets inserted for field: " + str(row.toros_field_id) + ".", "info")
