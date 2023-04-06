INSERT INTO torosdb.%(gaia_table_name)s
(toros_field_id, source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, teff_val,
parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error)
VALUES ('%(toros_field_id)s', %(source_id)s, %(ra)s, %(dec)s, %(phot_g_mean_mag)s, %(phot_bp_mean_mag)s,
%(phot_rp_mean_mag)s, %(teff_val)s, %(parallax)s, %(parallax_error)s, %(pmra)s, %(pmra_error)s, %(pmdec)s, %(
pmdec_error)s);
