SELECT ticlink.pk, ticlink.tessmag, ticlink.ra, ticlink.dec,
    q3c_dist(ticlink.ra, ticlink.dec, %(cen_ra)s, %(cen_de)s) AS rdist
FROM tic.ticlink
WHERE q3c_radial_query(ticlink.ra, ticlink.dec, %(cen_ra)s, %(cen_de)s, %(mx_dist)s)
ORDER BY q3c_dist(ticlink.ra, ticlink.dec, %(cen_ra)s, %(cen_de)s)