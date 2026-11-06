# Test data

| File | Contents | Source |
| --- | --- | --- |
| `cat_*.tsv`, `cat_atlas.txt` | minimal samples in VizieR TSV format (UCAC4, CMC-15, PPMXL, Gaia, ATLAS) for the catalog reader | synthetic |
| `horizons_hindemith.txt` | ephemerides of (5157) Hindemith at the mid-exposure instants of the example (`cpp/gui/assets/example`) | JPL Horizons API (query in the file header) |
| `orbits_mpcorb.dat` | elements of (1) Ceres, (3) Juno, (433) Eros and (5157) Hindemith in the MPCORB.DAT column layout | JPL Horizons API |
| `orbits_comets.dat` | elements of 2P/Encke and 9P/Tempel 1 in the MPC COMET.DAT column layout | JPL Horizons API |

## Orbital elements (`orbits_*.dat`)
Heliocentric osculating elements, ecliptic and equinox J2000 (ICRF), epoch 2020-12-17.0 TDB (JD 2459200.5,
`K20CH`), obtained on 2026-11-06 from https://ssd.jpl.nasa.gov/api/horizons.api with:

    format=text COMMAND='<target>' OBJ_DATA=YES MAKE_EPHEM=YES EPHEM_TYPE=ELEMENTS CENTER='500@10'
    REF_PLANE=ECLIPTIC REF_SYSTEM=ICRF TLIST=2459200.5 TIME_TYPE=TDB CSV_FORMAT=YES OUT_UNITS=AU-D

with `<target>` = `5157;`, `1;`, `3;`, `433;`, `DES=2P;CAP;NOFRAG`, `DES=9P;CAP;NOFRAG`. H and G (asteroids) and M1, k1
(comets; K = 2.5·k1) come from the OBJ_DATA of the same query; solutions JPL#71 (5157), JPL#48 (1), JPL#144 (3), JPL#659
(433), K273/27 (2P), K284/2 (9P). For comets the perihelion date Tp is in TDB. The `test_orbits` reference positions
(EPHEM_TYPE=OBSERVER, QUANTITIES='1,9', CENTER '500' or 'I41@399') are in the test itself.

Credit: NASA/JPL Solar System Dynamics Group, Horizons On-Line Ephemeris System (public data).
