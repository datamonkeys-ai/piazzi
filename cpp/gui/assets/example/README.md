# Piazzi example data

Three exposures of one ZTF quadrant (field 565, CCD 15, quadrant 2) taken on 2021-02-06 between 05:10 and 06:08 UTC
(filters r, i, g; 30 s each), with the main-belt asteroid **(5157) Hindemith** (V ≈ 17.6, 0.34″/min) near the
centre, plus the matching Gaia DR3 reference stars. "Open example" in the Piazzi apps reduces this set fully offline.

**This is example data — not your observation.** The observatory code is `XXX` and the apps block MPC submission while
the example is loaded. Never submit these measurements to the Minor Planet Center.

| File | Content |
| --- | --- |
| `ztf_hindemith_0.fits` … `_2.fits` | 356 × 356 px (6′ × 6′) central crops of the IRSA ZTF science-image cutouts (`ztf_20210206215787_000565_zr_c15_o`, `…232836_…_zi_…`, `…255880_…_zg_…`); pixel values unchanged, `NAXIS`/`CRPIX` adjusted |
| `gaia_dr3.tsv` | Gaia DR3 sources (RP < 19.1) in a 17′ box around RA 07h11m43.5s, Dec +22°24′18″ (from the NOIRLab Astro Data Lab TAP mirror of `gaia_dr3.gaia_source`, in VizieR TSV layout) |
| `example.cfg` | Session settings: focal length 3063 mm, 15 µm pixels, field centre, MPC code `XXX` |

Command line equivalent:

```sh
piazzi --cfg example.cfg --catalog gaia_dr3.tsv --center "07 11 43.534" "+22 24 17.78" ztf_hindemith_*.fits
```

## Credits

**ZTF.** Based on observations obtained with the Samuel Oschin Telescope 48-inch and the 60-inch Telescope at the
Palomar Observatory as part of the Zwicky Transient Facility project. ZTF is supported by the National Science
Foundation under Grants No. AST-1440341 and AST-2034437 and a collaboration including current partners Caltech, IPAC,
the Oskar Klein Center at Stockholm University, the University of Maryland, University of California, Berkeley, the
University of Wisconsin at Milwaukee, University of Warwick, Ruhr University, Cornell University, Northwestern
University and Drexel University. Operations are conducted by COO, IPAC, and UW.

The images were obtained from the NASA/IPAC Infrared Science Archive (IRSA), which is funded by the National
Aeronautics and Space Administration and operated by the California Institute of Technology. ZTF public data
(DOI: 10.26131/IRSA539).

**Gaia DR3.** This work has made use of data from the European Space Agency (ESA) mission Gaia
(https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC,
https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions,
in particular the institutions participating in the Gaia Multilateral Agreement. Credit: ESA/Gaia/DPAC.
Catalogue query served by the NSF NOIRLab Astro Data Lab.
