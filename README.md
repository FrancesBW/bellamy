# BELLAMY: A cross-matching package for the cynical astronomer

BELLAMY is a cross-matching algorithm that aims to match all sources in the supplied target catalogue to sources in a reference catalogue. BELLAMY utilises not only the position of a source on the sky, but also the flux data, to determine the most probable match in the reference catalog to the target source. This two variable approach means BELLAMY is a more rigorous cross-matcher than the matching tools supplied by TOPCAT or astropy alone. Additionally, BELLAMY attempts to undo any spatial distortion that may be affecting the target catalogue, by creating a model of the offsets of matched sources which is then applied to unmatched sources. Furthermore, BELLAMY provides the option for a simple flux correction across the target catalogue, in the case that there is some large scale calibration error affecting the target catalogue.

## Version support and required packages:

BELLAMY is compatible with Python 2.7, but currently untested with Python 3. 

The Python 2.7 version requires the following packages:
* numpy
* astropy (version 3.0 and above only support Python 3, so use versions 2.x)
* scipy
* argparse
* matplotlib

Note: Python 2.7 will become deprecated in late 2019/early 2020 so upgrade to Python 3.6 or greater before this time would be the goal.

## Installation:

To install BELLAMY you can use one of the following methods:

* Install directly from github
```
$ pip install git+https://github.com/FrancesBW/bellamy.git
```
* Download from github and install via the setup file
```
$ git clone https://github.com/FrancesBW/bellamy.git
$ cd bellamy
$ python setup.py install
```

## Usage:

BELLAMY only requires 3 inputs to run a cross match: a target catalogue, a reference catalogue and the reference catalogue format. The simplest execution on the command line would be:
```
$ bellamy --target TARGET_FILENAME --reference REFERENCE_FILENAME --refsurvey REFERENCE_SURVEY_NAME
```

If the reference catalogue is from a supported survey, calling the survey name is sufficient (e.g. if using the GLEAM catalogue enter '--refsurvey GLEAM' on the command line). If you are using a reference catalogue from a non-supported survey, or of your own design, you can edit the 'reference_catalogue_format.txt' file to reflect the column names in your reference catalogue for necessary data used by BELLAMY. Once you have edited the file you can call BELLAMY on the command line with '--refsurvey reference_catalogue_format.txt' and BELLAMY will read this file.

Currently supported catalogues are:
* GLEAM (see http://www.mwatelescope.org/gleam for instructions on downloading GLEAM catalogue)
* TGSS (see http://tgssadr.strw.leidenuniv.nl/doku.php for instructions on downloading TGSS catalogue)

### Required Inputs:

* ```--target TARGET_FILENAME``` filename/path for the target catalogue.

* ```--reference REFERENCE_FILENAME``` filename/path for the reference catalogue.

* ```--refsurvey REFERENCE_SURVEY_NAME``` provide the format file for the target catalogue (see Formats section for guide). If target catalogue is a supported catalogue, provide survey name.

### Settings:

* ```--norestrict``` will turn off the signal-to-noise restriction that is applied on the first iteration of cross-matching. Only target sources above an SNR cutoff (where SNR is peak_flux/local_rms) are passed to the matching function. This restriction is applied by default to provide the most accurate initial offset model. 

* ```--overwrite``` will allow output files to overwrite files of the same name. If those files already exist and overwrite is active, the logger will also display a warning on which files will be overwritten. 

* ```--nofluxmodel``` will turn off the modelling of a flux adjustment factor across the target catalogue field of view. 

* ```--plot``` will turn on plotting of measured and modelled offsets of sources matched between the target and reference catalogues. It also turns on plotting for the flux model applied to the target catalogue.

* ```--nofluxmatch``` will stop peak flux of target and reference sources being compared for matching. If flux match is turned off, BELLAMY returns a nearest neighbour match with modelled adjustments.

* ```--debug``` will turn on debug mode.

* ```--writelog``` will write out the log shown in the terminal to 'bellamy_log_(date)\_(time).log'. 

### Optional Inputs:

* ```--tarfreq FREQUENCY ``` provide the frequency of the target catalogue. This is only required if both flux-matching is on (as is the default), and the reference catalogue has flux measurements for more than one frequency. If the reference catalogue has flux measurements for only one fequency or has no specified frequency, then BELLAMY proceeds by assuming it is comparable to the target catalogue. 

* ```--outformat FILE_EXTENSION``` provide the desired file extension/format for the output tables. See astropy.table.Table docuemntation for supported writing formats (http://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers). Default format is 'fits'.

* ```--singlepercentile PERCENTILE``` provide the desired percentile cutoff for target sources with single reference match candidates within the search radius. Only target sources that meet this confidence level will be accepted as matches and used in the offset model. Default value is 0.95.

* ```--multipercentile PERCENTILE``` provide the desired normalised percentile cutoff for target sources with multiple reference match candidates within the search radius. Only target sources that meet this confidence level will be accepted as matches and used in the offset model. Default value is 0.60.

* ```--snr SNR_CUTOFF``` provide the desired signal-to-noise (SNR, defined as peak_flux/local_rms) cutoff for target sources matched in the first iteration. To remove this cutoff use '--norestrict'. Default value is 10.

* ```--modeldeg DEGREE``` provide the desired degree of model used in flux model (integer between 1 and 5). Modelled using scipy.interpolate.LSQBivariateSpline (see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.LSQBivariateSpline.html). Default value is 1.

* ```--tarformat FORMAT_FILENAME``` provide the format file for the target catalogue (see Formats section for guide). Default is detailed in Formats section.

* ```--filesuffix SUFFIX``` provide desired suffix to be appended for all output files (except the log file). Note a '\_' will be added between default filename and suffix, so there is no need to include this in the suffix. Default is None.

* ```--tarres TAR_RES``` provide a single value for the image resolution of the target catalogue 


## Formats:

As catalogue formats are wildly varied, it becomes necessary to have some sort of formatting guide for BELLAMY, to allow it to read in any catalogue. 

**Target catalogue**

By default it is assumed that the target catalogue is constructed by the source finding program Aegean (see https://github.com/PaulHancock/Aegean for more information). This points BELLAMY to search for columns called 'ra', 'dec' 'err_ra' etc. in the target catalogue (full list shown in target catalogue defaults column below). If you target catalogue was constructed by Aegean, or happens to have exactly the same column names, then you do not need to specify a format file for the target catalogue. If your target catalogue was not constructed by Aegean, or has different column names, edit 'target_catalogue_format.txt' (found in the format_templates directory) with the exact column names for the data required. Note that this file does not ask for the catalogue frequency (you can specify this using ```--tarfreq```). This means that BELLAMY cannot support target catalogues with multiple frequency data, like it can with reference catalogues. If you have data for multiple frequencies in your target catalogue, you must specify the column names for the specific frequency you want to match. 

**Reference catalogue**

BELLAMY does not have a default reference catalogue format, but it does have supported catalogues, for which the formats are already known. So far only the GLEAM catalogue (76-227 MHz) and TGSS catalogue (150 MHz) are supported but it is possible to add support for more large survey catalogues. If your reference catalogue is a supported one, you need only specify the survey acronym/name for the ```--refsurvey``` argument. If your reference catalogue is not supported, edit 'reference_catalogue_format.txt' (found in the format_templates directory) with the exact column names for the data required. 

This file (unlike the one for the target catalogue), asks you to specify the catalogue frequency. You can specify a single frequency, or a frequency range. If your reference catalogue has data for mulitple frequencies, and has data for the exact same frequency as your target catalogue, then you have to specify that as a single frequency. **For example**, if my target catalogue is measured at 200 MHz, and my reference catalogue has data for 150, 200, 250 and 300 MHz, I would specify in the formatting file that the catalogue frequency is 200. 

If the frequency of your target catalogue falls in between measured frequencies in your reference catalogue, then you can specify a frequency range. If you specify a frequency range, you MUST specify it as the nearest frequencies in the reference catalogue, to the target catalogue frequency. **For example**, if my target catalogue is measured at 145 MHz, and my reference catalogue has data for 130, 140, 150, 160 and 170 MHz, I would specify in the formatting file that the catalogue frequency is 140-150. 

The table below shows the default target catalogue column names and the supported catalogue column names. 

| Data type (typical units)        | target catalogue defaults | GLEAM         | TGSS        |
| -------------------------------- | ------------------------- | ------------- | ----------- |
| RA (J2000, deg)                  | ra                        | RAJ2000       | RA          |
| DEC (J2000, deg)                 | dec                       | DEJ2000       | DEC         |
| error in RA (deg)                | err_ra                    | err_RAJ2000   | E_RA        |
| error in DEC (deg)               | err_dec                   | err_DEJ2000   | E_DEC       |
| Beam/psf semimajor axis (arcsec) | psf_a                     | psf_a         |             |
| Beam/psf semiminor axis (arcsec) | psf_b                     | psf_b         |             |
| Source semimajor axis (arcsec)   | a                         | a             | Maj         |
| Source semiminor axis (arcsec)   | b                         | b             | Min         |
| Source position angle (deg)      | pa                        | pa            | PA          |
| Peak flux (Jy/Beam)              | peak_flux                 | peak_flux     | Peak_flux   |
| error in Peak flux (Jy/Beam)     | err_peak_flux             | err_peak_flux | E_Peak_flux |
| Integrated flux (Jy)             | int_flux                  | int_flux      | Total_flux  |
| RMS noise (Jy/Beam)              | local_rms                 | local_rms     |             |
| Unique source name               | uuid                      | Name          | Source_name |
| Catalogue frequency (MHz)        |                           | 076-227       | 150         |
| Frequency prefix or suffix       |                           | suffix        |             |

Note: The last row only applies to reference catalogues with flux data for multiple frequencies. In this case, columns pertaining to flux or psf etc. may have a prefix or suffix indicating the frequency the column refers to. If this is the case for your reference catalogue, you can edit the 'reference_catalog_format.txt' file and add either prefix or suffix at the end of the 'frequency_prefix_or_suffix=' line. 

**For example:**

GLEAM does not have a column called 'peak_flux' like its format suggests, but rather multiple columns called 'peak_flux_076', 'peak_flux_084', 'peak_flux_092' and so on. By specifying that the peak flux data is stored in a column generally called 'peak_flux' and specifying the frequency is a suffix on these column names, allows the algorithm to search for columns called 'peak_flux' with the frequency appended on the end (with or without a '\_'). 

## How it works

The first step taken by BELLAMY is to gauge how well the fluxes match between the reference and target catalogues. It uses sources with a signal to noise greater than 10, does a quick and dirty nearest neighbour match and compares the fluxes of the matched pairs. It then fits a model (2D bivariate spline) of the adjustment factor needed to bring the target catalogue fluxes in line with the reference catalogue and applies it (a rough form of calibration).

Then, the adjusted target catalogue and reference catalogued are passed to a cross matching function, which searches within a 10' radius of a target source to find any possible matching refernce sources. Then using the difference in position and flux between the target source and each possible reference match, each match is assigned a value (between 0 and 1) to determine how 'likely' it is to be correct. The likelihood is calculated from a Gaussian distribution like so:

![screenshot_2019-02-13 icrar summer project](https://user-images.githubusercontent.com/38412724/52688347-132d2580-2f91-11e9-81e8-c128f06e7cad.png)

The uncertainties used in these likelihood calculations are calculates like so:

![screenshot_2019-02-13 icrar summer project 1](https://user-images.githubusercontent.com/38412724/52688693-6f447980-2f92-11e9-88fc-7867fd32d5b1.png)

Target sources with only one possible match will be accepted as good matches if their raw likelihood meets the 'singlepercentile' threshold. For target sources with multiple possible matches, a normalised likelihood will be calculated (so that all likelihoods add up to 1) and the most likely reference source will be accepted as a good match if the normalised likelihood meets the 'multipercentile' threshold. Target sources with accepted matches are removed from the target list left to be matched. Target sources with no accepted matches, remain in the target list left to be matched. 

Accepted matches are then used to model the offsets of the target and reference sources. A Rbf model is used, which is then applied to the remaining target sources. The measured offsets and the model based on these, can then be plotted and output. 

This process is repeated until there are no additional accepted matches to be made. Note that by default, the first iteration of cross-matching runs with only target sources with a SNR of 10 or more, to ensure the initial offset model is as accurate as can be. Once no new matches are made, the algorithm returns all target sources with their most likely match and the raw likelihood of that match.

The algorithm will output 4 tables of the specified file format (fits by default) with the following names (with a suffix if the user has specified one):

* 'cross_matched_table': each row contains the positions, dimensions and unique identifiers of each target source and its matched reference source. It also contains the raw and normalised (where applicable) likelihoods of the match, as well as the number of possible matches that were considered for the target source within the search radius.

* 'leftover_reference_catalogue': contains the data for the reference sources that weren't matched or deemed the most likely to match a target source.

* 'leftover_target_catalogue': contains the data for the target sources that were unable to be matched, with their positions altered in accordance with the Rbf offset model.

* 'leftover_unwarped_target_catalogue': contains the data for the target sources that were unable to be matched, with their original positions.

