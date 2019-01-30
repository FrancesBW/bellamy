# MACCAS: Multi Attribute Cross-matcher with Correction for wArped Sky

MACCAS is a cross-matching algorithm that aims to match all sources in the supplied target catalogue to sources in a reference catalogue. MACCAS utilises not only the position of a source on the sky, but also the flux data, to determine the most probable match in the reference catalog to the target source. This two variable approach means MACCAS is a more rigorous cross-matcher than the matching tools supplied by TOPCAT or astropy alone. Additionally, MACCAS attempts to undo any spatial distortion that may be affecting the target catalogue, by creating a model of the offsets of matched sources which is then applied to unmatched sources. Furthermore, MACCAS provides the option for a simple flux correction across the target catalogue, in the case that there is some large scale calibration error affecting the target catalogue.

## Version support and required packages:

MACCAS is compatible with Python 2.7, but currently untested with Python 3. 

The Python 2.7 version requires the following packages:
* numpy
* astropy (version 3.0 and above only support Python 3, so use versions 2.x)
* scipy
* argparse
* matplotlib

## Installation:

To install MACCAS you can use one of the following methods:

* Install directly from github
```
$ pip install git+https://github.com/FrancesBW/maccas.git
```
* Download from github and install via the setup file
```
$ git clone https://github.com/FrancesBW/maccas.git
$ cd maccas
$ python setup.py install
```

## Usage:

MACCAS only requires 3 inputs to run a cross match: a target catalogue, a reference catalogue and the reference catalogue format. The simplest execution on the command line would be:
```
$ maccas --target TARGET_FILENAME --reference REFERENCE_FILENAME --refsurvey REFERENCE_SURVEY_NAME
```

If the reference catalogue is from a supported survey, calling the survey name is sufficient (e.g. if using the GLEAM catalogue enter '--refsurvey GLEAM' on the command line). If you are using a reference catalogue from a non-supported survey, or of your own design, you can edit the 'reference_catalogue_format.txt' file to reflect the column names in your reference catalogue for necessary data used by MACCAS. Once you have edited the file you can call MACCAS on the command line with '--refsurvey reference_catalogue_format.txt' and MACCAS will read this file.

Currently supported catalogues are:
* GLEAM (see http://www.mwatelescope.org/gleam for instructions on downloading GLEAM catalogue)
* TGSS (see http://tgssadr.strw.leidenuniv.nl/doku.php for instructions on downloading TGSS catalogue)

### Required Inputs:

* ```--target TARGET_FILENAME```

* ```--reference REFERENCE_FILENAME```

* ```--refsurvey REFERENCE_SURVEY_NAME```

### Settings:

* ```--norestrict``` will turn off the signal-to-noise restriction that is applied on the first iteration of cross-matching. Only target sources above an SNR cutoff (where SNR is peak_flux/local_rms) are passed to the matching function. This restriction is applied by default to provide the most accurate initial offset model. 

* ```--overwrite``` will allow output files to overwrite files of the same name. If those files already exist and overwrite is active, the logger will also display a warning on which files will be overwritten. 

* ```--nofluxmodel``` will turn off the modelling of a flux adjustment factor across the target catalogue field of view. 

* ```--plot``` will turn on plotting of measured and modelled offsets of sources matched between the target and reference catalogues. It also turns on plotting for the flux model applied to the target catalogue.

* ```--nofluxmatch``` will stop peak flux of target and reference sources being compared for matching. If flux match is turned off, MACCAS returns a nearest neighbour match with modelled adjustments.

* ```--debug``` will turn on debug mode.

* ```--writelog``` will write out the log shown in the terminal to 'maccas_log_(date)\_(time).log'. 

### Optional Inputs:

* ```--tarfreq FREQUENCY ``` provide the frequency of the target catalogue. This is only required if both flux-matching is on (as is the default), and the reference catalogue has flux measurements for more than one frequency. If the reference catalogue has flux measurements for only one fequency or has no specified frequency, then MACCAS proceeds by assuming it is comparable to the target catalogue. 

* ```--outformat FILE_EXTENSION``` provide the desired file extension/format for the output tables. See astropy.table.Table docuemntation for supported writing formats (http://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers). Default format is 'fits'.

* ```--singlepercentile PERCENTILE``` provide the desired percentile cutoff for target sources with single reference match candidates within the search radius. Only target sources that meet this confidence level will be accepted as matches and used in the offset model. Default value is 0.95.

* ```--multipercentile PERCENTILE``` provide the desired normalised percentile cutoff for target sources with multiple reference match candidates within the search radius. Only target sources that meet this confidence level will be accepted as matches and used in the offset model. Default value is 0.60.

* ```--snr SNR_CUTOFF``` provide the desired signal-to-noise (SNR, defined as peak_flux/local_rms) cutoff for target sources matched in the first iteration. To remove this cutoff use '--norestrict'. Default value is 10.

* ```--modeldeg DEGREE``` provide the desired degree of model used in flux model (integer between 1 and 5). Modelled using scipy.interpolate.LSQBivariateSpline (see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.LSQBivariateSpline.html). Default value is 1.

* ```--tarformat FORMAT_FILENAME``` provide the format file for the target catalogue (see Formats section for guide). Default is detailed in Formats section.

* ```--filesuffix SUFFIX``` provide desired suffix to be appended for all output files (except log). Note a '\_' will be added between default filename and suffix, so there is no need to include this in the suffix. Default is None.


## Formats:
| foo | bar |
| --- | --- |
| baz | bim |
