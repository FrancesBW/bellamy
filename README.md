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
$ git install https://github.com/FrancesBW/maccas.git
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
$ maccas --target target_filename --reference reference_filename --refsurvey reference_survey_name
```

If the reference catalogue is from a supported survey, calling the survey name is sufficient (e.g. if using the GLEAM catalogue enter '--refsurvey GLEAM' on the command line). If you are using a reference catalogue from a non-supported survey, or of your own design, you can edit the 'reference_catalogue_format.txt' file to reflect the column names in your reference catalogue for necessary data used by MACCAS. Once you have edited the file you can call MACCAS on the command line with '--refsurvey reference_catalogue_format.txt' and MACCAS will read this file.

Currently supported catalogues are:
* GLEAM (see http://www.mwatelescope.org/gleam for instructions on downloading GLEAM catalogue)
* TGSS (see http://tgssadr.strw.leidenuniv.nl/doku.php for instructions on downloading TGSS catalogue)

### Settings:

* ```--norestrict``` will turn off the signal-to-noise restriction that is applied on the first iteration of cross-matching. Only target sources above an SNR cutoff (where SNR is peak_flux/local_rms) are passed to the matching function. This restriction is applied by default to provide the most accurate initial offset model. 



