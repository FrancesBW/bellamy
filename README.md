MACCAS: Multi Attribute Cross-matcher with Correction for wArped Sky
======

MACCAS is a cross-matching algorithm that aims to match all sources in the supplied target catalogue to sources in a reference catalogue. MACCAS utilises not only the position of a source on the sky, but also the flux data, to determine the most probable match in the reference catalog to the target source. This two variable approach means MACCAS is a more rigorous cross-matcher than the matching tools supplied by TOPCAT or astropy alone. Additionally, MACCAS attempts to undo any spatial distortion that may be affecting the target catalogue, by creating a model of the offsets of matched sources which is then applied to unmatched sources. Furthermore, MACCAS provides the option for a simple flux correction across the target catalogue, in the case that there is some large scale calibration error affecting the target catalogue.

Version support and required packages:
======

MACCAS is compatible with Python 2.7, but currently untested with Python 3. 

The Python 2.7 version requires the following packages:
* numpy
* astropy (version 3.0 and above only support Python 3, so use versions 2.x)
* scipy
* argparse
* matplotlib

Installation:
======

To install MACCAS you can use one of the following methods:

* Install directly from github
```
git install https://github.com/FrancesBW/maccas.git
```
* Download from github and install via the setup file
```
git clone https://github.com/FrancesBW/maccas.git
cd maccas
python setup.py install
```

