MACCAS: Multi Attribute Cross-matcher with Correction for wArped Sky
======

MACCAS is a cross-matching algorithm that aims to match all sources in the supplied target catalogue to sources in a reference catalogue. MACCAS utilises not only the position of a source on the sky, but also the flux data, to determine the most probable match in the reference catalog to the target source. This two variable approach means MACCAS is a more rigorous cross-matcher than the matching tools supplied by TOPCAT or astropy alone. Additionally, MACCAS uses an iterative matching approach, where a model of the offsets of matched sources is applied to unmatched sources to correct for any image distortion 

MACCAS also provides the option for a simple flux correction across the target catalogue
