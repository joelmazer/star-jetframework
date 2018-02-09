# corrections
*this will be updated...*
This directory contains correction files and header for:
* event plane re-centering
* event plane shifting

## methods will include use of:
* BBC detector
* ZDC detector
* TPC detector

### TPC methods
When using the TPC, various methods will be employed.  When examining the effect of the leading jet, checks will be made by
* removing the jet constituents from the event plane calculation in the eta-phi cone around the jet axis
* removing ALL tracks from the event plane calculation in the eta-phi cone around the jet axis
* removing the eta strip around the axis of the jet

pt cuts will be made to the tracks used in reconstruction of
* 0.2 < pt < 2.0 GeV/c
* 0.2 < pt < 5.0 GeV/c

The weight of such tracks is the pt of the corresponding track.  One can use a linear weight up to max pt, or linear up to 2.0 GeV and a constant of 2 when tracks are > 2.0 GeV.

## Author
**Joel Mazer**

## Notes
up-to-date file usage
for BBC and ZDC recentering corrections:
* Method1: recenter_calib_file_bin0_STEP1_Feb1.root
* Std: recenter_calib_file.root

for BBC and ZDC shifting corrections:
* Method1: shift_calib_file_bin0_STEP2_Feb5.root
* Std: shift_calib_file.root

