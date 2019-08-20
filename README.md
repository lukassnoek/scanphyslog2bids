# PhilipsPhysio2BIDS
Code to convert Philips physiology files ("SCANPHYSLOG") to the BIDS-format, including volume triggers.
This is very ugly code, but it seems to work. It writes out BIDSified physio-files (as \*.tsv.gz and \*.json files).
I recommend using the [PhysIO toolbox](https://github.com/translationalneuromodeling/tapas/tree/master/PhysIO) to convert the BIDSified files to RETROICOR/HRV/RVT regressors
(worked really well for me in the past).

Feel free to submit a PR with proposed changes/enhancements/fixes.

