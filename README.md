# Mach wave analysis
Mach wave analysis in Python. The Obspy modulus are required. 

The objective of this script is to identify the Mach wave for supershear earthquakes and measure the cross-correlation coefficient between mainshock Rayleigh wave and EGF Rayleight wave. 

Mach waves, a unique signature of supershear ruptures, can be effectively identified in the far-field surface wavefield. For earthquakes rupturing at velocities slower than Rayleigh wave speed, wavefronts from different sub-sources arrive at any far-field receiver at different times. For supershear earthquakes, only at stations located on the Rayleigh Mach cone do waves from different parts of a supershear rupture arrive simultaneously and interfere constructively, resulting in an apparent point-source. As first proposed by Vallée and Dunham, along the Mach cone, the waveform of a supershear mainshock should be identical to that of a smaller reference event, referred to as an empirical Green’s function (EGF) event, with a similar focal mechanism in the vicinity of the mainshock. The Mach wave identification is conducted by evaluating the waveform similarity between a target event and its collocated EGF.

The code is contributed and maintained by Liuwei Xu (xuliuw1997@ucla.edu), Han Bao (hbrandon@ucla.edu), and Lingsen Meng (meng@epss.ucla.edu). 

- [x] Code and Data: MUSICBP.zip

- [x] Related Paper: Kokoxili 2001.pdf by M. Vallée and Eric M. Dunham, 2012
