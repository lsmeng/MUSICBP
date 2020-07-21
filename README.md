# MUSICBP
Beamforming and MUSIC back-projection (BP) of teleseismic P waves implemented in Matlab.

The objective of this package is to perform the Back-Projection Imaging on the seismograms of large earthquakes recorded by large-scale dense arrays. The instruction is performed on the 2018 Mw 7.5 Palu Earthquake, and you can also practise Back-Projection Imaging on the 2011 Mw 9.0 Tohoku Earthquake on your own. 

Back-Projection is an earthquake-rupture imaging technique utilizing the coherent teleseismic P wavefield based on seismic array processing. Back-tracking of seismic waves recorded by dense arrays allows Back-Projection to determine the spatio-temporal properties of the rupture (length, direction, speed, and segmentation). Over recent decades, the development of large-scale dense seismic networks has enabled the Back-Projection imaging of the rupture process of major large earthquakes. 

This code package include two versions of Back-Projection: Beamforming and Multiple Signal Classification (MUSIC). Beamforming stacks the seismograms directly in the time domain, while MUSIC is performed in the frequency domain based on the orthogonality between the noise and signal subspace of the data covariance matrix. Compared with Beamforming, MUSIC has the advantage of detecting multiple closeby sources simutaneously. More details about Back-Projection and MUSIC could be found in the following papers:

Meng, L., A. Inbal, and J.-P. Ampuero. 2011. “A window into the complexity of the dynamic rupture of the 2011 Mw 9 Tohoku-Oki earthquake”, Geophys. Res. Lett., 38, L00G07, doi:10.1029/2011GL048118.

Bao, H., Ampuero, J. P., Meng, L., Fielding, E. J., Liang, C., Milliner, C. W., ... & Huang, H. (2019). Early and persistent supershear rupture of the 2018 magnitude 7.5 Palu earthquake. Nature Geoscience, 12(3), 200-205.

Kiser, E., & Ishii, M. (2017). Back-projection imaging of earthquakes. Annual Review of Earth and Planetary Sciences, 45, 271-299.

Backprojection imaging is also performed routinely by IRIS for all new large earthquakes: https://ds.iris.edu/ds/products/backprojection/

The MUSICBP code is contributed and maintained by Han Bao (hbrandon@ucla.edu), Tian Feng (tianfengseis@gmail.com) and Lingsen Meng (meng@epss.ucla.edu). 

- [x] Instruction: MUSICBP.pdf

- [x] Code and Data: MUSICBP.zip

- [x] Related Paper: Bao_NatGeo2019.pdf and Meng_et_al-2011-Geophysical_Research_Letters.pdf
