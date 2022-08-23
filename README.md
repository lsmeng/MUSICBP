# Slowness-Enhanced Back-Projection
Slowness-Enhanced MUSIC back-projection (SEBP) of teleseismic P waves implemented in Matlab. The Matlab Signal processing, Mapping and Imaging tool boxes are required. 

The objective of this package is to perform the Back-Projection Imaging on the seismograms of large earthquakes recorded by large-scale dense arrays and calibrate the spatial errors by aftershock. The instruction is performed on the 2021 Mw 7.4 Maduo Earthquake, and you can also practise Back-Projection Imaging on the 2011 Mw 9.0 Tohoku Earthquake on your own. 

Back-Projection is an earthquake-rupture imaging technique utilizing the coherent teleseismic P wavefield based on seismic array processing. Back-tracking of seismic waves recorded by dense arrays allows Back-Projection to determine the spatio-temporal properties of the rupture (length, direction, speed, and segmentation). Over recent decades, the development of large-scale dense seismic networks has enabled the Back-Projection imaging of the rupture process of major large earthquakes. 

This code package is developed based on MUSIC BP. Please see the "MUSICBP" branch for more details about MUSIC BP.

Meng, L., A. Inbal, and J.-P. Ampuero. 2011. “A window into the complexity of the dynamic rupture of the 2011 Mw 9 Tohoku-Oki earthquake”, Geophys. Res. Lett., 38, L00G07, doi:10.1029/2011GL048118.

Bao, H., Ampuero, J. P., Meng, L., Fielding, E. J., Liang, C., Milliner, C. W., ... & Huang, H. (2019). Early and persistent supershear rupture of the 2018 magnitude 7.5 Palu earthquake. Nature Geoscience, 12(3), 200-205.

Kiser, E., & Ishii, M. (2017). Back-projection imaging of earthquakes. Annual Review of Earth and Planetary Sciences, 45, 271-299.

Backprojection imaging is also performed routinely by IRIS for all new large earthquakes: https://ds.iris.edu/ds/products/backprojection/

The MUSICBP code is contributed and maintained by Han Bao (hbrandon@ucla.edu), Tian Feng (tianfengseis@gmail.com) and Lingsen Meng (meng@epss.ucla.edu). 

- [x] Instruction: SEBP.pdf

- [x] Code and Data: AU*.m

- [x] Related Paper: BaoNG_2019.pdf
