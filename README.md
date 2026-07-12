# MUSICBP

**Beamforming and MUSIC Back-Projection (BP) of teleseismic P waves for earthquake rupture imaging, implemented in MATLAB.**

Back-Projection is an earthquake-rupture imaging technique utilizing the coherent teleseismic P wavefield based on seismic array processing. Back-tracking of seismic waves recorded by dense arrays allows Back-Projection to determine the spatio-temporal properties of the rupture (length, direction, speed, and segmentation). Over recent decades, the development of large-scale dense seismic networks has enabled the Back-Projection imaging of the rupture process of major large earthquakes.

This package includes two versions of Back-Projection:

- **Beamforming** — stacks the seismograms directly in the time domain.
- **MUSIC (Multiple Signal Classification)** — performed in the frequency domain based on the orthogonality between the noise and signal subspaces of the data covariance matrix. Compared with Beamforming, MUSIC has the advantage of resolving multiple closely spaced simultaneous sources.

The tutorial walks through the 2018 Mw 7.5 Palu earthquake, and sample data for the 2011 Mw 9.0 Tohoku earthquake are included for further practice. Back-Projection imaging is also performed routinely by IRIS for all new large earthquakes: https://ds.iris.edu/ds/products/backprojection/

## Repository layout

| Path | Content |
|------|---------|
| [`code/`](code/) | MATLAB source code (`General_BP.m` is the main driver), plate-boundary/coastline data files, and sample waveform datasets for the 2018 Palu (`PaluAUData/`) and 2011 Tohoku (`TohokuTAData/`) earthquakes |
| [`code/README.md`](code/README.md) | Step-by-step usage instructions |
| [`MUSICBP.pdf`](MUSICBP.pdf) | Full tutorial document |
| `Bao_NatGeo2019.pdf`, `Meng_et_al-2011-Geophysical_Research_Letters.pdf` | Reference papers |

## Requirements

MATLAB with the Signal Processing, Mapping, and Image Processing toolboxes.

## Quick start

1. Open [`code/General_BP.m`](code/General_BP.m), set `Initial_flag = 1` and a project name (e.g. `Palu_2018`), and run it to initialize the project folder.
2. Copy your SAC waveform files into the project's `Data/` folder (or use the included Palu/Tohoku sample data).
3. Follow the step-by-step flags (`readBP_flag` → `alignBP_flag` → `runBPbmfm_flag` / `runBPmusic_flag`) described in [`code/README.md`](code/README.md) and `MUSICBP.pdf` to read, align, and back-project the data. Outputs include BP movies (`movie.gif`), summary plots (`summary.pdf`), and radiator coordinates (`HFdots`).

## Citation

If you use this code, please cite:

- Meng, L., A. Inbal, and J.-P. Ampuero (2011). A window into the complexity of the dynamic rupture of the 2011 Mw 9 Tohoku-Oki earthquake. *Geophys. Res. Lett.*, 38, L00G07. doi:10.1029/2011GL048118
- Bao, H., Ampuero, J.-P., Meng, L., Fielding, E. J., Liang, C., Milliner, C. W., Feng, T., & Huang, H. (2019). Early and persistent supershear rupture of the 2018 magnitude 7.5 Palu earthquake. *Nature Geoscience*, 12(3), 200–205.

For a review of the method: Kiser, E., & Ishii, M. (2017). Back-projection imaging of earthquakes. *Annual Review of Earth and Planetary Sciences*, 45, 271–299.

## Selected publications using MUSICBP

The MUSIC back-projection method implemented in this package has underpinned rupture-imaging studies of many of the most significant recent earthquakes, including:

- Xu, L., Meng, L., Yunjun, Z., et al. (2025). Ultralong, supershear rupture of the 2025 Mw 7.7 Mandalay earthquake reveals unaccounted risk. *Science*.
- Xu, L., Ji, C., Meng, L., Ampuero, J.-P., et al. (2024). Dual-initiation ruptures in the 2024 Noto earthquake encircling a fault asperity at a swarm edge. *Science*, 385(6711), 871–876. (cover article)
- Bao, H., Xu, L., Meng, L., Ampuero, J.-P., et al. (2022). Global frequency of oceanic and continental supershear earthquakes. *Nature Geoscience*, 15, 942–949.
- Bao, H., Ampuero, J.-P., Meng, L., et al. (2019). Early and persistent supershear rupture of the 2018 magnitude 7.5 Palu earthquake. *Nature Geoscience*, 12(3), 200–205.
- Avouac, J.-P., Meng, L., Wei, S., Wang, T., & Ampuero, J.-P. (2015). Lower edge of locked Main Himalayan Thrust unzipped by the 2015 Gorkha earthquake. *Nature Geoscience*, 8, 708–711.
- Simons, M., Minson, S. E., Sladen, A., Ortega, F., Jiang, J., Owen, S. E., Meng, L., Ampuero, J.-P., et al. (2011). The 2011 magnitude 9.0 Tohoku-Oki earthquake: Mosaicking the megathrust from seconds to centuries. *Science*, 332(6036), 1421–1425. doi:10.1126/science.1206731

## Maintainers

The MUSICBP code is contributed and maintained by Han Bao (hbrandon@ucla.edu), Tian Feng (tianfengseis@gmail.com), and Lingsen Meng (meng@epss.ucla.edu), UCLA Earth, Planetary, and Space Sciences.

## License

MIT — see [LICENSE](LICENSE).
