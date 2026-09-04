# get_metadata

Returns a spatial grid (as an sf object) containing metadata and
download URLs for global 3D building footprint tiles (3D-GloBFP).

## Usage

``` r
get_metadata(test = FALSE, quiet = TRUE)
```

## Arguments

- test:

  logic, Ignored during normal use; included for internal testing
  purposes. Defaults to `FALSE`.

- quiet:

  logical. Accepted for consistency with other package functions;
  currently unused because this function does not emit cli messages.

## Value

sf a spatial polygon grid with attributes: `id`, `gridID`, bounding box
coordinates, and `download_url`.

## Details

The metadata of 3D Global Building Footprints (3D-GloBFP) dataset is
uploaded on zenodo. More detials about this dataset can to found
[here](https://zenodo.org/records/15487037).

The data is detailed in the following article

## References

Che, Y., Li, X., Liu, X., Wang, Y., Liao, W., Zheng, X., Zhang, X., Xu,
X., Shi, Q., Zhu, J., Zhang, H., Yuan, H., & Dai, Y. (2025). 3D-GloBFP:
the first global three-dimensional building footprint dataset. Zenodo.
https://doi.org/10.5281/zenodo.15487037

Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng
Xianwei, Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui,
Yuan Hua, & Dai Yongjiu (2024). 3D-GloBFP: the first global
three-dimensional building footprint dataset. Earth Syst. Sci. Data, 16,
5357-5374

## Examples

``` r
meta <- gloBFPr::get_metadata(test=TRUE)
```
