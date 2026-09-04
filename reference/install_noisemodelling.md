# Install the headless NoiseModelling runner

Downloads and unzips the official headless NoiseModelling release into
an R user cache directory. This is used by
[`get_noise_map()`](https://billbillbilly.github.io/gloBFPr/reference/get_noise_map.md)
when `run = TRUE` and no `nm_path` is supplied.

## Usage

``` r
install_noisemodelling(
  version = "5.0.1",
  destdir = noisemodelling_cache_dir(),
  quiet = TRUE
)
```

## Arguments

- version:

  NoiseModelling release version. Defaults to `"5.0.1"`.

- destdir:

  Destination directory. Defaults to the `gloBFPr` user cache.

- quiet:

  Logical. If `TRUE`, suppress download messages.

## Value

Path to the installed NoiseModelling directory.
