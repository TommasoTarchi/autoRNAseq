# Third-party software in the container images

The container images attached to the autoRNAseq releases contain software written by other authors. This software is **not** covered by the MIT License of this repository: each component is distributed under its own license, listed below.

The pipeline code in this repository only calls these tools as separate programs and does not include any of their code.

## Contents

- [Main tools](#main-tools)
- [Usage restriction: rMATS-turbo](#usage-restriction-rmats-turbo)
- [Source code](#source-code)
- [Other components of the images](#other-components-of-the-images)
- [License texts](#license-texts)

## Main tools

| Image file | Software | Version | License | Source code |
|---|---|---|---|---|
| `STAR-v2.7.11b.sif` | STAR | 2.7.10b | MIT | https://github.com/alexdobin/STAR/tree/2.7.10b |
| `trim_galore-v0.6.7.sif` | Trim Galore! | 0.6.7 | GPL-3.0 | https://github.com/FelixKrueger/TrimGalore/tree/0.6.7 |
| `trim_galore-v0.6.7.sif` | Cutadapt | 3.4 | MIT | https://github.com/marcelm/cutadapt/tree/v3.4 |
| `trim_galore-v0.6.7.sif` | FastQC | 0.11.9 | GPL-3.0-or-later (bundled Java libraries under their own licenses, see note) | https://github.com/s-andrews/FastQC/tree/v0.11.9 |
| `SAMtools-v1.3.1.sif` | SAMtools | 1.11 | MIT/Expat | https://github.com/samtools/samtools/tree/1.11 |
| `SAMtools-v1.3.1.sif` | HTSlib | 1.11 | MIT/Expat; Modified BSD for files in `cram/` | https://github.com/samtools/htslib/tree/1.11 |
| `picard-v3.1.1.sif` | Picard (picard-slim) | 2.27.4 | MIT | https://github.com/broadinstitute/picard/tree/2.27.4 |
| `featureCounts-v2.0.6.sif` | Subread (featureCounts) | 2.0.1 | GPL-3.0 | https://sourceforge.net/projects/subread/files/subread-2.0.1/ (subread-2.0.1-source.tar.gz) |
| `HTSeq-v2.0.2.sif` | HTSeq | 2.0.2 | GPL-3.0 | https://github.com/htseq/htseq/tree/release_2.0.2 |
| `rMATS-v4.3.0.sif` | rMATS-turbo | 4.3.0 | BSD-2-Clause, **academic and non-profit research use only** | https://github.com/Xinglab/rmats-turbo/tree/v4.3.0 |
| `rMATS-v4.3.0.sif` | STAR | 2.7.10b | MIT | https://github.com/alexdobin/STAR/tree/2.7.10b |
| `rMATS-v4.3.0.sif` | GNU Scientific Library (GSL) | 2.7.1 | GPL-3.0-or-later | https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz |
| `multiQC-v1.18.sif` | MultiQC | 1.18 | GPL-3.0 | https://github.com/MultiQC/MultiQC/tree/v1.18 |

The version column gives the version actually installed in each image.

**Note on FastQC:** FastQC includes third-party Java libraries (for example JHDF5, under the Apache License 2.0). Their license files are included in the image, in `/usr/local/opt/fastqc-0.11.9/`.

## Usage restriction: rMATS-turbo

rMATS-turbo, used in the splicing analysis step (image `rMATS-v4.3.0.sif`), may be used by academic and non-profit institutions for research purposes. For use by or transfer to commercial entities, contact Dr. Yi Xing (xingyi@chop.edu). The full license is reproduced [below](#rmats-turbo).

## Source code

The source code of each main tool, at the exact version included in the images, is available at the links in the [Main tools](#main-tools) table.

The other components of the images were installed as conda packages. The source code of each package can be obtained from the recipe that built it:

- packages from `conda-forge`: `https://github.com/conda-forge/<package>-feedstock`
- packages from `bioconda`: `https://github.com/bioconda/bioconda-recipes/tree/master/recipes/<package>`
- packages from `anaconda`: `https://github.com/AnacondaRecipes/<package>-feedstock`

Each recipe points to the upstream source archive for the packaged version. The source code of the Debian packages in the base systems can be obtained from the Debian archive (https://snapshot.debian.org).

## Other components of the images

Besides the main tools, each image contains a Linux base system built from Debian packages, and supporting packages installed with conda (Python, Perl, a Java runtime, compression libraries, and others).

- **Debian packages:** the copyright and license of each Debian package included are in `/usr/share/doc/<package>/copyright` inside the image.
- **Conda packages:** listed below for each image, with the license recorded in the package metadata. Where a package provides license files, they are inside the image under `<conda prefix>/pkgs/<package>/info/licenses/` or `<conda prefix>/conda-meta/`.

<details>
<summary><code>STAR-v2.7.11b.sif</code>: 36 conda packages, base image continuumio/miniconda3:4.8.3 (Debian 10)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | anaconda/main |
| ca-certificates | 2024.3.11 | MPL-2.0 | anaconda/main |
| certifi | 2024.2.2 | MPL-2.0 | anaconda/main |
| cffi | 1.14.0 | MIT | anaconda/main |
| chardet | 3.0.4 | LGPL2 | anaconda/main |
| conda | 4.8.3 | BSD 3-Clause | anaconda/main |
| conda-package-handling | 1.6.1 | BSD-3-Clause | anaconda/main |
| cryptography | 2.9.2 | Apache-2.0 AND BSD-3-Clause AND PSF-2.0 | anaconda/main |
| idna | 2.9 | BSD Like | anaconda/main |
| ld_impl_linux-64 | 2.33.1 | GPL 3 | anaconda/main |
| libedit | 3.1.20181209 | NetBSD | anaconda/main |
| libffi | 3.3 | Custom | anaconda/main |
| libgcc-ng | 9.1.0 | GPL | anaconda/main |
| libstdcxx-ng | 9.1.0 | GPL3 with runtime exception | anaconda/main |
| ncurses | 6.2 | Free software (MIT-like) | anaconda/main |
| openssl | 1.1.1w | OpenSSL | anaconda/main |
| pip | 20.0.2 | MIT | anaconda/main |
| pycosat | 0.6.3 | MIT | anaconda/main |
| pycparser | 2.20 | BSD 3-clause | anaconda/main |
| pyopenssl | 19.1.0 | Apache 2.0 | anaconda/main |
| pysocks | 1.7.1 | BSD 3-Clause | anaconda/main |
| python | 3.8.3 | PSF | anaconda/main |
| readline | 8.0 | GPL-3.0 | anaconda/main |
| requests | 2.23.0 | Apache-2.0 | anaconda/main |
| ruamel_yaml | 0.15.87 | MIT | anaconda/main |
| setuptools | 46.4.0 | MIT | anaconda/main |
| six | 1.14.0 | MIT | anaconda/main |
| sqlite | 3.31.1 | Public-Domain (http://www.sqlite.org/copyright.html) | anaconda/main |
| star | 2.7.10b | GPLv3 | bioconda |
| tk | 8.6.8 | Tcl/Tk | anaconda/main |
| tqdm | 4.46.0 | MPLv2.0, MIT | anaconda/main |
| urllib3 | 1.25.8 | MIT | anaconda/main |
| wheel | 0.34.2 | MIT | anaconda/main |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | anaconda/main |
| yaml | 0.1.7 | MIT | anaconda/main |
| zlib | 1.2.11 | zlib | anaconda/main |

</details>

<details>
<summary><code>trim_galore-v0.6.7.sif</code>: 77 conda packages, base image trim-galore:0.6.7--hdfd78af_0 (BioContainers, minimal Debian 10 base)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | conda-forge |
| _openmp_mutex | 4.5 | BSD-3-Clause | conda-forge |
| alsa-lib | 1.2.3 | LGPL-2.1-or-later | conda-forge |
| ca-certificates | 2021.5.30 | ISC | conda-forge |
| cairo | 1.16.0 | LGPL 2.1 or MPL 1.1 | conda-forge |
| certifi | 2021.5.30 | ISC | conda-forge |
| cutadapt | 3.4 | MIT | bioconda |
| dnaio | 0.5.1 | MIT | bioconda |
| fastqc | 0.11.9 | GPL >=3 | bioconda |
| font-ttf-dejavu-sans-mono | 2.37 | BSD-3-Clause | conda-forge |
| fontconfig | 2.13.1 | MIT | conda-forge |
| freetype | 2.10.4 | GPL-2.0-only and LicenseRef-FreeType | conda-forge |
| gettext | 0.19.8.1 | LGPL-2.1-or-later AND GPL-3.0-or-later | conda-forge |
| giflib | 5.2.1 | MIT | conda-forge |
| graphite2 | 1.3.13 | LGPLv2 | conda-forge |
| harfbuzz | 2.8.2 | MIT | conda-forge |
| icu | 68.1 | MIT | conda-forge |
| isa-l | 2.30.0 | BSD-3-Clause | conda-forge |
| jbig | 2.1 | GPL-2.0-or-later | conda-forge |
| jpeg | 9d | IJG | conda-forge |
| lcms2 | 2.12 | MIT | conda-forge |
| ld_impl_linux-64 | 2.36.1 | GPL-3.0-only | conda-forge |
| lerc | 2.2.1 | Apache-2.0 | conda-forge |
| libdeflate | 1.7 | MIT | conda-forge |
| libffi | 3.3 | MIT | conda-forge |
| libgcc-ng | 11.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libglib | 2.68.3 | LGPL-2.1-or-later | conda-forge |
| libgomp | 11.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libiconv | 1.16 | GPL and LGPL | conda-forge |
| libpng | 1.6.37 | zlib/libpng | conda-forge |
| libstdcxx-ng | 11.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libtiff | 4.3.0 | HPND | conda-forge |
| libuuid | 2.32.1 | BSD-3-Clause | conda-forge |
| libwebp-base | 1.2.0 | BSD-3-Clause | conda-forge |
| libxcb | 1.13 | MIT | conda-forge |
| libxml2 | 2.9.12 | MIT | conda-forge |
| lz4-c | 1.9.3 | BSD-2-Clause | conda-forge |
| ncurses | 6.2 | X11 AND BSD-3-Clause | conda-forge |
| openjdk | 11.0.9.1 | GPL-2.0-or-later WITH Classpath-exception-2.0 | conda-forge |
| openssl | 1.1.1k | OpenSSL | conda-forge |
| pcre | 8.45 | BSD-3-Clause | conda-forge |
| perl | 5.32.1 | GPL-1.0-or-later OR Artistic-1.0-Perl | conda-forge |
| pigz | 2.6 | zlib | conda-forge |
| pip | 21.1.3 | MIT | conda-forge |
| pixman | 0.40.0 | MIT | conda-forge |
| pthread-stubs | 0.4 | MIT | conda-forge |
| python | 3.9.6 | Python-2.0 | conda-forge |
| python-isal | 0.10.0 | MIT | conda-forge |
| python_abi | 3.9 | BSD-3-Clause | conda-forge |
| readline | 8.1 | GPL-3.0-only | conda-forge |
| setuptools | 49.6.0 | MIT | conda-forge |
| sqlite | 3.36.0 | Unlicense | conda-forge |
| tk | 8.6.10 | TCL | conda-forge |
| trim-galore | 0.6.7 | GPL | bioconda (local build) |
| tzdata | 2021a | LicenseRef-Public-Domain | conda-forge |
| wheel | 0.36.2 | MIT | conda-forge |
| xopen | 1.1.0 | MIT | conda-forge |
| xorg-fixesproto | 5.0 | MIT | conda-forge |
| xorg-inputproto | 2.3.2 | MIT | conda-forge |
| xorg-kbproto | 1.0.7 | MIT | conda-forge |
| xorg-libice | 1.0.10 | MIT | conda-forge |
| xorg-libsm | 1.2.3 | MIT | conda-forge |
| xorg-libx11 | 1.7.2 | MIT | conda-forge |
| xorg-libxau | 1.0.9 | MIT | conda-forge |
| xorg-libxdmcp | 1.1.3 | MIT | conda-forge |
| xorg-libxext | 1.3.4 | MIT | conda-forge |
| xorg-libxfixes | 5.0.3 | MIT | conda-forge |
| xorg-libxi | 1.7.10 | MIT | conda-forge |
| xorg-libxrender | 0.9.10 | MIT | conda-forge |
| xorg-libxtst | 1.2.3 | MIT | conda-forge |
| xorg-recordproto | 1.14.2 | MIT | conda-forge |
| xorg-renderproto | 0.11.1 | MIT | conda-forge |
| xorg-xextproto | 7.3.0 | MIT | conda-forge |
| xorg-xproto | 7.0.31 | MIT | conda-forge |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | conda-forge |
| zlib | 1.2.11 | Zlib | conda-forge |
| zstd | 1.5.0 | BSD-3-Clause | conda-forge |

</details>

<details>
<summary><code>SAMtools-v1.3.1.sif</code>: 42 conda packages, base image continuumio/miniconda3:4.8.3 (Debian 10)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | anaconda/main |
| bzip2 | 1.0.8 | bzip2 | anaconda/main |
| ca-certificates | 2024.3.11 | MPL-2.0 | anaconda/main |
| certifi | 2024.6.2 | MPL-2.0 | anaconda/main |
| cffi | 1.14.0 | MIT | anaconda/main |
| chardet | 3.0.4 | LGPL2 | anaconda/main |
| conda | 4.8.3 | BSD 3-Clause | anaconda/main |
| conda-package-handling | 1.6.1 | BSD-3-Clause | anaconda/main |
| cryptography | 2.9.2 | Apache-2.0 AND BSD-3-Clause AND PSF-2.0 | anaconda/main |
| htslib | 1.11 | MIT | bioconda |
| idna | 2.9 | BSD Like | anaconda/main |
| krb5 | 1.17.1 | MIT | anaconda/main |
| ld_impl_linux-64 | 2.33.1 | GPL 3 | anaconda/main |
| libcurl | 7.71.1 | MIT/X derivate (http://curl.haxx.se/docs/copyright.html) | conda-forge |
| libdeflate | 1.8 | MIT | anaconda/main |
| libedit | 3.1.20181209 | NetBSD | anaconda/main |
| libffi | 3.3 | Custom | anaconda/main |
| libgcc-ng | 9.1.0 | GPL | anaconda/main |
| libssh2 | 1.10.0 | BSD-3-Clause | anaconda/main |
| libstdcxx-ng | 9.1.0 | GPL3 with runtime exception | anaconda/main |
| ncurses | 6.2 | Free software (MIT-like) | anaconda/main |
| openssl | 1.1.1w | OpenSSL | anaconda/main |
| pip | 20.0.2 | MIT | anaconda/main |
| pycosat | 0.6.3 | MIT | anaconda/main |
| pycparser | 2.20 | BSD 3-clause | anaconda/main |
| pyopenssl | 19.1.0 | Apache 2.0 | anaconda/main |
| pysocks | 1.7.1 | BSD 3-Clause | anaconda/main |
| python | 3.8.3 | PSF | anaconda/main |
| readline | 8.0 | GPL-3.0 | anaconda/main |
| requests | 2.23.0 | Apache-2.0 | anaconda/main |
| ruamel_yaml | 0.15.87 | MIT | anaconda/main |
| samtools | 1.11 | MIT | bioconda |
| setuptools | 46.4.0 | MIT | anaconda/main |
| six | 1.14.0 | MIT | anaconda/main |
| sqlite | 3.31.1 | Public-Domain (http://www.sqlite.org/copyright.html) | anaconda/main |
| tk | 8.6.8 | Tcl/Tk | anaconda/main |
| tqdm | 4.46.0 | MPLv2.0, MIT | anaconda/main |
| urllib3 | 1.25.8 | MIT | anaconda/main |
| wheel | 0.34.2 | MIT | anaconda/main |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | anaconda/main |
| yaml | 0.1.7 | MIT | anaconda/main |
| zlib | 1.2.11 | zlib | anaconda/main |

</details>

<details>
<summary><code>picard-v3.1.1.sif</code>: 37 conda packages, base image continuumio/miniconda3:4.8.3 (Debian 10)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | anaconda/main |
| ca-certificates | 2024.3.11 | MPL-2.0 | anaconda/main |
| certifi | 2024.6.2 | MPL-2.0 | anaconda/main |
| cffi | 1.14.0 | MIT | anaconda/main |
| chardet | 3.0.4 | LGPL2 | anaconda/main |
| conda | 4.8.3 | BSD 3-Clause | anaconda/main |
| conda-package-handling | 1.6.1 | BSD-3-Clause | anaconda/main |
| cryptography | 2.9.2 | Apache-2.0 AND BSD-3-Clause AND PSF-2.0 | anaconda/main |
| idna | 2.9 | BSD Like | anaconda/main |
| ld_impl_linux-64 | 2.33.1 | GPL 3 | anaconda/main |
| libedit | 3.1.20181209 | NetBSD | anaconda/main |
| libffi | 3.3 | Custom | anaconda/main |
| libgcc-ng | 9.1.0 | GPL | anaconda/main |
| libstdcxx-ng | 9.1.0 | GPL3 with runtime exception | anaconda/main |
| ncurses | 6.2 | Free software (MIT-like) | anaconda/main |
| openjdk | 8.0.152 | GPLv2 | anaconda/main |
| openssl | 1.1.1w | OpenSSL | anaconda/main |
| picard-slim | 2.27.4 | MIT | bioconda |
| pip | 20.0.2 | MIT | anaconda/main |
| pycosat | 0.6.3 | MIT | anaconda/main |
| pycparser | 2.20 | BSD 3-clause | anaconda/main |
| pyopenssl | 19.1.0 | Apache 2.0 | anaconda/main |
| pysocks | 1.7.1 | BSD 3-Clause | anaconda/main |
| python | 3.8.3 | PSF | anaconda/main |
| readline | 8.0 | GPL-3.0 | anaconda/main |
| requests | 2.23.0 | Apache-2.0 | anaconda/main |
| ruamel_yaml | 0.15.87 | MIT | anaconda/main |
| setuptools | 46.4.0 | MIT | anaconda/main |
| six | 1.14.0 | MIT | anaconda/main |
| sqlite | 3.31.1 | Public-Domain (http://www.sqlite.org/copyright.html) | anaconda/main |
| tk | 8.6.8 | Tcl/Tk | anaconda/main |
| tqdm | 4.46.0 | MPLv2.0, MIT | anaconda/main |
| urllib3 | 1.25.8 | MIT | anaconda/main |
| wheel | 0.34.2 | MIT | anaconda/main |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | anaconda/main |
| yaml | 0.1.7 | MIT | anaconda/main |
| zlib | 1.2.11 | zlib | anaconda/main |

</details>

<details>
<summary><code>featureCounts-v2.0.6.sif</code>: 36 conda packages, base image continuumio/miniconda3:4.8.3 (Debian 10)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | anaconda/main |
| ca-certificates | 2024.3.11 | MPL-2.0 | anaconda/main |
| certifi | 2024.2.2 | MPL-2.0 | anaconda/main |
| cffi | 1.14.0 | MIT | anaconda/main |
| chardet | 3.0.4 | LGPL2 | anaconda/main |
| conda | 4.8.3 | BSD 3-Clause | anaconda/main |
| conda-package-handling | 1.6.1 | BSD-3-Clause | anaconda/main |
| cryptography | 2.9.2 | Apache-2.0 AND BSD-3-Clause AND PSF-2.0 | anaconda/main |
| idna | 2.9 | BSD Like | anaconda/main |
| ld_impl_linux-64 | 2.33.1 | GPL 3 | anaconda/main |
| libedit | 3.1.20181209 | NetBSD | anaconda/main |
| libffi | 3.3 | Custom | anaconda/main |
| libgcc-ng | 9.1.0 | GPL | anaconda/main |
| libstdcxx-ng | 9.1.0 | GPL3 with runtime exception | anaconda/main |
| ncurses | 6.2 | Free software (MIT-like) | anaconda/main |
| openssl | 1.1.1w | OpenSSL | anaconda/main |
| pip | 20.0.2 | MIT | anaconda/main |
| pycosat | 0.6.3 | MIT | anaconda/main |
| pycparser | 2.20 | BSD 3-clause | anaconda/main |
| pyopenssl | 19.1.0 | Apache 2.0 | anaconda/main |
| pysocks | 1.7.1 | BSD 3-Clause | anaconda/main |
| python | 3.8.3 | PSF | anaconda/main |
| readline | 8.0 | GPL-3.0 | anaconda/main |
| requests | 2.23.0 | Apache-2.0 | anaconda/main |
| ruamel_yaml | 0.15.87 | MIT | anaconda/main |
| setuptools | 46.4.0 | MIT | anaconda/main |
| six | 1.14.0 | MIT | anaconda/main |
| sqlite | 3.31.1 | Public-Domain (http://www.sqlite.org/copyright.html) | anaconda/main |
| subread | 2.0.1 | GPLv3 | bioconda |
| tk | 8.6.8 | Tcl/Tk | anaconda/main |
| tqdm | 4.46.0 | MPLv2.0, MIT | anaconda/main |
| urllib3 | 1.25.8 | MIT | anaconda/main |
| wheel | 0.34.2 | MIT | anaconda/main |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | anaconda/main |
| yaml | 0.1.7 | MIT | anaconda/main |
| zlib | 1.2.11 | zlib | anaconda/main |

</details>

<details>
<summary><code>HTSeq-v2.0.2.sif</code>: 76 conda packages, base image htseq:2.0.2--py310ha14a713_0 (BioContainers, minimal Debian 10 base)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | conda-forge |
| _openmp_mutex | 4.5 | BSD-3-Clause | conda-forge |
| brotli | 1.0.9 | MIT | conda-forge |
| brotli-bin | 1.0.9 | MIT | conda-forge |
| bzip2 | 1.0.8 | bzip2-1.0.6 | conda-forge |
| c-ares | 1.18.1 | MIT | conda-forge |
| ca-certificates | 2022.6.15 | ISC | conda-forge |
| certifi | 2022.6.15 | ISC | conda-forge |
| cycler | 0.11.0 | BSD-3-Clause | conda-forge |
| fonttools | 4.34.4 | MIT | conda-forge |
| freetype | 2.10.4 | GPL-2.0-only and LicenseRef-FreeType | conda-forge |
| giflib | 5.2.1 | MIT | conda-forge |
| htseq | 2.0.2 | GPL-3.0-only | bioconda (local build) |
| jpeg | 9e | IJG | conda-forge |
| keyutils | 1.6.1 | LGPL-2.1-or-later | conda-forge |
| kiwisolver | 1.4.4 | BSD-3-Clause | conda-forge |
| krb5 | 1.19.3 | MIT | conda-forge |
| lcms2 | 2.12 | MIT | conda-forge |
| ld_impl_linux-64 | 2.36.1 | GPL-3.0-only | conda-forge |
| lerc | 3.0 | Apache-2.0 | conda-forge |
| libblas | 3.9.0 | BSD-3-Clause | conda-forge |
| libbrotlicommon | 1.0.9 | MIT | conda-forge |
| libbrotlidec | 1.0.9 | MIT | conda-forge |
| libbrotlienc | 1.0.9 | MIT | conda-forge |
| libcblas | 3.9.0 | BSD-3-Clause | conda-forge |
| libcurl | 7.83.1 | curl | conda-forge |
| libdeflate | 1.10 | MIT | conda-forge |
| libedit | 3.1.20191231 | BSD-2-Clause | conda-forge |
| libev | 4.33 | BSD-2-Clause | conda-forge |
| libffi | 3.4.2 | MIT | conda-forge |
| libgcc-ng | 12.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgfortran-ng | 12.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgfortran5 | 12.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgomp | 12.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| liblapack | 3.9.0 | BSD-3-Clause | conda-forge |
| libnghttp2 | 1.47.0 | MIT | conda-forge |
| libnsl | 2.0.0 | GPL-2.0-only | conda-forge |
| libopenblas | 0.3.20 | BSD-3-Clause | conda-forge |
| libpng | 1.6.37 | zlib/libpng | conda-forge |
| libssh2 | 1.10.0 | BSD-3-Clause | conda-forge |
| libstdcxx-ng | 12.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libtiff | 4.4.0 | HPND | conda-forge |
| libuuid | 2.32.1 | BSD-3-Clause | conda-forge |
| libwebp | 1.2.3 | BSD-3-Clause | conda-forge |
| libwebp-base | 1.2.3 | BSD-3-Clause | conda-forge |
| libxcb | 1.13 | MIT | conda-forge |
| libzlib | 1.2.12 | Zlib | conda-forge |
| lz4-c | 1.9.3 | BSD-2-Clause | conda-forge |
| matplotlib-base | 3.5.2 | LicenseRef-PSF-based | conda-forge |
| munkres | 1.1.4 | Apache-2.0 | conda-forge |
| ncurses | 6.3 | X11 AND BSD-3-Clause | conda-forge |
| numpy | 1.23.1 | BSD-3-Clause | conda-forge |
| openjpeg | 2.4.0 | BSD-2-Clause | conda-forge |
| openssl | 1.1.1q | OpenSSL | conda-forge |
| packaging | 21.3 | Apache-2.0 | conda-forge |
| pillow | 9.2.0 | LicenseRef-PIL | conda-forge |
| pip | 22.2.1 | MIT | conda-forge |
| pthread-stubs | 0.4 | MIT | conda-forge |
| pyparsing | 3.0.9 | MIT | conda-forge |
| pysam | 0.19.1 | MIT | bioconda |
| python | 3.10.5 | Python-2.0 | conda-forge |
| python-dateutil | 2.8.2 | Apache-2.0 | conda-forge |
| python_abi | 3.10 | BSD-3-Clause | conda-forge |
| readline | 8.1.2 | GPL-3.0-only | conda-forge |
| setuptools | 63.3.0 | MIT | conda-forge |
| six | 1.16.0 | MIT | conda-forge |
| sqlite | 3.39.2 | Unlicense | conda-forge |
| tk | 8.6.12 | TCL | conda-forge |
| tzdata | 2022a | LicenseRef-Public-Domain | conda-forge |
| unicodedata2 | 14.0.0 | Apache-2.0 | conda-forge |
| wheel | 0.37.1 | MIT | conda-forge |
| xorg-libxau | 1.0.9 | MIT | conda-forge |
| xorg-libxdmcp | 1.1.3 | MIT | conda-forge |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | conda-forge |
| zlib | 1.2.12 | Zlib | conda-forge |
| zstd | 1.5.2 | BSD-3-Clause | conda-forge |

</details>

<details>
<summary><code>rMATS-v4.3.0.sif</code>: 49 conda packages, base image continuumio/miniconda3:4.8.3 (Debian 10)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | conda-forge |
| _openmp_mutex | 4.5 | BSD-3-Clause | conda-forge |
| ca-certificates | 2024.7.2 | MPL-2.0 | anaconda/main |
| certifi | 2024.7.4 | MPL-2.0 | anaconda/main |
| cffi | 1.14.0 | MIT | anaconda/main |
| chardet | 3.0.4 | LGPL2 | anaconda/main |
| conda | 4.8.3 | BSD 3-Clause | anaconda/main |
| conda-package-handling | 1.6.1 | BSD-3-Clause | anaconda/main |
| cryptography | 2.9.2 | Apache-2.0 AND BSD-3-Clause AND PSF-2.0 | anaconda/main |
| gsl | 2.7.1 | GPL-3.0-or-later | anaconda/main |
| idna | 2.9 | BSD Like | anaconda/main |
| ld_impl_linux-64 | 2.33.1 | GPL 3 | anaconda/main |
| libblas | 3.9.0 | BSD-3-Clause | conda-forge |
| libcblas | 3.9.0 | BSD-3-Clause | conda-forge |
| libedit | 3.1.20181209 | NetBSD | anaconda/main |
| libffi | 3.3 | Custom | anaconda/main |
| libgcc-ng | 14.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgfortran-ng | 14.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgfortran5 | 14.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgomp | 14.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| liblapack | 3.9.0 | BSD-3-Clause | conda-forge |
| libopenblas | 0.3.21 | BSD-3-Clause | anaconda/main |
| libstdcxx-ng | 14.1.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libzlib | 1.2.13 | Zlib | conda-forge |
| ncurses | 6.2 | Free software (MIT-like) | anaconda/main |
| numpy | 1.24.4 | BSD-3-Clause | conda-forge |
| openssl | 1.1.1w | OpenSSL | anaconda/main |
| pip | 20.0.2 | MIT | anaconda/main |
| pycosat | 0.6.3 | MIT | anaconda/main |
| pycparser | 2.20 | BSD 3-clause | anaconda/main |
| pyopenssl | 19.1.0 | Apache 2.0 | anaconda/main |
| pysocks | 1.7.1 | BSD 3-Clause | anaconda/main |
| python | 3.8.3 | PSF | anaconda/main |
| python_abi | 3.8 | BSD-3-Clause | conda-forge |
| readline | 8.0 | GPL-3.0 | anaconda/main |
| requests | 2.23.0 | Apache-2.0 | anaconda/main |
| rmats | 4.3.0 | Free for non-commercial use, see LICENSE file | bioconda |
| ruamel_yaml | 0.15.87 | MIT | anaconda/main |
| setuptools | 46.4.0 | MIT | anaconda/main |
| six | 1.14.0 | MIT | anaconda/main |
| sqlite | 3.31.1 | Public-Domain (http://www.sqlite.org/copyright.html) | anaconda/main |
| star | 2.7.10b | GPLv3 | bioconda |
| tk | 8.6.8 | Tcl/Tk | anaconda/main |
| tqdm | 4.46.0 | MPLv2.0, MIT | anaconda/main |
| urllib3 | 1.25.8 | MIT | anaconda/main |
| wheel | 0.34.2 | MIT | anaconda/main |
| xz | 5.2.5 | LGPL-2.1 and GPL-2.0 | anaconda/main |
| yaml | 0.1.7 | MIT | anaconda/main |
| zlib | 1.2.13 | Zlib | conda-forge |

</details>

<details>
<summary><code>multiQC-v1.18.sif</code>: 94 conda packages, base image multiqc:1.18--pyhdfd78af_0 (BioContainers, Debian 12)</summary>

| Package | Version | License (from package metadata) | Channel |
|---|---|---|---|
| _libgcc_mutex | 0.1 | not recorded | conda-forge |
| _openmp_mutex | 4.5 | BSD-3-Clause | conda-forge |
| brotli | 1.1.0 | MIT | conda-forge |
| brotli-bin | 1.1.0 | MIT | conda-forge |
| brotli-python | 1.1.0 | MIT | conda-forge |
| bzip2 | 1.0.8 | bzip2-1.0.6 | conda-forge |
| ca-certificates | 2023.7.22 | ISC | conda-forge |
| certifi | 2023.7.22 | ISC | conda-forge |
| charset-normalizer | 3.3.2 | MIT | conda-forge |
| click | 8.1.7 | BSD-3-Clause | conda-forge |
| coloredlogs | 15.0.1 | MIT | conda-forge |
| colormath | 3.0.0 | BSD | conda-forge |
| contourpy | 1.2.0 | BSD-3-Clause | conda-forge |
| cycler | 0.12.1 | BSD-3-Clause | conda-forge |
| fonttools | 4.44.3 | MIT | conda-forge |
| freetype | 2.12.1 | GPL-2.0-only OR FTL | conda-forge |
| future | 0.18.3 | MIT | conda-forge |
| humanfriendly | 10.0 | MIT | conda-forge |
| humanize | 4.8.0 | MIT | conda-forge |
| idna | 3.4 | BSD-3-Clause | conda-forge |
| importlib-metadata | 6.8.0 | Apache-2.0 | conda-forge |
| jinja2 | 3.1.2 | BSD-3-Clause | conda-forge |
| kiwisolver | 1.4.5 | BSD-3-Clause | conda-forge |
| lcms2 | 2.15 | MIT | conda-forge |
| ld_impl_linux-64 | 2.40 | GPL-3.0-only | conda-forge |
| lerc | 4.0.0 | Apache-2.0 | conda-forge |
| libblas | 3.9.0 | BSD-3-Clause | conda-forge |
| libbrotlicommon | 1.1.0 | MIT | conda-forge |
| libbrotlidec | 1.1.0 | MIT | conda-forge |
| libbrotlienc | 1.1.0 | MIT | conda-forge |
| libcblas | 3.9.0 | BSD-3-Clause | conda-forge |
| libdeflate | 1.19 | MIT | conda-forge |
| libexpat | 2.5.0 | MIT | conda-forge |
| libffi | 3.4.2 | MIT | conda-forge |
| libgcc-ng | 13.2.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgfortran-ng | 13.2.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgfortran5 | 13.2.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libgomp | 13.2.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libjpeg-turbo | 3.0.0 | IJG AND BSD-3-Clause AND Zlib | conda-forge |
| liblapack | 3.9.0 | BSD-3-Clause | conda-forge |
| libnsl | 2.0.1 | LGPL-2.1-only | conda-forge |
| libopenblas | 0.3.24 | BSD-3-Clause | conda-forge |
| libpng | 1.6.39 | zlib-acknowledgement | conda-forge |
| libsqlite | 3.44.0 | Unlicense | conda-forge |
| libstdcxx-ng | 13.2.0 | GPL-3.0-only WITH GCC-exception-3.1 | conda-forge |
| libtiff | 4.6.0 | HPND | conda-forge |
| libuuid | 2.38.1 | BSD-3-Clause | conda-forge |
| libwebp-base | 1.3.2 | BSD-3-Clause | conda-forge |
| libxcb | 1.15 | MIT | conda-forge |
| libzlib | 1.2.13 | Zlib | conda-forge |
| lzstring | 1.0.4 | MIT | conda-forge |
| markdown | 3.5.1 | BSD-3-Clause | conda-forge |
| markdown-it-py | 3.0.0 | MIT | conda-forge |
| markupsafe | 2.1.3 | BSD-3-Clause | conda-forge |
| matplotlib-base | 3.8.1 | PSF-2.0 | conda-forge |
| mdurl | 0.1.0 | MIT | conda-forge |
| multiqc | 1.18 | GNU General Public License v3 (GPLv3) | bioconda (local build) |
| munkres | 1.1.4 | Apache-2.0 | conda-forge |
| ncurses | 6.4 | X11 AND BSD-3-Clause | conda-forge |
| networkx | 3.2.1 | BSD-3-Clause | conda-forge |
| numpy | 1.26.0 | BSD-3-Clause | conda-forge |
| openjpeg | 2.5.0 | BSD-2-Clause | conda-forge |
| openssl | 3.1.4 | Apache-2.0 | conda-forge |
| packaging | 23.2 | Apache-2.0 | conda-forge |
| pillow | 10.1.0 | HPND | conda-forge |
| pip | 23.3.1 | MIT | conda-forge |
| pthread-stubs | 0.4 | MIT | conda-forge |
| pyaml-env | 1.1.0 | MIT | conda-forge |
| pygments | 2.16.1 | BSD-2-Clause | conda-forge |
| pyparsing | 3.1.1 | MIT | conda-forge |
| pysocks | 1.7.1 | BSD-3-Clause | conda-forge |
| python | 3.11.6 | Python-2.0 | conda-forge |
| python-dateutil | 2.8.2 | Apache-2.0 | conda-forge |
| python_abi | 3.11 | BSD-3-Clause | conda-forge |
| pyyaml | 5.4.1 | MIT | conda-forge |
| readline | 8.2 | GPL-3.0-only | conda-forge |
| requests | 2.31.0 | Apache-2.0 | conda-forge |
| rich | 13.7.0 | MIT | conda-forge |
| rich-click | 1.7.1 | MIT | conda-forge |
| setuptools | 68.2.2 | MIT | conda-forge |
| simplejson | 3.19.2 | MIT | conda-forge |
| six | 1.16.0 | MIT | conda-forge |
| spectra | 0.0.11 | MIT | conda-forge |
| tk | 8.6.13 | TCL | conda-forge |
| typing_extensions | 4.8.0 | PSF-2.0 | conda-forge |
| tzdata | 2023c | LicenseRef-Public-Domain | conda-forge |
| urllib3 | 2.1.0 | MIT | conda-forge |
| wheel | 0.41.3 | MIT | conda-forge |
| xorg-libxau | 1.0.11 | MIT | conda-forge |
| xorg-libxdmcp | 1.1.3 | MIT | conda-forge |
| xz | 5.2.6 | LGPL-2.1 and GPL-2.0 | conda-forge |
| yaml | 0.2.5 | MIT | conda-forge |
| zipp | 3.17.0 | MIT | conda-forge |
| zstd | 1.5.5 | BSD-3-Clause | conda-forge |

</details>

## License texts

### STAR

Applies to STAR 2.7.10b (images `STAR-v2.7.11b.sif` and `rMATS-v4.3.0.sif`).

```text
MIT License

Copyright (c) 2019 Alexander Dobin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Cutadapt

Applies to Cutadapt 3.4 (image `trim_galore-v0.6.7.sif`).

```text
Copyright (c) 2010-2020 Marcel Martin <marcel.martin@scilifelab.se>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

### SAMtools

Applies to SAMtools 1.11 (image `SAMtools-v1.3.1.sif`).

```text
The MIT/Expat License

Copyright (C) 2008-2020 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.


[The use of a range of years within a copyright notice in this distribution
should be interpreted as being equivalent to a list of years including the
first and last year specified and all consecutive years between them.

For example, a copyright notice that reads "Copyright (C) 2005, 2007-2009,
2011-2012" should be interpreted as being identical to a notice that reads
"Copyright (C) 2005, 2007, 2008, 2009, 2011, 2012" and a copyright notice
that reads "Copyright (C) 2005-2012" should be interpreted as being identical
to a notice that reads "Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010,
2011, 2012".]
```

### HTSlib

Applies to HTSlib 1.11 (image `SAMtools-v1.3.1.sif`).

```text
[Files in this distribution outwith the cram/ subdirectory are distributed
according to the terms of the following MIT/Expat license.]

The MIT/Expat License

Copyright (C) 2012-2020 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.


[Files within the cram/ subdirectory in this distribution are distributed
according to the terms of the following Modified 3-Clause BSD license.]

The Modified-BSD License

Copyright (C) 2012-2020 Genome Research Ltd.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
   nor the names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


[The use of a range of years within a copyright notice in this distribution
should be interpreted as being equivalent to a list of years including the
first and last year specified and all consecutive years between them.

For example, a copyright notice that reads "Copyright (C) 2005, 2007-2009,
2011-2012" should be interpreted as being identical to a notice that reads
"Copyright (C) 2005, 2007, 2008, 2009, 2011, 2012" and a copyright notice
that reads "Copyright (C) 2005-2012" should be interpreted as being identical
to a notice that reads "Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010,
2011, 2012".]
```

### Picard

Applies to Picard 2.27.4 (image `picard-v3.1.1.sif`).

```text
MIT License

Copyright (c) 2017-2021 Broad Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### rMATS-turbo

Applies to rMATS-turbo 4.3.0 (image `rMATS-v4.3.0.sif`).

```text
Use of this software is available to academic and non-profit institutions for research purposes subject to the terms of the 2-Clause BSD License (see copy below).  For use or transfers of the software to commercial entities, please inquire with Dr. Yi Xing at xingyi@chop.edu. © 2022 CHOP.

The FreeBSD Copyright

Copyright 1992-2019 The FreeBSD Project.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of the FreeBSD Project.
```

### GNU General Public License, version 3

Applies to Trim Galore! 0.6.7, FastQC 0.11.9, Subread 2.0.1, HTSeq 2.0.2, MultiQC 1.18 and GSL 2.7.1.

```text
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The GNU General Public License is a free, copyleft license for
software and other kinds of works.

  The licenses for most software and other practical works are designed
to take away your freedom to share and change the works.  By contrast,
the GNU General Public License is intended to guarantee your freedom to
share and change all versions of a program--to make sure it remains free
software for all its users.  We, the Free Software Foundation, use the
GNU General Public License for most of our software; it applies also to
any other work released this way by its authors.  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
them if you wish), that you receive source code or can get it if you
want it, that you can change the software or use pieces of it in new
free programs, and that you know you can do these things.

  To protect your rights, we need to prevent others from denying you
these rights or asking you to surrender the rights.  Therefore, you have
certain responsibilities if you distribute copies of the software, or if
you modify it: responsibilities to respect the freedom of others.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must pass on to the recipients the same
freedoms that you received.  You must make sure that they, too, receive
or can get the source code.  And you must show them these terms so they
know their rights.

  Developers that use the GNU GPL protect your rights with two steps:
(1) assert copyright on the software, and (2) offer you this License
giving you legal permission to copy, distribute and/or modify it.

  For the developers' and authors' protection, the GPL clearly explains
that there is no warranty for this free software.  For both users' and
authors' sake, the GPL requires that modified versions be marked as
changed, so that their problems will not be attributed erroneously to
authors of previous versions.

  Some devices are designed to deny users access to install or run
modified versions of the software inside them, although the manufacturer
can do so.  This is fundamentally incompatible with the aim of
protecting users' freedom to change the software.  The systematic
pattern of such abuse occurs in the area of products for individuals to
use, which is precisely where it is most unacceptable.  Therefore, we
have designed this version of the GPL to prohibit the practice for those
products.  If such problems arise substantially in other domains, we
stand ready to extend this provision to those domains in future versions
of the GPL, as needed to protect the freedom of users.

  Finally, every program is threatened constantly by software patents.
States should not allow patents to restrict development and use of
software on general-purpose computers, but in those that do, we wish to
avoid the special danger that patents applied to a free program could
make it effectively proprietary.  To prevent this, the GPL assures that
patents cannot be used to render the program non-free.

  The precise terms and conditions for copying, distribution and
modification follow.

                       TERMS AND CONDITIONS

  0. Definitions.

  "This License" refers to version 3 of the GNU General Public License.

  "Copyright" also means copyright-like laws that apply to other kinds of
works, such as semiconductor masks.

  "The Program" refers to any copyrightable work licensed under this
License.  Each licensee is addressed as "you".  "Licensees" and
"recipients" may be individuals or organizations.

  To "modify" a work means to copy from or adapt all or part of the work
in a fashion requiring copyright permission, other than the making of an
exact copy.  The resulting work is called a "modified version" of the
earlier work or a work "based on" the earlier work.

  A "covered work" means either the unmodified Program or a work based
on the Program.

  To "propagate" a work means to do anything with it that, without
permission, would make you directly or secondarily liable for
infringement under applicable copyright law, except executing it on a
computer or modifying a private copy.  Propagation includes copying,
distribution (with or without modification), making available to the
public, and in some countries other activities as well.

  To "convey" a work means any kind of propagation that enables other
parties to make or receive copies.  Mere interaction with a user through
a computer network, with no transfer of a copy, is not conveying.

  An interactive user interface displays "Appropriate Legal Notices"
to the extent that it includes a convenient and prominently visible
feature that (1) displays an appropriate copyright notice, and (2)
tells the user that there is no warranty for the work (except to the
extent that warranties are provided), that licensees may convey the
work under this License, and how to view a copy of this License.  If
the interface presents a list of user commands or options, such as a
menu, a prominent item in the list meets this criterion.

  1. Source Code.

  The "source code" for a work means the preferred form of the work
for making modifications to it.  "Object code" means any non-source
form of a work.

  A "Standard Interface" means an interface that either is an official
standard defined by a recognized standards body, or, in the case of
interfaces specified for a particular programming language, one that
is widely used among developers working in that language.

  The "System Libraries" of an executable work include anything, other
than the work as a whole, that (a) is included in the normal form of
packaging a Major Component, but which is not part of that Major
Component, and (b) serves only to enable use of the work with that
Major Component, or to implement a Standard Interface for which an
implementation is available to the public in source code form.  A
"Major Component", in this context, means a major essential component
(kernel, window system, and so on) of the specific operating system
(if any) on which the executable work runs, or a compiler used to
produce the work, or an object code interpreter used to run it.

  The "Corresponding Source" for a work in object code form means all
the source code needed to generate, install, and (for an executable
work) run the object code and to modify the work, including scripts to
control those activities.  However, it does not include the work's
System Libraries, or general-purpose tools or generally available free
programs which are used unmodified in performing those activities but
which are not part of the work.  For example, Corresponding Source
includes interface definition files associated with source files for
the work, and the source code for shared libraries and dynamically
linked subprograms that the work is specifically designed to require,
such as by intimate data communication or control flow between those
subprograms and other parts of the work.

  The Corresponding Source need not include anything that users
can regenerate automatically from other parts of the Corresponding
Source.

  The Corresponding Source for a work in source code form is that
same work.

  2. Basic Permissions.

  All rights granted under this License are granted for the term of
copyright on the Program, and are irrevocable provided the stated
conditions are met.  This License explicitly affirms your unlimited
permission to run the unmodified Program.  The output from running a
covered work is covered by this License only if the output, given its
content, constitutes a covered work.  This License acknowledges your
rights of fair use or other equivalent, as provided by copyright law.

  You may make, run and propagate covered works that you do not
convey, without conditions so long as your license otherwise remains
in force.  You may convey covered works to others for the sole purpose
of having them make modifications exclusively for you, or provide you
with facilities for running those works, provided that you comply with
the terms of this License in conveying all material for which you do
not control copyright.  Those thus making or running the covered works
for you must do so exclusively on your behalf, under your direction
and control, on terms that prohibit them from making any copies of
your copyrighted material outside their relationship with you.

  Conveying under any other circumstances is permitted solely under
the conditions stated below.  Sublicensing is not allowed; section 10
makes it unnecessary.

  3. Protecting Users' Legal Rights From Anti-Circumvention Law.

  No covered work shall be deemed part of an effective technological
measure under any applicable law fulfilling obligations under article
11 of the WIPO copyright treaty adopted on 20 December 1996, or
similar laws prohibiting or restricting circumvention of such
measures.

  When you convey a covered work, you waive any legal power to forbid
circumvention of technological measures to the extent such circumvention
is effected by exercising rights under this License with respect to
the covered work, and you disclaim any intention to limit operation or
modification of the work as a means of enforcing, against the work's
users, your or third parties' legal rights to forbid circumvention of
technological measures.

  4. Conveying Verbatim Copies.

  You may convey verbatim copies of the Program's source code as you
receive it, in any medium, provided that you conspicuously and
appropriately publish on each copy an appropriate copyright notice;
keep intact all notices stating that this License and any
non-permissive terms added in accord with section 7 apply to the code;
keep intact all notices of the absence of any warranty; and give all
recipients a copy of this License along with the Program.

  You may charge any price or no price for each copy that you convey,
and you may offer support or warranty protection for a fee.

  5. Conveying Modified Source Versions.

  You may convey a work based on the Program, or the modifications to
produce it from the Program, in the form of source code under the
terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified
    it, and giving a relevant date.

    b) The work must carry prominent notices stating that it is
    released under this License and any conditions added under section
    7.  This requirement modifies the requirement in section 4 to
    "keep intact all notices".

    c) You must license the entire work, as a whole, under this
    License to anyone who comes into possession of a copy.  This
    License will therefore apply, along with any applicable section 7
    additional terms, to the whole of the work, and all its parts,
    regardless of how they are packaged.  This License gives no
    permission to license the work in any other way, but it does not
    invalidate such permission if you have separately received it.

    d) If the work has interactive user interfaces, each must display
    Appropriate Legal Notices; however, if the Program has interactive
    interfaces that do not display Appropriate Legal Notices, your
    work need not make them do so.

  A compilation of a covered work with other separate and independent
works, which are not by their nature extensions of the covered work,
and which are not combined with it such as to form a larger program,
in or on a volume of a storage or distribution medium, is called an
"aggregate" if the compilation and its resulting copyright are not
used to limit the access or legal rights of the compilation's users
beyond what the individual works permit.  Inclusion of a covered work
in an aggregate does not cause this License to apply to the other
parts of the aggregate.

  6. Conveying Non-Source Forms.

  You may convey a covered work in object code form under the terms
of sections 4 and 5, provided that you also convey the
machine-readable Corresponding Source under the terms of this License,
in one of these ways:

    a) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by the
    Corresponding Source fixed on a durable physical medium
    customarily used for software interchange.

    b) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by a
    written offer, valid for at least three years and valid for as
    long as you offer spare parts or customer support for that product
    model, to give anyone who possesses the object code either (1) a
    copy of the Corresponding Source for all the software in the
    product that is covered by this License, on a durable physical
    medium customarily used for software interchange, for a price no
    more than your reasonable cost of physically performing this
    conveying of source, or (2) access to copy the
    Corresponding Source from a network server at no charge.

    c) Convey individual copies of the object code with a copy of the
    written offer to provide the Corresponding Source.  This
    alternative is allowed only occasionally and noncommercially, and
    only if you received the object code with such an offer, in accord
    with subsection 6b.

    d) Convey the object code by offering access from a designated
    place (gratis or for a charge), and offer equivalent access to the
    Corresponding Source in the same way through the same place at no
    further charge.  You need not require recipients to copy the
    Corresponding Source along with the object code.  If the place to
    copy the object code is a network server, the Corresponding Source
    may be on a different server (operated by you or a third party)
    that supports equivalent copying facilities, provided you maintain
    clear directions next to the object code saying where to find the
    Corresponding Source.  Regardless of what server hosts the
    Corresponding Source, you remain obligated to ensure that it is
    available for as long as needed to satisfy these requirements.

    e) Convey the object code using peer-to-peer transmission, provided
    you inform other peers where the object code and Corresponding
    Source of the work are being offered to the general public at no
    charge under subsection 6d.

  A separable portion of the object code, whose source code is excluded
from the Corresponding Source as a System Library, need not be
included in conveying the object code work.

  A "User Product" is either (1) a "consumer product", which means any
tangible personal property which is normally used for personal, family,
or household purposes, or (2) anything designed or sold for incorporation
into a dwelling.  In determining whether a product is a consumer product,
doubtful cases shall be resolved in favor of coverage.  For a particular
product received by a particular user, "normally used" refers to a
typical or common use of that class of product, regardless of the status
of the particular user or of the way in which the particular user
actually uses, or expects or is expected to use, the product.  A product
is a consumer product regardless of whether the product has substantial
commercial, industrial or non-consumer uses, unless such uses represent
the only significant mode of use of the product.

  "Installation Information" for a User Product means any methods,
procedures, authorization keys, or other information required to install
and execute modified versions of a covered work in that User Product from
a modified version of its Corresponding Source.  The information must
suffice to ensure that the continued functioning of the modified object
code is in no case prevented or interfered with solely because
modification has been made.

  If you convey an object code work under this section in, or with, or
specifically for use in, a User Product, and the conveying occurs as
part of a transaction in which the right of possession and use of the
User Product is transferred to the recipient in perpetuity or for a
fixed term (regardless of how the transaction is characterized), the
Corresponding Source conveyed under this section must be accompanied
by the Installation Information.  But this requirement does not apply
if neither you nor any third party retains the ability to install
modified object code on the User Product (for example, the work has
been installed in ROM).

  The requirement to provide Installation Information does not include a
requirement to continue to provide support service, warranty, or updates
for a work that has been modified or installed by the recipient, or for
the User Product in which it has been modified or installed.  Access to a
network may be denied when the modification itself materially and
adversely affects the operation of the network or violates the rules and
protocols for communication across the network.

  Corresponding Source conveyed, and Installation Information provided,
in accord with this section must be in a format that is publicly
documented (and with an implementation available to the public in
source code form), and must require no special password or key for
unpacking, reading or copying.

  7. Additional Terms.

  "Additional permissions" are terms that supplement the terms of this
License by making exceptions from one or more of its conditions.
Additional permissions that are applicable to the entire Program shall
be treated as though they were included in this License, to the extent
that they are valid under applicable law.  If additional permissions
apply only to part of the Program, that part may be used separately
under those permissions, but the entire Program remains governed by
this License without regard to the additional permissions.

  When you convey a copy of a covered work, you may at your option
remove any additional permissions from that copy, or from any part of
it.  (Additional permissions may be written to require their own
removal in certain cases when you modify the work.)  You may place
additional permissions on material, added by you to a covered work,
for which you have or can give appropriate copyright permission.

  Notwithstanding any other provision of this License, for material you
add to a covered work, you may (if authorized by the copyright holders of
that material) supplement the terms of this License with terms:

    a) Disclaiming warranty or limiting liability differently from the
    terms of sections 15 and 16 of this License; or

    b) Requiring preservation of specified reasonable legal notices or
    author attributions in that material or in the Appropriate Legal
    Notices displayed by works containing it; or

    c) Prohibiting misrepresentation of the origin of that material, or
    requiring that modified versions of such material be marked in
    reasonable ways as different from the original version; or

    d) Limiting the use for publicity purposes of names of licensors or
    authors of the material; or

    e) Declining to grant rights under trademark law for use of some
    trade names, trademarks, or service marks; or

    f) Requiring indemnification of licensors and authors of that
    material by anyone who conveys the material (or modified versions of
    it) with contractual assumptions of liability to the recipient, for
    any liability that these contractual assumptions directly impose on
    those licensors and authors.

  All other non-permissive additional terms are considered "further
restrictions" within the meaning of section 10.  If the Program as you
received it, or any part of it, contains a notice stating that it is
governed by this License along with a term that is a further
restriction, you may remove that term.  If a license document contains
a further restriction but permits relicensing or conveying under this
License, you may add to a covered work material governed by the terms
of that license document, provided that the further restriction does
not survive such relicensing or conveying.

  If you add terms to a covered work in accord with this section, you
must place, in the relevant source files, a statement of the
additional terms that apply to those files, or a notice indicating
where to find the applicable terms.

  Additional terms, permissive or non-permissive, may be stated in the
form of a separately written license, or stated as exceptions;
the above requirements apply either way.

  8. Termination.

  You may not propagate or modify a covered work except as expressly
provided under this License.  Any attempt otherwise to propagate or
modify it is void, and will automatically terminate your rights under
this License (including any patent licenses granted under the third
paragraph of section 11).

  However, if you cease all violation of this License, then your
license from a particular copyright holder is reinstated (a)
provisionally, unless and until the copyright holder explicitly and
finally terminates your license, and (b) permanently, if the copyright
holder fails to notify you of the violation by some reasonable means
prior to 60 days after the cessation.

  Moreover, your license from a particular copyright holder is
reinstated permanently if the copyright holder notifies you of the
violation by some reasonable means, this is the first time you have
received notice of violation of this License (for any work) from that
copyright holder, and you cure the violation prior to 30 days after
your receipt of the notice.

  Termination of your rights under this section does not terminate the
licenses of parties who have received copies or rights from you under
this License.  If your rights have been terminated and not permanently
reinstated, you do not qualify to receive new licenses for the same
material under section 10.

  9. Acceptance Not Required for Having Copies.

  You are not required to accept this License in order to receive or
run a copy of the Program.  Ancillary propagation of a covered work
occurring solely as a consequence of using peer-to-peer transmission
to receive a copy likewise does not require acceptance.  However,
nothing other than this License grants you permission to propagate or
modify any covered work.  These actions infringe copyright if you do
not accept this License.  Therefore, by modifying or propagating a
covered work, you indicate your acceptance of this License to do so.

  10. Automatic Licensing of Downstream Recipients.

  Each time you convey a covered work, the recipient automatically
receives a license from the original licensors, to run, modify and
propagate that work, subject to this License.  You are not responsible
for enforcing compliance by third parties with this License.

  An "entity transaction" is a transaction transferring control of an
organization, or substantially all assets of one, or subdividing an
organization, or merging organizations.  If propagation of a covered
work results from an entity transaction, each party to that
transaction who receives a copy of the work also receives whatever
licenses to the work the party's predecessor in interest had or could
give under the previous paragraph, plus a right to possession of the
Corresponding Source of the work from the predecessor in interest, if
the predecessor has it or can get it with reasonable efforts.

  You may not impose any further restrictions on the exercise of the
rights granted or affirmed under this License.  For example, you may
not impose a license fee, royalty, or other charge for exercise of
rights granted under this License, and you may not initiate litigation
(including a cross-claim or counterclaim in a lawsuit) alleging that
any patent claim is infringed by making, using, selling, offering for
sale, or importing the Program or any portion of it.

  11. Patents.

  A "contributor" is a copyright holder who authorizes use under this
License of the Program or a work on which the Program is based.  The
work thus licensed is called the contributor's "contributor version".

  A contributor's "essential patent claims" are all patent claims
owned or controlled by the contributor, whether already acquired or
hereafter acquired, that would be infringed by some manner, permitted
by this License, of making, using, or selling its contributor version,
but do not include claims that would be infringed only as a
consequence of further modification of the contributor version.  For
purposes of this definition, "control" includes the right to grant
patent sublicenses in a manner consistent with the requirements of
this License.

  Each contributor grants you a non-exclusive, worldwide, royalty-free
patent license under the contributor's essential patent claims, to
make, use, sell, offer for sale, import and otherwise run, modify and
propagate the contents of its contributor version.

  In the following three paragraphs, a "patent license" is any express
agreement or commitment, however denominated, not to enforce a patent
(such as an express permission to practice a patent or covenant not to
sue for patent infringement).  To "grant" such a patent license to a
party means to make such an agreement or commitment not to enforce a
patent against the party.

  If you convey a covered work, knowingly relying on a patent license,
and the Corresponding Source of the work is not available for anyone
to copy, free of charge and under the terms of this License, through a
publicly available network server or other readily accessible means,
then you must either (1) cause the Corresponding Source to be so
available, or (2) arrange to deprive yourself of the benefit of the
patent license for this particular work, or (3) arrange, in a manner
consistent with the requirements of this License, to extend the patent
license to downstream recipients.  "Knowingly relying" means you have
actual knowledge that, but for the patent license, your conveying the
covered work in a country, or your recipient's use of the covered work
in a country, would infringe one or more identifiable patents in that
country that you have reason to believe are valid.

  If, pursuant to or in connection with a single transaction or
arrangement, you convey, or propagate by procuring conveyance of, a
covered work, and grant a patent license to some of the parties
receiving the covered work authorizing them to use, propagate, modify
or convey a specific copy of the covered work, then the patent license
you grant is automatically extended to all recipients of the covered
work and works based on it.

  A patent license is "discriminatory" if it does not include within
the scope of its coverage, prohibits the exercise of, or is
conditioned on the non-exercise of one or more of the rights that are
specifically granted under this License.  You may not convey a covered
work if you are a party to an arrangement with a third party that is
in the business of distributing software, under which you make payment
to the third party based on the extent of your activity of conveying
the work, and under which the third party grants, to any of the
parties who would receive the covered work from you, a discriminatory
patent license (a) in connection with copies of the covered work
conveyed by you (or copies made from those copies), or (b) primarily
for and in connection with specific products or compilations that
contain the covered work, unless you entered into that arrangement,
or that patent license was granted, prior to 28 March 2007.

  Nothing in this License shall be construed as excluding or limiting
any implied license or other defenses to infringement that may
otherwise be available to you under applicable patent law.

  12. No Surrender of Others' Freedom.

  If conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot convey a
covered work so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you may
not convey it at all.  For example, if you agree to terms that obligate you
to collect a royalty for further conveying from those to whom you convey
the Program, the only way you could satisfy both those terms and this
License would be to refrain entirely from conveying the Program.

  13. Use with the GNU Affero General Public License.

  Notwithstanding any other provision of this License, you have
permission to link or combine any covered work with a work licensed
under version 3 of the GNU Affero General Public License into a single
combined work, and to convey the resulting work.  The terms of this
License will continue to apply to the part which is the covered work,
but the special requirements of the GNU Affero General Public License,
section 13, concerning interaction through a network will apply to the
combination as such.

  14. Revised Versions of this License.

  The Free Software Foundation may publish revised and/or new versions of
the GNU General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

  Each version is given a distinguishing version number.  If the
Program specifies that a certain numbered version of the GNU General
Public License "or any later version" applies to it, you have the
option of following the terms and conditions either of that numbered
version or of any later version published by the Free Software
Foundation.  If the Program does not specify a version number of the
GNU General Public License, you may choose any version ever published
by the Free Software Foundation.

  If the Program specifies that a proxy can decide which future
versions of the GNU General Public License can be used, that proxy's
public statement of acceptance of a version permanently authorizes you
to choose that version for the Program.

  Later license versions may give you additional or different
permissions.  However, no additional obligations are imposed on any
author or copyright holder as a result of your choosing to follow a
later version.

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

  17. Interpretation of Sections 15 and 16.

  If the disclaimer of warranty and limitation of liability provided
above cannot be given local legal effect according to their terms,
reviewing courts shall apply local law that most closely approximates
an absolute waiver of all civil liability in connection with the
Program, unless a warranty or assumption of liability accompanies a
copy of the Program in return for a fee.

                     END OF TERMS AND CONDITIONS

            How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to the public, the best way to achieve this is to make it
free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest
to attach them to the start of each source file to most effectively
state the exclusion of warranty; and each file should have at least
the "copyright" line and a pointer to where the full notice is found.

    {one line to give the program's name and a brief idea of what it does.}
    Copyright (C) {year}  {name of author}

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Also add information on how to contact you by electronic and paper mail.

  If the program does terminal interaction, make it output a short
notice like this when it starts in an interactive mode:

    {project}  Copyright (C) {year}  {fullname}
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate
parts of the General Public License.  Of course, your program's commands
might be different; for a GUI interface, you would use an "about box".

  You should also get your employer (if you work as a programmer) or school,
if any, to sign a "copyright disclaimer" for the program, if necessary.
For more information on this, and how to apply and follow the GNU GPL, see
<http://www.gnu.org/licenses/>.

  The GNU General Public License does not permit incorporating your program
into proprietary programs.  If your program is a subroutine library, you
may consider it more useful to permit linking proprietary applications with
the library.  If this is what you want to do, use the GNU Lesser General
Public License instead of this License.  But first, please read
<http://www.gnu.org/philosophy/why-not-lgpl.html>.
```
