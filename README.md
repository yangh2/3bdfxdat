
# Table of Contents

1.  [Purpose](#org80d7b85)
2.  [List of tools](#org2b01f78)
    1.  [init\_bins.awk [IMPLICIT]](#orge070d83)
    2.  [init\_bins.sh [IMPLICIT]](#org8a05acf)
    3.  [split\_xdat.awk [IMPLICIT]](#org8f21fe6)
    4.  [3bdfxdat\_bundle.sh [IMPLICIT]](#orgf8498a9)
    5.  [slurm\_3bdfxdat](#org78bcbe6)
    6.  [3bdfxdat\_sum](#org96b34d4)
    7.  [3bdf2s3](#org8bd5129)
    8.  [3bdf2s3\_eps](#org46a344e)
    9.  [Post-analysis tools](#org54a1789)
3.  [How to](#org7c54760)
    1.  [How to install](#org35fdf95)
    2.  [How to run](#orgca34307)
    3.  [How to analyze](#orgf8b58c2)
    4.  [Example: liquid Boron](#orgb7a79b6)
4.  [Usage in details](#orgd2897af)
    1.  [slurm\_3bdfxdat](#org419870f)
        1.  [Required files](#orgf07e1e0)
        2.  [Modify the slurm script at your needs](#org33f86db)
        3.  [Submit jobs](#org0975a8a)
        4.  [Output](#org043243c)
    2.  [3bdfxdat\_sum](#org3ccb181)
    3.  [3bdf2s3](#org3fd992f)
    4.  [3bdf2s3\_eps [OBSOLETE]](#org90c29b0)



<a id="org80d7b85"></a>

# Purpose

The purpose of this document is to explain how to calculate 3-body distribution function
from a XDATCAR file (VASP output). You are free to modify and
redistribute the code for general purpose.

The basic idea of this program is to divide a BIG XDATCAR file into a
reasonable amount (2 or 3 times of cpus on a cluster) of samller
XDATCAR files and calculate 3-body distribution functions for each
individual job. The final distribution function is obtained by averaging
over distribution functions from these jobs.

This program works best for small systems with long MD
simulations and worst for very large systems. Meanwhile it requires a
great amount of memory to allocate large arrays (depending on number
of bins in your 3-body distribution function histogram) for each job. 


<a id="org2b01f78"></a>

# List of tools

Tools that are implicitly required by others are labeled with [IMPLICIT].


<a id="orge070d83"></a>

## init\_bins.awk [IMPLICIT]

Generate input files. Required by init\_bin.sh.


<a id="org8a05acf"></a>

## init\_bins.sh [IMPLICIT]

Initialize basic settings for slurm\_3bdfxdat runs. Required by slurm\_3bdfxdat.


<a id="org8f21fe6"></a>

## split\_xdat.awk [IMPLICIT]

Split XDATCAR file into individual split files. Required by slurm\_3bdfxdat.


<a id="orgf8498a9"></a>

## 3bdfxdat\_bundle.sh [IMPLICIT]

Bundle a group of split files and calculate their 3-body distribution function. Required by slurm\_3bdfxdat.


<a id="org78bcbe6"></a>

## slurm\_3bdfxdat

Submit a 3-body distribution function calculation job. 


<a id="org96b34d4"></a>

## 3bdfxdat\_sum

Summarize all statistical results from a slurm\_3bdfxdat sjob.


<a id="org8bd5129"></a>

## 3bdf2s3

Calculate 3-body entropy based on results from 3bdfxdat\_sum and
generate supplementary files.


<a id="org46a344e"></a>

## 3bdf2s3\_eps

Calculate 3-body entropy based on results from 3bdfxdat\_sum and
generate supplementary files. Leave out all the ambiguity about
numerical integrals and interpolations caused by the shape of
histograms. It is expected to give an upper boundary of S3.


<a id="org54a1789"></a>

## Post-analysis tools

1.  3bdf\_view\_theta
    Print out angular distribution of bond R1 and bond R2.
    Usage:
    ./3bdf\_view\_theta fname R1 R2
    output
    theta[0,PI] R3[A] f
2.  3bdfMsup [OBSOLETE] 
    Print out g3-g2g2g2
3.  3bdf\_print\_sup [OBSOLETE]
    Print out g2g2g2
4.  3bdf2s3\_interp.py
    Interpolate S3info and S3fluct and write the tot S3.
    output four columns to the STDOUT which are R[A], S3tot, S3Fluct, S3info.


<a id="org7c54760"></a>

# How to


<a id="org35fdf95"></a>

## How to install

To compile, type "make cleanall" and "make" in the main directory.
Then type "./install.sh" to link binaries to your "~/bin".
All tools except slurm\_3bdfxdat are now in your default bin directory.


<a id="orgca34307"></a>

## How to run

Copy "slurm\_3bdfxdat" to your run directory and modify the slurm
script at your needs. Prepare all input files and submit the slurm
script e.g. "sbatch -p n12 -n 12 slurm\_3bdfxdat".


<a id="orgf8b58c2"></a>

## How to analyze

The slurm\_3bdfxdat will provide separate "3bdf\_bins.dat" files for a
certain number samples. Use "3bdfxdat\_sum" to gather all histograms
from each "3bdf\_bins.dat" file and write the overall histogram to an
output file "3bdf\_bins.dat-sum". Use "3bdf2s3" to normalize the
histogram from " 3bdf\_bins.dat-sum" which gives a 3-body distribution
function and compute S3 from 3bdf. Use "3bdf\_view\_theta" and "xmgrace"
to visualize your calculation.


<a id="orgb7a79b6"></a>

## Example: liquid Boron

Go to examples/data and submit slurm\_3bdfxdat


<a id="orgd2897af"></a>

# Usage in details


<a id="org419870f"></a>

## slurm\_3bdfxdat


<a id="orgf07e1e0"></a>

### Required files

-   Trun
    This file contains in value which is the temperature of MD simulation.
-   Mass
    This file contains three rows:
    1.  type of species.
    2.  atomic numbers.
    3.  atomic masses [a.u.] in atomic units.
-   POSCAR
    The POSCAR should be consistent with the XDATCAR file but not
    necessary to be one of the structures in the XDATCAR file.
-   pdf.HH
    Pair distribution functions. You can generate those using "pdfxdat"
    command.
-   XDATCAR 
    This has to be a XDATCAR file from AIMD NVT simulation. Be careful
    about the different in format of XDATCAR from a NVT simulation and a
    NPT simulation.


<a id="org33f86db"></a>

### Modify the slurm script at your needs

-   Cut-off radius
    We use R1, R2 and THETA to control the shape of our histogram. By
    default, R1, R2 have the same range [rmin, rmax] and the range of
    theta depends on R1 and R2. By changing "rmin" and "rmax" in the slurm template, you can
    balance your calculation efficiency and accuracy.
    It is recommended to have "rmax" no longer than a quarter of the
    size of your box and use "rmin" to remove empty regime.
-   Number of bins
    The number of bins for R1 and R2 are specified by dR where
    NR=CEIL[(rmax-rmin)/dR]. For THETA, it is more complex. See the
    [add a link to the ref] pdf for details.
-   Distribute jobs
    In seeking of efficiency, you have to tell the script the number of
    individual 3bdfxdat runs depending on the number of cpus on your
    cluster. To do that, you can change "nsample", "njobs" and "nxdat"
    in the script which follow "nsamples=njobs\*nxdat".
    -   "nsampels": the number samples you want to gather from a XDATCAR
        file to calculate 3-body distribute function.
    -   "njobs" : "nsampels" structures will be divided evenly into $njobs
        subroutines.
    -   "nxdat" : each subroutines have to deal with "nxdat"
        single-structure XDATCAR files.


<a id="org0975a8a"></a>

### Submit jobs

Now you can submit the slurm\_3bdfxdat. On euler, it is "sbatch -p n12
-n slurm\_3bdfxdat".


<a id="org043243c"></a>

### Output

It will return a series of 3bdf\_bins.dat\_??? files where "???" is the
index of first structures in the XDATCAR file. Each 3bdf\_bins.dat\_???
contains a statistical histogram after analyzing "nxdat" amount of files.
In 3bdf\_bins.dat\_???, the first line is the number of analyzed xdatcar
files "nxdat", the second line is
the number of species, the third line is "rmin", "rmax" and "dR". The
remaining of the file contains NRxNR lines and each line is an angular
distribution histogram of bond R1 and bond R2 at their length intervals.


<a id="org3ccb181"></a>

## 3bdfxdat\_sum

To summarize these 3bdf\_bins.dat\_??? files, try
"3bdfxdat\_sum fout 3bdf\_bins.dat\_\*". It will gather information from
all 3bdf\_bins.dat\_??? files and write to the "fout". These
3bdf\_bins.dat\_??? files can have different number of samples. For
example, the "3bdf\_bins.dat\_1" may have analyzed 100 samples while the
"3bdf\_bins.dat\_2" only contains 20 samples. In this way, it is
convenient to add additional information to existing results.


<a id="org3fd992f"></a>

## 3bdf2s3

The statistical work is quite standard and should have less ambiguity
in definitions. Yet the normalization  of the 3-body histogram and
3-body entropy calculation are more subtle and involves detailed
numerical realization and terminology definitions. It is recommended
to check the code and test with a few examples before you make any
conclusions. See the [add a link to the ref] pdf file for more description of
implementation of normalization and integrals.

To run "3bdf2s3", try "3bdf2s3 fin dim", for instance, "3bdf2s3
3bdf\_bins.dat-sum 3". This will normalize the histogram in
"3bdf\_bins.dat-sum" and compute the "S3Fluct" and "S3Info" based on
given "dim". Because pair distribution function and 3-body
distribution function are obtained with different programs and hence
with different shapes of histograms (different rmin, rmax and dr), we
use "dim" to navigate the integration based on a dimxdimxdim mesh grid
where functions at each grid are interpolated from pair correlation functions.

The program generates three files: S3B.dat, S3B-ext.dat S3Fluct.dat and
S3Info.dat. Each contains two columns: the Rc and fval.  The
mathematical expression are shown in the [add a link to the ref].

S3B-ext.dat only depends on two-body correlation function while
S3B.dat is interpolated to the shape of 3bdf bins. In principle, these
two agree when sizes of 3bdf bins equal to sizes of pdf bins. In
general, S3B-ext is more accurate because pdf usually has finer bin settings.


<a id="org90c29b0"></a>

## 3bdf2s3\_eps [OBSOLETE]

This command tends to minimize the differences between superposition
g2g2g2 and 3-body distribution histogram g3. In realization, the
command defines a difference function (Eps=g3-g2g2g2), and  for each
binned g3, it will find a possible range [supp\_min, supp\_max] of g2g2g2  within the bin's
volume. Three conditions may occur.

1.  g3<supp\_min: Eps=g3-supp\_min;
2.  g3>supp\_max: Eps=g3-supp\_max;
3.  supp\_min<g3<supp\_max: Eps=0;

In the way, it estimates the upper bound of S3 and reduces the errors
from numerical solutions and simulation noises.

To run "3bdf2s3\_eps", try "3bdf2s3 fin", for instance, "3bdf2s3
3bdf\_bins.dat-sum". This will normalize the histogram in
"3bdf\_bins.dat-sum" and compute the "S3Fluct" and "S3Info".

The program generates three files: S3B.dat, S3Fluct.dat and
S3Info.dat. Each contains two columns: the Rc and fval.  

