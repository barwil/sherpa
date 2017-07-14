# SHERPA: Simple HiErarchical Profile Aggregation

Written by Krzysiek Królak with modifications by Bartek Wilczyński 

## Quick summary ##
-------------------

This software is intended to identify hierarchically nested domains based on normalized hi-c matrices. The process is based on adaptive aggregation of  most similar correlation profiles of neighboring hi-c bins. It works on normalized hi-c matrices in numpy format and  generates text files with domain boundaries and provides basic visualisation capabilities by outputting png or pdf images.

## Prerequisites ##
-------------------
The code is written in python (version 2.7) so you need to have python installed. Additional packages are:
* matplotlib (tested on version 1.5.1)
* numpy (tested on 1.11.0)
* scipy (tested on 0.17.0)
* biopython (tested on 1.66)

We recommend using it on linux systems. In particular, it is working on the current long term support Ubuntu version (16.04), that will be supported until 2021. If you have that system, you can ensure that all prerequisites are met, by running:
``` sudo apt-get install python-numpy python-biopython python-matplotlib python-scipy ```

If you are running SHERPA on other platforms, we recommend using some scientific-related distriubution like Enthought Canopy or Continuum Analytics Anaconda, but you should be able to also install all the packages in the standard python install on any platform.

## Installation  ##
-------------------
Once you have all the prerequisites installed, you just need to download our scripts. This can be done by downloading and unpacking the zip file:

https://github.com/regulomics/sherpa/archive/master.zip

or by cloning our git repository:

```git clone https://github.com/regulomics/sherpa.git```

This should leave you with all necessary files in the "sherpa" folder in the location you chose.

## Quick Manual  ##
-------------------

### Input ###

Hi-C data are taken as normalized numpy matrices (.npy) for each chromosome separately. There is an example file called test_data.npy in the package you downloaded that represents a small toy chromosome. If you have your matrices in another format, please refer to the numpy manual 
https://docs.scipy.org/doc/numpy/reference/generated/numpy.load.html

### Finding Domains ###

The script "find_domains.py" is used to generate hierarchical domains from a hi-c matrix. It's really simple, you can see how it works by running:

```python find_domains.py test_data.npy > output.doms```

This will result in creation of a new file called "output.doms" in your sherpa directory.

### Testing if SHERPA works correctly ###

We have attached a pre-computed version of the output file called test_domains.doms in the distribution. If you want to make sure that SHERPA gives the same output on your machine, please compare the test_domains.doms file with the output.doms file you generated in the previous step. For example on linux systems running the command:
```diff output.doms test_domains.doms```

should give you an empty output indicating that the files are indeed the same.

### Output format ###
The output format is really simple. Each line contains the name of the chromosome (in our example "chr") and 3 numbers: the domain number, the first bin included in the domain and the last bin included in the domain. For example, the following line in the test_domains.doms file describes a domain number 209 that is covering bins from 101 to 115 (inclusive):

```
chr 209: 101 115
```

But the same area is covered earlier in smaller domains:

```
chr 146: 101 109
chr 147: 110 115 
```

The output from SHERPA will always only contain domains that are hierarchical, i.e. completely contained in one another or disjoint.
 
### Visualizing Domains ###
To show the results in a graphical form, you can use the "show_domains.py" script. For example, the command:

```python show_domains.py test_data.npy test_doms.doms test_data.npy output.doms```

Should result in opening of the following window:

![alt text](http://regulomics.mimuw.edu.pl/~bartek/sherpa-example.png "sherpa show_domains interface")

You should be able to save this image in a pdf or png format by using the floppy disk icon in the lower right corner.
