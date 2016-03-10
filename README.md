# sherpa
Simple HiErarchical Profile Aggregation

Written by Krzysiek Królak, some modifications by Bartek Wilczyński

You need to have matplotlib and  numpy installed.

Hi-C data are taken as normalized numpy matrices (.npy)

To find domains run:

python find_domains.py test_data.npy > output.doms

To show them on your data (you can show the same files on both halves of the matrix or show two different datasets for contrast):

python show_domains.py data_upper.npy doms_upper.doms data_lower.npy data_lower.doms



