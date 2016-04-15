'''
This program I have written for my work as a molecular biologist. In a given DNA sequence it looks for regions which
are rich in AT pairs and returns them in a .txt file and on a plot. It is useful for searching chromosomes regions
susceptible for DNA helix distortions.

Written for Anaconda3 - a scientific Python distribution.
'''

from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as np
import math
import statistics
import numpy
import re


# opening fasta file as DNA sequence
your_file = input("Type in the FASTA file (with an extension): ")
match_file = re.search(r'\w+\.fasta', your_file)
while not match_file:
    your_file = input("Wrong file format. Type the name of file with extension again: ")
    match_file = re.search(r'\w+\.fasta', your_file)

seqq = SeqIO.read(your_file, "fasta")
seqq = seqq.seq
print(len(seqq))


# finds all divisors of the integer
def divisor_generator(n):
    divisors = []
    for i in range(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                divisors.append(int(n / i))
    for divisor in reversed(divisors):
        yield divisor

div_list = divisor_generator(len(seqq))

# asks for preferred region's length
proposed_reg_len = input("Type in the proposed region length:")
while not proposed_reg_len.isdigit():
    proposed_reg_len = input("Wrong length format. Type in an integer: ")


# finds an integer that will be smallest divisor closest to the proposed region length
region_len_list = [i for i in div_list if i <= int(proposed_reg_len)]
region_len = max(region_len_list)  # length of the region that will be in further use
print(region_len_list)
print("Length of the region:", region_len)

regions = []  # regions' starting points
regions_at = {}  # dict: keys - starting position of region, values - AT% of the region
at_values = []
start = 0

# makes a list of the regions' starting position
while start < len(seqq):
    regions.append(start)
    start += region_len

print(regions)

# computes %AT of every region in DNA sequence
for i in regions:
    at_percentage = 100 - GC(seqq[i:i + region_len+1])
    regions_at[i] = at_percentage
    at_values.append(at_percentage)
    start += region_len

print(regions_at)
print(at_values)
print("Total AT%:", 100-GC(seqq))
print("Number of regions:", len(regions))

at_median = statistics.median(at_values)  # median value of AT% in regions
print("AT% median:", at_median)

at_75percentile = numpy.percentile(at_values, 75)  # calculates the 3rd quartile (0.75 percentile) of AT regions
print("3rd quartile of AT%:", at_75percentile)


# takes regions with AT% only in 3rd quartile
regions_at_up_75percentile = {}  # dict as regions_at but regions with higher AT than 3rd quartile
at_values_up_75percentile = []
for key in regions_at:
    if regions_at[key] > at_75percentile:
        regions_at_up_75percentile[(key, key+region_len)] = regions_at[key]
        at_values_up_75percentile.append(regions_at[key])
print(regions_at_up_75percentile)
print(at_values_up_75percentile)
print(len(at_values_up_75percentile))

# writes regions (upper 25 percentile) with their %AT to the file
file_regions_up_75 = open("result.txt", "w")
file_regions_up_75.writelines(str((key, regions_at_up_75percentile[key])) + "\n" for key in sorted(regions_at_up_75percentile))


# draws a histogram with number of regions of x% of AT pairs
np.ylabel("number of regions")
np.xlabel("%AT")
np.hist(at_values, bins=20)
np.xticks(range(0, 105, 5))
np.show()

# draws vlines-plot with DNA sequence length (x) and AT% of regions with AT% higher than 3rd quartile
x_points = numpy.arange(0, len(seqq), region_len)
np.vlines(x_points, [0], at_values)
np.axhline(y=at_75percentile, color="r", alpha=0.5)
np.ylabel("%AT")
np.xlabel("bp")
np.show()


