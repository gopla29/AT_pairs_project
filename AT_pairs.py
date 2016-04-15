from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import matplotlib.pyplot as np
import math
import statistics
import numpy


# opening txt as DNA sequence
# your_file = input("Name of the file with DNA sequence (with .txt extension):")
seq = open("mc2.txt", "r")
seqq = Seq(seq.read(), IUPAC.ambiguous_dna)
print("Sequence length: %d" % len(str(seqq)))


# finds all divisors of the integer
def find_divisors(x):
    divisors_list = []
    y = 1
    while y <= math.sqrt(x):
        if x % y == 0:
            divisors_list.append(y)
            divisors_list.append(int(x / y))
        y += 1
    return divisors_list

div_list = find_divisors(len(str(seqq)))

# asks for region's length
left = int(input("Type in the proposed region length:"))

region_len_list = [i for i in div_list if i >= left]
region_len = min(region_len_list)  # length of the region that will be in further use
print("Region length in usage:", region_len)

print(region_len_list)

regions = []  # regions starting points
regions_at = {}  # dict: keys - starting position of region, values: AT% of the region
at_values = []
start = 0

# makes a list of the regions' starting position
while start < len(seqq):
    regions.append(start)
    start += region_len

print(regions)

seq.close()

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
        regions_at_up_75percentile[key] = regions_at[key]
        at_values_up_75percentile.append(regions_at[key])
print(regions_at_up_75percentile)
print(at_values_up_75percentile)
print(len(at_values_up_75percentile))

# saves in a file starting and ending points of regions with AT% higher 3rd quartile
starts_stops_up_75 = []
for key in regions_at_up_75percentile.keys():
    starts_stops_up_75.append([key+1, key+region_len+1])

starts_stops_up_75.sort()
file_regions_up_75 = open("result.txt", "w")
file_regions_up_75.writelines(str(line)+"\n" for line in starts_stops_up_75)

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