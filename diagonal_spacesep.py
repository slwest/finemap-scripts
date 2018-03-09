
#this script changes the dilimiter of the symmetrical matrix to the space separated and outputs the result to a different directory - the file is altered in place
import numpy as np
import sys
infile = sys.argv[1] 
outfile = sys.argv[2]

data=np.loadtxt(infile , delimiter="\t")
print  data.shape
##print with no header
#print the number in decimal rather than scientific notation

np.savetxt(outfile , data, fmt='%.16g', delimiter=" ", comments='')
