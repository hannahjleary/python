import h5py
import numpy as np
import os
from csv import writer

# path = '../../../../../ix/eschneider/hjl28/data/cooling_tables/z_0.000.hdf5'
# out_path = '../../../../../ix/eschneider/hjl28/data/cooling_tables/'

# f = h5py.File(path, 'r') 
# T = np.array(f['Solar/Temperature_bins/'])
# L = np.array(f['Solar/Net_cooling/'])
# f.close()

# with open(out_path + "my_file.txt", "w") as f:
#     for i in range(len(T)):
#         f.write(str(T[i]) + '\t' + str(L[i]) + '\n')

with open(os.path.join("../../../../../ix/eschneider/hjl28/data/adiabatic_PPMP/super/plotting_data/dm_3.csv"), "r") as f_in:
    with open(os.path.join("../../../../../ix/eschneider/hjl28/data/adiabatic_PPMP/super/plotting_data/dm_3_new.csv"), "a") as f_out:
        writer_obj = writer(f_out)
        i=0
        for line in f_in:
            if line == '\n':
                writer_obj.writerow([])
                i=0
            else:
                line = line.split(",")
                if i == 0:
                    init = line[2]
                writer_obj.writerow([line[0],line[1],float(line[2])/float(init)])
                i+=1

