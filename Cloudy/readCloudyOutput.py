import numpy as np
from astropy.table import Table

if __name__ == '__main__':
    # get path to CLOUDY output file
    outPath = str(input('output file path (w/ file extension): '))

    # hden input
    hden = 10**float(input(f'log(n_H): '))

    # store vals here
    heatVals = []
    coolVals = []
    temps = []

    # open output file and store heat/cool vals
    with open(outPath, 'r') as f:
        for line in f:
            if 'HTOT=' in line:
                hIdx = line.find('HTOT= ') + 7
                cIdx = line.find('CTOT= ') + 7
                tIdx = line.find('temperature')+16
                heat = line[hIdx:hIdx+8]
                cool = line[cIdx:cIdx+8]
                temp = line[tIdx:tIdx+9]
                heatVals.append(float(heat))
                coolVals.append(float(cool))
                temps.append(float(temp))
    f.close()

    # make table for plotting
    data = Table()

    data['log(T)']=np.log10(np.array(temps))
    data['log(cool/n^2)']=np.log10(np.array(coolVals)/hden**2)
    data['log(heat/n^2)']=np.log10(np.array(heatVals)/hden**2)

    # table name input
    tname = str(input('output table name (w/o file extension): '))

    # save table
    data.write(tname+'.txt', format='ascii')