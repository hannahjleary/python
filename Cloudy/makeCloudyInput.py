if __name__ == '__main__':
    # CLOUDY file title
    title = str(input('file name: '))

    # temperature grid
    Tmin = float(input('Temp. min (log(T)): '))
    Tmax = float(input('Temp. max (log(T)): '))
    delT = float(input('log(\u0394T): '))
    # string for input file
    tempStr = 'grid '+str(Tmin)+' '+str(Tmax)+' '+str(delT)+' log'

    # hydrogen density
    hden = float(input(f'log(n_H): '))
    # string for input file
    hStr = 'hden '+str(hden)

    # metalicity scaling
    metals = float(input('metalicity (scaled-solar): '))
    # string for input file
    metalStr = 'metals '+str(metals)

    # redshift
    z = float(input('redshift: '))
    # background choice
    bg = str(input('background (HM12 or HM05): '))
    # string for input file
    bgStr = 'table '+str(bg)+' redshift '+str(z)

    # ask if user needs heating/cooling tables
    hcQ = str(input('need heating/cooling tables (Y/N): '))

    # write the input file
    filename = 'cool_n='+str(hden)+'_Z='+str(metals)
    f = open(filename+'.in', 'w')
    f.write('title '+title+'\n')
    f.write('coronal T=1e5K vary\n')
    f.write(tempStr+'\n')
    f.write(hStr+'\n')
    f.write(bgStr+'\n')
    f.write(metalStr + '\n')
    f.write('stop zone 1'+'\n')
    f.write('print cooling' + '\n')
    if hcQ == 'Y':
        f.write('punch cooling "cooling.col" last no hash no clobber'+'\n')
        f.write('punch heating "heating.het" last no hash no clobber'+'\n')
    f.write('save overview "'+str(filename)+'.txt"'+'\n')

    # show user filename
    print()
    print('input file generated')
    print('file name: ' + filename)
    print()