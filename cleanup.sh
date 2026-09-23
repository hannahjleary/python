for n in $(seq 0 10 500); do
    if [ $n -ne 0 ] && [ $n -ne 250 ] && [ $n -ne 400 ]; then
        # rm -f -r ../../../../../ix/eschneider/hjl28/data/radiative/super/4/hdf5/$n.h5.0
        rm -f -r ../../../../../ix/eschneider/hjl28/data/radiative/super/48_raw/
        # rm -f -r ../../../../../ix/eschneider/hjl28/data/radiative/super/8/hdf5/$n.h5
        # rm -f -r ../../../../../ix/eschneider/hjl28/data/radiative/super/8/hdf5/raw/
    fi
done