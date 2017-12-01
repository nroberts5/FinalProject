import numpy as np
import h5py


def main():
    mat = np.random.randint(0,100,(10,8))

    h5f = h5py.File('hdf5matrix.h5', 'w')
    h5f.create_dataset('dataset_1', data=mat)
    h5f.close()



if __name__ == '__main__':
    main()