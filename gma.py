import numpy as np
import os

def read(file_gma, dtype='float32'):
    if type(file_gma) == str: #filename was given
        with open(file_gma, 'rb') as fin:
            dim_arr = np.fromfile(fin, count=2, dtype='int32')
            arr_in = np.fromfile(fin,dtype=dtype).reshape(tuple(dim_arr))
    else: #file IO was provided:
        try:
            dim_arr = np.fromfile(file_gma, count=2, dtype='int32')
            arr_in = np.fromfile(file_gma, dtype=dtype).reshape(tuple(dim_arr))
        except: #likely object from tarifle module
            try:
                stream_in = file_gma.read()
                dim_arr = np.frombuffer(stream_in, count=2, dtype='int32')
                arr_in = np.frombuffer(stream_in, dtype=dtype, offset=8).reshape(tuple(dim_arr))
            except:
                print('ERROR: Cannot recognize the input object for GMA')

    return arr_in


def write(array, gma_path, dtype='float32'):
    if os.path.exists(gma_path):
        os.remove(gma_path)

    with open(gma_path, 'ab') as gma_file:
        arr_shape = array.shape
        np.array(arr_shape, dtype='int32').tofile(gma_file)
        array.astype(dtype).tofile(gma_file)

    return gma_path

