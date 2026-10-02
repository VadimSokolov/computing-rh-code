# Regenerate the zero lists g_1_400.npy, g_401_1000.npy, g_1001_1700.npy (ordinates as float64).
import numpy as np, mpmath as mp
from multiprocessing import Pool
def z(n):
    mp.mp.dps = 30
    return float(mp.zetazero(n).imag)
if __name__ == '__main__':
    with Pool(16) as p:
        g = np.array(p.map(z, range(1, 1701), chunksize=10))
    np.save('g_1_400.npy', g[:400]); np.save('g_401_1000.npy', g[400:1000]); np.save('g_1001_1700.npy', g[1000:1700])
    print('zeros', len(g), g[0], g[-1])
