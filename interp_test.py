import numpy as np
from scipy.interpolate import interp2d

x = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
y = [9, 8, 7, 6, 5, 4, 3, 2, 1]

# bortle  9     8     7     6    5    4    3     2    1   
z = [[175.0, 62.0, 22.0, 10.8, 5.3, 2.5, 1.2, 0.98, 0.8],  # f4.0
     [112.0, 40.0, 14.0, 7.2, 3.7, 1.7, 0.81, 0.64, 0.51], # f5.0
     [78.0,  27.0,  9.6, 5.0, 2.6, 1.2, 0.56, 0.45, 0.36], # f6.0
     [57.0,  20.0,  7.1, 3.7, 1.9, 0.9, 0.41, 0.33, 0.26], # f7.0
     [47.0,  17.0,  5.9, 3.0, 1.6, 0.7, 0.34, 0.27, 0.22], # f8.0
     [38.0,  13.0,  4.6, 2.4, 1.2, 0.6, 0.26, 0.21, 0.17], # f9.0
     [28.0,  10.0,  3.4, 1.7, 0.85,0.4, 0.19, 0.16, 0.13]] # f10.0

z = np.array(z).T

f = interp2d(x, y, z)

fstop_ff = (952.0 / 127.0) * 0.7
white_mtn = 3.5

print(f"Fstop: {fstop_ff}, white_mtn_bortle: {white_mtn}")

            # fstop, bortle
value = f(x = fstop_ff, y = white_mtn) 
print(f" Value = {value}")