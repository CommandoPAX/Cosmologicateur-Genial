import numpy as np

from pyMIN import*
import matplotlib.pyplot as plt
from math import*
import sys

print(sys.path)

data = np.random.normal(size=(64,64,64))


v0,v1,v2,v3 = calculateMFs(data)


plt.plot(v0,label="v0")
plt.legend()


plt.plot(v1,label="v1")
plt.legend()


plt.plot(v2,label="v2")
plt.legend()


plt.plot(v3,label="v3")

plt.legend()
plt.savefig("test_min.png")
