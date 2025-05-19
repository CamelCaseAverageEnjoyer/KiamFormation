"""
Импорт пакета kiam-formation. 
Задание параметров окружения и параметров КА в символьном виде
"""
import sys  
sys.path.insert(0, "../")
import kiamformation as kf
import matplotlib.pyplot as plt
import numpy as np
import quaternion
from sympy import *


def get_state_vector(func, obj: str, n: int = 1):
    kw = {'n': 3, 'numb': False}
    if func == kf.get_func:
        kw['t'] = t
    return ([func(name=f'r_{i}^{obj}', **kw) for i in range(n)],
            [func(name=f'v_{i}^{obj}', **kw) for i in range(n)],
            [kf.vec2quat(func(name=f'q_{i}^{obj}', **kw)) for i in range(n)],
            [func(name=f'ω_{i}^{obj}', **kw) for i in range(n)])
