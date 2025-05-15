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


def init_symbol_params(info: bool = False):
    o = kf.init()
    if info:
        print(f"Высота орбиты: {int(o.HEIGHT // 1e3)} км")
        print(f"Период орбиты: {round((2*np.pi/o.W_ORB) / 3600, 2)} часов")
        print(f"Плотность атмосферы: {o.RHO} кг/м³")

    # Задание символьных параметров окружения
    t, ω, μ, ρ, r_orb, v_orb = var("t w_0 μ ρ r_0 v_0", nonzero=True)
    o.dT = var('dt')
    o.p.t = t
    o.ORBIT_RADIUS = r_orb
    o.V_ORB = v_orb
    o.W_ORB = ω
    o.MU = μ
    o.RHO = ρ
    o.P = r_orb  # * (1 - self.ECCENTRICITY**2)  # Фокальный параметр
    o.a.r_irf[0] = Matrix([r_orb, 0, 0])

    # Задание символьных параметров КА
    o.f.J = diag(*kf.get_vars("J^d", 3, numb=False))
    o.c.J = diag(*kf.get_vars("J^c", 3, numb=False))
    o.f.mass, o.c.mass = var('m_d m_c')
    o.f.c_resist, o.c.c_resist = var('C_d C_c')
    o.f.size = kf.get_vars(name='s^d', n=3)
    o.c.size = kf.get_vars(name='s^c', n=3)

    return o, t, ω, μ, ρ
