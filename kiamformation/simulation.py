"""Здесь много всего так-то"""
import numpy as np
from kiamformation.config import Objects
# from kiamformation.static.cosmetic import my_print


def timer(func):
    import time

    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f'Выполнено "{func.__name__}" за {run_time:.4f} секунд')
        return value
    return wrapper_timer


def find_close_solution(o: Objects):
    """Функция находит близкие траектории, анализирует историю измерений.
    Args:
        o : Objects
            Класс, содержащий основные классы моделирования
    """
    import matplotlib.pyplot as plt

    o.v.IF_NAVIGATION = False
    o.v.IF_ANY_PRINT = False
    qf = o.f.q
    wf = o.f.w_brf
    rf = o.f.r_orf
    vf = o.f.v_orf
    qc = o.c.q
    wc = o.c.w_brf
    rc = o.c.r_orf
    vc = o.c.v_orf

    dr_list = [0, -0.1, 0.1]
    dv_list = [0, -0.01, 0.01]
    measurements = []

    for dr_x in dr_list:
        for dr_y in dr_list:
            for dr_z in dr_list:
                for dv_x in dv_list:
                    for dv_y in dv_list:
                        for dv_z in dv_list:
                            for i in range(o.f.n):
                                o.init_classes()
                                o.f.q[i] = qf[i]
                                o.f.w_brf[i] = wf[i]
                                o.f.r_orf[i] = rf[i] + np.array([dr_x, dr_y, dr_z])
                                o.f.v_orf[i] = vf[i] + np.array([dv_x, dv_y, dv_z])
                                o.integrate(t=o.v.TIME)
                                # measurements.append(o.v.MEASURES_VECTOR)
                                measurements.append(o.p.record[f'MEASURES_VECTOR {0}'])
    for i, m in enumerate(measurements):
        plt.plot(np.abs(m - measurements[0]), color='k' if i == 0 else 'g', ls=':')
    plt.legend()
    plt.grid()
    plt.show()
