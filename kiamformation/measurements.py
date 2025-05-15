"""Комплекс первичной информации.

Основные функции:
- measure_antennas_power: Расчёт измерений типа RSSI.
- measure_magnetic_field: Расчёт измерений магнитного поля.
- measure_gps: Расчёт измерений ГНСС.
"""
from kiamformation.spacecrafts import Apparatus
from kiamformation.config import Objects
from kiamformation.dynamics import PhysicModel
import numpy as np


def measure_antennas_power(c: Apparatus, f: Apparatus, p: PhysicModel, o: Objects, noise: float = None,
                           produce: bool = False, j: int = None, state=None):
    """Функция рассчитывает измерения типа RSSI между КА в группе.

    Parameters
    ----------
    c : Apparatus
        Класс материнских КА
    f : Apparatus
        Класс дочерних КА
    p : PhysicModel
        Класс численного моделирования
    vrs : Variables
        Класс параметров моделирования
    produce : bool
        Флаг, указывающий, надо ли записывать полученные величины в PhysicModel.record
    noise : float
        Стандартное отклонение белого шума в возвращаемых значениях (передаётся при produce)
    j : int
        Количество параметров на 1 дочерний КА (передаётся при НЕ produce, из KalmanFilter.calc)
    state
        Вектор состояния всех дочерних КА в одном векторе (передаётся при НЕ produce, из KalmanFilter.calc)

    Returns
    -------
    any
        None (при produce) | Вектор измерений (при НЕ produce).

    Examples
    --------
    Численное моделирование - запись измерений в в vrs.MEASURES_VECTOR \n
    >>> measure_antennas_power(c=c, f=f, p=p, o=o, noise=noise, produce=True)

    Расчёт измерений для данного вектора состояний \n
    >>> z = measure_antennas_power(c=c, f=f, p=p, o=o, j=j, state=state)
    """
    from kiamformation.flexmath import setvectype, norm, float2rational
    from kiamformation.math_tools import quart2dcm, vec2quat
    from kiamformation.spacecrafts import get_gain

    S_1, S_2, dr, distance = None, None, None, None
    g_all, anw = [], []

    def get_U(obj, ind):
        from kiamformation.dynamics import get_matrices
        U, S, A, R_orb = get_matrices(o=o, t=p.t, obj=obj, n=ind)
        return U  # float2rational(U, (1, 2), (1, 1))

    for obj1 in [c, f]:
        for obj2 in [f]:
            for i_1 in range(obj1.n):
                for i_2 in range(obj2.n if obj1 == c else i_1):
                    # >>>>>>>>>>>> Расчёт положений и ориентаций <<<<<<<<<<<<
                    
                    # Relative position
                    if produce:
                        dr = obj1.r_orf[i_1] - obj2.r_orf[i_2]
                        if isinstance(dr, np.ndarray):
                            p.record.loc[p.iter, f'{obj1.name}-{obj2.name} RealDistance {i_1} {i_2}'] = norm(dr)
                    else:
                        r1 = setvectype(state[i_1 * j + 0: i_1 * j + 3]) if obj1 == f else obj1.r_orf[i_1]
                        r2 = setvectype(state[i_2 * j + 0: i_2 * j + 3])
                        dr = r1 - r2
                    d2 = dr[0]**2 + dr[1]**2 + dr[2]**2
                            
                    # Relative orientation    
                    U_1 = get_U(obj1, i_1)
                    U_2 = get_U(obj2, i_2) 
                    if produce or not o.NAVIGATION_ANGLES:
                        A_1 = quart2dcm(obj1.q[i_1])
                        A_2 = quart2dcm(obj2.q[i_2])
                    else:
                        q1 = vec2quat(setvectype(state[i_1 * j + 3: i_1 * j + 6])) \
                            if obj1 == f else obj1.q[i_1]
                        q2 = vec2quat(setvectype(state[i_2 * j + 3: i_2 * j + 6]))
                        A_1 = quart2dcm(q1)
                        A_2 = quart2dcm(q2)   
                    S_1 = A_1 @ U_1.T
                    S_2 = A_2 @ U_2.T

                    # >>>>>>>>>>>> Расчёт G и сигнала <<<<<<<<<<<<
                    G1 = [float2rational(g, (1, 2), (1, 1)) for g in get_gain(o=o, obj=obj1, r=S_1 @ dr if obj1==f else dr)]
                    G2 = [float2rational(g, (1, 2), (1, 1)) for g in get_gain(o=o, obj=obj2, r=S_2 @ dr)]
                    g_vec = [g1 * g2 for g1 in G1 for g2 in G2]
                    g_all.extend(g_vec)

                    estimates = [(gg / d2) + np.random.normal(0, noise) for gg in g_vec] if produce else \
                                [gg / d2 for gg in g_vec]
                    anw.extend(estimates)

                    # >>>>>>>>>>>> Расчёт производных по λ <<<<<<<<<<<<
                    '''if not produce and vrs.NAVIGATION_ANGLES and obj2.gain_mode != 'isotropic':
                        screw = get_antisymmetric_matrix
                        a = [np.array(i) for i in
                             get_gain(vrs=vrs, obj=obj2, r=S_2 @ dr, return_dir=True)]
                        a = [a[i] for _ in G1 for i in range(len(G2))]
                        g1_vec = [g1 * 1 for g1 in G1 for _ in G2]
                        g2_vec = [1 * g2 for _ in G1 for g2 in G2]
                        e = dr / np.linalg.norm(dr)

                        tmp = [g1_vec[i] / d2 *
                               (6 * (
                                   (screw(a[i]) @ S_2 @ e).T @ (-screw(a[i]) @ A_2 @ screw(U_2.T @ e))
                               ) / g2_vec[i]**(2/3)) for i in range(len(g_vec))]
                        ...'''

                    # >>>>>>>>>>>> Запись <<<<<<<<<<<<
                    if produce and isinstance(dr, np.ndarray):
                        p.record.loc[p.iter, f'{obj1.name}-{obj2.name} EstimateDistance {i_1} {i_2}'] = \
                            np.mean(estimates)
                        p.record.loc[p.iter, f'{obj1.name}-{obj2.name} ErrorEstimateDistance {i_1} {i_2}'] = \
                            abs(np.mean(estimates) - norm(dr))
                        p.record.loc[p.iter, f'{obj1.name}-{obj2.name} ErrorEstimateDistance 1 {i_1} {i_2}'] = \
                            abs(np.min(estimates) - norm(dr))
                        p.record.loc[p.iter, f'{obj1.name}-{obj2.name} ErrorEstimateDistance 2 {i_1} {i_2}'] = \
                            abs(np.max(estimates) - norm(dr))

    if produce:
        o.MEASURES_VECTOR = setvectype(anw)
        if isinstance(dr, np.ndarray):
            p.record.loc[p.iter, f'G N'] = len(g_all)
            for i in range(len(g_all)):
                p.record.loc[p.iter, f'G {i}'] = g_all[i]
    else:
        return setvectype(anw)

def measure_magnetic_field(c: Apparatus, f: Apparatus, o: Objects, noise: float = 0.) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    pass
    '''
    for obj in [c, f]:
        for i in range(obj.n):
            obj.b_env[i] = np.zeros(3) + np.random.normal(0, noise, 3)'''

def measure_gps(f: Apparatus, noise: float) -> None:
    """Функция обновляет для объектов FemtoSat параметры _не_введено_"""
    pass
