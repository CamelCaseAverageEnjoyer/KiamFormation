"""Функции для моделирования динамики КА"""
import numpy as np
from kiamformation.config import Objects

# >>>>>>>>>>>> Физика орбитального движения <<<<<<<<<<<<
def get_atm_params(o: Objects, h: float, atm_model: str = None) -> tuple:
    """Функция рассчитывает параметры атмосферы {ρ, T, P} - плотность, температура, давление.
    В kiam-formation используется только плотность для расчёта атмосферного торможения.

    Parameters
    ----------
    o : Objects
        Конфигурационный класс моделирования
    h : float
        Высота над уровнем Земли.
    atm_model : str, optional
        Модель атмосферы. Можно выбрать из: {NASA, ПНБО, COESA62, COESA76}

    Returns
    -------
    tuple
        (ρ, T, P): плотность, температура, давление (Внимание! Для ПНБО (ρ, None, None)

    Notes
    -----
    NASA модель: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html (якорная точка: 25 км).
    ПНБО модель: https://www.energia.ru/ktt/archive/2022/04-2022/101-111.pdf (120-600 км).
    COESA62, COESA76 модели: библиотека poliastro.
    """
    atm_model = o.ATMOSPHERE_MODEL if atm_model is None else atm_model
    rho, T, p = None, None, None
    if atm_model == 'NASA':
        if h > 25e3:
            T = -131.21 + 0.00299 * h
            p = 2.488 * ((T + 273.1) / 216.6) ** -11.388
        elif h > 11e3:
            T = -56.46
            p = 22.65 * np.exp(1.73 - 0.000157 * h)
        else:
            T = 15.04 - 0.00649 * h
            p = 101.29 * (T + 273.1) / 288.08
        rho = p / (0.2869 * (T + 273.1))
    if atm_model == 'ПНБО':
        A, H_m, H_0, K_0, a_1, a_2, n_0, n_01, n_02 = 2.519e-10, 200e3, 290e3, 0.26, 100e3, 141.13e3, 6.34, 4.25, 4.37
        n = n_0 + K_0 * ((h - H_0) / a_1) ** n_01 - ((h - H_0) / a_2) ** n_02 if h > 290e3 else \
            n_0 + K_0 * ((H_0 - h) / a_1) ** n_01
        rho = A * (H_m / h) ** n
    if atm_model in ['COESA62', 'COESA76']:
        from astropy import units as u
        from poliastro.earth.atmosphere import COESA62, COESA76
        coesa = COESA62() if atm_model == 'COESA62' else COESA76()
        T, p, rho = coesa.properties(h * u.m)
        T = T.value
        p = p.value
        rho = rho.value
    return rho, T, p

def get_geopotential_acceleration(r, v, w0):
    """Расчёт ускорения КА от притяжения Земли.

    Parameters
    ----------
    r
        Радиус-вектор КА в ОСК
    v
        Вектор скорости КА в ОСК
    w0
        Модуль угловой скорости вращения КА вокруг Земли
    """
    from kiamformation.flexmath import setvectype
    return setvectype([-2 * w0 * v[2],
                       -w0 ** 2 * r[1],
                       3 * w0 ** 2 * r[2] + 2 * w0 * v[0]])

def get_aero_drag_acceleration(o: Objects, obj, i: int, r, v, rho=None):
    """Возвращает ускорение КА от сопротивления атмосферы."""
    from kiamformation.math_tools import quart2dcm, matrix2angle
    from kiamformation.flexmath import setvectype
    S = quart2dcm(obj.q[i])
    cos_alpha = matrix2angle(S) if obj.name == "FemtoSat" else 1

    v2 = setvectype([(0 + o.V_ORB)**2, 0, 0])  # v[0] +
    o.RHO = get_atm_params(o=o, h=r[2] + o.HEIGHT)[0] if rho is None else rho
    return - v2 * abs(obj.get_blown_surface(cos_alpha)) * o.RHO / obj.mass

def get_full_acceleration(o: Objects, obj, i: int, r, v):
    """Возвращает вектор силы в ОСК, принимает параметры в ОСК"""
    from kiamformation.flexmath import setvectype
    force = get_geopotential_acceleration(r=r, v=v, w0=o.W_ORB)
    if o.physics_model['aero drag']:
        force += setvectype(get_aero_drag_acceleration(o=o, r=r, v=v, obj=obj, i=i))
    return force


# >>>>>>>>>>>> Физика вращательного движения <<<<<<<<<<<<
def get_torque(o: Objects, obj, q, w, t, i):
    """Вектор внешнего углового ускорения"""
    from kiamformation.flexmath import norm, cross, inv
    q = obj.q if q is None else q
    w = obj.w_brf if w is None else w
    J = obj.J
    U, S, A, R_orb = get_matrices(o=o, t=t, obj=obj, n=i, q=q)
    R = A @ R_orb
    # m_grav = np.zeros(3)
    m_grav = 3*o.MU/norm(R)**5 * cross(R, J @ R)
    return inv(J) @ (m_grav - cross(w, J @ w))


# >>>>>>>>>>>> Интегрирование движения <<<<<<<<<<<<
def rhs(o: Objects, obj, t: float, i: int, rqvw) -> tuple:
    from kiamformation.flexmath import quat, float2rational
    from kiamformation.math_tools import q_dot

    r, q, v, w = rqvw if isinstance(rqvw, tuple) else \
        (rqvw[[0, 1, 2]], quat(rqvw[[3, 4, 5, 6]]), rqvw[[7, 8, 9]], rqvw[[10, 11, 12]])

    dr = v
    dq = 1 / 2 * q_dot(q, quat(w))
    dq = float2rational(dq, (1, 2))
    dv = get_full_acceleration(o=o, obj=obj, i=i, r=r, v=v)
    dw = get_torque(o=o, obj=obj, q=q, w=w, t=t, i=i)
    return (dr, dq, dv, dw) if isinstance(rqvw, tuple) else np.append(np.append(dr, dq.components), np.append(dv, dw))


def ode4(o: Objects, obj, t: float, i: int, dt: float = None, r=None, v=None, q=None, w=None) -> tuple:
    """Реализация схемы Рунге-Кутты 4-го порядка для поступательного и вращательного движения КА."""
    import quaternion
    from kiamformation.math_tools import vec2quat

    dt = o.dt if dt is None else dt
    r = obj.r_orf[i] if r is None else r
    v = obj.v_orf[i] if v is None else v
    q = obj.q[i] if q is None else q
    w = obj.w_brf[i] if w is None else w
    q4 = q if isinstance(q, np.quaternion) else vec2quat(q)

    rqvw = np.append(np.append(r, quaternion.as_float_array(q4)), np.append(v, w))
    k1 = rhs(o=o, obj=obj, i=i, t=t, rqvw=rqvw)
    k2 = rhs(o=o, obj=obj, i=i, t=t, rqvw=rqvw + k1 * dt / 2)
    k3 = rhs(o=o, obj=obj, i=i, t=t, rqvw=rqvw + k2 * dt / 2)
    k4 = rhs(o=o, obj=obj, i=i, t=t, rqvw=rqvw + k3 * dt)
    rqvw = dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    r_anw = r + rqvw[0:3]
    q_anw = (q4 + np.quaternion(*rqvw[3:7])).normalized()
    q_anw = (-1 * q_anw) if q_anw.w < 0 else q_anw
    v_anw = v + rqvw[7:10]
    w_anw = w + rqvw[10:13]

    return (r_anw, q_anw, v_anw, w_anw) if isinstance(q, np.quaternion) else (r_anw, q_anw.vec, v_anw, w_anw)


# >>>>>>>>>>>> Перевод между системами координат <<<<<<<<<<<<
def get_matrices(o: Objects, t, obj=None, n: int = None, first_init: bool = False, q=None):
    """Функция возвращает матрицы поворота.
    Инициализируется в dymancis.py, используется в spacecrafts, dynamics"""
    from kiamformation.flexmath import arctan, sqrt, setvectype, sin, cos, tan, cross, norm, float2rational, \
        get_same_type_conversion
    from kiamformation.math_tools import quart2dcm, vec2unit
    q = obj.q[n] if q is None else q
    E = t * o.W_ORB  # Эксцентрическая аномалия
    if o.ECCENTRICITY == 0:
        f = E
    else:
        f = 2 * arctan(sqrt((1 + o.ECCENTRICITY) / (1 - o.ECCENTRICITY)) * tan(E / 2))  # Истинная аномалия
    A = quart2dcm(q)
    if first_init:
        U = setvectype([[0, 1, 0],  # Поворот к экваториальной плоскости
                        [0, 0, 1],
                        [1, 0, 0]]) @ \
            setvectype([[cos(f), sin(f), 0],  # Разница между истинной аномалией и местной
                        [-sin(f), cos(f), 0],
                        [0, 0, 1]]) @ \
            setvectype([[1, 0, 0],  # Поворот к плоскости орбиты
                        [0, cos(o.INCLINATION), sin(o.INCLINATION)],
                        [0, -sin(o.INCLINATION), cos(o.INCLINATION)]])
        translation = o.P / (1 + o.ECCENTRICITY * cos(f))
        U = float2rational(U, (2, 1), (1, 2), (1, 1))
    else:
        e_z = vec2unit(o.a.r_irf[0])
        e_x = vec2unit(o.a.v_irf[0])
        e_y = vec2unit(cross(e_z, e_x))
        e_x = vec2unit(cross(e_y, e_z))
        U = setvectype([e_x, e_y, e_z])
        translation = norm(o.a.r_irf[0])
    S = A @ U.T
    tmp = get_same_type_conversion(U)
    R_orb = U.T @ tmp([0, 0, translation])
    return U, S, A, R_orb

def i_o(a, o: Objects, U, vec_type: str):
    """Инерциальная -> Орбитальная"""
    from kiamformation.flexmath import get_same_type_conversion

    t = get_same_type_conversion(a)
    if len(a.shape) == 1 and vec_type == "r":
        return U @ a - t([0, 0, o.ORBIT_RADIUS])
    if len(a.shape) == 1 and vec_type == "v":
        return U @ a - t([o.V_ORB, 0, 0])
    if len(a.shape) == 1 and vec_type == "w":
        return U @ (a - o.W_ORB_VEC_IRF)
    if len(a.shape) == 2:
        return U @ a @ U.T
    raise ValueError(f"Необходимо подать вектор или матрицу! Тип вектора {vec_type} должен быть из [r, v, w]")

def o_i(a, o: Objects, U, vec_type: str):
    """Орбитальная -> Инерциальная"""
    from kiamformation.flexmath import get_same_type_conversion

    t = get_same_type_conversion(a)
    if len(a.shape) == 1 and vec_type == "r":
        return U.T @ (a + t([0, 0, o.ORBIT_RADIUS]))
    if len(a.shape) == 1 and vec_type == "v":
        return U.T @ (a + t([o.V_ORB, 0, 0]))
    if len(a.shape) == 1 and vec_type == "w":
        return U.T @ a + o.W_ORB_VEC_IRF
    if len(a.shape) == 2:
        return U.T @ a @ U
    raise ValueError(f"Необходимо подать вектор или матрицу! Тип вектора {vec_type} должен быть из [r, v, w]")


# >>>>>>>>>>>> Класс динамики кубсатов и чипсатов <<<<<<<<<<<<
class PhysicModel:
    from kiamformation.spacecrafts import FemtoSat, CubeSat, Anchor

    def __init__(self, f: FemtoSat, c: CubeSat, a: Anchor, o: Objects):
        from pandas import DataFrame
        from kiamformation.navigation import KalmanFilter

        # Неизменные параметры
        self.t = 0.
        self.iter = 0
        self.o = o
        self.c = c
        self.f = f
        self.a = a
        self.spacecrafts_cd = [self.c, self.f]
        self.spacecrafts_all = [self.a, self.c, self.f]
        self.time2nav = 0.  # Если навигация включается не сразу

        # Инициализация фильтра
        self.k = KalmanFilter(f=f, c=c, p=self)

        # Инициализация траектории kiam-astro
        self.jd0, self.tr = None, None

        # Запись параметров
        self.record = DataFrame()

    # Шаг по времени
    def time_step(self):
        from scipy.spatial.transform import Rotation
        from kiamformation.navigation import navigate
        from kiamformation.guidance import guide
        from kiamformation.measurements import measure_antennas_power, measure_magnetic_field

        self.iter += 1
        self.t += self.o.dt
        
        if self.iter == 9900:
            import matplotlib.pyplot as plt
            plt.grid()
            plt.show()

        # Движение системы
        for j, obj in enumerate(self.spacecrafts_all):
            for i in range(obj.n):
                U, S, A, _ = get_matrices(o=self.o, t=self.t, obj=obj, n=i)

                # Интегрирование движения
                obj.r_orf[i], obj.q[i], obj.v_orf[i], obj.w_brf[i] = ode4(o=self.o, obj=obj, i=i, t=self.t)

                # Костыль
                if obj == self.c:
                    obj.q[i] = np.quaternion(*np.roll(Rotation.from_matrix(U).as_quat(), 1))
                    obj.w_orf[i] = np.zeros(3)
                    obj.update_irf_w(o=self.o, t=self.t)

                # Расчёт зависимых параметров
                U, _, A, _ = get_matrices(o=self.o, t=self.t, obj=obj, n=i)
                obj.w_irf[i] = A.T @ obj.w_brf[i]
                obj.w_orf[i] = i_o(o=self.o, a=obj.w_irf[i], U=U, vec_type='w')
                obj.r_irf[i] = o_i(o=self.o, a=obj.r_orf[i], U=U, vec_type='r')
                obj.v_irf[i] = o_i(o=self.o, a=obj.v_orf[i], U=U, vec_type='v')

        # Комплекс первичной информации
        noise = np.sqrt(self.o.kalman_coef['r'])
        measure_antennas_power(c=self.c, f=self.f, o=self.o, noise=noise, produce=True, p=self)
        measure_magnetic_field(c=self.c, f=self.f, o=self.o, noise=noise)

        # Изменение режимов работы
        guide(o=self.o, c=self.c, f=self.f, earth_turn=self.t * self.o.W_ORB / 2 / np.pi)

        # Навигация чипсатов
        if self.o.if_navigation:
            navigate(k=self.k, if_correction=self.time2nav <= 0.)
            self.time2nav = self.o.dt - self.o.dt if self.time2nav <= 0 else self.time2nav - self.o.dt

        # Запись параметров
        self.do_report()

    def do_report(self):
        i_t = self.iter
        d = self.record
        d.loc[i_t, f'i'] = self.iter
        d.loc[i_t, f't'] = self.t
        n_tmp = len(self.o.MEASURES_VECTOR)
        d.loc[i_t, f'MEASURES_VECTOR N'] = n_tmp
        d.loc[i_t, [f'MEASURES_VECTOR {i}' for i in range(n_tmp)]] = self.o.MEASURES_VECTOR
        for obj in self.spacecrafts_all:
            d.loc[i_t, f'{obj.name} n'] = obj.n
            for i_n in range(obj.n):
                d.loc[i_t, f'{obj.name} a x orf {i_n}'] = 0 if i_t == 1 else \
                    (obj.v_orf[i_n][0] - d.loc[i_t-1, f'{obj.name} v x orf {i_n}'])/self.o.dt
                for v in ['r', 'q', 'v', 'w']:
                    tmp = {'r': [obj.r_irf[i_n], obj.r_orf[i_n]],
                           'v': [obj.v_irf[i_n], obj.v_orf[i_n]],
                           'q': [obj.q[i_n].vec],
                           'w': [obj.w_irf[i_n], obj.w_orf[i_n], obj.w_brf[i_n]]}[v]
                    for i_fr, frame in enumerate(['irf', 'orf'] if v not in 'qw' else
                                                 (['irf'] if v != 'w' else ['irf', 'orf', 'brf'])):
                        for i_r, c in enumerate('xyz'):
                            d.loc[i_t, f'{obj.name} {v} {c} {frame} {i_n}'] = tmp[i_fr][i_r]

        for obj in [self.f]:
            for i_n in range(obj.n):
                if obj.operating_mode[i_n] != "lost":  # Иначе заполняется Null (в plot в self.v.NO_LINE_FLAG)
                    r_orf_estimation = self.k.get_estimation(i_f=i_n, v='r orf')
                    w_brf_estimation = self.k.get_estimation(i_f=i_n, v='w brf')
                    q_irf_estimation = self.k.get_estimation(i_f=i_n, v='q-3 irf')
                    r_orf = self.f.r_orf[i_n]
                    w_brf = self.f.w_brf[i_n]
                    q_irf = self.f.q[i_n].vec

                    w_orf, w_orf_estimation = [], []  # Чтобы PyCharm не ругался
                    if self.o.rotational_motion_navigate:
                        U, _, A, _ = get_matrices(o=self.o, t=self.t, obj=obj, n=i_n)
                        w_irf = A.T @ w_brf
                        w_orf = i_o(a=w_irf, o=self.o, vec_type='w', U=U)
                        w_irf_estimation = A.T @ w_brf_estimation
                        w_orf_estimation = i_o(a=w_irf_estimation, o=self.o, vec_type='w', U=U)

                    d.loc[i_t, f'{obj.name} KalmanPosEstimation r {i_n}'] = np.linalg.norm(r_orf_estimation)
                    d.loc[i_t, f'{obj.name} KalmanPosError r {i_n}'] = np.linalg.norm(r_orf_estimation - r_orf)
                    if self.o.rotational_motion_navigate:
                        d.loc[i_t, f'{obj.name} KalmanSpinError w {i_n}'] = np.linalg.norm(w_brf_estimation - w_brf)
                        d.loc[i_t, f'{obj.name} KalmanQuatError q {i_n}'] = np.linalg.norm(q_irf_estimation - q_irf)
                    for i_r, c in enumerate('xyz'):
                        d.loc[i_t, f'{obj.name} KalmanPosEstimation {c} {i_n}'] = r_orf_estimation[i_r]
                        d.loc[i_t, f'{obj.name} KalmanPosError {c} {i_n}'] = r_orf_estimation[i_r] - r_orf[i_r]
                        if self.o.rotational_motion_navigate:
                            d.loc[i_t, f'{obj.name} RealSpin BRF {c} {i_n}'] = w_brf[i_r]
                            d.loc[i_t, f'{obj.name} RealSpin ORF {c} {i_n}'] = w_orf[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinEstimation BRF {c} {i_n}'] = w_brf_estimation[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinEstimation ORF {c} {i_n}'] = w_orf_estimation[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinError BRF {c} {i_n}'] = w_brf_estimation[i_r] - w_brf[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinError ORF {c} {i_n}'] = w_orf_estimation[i_r] - w_orf[i_r]
                            d.loc[i_t, f'{obj.name} RealQuat {c} {i_n}'] = q_irf[i_r]
                            d.loc[i_t, f'{obj.name} KalmanQuatEstimation {c} {i_n}'] = q_irf_estimation[i_r]
                            d.loc[i_t, f'{obj.name} KalmanQuatError {c} {i_n}'] = q_irf_estimation[i_r] - q_irf[i_r]
