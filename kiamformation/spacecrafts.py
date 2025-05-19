"""Функции, связанные с архитектурой КА

Основные функции:
- local_dipole: Расчёт углового коэффициента усиления антенн
- get_gain: Расчёт угловых коэффициентов для определённого КА по определённому направлению

Классы:
- Apparatus: Космический аппарат
"""
import numpy as np
from kiamformation.config import Objects

# >>>>>>>>>>>> Диаграмма направленности антенн связи <<<<<<<<<<<<
def r_a(ind: str):
    """Единичный вектор по направлениям x, y или z"""
    return [int(ind == 'x'), int(ind == 'y'), int(ind == 'z')]

def local_dipole(o: Objects, r, ind: str = 'x', model: str = 'quarter-wave monopole'):
    """Возвращает коэффициент усиления от 1-й антенны

    Parameters
    ----------
    vrs : Variables
        Класс параметров моделирования
    r
        Направление к другому КА в ССК КА
    ind : str
        Координата направления антенны в ССК КА
    model : str
        Математическая модель антенны. Реализованы: [half-wave dipole, short dipole, quarter-wave monopole]

    Returns
    -------
    float
        Угловой коэффициент усиления в диапазоне [-1, 1]
    """
    from kiamformation.flexmath import setvectype, norm, cos, dot, cross, pi
    from kiamformation.math_tools import vec2unit
    
    r_12 = vec2unit(r) + setvectype([0.03 * o.G_distortion, 0.05 * o.G_distortion, 0])   # Искажение диаграммы
    r_antenna_brf = setvectype(r_a(ind))
    
    sin_theta = norm(cross(r_antenna_brf, r_12))
    cos_theta = dot(r_antenna_brf, r_12)  # HERE ABS VALUE IS NEED !!!

    if model == 'half-wave dipole':
        return cos(pi(r) / 2 * cos_theta) / sin_theta
    if model == 'short dipole':
        return sin_theta**2
    if model == 'quarter-wave monopole':
        return sin_theta**3


def get_gain(o: Objects, obj, r, gm=None, return_dir: bool = False):
    """Возвращает вектор коэффициентов усиления от каждой антенны КА

    Parameters
    ----------
    vrs : Variables
        Класс параметров моделирования
    obj : Apparatus
        Класс КА
    r
        Направление сигнала в ССК КА
    return_dir : bool, optional
        Флаг возврата вектора направления антенны в ССК КА
    gm : str, optional
        Количество и тип антенн

    Returns
    -------
    list
        Набор угловых коэффициентов усиления от всех антенн по направлению r
    """
    gm = obj.gain_mode if gm is None else gm
    if gm in o.GAIN_MODES and gm != 'isotropic':
        d = ['x', 'xy', 'xyz'][o.GAIN_MODES.index(gm) - 1]
        return [local_dipole(o, r, i) for i in d] if not return_dir else [r_a(i) for i in d]
    return [1]


# >>>>>>>>>>>> Классы аппаратов <<<<<<<<<<<<
class Apparatus:
    """Класс КА для численного моделирования в

    Note:
        gain_mode копируется из GAIN_MODEL_C или GAIN_MODEL_F

    Attributes
    ----------
    name : str
        Тип КА
    n : int
        Количество КА данного типа в группе
    mass : float
        Масса КА
    size : list
        Линейные размеры
    c_resist : float
        Коэффициент лобового сопротивления
    J
        Тензор инерции
    gain_mode : str
        Количество и тип антенн у КА
    b_env
        Вектор магнитного поля Земли в положении КА

    Methods
    -------
    get_blown_surface(cos_alpha)
        Возвращает характерную площадь обдувания
    update_irf_rv(v, t)
        Положения и скорости КА данного типа в ИСК рассчитываются из себя в ОСК
    """
    def __init__(self, o: Objects, n: int):
        from kiamformation.physics import get_matrices, o_i

        # Общие параметры
        self.name = "No exist"
        self.n = n
        self.mass = 1e10
        self.size = [1., 1., 1.]
        self.c_resist = 0
        self.J = np.diag([1, 1, 1])
        self.gain_mode = 'isotropic'
        self.b_env = None
        self.apriori_params = None

        # Индивидуальные параметры движения
        self.w_irf = [np.zeros(3) for _ in range(self.n)]
        self.w_orf = [np.zeros(3) for _ in range(self.n)]
        self.w_brf = [np.zeros(3) for _ in range(self.n)]
        self.q = [np.quaternion(1, 0, 0, 0) for _ in range(self.n)]
        self.r_orf = [np.zeros(3) for _ in range(self.n)]
        self.v_orf = [np.zeros(3) for _ in range(self.n)]
        U, _, _, _ = get_matrices(o=o, t=0, obj=self, n=0, first_init=True)
        self.r_irf = [o_i(o=o, a=self.r_orf[0], U=U, vec_type='r')]
        self.v_irf = [o_i(o=o, a=self.v_orf[0], U=U, vec_type='v')]

        # Индивидуальные параметры режимов работы
        self.operating_mode = [o.OPERATING_MODES[0] for _ in range(self.n)]

    def get_blown_surface(self, cos_alpha):
        return abs(self.size[0] * self.size[1] * abs(cos_alpha) * self.c_resist)

    def update_irf_rv(self, o: Objects, t: float = 0):
        from kiamformation.physics import o_i, get_matrices
        for i in range(self.n):
            U, _, _, _ = get_matrices(o=o, t=t, obj=self, n=i)
            self.r_irf[i] = o_i(o=o, a=self.r_orf[i], U=U, vec_type='r')
            self.v_irf[i] = o_i(o=o, a=self.v_orf[i], U=U, vec_type='v')

    def update_irf_w(self, o: Objects, t: float = 0, w_irf: list = None, w_orf: list = None, w_brf: list = None):
        from kiamformation.physics import o_i, get_matrices
        w_irf = self.w_irf if w_irf is None else w_irf
        w_orf = self.w_orf if w_orf is None else w_orf
        w_brf = self.w_brf if w_brf is None else w_brf
        for i in range(self.n):
            U, _, A, _ = get_matrices(o=o, t=t, obj=self, n=i)
            w_irf[i] = o_i(a=w_orf[i], o=o, U=U, vec_type='w')
            w_brf[i] = A @ w_irf[i]

    def init_correct_q_v(self, o, q: list = None):
        q = self.q if q is None else q
        for i in range(self.n):
            if o.specific_initial["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * o.W_ORB
            q[i] = q[i].normalized()
            if q[i].w < 0:
                q[i] *= -1


class Anchor(Apparatus):
    def __init__(self, o: Objects):
        """Класс мнимого КА, центр которого совпадает с центром ОСК"""
        super().__init__(o=o, n=1)
        self.name = "Anchor"

class CubeSat(Apparatus):
    """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.
    Все величны представлены в СИ."""
    def __init__(self, o: Objects):
        super().__init__(o=o, n=o.cubesat_amount)

        # Предопределённые параметры
        cubesat_property = {'1U': {'mass': 2.,
                                   'mass_center_error': [0.02, 0.02, 0.02],
                                   'dims': [0.1, 0.1, 0.1135]},
                            '1.5U': {'mass': 3.,
                                     'mass_center_error': [0.02, 0.02, 0.03],
                                     'dims': [0.1, 0.1, 0.1702]},
                            '2U': {'mass': 4.,
                                   'mass_center_error': [0.02, 0.02, 0.045],
                                   'dims': [0.1, 0.1, 0.227]},
                            '3U': {'mass': 6.,
                                   'mass_center_error': [0.02, 0.02, 0.07],
                                   'dims': [0.1, 0.1, 0.3405]},
                            '6U': {'mass': 12.,
                                   'mass_center_error': [4.5, 2., 7.],
                                   'dims': [0.2263, 0.1, 0.366]},
                            '12U': {'mass': 24.,
                                    'mass_center_error': [4.5, 4.5, 7.],
                                    'dims': [0.2263, 0.2263, 0.366]}}

        # Общие параметры
        self.name = "CubeSat"
        self.gain_mode = o.GAIN_MODEL_C
        self.mass = cubesat_property[o.CUBESAT_MODEL]['mass']
        self.size = cubesat_property[o.CUBESAT_MODEL]['dims']
        self.mass_center_error = cubesat_property[o.CUBESAT_MODEL]['mass_center_error']
        self.r_mass_center = np.array([np.random.uniform(-i, i) for i in self.mass_center_error])
        self.c_resist = 1.05
        # Пока что J диагонален
        J = np.array([self.size[1]**2 + self.size[2]**2,
                      self.size[0]**2 + self.size[2]**2,
                      self.size[0]**2 + self.size[1]**2]) * self.mass / 12
        self.J = np.diag(J)

        # Индивидуальные параметры движения
        self.r_orf = [o.spread('r', name=self.name) for _ in range(self.n)]
        self.v_orf = [o.spread('v', name=self.name) for _ in range(self.n)]
        self.w_orf = [np.zeros(3) for _ in range(self.n)]
        self.q = [np.quaternion(1, -1, -1, -1) for _ in range(self.n)]

        # СПЕЦИАЛЬНЫЕ НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ УДОВЛЕТВОРЕНИЯ ТРЕБОВАНИЯМ СТАТЬИ
        self.r_orf[0], self.v_orf[0], self.w_orf[0] = [np.array(o.specific_initial[f"CubeSat {i}"]) for i in "rvw"]
        self.q[0] = np.quaternion(*o.specific_initial[f"CubeSat q"])

        # Инициализируется автоматически
        self.init_correct_q_v(o=o)
        self.r_irf, self.v_irf, self.w_irf, self.w_irf = [[np.zeros(3) for _ in range(self.n)] for _ in range(4)]
        self.update_irf_rv(o=o, t=0)
        self.update_irf_w(o=o, t=0)

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

class FemtoSat(Apparatus):
    def __init__(self, o: Objects, c: CubeSat):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        super().__init__(o=o, n=o.chipsat_amount)

        # Предопределённые параметры
        chipsat_property = {'KickSat': {'mass': 0.005,
                                        'mass_center_error': [0.001, -0.001],
                                        'dims': [0.035, 0.035, 0.001]},
                            'Трисат': {'mass': 0.1,
                                       'mass_center_error': [0.005, 0.003],
                                       'dims': [0.4, 0.15, 0.001]}}

        # Общие параметры
        self.name = "FemtoSat"
        self.gain_mode = o.GAIN_MODEL_F
        self.mass = chipsat_property[o.CHIPSAT_MODEL]['mass']
        self.mass_center_error = chipsat_property[o.CHIPSAT_MODEL]['mass_center_error']
        self.size = chipsat_property[o.CHIPSAT_MODEL]['dims']
        self.c_resist = 1.17
        # Пока что J диагонален
        J = np.array([self.size[1]**2, self.size[0]**2, self.size[0]**2 + self.size[1]**2]) * self.mass / 12
        self.J = np.diag(J)
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001

        # Индивидуальные параметры движения
        self.deploy(o=o, c=c, i_c=0)
        self.w_orf_ = [o.spread('w', name=self.name) for _ in range(self.n)]
        self.r_irf, self.v_irf, self.w_irf, self.w_irf_, self.w_brf, self.w_brf_ = \
            [[np.zeros(3) for _ in range(self.n)] for _ in range(6)]
        self.q, self.q_ = [[np.quaternion(*np.random.uniform(-1, 1, 4)) for _ in range(self.n)] for _ in range(2)]

        # СПЕЦИАЛЬНЫЕ НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ УДОВЛЕТВОРЕНИЯ ТРЕБОВАНИЯМ СТАТЬИ
        self.r_orf[0], self.v_orf[0], self.w_orf[0] = [np.array(o.specific_initial[f"ChipSat {i}"]) for i in "rvw"]
        dr, dv = [np.array(o.specific_initial[f"ChipSat d{i}"]) for i in "rv"]
        self.q[0] = np.quaternion(*o.specific_initial[f"ChipSat q"])

        # Инициализируется автоматически
        self.init_correct_q_v(o=o)
        self.init_correct_q_v(o=o, q=self.q_)
        self.update_irf_rv(o=o, t=0)
        self.update_irf_w(o=o, t=0)
        self.update_irf_w(o=o, t=0, w_irf=self.w_irf_, w_orf=self.w_orf_, w_brf=self.w_brf_)

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        tol = 1 if o.START_NAVIGATION == o.NAVIGATIONS[0] else o.start_nav_tolerance
        tol = 0 if o.START_NAVIGATION == o.NAVIGATIONS[2] else tol

        if o.rotational_motion_navigate:
            self.apriori_params = {'r orf': [self.r_orf[i] * tol + o.spread('r', name=self.name) * (1 - tol)
                                             for i in range(self.n)],
                                   'v orf': [self.v_orf[i] * tol + o.spread('v', name=self.name) * (1 - tol)
                                             for i in range(self.n)],
                                   'w brf': [self.w_brf[i] * tol + self.w_brf_[i] * (1 - tol) for i in range(self.n)],
                                   'q-3 irf': [self.q[i].vec * tol + self.q_[i].vec * (1 - tol) for i in range(self.n)]}
        else:
            if tol in [0, 1]:
                self.apriori_params = {'r orf': [self.r_orf[i] * tol + o.spread('r', name=self.name) * (1 - tol)
                                                 for i in range(self.n)],
                                       'v orf': [self.v_orf[i] * tol + o.spread('v', name=self.name) * (1 - tol)
                                                 for i in range(self.n)],
                                       'w brf': [self.w_brf[i] for i in range(self.n)],
                                       'q-3 irf': [self.q[i].vec for i in range(self.n)]}
            else:
                self.apriori_params = {'r orf': [self.r_orf[i] + dr for i in range(self.n)],
                                       'v orf': [self.v_orf[i] + dv for i in range(self.n)],
                                       'w brf': [self.w_brf[i] for i in range(self.n)],
                                       'q-3 irf': [self.q[i].vec for i in range(self.n)]}

    def deploy(self, o: Objects, c: CubeSat, i_c: int) -> None:
        """Функция отделения задаёт начальные условия для дочерних КА из материнских КА
        :param v: объект Variables
        :param c: объект CubeSat
        :param i_c: id-номер материнского КА, от которого отделяются дочерние КА
        :return: {'r orf': ..., 'v orf': ..., 'q-3 irf': ..., 'w irf': ...}, где значения - list of np.ndarray
        """

        if o.DEPLOYMENT == o.DEPLOYMENTS[0]:  # Deploy: "No"
            self.r_orf = [o.spread('r', name=self.name) for _ in range(self.n)]
            self.v_orf = [o.spread('v', name=self.name) for _ in range(self.n)]
            self.w_orf = [o.spread('w', name=self.name) for _ in range(self.n)]
        elif o.DEPLOYMENT == o.DEPLOYMENTS[1]:  # Deploy: "Specific"
            r_before = c.r_orf[i_c]
            v_before = c.v_orf[i_c]
            dv = 1e-2
            v_deploy = o.RVW_ChipSat_SPREAD[1]
            self.r_orf = [r_before.copy() for _ in range(self.n)]
            self.v_orf = [v_before + np.array([0, 0, v_deploy]) + np.random.uniform(-dv, dv, 3) for _ in range(self.n)]
            self.w_orf = [o.spread('w', name=self.name) for _ in range(self.n)]
           

