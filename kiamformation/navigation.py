from kiamformation.spacecrafts import Apparatus
from kiamformation.static.cosmetic import my_print
import numpy as np


class KalmanFilter:
    """Оцениваемые параметры: ['r orf', 'q-3 irf', 'v orf', 'w brf']; согласовано с spacecrafts.py"""
    def __init__(self, f: Apparatus, c: Apparatus, p: any):
        from kiamformation.flexmath import block_diag

        # Общие параметры
        self.f = f  # Дочерние КА
        self.c = c  # Материнские КА
        self.p = p  # Динамическая модель
        self.o = p.o

        self.estimation_params = self.params_dict2vec(d=f.apriori_params, separate_spacecraft=True)
        self.j = len(self.estimation_params[0])  # Вектор состояния 1 чипсата
        self.estimation_params_d = self.params_vec2dict(self.estimation_params)
        self.STM = None
        self.observability_gramian = None

        # Матрицы фильтра в начальный момент времени
        if not self.o.rotational_motion_navigate:  # Вектор состояния содержит только положение и скорость
            self.D = np.vstack([np.zeros((3, 3)), np.eye(3)])
            self.P = [np.diag([self.o.kalman_coef['p'][0]] * 3 + [self.o.kalman_coef['p'][1]] * 3) for _ in range(f.n)]
            self.Q = np.diag([self.o.kalman_coef['q'][0]] * 3)
        else:  # Вектор состояния содержит ещё и угловые переменные
            self.D = np.vstack([np.zeros((6, 6)), np.eye(6)])
            self.P = [np.diag([self.o.kalman_coef['p'][0]] * 3 + [self.o.kalman_coef['p'][1]] * 3 +
                              [self.o.kalman_coef['p'][2]] * 3 + [self.o.kalman_coef['p'][3]] * 3) for _ in range(f.n)]
            self.Q = np.diag([self.o.kalman_coef['q'][0]] * 3 + [self.o.kalman_coef['q'][1]] * 3)
        self.Phi = self.get_Phi(w=None, w0=None, q=None)

        # Расширешние на учёт несколько аппаратов в фильтре
        self.D = block_diag(*[self.D for _ in range(self.f.n)])
        self.P = block_diag(*[self.P[i] for i in range(self.f.n)])
        self.Q = block_diag(*[self.Q for _ in range(self.f.n)])

    def get_Phi_1(self, i: int, w=None, w0=None, q=None):
        """Внимание: переделать секцию (if self.o.NAVIGATION_ANGLES)"""
        from kiamformation.math_tools import get_antisymmetric_matrix, quart2dcm, vec2quat, matrix2angle
        from kiamformation.flexmath import inv, sqrt, cos, sin, bmat, setvectype, zeros, eye
        from kiamformation.physics import get_matrices

        w0 = self.o.W_ORB if w0 is None else w0
        q = self.f.q if q is None else q

        O = zeros((3, 3), template=w0)
        E = eye(3, template=w0)
        A1 = setvectype([[0, 0, 0],
                         [0, -w0**2, 0],
                         [0, 0, 3*w0**2]], template=w0)
        A2 = setvectype([[0, 0, -2 * w0],
                         [0, 0, 0],
                         [2 * w0, 0, 0]], template=w0)

        if self.o.physics_model['aero drag']:
            cos_alpha = matrix2angle(quart2dcm(self.f.q[i]))
            c = self.f.get_blown_surface(cos_alpha=cos_alpha) * self.o.RHO / self.f.mass
            A2[0, 0] = -2 * c * (self.o.V_ORB + self.estimation_params_d['v orf'][i][0])

        if self.o.rotational_motion_navigate:  # Оценка орбитального и углового движения
            J = self.f.J
            Wj = inv(J) @ (-get_antisymmetric_matrix(w) @ J + get_antisymmetric_matrix(J @ w))
            Ww = get_antisymmetric_matrix(w)

            U, S, A, R_orb = get_matrices(o=self.o, t=self.p.t, obj=self.f, n=i)
            R = quart2dcm(vec2quat(*np.array(self.estimation_params_d['q-3 irf']))) @ R_orb
            eta = R / np.linalg.norm(R)
            Wr = get_antisymmetric_matrix(eta)
            Wq = inv(self.f.J) @ (6*self.o.W_ORB**2 * (Wr @ J @ Wr - get_antisymmetric_matrix(J @ eta) @ Wr))

            if self.o.physics_model['aero drag']:
                q_x, q_y, q_z = q
                t = self.p.t
                w_0 = self.o.W_ORB

                s1 = (2*q_x**2/sqrt(-q_x**2 - q_y**2 - q_z**2 + 1) - 2*sqrt(-q_x**2 - q_y**2 - q_z**2 + 1))*cos(t*w_0) - (2*q_x*q_y/sqrt(-q_x**2 - q_y**2 - q_z**2 + 1) + 2*q_z)*sin(t*w_0)
                s2 = -(2*q_y**2/sqrt(-q_x**2 - q_y**2 - q_z**2 + 1) - 2*sqrt(-q_x**2 - q_y**2 - q_z**2 + 1))*sin(t*w_0) + (2*q_x*q_y/sqrt(-q_x**2 - q_y**2 - q_z**2 + 1) + 2*q_z)*cos(t*w_0)
                s3 = -(2*q_x + 2*q_y*q_z/sqrt(-q_x**2 - q_y**2 - q_z**2 + 1))*sin(t*w_0) + (2*q_x*q_z/sqrt(-q_x**2 - q_y**2 - q_z**2 + 1) + 2*q_y)*cos(t*w_0)

                cos_alpha = matrix2angle(quart2dcm(self.estimation_params_d['q-3 irf'][i]))
                c = self.f.get_blown_surface(cos_alpha=cos_alpha) * self.o.RHO / self.f.mass
                v2 = self.o.V_ORB**2
                s1, s2, s3 = c*s1*v2, c*s2*v2, c*s3*v2
            else:
                s1, s2, s3 = 0, 0, 0

            A3 = setvectype([[s1, s2, s3],
                             [0, 0, 0],
                             [0, 0, 0]])
            F = bmat([[O, O, E, O],
                      [O, -Ww, O, E/2],
                      [A1, A3, A2, O],
                      [O, Wq, O, Wj]])
        else:
            F = bmat([[O, E],
                      [A1, A2]], template=w0)

        return F * self.o.dt + np.eye(self.j)

    def get_Phi(self, w=None, w0=None, q=None):
        from kiamformation.flexmath import block_diag
        w0 = self.o.W_ORB if w0 is None else w0
        w = [self.estimation_params_d['w brf'][i] for i in range(self.f.n)] if w is None else w
        q = [self.estimation_params_d['q-3 irf'][i] for i in range(self.f.n)] if q is None else q
        return block_diag(*[self.get_Phi_1(w=w[i], w0=w0, i=i, q=q[i]) for i in range(self.f.n)])

    def params_vec2dict(self, params: list = None, j: int = None, separate_spacecraft: bool = True):
        p = self.estimation_params if params is None else params
        j = self.j if j is None else j
        if self.o.rotational_motion_navigate:
            if separate_spacecraft:
                r_orf = [p[i][0: 3] for i in range(self.f.n)]
                q_irf = [p[i][3: 6] for i in range(self.f.n)]
                v_orf = [p[i][6: 9] for i in range(self.f.n)]
                w_orf = [p[i][9: 12] for i in range(self.f.n)]
            else:
                r_orf = [p[i*j + 0: i*j + 3] for i in range(self.f.n)]
                q_irf = [p[i*j + 3: i*j + 6] for i in range(self.f.n)]
                v_orf = [p[i*j + 6: i*j + 9] for i in range(self.f.n)]
                w_orf = [p[i*j + 9: i*j + 12] for i in range(self.f.n)]
        else:
            if separate_spacecraft:
                r_orf = [p[i][0: 3] for i in range(self.f.n)]
                v_orf = [p[i][3: 6] for i in range(self.f.n)]
            else:
                r_orf = [p[i*j + 0: i*j + 3] for i in range(self.f.n)]
                v_orf = [p[i*j + 3: i*j + 6] for i in range(self.f.n)]
            q_irf, w_orf = [[None for _ in range(self.f.n)] for _ in range(2)]
        return {'r orf': r_orf, 'v orf': v_orf, 'w brf': w_orf, 'q-3 irf': q_irf}

    def params_dict2vec(self, d: dict, separate_spacecraft: bool = True):
        variables = ['r orf', 'q-3 irf', 'v orf', 'w brf'] if self.o.rotational_motion_navigate else ['r orf', 'v orf']
        if separate_spacecraft:
            return [np.array([d[v][i][j] for v in variables for j in range(3)]) for i in range(self.f.n)]
        else:
            return np.array([d[v][i][j] for i in range(self.f.n) for v in variables for j in range(3)])

    def get_estimation(self, i_f: int, v: str):
        d = self.params_vec2dict()
        return d[v][i_f]

    def calc(self, if_correction: bool) -> None:
        from kiamformation.measurements import measure_antennas_power
        from kiamformation.physics import ode4
        from kiamformation.math_tools import vec2quat
        from kiamformation.spacecrafts import get_gain
        from control import obsv

        # >>>>>>>>>>>> Предварительный расчёт <<<<<<<<<<<<
        if True:
            f = self.f
            c = self.c
            o = self.o
            p = self.p
            j = self.j
            c_len = len(get_gain(o=o, obj=c, r=np.ones(3)))
            f_len = len(get_gain(o=o, obj=f, r=np.ones(3)))
            z_len = int(f.n * c.n * (c_len * f_len) + f.n * (f.n - 1) * f_len**2 // 2)
            my_print(f"Количество измерений:", color='b', end=' ', if_print=p.iter == 1 and o.IF_TEST_PRINT)
            my_print(f"{z_len}", color='g', if_print=p.iter == 1 and o.IF_TEST_PRINT)

        # >>>>>>>>>>>> Этап экстраполяции <<<<<<<<<<<<
        d = self.params_vec2dict()

        # Интегрирование движения -> вектор состояния x_m
        rqvw = [ode4(o=o, obj=f, i=i, t=self.p.t, r=d['r orf'][i], v=d['v orf'][i], q=d['q-3 irf'][i], w=d['w brf'][i])
                for i in range(f.n)]
        x = self.params_dict2vec(d={'r orf': [rqvw[i][0] for i in range(f.n)],
                                    'q-3 irf': [rqvw[i][1] for i in range(f.n)],
                                    'v orf': [rqvw[i][2] for i in range(f.n)],
                                    'w brf': [rqvw[i][3] for i in range(f.n)]}, separate_spacecraft=False)
        d = self.params_vec2dict(params=x, separate_spacecraft=False)
        self.estimation_params_d = d

        # Измерения с поправкой на угловой коэффициент усиления G (signal_rate)
        y = o.MEASURES_VECTOR

        # Измерения согласно модели
        y_model = measure_antennas_power(c=c, f=f, o=o, p=p, j=j, state=x)
        
        tmp = np.abs(y_model - y)
        p.record.loc[p.iter, f'ZModel&RealDifference'] = tmp.mean()
        p.record.loc[p.iter, f'ZModel&RealDifference min'] = tmp.min()
        p.record.loc[p.iter, f'ZModel&RealDifference max'] = tmp.max()
        p.record.loc[p.iter, f'ZModel&RealDifference N'] = len(y_model)
        p.record.loc[p.iter, f'ZReal N'] = len(y)
        p.record.loc[p.iter, f'ZModel N'] = len(y_model)
        for i in range(len(y_model)):
            p.record.loc[p.iter, f'ZModel&RealDifference {i}'] = tmp[i]
            p.record.loc[p.iter, f'ZReal {i}'] = abs(y[i])
            p.record.loc[p.iter, f'ZModel {i}'] = abs(y_model[i])

        # >>>>>>>>>>>> Этап коррекции <<<<<<<<<<<<
        if if_correction:
            from kiamformation.H_matrix import h_matrix

            self.Phi = self.get_Phi(w=None, w0=None, q=None)
            Q_tilda = self.Phi @ self.D @ self.Q @ self.D.T @ self.Phi.T  # * o.dT
            P_m = self.Phi @ self.P @ self.Phi.T + Q_tilda
            q_f = d['q-3 irf'] if o.rotational_motion_navigate else [f.q[i].vec for i in range(f.n)]
            H = h_matrix(t=p.t, o=o, f=f, c=c, r_f=d['r orf'], r_c=c.r_orf, q_f=q_f, q_c=[c.q[i].vec for i in range(c.n)])

            R = np.eye(z_len) * o.kalman_coef['r']
            K = P_m @ H.T @ np.linalg.inv(H @ P_m @ H.T + R)
            self.P = (np.eye(j * f.n) - K @ H) @ P_m
            raw_estimation_params = np.array(np.matrix(x) + K @ (y - y_model))[0]

            # Численный расчёт STM
            self.STM = self.Phi if self.STM is None else self.Phi @ self.STM
            tmp = self.STM.T @ H.T @ H @ self.STM
            self.observability_gramian = tmp if self.observability_gramian is None else self.observability_gramian + tmp
            _, simgas, _ = np.linalg.svd(self.observability_gramian)
            p.record.loc[p.iter, f'gramian sigma criteria'] = np.min(simgas)/np.max(simgas)

            tmp = obsv((self.Phi - np.eye(self.Phi.shape[0])) / self.o.dt, H)
            _, simgas, _ = np.linalg.svd(tmp)
            p.record.loc[p.iter, f'linear rank criteria'] = np.linalg.matrix_rank(tmp)
            p.record.loc[p.iter, f'linear sigma criteria'] = np.min(simgas)/np.max(simgas)
        else:
            raw_estimation_params = x

        # >>>>>>>>>>>> Обновление оценки <<<<<<<<<<<<
        for i in range(f.n):
            tmp = raw_estimation_params[(0 + i) * j: (1 + i) * j]

            # Нормировка кватерниона
            if o.rotational_motion_navigate:
                tmp[3:6] = vec2quat(tmp[3:6]).normalized().vec

            # Запись
            self.estimation_params[i] = tmp


def navigate(k: KalmanFilter, if_correction: bool = True):
    k.calc(if_correction=if_correction)  # Пока что при любом OPERATING_MODES (если все КА включены и не потеряны)
