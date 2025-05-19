"""2D and 3D display of numerical modelling results

Attributes:
    FEMTO_RATE : Variation of 3D FemtoSat dimensions
    CUBE_RATE : Variation of 3D CubeSat dimensions
    TITLE_SIZE : Size of 2D plot titles
    CAPTION_SIZE : Size of 2D plot x- and y-labels
"""
import numpy as np
from kiamformation.config import Objects
from kiamformation.static.cosmetic import my_print


FEMTO_RATE = 1
CUBE_RATE = 1/3
TITLE_SIZE = 15
CAPTION_SIZE = 13

# >>>>>>>>>>>> 2D графики <<<<<<<<<<<<
def plot_observability_criteria(o):
    global TITLE_SIZE, CAPTION_SIZE
    import matplotlib.pyplot as plt
    if 'linear sigma criteria' in o.p.record.keys():
        x = o.p.record['t'].to_list()
        label_time = {"рус": f"Время, с", "eng": f"Time, s"}[o.LANGUAGE]
        labels = {"рус": ["", "", ""],
                  "eng": ["Linear singular values ratio σᴸ", "Linear rank",
                          "Gramian singular values ratio σᵂ"]}[o.LANGUAGE]
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))  # , gridspec_kw={'height_ratios': [3, 1]})
        for i, s in enumerate(['linear sigma criteria']):  # , 'linear rank criteria'
            ax.plot(x, o.p.record[s].to_list(), c='blue', label=labels[i])
            ax.grid(True)
            ax.set_xlabel(label_time, fontsize=CAPTION_SIZE)

        ax2 = ax.twinx()
        ax.set_ylabel(labels[0], fontsize=CAPTION_SIZE).set_color('blue')
        ax2.set_ylabel(labels[2], fontsize=CAPTION_SIZE).set_color('red')
        ax2.plot(x, o.p.record['gramian sigma criteria'].to_list(), c='red')

        ax.set_yscale('log')
        ax2.set_yscale('log')
        plt.show()

        '''
        for p, param in enumerate(['qλ', 'wω']):
            fig, ax = plt.subplots(1, 1, figsize=(8, 5))
            for i, j in zip('xyz', 'ˣʸᶻ'):
                ax.plot(x, o.p.record[f"{o.f.name} {param[0]} {i} {'ib'[p]}rf 0"].to_list(), label=f"{param[1]}{j}")
            ax.legend(fontsize=CAPTION_SIZE)
            ax.grid(True)
            ax.set_xlabel(label_time, fontsize=CAPTION_SIZE)
            ax.set_ylabel(["Quaternion components", "Angular velocity components"][p], fontsize=CAPTION_SIZE)
            plt.show()
        # '''


def plot_distance(o):
    global TITLE_SIZE, CAPTION_SIZE
    import matplotlib.pyplot as plt

    plot_observability_criteria(o)

    fig, ax = plt.subplots(2 if o.rotational_motion_navigate else 1,
                           2 if o.rotational_motion_navigate else 2,
                           figsize=(12, 20 if o.rotational_motion_navigate else 5))
    axes = ax[0] if o.rotational_motion_navigate else ax
    title = {"рус": f"Неточности в навигации", "eng": f"Navigation Errors"}[o.LANGUAGE]
    label_time = {"рус": f"Время, с", "eng": f"Time, s"}[o.LANGUAGE]
    fig.suptitle(title, fontsize=TITLE_SIZE)

    x = o.p.record['t'].to_list()
    for i_c in range(o.c.n):
        for i_f in range(o.f.n):
            labels = {"рус": ["Разница измерений модельных и полученных",
                              "Ошибка определения положения Δr"],
                      "eng": ["Residuals yᵖʳᵉᵈⁱᶜᵗ-y",
                              "Real measurement y",
                              "Model measurement yᵖʳᵉᵈⁱᶜᵗ"]}[o.LANGUAGE] \
                if i_f == 0 else [None for _ in range(100)]
            for jj in range(int(o.p.record[f'ZModel&RealDifference N'][1])):
                y2 = o.p.record[f'ZModel&RealDifference {jj}'].to_list()
                # axes[0].plot(x, y2, c='blue', label=labels[0] if i_c == 0 and jj == 0 else None, lw=1)
                axes[0].plot(x, y2, c='blue', label=labels[0] if i_c == 0 and jj == 0 else None, lw=1)
            '''for jj in range(int(o.p.record[f'ZReal N'][1])):
                y1 = o.p.record[f'ZReal {jj}'].to_list()
                axes[0].plot(x, y1, label=labels[1] if i_c == 0 and jj == 0 else None, lw=1, ls="-")
            for jj in range(int(o.p.record[f'ZModel N'][1])):
                y1 = o.p.record[f'ZModel {jj}'].to_list()
                axes[0].plot(x, y1, label=labels[2] if i_c == 0 and jj == 0 else None, lw=1, ls="-")'''
    axes[0].set_xlabel(label_time, fontsize=CAPTION_SIZE)
    axes[0].set_ylabel({"рус": f"Невязка", "eng": f"Residuals, 1/m²"}[o.LANGUAGE], fontsize=CAPTION_SIZE)
    axes[0].legend(fontsize=CAPTION_SIZE)
    axes[0].grid(True)

    for i_f in range(o.f.n):
        labels = ["ΔX", "ΔY", "ΔZ", "Error of estimated position Δr"]
        for j, c in enumerate('xyz'):
            y = o.p.record[f'{o.f.name} KalmanPosError {c} {i_f}'].to_list()
            axes[1].plot(x, y, label=labels[j] if i_f == 0 else None)
        # y3 = o.p.record[f'{o.f.name} KalmanPosError r {i_f}'].to_list()
        # axes[1].plot(x, y3, c='k', label=labels[3] if i_f == 0 else None)
    axes[1].set_xlabel(label_time, fontsize=CAPTION_SIZE)
    axes[1].set_ylabel({"рус": f"Δr компоненты, м", "eng": f"Estimation error, m"}[o.LANGUAGE], fontsize=CAPTION_SIZE)
    axes[1].legend(fontsize=CAPTION_SIZE)
    axes[1].grid(True)

    if o.rotational_motion_navigate:
        for i_f in range(o.f.n):
            labels_dq = ["Δλˣ", "Δλʸ", "Δλᶻ"]
            labels_dw = ["Δωˣ", "Δωʸ", "Δωᶻ"]
            for j, c in enumerate('xyz'):
                y1 = o.p.record[f'{o.f.name} KalmanQuatError {c} {i_f}'].to_list()
                y2 = o.p.record[f'{o.f.name} KalmanSpinError ORF {c} {i_f}'].to_list()
                ax[1][0].plot(x, y1, label=labels_dq[j] if i_f == 0 else None)
                ax[1][1].plot(x, y2, label=labels_dw[j] if i_f == 0 else None)
        for ii in [1]:
            ax[ii][0].set_ylabel({"рус": "Ошибки λ",
                                  "eng": "Quaternion components errors"}[o.LANGUAGE], fontsize=CAPTION_SIZE)
            ax[ii][1].set_ylabel({"рус": "Ошибки ω (ORF)",
                                  "eng": "Angular velocity errors, rad/s²"}[o.LANGUAGE], fontsize=CAPTION_SIZE)
            for j in range(2):
                ax[ii][j].legend(fontsize=CAPTION_SIZE)
                ax[ii][j].set_xlabel(label_time, fontsize=CAPTION_SIZE)
                ax[ii][j].grid(True)
    plt.show()

    '''fig, ax = plt.subplots()
    for i in range(int(o.p.record[f'G N'][1])):
        ax.plot(x, o.p.record[f'G {i}'].to_list())
    plt.title("gains")
    plt.show()'''

def plot_atmosphere_models(n: int = 100):
    import matplotlib.pyplot as plt
    from kiamformation.physics import get_atm_params

    o = Objects()
    range_km = [300, 500]
    fig, axs = plt.subplots(1, 2, figsize=(19, 5))
    fig.suptitle('Модели атмосферы')
    for m in range(len(o.ATMOSPHERE_MODELS)):
        for j in range(2):
            z = np.linspace(100e3 if j == 0 else range_km[0]*1e3, range_km[1]*1e3, n)
            rho = [get_atm_params(h=z[i], o=o, atm_model=o.ATMOSPHERE_MODELS[m])[0] for i in range(n)]
            tmp = ", используемая" if o.ATMOSPHERE_MODELS[m] == o.ATMOSPHERE_MODEL else ""
            axs[j].plot(z, rho, ls="-" if o.ATMOSPHERE_MODELS[m] == o.ATMOSPHERE_MODEL else ":",
                        label=f"Модель: {o.ATMOSPHERE_MODELS[m]}{tmp}")
            axs[j].set_ylabel(f"Плотность ρ, кг/м³")
            axs[j].set_xlabel(f"Высота z, м")
            axs[j].legend()
            axs[j].grid()
    axs[0].set_title(f"От линии Кармана до {range_km[1]} км")
    axs[1].set_title(f"От {range_km[0]} до {range_km[1]} км")
    plt.show()


# >>>>>>>>>>>> 3D отображение в ОСК <<<<<<<<<<<<
def show_chipsat(o, j, clr, opacity, reference_frame: str, return_go: bool = True, ax=None, xyz=None) -> list:
    """Функция отображения дочернего КА (квадратная пластина)
    :param o: Objects
    :param j: номер чипсата
    :param clr:
    :param opacity:
    :param reference_frame:
    :param return_go: По умолчанию, при return_go=False считается, что reference_frame="BRF"
    :param ax:
    :param xyz: Размеры КА
    :return:
    """
    import plotly.graph_objs as go
    from kiamformation.physics import get_matrices
    global FEMTO_RATE

    rate = FEMTO_RATE if reference_frame != "BRF" else 2 / min(o.f.size)
    x, y, z = ([], [], [0 for _ in range(4)])
    for x_shift in [-o.f.size[0] * rate, o.f.size[0] * rate]:
        for y_shift in [-o.f.size[1] * rate, o.f.size[1] * rate]:
            x += [x_shift]
            y += [y_shift]
    U, S, A, _ = get_matrices(o=o, t=o.p.t, obj=o.f, n=j)

    if reference_frame != "BRF":
        for i in range(4):
            r = S.T @ np.array([x[i], y[i], z[i]])
            if reference_frame == 'ORF' and xyz is not None:
                r[0] *= xyz[0] if o.RELATIVE_SIDES else max(xyz)
                r[1] *= xyz[1] if o.RELATIVE_SIDES else max(xyz)
                r[2] *= xyz[2] if o.RELATIVE_SIDES else max(xyz)
            x[i] = r[0] + o.p.record[f'{o.f.name} r x {reference_frame.lower()} {j}'].to_list()[-1]
            y[i] = r[1] + o.p.record[f'{o.f.name} r y {reference_frame.lower()} {j}'].to_list()[-1]
            z[i] = r[2] + o.p.record[f'{o.f.name} r z {reference_frame.lower()} {j}'].to_list()[-1]
        r_real, r_estimation = [], []
        for c in "xyz":
            r_real.append(o.p.record[f'{o.f.name} r {c} {reference_frame.lower()} {j}'].to_list())
            r_estimation.append(o.p.record[f'{o.f.name} KalmanPosEstimation {c} {j}'].to_list())
        if return_go:
            if j == 0:
                name_1 = 'Modelled ChipSat trajectory' if o.f.n == 1 else 'Modelled ChipSat trajectories'
                name_2 = 'Estimated ChipSat trajectory' if o.f.n == 1 else 'Estimated ChipSat trajectories'
            else:
                name_1 = name_2 = None
            return [go.Mesh3d(x=x, y=y, z=z, color=clr, opacity=opacity),
                    go.Scatter3d(x=r_real[0], y=r_real[1], z=r_real[2], mode='lines',
                                 line=dict(color='blue', width=5), name=name_1,
                                 legend="legend2" if name_1 is None else "legend"),
                    go.Scatter3d(x=r_estimation[0], y=r_estimation[1], z=r_estimation[2], mode='lines',
                                 line=dict(color='red', width=5), name=name_2, 
                                 legend="legend2" if name_2 is None else "legend")]
    ax.plot([[c[i] for i in [0, 2, 3, 1, 0]] for c in [x, y, z]], c='gray', linewidth=3)


def show_cubesat(o, j, reference_frame: str, xyz=None) -> list:
    import plotly.graph_objs as go
    from kiamformation.physics import get_matrices
    global CUBE_RATE

    n_legs = 4
    total_cubes = 6 * (n_legs + 1)
    r = [[[] for _ in range(total_cubes)] for _ in range(3)]
    # Памятка: r[x/y/z][0..5 - yellow, 6..29 - gray][1..4 - sides of square]
    sequence = [[0, 0], [0, 1], [1, 1], [1, 0]]

    shift = [[-o.c.size[i] * CUBE_RATE, o.c.size[i] * CUBE_RATE] for i in range(3)]
    legs = [o.c.legs_x, o.c.legs_x, o.c.legs_z]
    bound_legs = [[[((-1)**sequence[s][0] * o.c.size[0] - legs[0]) * CUBE_RATE,
                    ((-1)**sequence[s][0] * o.c.size[0] + legs[0]) * CUBE_RATE],
                   [((-1)**sequence[s][1] * o.c.size[1] - legs[1]) * CUBE_RATE,
                    ((-1)**sequence[s][1] * o.c.size[1] + legs[1]) * CUBE_RATE],
                   [(-o.c.size[2] - legs[2]) * CUBE_RATE, (o.c.size[2] + legs[2]) * CUBE_RATE]] for s in range(n_legs)]

    for i in range(3):
        for k in range(2):
            r[i][k + i * 2] = [(-1)**(k+1) * o.c.size[i] * CUBE_RATE for _ in range(4)]
            tmp = k + i * 2
            ind1 = 0 + int(i < 1)
            ind2 = 1 + int(i < 2)
            r[ind1][tmp].extend([shift[ind1][sequence[m][0]] for m in range(4)])
            r[ind2][tmp].extend([shift[ind2][sequence[m][1]] for m in range(4)])

    for s in range(n_legs):
        r[0][(s + 1) * 6 + 0] = [bound_legs[s][0][0], bound_legs[s][0][0], bound_legs[s][0][0], bound_legs[s][0][0]]
        r[1][(s + 1) * 6 + 0] = [bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 0] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 1] = [bound_legs[s][0][1], bound_legs[s][0][1], bound_legs[s][0][1], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 1] = [bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 1] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 2] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 2] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 2] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 3] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 3] = [bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 3] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 4] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 4] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 4] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][0]]

        r[0][(s + 1) * 6 + 5] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 5] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 5] = [bound_legs[s][2][1], bound_legs[s][2][1], bound_legs[s][2][1], bound_legs[s][2][1]]

    # Костыль: не отображаются ровно вертикальные грани
    U, S, A, _ = get_matrices(o=o, t=o.p.t, obj=o.c, n=j, q=o.c.q[j] + np.quaternion(0.01, 0, 0, 0))
    # S = quart2dcm(np.quaternion(1/2+0.01, -1/2, -1/2, -1/2).normalized())
    for k in range(total_cubes):
        for i in range(4):
            r1 = S.T @ np.array([r[0][k][i], r[1][k][i], r[2][k][i]])
            if reference_frame == 'ORF' and xyz is not None:
                r1[0] *= xyz[0] if o.RELATIVE_SIDES else max(xyz)
                r1[1] *= xyz[1] if o.RELATIVE_SIDES else max(xyz)
                r1[2] *= xyz[2] if o.RELATIVE_SIDES else max(xyz)
            r[0][k][i] = r1[0] + o.p.record[f'{o.c.name} r x {reference_frame.lower()} {j}'].to_list()[-1]
            r[1][k][i] = r1[1] + o.p.record[f'{o.c.name} r y {reference_frame.lower()} {j}'].to_list()[-1]
            r[2][k][i] = r1[2] + o.p.record[f'{o.c.name} r z {reference_frame.lower()} {j}'].to_list()[-1]
    anw = []
    for i in range(total_cubes):
        color = 'yellow' if i < 6 else 'gray'
        anw.append(go.Mesh3d(x=r[0][i], y=r[1][i], z=r[2][i], color=color, opacity=1))
    r = []
    for c in 'xyz':
        r.append(o.p.record[f'{o.c.name} r {c} {reference_frame.lower()} {j}'].to_list())
    
    name_1 = 'CubeSat trajectory' if o.c.n == 1 else ('CubeSat trajectories' if j == 0 else None)
    anw.append(go.Scatter3d(x=r[0], y=r[1], z=r[2], mode='lines', line=dict(color='tan', width=5), 
                            name=name_1, legend="legend2" if name_1 is None else "legend"))
    return anw


def plot_model_gain(o: Objects, n: int = 20):
    import matplotlib.pyplot as plt
    from kiamformation.math_tools import pol2dec
    from kiamformation.spacecrafts import get_gain

    fig = plt.figure(figsize=(15, 10))

    for i in range(2):
        for j in range(2):
            ax = fig.add_subplot(2, 2, i + 1 + 2*j, projection='3d')
            obj = [o.f, o.c][i]
            my_print(f"Диаграмма направленностей для {obj.name}: {o.GAIN_MODEL_F}", color="b")

            u = np.linspace(0, 2 * np.pi, n)
            v = np.linspace(-np.pi / 2, np.pi / 2, n)
            U, V = np.meshgrid(u, v)

            max_g = 0
            for k in range(len(get_gain(o=o, obj=obj, r=pol2dec(1, u[0], v[0])))):
                g = np.array([[get_gain(o=o, obj=obj, r=pol2dec(1, u[ii], v[jj]))[k]
                               for ii in range(n)] for jj in range(n)])
                X, Y, Z = pol2dec(g, U, V)
                ax.plot_surface(X, Y, Z, cmap=o.MY_COLORMAPS[k])
                max_g = max(max_g, np.max(g.flatten()))

            # if obj.name == "FemtoSat":
            #     max_g = max(max_g, 2 * np.max(o.f.size))
            #     show_chipsat(o=o, j=0, reference_frame="BRF", return_go=False, ax=ax, clr=None, opacity=None)

            ax.set_xlim(-max_g, max_g)
            ax.set_ylim(-max_g, max_g)
            ax.set_zlim(-max_g, max_g)
            ax.set_box_aspect([1, 1, 1])
            title_text = f"Диаграмма направленностей для {obj.name} | " \
                         f"GAIN_MODEL = {o.GAIN_MODEL_F if  obj.name == 'FemtoSat' else o.GAIN_MODEL_C}\n" \
                         f"Отправленный сигнал"
            ax.set_title(title_text if j == 0 else "Принятый сигнал")
    plt.show()


# >>>>>>>>>>>> 3D отображение в ИСК <<<<<<<<<<<<
def arrows3d(ends, starts=None, ax=None, label: str = None, **kwargs):
    """Построение 3D стрелок
    GitHub: https://github.com/matplotlib/matplotlib/issues/22571
    :param ends: (N, 3) size array of arrow end coordinates
    :param starts: (N, 3) size array of arrow start coordinates.
    :param ax: (Axes3DSubplot) existing axes to add to
    :param label: legend label to apply to this group of arrows
    :param kwargs: additional arrow properties"""
    from matplotlib.patches import FancyArrowPatch
    import matplotlib.pyplot as plt

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            super().__init__((0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def do_3d_projection(self, renderer=None):
            from mpl_toolkits.mplot3d import proj3d

            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

            return np.min(zs)

    if starts is None:
        starts = np.zeros_like(ends)

    assert starts.shape == ends.shape, "`starts` and `ends` shape must match"
    assert len(ends.shape) == 2 and ends.shape[1] == 3, \
        "`starts` and `ends` must be shape (N, 3)"


    # create new axes if none given
    if ax is None:
        ax = plt.figure().add_subplot(111, projection='3d')
    arrow_prop_dict = dict(mutation_scale=20, arrowstyle='-|>', color='k', shrinkA=0, shrinkB=0)
    arrow_prop_dict.update(kwargs)
    for ind, (s, e) in enumerate(np.stack((starts, ends), axis=1)):
        a = Arrow3D(
            [s[0], e[0]], [s[1], e[1]], [s[2], e[2]],
            # only give label to first arrow
            label=label if ind == 0 else None,
            **arrow_prop_dict)
        ax.add_artist(a)
    ax.points = np.vstack((starts, ends, getattr(ax, 'points', np.empty((0, 3)))))
    return ax


def plot_the_earth_mpl(ax, o: Objects, res: int = 1, pth: str = "./", earth_image=None):
    """Отрисовка слева красивой Земли из одной линии"""
    from PIL import Image

    x_points = 180 * res
    y_points = 90 * res

    if earth_image is None:
        bm = Image.open(f'{pth}source/skins/{o.EARTH_FILE_NAME}')
        bm = np.array(bm.resize((x_points, y_points))) / 256.
    else:
        bm = earth_image

    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi / 180
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi / 180

    x = o.EARTH_RADIUS * np.outer(np.cos(lons), np.cos(lats)).T
    y = o.EARTH_RADIUS * np.outer(np.sin(lons), np.cos(lats)).T
    z = o.EARTH_RADIUS * np.outer(np.ones(np.size(lons)), np.sin(lats)).T
    ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=bm, alpha=1)
    ax.set_xlabel("x, тыс. км")
    ax.set_ylabel("y, тыс. км")
    ax.set_zlabel("z, тыс. км")
    return ax


def plot_the_earth_go(o: Objects):
    import plotly.graph_objs as go

    spherical_earth_map = np.load('kiamformation/static/map_sphere.npy')
    xm, ym, zm = spherical_earth_map.T * o.EARTH_RADIUS
    return go.Scatter3d(x=xm, y=ym, z=zm, mode='lines')


def plot_reference_frames(ax, o, txt: str, color: str = "gray", t: float = None):
    from kiamformation.physics import get_matrices

    x, y, z = np.eye(3)
    arrows, start = np.array([x, y, z]) * o.ORBIT_RADIUS, None
    if txt == "ОСК":
        U, _, _, R_orb = get_matrices(obj=o.a, n=0, t=t, o=o)
        arrows = np.array([U.T @ x, U.T @ y, U.T @ z]) * o.ORBIT_RADIUS / 2
        start = R_orb
    if txt == "ИСК":
        # Отрисовка кружочка
        x, y, z = [o.p.record[f'Anchor r {c} irf 0'].to_list() for c in 'xyz']
        ax.plot(x, y, z, color)
        start = np.zeros(3)
    ax = arrows3d(starts=np.array([start for _ in range(3)]), ends=np.array([start + arrows[i] for i in range(3)]),
                  ax=ax, color=color, label=txt)
    for i in range(3):
        label = ["x, m", "y, m", "z, m"][i]
        a = start + arrows[i] + arrows[i] / np.linalg.norm(arrows[i]) * 0.2
        ax.text(a[0], a[1], a[2], c=color, s=label)
    return ax


# >>>>>>>>>>>> Анимация <<<<<<<<<<<<
def animate_reference_frames(resolution: int = 3, n: int = 50):
    import matplotlib.pyplot as plt
    from PIL import Image
    from os import remove

    o = Objects()
    # o.dT = o.SEC_IN_TURN / (n - 3)
    TIME = 2*np.pi / o.W_ORB
    o.dt = TIME / (n - 3)
    o.if_navigation = False

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim3d([-1e7, 1e7])
    ax.set_ylim3d([-1e7, 1e7])
    ax.set_zlim3d([-1e7, 1e7])
    x_points = 180 * resolution
    y_points = 90 * resolution
    earth_image = Image.open(f'kiamformation/static/{o.EARTH_FILE_NAME}')
    earth_image = np.array(earth_image.resize((x_points, y_points))) / 256.
    for i in range(n):
        o.p.time_step()
        ax = plot_the_earth_mpl(ax, o=o, earth_image=earth_image)
        ax = plot_reference_frames(ax, o, t=o.p.t, txt="ИСК", color="lime")
        ax = plot_reference_frames(ax, o, t=o.p.t, txt="ОСК", color="red")
        ax.view_init(azim=20, elev=30, roll=0)
        ax.axis('equal')
        plt.title(f"Наклонение: {o.INCLINATION}°, эксцентриситет: {round(o.ECCENTRICITY, 2)}, "
                  f"апогей: {round(o.APOGEE / 1e3)} км, перигей: {round(o.PERIGEE / 1e3)} км")
        plt.legend()
        plt.savefig(f"kiamformation/temporary/to_delete_{'{:04}'.format(i)}.jpg")
        ax.clear()
    plt.close()

    images = [Image.open(f"kiamformation/temporary/to_delete_{'{:04}'.format(i)}.jpg") for i in range(n)]
    images[0].save('kiamformation/temporary/res.gif', save_all=True, append_images=images[1:], duration=20, loop=0)
    for i in range(n):
        remove(f"kiamformation/temporary/to_delete_{'{:04}'.format(i)}.jpg")


# >>>>>>>>>>>> Функции отображения, готовые к использованию <<<<<<<<<<<<
def show_chipsats_and_cubesats(o, reference_frame: str, clr: str = 'lightpink', opacity: float = 1):
    xyz_max, xyz_min = [-1e4] * 3, [1e4] * 3
    xyz = np.zeros(3)
    for j, c in enumerate("xyz"):
        for i in range(o.f.n):
            for l in [o.p.record[f'{o.f.name} r {c} {reference_frame.lower()} {i}'].to_list(),
                      o.p.record[f'{o.f.name} KalmanPosEstimation {c} {i}'].to_list()]:
                xyz_max[j] = max(xyz_max[j], max(l))
                xyz_min[j] = min(xyz_min[j], min(l))
        for i in range(o.c.n):
            l = o.p.record[f'{o.c.name} r {c} {reference_frame.lower()} {i}'].to_list()
            xyz_max[j] = max(xyz_max[j], max(l))
            xyz_min[j] = min(xyz_min[j], min(l))
        xyz[j] = xyz_max[j] - xyz_min[j]
    data = []
    for i in range(o.f.n):
        data += show_chipsat(o, i, clr, opacity, reference_frame, xyz=xyz)
    for i in range(o.c.n):
        data += show_cubesat(o, i, reference_frame, xyz=xyz)
    return data


def plot_all(o, save: bool = False, count: int = None):
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'surface'}, ]], # {'type': 'surface'}
                        subplot_titles=('Орбитальная СК',))  # 'Инерциальная СК'
    for i in range(1):
        tmp = show_chipsats_and_cubesats(o, ['ORF', 'IRF'][i], clr='black')
        for surf in tmp:
            fig.add_trace(surf, row=1, col=i+1)
    # fig.add_trace(plot_the_earth_go(v=o.v), row=1, col=1)
    fig.update_layout(title=dict(text=f"Солвер: {o.SOLVER} "
                                      f"{'(' if o.physics_model['aero drag'] or o.physics_model['j2'] else ''}"
                                      f"{' +Лобовое сопротивление' if o.physics_model['aero drag'] else ''}"
                                      f"{' +Вторая гармоника' if o.physics_model['j2'] else ''}"
                                      f"{' )' if o.physics_model['aero drag'] or o.physics_model['j2'] else ''}"
                                      f" | Время: {o.time} ({round(o.time / (3600 * 24), 2)} дней)  |  "
                                      f"i={o.INCLINATION}°, e={o.ECCENTRICITY}"),
                      legend=dict(yanchor='top', xanchor='left', y=0.7, x=0.3))
    fig.update_scenes(xaxis_title_text='x, m',
                      yaxis_title_text='y, m',
                      zaxis_title_text='z, m')
    if save:
        fig.write_image('../../img/' + str('{:04}'.format(count)) + '.jpg')
    else:
        fig.show()
