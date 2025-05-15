"""
Файл сгененрирован программой OpenCompendium/RnD-PhD/notebooks/H matrix generation.ipynb
"""

def h_element(gm_1: str, gm_2: str, fn, angles_navigation, relation, r1, r2, w_0: float, t: float, q1=None, q2=None):
    """
    :param fn: Количество дочерних КА
    :param angles_navigation: Оценивается ли вращательное движение
    :param r1: Положение 1-го КА
    :param r2: Положение 2-го КА
    :param w_0: Угловая скорость вращения ОСК относительно ИСК
    :param t: Текущее время
    :param q1: Вектор-часть кватерниона 1-го КА (при angles_navigation=True)
    :param q2: Вектор-часть кватерниона 2-го КА (при angles_navigation=True)
    relation : str
        Chief-deputy (cd) or deputy-deputy (dd)
    """
    from kiamformation.flexmath import pi, setvectype, sqrt, cos, sin

    ff_sequence = []  # Последовательность номеров непустых столбцов, длина ff_sequence - кол-во строк нижней подматицы
    for i_f1 in range(fn):
        for i_f2 in range(i_f1):
            if i_f1 != i_f2:
                ff_sequence += [[i_f1, i_f2]]

    r1_x, r1_y, r1_z = r1
    r2_x, r2_y, r2_z = r2
    r12x = r1_x - r2_x
    r12y = r1_y - r2_y
    r12z = r1_z - r2_z
    r12 = sqrt(r12x**2 + r12y**2 + r12z**2)
    if q1 is not None:
        q1_x, q1_y, q1_z = q1
        q1_0 = sqrt(1 - q1_x**2 - q1_y**2 - q1_z**2)
    if q2 is not None:
        q2_x, q2_y, q2_z = q2
        q2_0 = sqrt(1 - q2_x**2 - q2_y**2 - q2_z**2)
    pi = pi(r1_x)

    <to_replace>
    
    # if angles_navigation:
    swt = sin(t*w_0)
    cwt = cos(t*w_0)
    s1_x_r12 = s1_11*(r12x) + s1_12*(r12y) + s1_13*(r12z)
    s1_y_r12 = s1_21*(r12x) + s1_22*(r12y) + s1_23*(r12z)
    s1_z_r12 = s1_31*(r12x) + s1_32*(r12y) + s1_33*(r12z)
    s2_x_r12 = s2_11*(r12x) + s2_12*(r12y) + s2_13*(r12z)
    s2_y_r12 = s2_21*(r12x) + s2_22*(r12y) + s2_23*(r12z)
    s2_z_r12 = s2_31*(r12x) + s2_32*(r12y) + s2_33*(r12z)
    s1_r12 = sqrt((s1_x_r12)**2 + (s1_y_r12)**2 + (s1_z_r12)**2)
    s2_r12 = sqrt((s2_x_r12)**2 + (s2_y_r12)**2 + (s2_z_r12)**2)
    s1_r12_2 = ((s1_x_r12)**2 + (s1_y_r12)**2 + (s1_z_r12)**2)
    s2_r12_2 = ((s2_x_r12)**2 + (s2_y_r12)**2 + (s2_z_r12)**2)

    s1_xyx = sqrt((s1_x_r12)**2 + (s1_y_r12)**2 + (s1_x_r12)**2)
    s1_xy = sqrt((s1_x_r12)**2 + (s1_y_r12)**2)
    s1_yz = sqrt((s1_y_r12)**2 + (s1_z_r12)**2)
    s1_xz = sqrt((s1_x_r12)**2 + (s1_z_r12)**2)
    s2_xyx = sqrt((s2_x_r12)**2 + (s2_y_r12)**2 + (s2_x_r12)**2)
    s2_xy = sqrt((s2_x_r12)**2 + (s2_y_r12)**2)
    s2_yz = sqrt((s2_y_r12)**2 + (s2_z_r12)**2)
    s2_xz = sqrt((s2_x_r12)**2 + (s2_z_r12)**2)

    s1_xyx_2 = ((s1_x_r12)**2 + (s1_y_r12)**2 + (s1_x_r12)**2)
    s1_xy_2 = ((s1_x_r12)**2 + (s1_y_r12)**2)
    s1_yz_2 = ((s1_y_r12)**2 + (s1_z_r12)**2)
    s1_xz_2 = ((s1_x_r12)**2 + (s1_z_r12)**2)
    s2_xyx_2 = ((s2_x_r12)**2 + (s2_y_r12)**2 + (s2_x_r12)**2)
    s2_xy_2 = ((s2_x_r12)**2 + (s2_y_r12)**2)
    s2_yz_2 = ((s2_y_r12)**2 + (s2_z_r12)**2)
    s2_xz_2 = ((s2_x_r12)**2 + (s2_z_r12)**2)

    c_s1_x = cos(pi*(s1_x_r12)/(2*s1_r12))
    c_s1_y = cos(pi*(s1_y_r12)/(2*s1_r12))
    c_s1_z = cos(pi*(s1_z_r12)/(2*s1_r12))
    c_s2_x = cos(pi*(s2_x_r12)/(2*s2_r12))
    c_s2_y = cos(pi*(s2_y_r12)/(2*s2_r12))
    c_s2_z = cos(pi*(s2_z_r12)/(2*s2_r12))

    q1_y_z2 = -2*q1_y**2 - 2*q1_z**2 + 1
    q1_x_z2 = -2*q1_x**2 - 2*q1_z**2 + 1
    q1_0z_xy1 = 2*q1_0*q1_z + 2*q1_x*q1_y
    q1_0z_xy2 = -2*q1_0*q1_z + 2*q1_x*q1_y
    q1_0y_xz1 = 2*q1_0*q1_y + 2*q1_x*q1_z
    q1_0x_yz2 = -2*q1_0*q1_x + 2*q1_y*q1_z

    q2_y_z2 = -2*q2_y**2 - 2*q2_z**2 + 1
    q2_x_z2 = -2*q2_x**2 - 2*q2_z**2 + 1
    q2_0z_xy1 = 2*q2_0*q2_z + 2*q2_x*q2_y
    q2_0z_xy2 = -2*q2_0*q2_z + 2*q2_x*q2_y
    q2_0y_xz1 = 2*q2_0*q2_y + 2*q2_x*q2_z
    q2_0x_yz2 = -2*q2_0*q2_x + 2*q2_y*q2_z

    i001 = 2*q1_x*q1_z + 2*q1_y*q1_0
    i002 = 2*q2_x*q2_z + 2*q2_y*q2_0
    i003 = 2*q1_x*q1_y + 2*q1_z*q1_0
    i004 = 2*q2_x*q2_y + 2*q2_z*q2_0
    i005 = 2*q1_y*q1_z + 2*q1_x*q1_0
    i006 = 2*q1_z**2/q1_0 - 2*q1_0
    i007 = 2*q2_z**2/q2_0 - 2*q2_0
    i008 = 2*q1_x*q1_z/q1_0 + 2*q1_y
    i009 = 2*q2_x*q2_z/q2_0 + 2*q2_y
    i010 = 2*q1_x*q1_y/q1_0 + 2*q1_z
    i011 = 2*q2_x*q2_y/q2_0 + 2*q2_z
    i012 = 2*q1_x**2/q1_0 - 2*q1_0
    i013 = 2*q2_x**2/q2_0 - 2*q2_0
    i014 = -2*q1_x*q1_0 + 2*q1_y*q1_z
    i015 = -2*q2_x*q2_0 + 2*q2_y*q2_z
    i016 = -2*q2_x*q2_z/q2_0 + 2*q2_y
    i017 = 2*q2_x + 2*q2_y*q2_z/q2_0
    i018 = 2*q2_x - 2*q2_y*q2_z/q2_0
    i019 = 2*q2_x*q2_y - 2*q2_z*q2_0
    i020 = -4.0*cwt*q2_x - 1.0*swt*(i009)
    i021 = 1.0*cwt*(i009) - 4.0*q2_x*swt
    i022 = -2.0*q2_x**2/q2_0 + 2.0*q2_0
    i023 = -2*q2_y**2/q2_0 + 2*q2_0
    i024 = -8.0*q2_x*(r12y) + 2*(r12x)*(1.0*cwt*(i013) - 1.0*swt*(-i011))
    i025 = (-2.0*cwt*(i002) - 2.0*swt*(i015))*(s2_z_r12)
    i026 = (-2.0*cwt*(q2_x_z2) + 2.0*swt*(i019))*(s2_y_r12)
    i027 = (-4.0*q2_x*q2_0 - 4.0*q2_y*q2_z)*(s2_y_r12)
    i028 = (s2_y_r12)**2/s2_r12_2 + (s2_z_r12)**2/s2_r12_2
    i029 = (s1_y_r12)**2/s1_xyx_2 + (s1_x_r12)**2/s1_xyx_2
    i030 = (r12x)**2 + (r12y)**2 + (r12z)**2
