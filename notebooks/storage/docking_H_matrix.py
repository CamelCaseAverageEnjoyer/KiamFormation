

def h_matrix(t, o, f, c, r_f, r_c, q_f, q_c: list, return_template: bool = False):
    """Возвращает матрицу отображения наблюдаемости H = dz/dx. Только измерения RSSI""" 
    from sympy import var, Matrix
    from kiamformation.spacecrafts import get_gain
    from kiamformation.flexmath import zeros, block_diag, vstack, bmat
    import numpy as np

    fn = f.n
    cn = c.n
    c_g = o.GAIN_MODEL_C
    f_g = o.GAIN_MODEL_F
    angles_navigation = o.rotational_motion_navigate
    w_0 = o.W_ORB

    f_take_len = f_send_len = len(get_gain(o=o, obj=f, r=np.ones(3)))
    
    ff_sequence = []  # Последовательность номеров непустых столбцов, длина ff_sequence - кол-во строк нижней подматицы
    for i_f1 in range(fn):
        for i_f2 in range(i_f1):
            if i_f1 != i_f2:
                ff_sequence += [[i_f1, i_f2]]

    H_cd = None
    for i_c in range(cn):
        row = []
        for i_f in range(fn):
            if return_template:
                row.append(Matrix([var('c_{' + str(i_c) + '}d_{' + str(i_f) + '}')]))
            else:
                row.append(h_element(gm_1=c_g, gm_2=f_g, fn=fn, angles_navigation=angles_navigation, r1=r_c[i_c], r2=r_f[i_f], q1=q_c[i_c], q2=q_f[i_f], w_0=w_0, t=t, relation='cd'))
                
        row = block_diag(*row)
        H_cd = row if H_cd is None else vstack(H_cd, row)

    # >>>>>>>>>>>> Нижняя подматрица <<<<<<<<<<<<
    H_dd = None
    for i_y in range(len(ff_sequence)):  # то же самое, что range(int(fn*(fn-1)/2)):
        row = []
        for i_x in range(fn):
            if i_x in ff_sequence[i_y]:
                i_1, i_2 = ff_sequence[i_y] if i_x == ff_sequence[i_y][1] else ff_sequence[i_y][::-1]
                if return_template:
                    row.append(Matrix([var('d_{' + str(i_1) + '}d_{' + str(i_2) + '}')]))
                else:
                    row.append(h_element(gm_1=f_g, gm_2=f_g, fn=fn, angles_navigation=angles_navigation, r1=r_f[i_1], r2=r_f[i_2], q1=q_f[i_1], q2=q_f[i_2], w_0=w_0, t=t, relation='dd'))
            else:
                if return_template:
                    row.append(Matrix([0]))
                else:
                    row.append(zeros((f_take_len * f_send_len, 12 if angles_navigation else 6), template=w_0))
        row = bmat(row)
        H_dd = row if H_dd is None else vstack(H_dd, row)   

    # >>>>>>>>>>>> Компановка <<<<<<<<<<<<
    return vstack(H_cd, H_dd) if H_dd is not None else H_cd
