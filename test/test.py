"""Тестирование проекта kiam-formation.

Добавить:
- Кватернион в rk4 принимается и отдаётся N-вектором (передали 3 - получили 3)
"""
import numpy as np
import sympy
import sys
sys.path.insert(1, f"../")
import kiamformation as kf


def test_save_load_delete():
    assert (1, 2, 3) == (1, 2, 3)


def test_1():
    assert (1, 2, 3) == (1, 2, 3)

def test_flexmath():
    a1, a2, a3, a4 = sympy.symbols("a_1 a_2 a_3 a_4")

    a = np.array([0, 1, 2, 3, 4])
    b = kf.append([0, 1], [2, 3, 4])
    assert (a == b).all()

    a = sympy.Matrix([a1, a2, a3, a4])
    b = kf.append(sympy.Matrix([a1, a2]), sympy.Matrix([a3, a4]))
    assert a == b

    a = (1 + 2 + 3 + 4) / 4
    b = kf.mean([1, 2, 3, 4])
    assert a == b

    a = (1 + 2 + 3 + 4) / 4
    b = kf.mean(np.array([1, 2, 3, 4]))
    assert a == b

    a = (a1 + a2 + a3 + a4) / 4
    b = kf.mean(sympy.Matrix([a1, a2, a3, a4]))
    assert a == b

def test_kiam_astro():
    from kiam_astro import trajectory, kiam
    print(kiam.units('earth')['TimeUnit'] * 24 * 2 * np.pi)
