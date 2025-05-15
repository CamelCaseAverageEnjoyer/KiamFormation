"""
Пакет kiam-formation для численного моделирования задачи навигации космических аппаратов в групповом полёте.
"""
import sys
sys.path.insert(1, f"{sys.path[0]}/kiamformation/.")
sys.path.insert(1, f"{sys.path[0]}/kiamformation")
# sys.path.insert(1, f"{sys.path[0]}/test")

from kiamformation.config import *
from kiamformation.static.cosmetic import *
from kiamformation.dynamics import *
from kiamformation.flexmath import *
from kiamformation.navigation import *
from kiamformation.H_matrix import *
# from .interface import *
from kiamformation.math_tools import *
from kiamformation.visualization import *
from kiamformation.measurements import *
from kiamformation.simulation import *
from kiamformation.spacecrafts import *
