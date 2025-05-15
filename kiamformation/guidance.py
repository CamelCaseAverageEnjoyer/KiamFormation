from kiamformation.spacecrafts import Apparatus
from kiamformation.config import Objects

# >>>>>>>>>>>> Guidance <<<<<<<<<<<<
def guide(c: Apparatus, f: Apparatus, o: Objects, earth_turn: float) -> None:
    """Функция включает и выключает КА. Моделирование режимов работы."""
    pass
    '''for obj in [c, f]:
        for i in range(obj.n):
            if obj.operating_modes[i] == v.OPERATING_MODES_CHANGE[1]:  # Отсутствие аккумулятора на чипсате
                if earth_turn % 1 < 0.5 and obj.operating_mode[i] == v.OPERATING_MODES[-1]:
                    obj.operating_mode[i] = v.OPERATING_MODES[0]
                if earth_turn % 1 > 0.5 and obj.operating_mode[i] != v.OPERATING_MODES[-1]:
                    obj.operating_mode[i] = v.OPERATING_MODES[-1]'''
