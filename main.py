"""Численное моделирование с использованием интерфейса."""
if __name__ == '__main__':
    import sys
    import numpy as np
    import pandas as pd
    from kiamformation.interface import interface_window
    from kiamformation.config import init
    from warnings import simplefilter

    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    # Инициализация объектов
    o = init()
    np.set_printoptions(linewidth=300)

    # Интерфейс
    app, window = interface_window(o=o)
    sys.exit(app.exec_())
