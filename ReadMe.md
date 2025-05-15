### KiamFormation
| <img src="kiamformation/static/robot1.png" alt="robot image" width="50"/> | "Помнишь, я тебе говорила про мусор, который стоит? Стоит и смердит? Так вот — это была метафора. Я имела в виду тебя." |
|---------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| 1) При запуске                                                            | Прямой запуск модуля `interface` (интерфейс PyQt5)                                                                      |
| 2) При импорте                                                            | Импорт всех подмодулей кроме `interface`                                                                                |

Пакет `kiam-formation` предназначен для численного моделирования навигации централизованного группового полёта космических аппаратов (КА). Навигация основана на RSSI.

#### <u>1. Запуск модуля</u>
``` console
python3 kiamformation
```
В открытом окне доступны настройка параметров и запуск численного моделирования.

| Иконка                                                     | Функция                              | Описание                                                                                                                                                                  |
|------------------------------------------------------------|--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| <img src="kiamformation/static/robot1.png" width="30"/>    | `static.cosmetic.talk()`             | Случайная фраза.                                                                                                                                                          |
| <img src="kiamformation/static/integral.png" width="30"/>  | `interface.main_run()`               | Запуск/продолжение численного моделирования (с отображением `my_plot.plot_distance()`)                                                                                    |
| <img src="kiamformation/static/plot.png" width="30"/>      | `my_plot.plot_distance()`            | 2D-отрисовка результатов численного моделирования (ошибки навигации, численные критерии наблюдаемости)                                                                    |
| <img src="kiamformation/static/orbit.png" width="30"/>     | `my_plot.plot_all()`                 | 3D-отрисовка результатов численного моделирования в браузере по умолчанию. Пример отрисовки: <img src="images/example.gif" width="300">                                   |
| <img src="kiamformation/static/param.png" width="30"/>     | `interface.plot_1_param()`           | Ручная отрисовка записанных параметров из таблицы `dynamic.PhysicModel.record` (после выбора каждого параметра "ок", для отображения выбранных параметров "cancel")       |
| <img src="kiamformation/static/antenna.png" width="30"/>   | `my_plot.plot_model_gain()`          | Отрисовка диаграмм направленностей для выбранных антенн материнских и дочерних КА (выбор модели в `spacecrafts.local_dipole`)                                             |
| <img src="kiamformation/static/air.png" width="30"/>       | `my_plot.plot_atmosphere_models()`   | Отрисовка доступных моделей расчёта плотности атмосферы (выбор модели в `config.Variables.ATMOSPHERE_MODEL`)                                                              |
| <img src="kiamformation/static/animation.png" width="30"/> | `my_plot.animate_reference_frames()` | Анимирование орбитального движения в `/localfiles/res.gif`. Нужна для валидации моделирования. Создаёт анимацию следующего вида: <img src="images/earth.gif" width="300"> |
| <img src="kiamformation/static/save.png" width="30"/>      | `interface.save_trajectories()`      | Сохраняет результаты численного моделирования в `/srs/kiamformation/data/trajectories`                                                                                    |
| <img src="kiamformation/static/load.png" width="30"/>      | `interface.load_trajectories()`      | Загружает результаты численного моделирования из `/srs/kiamformation/data/trajectories`                                                                                   |
| <img src="kiamformation/static/eraser.png" width="30"/>    | `interface.remove_trajectories()`    | Удаляет результаты численного моделирования в `/srs/kiamformation/data/trajectories`                                                                                      |


#### <u>2. Импорт модуля</u>
``` Python
import kiamformation as kf
o = kf.init()
```

Переменная `o` класса `config.Objects` содержит:
1) `v` - переменная класса `config.Variables` с параметрами моделирования
2) `a` - переменная класса `spacecrafts.Anchor(Apparatus)` (мнимый КА в центре ОСК)
3) `c` - переменная класса `spacecrafts.CubeSat(Apparatus)` (материнские КА одного типа)
4) `f` - переменная класса `spacecrafts.FemtoSat(Apparatus)` (дочерние КА одного типа)
5) `p` - переменная класса `dynamic.PhysicModel` с функциями численного моделирования

Автор пакета использует импорт модуля в целях символьного анализа проблемы навигации. Для этого следует задать параметры окружения и КА в символьном виде:
``` Python
from sympy import *
```

#### <u>3. Содержание пакета kiam-formation</u>
