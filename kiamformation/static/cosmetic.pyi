# cosmetic.pyi
def my_print(txt: any, color: str = None, if_print: bool = True, bold: bool = False, if_return: bool = False,
             end: str = '\n') -> None:
    """Функция вывода цветного текста
    :param txt: Выводимый текст
    :param color: Цвет текста {b, g, y, r, c, m}
    :param if_print: Флаг вывода для экономии места
    :param bold: Жирный текст
    :param if_return: Надо ли возвращать строку"""
    pass


def ending(n: int) -> str:
    """Окончание слова в зависимости от слова."""
    pass


def talk(aloud=True):
    """Хорошая штучка."""
    pass


