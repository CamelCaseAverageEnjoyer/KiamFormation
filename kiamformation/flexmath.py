import numpy
import sympy


def _numcheck(a):
    return isinstance(a, numpy.ndarray) or isinstance(a, int) or isinstance(a, float) or isinstance(a, complex) or \
        isinstance(a, numpy.quaternion)


def setvectype(a, template=None):
    """
    setvectype(a) create numpy.array() or sympy.Matrix() from array a.
    Note: 'a' array should contain at least 1 element.
    Examples:
        setvectype([1, 2, 3])
        setvectype(sympy.var("a_1 a_2 a_3"))
    """
    if template is not None:
        return numpy.array(a) if _numcheck(template) else sympy.Matrix(a)
    try:
        return numpy.array(a) if _numcheck(a[0][0]) else sympy.Matrix(a)
    except:
        return numpy.array(a) if _numcheck(a[0]) else sympy.Matrix(a)


def pi(a):
    """
    pi(a) return π value the same type as 'a' parameter.
    For numeric π=3.1415... return parameter 'a' should be (int | float | complex | numpy.ndarray).
    Examples:
        pi(1)
        pi(sympy.abc.x)
    """
    return numpy.pi if _numcheck(a) else sympy.pi


def quat(a):
    """
    quat(a) return quaternion by 'a' array.
    Note: 'a' array should be 4-(or less)-dimentional.
    Note: for numerical return, numpy-quaternion library should be installed.
    Examples:
        quat([1, 2, 3])
        quat([0, 1, 2, 3])
        quat(sympy.var("a_1 a_2 a_3"))
        quat(sympy.var("a_0 a_1 a_2 a_3"))
    """
    import quaternion
    if _numcheck(a[0]):
        import quaternion
        return numpy.quaternion(*a)
    else:
        return sympy.Matrix([0]*(4 - len(a)) + [i for i in a])


def zeros(dims, template):
    """
    zeros() return zeros array with dimentions=dims and type of template variable.
    Examples:
        zeros((3, 3), 1)
        zeros((3, 3), sympy.abc.x)
    """
    return numpy.zeros(dims) if _numcheck(template) else sympy.zeros(*dims)


def eye(dims, template):
    """
    eye() return eye matrix with dimentions=dims and type of template variable.
    Examples:
        eye(3, 1)
        eye(3, sympy.abc.x)
    """
    return numpy.eye(dims) if _numcheck(template) else sympy.eye(dims)


def cross(a, b):
    """
    cross(a, b) returns cross product of two 3-dimentional vectors.
    Note: numpy.cross() realization is not used because of slow implementation.
    Examples:
        cross([1, 2, 3], [4, 5, 6])
        cross(sympy.var("a_1 a_2 a_3"), sympy.var("b_1 b_2 b_3"))
    """
    return setvectype([a[1] * b[2] - a[2] * b[1],
                       a[2] * b[0] - a[0] * b[2],
                       a[0] * b[1] - a[1] * b[0]])


def dot(a, b):
    """
    dot(a, b) returns dot product of two N-dimentional vectors.
    Examples:
        dot([1, 2, 3, 4], [5, 6, 7, 8])
        dot(sympy.var("a_1 a_2 a_3"), sympy.var("b_1 b_2 b_3"))
    """
    if _numcheck(a) or _numcheck(a[0]):
        return numpy.dot(a, b)
    else:
        _a = a if isinstance(a, sympy.Matrix) else sympy.Matrix(a)
        _b = b if isinstance(b, sympy.Matrix) else sympy.Matrix(b)
        return (_a.T @ _b)[0]


def sin(a):
    return numpy.sin(a) if _numcheck(a) else sympy.sin(a)


def cos(a):
    return numpy.cos(a) if _numcheck(a) else sympy.cos(a)


def sqrt(a):
    return numpy.sqrt(a) if _numcheck(a) else sympy.sqrt(a)


def tan(a):
    return numpy.tan(a) if _numcheck(a) else sympy.tan(a)


def arctan(a):
    return numpy.arctan(a) if _numcheck(a) else sympy.atan(a)


def mean(a):
    """
    mean(a) return mean value of array a.
    Note: 'a' array should contain at least 1 element.
    Examples:
        mean([1, 2, 3])
        mean(sympy.var("a_1 a_2 a_3"))
    """
    return numpy.mean(a) if _numcheck(a[0]) else sum(list(a)) / len(a)


def norm(a):
    """
    mean(a) return norm value of array a.
    Examples:
        norm([1, 2, 3])
        norm(sympy.var("a_1 a_2 a_3"))
    """
    return numpy.linalg.norm(a) if (_numcheck(a) or isinstance(a, list)) else sympy.sqrt(dot(a, a))


def inv(a):
    """
    inv(a) return inverse matrix of matrix a.
    Examples:
        inv(numpy.eye(3))
        inv(sympy.Matrix([[sympy.var(f'a_{i}^{j}') for i in range(3)] for j in range(3)]))
    """
    return numpy.linalg.inv(a) if _numcheck(a) else a.inv()


def append(*args):
    """
    append() returns merged array.
    Note: for numerical return, only 2 arrays in input are expected.
    Examples:
        append([1, 2, 3], [4, 5])
        append(sympy.var("a_1 a_2 a_3"), sympy.var("b_1 b_2"), sympy.var("c_1 c_2"))
    """
    return numpy.append(*args) if _numcheck(args[0][0]) else sympy.Matrix(sympy.BlockMatrix([sympy.Matrix(i).T for i in args]).T)


def vstack(*args):
    """
    vstack() return stack arrays in sequence vertically.
    Examples:
        vstack(numpy.array([[1, 2], [3, 4]]), numpy.array([[5, 6], [7, 8]]))
        vstack(sympy.eye(3), sympy.zeros(1, 3))
    """
    return numpy.vstack(args) if _numcheck(args[0][0]) else sympy.Matrix.vstack(*args)


def bmat(*args, template=None):
    """
    bmat() return block matrix from given ones.
    Examples:
        bmat([[numpy.zeros((3, 1)), numpy.eye(3)], [numpy.zeros((1, 4))]])
        bmat([[sympy.zeros(3, 1), sympy.eye(3)], [sympy.eye(1, 1), sympy.zeros(1, 3)]])
    """
    if template is not None:
        return numpy.bmat(*args) if _numcheck(template) else sympy.Matrix(sympy.BlockMatrix(*args))
    try:
        return numpy.bmat(*args) if _numcheck(args[0][0][0]) else sympy.Matrix(sympy.BlockMatrix(*args))
    except:
        return numpy.bmat(*args) if _numcheck(args[0][0]) else sympy.Matrix(sympy.BlockMatrix(*args))


def block_diag(*args):
    """
    block_diag() return block-diagonal matrix from given ones.
    Examples:
        block_diag(numpy.eye(2), numpy.zeros((1,2)), numpy.eye(2))
        block_diag(sympy.eye(2), sympy.zeros(1,2), sympy.eye(2))
    """
    if _numcheck(args[0]):
        from scipy.linalg import block_diag
        return block_diag(*args)
    else:
        return sympy.Matrix(sympy.BlockDiagMatrix(*args))


def get_same_type_conversion(a):
    if isinstance(a, list):
        return list
    if isinstance(a, numpy.ndarray):
        return numpy.array
    if isinstance(a, sympy.Matrix):
        return sympy.Matrix


def float2rational(a, *args):
    """
    float2rational() replaces float number in sympy expression to fraction.
    Example:
        a = 1/2 * sympy.abc.x + 1/3 * sympy.abc.y
        a = float2rational(a, (1, 2), (1, 3))
    """
    return a if _numcheck(a) else a.subs([(sympy.Float(i[0] / i[1]), sympy.Rational(i[0], i[1])) for i in args])
