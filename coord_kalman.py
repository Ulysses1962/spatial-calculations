# coding: cp1251
__author__ = 'Fliegende Hollander'

from math import sin, cos, tan, radians, sqrt, pi, log, exp, atan, fabs, pow, hypot, floor
from numpy import matrix, linalg
from random import uniform, randint
import winsound as ws

#========================================================================================
# Процедуры преобразования координат точки между геодезическими WGS84 и проекцией Меркатора
#========================================================================================
def WGS84_to_Mercator(lat, lon):
    """
    Производит перерасчет геодезических координат на геоиде WGS84 в координаты проекции Меркатора
    :param lat:     - Широта заданной точки, градусы
    :param lon:     - Долгота заданной точки. градусы
    :return:        - Координаты проекции Меркатора для заданной точки
    """
    a = 6378137.0
    b = 6356752.3142
    f = (a - b) / a
    e = sqrt(f * (2 - f))
    x = a * radians(lon)
    i_1 = tan(pi / 4 + radians(lat) / 2)
    i_2 = ((1 - e * sin(radians(lat))) / (1 + e * sin(radians(lat)))) ** (e / 2)
    y = a * log(i_1 * i_2)

    return x, y
def Mercator_to_WGS84(x, y, eps = 0.000000001):
    """
    Производит перерасчет координат проекции Меркатора в геодезическую долготу и широту
    :param x:   - X проекции Меркатора заданной точки
    :param y:   - Y проекции Меркатора заданной точки
    :param eps: - Точность вычисления широты
    :return:    - Геодезические координаты на геоиде WGS84, градусы
    """
    a = 6378137.0
    b = 6356752.3142
    f = (a - b) / a
    e = sqrt(f * (2 - f))
    lon = (180 * x / a) / pi
    # Итерационный процесс вычисления широты
    eh  = e / 2
    pih = pi / 2
    ts  = exp(-y / a)
    PHI = pih - 2 * atan(ts)
    i    = 0
    dPHI = 1
    while fabs(dPHI) > eps and i < 15:
        con = e * sin(PHI)
        dPHI = pih - 2 * atan(ts * pow(((1 - con) / (1 + con)), eh)) - PHI
        PHI += dPHI
        i += 1

    lat = (180 * PHI) / pi

    return lat, lon

#========================================================================================
# Процедуры преобразования координат точки между геодезическими WGS84 и проекцией ECEF
#========================================================================================
def WGS84_to_ECEF(lat, lon, alt = 0):
    """
    Производит перерасчет геодезических координат на геоиде WGS84 в координаты ECEF (геоцентрическая система координат)
    :param lat:     - Широта заданной точки
    :param lon:     - Долгота заданной точки
    :param alt:     - Высота над уровнем моря заданной точки, м
    :return:        - ECEF-координаты заданной точки
    """
    a = 6378137.0
    b = 6356752.3142
    f = (a - b) / a
    e = sqrt(f * (2 - f))
    Nfi = a / sqrt(1 - (e * sin(radians(lat))) ** 2.0)
    x = (Nfi +alt) * cos(radians(lat)) * cos(radians(lon))
    y = (Nfi + alt) * cos(radians(lat)) * sin(radians(lon))
    z = (Nfi * (1 - e ** 2.0) + alt) * sin(radians(lat))

    return  x, y, z
def ECEF_to_WGS84(x, y, z, eps=0.000000001):
    """
    Производит перерасчет координтат в системе ECEF в геодезические на геоиде WGS84
    :param x:   - X ECEF заданной точки
    :param y:   - Y ECEF заданной точки
    :param z:   - Z ECEF заданной точки
    :param eps: - Точность вычисления широты и высоты над уровнем моря
    :return:    - Геодезические координаты заданной точки на геоиде  WGS84
    """
    a = 6378137.0
    b = 6356752.3142
    f = (a - b) / a
    e = sqrt(f * (2 - f))
    lon = 180 * atan(y / x) / pi
    # Итерационный процесс вычисления широты и высоты над уровнем моря
    p = sqrt(x ** 2.0 + y ** 2.0)
    dfi = 1.0
    h   = 0
    fi  = atan(z /(p * (1 - e ** 2.0)))
    while fabs(dfi) > eps:
        N = a / sqrt(1 - (e * sin(fi)) ** 2.0)
        h = p / cos(fi) - N
        fi_N = atan(z / (p * (1 - (N * e ** 2.0) / (N + h))))
        dfi = fi_N -fi
        fi = fi_N

    lat = 180 * fi / pi
    alt = h

    return lat, lon, alt

#========================================================================================
# Вспомогательные процедуры пересчета координат WGS84 - Universal Transverse Mercator
#========================================================================================
# Расчет параметров геоида
def geoid_params():
    """
    Процедура вычисления характеристик геоида
    :return:    Набор характеристик геоида
    """
    a = 6378137.0
    b = 6356752.3142
    f    = (a - b) / a
    e_2  = f * (2 - f)
    e_21 = e_2 / (1 - e_2)
    n    = f / (2 - f)
    A = a * (1 - n + 5 * (n ** 2 - n ** 3) / 4. + 81 * (n ** 4 - n ** 5) / 64.)
    B = 3 * a * (n - n ** 2 + 7 * (n ** 3 - n ** 4) / 8. + 55 * (n ** 5) / 64.) / 2.
    C = 15 * a * (n ** 2 - n ** 3 + 3 * (n ** 4 - n ** 5) / 4.) / 16.
    D = 35 * a * (n ** 3 - n ** 4 + 11 * (n ** 5) / 16.) /48.
    E = 315 * a * (n ** 4 - n ** 5) / 512.

    return a, b, f, e_2, e_21, n, A, B, C, D, E
# Расчет степение синуса угла (для ускорения расчета
def sins(fi):
    """
    Процедура вычисления степеней синуса широты (для ускорения вычислений)
    :param fi:    Широта в радианах
    :return:      Массив степеней синуса
    """
    s = []
    for i in range(1, 9):
        s.append(sin(fi) ** i)

    return s
# Расчет степеней косинуса угла
def coss(fi):
    """
    Процедура вычисления степеней косинуса широты (для ускорения вычислений)
    :param fi:    Широта в радианах
    :return:      Массив степеней косинуса
    """
    s = []
    for i in range(1, 9):
        s.append(cos(fi) ** i)

    return s
# Расчет степеней тангенса угла
def tans(fi):
    """
    Процедура вычисления степеней тангенса широты (для ускорения вычислений)
    :param fi:    Широта в радианах
    :return:      Массив степеней тангенса
    """
    s = []
    for i in range(1, 9):
        s.append(tan(fi) ** i)

    return s
# Расчет длины меридианной дуги для заданной широты
def S(fi):
    """
    Расчет длины меридианной дуги
    :param fi:  Широта точки в радианах
    :return:    Длина меридиональной дуги
    """
    geoid = geoid_params()
    A = geoid[6]
    B = geoid[7]
    C = geoid[8]
    D = geoid[9]
    E = geoid[10]
    S = A * fi - B * sin(2 * fi) + C * sin(4 * fi) - D * sin(6 * fi) + E * sin(8 * fi)

    return S
# Расчет широты основаня перпендикуляра от точки к центральному меридиану зоны
def UTM_footpoint(N, hemisphere, eps = 0.000000001):
    """
    Вычисление широты основания перпендикуляра от точки до центрального меридиана зоны
    :param N:           Координата N точки
    :param eps:         Точность вычислений
    :param hemisphere:  Полушарие - северное 'N', южное - 'S'
    :return:            Footpoint-широта заданной точки
    """
    k0 = 0.9996
    fi = 0
    dN  = 1
    # Максимальное значение дуги
    Nm  = S(radians(84.)) * k0 if hemisphere == 'N' else S(radians(80.)) * k0
    # Итерационной вычисление искомой широты
    while fabs(dN) > eps:
        N0 = S(fi) * k0
        dN = N - N0
        fi += dN / Nm

    return fi
# Вычисление термов для прямого преобразования WGS84 -> UTM
def T_WGS84_UTM(fi):
    """
    Вычисляет термы для прямого преобразования WGS84 в UTM
    :param fi:      Широта в радианах
    :return:        Массив термов
    """
    k0 = 0.9996
    g = geoid_params()
    s = sins(fi)
    c = coss(fi)
    t = tans(fi)
    nu  = g[0] / sqrt(1 - g[3] * s[1])
    T1  = S(fi) * k0
    T2  = nu * s[0] * c[0] * k0 / 2.
    T3  = nu * s[0] * c[2] * k0 * (5 - t[1] + 9 * g[4] * c[1] + 4 * g[4] ** 2 * c[3]) / 24.
    T4  = nu * s[0] * c[4] * k0 * (61 - 58 * t[1] + t[3] * 270 * g[4] * c[1] - 330 * t[1] * g[4] * c[1] + \
          445 * g[4] ** 2 *c[3] + 324 * g[4] ** 3 * c[5] - 680 * t[1] * g[4] ** 2 * c[3] + \
          88 * g[4] ** 4 * c[7] - 600 * t[1] * g[4] ** 3 * c[5] - 192 *t[1] * g[4] ** 4 * c[7]) / 720.
    T5  = nu * s[0] * c[6] * k0 * (1385 - 3111 * t[1] + 543 * t[3] - t[5]) / 40320.
    T6  = nu * c[0] * k0
    T7  = nu * c[2] * k0 * (1 - t[1] + g[4] * c[1]) / 6.
    T8  = nu * c[4] * k0 * (5 - 18 * t[1] + t[3] + 14 * g[4] * c[1] - 58 *t[1] * g[4] * c[1] + \
          13 * g[4] ** 2 * c[3] + 4 * g[4] ** 3 * c[5] - 64 * t[1] * g[4] ** 2 * c[3] - 24 * t[1] * g[4] ** 3 * c[5]) / 120.
    T9  = nu * c[6] * k0 * (61 - 479 * t[1] + 179 * t[3] - t[5]) / 5040.

    return T1, T2, T3, T4, T5, T6, T7, T8, T9
# Вычисление термов для обратного преобразования UTM -> WGS84
def T_UTM_WGS84(N, hemisphere):
    """
    Вычисляет термы для обратного преобразования UTM и WGS84
    :param N:           Северное склонение точки
    :param hemisphere:  Полушарие
    :return:            Массив термов
    """
    k0 = 0.9996
    fi1 = UTM_footpoint(N, hemisphere)

    g = geoid_params()
    s = sins(fi1)
    c = coss(fi1)
    t = tans(fi1)
    nu  = g[0] / sqrt(1 - g[3] * s[1])
    ro  = g[0] * (1 - g[3]) / sqrt((1 - g[3] * s[1]) ** 3)
    T10 = t[0] / (2 * ro * nu * k0 ** 2)
    T11 = t[0] * (5 + 3 * t[1] + g[4] * c[1] - 4 * g[4] ** 2 * c[3] - 9 * t[1] * g[4] * c[1]) / (24 * ro * nu ** 3 * k0 ** 4)
    T12 = t[0] * (61 + 90 * t[1] + 46 * g[4] * c[1] + 45 * t[3] - 252 * t[1] * g[4] * c[1] - \
          3 * g[4] ** 2 * c[3] + 100 * g[4] ** 3 * c[5] - 66 * t[1] * g[4] ** 2 * c[3] - \
          90 * t[3] * g[4] * c[1] + 88 * g[4] ** 4 * c[7] + 225 * t[3] * g[4] ** 2 * c[3] + \
          84 * t[1] * g[4] ** 3 * c[5] - 192 * t[1] * g[4] ** 4 * c[7]) / (720 * ro * nu ** 5 * k0 ** 6)
    T13 = t[1] * (1385 + 3633 * t[1] + 4095 * t[3] + 1575 * t[5]) / (40320 * ro * nu ** 7 * k0 ** 8)
    T14 = 1. / (nu * c[0] * k0)
    T15 = (1 + 2 * t[1] + g[4] * c[1]) / (6 * nu ** 3 * c[0] * k0 ** 3)
    T16 = (5 + 6 * g[4] * c[1] + 28 * t[1] - 3 * g[4] ** 2 * c[3] + 8 * t[1] * g[4] * c[1] + \
           24 * t[3] - 4 * g[4] ** 3 * c[5] + 4 * t[1] * g[4] ** 2 * c[3] + 24 * t[1] * g[4] ** 3 * c[5]) / (120 * nu ** 5 * c[0] * k0 ** 5)
    T17 = (61 + 6621 * t[1] + 1320 * t[3] + 720 * t[5]) / (5040 * nu ** 7 * c[0] * k0 ** 7)

    return T10, T11, T12, T13, T14, T15, T16, T17, fi1

#========================================================================================
# Процедуры пересчета координат WGS84 <-> UTM
#========================================================================================
def WGS84_UTM(fi, hemisphere, lo, easting):
    """
    Процедура прямого преобразования координат WGS84 -> UTM
    :param fi:          Широта в градусах
    :param hemisphere:  Полушарие N - северное, S - южное
    :param lo:          Долгота в градусах
    :param easting:     Положение относительно нулевого меридиана E - долгота восточная, W - западная
    :return:            Координаты точки в проекции Universal Transverse Mercator, номер зоны
    """
    FN = 0 if hemisphere == 'N' else 10000000
    FE = 500000
    zone = (30 + floor(lo / 6) + 1) if easting == "E" else (floor((180 - lo) / 6) + 1)
    lo_0 = radians((zone - 1) * 6 - 180 + 3) if easting == 'E' else radians(180 - (zone - 1) * 6 + 3)
    dlo  = radians(lo) - lo_0

    T = T_WGS84_UTM(radians(fi))

    N = FN + (T[0] + dlo ** 2 * T[1] + dlo ** 4 * T[2] + dlo ** 6 * T[3] + dlo ** 8 * T[4])
    E = FE + (dlo * T[5] + dlo ** 3 * T[6] + dlo ** 5 * T[7] + dlo ** 7 * T[8])

    return N, E, zone
def UTM_WGS84(N, E, zone, hemisphere):
    """
    Процедура обратного преобразования UTM -> WGS84
    :param N:           Северное склонение в метрах
    :param E:           Восточное склонение в метрах
    :param zone:        Номер зоны UTM
    :param hemisphere:  Полушарие - северное или южное
    :return:            Геодезические широта и долгота в радианах
    """
    lo_0 = radians((zone - 1) * 6 - 180 + 3) if zone > 30 else radians(180 - (zone - 1) * 6 + 3)
    FE = 500000
    DE = E - FE
    T = T_UTM_WGS84(N, hemisphere)

    lat = T[8] - DE ** 2 * T[0] + DE ** 4 * T[1] - DE ** 6 * T[2] + DE ** 8 * T[3]
    lon = lo_0 + DE * T[4] - DE ** 3 * T[5] + DE ** 5 * T[6] - DE ** 7 * T[7]

    return lat, lon

#========================================================================================
# Процедуры вычисления радиуса наименьшей охватывающей окружности множества точек на плоскости
#========================================================================================
# Построение окружности по 2-м точкам на плоскости
def circle_2p(p1, p2):
    """
    :param p1:  Первая диаметральная точка окружности
    :param p2:  Вторая диаметральная точка окружности
    :return:    x, y, d - координаты центра и диаметр окружности
    """
    d = hypot(p2[0] - p1[0], p2[1] - p1[1])
    x = (p2[0] + p1[0]) / 2.
    y = (p2[1] + p1[1]) / 2.

    return x, y, d
# Построение окружности по 3-м точкам на поскости
def circle_3p(p1, p2, p3):
    """
    :param p1:  Первая точка окружности
    :param p2:  Вторая точка окружности
    :param p3:  Третья точка окружности
    :return:    x, y, d - координаты центра и диаметр окружности
    """
    x1 = p1[0]; y1 = p1[1]
    x2 = p2[0]; y2 = p2[1]
    x3 = p3[0]; y3 = p3[1]

    if x1 == x2 and x2 != x3:
        mb = (y3 - y2) / (x3 - x2)
        Xc = (mb * (y3 - y1) + x2 + x3) / 2.
        Yc = -1. *(Xc - (x2 + x3) / 2.) / mb + (y2 + y3) / 2.
    elif x1 != x2 and x2 == x3:
        ma = (y2 - y1) / (x2 - x1)
        Xc = (ma * (y1 - y3) + x2 + x1) / 2.
        Yc = -1. *(Xc - (x2 + x1) / 2.) / ma + (y2 + y1) / 2.
    else:
        ma = (y2 - y1) / (x2 - x1)
        mb = (y3 - y2) / (x3 - x2)
        Xc = (ma * mb * (y1 - y3) + mb * (x1 + x2) - ma * (x2 + x3)) / (2. * (mb - ma))
        Yc = - 1. * (Xc - (x1 + x2) / 2.) / ma + (y1 + y2) / 2.

    D = hypot(Xc - x1, Yc - y1) * 2.

    return Xc, Yc, D
# Определение попадания точки внутрь заданной окружности
def is_in_circle(p, c, eps = 1e-12):
    """
    :param p:   Точка, для которой нужно произвести проверку попадания внутрь окружности
    :param C:   Окружность, относительно которой проверяется точка
    :param eps: Точность определения попадания точки внутрь окружности
    :return:    true, если попали, в противном случае - false
    """
    return hypot(p[0] - c[0], p[1] - c[1]) < (c[2] / 2. + eps)
# Вычисление минимальной охватывающей окружности для множества точек на плоскости
def min_enclosing_circle(points):
    """
    :param points:  Массив отсчетов координат в системе ECEF или Меркатора
    :return:        диаметр наименьшей охватывающей окружности, м
    """
    c = circle_2p(points[0], points[1])
    for i in range(2, len(points)):
        if not is_in_circle(points[i], c):
            c = circle_2p(points[0], points[i])
            for j in range(1, i):
                if not is_in_circle(points[j], c):
                    c = circle_2p(points[i], points[j])
                    for k in range(0, j):
                        if not is_in_circle(points[k], c):
                            c = circle_3p(points[i], points[j], points[k])

    return c

#========================================================================================
# Фильтр Калмана
#========================================================================================
class KF(object):
    """
    Реализация многомерного фильтра Калмана для обработки измерений координат
    """
    # Filter estimate state vector
    X_ = matrix([[0.], [0.], [0.], [0.]])
    # Filter a posteriori state vector
    X = matrix([[0.], [0.], [0.], [0.]])
    # Filter prediction matrix
    F = matrix([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])
    # Filter observation model matrix
    H = matrix([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])
    # Filter process noise covariance matrix
    Q = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
    # Filter measurement noise covariance matrix
    R = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
    # System a priori error covariance matrix
    P_ = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
    # System a posteriori error covariance matrix
    P = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
    # Measurements vector
    Z = matrix([[0., 0., 0., 0.]])
    # Identity matrix
    I = matrix([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])

    # Filter constructor
    def __init__(self):
        # Filter time quantum
        self.DT = 0.

    # Set noise covariance matricies
    def SetQR(self, Q, R):
        """
        :param Q:   Матрица комариации шума процесса
        :param R:   Матрица ковариации шума измерений
        """
        self.Q = Q
        self.R = R

    # Initial data and system covariance matricies setting
    def SetInitialData(self, X0, P0):
        """
        Установка начальных параметров для фильтрации
        :param X0: Матрица начальных значений состояния
        :param P0: Матрица начальных ковариаций ошибок
        """
        self.X = X0
        self.P = P0

    # Prediction matrix calculation
    def CalcF(self, dt):
        """
        Вычисление матрицы предсказания состояния системы
        :param dt: Временной интервал
        """
        self.DT = dt
        A = self.DT * cos(self.X.item(2))
        B = self.DT * sin(self.X.item(2))
        self.F.itemset((0, 3), A)
        self.F.itemset((1, 3), B)

    # Prediction stage
    def Predict(self, dt):
        """
        Предсказание состояния системы
        :param dt: Временной интервал
        """
        self.CalcF(dt)
        # Estimate state
        self.X_ = self.F * self.X
        # A priori error covariance
        FT = self.F.getT()
        self.P_ = self.F * self.P * FT + self.Q

    # State update stage
    def Correct(self, Z, dt):
        """
        Корректировка вектора состояния фильтра
        :param Z:   Вектор измерений (МАТРИЦА-СТОЛБЕЦ!)
        :param dt:  Временной интервал
        :return:    Текущий массив состояния системы
        """
        # Predict state
        self.Predict(dt)
        # Gain matrix calculation
        HT = self.H.getT()
        MM = self.H * self.P_ * HT - self.R
        S = linalg.inv(MM)
        K = self.P_ * HT * S
        # A posteriori estimate of system state
        self.X = self.X_ + K * (Z - self.H * self.X_)
        #A posteriori error covariance
        self.P = (self.I - K * self.H) * self.P_

        return self.X

#========================================================================================
# Вспомогательные функции калибровки фильтра Калмана
#========================================================================================
# Вычисление массива тестовых данных
def TestDataCalc(Xc, Yc):
    """
    Вычисляет массив тестовых данных, по которым калибруется фильтр
    :return: Массив из 18 точек, отстоящих от центра на равном расстоянии
             под разными азимутами. Координаты центра тестовой фигуры - (0, 0)
    """
    TM = []
    distance = 30.0
    delta    = radians(20.0)
    for i in range(18):
        teta = i * delta
        x = Xc + round(distance * sin(teta), 10)
        y = Yc + round(distance * cos(teta), 10)
        TM.append(matrix([[x], [y], [teta], [0.]]))

    return TM
# Процедура калибровки фильтра
def Particulate(Xc, Yc):
    """
    Производит калибровку фильтра Калмана
    :return: Оптимальные матрицы ковариаций Q и R
    """
    #ta = TestDataCalc(Xc, Yc)
    TEST_DATA = []
    tf = open("RTCM0.txt")
    i = 0
    for line in tf:
        line = line.rstrip()
        lst = line.split("\t")
        N, E, ZONE = WGS84_UTM(float(lst[0]), 'N', float(lst[1]), 'E')
        TEST_DATA.append(matrix([[E], [N], [0.], [0.]]))
        i += 1
    tf.close()

    print " Starting filter particulation. Velocity is 2 m/s"
    min_RMSE = 0
    for j in range(30000):
        Q  = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
        R  = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
        P0 = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
        for i in range(4):
            Q.itemset((i, i), uniform(0, 6))
            R.itemset((i, i), uniform(0, 6))
        # Вычисление фильтрованных тестовых данных
        KFL = KF()
        KFL.SetQR(Q, R)
        X0 = matrix([[Xc], [Yc], [0.], [0.]])
        KFL.SetInitialData(X0, P0)
        RES = []
        for k in range(18):
            X = KFL.Correct(TEST_DATA[k], 0.0)
            RES.append(X)
        # Вычисляем расхождение по методу наименьших квадратов
        s = 0
        for k in range(len(RES)):
            s += ((RES[k].item(0) - Xc) ** 2 + (RES[k].item(1) - Yc) ** 2)
        RMSE = sqrt(s / len(RES))
        if j == 0:
            min_RMSE = RMSE
            min_Q    = Q
            min_R    = R
        else:
            if RMSE < min_RMSE:
                min_RMSE = RMSE
                min_Q    = Q
                min_R    = R

    # Вывод оптимальных матриц ковариаций
    print "================================================================="
    print "    OPTIMAL COVARIANCE MATRICIES"
    print "================================================================="
    print "    Q matrix is "
    print "-----------------------------------------------------------------"
    print min_Q
    print "-----------------------------------------------------------------"
    print "    R matrix is "
    print "-----------------------------------------------------------------"
    print min_R
    print "-----------------------------------------------------------------"
    print "    Minimal RMSE is ", min_RMSE

    return min_Q, min_R

#========================================================================================
# Вспомогательные функции подготовки тестовых наборов данных
#========================================================================================
def GetTESTDATA(Xc, Yc):
    TA = TestDataCalc(Xc, Yc)
    TEST_DATA = []
    tf = open("7082015.txt")
    i = 0
    for line in tf:
        line = line.rstrip()
        lst = line.split("\t")
        N, E, ZONE = WGS84_UTM(float(lst[0]), 'N', float(lst[1]), 'E')
        TEST_DATA.append(matrix([[E], [N], [TA[i].item(2)], [2.]]))
        i += 1
    tf.close()

    return TA, TEST_DATA


if __name__ == "__main__":
    f = open("RTCM0.txt")
    coordinates = []
    for line in f:
        line = line.rstrip()
        lst = line.split("\t")
        N, E, ZONE = WGS84_UTM(float(lst[0]), 'N', float(lst[1]), 'E')
        coordinates.append((E, N))
    f.close()

    c = min_enclosing_circle(coordinates)

    print "----------------------------------------------------------------"
    print " ENCLOSING CIRCLE PARAMETERS"
    print "----------------------------------------------------------------"
    print " Xc= ", c[0], "\t", "Yc= ", c[1], "\t", "D= ", c[2]
    print "----------------------------------------------------------------"

    Q  = matrix([[0.00443241, 0., 0., 0.], [0., 4.32669338, 0., 0.], [0., 0., 4.71483075, 0.], [0., 0., 0., 4.58424025]])
    R  = matrix([[3.57602215, 0., 0., 0.], [0., 0.19994744, 0., 0.], [0., 0., 5.84132996, 0.], [0., 0., 0., 2.99903718]])
    #Q  = matrix([[1.96802013, 0., 0., 0.], [0., 4.94848626, 0., 0.], [0., 0., 2.61295058, 0.], [0., 0., 0., 5.62512631]])
    #R  = matrix([[0.09670168, 0., 0., 0.], [0., 0.00017878, 0., 0.], [0., 0., 2.30767356, 0.], [0., 0., 0., 1.25258675]])
    P0 = matrix([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
    # Вычисление фильтрованных тестовых данных
    KL = KF()
    KL.SetQR(Q, R)
    X0 = matrix([[c[0]], [c[1]], [0.], [0.]])
    KL.SetInitialData(X0, P0)
    RES = []
    for k in range(len(coordinates)):
        Z = matrix([[coordinates[k][0]], [coordinates[k][1]], [0.], [0.]])
        X = KL.Correct(Z, 1.0)
        RES.append(X)

    flt_coord = []
    for i in range(len(RES)):
        flt_coord.append((RES[i].item(0), RES[i].item(1)))

    c1 = min_enclosing_circle(flt_coord)

    print "----------------------------------------------------------------"
    print " FILTERED SET ENCLOSING CIRCLE PARAMETERS"
    print "----------------------------------------------------------------"
    print " Xc= ", c1[0], "\t", "Yc= ", c1[1], "\t", "D= ", c1[2]
    print "----------------------------------------------------------------"
    '''
    Xc, Yc, ZONE = WGS84_UTM(49.4345794666, 'N', 27.0021217333, 'E')
    Particulate(Xc, Yc)
    '''
    ws.PlaySound('SystemExit', ws.SND_ALIAS)



