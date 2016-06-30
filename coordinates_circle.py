# coding: cp1251
__author__ = 'Fliegende Hollander'

import math
import matplotlib.pyplot as plt

lat_E = 49.4345794666
lon_E = 27.0021217333


def read_coordinates(name):
    """
    Reads GPS samples from file
    :param name:    File name to process
    :return:        Coordinates points list
    """
    f = open(name)
    coordinates = []
    for line in f:
        line = line.rstrip()
        lst = line.split("\t")
        coordinates.append((float(lst[0]), float(lst[1])))

    f.close()
    return coordinates

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
    e = math.sqrt(f * (2 - f))
    x = a * math.radians(lon)
    i_1 = math.tan(math.pi / 4 + math.radians(lat) / 2)
    i_2 = ((1 - e * math.sin(math.radians(lat))) / (1 + e * math.sin(math.radians(lat)))) ** (e / 2)
    y = a * math.log(i_1 * i_2)

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
    e = math.sqrt(f * (2 - f))
    lon = (180 * x / a) / math.pi
    # Итерационный процесс вычисления широты
    eh  = e / 2
    pih = math.pi / 2
    ts  = math.exp(-y / a)
    PHI = pih - 2 * math.atan(ts)
    i    = 0
    dPHI = 1
    while math.fabs(dPHI) > eps and i < 15:
        con = e * math.sin(PHI)
        dPHI = pih - 2 * math.atan(ts * math.pow(((1 - con) / (1 + con)), eh)) - PHI
        PHI += dPHI
        i += 1

    lat = (180 * PHI) / math.pi

    return lat, lon

#========================================================================================
# Процедуры преобразования координат точки между геодезическими WGS84 и проекцией ECEF
#========================================================================================
def WGS84_to_ECEF(lat, lon, alt = 0):
    """
    Производит перерасчет геодезических координат на геоиде WGS84 в координаты ECEF (геоцентрическая система координат)
    :param lat:     - Широта заданной точки
    :param lon:     - Долгота заданной точки
    :param lat:     - Высота над уровнем моря заданной точки, м
    :return:        - ECEF-координаты заданной точки
    """
    a = 6378137.0
    b = 6356752.3142
    f = (a - b) / a
    e = math.sqrt(f * (2 - f))
    Nfi = a / math.sqrt(1 - (e * math.sin(math.radians(lat))) ** 2.0)
    x = (Nfi +alt) * math.cos(math.radians(lat)) * math.cos(math.radians(lon))
    y = (Nfi + alt) * math.cos(math.radians(lat)) * math.sin(math.radians(lon))
    z = (Nfi * (1 - e ** 2.0) + alt) * math.sin(math.radians(lat))

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
    e = math.sqrt(f * (2 - f))
    lon = 180 * math.atan(y / x) / math.pi
    # Итерационный процесс вычисления широты и высоты над уровнем моря
    p = math.sqrt(x ** 2.0 + y ** 2.0)
    dfi = 1.0
    h   = 0
    fi  = math.atan(z /(p * (1 - e ** 2.0)))
    while math.fabs(dfi) > eps:
        N = a / math.sqrt(1 - (e * math.sin(fi)) ** 2.0)
        h = p / math.cos(fi) - N
        fi_N = math.atan(z / (p * (1 - (N * e ** 2.0) / (N + h))))
        dfi = fi_N -fi
        fi = fi_N

    lat = 180 * fi / math.pi
    alt = h

    return lat, lon, alt

#========================================================================================
# Процедуры вычисления радиуса наименьшей охватывающей окружности
#========================================================================================
def prepare_ECEF_points(points_WGS):
    """
    :param points_WGS:  Массив отсчетов в координатах WGS84
    :return:            Массив точек в координатах ECEF
    """
    points_ECEF = []
    for i in range(len(points_WGS)):
        x, y, z = WGS84_to_ECEF(points_WGS[i][0], points_WGS[i][1])
        points_ECEF.append((x, y, z))

    return points_ECEF

def prepare_Mercator_points(points_WGS):
    """
    :param points_WGS:  Массив отсчетов в координатах WGS84
    :return:            Массив точек в координатах проекции Меркатора
    """
    points_Mercator= []
    for i in range(len(points_WGS)):
        x, y = WGS84_to_Mercator(points_WGS[i][0], points_WGS[i][1])
        points_Mercator.append((x, y))

    return points_Mercator

def circle_2p(p1, p2):
    """
    :param p1:  Первая диаметральная точка окружности
    :param p2:  Вторая диаметральная точка окружности
    :return:    x, y, d - координаты центра и диаметр окружности
    """
    d = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
    x = (p2[0] + p1[0]) / 2.
    y = (p2[1] + p1[1]) / 2.

    return x, y, d

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

    D = math.hypot(Xc - x1, Yc - y1) * 2.

    return Xc, Yc, D

def is_in_circle(p, c, eps = 1e-12):
    """
    :param p:   Точка, для которой нужно произвести проверку попадания внутрь окружности
    :param C:   Окружность, относительно которой проверяется точка
    :param eps: Точность определения попадания точки внутрь окружности
    :return:    true, если попали, в противном случае - false
    """
    return math.hypot(p[0] - c[0], p[1] - c[1]) < (c[2] / 2. + eps)

def is_on_circle(p, c, eps = 0.1):
    """
    :param p:   Точка, для которой нужно произвести проверку попадания на окружность
    :param c:   Окружность, относительно которой проверяется точка
    :param eps: Точность определения попадания точки на окружность
    :return:    true, если попали, в противном случае - false
    """
    return math.fabs((math.hypot(p[0] - c[0], p[1] - c[1]) - c[2] / 2.)) <= eps


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

def min_enclosing_circle_wt(points, treshold = 0.9):
    """
    :param points:      Массив точек для фильтрации
    :param treshold:    Порог фильтрации
    :return:            Окружность, охватывающая заданный процент точек заданного множества
    """
    src_len = len(points)
    # Расчитываем минимальную охватывающую окружность для заданного множества
    c = min_enclosing_circle(points)
    # Находим точки, лежащин на окружности и исключаем их из рассмотрения
    points_new = []
    disp = 1.0
    while disp >= treshold:
        print c[0], "\t", c[1], "\t", c[2]
        for i in range(len(points)):
            if not is_on_circle(points[i], c):
                points_new.append(points[i])

        disp = len(points_new) / float(src_len)
        print "disp = ", disp, " new set len = ", len(points_new)
        c = min_enclosing_circle(points_new)

        points = points_new
        points_new = []

    print src_len, len(points_new)
    return c

if __name__ == "__main__":
    points = read_coordinates("bezP.txt")
    Mercator_points = prepare_Mercator_points(points)

    x_draw = []
    y_draw = []

    for i in range(len(Mercator_points)):
        x_draw.append(Mercator_points[i][0])
        y_draw.append(Mercator_points[i][1])

    c = min_enclosing_circle(Mercator_points)

    fig, ax = plt.subplots(figsize=(6,6))
    circle = plt.Circle((c[0], c[1]), c[2] * 0.5, facecolor='none', edgecolor='g', linewidth=1, alpha= 0.5)
    ax.add_patch(circle)

    plt.plot(x_draw, y_draw, 'yx:', alpha=0.8)
    plt.axis('equal')
    plt.show()

    lat, lon = Mercator_to_WGS84(c[0], c[1])
    print "======================================================================="
    print "Minimal enclosing circle for given points set is"
    print "======================================================================="
    print " center X   = ", c[0]
    print " center Y   = ", c[1]
    print " center LAT = ", lat
    print " center LON = ", lon
    print " diameter   = ", c[2]
