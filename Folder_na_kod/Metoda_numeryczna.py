import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np


def modifiedEulerMethod(function, initial_cond, parting, division=10 ** 5):
    """
    Funkcja rozwiązuje równanie różniczkowe zadane w postaci: x''(t) = f(t, x, x') z określonymi warunkami
    początkowymi na danym przedziale np. na przedziale (0, 12.5) z następującymi warunkami: x(0) = 0,
    x'(0) = 2 :param function: funkcja f(t, x, x') :param initial_cond: warunek początkowy :param parting: przedział,
    dla którego otrzymamy wartości funkcji po rozwiązaniu równania różniczkowego :param division: podział przedziału na
    division części:return: funkcja zwraca dwie tablice - pierwsza z nich to podział odcinka parting na division części,
    druga tablica z kolei zwraca wartości szukanej funkcji dla odpowiadających jej punktów z pierwszej tablicy
    """
    # Sprawdzenie, czy warunek początkowy jest zadany dla punktu z przedziału parting
    if parting[0] <= initial_cond[0][0] <= parting[1]:
        # podział odcinka parting względem punktu
        left = initial_cond[0][0] - parting[0]
        right = parting[1] - initial_cond[0][0]
        # podział na równe części
        left_div = int(left / (left + right) * division)
        right_div = int(right / (left + right) * division)
        # tworzenie pustych tablicy na wartości. y_list to tablica przechowująca wartości x', potrzebne do obliczenia x
        x_list = left_div * [0]
        y_list = left_div * [0]
        w_list = left_div * [0]
        w = initial_cond[0][0]  # ustawienie wartości w na zadaną w warunku początkowym
        c = left_div - 1        # ustawienie licznika
        x_list.append(initial_cond[0][1])
        y_list.append(initial_cond[1][1])
        w_list.append(w)
        if left_div != 0:       # sprawdzenie, czy warunek początkowy nie był zadany na lewym krańcu przedziału
            left = left / left_div      # ustalenie podziału odcinka
            while c >= 0:
                # korzystamy z odwróconej metody Eulera
                x_list[c] = x_list[c + 1] - left * y_list[c + 1]
                y_list[c] = y_list[c + 1] - left * function(w, x_list[c + 1], y_list[c + 1])
                w -= left
                w_list[c] = w
                c -= 1
        if right_div != 0:      # sprawdzenie, czy warunek początkowy nie był zadany na prawym krańcu przedziału
            w = initial_cond[0][0]
            right = right / right_div   # ustalenie podziału odcinka
            while w <= parting[1] - right:
                # korzystamy z jawnej metody Eulera
                x_list.append(x_list[-1] + right * y_list[-1])
                y_list.append(y_list[-1] + right * function(w, x_list[-2], y_list[-1]))
                w += right
                w_list.append(w)
        return w_list, x_list


def shootingMethod(function, boundary_cond, parting, a_part, division=10 ** 5 +1, epsilon=0.001):
    """
    Funkcja rozwiązuje zagadnienie brzegowe dla równania różniczkowego zadanego w postaci: x''(t) = f(t, x, x')
    :param function: funkcja f(t, x, x')
    :param boundary_cond: warunkek brzegowy zadany w postaci: [(t_0, x_0), (t_1, x_1)]
    :param parting: przedział, dla jakiego chcemy mieć tablicę wartości funkcji x
    :param a_part: przedział, z którego będziemy szukać odpowiedniego parametru c
    :param division: podział przedziału parting
    :param epsilon: parametr, który określa nam dokładność, z jaką chcemy rozwiązać warunek brzegowy
    :return: c, gdzie c jest wartością x'(t_0)
    """
    # sprawdzanie wartości na krańcach przedziału a_part. Szukanie c metodą bisekcji
    a = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], a_part[0])], parting, division)[1][-1] - \
        boundary_cond[1][1]
    b = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], a_part[1])], parting, division)[1][-1] - \
        boundary_cond[1][1]
    if a * b >= 0:
        raise Exception("Dobierz inny przedział a_part")
    while True:
        c = (a_part[1] + a_part[0]) / 2
        p = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], c)], parting, division)
        q = p[1][-1] - boundary_cond[1][1]
        if abs(q) <= epsilon:
            return c
        else:
            if a * q > 0:
                a_part = (c, a_part[1])
            elif a * q < 0:
                a_part = (a_part[0], c)


def analyticMethodSODE(list_of_coefs, initial_cond, parting, div=10**5):
    """
    Funkcja rozwiązuje analitycznie równanie różniczkowe jednorodne zadane w postaci: a * y''(t) + b * y'(t) + c * y(t) = 0
    :param list_of_coefs: współczynniki w kolejności: b, c
    :param initial_cond: warunki początkowe zadane w postaci: [(t0, y(t0) = y1), (t0, y'(t0) = y2)]
    :param parting: przedział, dla jakiego chcemy otrzymać wartości funkcji
    :return:  arguments_list, values_list - listy argumentów i wartości
    """
    if initial_cond[0][0] != initial_cond[1][0]:
        raise ValueError('Pochodna jest w innym punkcie!')
    a, b, c = list_of_coefs
    if a == 0:
        raise ValueError('Źle dobrane współczynniki! Współczynnik przy drugiej pochodnej musi być różny od zera!')
    values_list = []
    arguments_list = []
    # obliczanie wielomianu charakterystycznego
    delta_r = b ** 2 - 4 * c
    if delta_r == 0:
        r = -b/(2*a)
        expon = math.exp(r*initial_cond[0][0])
        # wyznaczanie współczynników z układu równań
        A = np.array([[expon, initial_cond[1][0] * expon], [r * expon, r * initial_cond[1][0] * expon + expon]])
        B = np.array([initial_cond[0][1], initial_cond[1][1]])
        C_1, C_2 = np.linalg.inv(A).dot(B)
        t = parting[0]
        delta_t = (parting[1] - parting[0])/div
        for i in range(div):
            values_list.append(C_1*math.exp(r*t) + C_2*math.exp(r*t)*t)
            arguments_list.append(t)
            t += delta_t
        return arguments_list, values_list
    elif delta_r > 0:
        r1 = (- math.sqrt(delta_r) - b)/(2 * a)
        r2 = (math.sqrt(delta_r) - b)/(2 * a)
        exp1 = math.exp(r1 * initial_cond[0][0])
        exp2 = math.exp(r2 * initial_cond[1][0])
        # wyznaczanie współczynników z układu równań
        A = np.array([[exp1, exp2], [r1 * exp1, r2 * exp2]])
        B = np.array([initial_cond[0][1], initial_cond[1][1]])
        C_1, C_2 = np.linalg.inv(A).dot(B)
        t = parting[0]
        delta_t = (parting[1] - parting[0]) / div
        for i in range(div):
            values_list.append(C_1 * math.exp(r1 * t) + C_2 * math.exp(r2 * t))
            arguments_list.append(t)
            t += delta_t
        return arguments_list, values_list
    else:
        re = -b/(2*a)
        im = math.sqrt(-delta_r)/(2*a)
        exp = math.exp(re * initial_cond[0][0])
        sin = math.sin(im * initial_cond[0][0])
        cos = math.cos(im * initial_cond[0][0])
        # wyznaczanie współczynników z układu równań
        A = np.array([[exp * sin, exp * cos], [re*exp*sin + exp*cos*im, re*exp*cos - exp*sin*im]])
        B = np.array([initial_cond[0][1], initial_cond[1][1]])
        C_1, C_2 = np.linalg.inv(A).dot(B)
        t = parting[0]
        delta_t = (parting[1] - parting[0]) / div
        for i in range(div):
            values_list.append(math.exp(re*t) * (C_1 * math.sin(im*t) + C_2 * math.cos(im*t)))
            arguments_list.append(t)
            t += delta_t
        return arguments_list, values_list


def funP(x, P, Y, n=4.0, l=10.0):
    return - (math.pi * n / l) ** 2 * P


def funQ(t, Q, Y, n=2.0, l=10.0, c=2.0):
    return - (math.pi * n / l) ** 2 * c ** 2 * Q


def f_func(x):
    return -0.1 * (x - 5) ** 2 + 2.5


def g_func(x):
    return -0.1 * (x - 5) ** 2 + 2.5


def initialConditionsToQ(x, f_func, g_func, fun_p_lists, eps=0.0001):
    """
    Funkcja wyznacza nam warunki początkowe dla fukcji Q(t) w punkcie t=0, korzystając z rozwiązania funkcji P.
    :param x: dowolny x, dla którego znamy wartość funkcji f, g i P. Jest on potrzebny do obliczania ilorazu.
    :param f_func: funkcja f(x)
    :param g_func: funkcja g(x)
    :param fun_p_lists: lista argumentów i wartości funkcji P
    :param eps: dokładność warunków początkowych
    :return: warunek początkowy zadany w postaci: [(0, Q(0)), (0, Q'(0))]
    """
    for i in range(len(fun_p_lists[0])):
        if abs(fun_p_lists[0][i] - x) <= eps:
            a = fun_p_lists[1][i]
            return [(0, f_func(x) / a), (0, g_func(x) / a)]
    raise 'No value for given argument'


if __name__ == '__main__':
    n = 3
    l = 10.0
    c = 1.0
    d = shootingMethod(lambda w, x, y: funP(w, x, y, n, l), [(0.0, 0.0), (l, 0.0)], (0.0, l), (-3.0, 2.0))
    P = modifiedEulerMethod(lambda w, x, y: funQ(w, x, y, n, l, 1.0), [(0.0, 0.0), (0.0, d)], (0.0, l))

    P_an = analyticMethodSODE([1, 0, (math.pi * n / l) ** 2], [(0.0, 0.0), (0.0, d)], (0.0, l))

    # plt.plot(P[0], P[1])
    # plt.plot(P_an[0], P_an[1])
    # plt.show()

    Q = modifiedEulerMethod(lambda w, x, y: funQ(w, x, y, n, l, 1.0), initialConditionsToQ(3, f_func, g_func, P),
                            (0.0, l))

    Q_an = analyticMethodSODE([1, 0, (math.pi * n / l) ** 2 * c ** 2], initialConditionsToQ(3, f_func, g_func, P), (0.0, l))

    X = np.array(P[1]) * Q[1][7000]
    X_an = np.array(P_an[1]) * Q_an[1][7000]

    plt.plot(P[0], X)
    plt.plot(P_an[0], X_an)
    plt.show()
