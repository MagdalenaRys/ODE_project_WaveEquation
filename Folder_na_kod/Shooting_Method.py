import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import math


def modifiedEulerMethod(function, initial_cond, parting, division=10**3):
    """
    Funkcja rozwiązuje równanie różniczkowe zadane w postaci: x''(t) = lamb * f(t, x, x') z określonymi warunkami początkowymi
    na danym przedziale np. na przedziale (0, 12.5) z następującymi warunkami: x(0) = 0, x'(0) = 2
    :param lamb: lambda
    :param function: funkcja zależąca od zmiennej; f(t, x, x')
    :param initial_cond: warunek początkowy; wartości funkcji w punkcie parting[0]
    :param parting: przedział, dla którego otrzymamy wartości funkcji po rozwiązaniu równania różniczkowego
    :param division: podział przedziału na division części
    :return: funkcja zwraca dwie tablice - pierwsza z nich to podział odcinka parting na division części, druga tablica
    z kolei zwraca wartości szukanej funkcji dla odpowiadających jej punktów z pierwszej tablicy
    """
    if parting[0] <= initial_cond[0][0] <= parting[1]:
        left = initial_cond[0][0]-parting[0]
        right = parting[1]-initial_cond[0][0]
        left_div = int(left/(left + right) * division)
        right_div = int(right/(left + right) * division)
        x_list = left_div*[0]
        y_list = left_div * [0]
        w_list = left_div*[0]
        w = initial_cond[0][0]
        c = left_div - 1
        x_list.append(initial_cond[0][1])
        y_list.append(initial_cond[1][1])
        w_list.append(w)
        if left_div > 0:
            left = left/left_div
            while c >= 0:
                x_list[c] = x_list[c+1] - left * y_list[c+1]
                y_list[c] = y_list[c+1] - left * function(w, x_list[c+1], y_list[c+1])
                w -= left
                w_list[c] = w
                c -= 1
        if right_div > 0:
            w = initial_cond[0][0]
            right = right / right_div
            while w <= parting[1] - right:
                x_list.append(x_list[-1] + right * y_list[-1])
                y_list.append(y_list[-1] + right * function(w, x_list[-2], y_list[-1]))
                w += right
                w_list.append(w)
        return w_list, x_list


def shootingMethod(function, boundary_cond, parting, a_part, division=10**3, epsilon=0.01):
    """
    :param function:
    :param boundary_cond:
    :param parting:
    :param division:
    :param epsilon:
    :return:
    """
    a = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], a_part[0])], parting, division)[1][-1] - boundary_cond[1][1]
    b = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], a_part[1])], parting, division)[1][-1] - boundary_cond[1][1]
    if a*b >= 0:
        raise Exception("coś nie tak z przedziałem XD")
    while True:
        c = (a_part[1]+a_part[0])/2
        p = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], c)], parting, division)
        q = p[1][-1] - boundary_cond[1][1]
        if abs(q) <= epsilon:
            return p
        else:
            if a*q > 0:
                a_part = (c, a_part[1])
            elif a*q < 0:
                a_part = (a_part[0], c)


def funP(x, P, Y, lamb=4.0):
    return lamb * P


def funQ(t, Q, Y, lamb=2.0, c=2.0):
    return -lamb**2 * c**2 * Q


def f_func(x):
    return -0.1*(x-5)**2 + 2.5


def g_func(x):
    return -0.1*(x-5)**2 + 2.5


P = shootingMethod(lambda w, x, y: funP(w, x, y, -math.pi/10), [(0.0, 0.0), (10.0, 0.0)], (-2.0, 10.0), (-5.0, 10.0))


def initialConditionsToQ(x, f_func, g_func, fun_p_lists, eps=0.01):
    a = 0
    for i in range(len(fun_p_lists[0])):
        if abs(fun_p_lists[0][i] - x) <= eps:
            a = fun_p_lists[1][i]
            break
    return [(0, f_func(x)/a), (0, g_func(x)/a)]


Q = modifiedEulerMethod(funQ, initialConditionsToQ(1, f_func, g_func, P), (-2.0, 10.0))


U = np.multiply(P[1], Q[1])
plt.plot(P[0], U)
plt.show()
