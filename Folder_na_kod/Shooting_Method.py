import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation
from IPython import display


def modifiedEulerMethod(function, initial_cond, parting, division=10 ** 5):
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
        left = initial_cond[0][0] - parting[0]
        right = parting[1] - initial_cond[0][0]
        left_div = int(left / (left + right) * division)
        right_div = int(right / (left + right) * division)
        x_list = left_div * [0]
        y_list = left_div * [0]
        w_list = left_div * [0]
        w = initial_cond[0][0]
        c = left_div - 1
        x_list.append(initial_cond[0][1])
        y_list.append(initial_cond[1][1])
        w_list.append(w)
        if left_div != 0:
            left = left / left_div
            while c >= 0:
                x_list[c] = x_list[c + 1] - left * y_list[c + 1]
                y_list[c] = y_list[c + 1] - left * function(w, x_list[c + 1], y_list[c + 1])
                w -= left
                w_list[c] = w
                c -= 1
        if right_div != 0:
            w = initial_cond[0][0]
            right = right / right_div
            while w <= parting[1] - right:
                x_list.append(x_list[-1] + right * y_list[-1])
                y_list.append(y_list[-1] + right * function(w, x_list[-2], y_list[-1]))
                w += right
                w_list.append(w)
        return w_list, x_list


def shootingMethod(function, boundary_cond, parting, a_part, division=10 ** 5, epsilon=0.001):
    """
    :param function:
    :param boundary_cond:
    :param parting:
    :param division:
    :param epsilon:
    :return:
    """
    a = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], a_part[0])], parting, division)[1][-1] - \
        boundary_cond[1][1]
    b = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], a_part[1])], parting, division)[1][-1] - \
        boundary_cond[1][1]
    if a * b >= 0:
        raise Exception("coś nie tak z przedziałem XD")
    while True:
        c = (a_part[1] + a_part[0]) / 2
        p = modifiedEulerMethod(function, [boundary_cond[0], (boundary_cond[0][0], c)], parting, division)
        q = p[1][-1] - boundary_cond[1][1]
        if abs(q) <= epsilon:
            return p
        else:
            if a * q > 0:
                a_part = (c, a_part[1])
            elif a * q < 0:
                a_part = (a_part[0], c)


def funP(x, P, Y, lamb=4.0):
    return - lamb ** 2 * P


def funQ(t, Q, Y, lamb=2.0, c=2.0):
    return - lamb ** 2 * c ** 2 * Q


def f_func(x):
    return -0.1 * (x - 5) ** 2 + 2.5


def g_func(x):
    return -0.1 * (x - 5) ** 2 + 2.5


def initialConditionsToQ(x, f_func, g_func, fun_p_lists, eps=0.0001):
    for i in range(len(fun_p_lists[0])):
        if abs(fun_p_lists[0][i] - x) <= eps:
            a = fun_p_lists[1][i]
            return [(0, f_func(x) / a), (0, g_func(x) / a)]
    raise 'No value for given argument'


if __name__ == '__main__':
    P = shootingMethod(lambda w, x, y: funP(w, x, y, math.pi / 10), [(0.0, 0.0), (10.0, 0.0)], (0.0, 10.0), (-3.0, 2.0))

    # p = modifiedEulerMethod(lambda w, x, y: funP(w, x, y, math.pi), [(0.0, 0.0), (0.0, 1.0)], (0.0, 10.0))

    Q = modifiedEulerMethod(lambda w, x, y: funQ(w, x, y, math.pi / 10, 1.0), initialConditionsToQ(3, f_func, g_func, P),
                            (0.0, 10.0))

    # U0 = np.array(P[1])*Q[1][0]
    # U1 = np.array(P[1])*Q[1][1000]
    # U2 = np.array(P[1])*Q[1][2000]
    # U3 = np.array(P[1])*Q[1][3000]
    # U4 = np.array(P[1])*Q[1][90000]
    #
    # plt.plot(P[0], U0)
    # plt.plot(P[0], U1)
    # plt.plot(P[0], U2)
    # plt.plot(P[0], U3)
    # plt.plot(P[0], U4)
    # plt.show()
    x = P[0]


    def AnimationFunction(frame):
        # setting y according to frame
        # number and + x. It's logic
        y = np.cos(x + 2 * np.pi * frame / 100)

        # line is set with new values of x and y
        line_plotted.set_data((x, y))


    Figure = plt.figure()
    # creating a plot
    lines_plotted = plt.plot([])

    # putting limits on x axis since
    # it is a trigonometry function
    # (0,2∏)
    line_plotted = lines_plotted[0]
    # plt.xlim(0.0, 10.0)
    anim_created = FuncAnimation(Figure, AnimationFunction, frames=100, interval=25)
    video = anim_created.to_html5_video()
    html = display.HTML(video)
    display.display(html)

    # good practice to close the plt object.
    plt.close()
