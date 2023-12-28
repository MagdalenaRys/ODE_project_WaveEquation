import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def modifiedEulerMethod(function, initial_cond, parting, division=10**2):
    """
    :param function:
    :param initial_cond:
    :param parting:
    :param division:
    :return:
    """

    delta_h = (parting[1]-parting[0])/division
    x_list = [initial_cond[0]]
    y_list = [initial_cond[1]]
    w_list = [parting[0]]
    w = parting[0]
    while w <= parting[1] - delta_h:
        x_list.append(x_list[-1] + delta_h * y_list[-1])
        y_list.append(y_list[-1] + delta_h * function(w, x_list[-2], y_list[-1]))
        w += delta_h
        w_list.append(w)
    return x_list, w_list


def shootingMethod(function, boundary_cond, parting, a_part, division=10**2, epsilon=0.01):
    p = modifiedEulerMethod(function, (boundary_cond[0], (a_part[1]-a_part[0])/2), parting, division)
    if abs(p[-1] - boundary_cond[1]) <= epsilon:
        return p
    else:
        pass


def funX(x, X, Y, lamb=2):
    return -lamb**2 * X


def funT(t, T, Y, lamb=2, c=2):
    return -lamb**2 * c**2 * T


X = modifiedEulerMethod(funX, (0, 1), (0, 3.1415926))

print(X[0])

plt.plot(X[1], X[0])
plt.show()
