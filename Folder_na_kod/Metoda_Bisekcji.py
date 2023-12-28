import math


def Value_Of_Polynom(list_of_coeffs, value):
    sum = 0
    for i in range(len(list_of_coeffs)):
        sum += list_of_coeffs[i]*value**i
    return sum


def DarbouxMethod(function, parting, accuracy):
    x1 = (parting[0]+parting[1])/2
    f_x1 = function(x1)
    if abs(f_x1) <= accuracy:
        return x1, f_x1
    elif f_x1 * function(parting[0]) < 0:
        return DarbouxMethod(function, [parting[0], x1], accuracy)
    elif f_x1 * function(parting[1]) < 0:
        return DarbouxMethod(function, [x1, parting[1]], accuracy)

