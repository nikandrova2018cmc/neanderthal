import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

def comparaison(reality,calculated):
    total_gen = sum(list(map(sum,reality)))
    errors = []
    not_inside = 0
    not_inside1 = 0
    for calc in calculated:
        flag = True
        for real in reality:
            start_c,end_c,start_r,end_r = calc[0],calc[1],real[0],real[1]
            maximum_margin = int(end_c-start_c)*0.75
            if abs(start_c-start_r) < maximum_margin and abs(end_c-end_r) < maximum_margin:
                error = abs(start_c-start_r) + abs(end_c-end_r)
                errors.append(error)
                flag = False
                break
            error = abs(start_c-end_c)
            errors.append(error)
        if flag:
            not_inside += 1
        not_inside1  +=1
    rate_I_NI = not_inside/not_inside1
    error_rate = (sum(errors)/total_gen)*100
    return error_rate,rate_I_NI

def graph(X,Y):
    x,y = np.array(X),np.array(Y)
    X_Y_Spline = make_interp_spline(x, y)


    X_ = np.linspace(x.min(), x.max(), 500)
    Y_ = X_Y_Spline(X_)

    plt.plot(X_, Y_)
    plt.title("Relation between the number of African genomes with the incorrect intervals ratio")
    plt.xlabel("Number of Afrian genomes")
    plt.ylabel("Incorrect intervals ratio")
    plt.show()
