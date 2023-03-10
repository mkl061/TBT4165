import matplotlib.pyplot as plt
import numpy as np

def multi_plot(plot_type, lst_of_xs, lst_of_ys):
    figure, axis = plt.subplots(1, len(lst_of_ys))
    
    plot_types = ["scatter", "plot"]
    if plot_type not in plot_types:
        raise Exception(f"Method {plot_type} not implemented for axis")
        
    for i in range(len(lst_of_ys)):
        eval(f"axis[{i}].{plot_type}(lst_of_xs[{i}],lst_of_ys[{i}])")
    
    return axis
    #plt.show()

y1 = np.random.normal(5,2,10)
y2 = np.random.normal(1,4,10)
y3 = np.random.normal(6,5,10)
x = range(len(y1))

ys = [y1,y2,y3]
#plt.scatter(x,y2)
#plt.show()

multi_plot("scatter", [x,x,x], ys)
