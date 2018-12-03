import matplotlib.pyplot as plt
def creat_plot(xaxis, yaxis, type, xlabel, ylabel, title):
    plt.plot(xaxis, yaxis, str(type))
    plt.xlabel(str(xlabel))
    plt.ylabel(str(ylabel))
    plt.title(str(title))
    return plt