import re,sys
import matplotlib.pyplot as plt
import numpy as np

def main():
    with open(sys.argv[1]) as fp:
        data = fp.readlines()
    xTitle,yTitle = data[0].rstrip().split(',')    
    X = []
    Y = []
    for line in data[1:]:
        x,y = line.rstrip().split(',') 
        X.append(float(x))
        Y.append(float(y))
    X = np.array(X)
    Y = np.array(Y)
    plt.plot(X,Y)
    plt.xlabel(xTitle)
    plt.ylabel(yTitle)
    plt.title("DOA estimation based on MUSIC algorithm")
    plt.show()
if __name__ == '__main__':
    main()
