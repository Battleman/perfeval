import numpy as np
import matplotlib.pyplot as plt
from helpers import ci_median, ci_mean, pi_on_mean, pi_order_stat
import pylab 
import scipy.stats as stats

def normal_plot(n):
    values = np.random.normal(size=n)
    plt.hist(values)
    plt.xlim(-3,3)
    plt.title("{} samples".format(n))

def plot_distributions(points1, points2, legend1, legend2, bins):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    
    ax1.hist(points1, bins=bins);
    ax1.set_title(legend1)
    
    ax2.hist(points2, bins=bins);
    ax2.set_title(legend2)
    
    plt.tight_layout()

    
def qq_plots(points1, points2, legend1, legend2):
    fig = plt.figure(figsize=(13,5))
    
    fig.add_subplot(1,2,1)
    stats.probplot(points1, dist="norm", plot=pylab)
    plt.title("QQ-plot " + legend1)
    
    fig.add_subplot(1,2,2)
    stats.probplot(points2, dist="norm", plot=pylab)
    plt.title("QQ-plot " + legend2)
    
    pylab.show()
    
    
def plot_boxplots_ci(X, Y, j, k, n):
    
    ci_median_X, ci_median_Y = ci_median(X,j,k), ci_median(Y,j,k)
    ci_mean_X, ci_mean_Y = ci_mean(X), ci_mean(Y)
    
    fig, ax= plt.subplots(1,2, figsize=(20,25)) 

    #################
    ### Event average
    #################
    ax[0].boxplot(X, 1,"")

    # plot mean and median
    ax[0].plot([0.85, 1.15], [np.mean(X)]*2, label="Mean", color='orange')
    ax[0].plot([0.85, 1.15], [np.median(X)]*2, label="Median", color='cyan')
    
    #plot CI for mean
    ax[0].plot([0.75, 1.25], [ci_mean_X[0], ci_mean_X[0]], 
               linestyle="dashed", color='red', linewidth=1, label="CI for mean")
    ax[0].plot([0.75, 1.25], [ci_mean_X[1], ci_mean_X[1]], 
               linestyle="dashed", color='red', linewidth=1)

    #plot CI for median
    ax[0].plot([0.75, 1.25], [ci_median_X[0], ci_median_X[0]], 
               linestyle="dashed", color='blue', linewidth=1, label="CI for median")
    ax[0].plot([0.75, 1.25], [ci_median_X[1], ci_median_X[1]], 
               linestyle="dashed", color='blue', linewidth=1)

    #annotate CI for mean
    ax[0].annotate("{:.4f}".format(ci_mean_X[0]), 
                   xy=(0.75, ci_mean_X[0]), xytext=(0.7, ci_mean_X[0]-0.05), 
                   arrowprops=dict(facecolor='red', shrink=0.05,shrinkA=0, width=0.5))
    ax[0].annotate("{:.4f}".format(ci_mean_X[1]), 
                   xy=(0.75, ci_mean_X[1]), xytext=(0.7, ci_mean_X[1]+0.05), 
                   arrowprops=dict(facecolor='red', shrink=0.05,shrinkA=0, width=0.5))

    #annotate CI for median
    ax[0].annotate("{:.4f}".format(ci_median_X[0]), 
                   xy=(0.9, ci_median_X[0]), xytext=(0.9, ci_median_X[0]-0.05), 
                   arrowprops=dict(facecolor='blue', shrink=0.05, shrinkA=0,width=0.5))
    ax[0].annotate("{:.4f}".format(ci_median_X[1]), 
                   xy=(0.9, ci_median_X[1]), xytext=(0.9, ci_median_X[1]+0.05), 
                   arrowprops=dict(facecolor='blue', shrink=0.05, shrinkA=0,width=0.5))

    #annotate median and mean
    ax[0].annotate("Median\n{:.4f}".format(np.median(X)), 
                   xy=(1, np.median(X)), xytext=(1.25, np.median(X)-0.05), 
                   arrowprops=dict(facecolor='cyan', shrink=0.05, shrinkA=0,width=0.5))
    ax[0].annotate("Mean:\n{:.4f}".format(np.mean(X)), 
                   xy=(1, np.mean(X)), xytext=(1.25, np.mean(X)+0.05), 
                   arrowprops=dict(facecolor='orange', shrink=0.05, shrinkA=0,width=0.5))

    ax[0].set_title("Event Average N={}".format(n))
    ax[0].legend()
    
    
    ###############
    #### TIME AVERAGE
    ###############
    
    #make boxplot
    ax[1].boxplot(Y, 1,"")
    
    ax[1].plot([0.75, 1.25], [np.mean(Y)]*2, label="Mean", c="orange")
    ax[1].plot([0.75, 1.25], [np.median(Y)]*2, label="Median", c="cyan")
    
    #plot CI for mean
    ax[1].plot([0.75, 1.25], [ci_mean_Y[0], ci_mean_Y[0]], 
               linestyle="dashed", color='red', linewidth=1, label="CI for mean")
    ax[1].plot([0.75, 1.25], [ci_mean_Y[1], ci_mean_Y[1]], 
               linestyle="dashed", color='red', linewidth=1)
    
    #plot CI for median
    ax[1].plot([0.85, 1.15], [ci_median_Y[0], ci_median_Y[0]], 
               linestyle="dashed", color='blue', linewidth=1, label="CI for median")
    ax[1].plot([0.85, 1.15], [ci_median_Y[1], ci_median_Y[1]], 
               linestyle="dashed", color='blue', linewidth=1)

    #annotate CI for mean
    ax[1].annotate("{:.4f}".format(ci_mean_Y[0]), 
                   xy=(0.75, ci_mean_Y[0]), xytext=(0.6, ci_mean_Y[0]-0.05), 
                   arrowprops=dict(facecolor='red', shrink=0.05, shrinkA=0, width=0.5))
    ax[1].annotate("{:.4f}".format(ci_mean_Y[1]), 
                   xy=(0.75, ci_mean_Y[1]), xytext=(0.6, ci_mean_Y[1]+0.05), 
                   arrowprops=dict(facecolor='red', shrink=0.05, shrinkA=0, width=0.5))
    #annotate  CI for median
    ax[1].annotate("{:.4f}".format(ci_median_Y[0]), 
                   xy=(0.9, ci_median_Y[0]), xytext=(0.6, ci_median_Y[0]-0.05), 
                   arrowprops=dict(facecolor='blue', shrink=0.05, shrinkA=0, width=0.5))
    ax[1].annotate("{:.4f}".format(ci_median_Y[1]), 
                   xy=(0.9, ci_median_Y[1]), xytext=(0.6, ci_median_Y[1]+0.05),
                   arrowprops=dict(facecolor='blue', shrink=0.05, shrinkA=0,width=0.5))

    #annotate median and mean
    ax[1].annotate("Median:\n{:.4f}".format(np.median(Y)), 
                   xy=(1, np.median(Y)), xytext=(1.25, np.median(Y)-0.05), 
                   arrowprops=dict(facecolor='cyan', shrink=0.05, shrinkA=0, width=0.5))
    ax[1].annotate("Mean:\n{:.4f}".format(np.mean(Y)), 
                   xy=(1, np.mean(Y)), xytext=(1.25, np.mean(Y)+0.05), 
                   arrowprops=dict(facecolor='orange', shrink=0.05, shrinkA=0, width=0.5))

    ax[1].set_title("Time Average N={}".format(n))
    ax[1].legend()

    plt.show()
    
def plot_boxplots_pi(Y):

    pi_mean = pi_on_mean(Y)
    pi_ostat = pi_order_stat(Y, 0.95)
       
    plt.figure(figsize=(15, 10))
    #make boxplot
    plt.boxplot(Y, 1,"")
    

    plt.plot([0.75, 1.25], [pi_mean[0]]*2, label="PI based on mean", c="blue")
    plt.plot([0.75, 1.25], [pi_mean[1]]*2, c="blue")
    plt.annotate("{:.4f}".format(pi_mean[0]), xy=(1.2,pi_mean[0]+0.01))
    plt.annotate("{:.4f}".format(pi_mean[1]), xy=(1.2,pi_mean[1]+0.01))
    
    if pi_ostat[0] > 0 and pi_ostat[1] > 0: # if method worked, plot
        plt.plot([0.75, 1.25], [pi_ostat[0]]*2, label="PI based on order stat", c="red")
        plt.plot([0.75, 1.25], [pi_ostat[1]]*2, c="red")
        plt.annotate("{:.4f}".format(pi_ostat[0]), xy=(1.2,pi_ostat[0]+0.01))
        plt.annotate("{:.4f}".format(pi_ostat[1]), xy=(1.2,pi_ostat[1]+0.01))

    plt.title("Time Average N={}".format(len(Y)))

    plt.legend()
    plt.show()