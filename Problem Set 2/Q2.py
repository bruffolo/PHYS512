import numpy as np
from matplotlib import pyplot as plt

plot = True

t,flux,u = np.loadtxt("229614158_PDCSAP_SC6.txt",delimiter=',',unpack = True)

# Identify the spike
spike = max(flux)
index = np.argwhere(flux == spike)[0][0]
t_spike = t[index]

#cutoff = index+50
#flux,t = flux[index-20:cutoff],t[index-20:cutoff]

# Model function
def model(t,A,tau,b):
    out = A * np.exp(-tau*(t-t_spike)) + b
    out[(t - t_spike) < 0] = b
    return out

# Derivative calculation
def calc_model(t,p):
    y = model(t,*p)
    grad=np.zeros([t.size,p.size])
    #now differentiate w.r.t. all the parameters
    for i in range(0,len(t)):
        if( (t[i]-t_spike) < 0):
            grad[i,0]= 0
            grad[i,1]= 0
            grad[i,2] = 1.0
        else:
            grad[i,0]= np.exp(-p[1]*(t[i]-t_spike))
            grad[i,1]= -p[0]*(t[i]-t_spike)*np.exp(-p[1]*(t[i]-t_spike))
            grad[i,2] = 1.0
    return y,grad


# Starting guess
p0 = np.array([max(flux)-1,50,1])

p=p0.copy()
y,grad=calc_model(t,p)

# Newton's method loop
for j in range(5):
    y,grad=calc_model(t,p)
    r=flux-y
    err=(r**2).sum()
    r=np.matrix(r).transpose()
    grad=np.matrix(grad)

    lhs=grad.transpose()*grad
    rhs=grad.transpose()*r
    dp=np.linalg.inv(lhs)*(rhs)
    for jj in range(p.size):
        p[jj]=p[jj]+dp[jj]
    #print(p,err)

print("Best Fit Parameters:\nA = ",p[0],"\ntau = ",p[1],"\nb = ",p[2])

if(plot):
    # Main figure
    fig1 = plt.figure(1) # Add main frame
    frame1=fig1.add_axes((.13,.3,.83,.6))

    plt.plot(t,flux,'o',color = "blue",markersize = 2, label = "Data")
    plt.plot(t,model(t,*p), color = 'green', label = "Optimized Fit")
    plt.plot(t,model(t,*p0), color = 'Orange', label = "Guessed Fit")
    frame1.set_xticklabels([]) # remove x-tic labels for the first frame
    frame1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    plt.ylabel("Flux")
    plt.grid()
    #plt.xlim([min(t),max(t)]) 
    plt.xlim([t_spike-.3,max(t)-3.2])
    plt.legend()

    # Residual figure
    frame2=fig1.add_axes((.13,.09,.83,.20))  # Add frame for Residuals
    res = flux-model(t,*p)
    plt.plot(t,res,'or', color = "Green",markersize = 4) # plot the residuals
    res = flux-model(t,*p0)
    plt.plot(t,res,'or', color = "Orange",markersize = 4) # plot the residuals
    plt.xlabel("Time")
    plt.ylabel("Residuals")
    #plt.xlim([min(t),max(t)])
    plt.xlim([t_spike-.3,max(t)-3.2])
    plt.grid()

    # Save the plot
    plt.savefig("Q2_a_b")
    plt.show()

