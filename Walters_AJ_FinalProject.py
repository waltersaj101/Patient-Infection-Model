import numpy as np
import matplotlib.pyplot as plt

nSens = 10**9              	# initial number of sensitive bacteria
nRes = 0                    # initial number of resistant bacteria
nImm = 1                    # initial number of recruited immune cells
cAB1 = 0                    # initial concentration of isonizid (INH) antibiotic
cAB2 = 0                    # initial concentration of pyrazinamide (PZA) antibiotic
cAB3 = 0                    # initial concentration of rifampicin (RIF) antibiotic
betaS = 0.8                 # growth rate of sensitive bacteria, (per day)
betaR = 0.4                 # growth rate of resistant bacteria, (per day: range 0.4-0.1)
k = 0.6                     # growth rate of immune cells, (per day)
eta = 0.3                   # rate of bacteria destroyed by immune cells, (per day)
omega = 1                   # carrying capacity of immune cells, ratio of bacteria present to immune cells (immune cells)
T = 10**9                   # carrying capacity of bacteria, (bacteria)
z = 1						# scaling factor for T-Cell therapy
ab = np.zeros((3, 4))       # antibiotics, rows = anitbiotics, columns = name, alpha, d, delta, mu
ab[0][0] = 10**-6           # alphaBar, mutation rate of isonaizid (INH) antibiotic (mutations*generation)
ab[0][1] = 0.0039           # dBar, elimination rate of sensitive bacteria due to INH (per day)
ab[0][2] = 5                # delta, daily dose of INH (mg/(kg*day))
ab[0][3] = 0.06             # mu, uptake rate of INH (per day)
ab[1][0] = 0                # alphaBar, mutation rate of pyrazinamide (PZA) antibiotic (mutations*generation)
ab[1][1] = 0.0001625        # dBar, elimination rate of sensitive bacteria due to PZA (per day)
ab[1][2] = 35               # delta, daily dose of PZA (mg/(kg*day): range 35-20)
ab[1][3] = 0.03             # mu, uptake rate of PZA (per day)
ab[2][0] = 9.8*10**-9       # alphaBar, mutation rate of rifampicin (RIF) antibiotic (mutations*generation)
ab[2][1] = 0.00283          # dBar, elimination rate of sensitive bacteria due to RIF (per day)
ab[2][2] = 50               # delta, daily dose of RIF (mg/(kg*day))
ab[2][3] = 0.1              # mu, uptake rate of RIF (per day)
dt = 1                      # time step (day)
tf = 100                    # final time (days)

"""
Normalization steps:
s = S/T
r = R/T
b = B/(omega*T)
ai = Ai/(deltai/mui)
"""

p = np.zeros((int(tf/dt)+1, 8))  	# mock patient, rows = time, columns = time (t), sensitive bacteria (s), resistant bacteria (r), immune cells (b), INH (a1), PZA (a2), RIF (a3)
p[0][1] = nSens/T 					# sensitive bacteria (s)					
p[0][2] = nRes/T 					# resistant bacteria (r)
p[0][3] = nImm/(omega*T)			# immune cells (b)
p[0][4] = cAB1/(ab[0][2]/ab[0][3])	# INH
p[0][5] = cAB2/(ab[1][2]/ab[1][3])  # PZA
p[0][6] = cAB3/(ab[2][2]/ab[2][3])  # RIF
p[0][7] = p[0][1] + p[0][2]			# Total cell population

for i in range(1, len(p)):	# Iterate through time
    p[i][0] = i*dt          # Save time
    deathSum, mutSum = 0,0
    for j in range(len(ab)):	# Iterate through the antibiotics
        alpha = ab[j][0]*(ab[j][2]/ab[j][3]) # Calculte alpha
        d = ab[j][1]*(ab[j][2]/ab[j][3]) # Calcualte d
        deathSum += (alpha + d)*p[i-1][j+4] # Find sum for ds/dt
        mutSum += alpha*p[i-1][j+4]	# Find sum for dr/dt
    # Calculate change in sensitive and resistant populations over time
    p[i][1] = p[i-1][1] + (betaS*p[i-1][1]*(1-(p[i-1][1]+p[i-1][2])) - eta*p[i-1][1]*p[i-1][3] - p[i-1][1]*deathSum)*dt
    p[i][2] = p[i-1][2] + (betaR*p[i-1][2]*(1-(p[i-1][1]+p[i-1][2])) - eta*p[i-1][2]*p[i-1][3] + p[i-1][1]*mutSum)*dt
    # Calculate change in immune cell population over time
    p[i][3] = p[i-1][3] + k*p[i-1][3]*(1-(p[i-1][3]/(p[i-1][1]+p[i-1][2])))*dt
    # Calculate change in antibiotic concentration over time 
    p[i][4] = p[i-1][4] + ab[0][3]*(1-p[i-1][4])*dt
    p[i][5] = p[i-1][5] + ab[1][3]*(1-p[i-1][5])*dt
    p[i][6] = p[i-1][6] + ab[2][3]*(1-p[i-1][6])*dt
    # Calculate total change in bacteria population over time
    p[i][7] = p[i][1] + p[i][2]

    # Tunable Parameters:
    # if p[i-1][0] < 75: # Change day for when T-Cell therapy is applied
    #     p[i][3] = p[i-1][3] + k*p[i-1][3]*(1-(p[i-1][3]/(p[i-1][1]+p[i-1][2])))*dt
    # else: # Apply on T-Cell therapy
    #     eta = 0.9
    #     p[i][3] = p[i-1][3] + k*p[i-1][3]*(1-(p[i-1][3]/(z*(p[i-1][1]+p[i-1][2]))))*dt
    # # Tune when antibiotics are applied. Above = all at once, this staggers it
    # if p[i-1][0] > 15 and p[i-1][0] < 100:
    #     p[i][4] = p[i-1][4] + ab[0][3]*(1-p[i-1][4])*dt
    # if p[i-1][0] > 30:
    #     p[i][5] = p[i-1][5] + ab[1][3]*(1-p[i-1][5])*dt
    # if p[i-1][0] > 45:
    #     p[i][6] = p[i-1][6] + ab[2][3]*(1-p[i-1][6])*dt
    # p[i][7] = p[i][1] + p[i][2]


# np.savetxt("FinalProj.csv", p, fmt='%1.3f', header='Time(hr),s,r,b,INH,PZA', delimiter=',', comments='')

# data = np.loadtxt('FinalProj.csv', delimiter=',', skiprows=1)                   # Import data from csv
data = p
ax = plt.subplot()
line1 = ax.plot(data[:,0], data[:,1], c='green', label='sensitive bacteria')
line2 = ax.plot(data[:,0], data[:,2], c='red', label='resistant bacteria')
line3 = ax.plot(data[:,0], data[:,3], c='blue', label='immune cells')
line4 = ax.plot(data[:,0], data[:,4], c='purple', label='INH')
line5 = ax.plot(data[:,0], data[:,5], c='black', label='PZA')
line6 = ax.plot(data[:,0], data[:,6], c='gray', label='RIF')
# line7 = ax.plot(data[:,0], data[:,7], c='gray', label='total bacteria population')
plt.xlabel('Time (Days)')
plt.ylabel('Normalized Concentration')
plt.title('Patient Model')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.show()
