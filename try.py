"""
from math import *
import random


landmarks = [[20.0, 20.0], [80.0, 80.0], [20.0, 80.0], [80.0, 20.0]]
world_size = 100.0

class robot:
    def __init__(self):
        self.x = random.random() * world_size
        self.y = random.random() * world_size
        self.orientation = random.random() * 2.0 * pi
        self.forward_noise = 0.0;
        self.turn_noise = 0.0;
        self.sense_noise = 0.0;

    def set(self, new_x, new_y, new_orientation):
        if new_x < 0 or new_x >= world_size:
            raise ValueError, 'X coordinate out of bound'
        if new_y < 0 or new_y >= world_size:
            raise ValueError, 'Y coordinate out of bound'
        if new_orientation < 0 or new_orientation >= 2 * pi:
            raise ValueError, 'Orientation must be in [0..2pi]'
        self.x = float(new_x)
        self.y = float(new_y)
        self.orientation = float(new_orientation)

    def set_noise(self, new_f_noise, new_t_noise, new_s_noise):
        # makes it possible to change the noise parameters
        # this is often useful in particle filters
        self.forward_noise = float(new_f_noise);
        self.turn_noise = float(new_t_noise);
        self.sense_noise = float(new_s_noise);

    def sense(self):
        Z = []
        for i in range(len(landmarks)):
            dist = sqrt((self.x - landmarks[i][0]) ** 2 + (self.y - landmarks[i][1]) ** 2)
            dist += random.gauss(0.0, self.sense_noise)
            Z.append(dist)
        return Z
            
    def move(self, turn, forward):
        if forward < 0:
            raise ValueError, 'Robot cant move backwards'
        # turn, and add randomness to the turning command
        orientation = self.orientation + float(turn) + random.gauss(0.0, self.turn_noise)
        orientation %= 2 * pi

        # move, and add randomness to the motion command
        dist = float(forward) + random.gauss(0.0, self.forward_noise)
        x = self.x + (cos(orientation) * dist)
        y = self.y + (sin(orientation) * dist)
        x %= world_size # cyclic truncate
        y %= world_size
        # set particle
        res = robot()
        res.set(x, y, orientation)
        res.set_noise(self.forward_noise, self.turn_noise, self.sense_noise)
        return res
            
    def Gaussian(self, mu, sigma, x):
        # calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma
        return exp(- ((mu - x) ** 2) / (sigma ** 2) / 2.0) / sqrt(2.0 * pi * (sigma ** 2))

    def measurement_prob(self, measurement):
        # calculates how likely a measurement should be
        prob = 1.0;
                
        for i in range(len(landmarks)):
            dist = sqrt((self.x - landmarks[i][0]) ** 2 + (self.y - landmarks[i][1]) ** 2)
            prob *= self.Gaussian(dist, self.sense_noise, measurement[i])
        return prob
 
    def __repr__(self):
        return '[x=%.6s y=%.6s orient=%.6s]' % (str(self.x), str(self.y), str(self.orientation))
            
def eval(r, p):
    sum = 0.0;
    for i in range(len(p)): # calculate mean error
        dx = (p[i].x - r.x + (world_size/2.0)) % world_size - (world_size/2.0)
        dy = (p[i].y - r.y + (world_size/2.0)) % world_size - (world_size/2.0)
        err = sqrt(dx * dx + dy * dy)
        sum += err
    return sum / float(len(p))
        
# --------

N = 1000
T = 10
myrobot = robot()

p = []
for i in range(N):
    r = robot()
    r.set_noise(0.05, 0.05, 5.0) # provided noise.
    p.append(r)
    
print 'Mean error at start', eval(myrobot, p)
# show particle's initial locations
print p    

for t in range(T):
#    print p
    myrobot= myrobot.move(0.1, 5.0)
    Z = myrobot.sense()
    
    p2 = []

    for i in range(N):
        p2.append(p[i].move(0.1, 5.0))
    p = p2
    
    w = []
    for i in range(N):
        w.append(p[i].measurement_prob(Z))
    p3 = []

# this is importance sampling code

    index = int(random.random() * N)
    beta = 0.0
    mw = max(w)
    for i in range(N):
        beta += random.random() * 2.0 * mw
        while beta > w[index]:
            beta -= w[index]
            index = (index + 1) % N
        p3.append(p[index])
    p = p3
    
    print 'Mean error',eval(myrobot, p)
    
print ' '
if eval(myrobot, p) > 0.0:
    for i in range(N/100):        
        print 'Final particle #', i*100, p[i*100]
    print ' '
    print 'Actual Robot Location', myrobot
"""






""" cal num of snps for each chr
import numpy as np


data = np.load("./list_snp_info.npy")
print len(data)
rep = {}
for i in range(1, 23):
	chr = str(i)
	rep[chr] = 0
for i in range(len(data)):
	chr = data[i][1]
	rep[chr] += 1

list_pos = []
for i in range(1, 23):
	chr = str(i)
	list_pos.append(rep[chr])
list_pos = np.array(list_pos)
np.save("./list_num_snp", list_pos)

for i in range(len(list_pos)):
	print i, list_pos[i]
print np.sum(list_pos)
"""





import numpy as np


list_gene = np.load("./data_prepared/Gene_list.npy")
print "there are # of genes:", len(list_gene)


##==== gene tss repo
rep_gene_tss = {}
file = open("../../GTEx_data/data_raw/gene_tss_gencode.v19.v6p.txt", 'r')
while 1:
	line = (file.readline()).strip()
	if not line:
		break

	line = line.split('\t')
	if line[1] == 'X' or line[1] == 'Y' or line[1] == 'MT':
		continue

	gene = line[0]
	chr = int(line[1])
	tss = int(line[2])
	rep_gene_tss[gene] = (chr, tss)
file.close()



rep = {}
for i in range(1, 23):
	chr = str(i)
	rep[chr] = 0
for gene in list_gene:
	chr = str(rep_gene_tss[gene][0])
	rep[chr] += 1

list_pos = []
for i in range(1, 23):
	chr = str(i)
	list_pos.append(rep[chr])
list_pos = np.array(list_pos)
np.save("./data_prepared/list_num_gene", list_pos)

for i in range(len(list_pos)):
	print i, list_pos[i]
print np.sum(list_pos)



