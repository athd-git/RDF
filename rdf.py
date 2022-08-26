# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def distance(a, b): #Distance between 2 arrays (x,y,z) coordinates of particles
    dx = abs(a[0] - b[0])
    x = min(dx, abs(A - dx)) #Computational parameters A,B,C
     
    dy = abs(a[1] - b[1])
    y = min(dy, abs(B - dy))
     
    dz = abs(a[2] - b[2])
    z = min(dz, abs(C - dz))
 
    return sp.sqrt(x**2 + y**2 + z**2)

def volume(z, r, z_bot, z_top):
    """ volume of a sphere of radius r located at height z """
    volume = 4.0 / 3.0 * sp.pi * r**3
     
    """ cut off a spherical cap if the sphere extends below z_bot """
    if z - r < z_bot:
        h = z_bot - (z - r)
        volume -= sp.pi * h**2 * (r - h / 3.0)
 
    """ cut off a spherical cap if the sphere extends above z_top """
    if z + r > z_top:
        h = z + r - z_top
        volume -= sp.pi * h**2 * (r - h / 3.0)
     
    return volume

class Trajectory: # For passing .xyz files
    def __init__(self, filename, skip, z_bot_interface, z_top_interface,
                 interface_offset=0.4, resolution=32):
 
        self.filename = filename
        self.skip = skip
        self.z_top = z_top_interface - interface_offset
        self.z_bot = z_bot_interface + interface_offset
        self.resolution = resolution
         
        self.parse_input()
        self.compute_volume_per_h2o()
        

    def parse_input(self):
        with open(self.filename, 'r') as f:
            data = f.readlines()
 
        self.n_atoms = int(data[0].split()[0])
        self.n_steps_total = int(len(data) / (self.n_atoms + 2))
        self.n_steps = self.n_steps_total // self.skip
         
        self.atom_list = [line.split()[0] for line in
                          data[2 : self.n_atoms + 2]]
 
        self.coordinates = np.zeros((self.n_steps, self.n_atoms, 3))
        for step in range(self.n_steps):
            coords = np.zeros((self.n_atoms, 3))
             
            i = step * self.skip * (self.n_atoms + 2)
            for j, line in enumerate(data[i + 2 : i + self.n_atoms + 2]):
                coords[j, :] = [float(value) for value in line.split()[1:]]
             
            self.coordinates[step] = coords
            
    
    def compute_volume_per_h2o(self): #Volume of computational cells
        n_h2o = 0
        for step in range(self.n_steps):
            for i, atom in enumerate(self.coordinates[step]):
                z = atom[2]
                if self.atom_list[i] == "O" and self.z_bot < z < self.z_top:
                    n_h2o += 1
         
        volume = A * B * (self.z_top - self.z_bot)
        average_n_h2o = n_h2o / self.n_steps
        self.volume_per_h2o = volume / average_n_h2o
        
        
    def compute_radial_distribution(self):
        """ Mirror images start to be double-counted so don't have 
        to go above half of the smallest lattice parameter """
        r_cutoff = min(A, B) / 2.0
 
        dr = r_cutoff / self.resolution
        self.radii = sp.linspace(0.0, self.resolution * dr, self.resolution)
 
        volumes = sp.zeros(self.resolution)
        self.g_of_r = sp.zeros(self.resolution)
         
        for step in range(self.n_steps):
            print('{:4d} : {:4d}'.format(step, self.n_steps))
             
            """ isolate all liquid water molecules based on their oxygens """
            data_oxygen = []
            for i, atom in enumerate(self.coordinates[step]):
                if self.atom_list[i] == "O":
                    if self.z_bot < atom[2] < self.z_top:
                        data_oxygen.append(atom)
            data_oxygen = sp.array(data_oxygen)
             
            for i, oxygen1 in enumerate(data_oxygen):
                """ calculate the volume of the cut spherical shells """
                for j in range(self.resolution):
                    r1 = j * dr
                    r2 = r1 + dr
                    v1 = volume(oxygen1[2], r1, self.z_bot, self.z_top)
                    v2 = volume(oxygen1[2], r2, self.z_bot, self.z_top)
                    volumes[j] += v2 - v1
 
                """ loop over pairs of O atoms, each pair counts as 2 """
                for oxygen2 in data_oxygen[i:]:
                    dist = distance(oxygen1, oxygen2)
                    index = int(dist / dr)
                    if 0 < index < self.resolution:
                        self.g_of_r[index] += 2.0
         
        for i, value in enumerate(self.g_of_r):
            self.g_of_r[i] = value * self.volume_per_h2o / volumes[i]
            
    
    def plot(self, filename=""):
        """ plots the radial distribution function
        if filename is specified, prints it to file as a png """
         
        if not self.g_of_r.any():
            print('compute the radial distribution function first\n')
            return
         
        plt.xlabel('r (Å)')
        plt.ylabel('g$_{OO}$(r)')
        plt.plot(self.radii, self.g_of_r)
         
        if filename:
            plt.savefig('radial_distribution_function.png', dpi=300,
                        bbox='tight', format='png')
         
        plt.show()
        

A = 18.991
B = 16.832
C = 30.577
bottom_interface = 8
top_interface = 10.0
 
H2O = Trajectory('water_traj.xyz', 10, bottom_interface, top_interface)
H2O.compute_radial_distribution()
H2O.plot()

