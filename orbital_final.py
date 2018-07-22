#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 21:48:44 2018
@author: Konrad
"""
import numpy as np
import tkinter as tk

class Orbital:
    def __init__(self, orbiting_mass, central_mass, initial_position, initial_velocity, time, step):
        self.o_mass = orbiting_mass # defines mass of object orbiting central mass
        self.c_mass = central_mass # defines mass of central object in 2-body system
        self.initial_pos = initial_position # defines starting postion of orbiting mass relative to central mass
        self.initial_vel = initial_velocity # defines starting velocity of orbiting mass
        self.step = step # defines time step 
        self.Grav_constant = 6.67408e-11 # Gravitational Constant
        self.end = time # defines when simulation will end
    def time(self, start = 0):
        """function that returns an array of times between the initial time and final time with selected time step."""
        return np.linspace(start, self.end, int(self.end / self.step) + 1)
    def xvelocity(self, pos_x, pos_y):
        """function that returns the calculated x-component of the orbiting mass's velocity affected by the central mass,
        using Newtons law of gravitation."""
        return -(self.Grav_constant * self.c_mass * pos_x) / ((pos_x ** 2 + pos_y ** 2) ** (3 / 2))
    def yvelocity(self, pos_x, pos_y):
        """function that returns the calculated y-component of the orbiting mass's velocity affected by the central mass,
        using Newtons law of gravitation."""
        return -(self.Grav_constant * self.c_mass * pos_y) / ((pos_x ** 2 + pos_y ** 2) ** (3 / 2))
    def k_pos_vel(self, pos_x, pos_y, vel_x, vel_y):
        """function that returns multiple arrays of Runge-Kutta 4th Order values for position and velocity of the orbiting mass
        in x and y components. These arrays are later used with the tstepping functions for velocity and position to determine
        how the orbiting mass moves around the central mass."""
        # calculates k1 for x and y component of position and velocity
        k1pos_x = vel_x
        k1pos_y = vel_y
        k1vel_x = self.xvelocity(pos_x, pos_y)
        k1vel_y = self.yvelocity(pos_x, pos_y)
        # calculates k2 for x and y component of position and velocity
        k2pos_x = vel_x + (self.step * k1vel_x) / 2
        k2pos_y = vel_y + (self.step * k1vel_y) / 2
        k2vel_x = self.xvelocity(pos_x + (self.step * k1pos_x) / 2, pos_y + (self.step * k1pos_y) / 2)
        k2vel_y = self.yvelocity(pos_x + (self.step * k1pos_x) / 2, pos_y + (self.step * k1pos_y) / 2)
        # calculates k3 for x and y component of position and velocity
        k3pos_x = vel_x + (self.step * k2vel_x) / 2
        k3pos_y = vel_y + (self.step * k2vel_y) / 2
        k3vel_x = self.xvelocity(pos_x + (self.step * k2pos_x) / 2, pos_y + (self.step * k2pos_y) / 2)
        k3vel_y = self.yvelocity(pos_x + (self.step * k2pos_x) / 2, pos_y + (self.step * k2pos_y) / 2) 
        # calculates k4 for x and y component of position and velocity
        k4pos_x = vel_x + self.step * k3vel_x
        k4pos_y = vel_y + self.step * k3vel_y
        k4vel_x = self.xvelocity(pos_x + self.step * k3pos_x, pos_y + self.step * k3pos_y)
        k4vel_y = self.yvelocity(pos_x + self.step * k3pos_x, pos_y + self.step * k3pos_y)
        return [k1pos_x, k2pos_x, k3pos_x, k4pos_x], [k1pos_y, k2pos_y, k3pos_y, k4pos_y], [k1vel_x, k2vel_x, k3vel_x, k4vel_x], [k1vel_y, k2vel_y, k3vel_y, k4vel_y]
    def tstepping_position(self, pos_x, pos_y, kpos_xlist, kpos_ylist):
        """function that returns the x and y component of the orbiting mass's position for the next time step."""
        pos_x = pos_x + (self.step / 6) * (kpos_xlist[0] + 2 * kpos_xlist[1] + 2 * kpos_xlist[2] + kpos_xlist[3])
        pos_y = pos_y + (self.step / 6) * (kpos_ylist[0] + 2 * kpos_ylist[1] + 2 * kpos_ylist[2] + kpos_ylist[3])
        return pos_x, pos_y
    def tstepping_velocity(self, vel_x, vel_y, kvel_xlist, kvel_ylist):
        """function that returns the x and y component of the orbiting mass's velocity for the next time step."""
        vel_x = vel_x + (self.step / 6) * (kvel_xlist[0] + 2 * kvel_xlist[1] + 2 * kvel_xlist[2] + kvel_xlist[3])
        vel_y = vel_y + (self.step / 6) * (kvel_ylist[0] + 2 * kvel_ylist[1] + 2 * kvel_ylist[2] + kvel_ylist[3])
        return vel_x, vel_y
    def simulate_orbit(self):
        """function that returns the x and y components of the orbiting mass's position and velocity throughout the entire simulation
        as arrays used to animate the simulation or used in matplotlib.pyplot to plot a graph of the mass's orbit."""
        pos_x, pos_y = self.initial_pos[0], self.initial_pos[1]
        vel_x, vel_y = self.initial_vel[0], self.initial_vel[1]
        time_list = self.time()
        rx, ry = np.zeros(len(time_list)), np.zeros(len(time_list))
        rdotx, rdoty = np.zeros(len(time_list)), np.zeros(len(time_list))
        for i in range(0, len(time_list)):
            kpos_xlist, kpos_ylist, kvel_xlist, kvel_ylist = self.k_pos_vel(pos_x, pos_y, vel_x, vel_y)
            pos_x, pos_y = self.tstepping_position(pos_x, pos_y, kpos_xlist, kpos_ylist)
            vel_x, vel_y = self.tstepping_velocity(vel_x, vel_y, kvel_xlist, kvel_ylist)
            rx[i], ry[i], rdotx[i], rdoty[i] = pos_x, pos_y, vel_x, vel_y
        return rx, ry, rdotx, rdoty
    def animation(rx, ry, rdotx, rdoty):
        """function that deals with the animation, initialising the tkinter canvas and orbiting objects to help visualise the animation."""
        root = tk.Tk()
        root.title("Orbitng Satellite Simulation")
        window_width = 500 # defines the window dimensions
        window_height = 500
        center_x = window_width / 2 #defines the centre of the window
        center_y = window_height / 2
        scale_x = (window_width - 100) / (2 * rx.max() - rx.min())  # values from the numerical method are too large for the canvas
        scale_y = (window_height - 100) / (2 * ry.max() - ry.min()) # and so are scaled so that the animation can be viewed properly
        rx *= scale_x # the following lines apply the scale to the simulated results
        ry *= -scale_y
        rx += center_x
        ry += center_y
        window = tk.Canvas(root, width = window_width, height = window_height, bg = "#02198c") # initialises the canvas with desired specifications
        window.pack()
        radius_e = 10 # radius for earth tkinter object
        radius_s = 5 # radius for satellite tkinter object
        x0e, y0e, x1e, y1e = center_x - radius_e, center_y - radius_e, center_x + radius_e, center_y + radius_e # next two lines position the objects
        x0s, y0s, x1s, y1s = center_x, center_y, center_x, center_y                                             # in the canvas so they are visible
        earth = window.create_oval(x0e, y0e, x1e, y1e, fill = "white", outline = "white")     # initalises the earth and moon object with desired
        satellite = window.create_oval(x0s, y0s, x1s, y1s, fill = "white", outline = "white") # positions and colour scheme
        radial_line = window.create_line(center_x, center_y, center_x, center_y, fill = "white") # creates radial line from central to orbital mass
        stats = window.create_rectangle(10, window_height - 5, 200, window_height - 45, fill = "white", outline = "") # creates a rectangular box used for showing stats
        earth_text = window.create_text(center_x + 10, center_y + 10, text = "EARTH", fill = "white", anchor = "nw")  # creates text that names the central mass on canvas
        sat_text = window.create_text(rx[0] + 10, ry[0] + 10, text = "SATELLITE", fill = "white", anchor = "s") # creates text that follows the orbiting mass
        polygon_coords = [] # creates and array that will hold orbiting mass position
        for i in range(0, len(rdotx)): # for loop appends the orbiting mass position to the polygon_coords array
            polygon_coords.append(rx[i])
            polygon_coords.append(ry[i])
        sat_path2 = window.create_polygon(polygon_coords, outline = "white", dash = "-.-", fill = "") # creates a polygon in the canvas that shows the trajectory of the orbiting mass
        counter = 0
        while counter <= len(rdotx) + 1:
            if counter == len(rdotx):
                counter = 0
            window.coords(satellite, rx[counter] - radius_s, ry[counter] - radius_s, rx[counter] + radius_s, ry[counter] + radius_s) # moves the orbiting mass accordingly
            window.coords(sat_text, rx[counter], ry[counter] - 10) # moves the sat_text to follow the orbiting mass text
            window.coords(radial_line, center_x, center_y, rx[counter], ry[counter]) # moves radial line from central to orbital mass
            distance = np.sqrt(((center_x - rx[counter]) ** 2) + ((center_y - ry[counter]) ** 2)) * 2 / (scale_x + scale_y) # determines the distance between the orbiting and central mass
            if distance < 6371e3:
                print("Unstable Orbit, Satellite has crashed!") # if distance between central and orbiting mass is smaller than the radius of the central mass, assume orbiting mass has crashed
                window.delete(stats)
                break # closes the tkinter window
            velocity = np.sqrt(rdotx[counter] ** 2 + rdoty[counter] ** 2) # next 7 lines of code display the stats of the orbiting mass in the stats box; distance from central mass and velocity
            velocity, distance = "%.4f" % velocity, "%.4f" % distance
            radial_dist = window.create_text(15, window_height - 25, text = "Radial Distance : " + distance + " m", fill = "#02198c", anchor = "sw")
            tang_velocity = window.create_text(15, window_height - 10, text = "Velocity : " + velocity + " m/s", fill = "#02198c", anchor = "sw")
            window.after(1, window.update())
            window.delete(radial_dist)
            window.delete(tang_velocity)
            counter += 1
        root.mainloop()
def main():
    """function that runs the program."""
    day = 24 * 60 * 60 # seconds
    step = 100 # seconds
    earth_radius = 6371e3 # metres
    orbital_radius = 35786e3 + earth_radius # metres
    earth_mass = 5.972e24 # kilograms
    satellite_mass = 768 # kilograms
    satellite_velocity = 3.07e3 # metres per second 
    multiplier = 1 # defines the amount of days the simulation should run through
    satellite = Orbital(satellite_mass, earth_mass, [0, orbital_radius], [satellite_velocity, 0], (multiplier * day), (10))
    rx, ry, rdotx, rdoty = satellite.simulate_orbit()
    Orbital.animation(rx, ry, rdotx, rdoty)
main() # calls the main function