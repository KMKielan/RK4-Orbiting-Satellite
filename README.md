# RK4-Orbiting-Satellite

Program that simulates the orbit of an orbiting mass around a central mass(satellite and earth) using the Runge-Kutta 4th Order numerical method. The initial conditions of the satellite are those find for a geostationary orbit. With an initial position above the Earth's surface of 35786 kilometres and an initial velocity of 3.07 kilometres per second resulting in a circular orbit. Modifying the initial speed slightly changes the eccentricity of the orbit, making it more elliptical the more the value is changed from 3.07km/s. The period of the simulation can be changed, however, currently is set for one day, as a geostationary takes around that much time for single orbit. The entire animation is written using the tkinter module and optional things like colour, canvas dimensions and even tracer and radial lines(tracer and radial lines being graphical features I have added) can be changed or disabled.

Below can be seen two images of the satellite's orbit with differing initial velocities:

<img src="https://github.com/KMKielan/RK4-Orbiting-Satellite/blob/master/CircularOrbit.gif" width="250" height="250" /> <img src="https://github.com/KMKielan/RK4-Orbiting-Satellite/blob/master/EllipticalOrbit.gif" width="250" height="250" />

The image on the left with an initial velocity of 3.07km/s and the image on the right with an initial velocity of 2.07km/s. Changing the velocity by 1km/s has drastically altered the trajectory of the satellite around the Earth. Changing the velocity even more can result in a hyperbolic orbit with the satellite never returning.
