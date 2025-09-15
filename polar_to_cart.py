import numpy as np

"""function that converts polar to cartesian coordinates based on 
   arrays containing necessary polar coordinates"""
def polar_to_cartesian(polar_coord):
    """input: polar_coord (array), array containing values for 
              radial position, theta angular position, radial velocity 
              and angular velocity
              
       returns: cartesian_vals (array), array containing cartesian 
                values position x, position y, velocity x and velocity 
                y in same order as they are given in 
    """
    
    polar_coord = np.atleast_2d(np.asarray(polar_coord, dtype=float)) #shaping for iterations
    
    cartesian_vals = []  #list for adding cartesian values
    
    for r , theta , vr , vtheta in polar_coord:
        
        x = r * np.cos(theta) #polar to x position conversion
        y = r * np.sin(theta) #polar to y position conversion
    
        vx = vr * np.cos(theta) - vtheta * np.sin(theta) #conv polar to x vel
        vy = vr * np.sin(theta) + vtheta * np.cos(theta) #conv polar to y vel
    
        cartesian_vals.append([x , y , vx , vy]) #adding cartesian coords to list
    cartesian_vals = np.array(cartesian_vals) #list into array
    return cartesian_vals
    