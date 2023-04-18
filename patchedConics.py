'''
A script to calculated the necessary Hohmann transfer parameters
for a interplanetary injection using patched conics
'''

import numpy as np

R1 = 149.6e6
R2 = 1.426e9

mu_sun = 132.71e9
mu_e   = 398600
mu_sat = 126686000
r_sat = 60270

# Parking orbit at departure (circular)
rpD = 185.2 + 6371

# Parking orbit at arrival 
e = 0.951

# --------- Departure ------------
v_infD = np.sqrt(mu_sun/R1) *(np.sqrt(2*R2/(R1+R2))-1)

v_circ = np.sqrt(mu_e/rpD)

deltaVd = v_circ*(np.sqrt(2+(v_infD/v_circ)**2)-1)
betaD = np.arccos(1/(1+(rpD*v_infD**2)/mu_e))

print('---------- Departure -----------')
print(f'V_inf = {v_infD:.2f}')
print(f'V_circ = {v_circ:.2f}')
print(f'deltaV = {deltaVd:.2f}')
print(f'beta = {betaD*180/np.pi:.2f}')

# --------- Arrival ------------
v_infA = np.sqrt(mu_sun/R2) *(1-np.sqrt(2*R1/(R1+R2)))

deltaVa = v_infA*np.sqrt((1-e)/2)

rpA = 2*mu_sat/v_infA**2 * (1-e)/(1+e)
raA = rpA*(1+e)/(1-e)

betaA = np.arccos(1/(1+(rpA*v_infA**2)/mu_sat))
DeltaA = rpA*np.sqrt(2/(1-e))

print('---------- Arrival -----------')
print(f'V_inf = {v_infA:.2f}')
print(f'deltaV = {deltaVa:.2f}')
print(f'beta = {betaA*180/np.pi:.2f}')
print(f'Aiming radius (\Delta) = {DeltaA:.2f}')
print(f'[{rpA-r_sat:.0f} km by {raA-r_sat:.0f} km] about Saturn (e = {e:.2f})')
