import spiceypy as spice
import numpy as np

def get_objects(fileName, display=False):
    objects = spice.spkobj(fileName)
    ids, names, tcs_sec, tcs_cal = [],[],[],[]
    n = 0
    if display:
        print(f'\nObjects in {fileName}')

    for o in objects:
        ids.append(o)

        # time coverage in seconds since J2000
        tc_sec = spice.wnfetd(spice.spkcov(fileName, ids[n]),n)

        # Human readable time format
        tc_cal = [spice.timout(f, "YYYY MON DD HR:MN:SC.### (TDB) :: (TBD)") for f in tc_sec]

        # Append time coverages to outputs
        tcs_sec.append(tc_sec)
        tcs_cal.append(tc_cal)

        # get Body name
        try:
            names.append(id2body(o))
        except:
            names.append('Unknown Name')

        if display:
            print(f'id:{ids[-1]}\t\tname: {names[-1]}\t\t tc: {tcs_sec[-1]} --> {tcs_cal[-1]}')

        return ids, names, tcs_sec, tcs_cal
