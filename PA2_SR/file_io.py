"""
Modules to extract fiducial data

Sidharth Raghavan
"""

import numpy as np
from typing import Tuple


def read_ct_fiducials(filename: str) -> Tuple[int, str, np.ndarray]:
    """
    Read CT fiducials data from file, return number of CT points, and coordinates of fiducials 
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    #header parising
    header_parts = [x.strip() for x in lines[0].strip().split(',')]
    Nb = int(header_parts[0])
    name = header_parts[1]
    
    #iterate through all data points
    data_points = []
    for line in lines[1:]:
        if line.strip():
            coords = [float(x.strip()) for x in line.split(',')]
            data_points.append(coords)
    
    #get first Nb rows
    bs = np.array(data_points[:Nb])
    
    return Nb, name, bs