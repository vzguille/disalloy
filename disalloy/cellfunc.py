import numpy as np
import ase


def rnd_ball(dimension=3):
    # Generate a vector with normally distributed components
    v = np.random.normal(0, 1, dimension)
    # Normalize the vector to make it a unit vector
    v /= np.linalg.norm(v)
    # Scale factor to ensure uniform distribution in the sphere
    r = np.random.uniform(0, 1) ** (1.0 / dimension)
    # Scale the unit vector to a random point within the unit sphere
    return v * r 
    

def jitter_ia(structure, ia=0.01, seed=None):
    """
    jitter atoms, based on cellcvrt.c++ by atat

    Parameters:
    structure: ase atom structure.
    ia: maximum amount of jitter in Angstroms.
    """
    if seed is not None:
        np.random.seed(seed)
    translation = ia*np.array([rnd_ball() for _ in range(len(structure))])
    structure.set_positions(structure.positions + translation)
    return structure


def jitter_ic(structure, ic=0.005, seed=None):
    """
    jitter cell, based on cellcvrt.c++ by atat

    Parameters:
    structure: ase atom structure.
    ic: maximum amount of jitter in Angstroms.
    """
    if seed is not None:
        np.random.seed(seed)
    transform = np.identity(3) + ic*np.array([rnd_ball() for _ in range(3)])
    structure.set_cell(transform@structure.get_cell(), scale_atoms=True)
    return structure


def compute_distortion(before, after, quiet=True):
    if type(before) is ase.atoms.Atoms:
        before = before.cell.array
    if type(after) is ase.atoms.Atoms:
        after = after.cell.array
        
    # Normalize the matrices by their cubic root of the determinant
    before /= np.power(np.linalg.det(before), 1.0 / 3.0)
    after /= np.power(np.linalg.det(after), 1.0 / 3.0)
    
    # Compute the difference matrix
    # diff = np.linalg.inv(before).dot(after)
    diff = np.linalg.inv(before) @ after
    
    # Symmetrize the difference matrix and subtract the identity matrix
    diff = (diff + diff.T) / 2.0 - np.eye(3)
    
    # Output the distortion measure
    if not quiet:
        print("Distorsion measure:")
        print(f"{np.linalg.norm(diff):.4f}")
    return np.linalg.norm(diff)