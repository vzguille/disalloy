import random
import time
from pymatgen.io import ase as pg_ase
from datetime import datetime


def write_log(message, log_file='application.log'):
    # Get the current time and format it
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Create the full log message
    log_entry = f"{timestamp} - {message}\n"
    
    # Open the log file and append the log entry
    print(log_entry)
    with open(log_file, 'a') as file:
        file.write(log_entry)





ase_to_pmg = pg_ase.AseAtomsAdaptor.get_structure
pmg_to_ase = pg_ase.AseAtomsAdaptor.get_atoms

def run_relaxer(structure, relaxer):
    b = relaxer(structure.copy())
    atoms_clean = b['trajectory'].atoms.copy()
    atoms_clean.calculator = None
    return atoms_clean, b['trajectory'].energies[-1]
    # return b['trajectory'].atoms.cell.array,
    # b['trajectory'].atoms.get_positions(), 
    # b['trajectory'].energies[-1]

def relax_row(row, calculator_label, relaxer):
    return run_relaxer(row['init_structure'], relaxer)


def run_packet(df_name, df_global, relaxer_dict, size_of_packet, structure_size, copy_calculated = None):
    write_log('starting packet run: attempting {} runs of size {}'.format(
        str(size_of_packet), str(structure_size)),
              log_file=df_name+'.log')
    calculator_label, relaxer = relaxer_dict['calculator_label'], relaxer_dict['relaxer']
    df = df_global #.copy()
    if not calculator_label in df.columns:
        df[calculator_label] = None
        df[calculator_label] = df[calculator_label].astype(object)
        df.attrs[calculator_label] = {}
        df.attrs[calculator_label]['calculated'] = []
        uncalculated = df.index
    else:
        uncalculated = df.attrs[calculator_label]['uncalculated']
    

    
    uncalculated_size = list(
        set(df[df['size'] == structure_size].index.to_list()) & 
        set(uncalculated)) 

    if copy_calculated is not None:
        if copy_calculated not in df.columns:
            write_log('this relaxer has not been applied yet and cannot be copied therefore',
              log_file=df_name+'.log')
            return
        if copy_calculated == calculator_label:
            write_log('same relaxer selected for copy,'
                'please select another relaxer calculated index, or appply a new relaxer',
              log_file=df_name+'.log')
            return
        uncalculated_size = list( set(uncalculated_size) & 
                                 set(df.attrs[copy_calculated]['calculated']))
        
    rest_number = len(uncalculated_size)
    if rest_number == 0:
        write_log('no more structures to run of this size', log_file=df_name+'.log')
        return
    if rest_number > size_of_packet:
        to_calculate = random.sample(uncalculated_size, size_of_packet)
    else:
        to_calculate = uncalculated_size
    ti = time.time()
    df.loc[to_calculate, calculator_label] = df.loc[to_calculate].apply(
        relax_row,args = (calculator_label, relaxer), axis=1)
    ti = time.time() - ti

    df.attrs[calculator_label]['uncalculated'] = list(set(uncalculated) - set(to_calculate))
    df.attrs[calculator_label]['calculated'] += to_calculate

    df.to_pickle(df_name+'.pkl')
    write_log('index of ran strctures:\n {}'.format(str(to_calculate)),
              log_file=df_name+'.log')
    write_log('Succesfuly ran packet of {} runs in {}s'.format(len(to_calculate), ti),
              log_file=df_name+'.log')