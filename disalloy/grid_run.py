
import random

def run_relaxer(structure, relaxer):
    b = relaxer(structure.copy())
    return b['trajectory'].atoms, b['trajectory'].energies[-1]

def relax_row(row, calculator_label, relaxer):
    return run_relaxer(row['init_structure'], relaxer)


def run_packet(df_global, relaxer_dict, size_of_packet, structure_size, copy_calculated = None):
    calculator_label = relaxer_dict['calculator_label']
    relaxer = relaxer_dict['relaxer']
    
    df = df_global#.copy()
    """ """
    
    if not calculator_label in df.columns:
        df[calculator_label] = None
        df[calculator_label] = df[calculator_label].astype(object)
        df.attrs[calculator_label] = {}
        df.attrs[calculator_label]['calculated'] = []
        uncalculated = df.index
        
    else:
        uncalculated = df.attrs[calculator_label]['uncalculated']
    #uncalculated_size = df[(df['size'] == structure_size) & uncalculated].index.to_list()

    
    uncalculated_size = list(
        set(df[df['size'] == structure_size].index.to_list()) & 
        set(uncalculated)) 

    if copy_calculated is not None:
        if copy_calculated not in df.columns:
            print('this relaxer has not been applied yet')
            return
        if copy_calculated == calculator_label:
            print('same relaxer selected for copy,'
                'please select another relaxer calculated index, or appply a new relaxer')
            return
        uncalculated_size = list( set(uncalculated_size) & 
                                 set(df.attrs[copy_calculated]['calculated']))
        
    rest_number = len(uncalculated_size)
    if rest_number == 0:
        print('no more structures to run of this size')
        return
    if rest_number > size_of_packet:
        to_calculate = random.sample(uncalculated_size, size_of_packet)
    else:
        to_calculate = uncalculated_size
        
    print(to_calculate)
    df.loc[to_calculate, calculator_label] = df.loc[to_calculate].apply(
        relax_row,args = (calculator_label, relaxer), axis=1)

    df.attrs[calculator_label]['uncalculated'] = list(set(uncalculated) - set(to_calculate))
    df.attrs[calculator_label]['calculated'] += to_calculate