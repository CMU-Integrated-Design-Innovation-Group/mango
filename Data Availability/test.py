import dill
import os

def edge_cutoff_constraint():
    pass

def porosity_objective():
    pass

def cylinder_volume():
    pass


times = []
for des in os.listdir('./data_output/InternalFeatures'):
    if '.DS_Store' not in des:
        fp = os.path.join('./data_output/InternalFeatures', des)
        with open(fp, 'rb') as f:
            data = dill.load(f)
        print(data.sim_time_seconds)