from mango.mango_features.preserved_regions import PreservedEdge, PreservedVertex
import numpy as np

if __name__ == '__main__':
    c1, c2 = np.array([0, 0, 0]), np.array([10, 10, 10])
    preserved_vertex = PreservedVertex(v1=c1)
    preserved_edge = PreservedEdge(v1=c1, v2=c2)
    