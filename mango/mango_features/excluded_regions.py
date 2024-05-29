"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This modules defines excluded region areas that will be used in defining the polyhedral_design_space to prevent material
from being added into these spaces. These regions are specifically useful if a designer has a pre-allocated space to be
occupied by a nanoparticle, and this conceptually prevents DNA from being added to these locations.
"""
# Importing Modules:
from dataclasses import dataclass
import trimesh  # NOTE: When installing trimesh use pip install trimesh\[easy\]
from mango.utils.mango_math import *
from mango.mango_features.mesh_face import MeshFace

@dataclass
class RectangularPrism(object):
    """
    This is an excluded region representing a simple box. This is more "efficient" than the sphere as there are fewer
    geometrical faces to check, although it takes up more space that you may not want.

    Parameters
    ------------
    c1 : Corner 1 defined by a numpy array in [X Y Z] (nm). Note C1 < C2 for all dimensions.
    c2 : Corner 2 defined by a numpy array in [X Y Z] (nm)
    numDecimals: Number of decimals to round positional values to
    """
    # We locate the Nanoparticle based on the bottom left corner and the top right corners 3D Coordinates as a 1D array
    # These MUST be OPPOSING corners to define out the "bounding box" of sorts.
    c1: np.array
    c2: np.array
    numDecimals: int = 5

    def verify_corners_opposing(self) -> None:
        """
        This function verifies that the corners passed in ARE in fact opposing corners
        :return: Nothing, raises an exception if they are invalid.
        """
        # Basically, all values of C1 MUST be less than C2:
        if any(C1 >= C2 for C1, C2 in zip(self.c1, self.c2)):
            raise ValueError('The input conditions for the opposing corners must be that all values in c1 are less than'
                             ' the values found in c2 for all of X, Y, Z.')


    def point_inside_space(self, pt: np.array) -> bool:
        """
        This function will determine if a passed in point is located inside this space, and is therefore an invalid
        design
        :param pt: XYZ of a point we are checking
        :return: True: Point inside of the excluded region; False: point not inside
        """
        if (self.c1[0] < pt[0] < self.c2[0] or self.c2[0] < pt[0] < self.c1[0]) and \
                (self.c1[1] < pt[1] < self.c2[1] or self.c2[1] < pt[1] < self.c1[1]) and \
                (self.c1[2] < pt[2] < self.c2[2] or self.c2[2] < pt[2] < self.c1[2]):
            return True
        # If we are not within the box, we return False
        return False


    def edge_intersecting_face(self, p1: np.array, p2: np.array) -> bool:
        """
        This function will determine if a passed in line segment intersects a face of an excluded Region.
        :param p1: Point 1 of a line segment definition
        :param p2: Point 2 of a line segment definition
        :return: True: Edge intersects the excluded space; False: edge does not intersect
        """
        face_intersections = []
        u = p2 - p1
        intersector = trimesh.ray.ray_triangle.RayMeshIntersector(self.cube_mesh)
        intersections, ignore1, ignore2 = intersector.intersects_location([p1], [u], multiple_hits=False)

        # Now we check to verify the intersections:
        for hit in intersections:
            if np.array_equal(p1, hit) or np.array_equal(p2, hit):
                pass
            else:
                # In this case, we also need to check to see if this intersection hit is inbetween our two points or
                # not, since a Ray does not have any bounds.
                d_12 = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2 + (p2[2] - p1[2]) ** 2)
                d_1_hit = np.sqrt((hit[0] - p1[0]) ** 2 + (hit[1] - p1[1]) ** 2 + (hit[2] - p1[2]) ** 2)
                d_2_hit = np.sqrt((p2[0] - hit[0]) ** 2 + (p2[1] - hit[1]) ** 2 + (p2[2] - hit[2]) ** 2)
                # Now, if the following check is true, then we have a face intersection:
                if d_1_hit < d_12 and d_2_hit < d_12:
                    face_intersections.append(True)

        if any(face_intersections):
            return True
        else:
            return False


    def __post_init__(self):
        self.relabel_dict = {}
        self.faces = {}
        self.verify_corners_opposing()  # Verify corners are opposing
        self.cube_mesh = trimesh.creation.box(bounds=np.array([self.c1, self.c2]))
        self.excluded_region_as_graph = nx.Graph()  # Initialize a graph that will store data of these points
        # Calculate the differences in X Y and Z:
        dx = round(self.c2[0] - self.c1[0], self.numDecimals)
        dy = round(self.c2[1] - self.c1[1], self.numDecimals)
        dz = round(self.c2[2] - self.c1[2], self.numDecimals)

        # Manually define the graph nodes and edges:
        graph_points = [(0, {'x': self.c1[0], 'y': self.c1[1], 'z': self.c1[2], 'terminal': True}),
                        (1, {'x': self.c1[0] + dx, 'y': self.c1[1], 'z': self.c1[2], 'terminal': True}),
                        (2, {'x': self.c1[0], 'y': self.c1[1] + dy, 'z': self.c1[2], 'terminal': True}),
                        (3, {'x': self.c1[0], 'y': self.c1[1], 'z': self.c1[2] + dz, 'terminal': True}),
                        (4, {'x': self.c2[0] - dx, 'y': self.c2[1], 'z': self.c2[2], 'terminal': True}),
                        (5, {'x': self.c2[0], 'y': self.c2[1] - dy, 'z': self.c2[2], 'terminal': True}),
                        (6, {'x': self.c2[0], 'y': self.c2[1], 'z': self.c2[2] - dz, 'terminal': True}),
                        (7, {'x': self.c2[0], 'y': self.c2[1], 'z': self.c2[2], 'terminal': True})
                        ]
        graph_edges = [(0, 1), (0, 2), (0, 3), (3, 4), (3, 5), (4, 2), (1, 5), (1, 6), (6, 2), (6, 7), (7, 5), (7, 4)]

        self.excluded_region_as_graph.add_nodes_from(graph_points)
        self.excluded_region_as_graph.add_edges_from(graph_edges)

        all_faces = [(2, 6, 1, 0), (2, 0, 3, 4), (4, 3, 5, 7), (7, 5, 1, 6), (2, 4, 7, 6), (0, 1, 5, 3)]
        for f in all_faces:
            # We can't determine face-directionality until the 3D structure is finalized. However, there are certain
            # functionalities we will want for the face:
            verts = []
            for v in list(f):
                # Grab the X, Y, and Z values:
                pt = xyz_from_graph(graph=self.excluded_region_as_graph, node=v)
                verts.append(pt)
            self.faces[f] = MeshFace(v1=verts[0], v2=verts[1], v3=verts[2], v4=verts[3], face_type='rect')


@dataclass
class Sphere(object):
    """
    This is an excluded region representing a simple box. This is more "efficient" than the sphere as there are fewer
    geometrical faces to check, although it takes up more space that you may not want.

    Parameters
    ------------
    center : Center of sphere defined by a numpy array in [X Y Z] (nm).
    diameter : floating point of the sphere (nm)
    numDecimals: Number of decimals to round positional values to
    """
    center: np.array
    diameter: float
    numDecimals: int = 5

    def recenter_mesh(self, new_center: np.array) -> None:
        """
        This function simply translates the mesh to a defined center point as needed due to trimesh requirements
        :param new_center: Numpy array defining the center point of the spherical excluded region
        """
        # Now translate the mesh:
        translation = new_center - self.mesh.centroid
        self.mesh.apply_translation(translation)

        # After re-centering the mesh, we then must update the faces:
        self.update_icosa_faces()

    def point_inside_space(self, pt: np.array) -> bool:
        """
        This function will determine if a passed in point is located inside this space, and is therefore an invalid
        design
        :param pt: XYZ of a point we are checking
        :return: True: Point inside of the excluded region; False: point not inside
        """
        # If the point is contained within the trimesh we return True
        # We make this a 2D array because trimesh.contains wants it this way:d
        reshaped_arr = np.reshape(pt, (1, 3))
        if self.mesh.contains(reshaped_arr):
            return True
        # If we are not within the box, we return False
        return False

    def edge_intersecting_face(self, p1: np.array, p2: np.array) -> bool:
        """
        This function will determine if a passed in line segment intersects a face of an excluded Region.
        :param p1: Point 1 of a line segment definition
        :param p2: Point 2 of a line segment definition
        :return: True: Edge intersects the excluded space; False: edge does not intersect
        """
        face_intersections = []
        u = p2 - p1
        intersector = trimesh.ray.ray_triangle.RayMeshIntersector(self.mesh)
        intersections, ignore1, ignore2 = intersector.intersects_location([p1], [u], multiple_hits=False)

        # Now we check to verify the intersections:
        for hit in intersections:
            if np.array_equal(p1, hit) or np.array_equal(p2, hit):
                pass
            else:
                # In this case, we also need to check to see if this intersection hit is inbetween our two points or
                # not, since a Ray does not have any bounds.
                d_12 = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2 + (p2[2] - p1[2]) ** 2)
                d_1_hit = np.sqrt((hit[0] - p1[0]) ** 2 + (hit[1] - p1[1]) ** 2 + (hit[2] - p1[2]) ** 2)
                d_2_hit = np.sqrt((p2[0] - hit[0]) ** 2 + (p2[1] - hit[1]) ** 2 + (p2[2] - hit[2]) ** 2)
                # Now, if the following check is true, then we have a face intersection:
                if d_1_hit < d_12 and d_2_hit < d_12:
                    face_intersections.append(True)

        if any(face_intersections):
            return True
        else:
            return False

    def update_icosa_faces(self) -> None:
        """ This function simply creates a list of MeshFace elements for use in constraint checking """
        verts = self.mesh.vertices
        for f in self.mesh.faces:
            # We can't determine face-directionality until the 3D structure is finalized. However, there are certain
            # functionalities we will want for the face:
            self.faces[f] = MeshFace(v1=verts[f[0]], v2=verts[f[1]], v3=verts[f[2]])

    def __post_init__(self) -> None:
        self.relabel_dict = {}
        self.faces = {}
        self.mesh = trimesh.creation.icosphere(subdivisions=1, radius=self.diameter/2)
        # Now recenter the mesh initially:
        self.recenter_mesh(new_center=self.center)

