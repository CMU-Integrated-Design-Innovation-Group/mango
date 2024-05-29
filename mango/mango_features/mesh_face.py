"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script contains a MeshFace which is simple a triangulated face. I offer the ability to use a square face but
I do not specifically use this functionality as it increases computational complexity slights (need more checks to
validate the face orientation).

This is generally not a user-facing class, it is really just used by internal functions.
"""

# Import Modules:
from dataclasses import dataclass
from mango.utils.mango_math import *
from typing import Optional
import plotly.graph_objs as go
from random import randint

@dataclass
class MeshFace(object):
    """
    The mesh face is a collection of vertices and edges. This class has some simple methods to plot / check various
    computational geometry attributes during the generative process.

    Parameters
    ------------
    v1 : Numpy array of X Y Z of a vertex defining a mesh face
    v2 : Numpy array of X Y Z of a second vertex defining a mesh face
    v3 : Numpy array of X Y Z of a third vertex defining a mesh face
    v4 : (OPTIONAL) Numpy array of X Y Z of a fourth vertex defining a mesh face
    face_type : String defining if the face is a triangle or a rectangle, not recommended to change
    numDecimals : integer that is automatically specified to 2. You can lower the precision if desired but consider DNA!
    """
    # Our faces can be either a triangle or a square, but mainly a triangle.
    v1: np.array
    v2: np.array
    v3: np.array
    v4: Optional[np.array] = None
    face_type: str = 'triangle'
    numDecimals: int = 2   # I use 2 decimals as that is the assumed bp rise in BDNA

    def __post_init__(self):
        # Validate typing...:
        self.v1 = np.array([self.v1[0], self.v1[1], self.v1[2]])
        self.v2 = np.array([self.v2[0], self.v2[1], self.v2[2]])
        self.v3 = np.array([self.v3[0], self.v3[1], self.v3[2]])
        self.area = 0.0

        if self.face_type != 'triangle' and self.face_type != 'rect':
            raise Exception(f'Invalid face type passed in: {self.face_type}')
        self.face_angles = [90, 90, 90, 90]
        if self.face_type == 'triangle':
            # Calculate the area:
            AB = self.v2 - self.v1
            AC = self.v3 - self.v1
            self.area = 0.5 * np.linalg.norm(np.cross(AB, AC))
        elif self.face_type == 'rect':
            # Calculating area of a rectangular face by adding two triangles:
            AB = self.v2 - self.v1
            AC = self.v3 - self.v1
            A1 = 0.5 * np.linalg.norm(np.cross(AB, AC))

            AD = self.v4 - self.v1
            A2 = 0.5 * np.linalg.norm(np.cross(AC, AD))
            self.area = A1 + A2


    def equation_of_face(self) -> tuple[float, float, float, float]:
        """
        This function will calculate the equation of the triangular face and store it as the plane values of a b c d

        :return: Values of a, b, c, and d to define a plane
        """
        AB = self.v2 - self.v1
        AC = self.v3 - self.v1
        a, b, c = np.cross(AB, AC)
        d = -1 * np.dot(np.cross(AB, AC), self.v1)
        return a, b, c, d


    def get_face_plot(self, color_type: str, override_opacity: bool = False, override_color: str = None) -> go.Mesh3d:
        """
        This function will be used to plot the mesh face simply via Plotly

        :param override_color: Used to color a face a specific color
        :param override_opacity: Determines opacity of face plot
        :param color_type: Dictates which color to use in the plot via the plotting dictionary
        :return: A Mesh3D value that will be used to plot
        """
        # Color dictionary used in plotting the faces below:
        color_key = {'NP': 'red', 'binding_region': 'green', 'other': 'grey', 'manual': override_color}

        if self.face_type == 'rect':
            all_verts = np.stack((self.v1, self.v2, self.v3, self.v4), axis=0)
            i = [0, 0]
            j = [1, 2]
            k = [2, 3]
        else:
            all_verts = np.stack((self.v1, self.v2, self.v3), axis=0)
            i = [0]
            j = [1]
            k = [2]
        # Now create the plot trace and return:
        if override_opacity:
            return go.Mesh3d(x=all_verts[:, 0], y=all_verts[:, 1], z=all_verts[:, 2], i=i, j=j, k=k,
                             color=color_key[color_type], flatshading=True, opacity=0.9)
        else:
            return go.Mesh3d(x=all_verts[:, 0], y=all_verts[:, 1], z=all_verts[:, 2], i=i, j=j, k=k,
                             color=color_key[color_type], flatshading=True, opacity=1)


    def divide_triangular_face(self):
        """ Divide face functionality for computational cost efficiency """
        if self.v4 is not None:
            raise Exception('A square / non-triangular face was selected which should never happen!')
        # We select a random vertex on this face and calculate the bisection to the other triangle:
        verts = [self.v1, self.v2, self.v3]
        randIndex = randint(0, 2)  # Select a random vertex
        A = verts.pop(randIndex)

        # Next calculate the midpoint of the two remaining vertices:
        bisecting_midpoint = (verts[0] + verts[1]) / 2

        # New triangles:
        new_triangle1 = (A, bisecting_midpoint, verts[0])
        new_triangle2 = (A, bisecting_midpoint, verts[1])

        # We will return this mid-point, the two new triangles we will be forming, as well as the divided_edge which is
        # just the remaining vertices in verts
        return bisecting_midpoint, new_triangle1, new_triangle2, verts, A


    def update_vertex_values(self, old_xyz: np.array, new_xyz: np.array) -> None:
        """
        This function is called whenever the vertices of a face change values. This function will also update teh abcd
        values.
        :param old_xyz: Old position to search for to change the correct vertex XYZ
        :param new_xyz: New position to update with
        :return: Nothing, just changing values of object
        """
        # We first ensure decimal place values:
        rounded_new_xyz = np.round(new_xyz, decimals=self.numDecimals)

        # I do not know if there is a better way to do this, but since the faces only have a few values I am just
        # going to use if statements. I tried a for loop but that didn't update the object value:
        if np.array_equal(old_xyz, self.v1):
            self.v1 = rounded_new_xyz
        elif np.array_equal(old_xyz, self.v2):
            self.v2 = rounded_new_xyz
        elif np.array_equal(old_xyz, self.v3):
            self.v3 = rounded_new_xyz
        elif np.array_equal(old_xyz, self.v4):
            self.v4 = rounded_new_xyz
        else:
            raise Exception('Unable to find original vertex value and therefore there is a logic error --> ERROR!!!')


    def calculate_face_centroid(self) -> np.array:
        """
        This function calculates and returns the centroid of the face to use specifically with the re-triangulate
        grammar
        :return: XYZ Coordinates of face centroid
        """
        if self.face_type == 'triangle':
            return np.round(((self.v1 + self.v2 + self.v3) / 3), decimals=self.numDecimals)
        else:
            return np.round(((self.v1 + self.v2 + self.v3 + self.v4) / 4), decimals=self.numDecimals)