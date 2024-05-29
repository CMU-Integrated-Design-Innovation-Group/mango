"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script contains the dataclasses used in the grammar sets for the initial mango design framework. Further grammars
can be added here as long as the top level "GrammarSet" object is used in the dataclass.
"""
from dataclasses import dataclass, field
import time
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
from random import choice, seed
from mango.utils.mango_math import *
from mango.mango_features.mesh_face import MeshFace
from mango.utils.DNA_property_constants import BDNA
from scipy.spatial.transform import Rotation
from mango.mango_features.preserved_regions import PreservedVertex, PreservedEdge

@dataclass
class GrammarSet(object):
    """
    This is a top level class that simply defines the number of decimals to use when rounding & control the random
    seed number of the selection

    Parameters
    ------------
    numDecimals: Integer
    """
    # Precision that can be manually defined. I just use "2" for default since 0.34 is the presumed bp rise in BDNA
    numDecimals: int = 2

    @staticmethod
    def set_seed(seed_number: int):
        """
        Simply seeds the random number for grammar selections
        :param seed_number: User specified seed number
        """
        seed(seed_number)

@dataclass
class TriangulationGrammars(GrammarSet):
    """
    This class controls the grammars used in the generative process. The methods of this class modify the design graph
    of the current design iteration.

    Parameters
    ------------
    grammars_to_not_use: list of grammars to not use if desired. Must match names in grammar_names
    """
    grammars_to_not_use: list = field(default_factory=list)  # This is a potential list of grammars to not use.

    def __post_init__(self):
        self.grammar_names = ['Divide Face', 'Merge Face', 'Extend Vertex', 'Retriangulate Face', 'Edge Rotation']
        # Remove any potential grammars:
        for grammar in self.grammars_to_not_use:
            if grammar not in self.grammar_names:
                raise Exception(f'The grammar {grammar} was not found in the grammar_names list. The only valid grammars '
                                f'you can remove from this set are: {self.grammar_names}')
            else:
                # Here we remove the grammar name from the list such that pick_random_grammar will therefore never
                # pick a value.
                self.grammar_names.remove(grammar)

    @staticmethod
    def set_seed(seed_number: int):
        """ Seeds the random seed """
        seed(seed_number)


    def pick_random_grammar(self) -> str:
        """ Randomly selects any grammar with equal probability """
        return choice(self.grammar_names)


    def call_grammar_function(self, grammar_selected: str, design_space: PolyhedralSpace, extension_value: float = None,
                              override=None) -> tuple[bool, float]:
        """
        This function calls the proper grammar to apply to the active design_space which is passed in

        :param extension_value: If an edge extension is being used we pass in the extension_value
        :param grammar_selected: Whichever grammar was selected to be applied
        :param design_space: The PolyhedralDesign space that is being manipulated
        :param override: Used for testing / image creation... manually apply rules to specific faces
        :return:
        """
        start_time = time.time()
        if grammar_selected == 'Divide Face':
            check = self.divide_face(design_space=design_space, override=override)

        elif grammar_selected == 'Merge Face':
            check = self.merge_face(design_space=design_space, override=override)

        elif grammar_selected == 'Extend Vertex':
            check = self.extend_vertex(design_space=design_space, extension_value=extension_value)

        elif grammar_selected == 'Retriangulate Face':
            check = self.retriangulate_face(design_space=design_space, override=override)

        elif grammar_selected == 'Edge Rotation':
            check = self.edge_rotation(design_space=design_space, rotation_value=extension_value)

        else:
            raise Exception(f'Invalid rule passed in: {grammar_selected}')

        # Now, if check is returned as False, we had an invalid rule applied so we need ot set design_constraint_failure
        end_time = time.time() - start_time
        return check, end_time

    ## ACTUAL GRAMMAR LOGIC BELOW
    def divide_face(self, design_space: PolyhedralSpace, override=None):
        """
        This grammar selects a random face and divides it into two.

        :param design_space: design space that is currently being generated / optimized
        :param override: A tuple of indices to override the random selection
        """
        all_faces_list = design_space.plotting_faces['other']
        face_verts_to_divide = choice(all_faces_list)
        new_faces = []  # Tracking what new faces are created
        if override is not None:
            face_verts_to_divide = override
        face_to_divide = design_space.all_faces[face_verts_to_divide]
        # We divide the face using the method from the mesh_face.py file in mango_features:
        new_vertex, new_triangle1, new_triangle2, divided_edge, og_vertex = face_to_divide.divide_triangular_face()

        # Now, if the new vertex we created already exists in the graph, we will simply return False meaning that the
        # rule has failed:
        faces = np.array(list(design_space.all_faces.keys()), dtype=object)
        vertex_mapping = design_space.renumber_vertices(faces=faces)
        all_verts = np.array(get_all_verts_in_graph(graph=design_space.design_graph, vertex_mapping=vertex_mapping))
        if any(np.all(row == new_vertex) for row in all_verts):
            # If this new_vertex matches a whole row in all_verts (meaning this new point already exists),
            # we simply return False meaning "INVALID RULE APPLICATION"
            ## Note to self: I have never seen this return False, I just thought it would make sense to error
            ##               catch for it in case this is a Fringe case.
            return False

        # Now, due to a requirement of the DAEDALUS and vHelix algorithms, we actually need to "double divide" to
        # maintain a valid mesh. Therefore, we will select a random face that is along the divided edge and also
        # divide it:
        potential_faces = []
        n1 = node_from_xyz(graph=design_space.design_graph, xyz=divided_edge[0])
        n2 = node_from_xyz(graph=design_space.design_graph, xyz=divided_edge[1])
        for f in design_space.plotting_faces['other']:
            face = list(f)
            if (n1 in face and n2 in face) and f != face_verts_to_divide:
                connected_node = next(node for node in set(face) - {n1, n2})
                potential_faces.append((face, connected_node))

        # If we make it here, we then add the new_vertex and begin creating the new faces and edges:
        design_space.add_new_node(x=new_vertex[0], y=new_vertex[1], z=new_vertex[2], terminal=False)
        new_node = node_from_xyz(graph=design_space.design_graph, xyz=new_vertex)
        reciprocal_face_to_divide = []
        if len(potential_faces) != 0:
            # If this is FALSE there are no other faces sharing this edge, we may simply remove the old and add the
            # new faces which we will do outside the if statement since this always happens.
            # IF this is TRUE then we also need to divide out a secondary face:
            reciprocal_face_to_divide, connect = choice(potential_faces)
            # Now we remove reciprocal_face_to_divide and add an edge between the new_vertex and the connect
            design_space.design_graph.add_edge(new_node, connect)
            all_faces_list.remove(tuple(reciprocal_face_to_divide))
            design_space.all_faces.pop(tuple(reciprocal_face_to_divide), None)
            # Add the new face:
            newF1 = (n1, new_node, connect)
            newF2 = (n2, new_node, connect)
            all_faces_list.extend([newF1, newF2])
            new_faces.extend([newF1, newF2])
            design_space.all_faces[newF1] = MeshFace(v1=xyz_from_graph(graph=design_space.design_graph, node=n1),
                                                     v2=xyz_from_graph(graph=design_space.design_graph, node=new_node),
                                                     v3=xyz_from_graph(graph=design_space.design_graph, node=connect),
                                                     numDecimals=self.numDecimals)
            design_space.all_faces[newF2] = MeshFace(v1=xyz_from_graph(graph=design_space.design_graph, node=n2),
                                                     v2=xyz_from_graph(graph=design_space.design_graph, node=new_node),
                                                     v3=xyz_from_graph(graph=design_space.design_graph, node=connect),
                                                     numDecimals=self.numDecimals)

        # Now we repeat the above but with the original face (since this always happens):
        design_space.design_graph.add_edge(new_node, node_from_xyz(graph=design_space.design_graph, xyz=og_vertex))
        # Remove the old edges:
        design_space.design_graph.add_edge(n1, new_node)
        design_space.design_graph.add_edge(n2, new_node)
        design_space.design_graph.remove_edge(n1, n2)
        design_space.edge_divide_nodes[new_node] = [n1, n2]  # Label the edge we divided for other rule uses.
        all_faces_list.remove(face_verts_to_divide)
        design_space.all_faces.pop(face_verts_to_divide, None)
        # Add the new face:
        newF1 = (
            node_from_xyz(graph=design_space.design_graph, xyz=new_triangle1[0]),
            node_from_xyz(graph=design_space.design_graph, xyz=new_triangle1[1]),
            node_from_xyz(graph=design_space.design_graph, xyz=new_triangle1[2])
        )
        newF2 = (
            node_from_xyz(graph=design_space.design_graph, xyz=new_triangle2[0]),
            node_from_xyz(graph=design_space.design_graph, xyz=new_triangle2[1]),
            node_from_xyz(graph=design_space.design_graph, xyz=new_triangle2[2])
        )
        all_faces_list.extend([newF1, newF2])
        design_space.all_faces[newF1] = MeshFace(v1=new_triangle1[0], v2=new_triangle1[1], v3=new_triangle1[2],
                                                 numDecimals=self.numDecimals)
        design_space.all_faces[newF2] = MeshFace(v1=new_triangle2[0], v2=new_triangle2[1], v3=new_triangle2[2],
                                                 numDecimals=self.numDecimals)

        # Finally, we label these faces:
        new_faces += [newF1, newF2]
        if reciprocal_face_to_divide:  # If we added the reciprocal face we update the value of face_verts_to_divide
            face_verts_to_divide = (face_verts_to_divide, tuple(reciprocal_face_to_divide))
        self.label_retri_or_divide_for_backout_rule(design_space=design_space,
                                                    new_faces=new_faces, old_face=face_verts_to_divide)

        # Finally we return True signalling the rule was applied successfully
        return True


    def merge_face(self, design_space: PolyhedralSpace, override=None):
        """
        This rule allows us to "back out" of edge divisions as well as face re-triangulations by checking the labelled
        dictionary and randomly selecting a face to back out of.
        NOTE: This rule currently can lead to very significant design changes as the merge face will merge ALL faces
              to back out. It does NOT require vertices to be coplanar and it will remove a significant amount of
              points and edges depending on which face is randomly selected to be merged back.

        :param design_space: design space that is currently being generated / optimized
        :param override: A tuple of indices to override the random selection
        """
        # First check to see if there are even any potential faces to merge:
        if design_space.merge_label_dict == {}:
            return False

        # Otherwise, we select the LAST key in the dictionary, as we can only back out of the most recent edge division
        # or triangulation. This is for 2 reasons: Computational complexity of checking to ensure a mesh stays water
        # -tight is too high and because a random choice can lead to very drastic design changes which is not ideal imo
        ####recreate_this_face = choice(list(self.merge_label_dict.keys()))
        if override is None:
            recreate_this_face = list(design_space.merge_label_dict.keys())[-1]
        else:
            recreate_this_face = override
        all_faces = design_space.plotting_faces['other']

        # Depending on if the face we are merging was from a division or a re-triangulation:
        if len(recreate_this_face) == 2:
            # If this is length 2, then we are recreating 2 faces from an edge division rule
            delete_these1 = self.get_children_nodes_to_merge(design_space.merge_label_graph,
                                                             start_node=recreate_this_face[0])
            delete_these2 = self.get_children_nodes_to_merge(design_space.merge_label_graph,
                                                             start_node=recreate_this_face[1])
            delete_these = list(set(delete_these1) & set(delete_these2))

            for f in recreate_this_face:
                v1 = xyz_from_graph(graph=design_space.design_graph, node=f[0])
                v2 = xyz_from_graph(graph=design_space.design_graph, node=f[1])
                v3 = xyz_from_graph(graph=design_space.design_graph, node=f[2])
                design_space.all_faces[f] = MeshFace(v1=v1, v2=v2, v3=v3, numDecimals=self.numDecimals)
                all_faces.append(f)
                keys_to_delete = self.is_f_in_keys(merge_label_dict=design_space.merge_label_dict, fa=f)
                for k in keys_to_delete:
                    design_space.merge_label_dict.pop(k, None)

        else:
            # Otherwise, we are only re-creating one face
            delete_these = self.get_children_nodes_to_merge(design_space.merge_label_graph, start_node=recreate_this_face)
            # Also add back the recreate_this_face:
            v1 = xyz_from_graph(graph=design_space.design_graph, node=recreate_this_face[0])
            v2 = xyz_from_graph(graph=design_space.design_graph, node=recreate_this_face[1])
            v3 = xyz_from_graph(graph=design_space.design_graph, node=recreate_this_face[2])

            design_space.all_faces[recreate_this_face] = MeshFace(v1=v1, v2=v2, v3=v3, numDecimals=self.numDecimals)
            all_faces.append(recreate_this_face)
            design_space.merge_label_dict.pop(recreate_this_face, None)

        # Next we need to remove the faces_to_delete
        ## First: Using the merge_label_graph we can find all of the downstream / "children" nodes and verify
        for face in delete_these:
            # First delete the node from the merge_label_graph and remove the key (if it exists)from merge_label_dict
            design_space.merge_label_graph.remove_node(face)
            keys_to_delete = self.is_f_in_keys(merge_label_dict=design_space.merge_label_dict, fa=face)
            for k in keys_to_delete:
                design_space.merge_label_dict.pop(k, None)
            # Next delete the face from all_faces
            design_space.all_faces.pop(face, None)
            if face in all_faces:
                all_faces.remove(face)

        # Next, we need to find any nodes that are no longer a part of the graph. We want to find the unique values of
        # all_faces to find this:
        all_nodes_in_graph = set(design_space.design_graph.nodes)
        unique_values = set(value for tpl in all_faces for value in tpl)
        nodes_to_remove = all_nodes_in_graph - unique_values
        for n in list(nodes_to_remove):
            # Remove any references to these nodes that are currently stored elsewhere:
            design_space.design_graph.remove_node(n)
            if n in design_space.edge_divide_nodes.keys():
                n1, n2 = design_space.edge_divide_nodes[n]
                design_space.design_graph.add_edge(n1, n2)
            design_space.edge_divide_nodes.pop(n, None)
            #design_space.normal_dir_label_dict.pop(n, None)

        # If we make it here, we signal True meaning rule applied successfully
        return True


    def extend_vertex(self, design_space: PolyhedralSpace, extension_value: float = None):
        """
        This grammar randomly selects a vertex that is non-terminal and moves it along an edge by the passed in
        extension_value

        :param design_space: design space that is currently being generated / optimized
        :param extension_value: A value to move the selected node along the edge by (in nm)
        """
        # First we update the list of nonterminal nodes, and if there are not any we return False
        design_space.find_nonterminal_nodes()
        if not design_space.nonterminal_nodes:
            # If there are no nonterminal nodes, we return False meaning no valid rule application
            return False
        # Otherwise, we start by getting all faces referencing this node:
        selected_node = choice(design_space.nonterminal_nodes)
        dir_to_move_towards = choice(list(design_space.design_graph.neighbors(selected_node)))
        # Next we grab the X Y Z values we need and determine the directionality to move in. However, we make a slightly
        # informed decision to move in:
        if design_space.design_graph.nodes[dir_to_move_towards]['terminal']:
            # if the direction we are moving towards is a terminal node, then we can only move in the direction of the
            # nonterminal (or selected) node:
            xyz_moving_towards = xyz_from_graph(graph=design_space.design_graph, node=dir_to_move_towards)
            og_xyz_to_search_for = xyz_from_graph(graph=design_space.design_graph, node=selected_node)
        else:
            # If the direction we are moving towards is also a non-terminal node, then we can choose either direction
            if choice([0, 1]) == 0:
                # If we choose 0, we move towards the dir_to_move_towards
                og_xyz_to_search_for = xyz_from_graph(graph=design_space.design_graph, node=selected_node)
                xyz_moving_towards = xyz_from_graph(graph=design_space.design_graph, node=dir_to_move_towards)
            else:
                # If we choose 1, we move towards the selected_node
                xyz_moving_towards = xyz_from_graph(graph=design_space.design_graph, node=selected_node)
                og_xyz_to_search_for = xyz_from_graph(graph=design_space.design_graph, node=dir_to_move_towards)
                # We also over-write the selected_node value so the positions update properly:
                selected_node = dir_to_move_towards

        # Find the encompassing unit vector to determine which "direction" to extend the edge in:
        u_v = unit_vector(P1=xyz_moving_towards, P2=og_xyz_to_search_for)
        # Update the edge extension by extend_rule_distance (which is passed in
        if extension_value is None:
            extension_value = BDNA.pitch_per_rise  # Set default value to 0.34
        new_value = og_xyz_to_search_for + (extension_value * u_v)
        update_nodal_position(graph=design_space.design_graph, node=selected_node, new_xyz=new_value,
                              numDecimals=self.numDecimals)
        self.update_face_values(design_space=design_space, node=selected_node, old_xyz=og_xyz_to_search_for,
                                new_xyz=new_value)

        # If we get here without issue, we return True meaning it applied successfully
        return True


    def retriangulate_face(self, design_space: PolyhedralSpace, override=None):
        """
        This grammar randomly selects a face and divides it into 3 triangles by placing a new node at the centroid
        and connecting that new node to all 3 edges of the selected triangle.

        :param design_space: design space that is currently being generated / optimized
        :param override: A tuple of 3 face index values to override the random selection
        """
        # First we need to select a face that is not part of the binding regions or NP faces:
        face_verts = choice(design_space.plotting_faces['other'])
        if override is not None:
            face_verts = override
        face_to_triangulate = design_space.all_faces[face_verts]

        # Finding the face centroid:
        centroid = face_to_triangulate.calculate_face_centroid()

        # Next we remove the original faces:
        self.remove_faces(design_space=design_space, face_verts=face_verts)

        # Then we simply add a new vertex here and create edges to the other face verts
        design_space.add_new_node(x=centroid[0], y=centroid[1], z=centroid[2], terminal=False)
        # NOTE: we add edges to node_count-1 because the add_new_node above will add 1 to the node_count.
        design_space.design_graph.add_edge(face_verts[0], design_space.node_count - 1)
        design_space.design_graph.add_edge(face_verts[1], design_space.node_count - 1)
        design_space.design_graph.add_edge(face_verts[2], design_space.node_count - 1)

        # Finally we add the 3 new faces to all_faces and plotting_faces['other']:
        all_faces = design_space.plotting_faces['other']
        new_faces = []
        for v1, v2 in combinations(list(face_verts), 2):
            all_faces.append((design_space.node_count - 1, v1, v2))
            new_faces.append((design_space.node_count - 1, v1, v2))
            v11 = xyz_from_graph(graph=design_space.design_graph, node=design_space.node_count - 1)
            v22 = xyz_from_graph(graph=design_space.design_graph, node=v1)
            v33 = xyz_from_graph(graph=design_space.design_graph, node=v2)
            design_space.all_faces[(design_space.node_count - 1, v1, v2)] = MeshFace(v1=v11, v2=v22, v3=v33,
                                                                                     numDecimals=self.numDecimals)
        # Next we label the new face:
        self.label_retri_or_divide_for_backout_rule(design_space=design_space, new_faces=new_faces, old_face=face_verts)

        # If we get here without issue, we return True meaning it applied successfully
        return True

    def edge_rotation(self, design_space: PolyhedralSpace, rotation_value: float = None):
        """
        This grammar randomly selects an edge containing a nonterminal node and rotates it along any face by which
        the edge exists. The face is randomly selected.

        :param design_space: design space that is currently being generated / optimized
        :param rotation_value: A value to rotate the edge along the face by (degrees)
        """
        # First we update the list of nonterminal nodes, and if there are not any we return False
        design_space.find_nonterminal_nodes()
        if not design_space.nonterminal_nodes:
            # If there are no nonterminal nodes, we return False meaning no valid rule application
            return False
        # Now, we can't move into the normal direction of a labelled edge division, so we remove those:
        selected_node = choice(design_space.nonterminal_nodes)
        selected_node_connection = []
        # Now we check if the selected_node has already had a labelled normal-direction move applied:
        if selected_node in design_space.edge_divide_nodes.keys():
            # If we select a node that was a part of an edge division, we need to select its non-collinear
            # points to create the plane label.
            connected_nodes = list(design_space.design_graph.neighbors(selected_node))
            n1, n2, n3 = find_non_collinear_nodes(graph=design_space.design_graph, node=selected_node,
                                                  potential_nodes=connected_nodes)
            if n1 == -1:
                # In this case, we could not find 3 non-collinear nodes forming the plane
                return False
            v = [xyz_from_graph(graph=design_space.design_graph, node=selected_node)]
            for i in [n2, n3]:  # We don't use n1 because that is the same as selected_node
                v.append(xyz_from_graph(graph=design_space.design_graph, node=i))
                selected_node_connection.append(xyz_from_graph(graph=design_space.design_graph, node=i))
            local_plane = LocalPlane(verts=np.array(v), og_point=xyz_from_graph(graph=design_space.design_graph,
                                                                                node=selected_node),
                                     divide_edge_rule_used=True)
        else:
            # Otherwise we must create a local plane label for the re-triangulation such:
            connected_nodes = list(design_space.design_graph.neighbors(selected_node))
            v = [xyz_from_graph(graph=design_space.design_graph, node=selected_node)]
            for i in connected_nodes:
                v.append(xyz_from_graph(graph=design_space.design_graph, node=i))
                selected_node_connection.append(xyz_from_graph(graph=design_space.design_graph, node=i))
            # Now with all faces that are connected to this node, we create a local plane and store in the dict:
            local_plane = LocalPlane(verts=np.array(v), og_point=xyz_from_graph(graph=design_space.design_graph,
                                                                                node=selected_node))

        # With the local_plane defined, we are going to move by extend_rule_distance into one of the normal dirs:
        normal_dir1 = np.array([local_plane.A, local_plane.B, local_plane.C])

        # Check both rotations about the face:
        rotation_matrix1 = Rotation.from_rotvec(normal_dir1 * rotation_value).as_matrix()
        rotation_matrix2 = Rotation.from_rotvec(normal_dir1 * (2*np.pi - rotation_value)).as_matrix()
        selected_rotation = choice([rotation_matrix1, rotation_matrix2])

        cur_xyz = xyz_from_graph(graph=design_space.design_graph, node=selected_node)
        connected_node_xyz = choice(selected_node_connection)  # Doesn't matter which we pick, just need the edge:
        edge_vector = connected_node_xyz - cur_xyz
        unit_v = edge_vector / np.linalg.norm(edge_vector)

        # Based on local_plan, we check if the unit vector is not none and overwrite otherwise:
        if local_plane.unit_bisector_vector is not None:
            unit_v1 = local_plane.unit_bisector_vector
            unit_v2 = -1 * local_plane.unit_bisector_vector
            unit_v = choice([unit_v1, unit_v2])

        new_value = cur_xyz + np.dot(selected_rotation, unit_v)
        # Then we update the positions for this node to all lists:
        update_nodal_position(graph=design_space.design_graph, node=selected_node, new_xyz=new_value,
                              numDecimals=self.numDecimals)
        self.update_face_values(design_space=design_space, node=selected_node, old_xyz=cur_xyz, new_xyz=new_value)

        # Signal that the rule has applied successfully:
        return True

    # HELPER FUNCTIONS FOR ABOVE GRAMMARS
    @staticmethod
    def label_retri_or_divide_for_backout_rule(design_space: PolyhedralSpace, new_faces: list, old_face: tuple) -> None:
        """
        This function is called whenever the retriangulation or division rule is applied such that we can
        back out of these geometries
        """
        # First we must check if the face we divided / merged
        design_space.merge_label_dict[old_face] = new_faces
        design_space.merge_label_graph.add_nodes_from(new_faces)  # Add new faces as nodes
        # Next we add to the directed graph:
        if len(old_face) == 2:
            # If old_face is length 2 that means two nodes are connected to "children" that we need to track:
            edges_to_add = [(old_face[0], target_node) for target_node in new_faces]
            edges_to_add2 = [(old_face[1], target_node) for target_node in new_faces]
            design_space.merge_label_graph.add_edges_from(edges_to_add)
            design_space.merge_label_graph.add_edges_from(edges_to_add2)
        else:
            edges_to_add = [(old_face, target_node) for target_node in new_faces]
            design_space.merge_label_graph.add_edges_from(edges_to_add)


    @staticmethod
    def get_children_nodes_to_merge(graph, start_node):
        downstream_nodes = set()
        stack = [start_node]

        while stack:
            current_node = stack.pop()
            downstream_nodes.add(current_node)

            # Get successors and add them to the stack
            successors = graph.successors(current_node)
            stack.extend(successors)
        downstream_nodes.discard(start_node)
        return list(downstream_nodes)


    @staticmethod
    def is_f_in_keys(merge_label_dict: dict, fa: tuple):
        """
        This function isn't really math but it searches a dictionary to find if a face is found within the dictionary

        :param merge_label_dict: Dictionary containing merge labels
        :param fa: Face we are looking for
        :return:
        """
        ke = []
        for v in merge_label_dict.keys():
            if len(v) == 2:
                if fa in list(v):
                    ke.append(v)
        return ke


    @staticmethod
    def update_face_values(design_space: PolyhedralSpace, node: int, old_xyz: np.array, new_xyz: np.array) -> None:
        """
        This function will update all of the face values in the current mesh with the new nodal values
        :param design_space: Design space we are updating
        :param node: Node number to update
        :param old_xyz: Old XYZ Position
        :param new_xyz: New XYZ Position
        :return: Nothing, just updates values
        """
        # We loop over all the faces which contain vertex information and update these values:
        for k, v in design_space.all_faces.items():
            # If the node is in the list of vertices for the face, we need to update it's mesh face:
            if node in list(k):
                v.update_vertex_values(old_xyz=old_xyz, new_xyz=new_xyz)


    @staticmethod
    def remove_faces(design_space: PolyhedralSpace, face_verts: tuple) -> None:
        """
        This function removes a face from the self.all_faces and self.plotting_faces['other']
        :param design_space: Design space being modified
        :param face_verts: Tuple of vertices to search for
        :return: Nothing, just updates values
        """
        design_space.all_faces.pop(face_verts, None)
        all_faces = design_space.plotting_faces['other']
        all_faces.remove(face_verts)
        design_space.plotting_faces['other'] = all_faces


@dataclass
class ParallelepipedGrammars(GrammarSet):
    """
    This class controls the grammars used in the generative process that modify the bounding box conditions.

    Parameters
    ------------
    cell_type: string name of a desired cell_type if input (e.g., triclinic, monoclinic, etc.).
    """
    cell_type: str = 'none'   # This determines which rules can be applied to which directions DesignConstraints used to
                              # hold cell walls constant.


    def __post_init__(self):
        # Which grammars we can call depends on the cell_type:
        if self.cell_type == 'triclinic':
            self.grammar_names = ['Vary_a', 'Vary_b', 'Vary_c',
                                  'Rotate_alpha', 'Rotate_beta', 'Rotate_gamma']

        elif self.cell_type == 'monoclinic':
            self.grammar_names = ['Vary_a', 'Vary_b', 'Vary_c',
                                  'Rotate_beta']

        elif self.cell_type == 'orthorhombic':
            self.grammar_names = ['Vary_a', 'Vary_b', 'Vary_c']

        elif self.cell_type == 'tetragonal':
            self.grammar_names = ['Vary_a', 'Vary_c']

        elif self.cell_type == 'rhombohedral' or self.cell_type == 'trigonal':
            self.grammar_names = ['Vary_a',
                                  'Rotate_alpha']

        elif self.cell_type == 'hexagonal':
            self.grammar_names = ['Vary_a', 'Vary_c']

        elif self.cell_type == 'isometric' or self.cell_type == 'cubic':
            self.grammar_names = ['Vary_a']
        else:
            raise Exception('Invalid parallelepiped volume calculation')


    @staticmethod
    def set_seed(seed_number: int):
        """ Sets the random seed """
        seed(seed_number)

    def pick_random_grammar(self) -> str:
        """ Randomly selects a grammar with equal probabilty """
        return choice(self.grammar_names)


    def call_grammar_function(self, grammar_selected: str, design_space: PolyhedralSpace, extension_value: float = None,
                              rotation_value: float = 1) -> tuple[bool, float]:
        """
        This function calls the proper grammar to apply to the active design_space which is passed in

        :param extension_value: How far to extend the grammar selected by
        :param grammar_selected: String of whatever grammar was selected to be applied
        :param design_space: The PolyhedralDesign space that is being manipulated
        :param rotation_value: How much to rotate an angle of the bounding box by
        :return:
        """
        start_time = time.time()
        grammar_rule, geometry = grammar_selected.split('_')

        if grammar_rule == 'Vary':
            check = self.extend_or_shorten_edge(design_space=design_space, extension_value=extension_value,
                                                param=geometry)

        elif grammar_rule == 'Rotate':
            check = self.rotate_edge(design_space=design_space, rotation_value=rotation_value, param=geometry)

        else:
            raise Exception(f'Invalid rule passed in: {grammar_selected}')

        # Now, if check is returned as False, we had an invalid rule applied so we need ot set design_constraint_failure
        end_time = time.time() - start_time
        return check, end_time

    # ACTUAL GRAMMAR LOGIC BELOW
    def extend_or_shorten_edge(self, design_space: PolyhedralSpace, extension_value: float, param: str) -> bool:
        """
        This function calls to extend (or shorten using a negative extension_value) one of the bounding box parameters

        :param extension_value: How far to extend the grammar selected by
        :param param: Which bounding box value to change (a, b, c)
        :param design_space: The PolyhedralDesign space that is being manipulated
        :return: Bool dictating if function ran properly
        """
        # This grammar will effectively change the value of 'a' which will modify the PreservedRegion vertex values
        # and therefore changes the size of the structure compliant to the constraints of the problem. The constraints
        # are very important for this type of optimization to narrow the design search to a reasonable search

        ## Step 1: Find the current value of the param:
        new_bounding_box = design_space.bounding_box
        curVal = round(getattr(new_bounding_box, param), self.numDecimals)

        # Randomly choose to either shorten or lengthen the selected param:
        random_choice = choice([-1, 1])  # -1 will shorten, 1 will lengthen the value
        curVal += round((random_choice * extension_value), self.numDecimals)

        # Set the new value to design_space.bounding_box.param
        for param_to_update in design_space.bounding_box.param_map[param]:
            setattr(new_bounding_box, param_to_update, curVal)

        # Now update:
        if self.update_design_space(design_space=design_space, new_bounding_box=new_bounding_box):
            return True  # If the design updated properly we send "True" back signalling rule aplied
        else:
            return False  # Otherwise signal the rule could not be applied


    def rotate_edge(self, design_space: PolyhedralSpace, rotation_value: float, param: str) -> bool:
        """
        This function calls to rotate one of the bounding box parameters

        :param rotation_value: How far to rotate the selected bounding box parameter
        :param param: Which bounding box value to change (alpha, beta, gamma)
        :param design_space: The PolyhedralDesign space that is being manipulated
        :return: Bool dictating if function ran properly
        """
        # This follows a very similar logic to extend_or_shorten where here we "open" or "close" an angle correspondin
        # to a specific alpha beta or gamma of a unit cell

        ## Step 1: Find the current value of the param:
        new_bounding_box = design_space.bounding_box
        curVal = getattr(new_bounding_box, param)

        # Randomly choose to either shorten or lengthen the selected param:
        random_choice = choice([-1, 1])  # -1 will shorten, 1 will lengthen the value
        curVal += (random_choice * rotation_value)  # Should be in RADIANS (adding or subtracting degrees)

        # Set the new value to design_space.bounding_box.param
        for param_to_update in design_space.bounding_box.param_map[param]:
            setattr(new_bounding_box, param_to_update, curVal)

        # Now update:
        if self.update_design_space(design_space=design_space, new_bounding_box=new_bounding_box):
            return True  # If the design updated properly we send "True" back signalling rule aplied
        else:
            return False  # Otherwise signal the rule could not be aplpied

    # HELPER FUNCTIONS FOR ABOVE GRAMMARS
    def update_design_space(self, design_space: PolyhedralSpace, new_bounding_box) -> bool:
        """
        This function actually updated the bounding box hyperparmeters so that the volume / other values can be
        properly updated
        :param new_bounding_box: New boudning box to be used in the design iteration
        :param design_space: The PolyhedralDesign space that is being manipulated
        :return: Bool dictating if function ran properly
        """
        # Update midpoint values to compare
        # NOTE (WIP): In the future I need to change this logic to allow for more "arbitrary" searches. Right now
        #             I am just updating the face midpoint
        # First we get the current midpoints:
        previous_vertices = new_bounding_box.shape.calculate_vertices_and_edges()
        current_midpoints = new_bounding_box.shape.get_midpoint_of_each_face(vertices=previous_vertices)

        new_bounding_box.update_shape()
        new_vertices = new_bounding_box.shape.calculate_vertices_and_edges()
        new_bounding_box.shape.update_face_equations()
        new_bounding_box.shape.update_volume()
        new_midpoints = new_bounding_box.shape.get_midpoint_of_each_face(vertices=new_vertices)
        midpoint_map = {}
        for i, j in zip(current_midpoints, new_midpoints):
            midpoint_map[i] = j

        # Update to new bounding box:
        design_space.bounding_box = new_bounding_box


        ## Also need to update the design graph preserved region values potentially (if they are on the bounding box):
        ## NOTE (WIP): Need to scale things up properly? Like I think everything need to "proportionally" scale
        for preserved_region in design_space.preserved:
            if isinstance(preserved_region, PreservedVertex):
                curPoint = preserved_region.v1
                if tuple(curPoint) in midpoint_map.keys():
                    newPoint = midpoint_map[tuple(curPoint)]
                    preserved_region.v1 = np.array([newPoint[0], newPoint[1], newPoint[2]])

                    # Now we need to update this nodal position in design graph:
                    curNode = node_from_xyz(graph=design_space.design_graph, xyz=curPoint)
                    update_nodal_position(graph=design_space.design_graph, node=curNode,
                                          new_xyz=preserved_region.v1, numDecimals=self.numDecimals)
                    update_vertex_in_all_faces(all_faces=design_space.all_faces, node=curNode,
                                               new_xyz=preserved_region.v1)


            elif isinstance(preserved_region, PreservedEdge):
                print('WIP: Need to implement this')
        # After this all occurs, we send True signalling 'grammar applied successfully'
        return True

@dataclass
class LocalPlane(object):
    """
    This class is not used directly by the user, it is used to check geometry of a generated design
    """
    verts: np.array
    og_point: np.array  # Original vertex / significant point we need to tract
    divide_edge_rule_used: bool = False

    def __post_init__(self):
        ### The plane we create for 3 coplanar points from a retriangulation rule is different:
        # Select the first 3 points and verify that they are NOT co-linear
        p1, p2, p3 = self.verts[:3, :]

        # Find two vectors lying on the plane
        v1 = p3 - p1
        v2 = p2 - p1

        # Calculate the normal vector using the cross product
        normal_vector = np.cross(v1, v2)
        self.A, self.B, self.C = np.round(normal_vector, 3)
        self.D = np.round(np.dot(normal_vector, p3), 3)

        # If we create a plan using the divide_edge rule, we also set the angle:
        if self.divide_edge_rule_used:
            bisector_vector = (v1 / np.linalg.norm(v1)) + (v2 / np.linalg.norm(v2))
            # Normalize the bisector vector to obtain the unit vector
            self.unit_bisector_vector = bisector_vector / np.linalg.norm(bisector_vector)

        else:
            self.unit_bisector_vector = None

@dataclass
class CustomGrammarSet(GrammarSet):
    """
    This class is not used directly by the user, but rather it is created if the user specifies more than 1 GrammarSets
    """
    grammar_sets: list = field(default_factory=list)

    def __post_init__(self):
        # First need to create a list of grammar_names:
        self.grammar_names = []
        self.grammar_map = {}
        for s in self.grammar_sets:
            # We add all grammar names from all grammar sets as potential rules to pick
            for name in s.grammar_names:
                self.grammar_names.append(name)
                # We also add this name to the grammar_map so we know which grammar functions to call:
                self.grammar_map[name] = s

    @staticmethod
    def set_seed(seed_number: int):
        seed(seed_number)


    def pick_random_grammar(self) -> str:
        # Simply select a random rule
        return choice(self.grammar_names)


    def call_grammar_function(self, grammar_selected: str, design_space: PolyhedralSpace,
                              extension_value: float = None) -> tuple[bool, float]:
        # First get the GrammarSet for the grammar_selected:
        grammar_set = self.grammar_map[grammar_selected]

        # Now call the call_grammar_function:
        grammar_applied_successfully, grammar_time = grammar_set.call_grammar_function(
            grammar_selected=grammar_selected,
            design_space=design_space,
            extension_value=extension_value)

        # Finally just "pass along" the grammar_applied_successfully and grammar_time functions:
        return grammar_applied_successfully, grammar_time

