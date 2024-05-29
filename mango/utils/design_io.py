"""
A.J. Vetturini
IDIG and MMBL
Carnegie Mellon University

This script is for reading in / outputting design files. For example, exporting to PLY file to convert a design to a
DNA file via DAEDALUS2.
"""
from mango.design_spaces.polyhedral_design_space import PolyhedralSpace
import os
from trimesh.exchange.ply import export_ply
from trimesh.exchange.obj import export_obj
from mango.utils.mango_math import *
from mango.utils.DNA_property_constants import BDNA
import subprocess
import dill
import shutil

def export_design_to_ply(design: PolyhedralSpace, savepath: str, savename_no_extension: str):
    """
    NOT INTENDED FOR END USER USE. PLEASE SEE export_DNA_design for the user-facing function

    This function exports a given PolyhedralSpace to a .ply file for the use in a automated scaffold routing algorithm

    :param design: PolyhedralSpace design space that is to be converted to a PLY file
    :param savepath: Filepath where to save the ply file
    :param savename_no_extension: Filename of the ply file to use
    :return: Directory path of where the file was saved
    """
    # First create the PLY (mesh) file:
    design.create_trimesh_mesh()  # Re-create the trimesh mesh just in case
    mesh = design.mesh_file
    mesh.fix_normals()
    export = export_ply(mesh=mesh, encoding='ascii')

    # Create the savepath:
    save_name = savename_no_extension + '.ply'
    # Check if the path exists:
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    writeDir = os.path.join(savepath, save_name)
    with open(writeDir, "wb") as f:
        f.write(export)
    return writeDir

def export_design_to_obj(design: PolyhedralSpace, savepath: str, savename_no_extension: str):
    """
    NOT INTENDED FOR END USER USE. PLEASE SEE export_DNA_design for the user-facing function

    This function exports a given PolyhedralSpace to a .obj file for the use in a automated scaffold routing algorithm

    :param design: PolyhedralSpace design space that is to be converted to a obj file
    :param savepath: Filepath where to save the obj file
    :param savename_no_extension: Filename of the obj file to use
    :return: Directory path of where the file was saved
    """
    # First create the PLY (mesh) file:
    design.create_trimesh_mesh()  # Re-create the trimesh mesh just in case
    mesh = design.mesh_file
    mesh.fix_normals()
    export = export_obj(mesh=mesh)

    # Create the savepath:
    save_name = savename_no_extension + '.obj'
    # Check if the path exists:
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    writeDir = os.path.join(savepath, save_name)
    with open(writeDir, "w") as f:
        f.write(export)
    return writeDir


def all_same(strings):
    """ Simple function to verify that the filename is unique when doing a mass export """
    # Check if the list is empty
    if not strings:
        return False

    # Take the first string as a reference
    reference_string = strings[0]

    # Compare each string to the reference string
    for string in strings[1:]:
        if string != reference_string:
            return False

    # If all strings are the same, return True
    return True


def export_list_of_designs_to_ply(design_list: list, savepath: str, list_of_savenames_no_extensions: list):
    """
    NOT INTENDED FOR END USER USE. PLEASE SEE mass_design_export for the user-facing function

    This function will export a whole list of designs such that large amounts of data can be gathered.

    :param design_list: List of PolyhedralSpace design space's that are to be converted to PLY files
    :param savepath: Filepath where to save the ply file
    :param list_of_savenames_no_extensions: List of filenames to save the files as (should be length of design_list)
    """
    if len(design_list) != len(list_of_savenames_no_extensions):
        raise Exception('Lengths of lists do NOT match -> Can NOT export')
    if all_same(list_of_savenames_no_extensions):
        raise Exception('Duplicate savenames found -> Can NOT export')

    for design, savename in zip(design_list, list_of_savenames_no_extensions):
        export_design_to_ply(design=design,
                             savepath=savepath,
                             savename_no_extension=savename)


def export_DNA_design(automated_scaffold_executable_path: str, design: PolyhedralSpace, export_path: str,
                      savename_no_extension: str, DNA_geometry=BDNA, scaffold_sequence_filepath: str = None,
                      ply_or_obj: str = 'ply') -> int:
    """
    This is the top level function a user should call which will export a single design. The export_design_to_ply
    is called by this function.

    :param design: PolyhedralSpace design space that is to be converted to a PLY file
    :param automated_scaffold_executable_path: Filepath to automated algorithm (e.g. DAEDALUS2.exe path)
    :param export_path: Filepath where to save the DNA design files to
    :param savename_no_extension: User defined name to save the design as
    :param DNA_geometry: Do not change
    :param scaffold_sequence_filepath: Path to sequence file to use, default is M13
    :param ply_or_obj: Ply will export a PLY file, obj will export an OBJ file
    """
    # First we need the minimal edge length in the design:
    all_edges = calculate_design_edge_lengths(graph=design.design_graph)
    min_edge_length = min(all_edges)  # min(all_edges) is in units of nm, need nucleobases:
    min_edge_length_nb = int(min_edge_length / DNA_geometry.pitch_per_rise)
    # With the min_edge_length found, we can create a PLY file in the specified path:
    if ply_or_obj == 'ply':
        ply_path = export_design_to_ply(design=design, savepath=export_path,
                                        savename_no_extension=savename_no_extension)
    else:
        ply_path = export_design_to_obj(design=design, savepath=export_path,
                                        savename_no_extension=savename_no_extension)
    # Then we will call the automated algorithm.
    ## First check if user specified a sequence to use:
    SEQ_PATH = './M13.txt'
    if scaffold_sequence_filepath is not None:
        SEQ_PATH = scaffold_sequence_filepath
    if 'DAEDALUS' in automated_scaffold_executable_path:
        #command = f"{automated_scaffold_executable_path} {export_path} {ply_path} {SEQ_PATH} 1 2 0 {min_edge_length_nb} 0.0 m > output.log 2>&1"
        ## NOTE TO SELF: 1 2 1 is outputting the correct scale in design... interesting, 0 is scaling up. FORTRAN code
        ## is a bit rather non-obvious as to why this scales up, but these designs simulated look good
        command = f"{automated_scaffold_executable_path} {export_path} {ply_path} {SEQ_PATH} 1 2 1 {min_edge_length_nb} 0.0 m > output.log 2>&1"
    elif 'TALOS' not in automated_scaffold_executable_path:
        command = f"{automated_scaffold_executable_path} {export_path} {ply_path} {SEQ_PATH} 3 2 0 {min_edge_length_nb} 0.0 m > output.log 2>&1"
    else:
        raise Exception('Only DAEDALUS and TALOS are supported for the executable path.')
    # Run:
    try:
        subprocess.run(command, shell=True, text=True, check=True, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError:
        pass
    # We return this value although most of the time it isn't really needed, just useful for mass export:
    return min_edge_length_nb


def mass_design_export(automated_scaffold_executable_path: str, export_path: str, aj1_filepath: str,
                       scaffold_sequence_filepath: str, DNA_geometry=BDNA, ply_or_obj: str = 'ply'):
    """
    This is the top level function a user should call which will export a single design. The export_design_to_ply
    is called by this function.

    :param aj1_filepath: Filepath to .aj1 file that was exported from a multi-objective optimization process
    :param automated_scaffold_executable_path:  Filepath to automated algorithm (e.g. DAEDALUS2.exe path)
    :param export_path: Filepath where to save the DNA design files to
    :param DNA_geometry: Do not change
    :param scaffold_sequence_filepath: Path to sequence file to use, default is M13
    """
    print('Beginning mass design export, note that this will take a while. There may be warnings printed below but '
          'the process is likely working properly!')
    # Check if path exists:
    if not os.path.exists(export_path):
        os.makedirs(export_path)
    save_json_path = os.path.join(export_path, 'JSON_files')
    foldername = 'PLY_Files' if ply_or_obj == 'ply' else 'OBJ_files'
    file_extension = '.ply' if ply_or_obj == 'ply' else '.obj'
    save_ply_path = os.path.join(export_path, foldername)
    if not os.path.exists(save_json_path):
        os.makedirs(save_json_path)
    if not os.path.exists(save_ply_path):
        os.makedirs(save_ply_path)

    with open(aj1_filepath, 'rb') as f:
        data = dill.load(f)
    all_files = []
    failed_conversion_designs = []
    for datapoint, design in data.MOSA_archive.items():
        # First we need to create a unique name identifier:
        fname = 'Design_' + str(datapoint[0]) + '_' + str(datapoint[1])
        ply_name = os.path.join(export_path, fname) + file_extension
        # Call the export function:
        nb_value = export_DNA_design(automated_scaffold_executable_path, design, export_path, fname,
                                     scaffold_sequence_filepath=scaffold_sequence_filepath)
        # After running the export we need to clean up the folder a bit:
        folder_path = fname + '_DX_' + str(nb_value) + 'bp'
        savename = folder_path + '.json'
        json_file = os.path.join(os.path.join(export_path, folder_path), savename)
        try:
            shutil.copy(json_file, save_json_path)
            shutil.move(ply_name, os.path.join(export_path, foldername))
            # Remove the various DAEDALUS files to preserve space:
            shutil.rmtree(os.path.join(export_path, folder_path))

            # Append tuple to all_files to write out the SNUPPI input file:
            all_files.append(('H', os.path.join('JSON_files', savename)))  # H is for honeycomb since DAEDALUS is always a honeycomb

        except FileNotFoundError:
            print(f'NOTE: DAEDALUS failed to convert the design: {fname}')
            failed_conversion_designs.append(savename)


    # Next we need to create the SNUPPI input file
    with open(os.path.join(export_path, 'Input.txt'), 'w') as f:
        for file in all_files:
            f.write(str(file[0]) + '\t' + str(file[1]) + '\n')
    with open(os.path.join(export_path, 'README'), 'w') as f:
        f.write('Place the JSON_files folder into the SNUPI root directory (shown below) and simply run SNUPI. \n'
                'The input file has been formatting to sequentially run all designs in the JSON_files filepath. \n'
                'Note that the default SNUPI input configuration is not ideal for this process as it will constantly try and populate windows and images, so please check that out.\n')
        f.write('SNUPI_(version #)\n')
        f.write('-- SNUPI.exe\n')
        f.write('-- OUTPUT\n')
        f.write('-- (Place JSON_files folder here)\n\n')
        f.write('Below is a list of all failed conversions for logging purposes. These designs failed to convert via DAEDALUS2: \n')
        for file in failed_conversion_designs:
            f.write(file + '\n')
