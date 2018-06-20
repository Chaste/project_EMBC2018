#!/usr/bin/env python

import subprocess
#mport vtk #uncommenting this causes an "Abort Trap 6"
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np

data_path = '/Users/jmosborne/workspace/ChasteJMOsborne/projects/EMBC2018/test/data/'
working_directory = '/Users/jmosborne/workspace/ChasteWorkingDirectory'

chaste_setup_exe = '/Users/jmosborne/workspace/ChasteJMOsborne/projects/EMBC2018/build/optimised/TestSetupFlowInPipeRunner'
chaste_run_exe = '/Users/jmosborne/workspace/ChasteJMOsborne/projects/EMBC2018/build/optimised/TestRunFlowInPipeRunner'
hemelb_setup_exe = 'env PYTHONPATH=hemelb/Tools/setuptool/:$PYTHONPATH hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui'

def generate_flow_vtus(timestep):
    print 'Turning hemelb results into vtus'
    command = 'python -m hemeTools.converters.ExtractedPropertyUnstructuredGridReader ' + working_directory + 'config.xml ' + working_directory + 'results/Extracted/surface-pressure.xtr ' + working_directory + 'results/Extracted/wholegeometry-velocity.xtr ' + working_directory + 'results/Extracted/surface-traction.xtr ' + working_directory + 'results/Extracted/surface-tangentialprojectiontraction.xtr'
    subprocess.call(command, shell=True)

    output_folders = glob.glob(working_directory + '/results_from_time_*')

    for folder in output_folders:
        folder_timestep = folder.split('results_from_time_')[-1]
        if float(folder_timestep) == timestep:
            shutil.copyfile(working_directory + 'results/Extracted/surface-pressure_2000.vtu', folder + '/surface-pressure.vtu')
            shutil.copyfile(working_directory + 'results/Extracted/surface-traction_2000.vtu', folder + '/surface-traction.vtu')
            shutil.copyfile(working_directory + 'results/Extracted/wholegeometry-velocity_2000.vtu', folder + '/wholegeometry-velocity.vtu')
            shutil.copyfile(working_directory + 'results/Extracted/surface-tangentialprojectiontraction_2000.vtu', folder + '/surface-tangentialprojectiontraction.vtu')

def vtu2stl():
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(working_directory + 'config.vtu')

    # vtkAppendFilter appends one or more datasets together into a single unstructured grid
    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(working_directory + 'config.stl')
    stl_writer.Write()

def run_hemelb_setup():
    # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you
    if  input_folder == 'embryo_plexus/': 
        voxel_size = 1.222e-6 
    else: 
        voxel_size = 0.15e-3

    heme_profile_filename = data_path + input_folder + 'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename + ' --stl ' + working_directory + 'config.stl' + ' --voxel ' + str(voxel_size) + ' --geometry ' + working_directory + 'config.gmy' + ' --xml ' + working_directory + 'config.xml'
    subprocess.call(command, shell=True)

def update_xml_file(iter_num, num_iters):
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    # Add monitoring of incompressibility and convergence
    monitoring = ElementTree.SubElement(root, 'monitoring')
    ElementTree.SubElement(monitoring, 'incompressibility') 
    convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': '1e-3', 'terminate': 'false'})
    ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': '0.01', 'units': 'm/s'})
    
    # Add definition of properties to be extracted
    extr = ElementTree.SubElement(root, 'properties')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(1000), 'file': 'surface-tangentialprojectiontraction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    # Save XML file to disk
    tree.write(filename)


def run_hemelb():
    shutil.rmtree(working_directory + '/results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 4 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)


if __name__=="__main__":

    existing_simulations = ['bifurcation_cut']
        
    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('-o', dest='overwrite_results', action='store_true', help='Allow overwriting of existing results folder.')
    parser.add_argument('-test', dest='test_chaste', action='store_true', help='Dont run HemeLB use a saved traction file to test Chaste code.')
    parser.add_argument('--hemelb_vtu', dest='flow_to_vtu', action='store_true', help='Create VTU files with flow results.')
    parser.add_argument('simulation_type', choices=existing_simulations, help='Type of simulation to be run.')
    parser.add_argument('num_iterations', nargs='?', default=5, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
    parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')
    parser.add_argument('--mesh_scale', dest='mesh_scale', default=1e-3, help='This specifies what to scale the mesh by so that all distances are in meters (defaults to mm).')
    
    # Parse arguments (this will create args.flow_to_vtu etc. variables)
    args = parser.parse_args()

    # Define input folder from simulation type chosen
    input_folder = args.simulation_type + '/'
    
    # Sort out working directory (create, overwrite if allowed, etc.)
    if len(args.output_postfix) == 0:
        working_directory = working_directory + '/'
    else:
        working_directory = working_directory + '-' + args.output_postfix + '/'

    if os.path.isdir(working_directory):
        if args.overwrite_results:
            shutil.rmtree(working_directory)
            os.mkdir(working_directory)
        else:
            raise Exception('Results folder {} exists. Enable results overwriting with -o. '.format(working_directory))
    else:
        os.mkdir(working_directory)

    os.environ['CHASTE_TEST_OUTPUT'] = working_directory

    print 'Working directory = ' + working_directory


    #
    # Step 1: Run preliminary Chaste setup
    #
    mesh_filename = data_path + input_folder + 'config.vtu' # Use proper path concatentation
    xml_filename = data_path + input_folder + 'config.xml' # Use proper path concatentation

    chaste_setup_call = chaste_setup_exe + ' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -mesh_scale ' +  str(args.mesh_scale)
    
    if args.simulation_type == 'bifurcation_cut':
        chaste_setup_call += ' -label_weak_region'

    subprocess.call(chaste_setup_call, shell=True)

    #
    # Step 1.a: Prepare mesh for HemeLB setup tool
    #
    if not args.test_chaste:
        vtu2stl()


    start_time = 0
    duration = 1
    for iter in range(args.num_iterations):
        print 'Iteration ' + str(iter)
        
        if not args.test_chaste:
            #
            # Step 2: Run HemeLB setup
            #
            run_hemelb_setup()

            #
            # Step 3: Use HemeLB to calculate tractions
            #
            update_xml_file(iter, args.num_iterations)
            if args.run_distributed_with_muscle:
                # Run the simulation described by /home/mobernabeu/workspace/ChasteWorkingDirectory/config.xml and /home/mobernabeu/workspace/ChasteWorkingDirectory/config.gmy
                raise NotImplementedError
                # Ensure that /home/mobernabeu/workspace/results/Extracted/surface-traction.xtr' exists for Chaste to carry on.
            else:
                run_hemelb()
        
        #
        # Step 4: Run Chaste with the tractions previously computed
        #
        traction_filename = working_directory + 'results/Extracted/surface-tractions.xtr'
        command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  traction_filename + ' -mesh_scale ' +  str(args.mesh_scale)
        subprocess.call(command, shell=True)
        start_time += duration

        if not args.test_chaste:
            # This has to be done after running Chaste, otherwise it will overwrite the folder
            if args.flow_to_vtu:
                generate_flow_vtus(iter*duration)
                    
            #
            # Step 5: Convert Chaste output (vtu file) into the input of the next iteration (stl file required by hemelb setup tool)
            #
            vtu2stl()

        
            
