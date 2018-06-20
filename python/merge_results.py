#!/usr/bin/env python

import re, glob, shutil, os
from argparse import ArgumentParser

parser = ArgumentParser(description='Prepare the results of a vascular remodelling simulation for Paraview viz')
parser.add_argument('--timestep', dest='timestep', type=float, default='1', help='Separation in time between folders results_from_time_*')
parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')

args = parser.parse_args()

working_directory = 'ChasteWorkingDirectory'
merged_directory = 'MergedResults'

# Sort out merged results directory
if len(args.output_postfix) == 0:
    working_directory = working_directory + '/'
    merged_directory = merged_directory + '/'

else:
    working_directory = working_directory + '-' + args.output_postfix + '/'
    merged_directory = merged_directory + '-' + args.output_postfix + '/'


if os.path.isdir(merged_directory):
    shutil.rmtree(merged_directory)
    os.mkdir(merged_directory)
else:
    os.mkdir(merged_directory)

vtu_files = glob.glob(working_directory + '/*/*.vtu')

filelist = glob.glob(merged_directory +"*.vtu")

for f in filelist:
    os.remove(f)


for pathname in vtu_files:
    filename= os.path.basename(pathname)
    print 'Processing: ' + pathname

    # Remove everything before the results_from_time_
    pattern = working_directory + 'results_from_time_'

    print(pattern)
    replacement  = ''
    tmp = re.sub(pattern, replacement, pathname)

    print(tmp)
    # Remove everything after the /
    seperator = '/'
    tmp2 = tmp.split(seperator, 1)[0]
    print(tmp2)
    file_id =  str(int(float(tmp2)/args.timestep))
    print 'file id', file_id
    
    padded_length = 3
    time = file_id.rjust(padded_length,'0')
    assert len(time) == padded_length # If the assertion trips increase padded_length
    print 'time', time

    #Pad the filename with zeros
    filename_no_extension = re.sub('.vtu', '' , filename)
    print 'filename_no_extension', filename_no_extension
    frame_number = filename_no_extension.split('_', 1)[-1]
    print 'frame_number', frame_number
    
    if (int(frame_number) == 0) and (int(time) > 0):
        continue

    padded_frame_number = frame_number.zfill(6)
    padded_filename = re.sub(frame_number, padded_frame_number, filename)
    
    print(padded_filename)
    
    pattern = '_'
    replacement  = '_' + time + ''
    new_filename = re.sub(pattern, replacement, padded_filename)

    print 'new filename ' + new_filename

    shutil.copy(pathname, merged_directory + new_filename)

    if filename.split('_')[0] == 'results':
        foldername= os.path.dirname(pathname)
        flow_filename = re.sub('results', 'velocity', new_filename)
        try:
            shutil.copy(foldername + '/wholegeometry-velocity.vtu', merged_directory + flow_filename)
        except:
            pass

        pressure_filename = re.sub('results', 'pressure', new_filename)
        try:
            shutil.copy(foldername + '/surface-pressure.vtu', merged_directory + pressure_filename)
        except:
            pass

        traction_filename = re.sub('results', 'tractions', new_filename)
        try:
            shutil.copy(foldername + '/surface-tractions.vtu', merged_directory + traction_filename)
        except:
            pass


