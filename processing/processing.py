import numpy as np
import pandas
import os
import glob
import json

def single_spheroid_process(spheroid_frame:pandas.DataFrame,
                            intensity_frame:pandas.DataFrame = None,
                            include_color:bool = False,
                            color_description:str = ''
                            ):
  
    spheroid = {}

    if include_color:

        spheroid['intensity'] = color_description
  
    cells = {}

    for ind in spheroid_frame.index:

        unique_cell = {}

        loc_frame = spheroid_frame.iloc[ind]
    
        unique_cell['x'] = loc_frame['x']
        unique_cell['y'] = loc_frame['y']
        unique_cell['z'] = loc_frame['z']

        # specific to cells where color described

        if include_color:

            assert intensity_frame is not None

            int_frame = intensity_frame.iloc[ind]

            if int_frame['GMM Color'] == -1:
        
                unique_cell['color'] = 'r'

            if int_frame['GMM Color'] == 1:
        
                unique_cell['color'] = 'g'

        cells[int(ind)] = unique_cell
    
    spheroid['cells'] = cells

    return spheroid


def generate_artificial_spheroid(n:int):

    columns = ['x', 'y', 'z']

    data = np.random.rand(n,3)

    Sf = pandas.DataFrame(data = data, columns = columns)

    return single_spheroid_process(Sf)

def batch_spheroid_process(spheroid_folder:str,
                           save_folder:str,
                           include_intensity:bool = False,
                           intensity_folder:str = '',
                           intensity_description:str = ''):

    for file_name in glob.glob(os.path.join(spheroid_folder, '*.csv')):

        spheroid_frame = pandas.read_csv(
            os.path.join(spheroid_folder, file_name))

        _, position, time = file_name.split('_')
        time, _ = time.split('.')

        save_name = position + '_' + time

        if include_intensity:

            list_files = glob.glob(os.path.join(intensity_folder, position + '_' + time))

            assert len(list_files) == 1

            file_name = list_files[0]

            intensity_frame = pandas.read_csv(
                os.path.join(intensity_folder, file_name))

            spheroid = single_spheroid_process(spheroid_frame,
                                        intensity_frame = intensity_frame,
                                        include_color = True,
                                        color_description = intensity_description)

        spheroid = single_spheroid_process(spheroid_frame)

    

    with open(os.path.join(save_folder, save_name + '.json', 'w') as fp:
        json.dump(spheroid, fp, sort_keys=True, indent=4)

    return