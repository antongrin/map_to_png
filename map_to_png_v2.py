# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 09:59:28 2017

@author: GrinevskiyAS
"""

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import png as png
import AGTools as agt


#----==== My Functions =====---- 
def AG_matrix_from_xyz(il, xl, amp, slice_stats):
    """returns a matrix from xyz data, IL in 0th dim, XL in 1st dim"""
    min_il_real = np.min(il)
    max_il_real = np.max(il)
    min_xl_real = np.min(xl)    
    max_xl_real = np.max(xl)
    ilmin=slice_stats['ilmin']
    ilmax=slice_stats['ilmax']
    xlmin=slice_stats['xlmin']
    xlmax=slice_stats['xlmax']
    dil=slice_stats['dil']
    dxl=slice_stats['dxl']
    il_num=slice_stats['nil']
    xl_num=slice_stats['nxl']
    
    print 'creating 2D array from XYZ.'
    print 'Determined data stats: \tIL {0} - {1}, XL {2} - {3}'.format(min_il_real, max_il_real, min_xl_real, max_xl_real)
    print 'Using data stats: \tIL {0} - {1}, XL {2} - {3}'.format(ilmin, ilmax, xlmin, xlmax)
    print 'Array size: \t{0}x{1}'.format(il_num, xl_num)
    
    datamat=np.zeros((il_num,xl_num))*np.nan

            
    ind_il=((il-ilmin)/dil).astype(int)
    ind_xl=((xl-xlmin)/dxl).astype(int)
    
    for ipoint in xrange(len(il)):
        datamat[ ind_il[ipoint] , ind_xl[ipoint] ] = amp[ipoint]

        
    return datamat    

    
def AG_xyz_from_matrix(ilmat, xlmat, data_out):
    """returns a matrix from xyz data, IL in 0th dim, XL in 1st dim"""
      
    xlm,ilm=np.meshgrid(xlmat,ilmat)
    data_col= data_out.reshape(np.size(data_out),1)     
    ind_notnan= ~(np.isnan(data_col))
    array_to_sv=np.column_stack((ilm.reshape(np.size(ilm),1)[ind_notnan], xlm.reshape(np.size(xlm),1)[ind_notnan],data_col[ind_notnan] ))

    return array_to_sv  


def scale_8bit(input_data, reference_data, direction = 'forward'):
    if direction == 'forward':
        out_data = (((input_data - np.nanmin(reference_data)) / (np.nanmax(reference_data) - np.nanmin(reference_data))) * 255).astype(np.uint8)
    elif direction == 'reverse':
        out_data = np.nanmin(reference_data) + input_data*(np.nanmax(reference_data) - np.nanmin(reference_data))/255.0
        
    return out_data, palette_text


def get_arr_from_png(path):
    print 'Reading file {0}'.format(path)    
    r = png.Reader(path)
    out_data = r.read()
    px_array1 = np.array( map( np.uint8, out_data[2] ))
    px_array = np.reshape(px_array1, (np.shape(px_array1)[0], int(np.shape(px_array1)[1]/4), 4))[:,:,0].astype(float)
    print 'Reading complete.'
    print 'Masking...'
    px_array[empty_points_ffm] = 0*np.nan
    print 'Masking complete.'
    #plt.imshow(abs(px_array.astype(float)-I8.astype(float)), cmap='Spectral_r')#, vmax = 256, vmin = 0)
    
    
    #out_map = np.nanmin(datamat_ffm) +  px_array*(np.nanmax(datamat_ffm) - np.nanmin(datamat_ffm))/255
    out_map = scale_8bit(px_array, datamat_ffm, direction='reverse')

    return out_map




def save_arr_as_png(datamat, filename):
    I8, palette_text = scale_8bit(datamat, datamat, direction='forward')
    I8_a = 255*np.ones_like(I8, dtype='int')
    I8_a[empty_points_ffm] = 0
    I8_exp = np.concatenate((I8[...,None], I8_a[...,None]), axis = 2)
    #I16 = rescale16b(datamat_ffm)
    #I16_a = 65535*np.ones_like(I16, dtype='int')
    #I16_a[empty_points_ffm] = 0
    #I16_exp = np.concatenate((I16[...,None], I16_a[...,None]), axis = 2)
    
    #with open('foo_gray16.png', 'wb') as f:
    #    w = png.Writer(width = I16.shape[1], height=I16.shape[0], bitdepth=16, greyscale=True, alpha = True)
    #    zgray2list = I16_exp.reshape((I16.shape[0], I16.shape[1] * 2))
    #    w.write(f, zgray2list)
    #
    with open(filename, 'wb') as f:
        writer = png.Writer(width = I8.shape[1], height=I8.shape[0], bitdepth=8, greyscale=True, alpha = True)
        zgray2list = I8_exp.reshape((I8.shape[0], I8.shape[1] * 2))
        writer.write(f, zgray2list)
        

 


fnm = r"D:\Projects\Denis\Misc\hrz_for_merging_v2.txt"


data = pd.read_csv(fnm, skiprows=9, delim_whitespace=True, header=None)
print 'Data imported successfully. Shape is {0}'.format(data.shape)
    
il = data.iloc[:, 0].values.astype(int)
xl = data.iloc[:, 1].values.astype(int)
hor_ffm = data.iloc[:, 2].values
hor_f1 = data.iloc[:, 3].values
hor_zd = data.iloc[:, 4].values
hor_reef = data.iloc[:, 5].values

#slice_stats={'ilmin': np.min(il), 'xlmin': np.min(xl), 
#             'xlmax': np.max(xl), 'ilmax': np.max(il),
#             'dil': 1, 'dxl': 1,
#             'nil':((np.max(il)-np.min(il))/1+1).astype(int),'nxl':((np.max(xl)-np.min(xl))/1+1).astype(int)}

slice_stats={'ilmin': 1000, 'ilmax': 5780,
             'xlmin': 5000, 'xlmax': 6800,
             'dil': 1, 'dxl': 1}
slice_stats.update({'nil': int((slice_stats['ilmax'] - slice_stats['ilmin'])/slice_stats['dil']),
                    'nxl': int((slice_stats['xlmax'] - slice_stats['xlmin'])/slice_stats['dxl'])})

datamat_ffm = AG_matrix_from_xyz(il, xl, hor_ffm, slice_stats)
empty_points_ffm = np.isnan(datamat_ffm)
ilmat = np.linspace(slice_stats['ilmin'], slice_stats['ilmax'], slice_stats['nil'])
xlmat = np.linspace(slice_stats['xlmin'], slice_stats['xlmax'], slice_stats['nxl'])

save_arr_as_png(datamat_ffm, 'foo_gray8_v1.png')



####Now converting png file to array and exporting to XYZ
out_map = get_arr_from_png(path = r"D:\Anton\WORK\My_Programs\map_to_png\foo_gray16_ed3.png")

data_to_save = AG_xyz_from_matrix(ilmat, xlmat, out_map)
np.savetxt(r'D:\Anton\WORK\My_Programs\map_to_png\out_matrix_2.txt', data_to_save, delimiter='\t', fmt = ['%6d','%6d','%0.2f'])

plt.imshow(abs(out_map-datamat_ffm), cmap='Spectral_r')#, vmax = 256, vmin = 0)




