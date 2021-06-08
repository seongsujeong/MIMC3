import scipy.io as sio
import gma
import glob
import numpy as np
import os
import matplotlib.pyplot as plt



class vmap:
    def __init__(self,path=None):
        self.path=path
        self._core_data={
            'x':None,
            'y':None,
            'vx':None,
            'vy':None,
            'spd':None,
            'ex':None,
            'ey':None,
            'qual':None,
            'flagcp':None,
            'meta':None
        }

        self.dict_dtype={
            'x':['float64','GMA'],
            'y':['float64','GMA'],
            'vx':['float32','GMA'],
            'vy':['float32','GMA'],
            'ex':['float32','GMA'],
            'ey':['float32','GMA'],
            'qual':['float32','GMA'],
            'flagcp':['uint8','GMA'],
            'meta':['text','txt']
        }
    
    def load_gma(self,path_vmap=None,field=None): #Load vmap as in GMA format
        if path_vmap is not None:
            self.path=path_vmap

        if field==None:
            list_field=self.dict_dtype.keys() #load all fields defined in the class
        else:
            list_field=field
        
        for field_of_interest in list_field:
            filename_in=glob.glob('{}/*_{}.{}'.format(self.path,field_of_interest,self.dict_dtype[field_of_interest][1]))
            if len(filename_in)==1:
                if self.dict_dtype[field_of_interest][1]=='txt':    #parse _meta.txt
                    with open(filename_in[0],'r') as fin:
                        lines_meta=fin.read().split('\n')
                        if lines_meta[-1]=='\n':
                            lines_meta=lines_meta[:-1]
                    self._core_data['meta']={}
                    for line in lines_meta:
                        seg_line=line.split('=')
                        if 'cp_offset_int' in seg_line[0]:
                            try:
                                self._core_data[field_of_interest][seg_line[0]]=int(seg_line[1])
                            except:
                                self._core_data[field_of_interest][seg_line[0]]=float(seg_line[1])
                        elif line=='': #blank line; do nothing
                            pass
                        else:
                            self._core_data[field_of_interest][seg_line[0]]=seg_line[1]
                else:
                    #treat the input file as GMA and load as such
                    self._core_data[field_of_interest]=gma.read(filename_in[0],dtype=self.dict_dtype[field_of_interest][0])
            else:
                print('ERROR: vmap.vmap.load() - Cannot find file: {}'.format(filename_in))
            

    @property
    def x(self):
        if self._core_data['x'] is None:
            try:
                self.load_gma(field=['x'])
            except:
                print('ERROR: vmap.x()')
        return self._core_data['x']
    
    @property
    def y(self):
        if self._core_data['y'] is None:
            try:
                self.load_gma(field=['y'])
            except:
                print('ERROR: vmap.y()')
        return self._core_data['y']

    @property
    def vx(self):
        if self._core_data['vx'] is None:
            try:
                self.load_gma(field=['vx'])
            except:
                print('ERROR: vmap.vx()')
        return self._core_data['vx']
    
    @property
    def vy(self):
        if self._core_data['vy'] is None:
            try:
                self.load_gma(field=['vy'])
            except:
                print('ERROR: vmap.vy()')
        return self._core_data['vy']

    @property
    def spd(self):
        if self._core_data['spd'] is None:
            if self._core_data['vx'] is None:
                try:
                    self.load_gma(field=['vx'])
                except:
                    print('ERROR: vmap.spd() when loading vx')

            if self._core_data['vy'] is None:
                try:
                    self.load_gma(field=['vy'])
                except:
                    print('ERROR: vmap.spd() when loading vy')

        self._core_data['spd']=np.sqrt(self.vx**2+self.vy**2)

        return self._core_data['spd']

    @property
    def ex(self):
        if self._core_data['ex'] is None:
            try:
                self.load_gma(field=['ex'])
            except:
                print('ERROR: vmap.ex()')
        return self._core_data['ex']
    
    @property
    def ey(self):
        if self._core_data['ey'] is None:
            try:
                self.load_gma(field=['ey'])
            except:
                print('ERROR: vmap.ey()')
        return self._core_data['ey']

    @property
    def qual(self):
        if self._core_data['qual'] is None:
            try:
                self.load_gma(field=['qual'])
            except:
                print('ERROR: vmap.qual()')
        return self._core_data['qual']

    @property
    def flagcp(self):
        if self._core_data['flagcp'] is None:
            try:
                self.load_gma(field=['flagcp'])
            except:
                print('ERROR: vmap.flagcp()')
        return self._core_data['flagcp']

    @property
    def meta(self):
        if self._core_data['meta'] is None:
            try:
                self.load_gma(field=['meta'])
            except:
                print('ERROR: vmap.meta()')
        return self._core_data['meta']

    def adjust(self):
        vxcp=self.vx[self.flagcp!=0]
        vycp=self.vy[self.flagcp!=0]
        
        mvxcp=np.nanmean(vxcp)
        mvycp=np.nanmean(vycp)

        print('bias on cp: vx={}, vy={}'.format(mvxcp,mvycp))

        self._core_data['vx']=self._core_data['vx']-mvxcp
        self._core_data['vy']=self._core_data['vy']-mvycp
        
