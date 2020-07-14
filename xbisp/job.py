# pipeline submit
import os, numpy as np
from quicksub import *

def set_pfile(tag,btype='ggk',zn=30,zm=1,ngal=20.,expname='PLK'):
  f = 'sub_'+tag+'.sh'
  p = 'sub_p_'+tag+'.ini'
  set_sbatch_params(f,tag,mem='8G',t='1-00:00')
  add('./exe '+p,f)
  add('mv snr_'+btype+'_zn'+str(zn)+'.dat snr_'+expname+'_'+btype+'_zn'+str(zn)+'_zm'+str(zm)+'_ngal'+str(ngal)+'.dat',f)
  add('Om = 0.31201550908',p,True)
  add('H0 = 67.51',p)
  add('w0 = -1',p)
  add('wa = 0',p)
  add('mnu = 0.06',p)
  add('alpha = 2',p)
  add('beta  = 1',p)
  add('zm    = '+str(zm),p)
  add('ngal  = '+str(ngal),p)
  add('imax  = 15',p)
  add('bisptype = '+btype,p)
  add('pklin = ../../../data/forecast/pk/pklin.dat',p)
  if expname=='S4':  add('nldd = ../../../data/forecast/nldd/S4_s1_t3_rlmax4000.dat',p)
  if expname=='PLK': add('nldd = ../../../data/forecast/nldd/plk_s40_t7_quad_rlmax4000.dat',p)
  if expname=='PB':  add('nldd = ../../../data/forecast/nldd/PB_s7.22_t3.5_quad_rlmax4000.dat',p)
  add('clfile = ../../../data/forecast/cl/ucls.dat',p)
  add('zspace = 1',p)
  add('zs     = 1088.69',p)
  add('zran   = 0.001, 40',p)
  add('zn     = '+str(zn),p)
  add('model  = GM12',p)
  os.system('sbatch '+f)

for btype in ['gkk','ggk']:
  for zn in [30,40,50]:
    for zm, ngal in [(.5,10.),(1.,20.),(1.5,50.)]:
      set_pfile(btype+'_S4_zn'+str(zn)+'_zm'+str(zm)+'_ngal'+str(ngal),btype,zn,zm,ngal,'S4')
      set_pfile(btype+'_PLK_zn'+str(zn)+'_zm'+str(zm)+'_ngal'+str(ngal),btype,zn,zm,ngal,'PLK')
      set_pfile(btype+'_PB_zn'+str(zn)+'_zm'+str(zm)+'_ngal'+str(ngal),btype,zn,zm,ngal,'PB')

