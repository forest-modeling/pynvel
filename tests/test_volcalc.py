'''
Created on Mar 5, 2016

@author: THAREN
'''
import unittest

import numpy as np
import pandas as pd

import pynvel

class Test(unittest.TestCase):


  def setUp(self):
    pass


  def tearDown(self):
    pass

  def test_products(self):
    mrule = {
      'evod':2, 'opt':23
      , 'maxlen':40.0, 'minlen':12.0
      , 'minlent':12.0, 'merchl': 20.0
      , 'mtopp':5.0, 'mtops':2.0, 'stump':1.0, 'trim':1.0
      , 'btr':0.0, 'dbtbh':0.0, 'minbfd':8.0, 'cor':'Y'
      }
    
    vol_eq = 'F01FW2W202'
    # vol_eq = 'NVBM240202'
    dbh_ob = 20
    total_ht = 130

    # 40 foot maximum log lenght
    mrule['maxlen'] = 40
    mrule['mtopp'] = 5.0
    _mrule = pynvel.init_merchrule(**mrule)

    vc = pynvel.VolumeCalculator(
      volume_eq=vol_eq, 
      merch_rule=_mrule,
      calc_products=True
      )
    e = vc.calc(dbh_ob=dbh_ob, total_ht=total_ht)

    print('* Products *')
    # print(pd.DataFrame([p.as_dict() for p in vc.products]))
    print(pd.DataFrame.from_dict(vc.products, orient='index').round(1))

  def test_loglen(self):
    """Ensure that changing mrule%maxlen has the expected effect
    """
    mrule = {
      'evod':2, 'opt':23
      , 'maxlen':40.0, 'minlen':12.0
      , 'minlent':12.0, 'merchl': 20.0
      , 'mtopp':5.0, 'mtops':2.0, 'stump':1.0, 'trim':1.0
      , 'btr':0.0, 'dbtbh':0.0, 'minbfd':8.0, 'cor':'Y'
      }
    
    vol_eq = 'F01FW2W202'
    # vol_eq = 'NVBM240202'
    dbh_ob = 20
    total_ht = 130

    # 40 foot maximum log lenght
    mrule['maxlen'] = 40
    mrule['mtopp'] = 5.0
    _mrule = pynvel.init_merchrule(**mrule)

    vc_40 = pynvel.VolumeCalculator(volume_eq=vol_eq, merch_rule=_mrule)
    e = vc_40.calc(dbh_ob=dbh_ob, total_ht=total_ht)

    bdft_40 = vc_40.volume['bdft_gross_prim']
    cuft_40 = vc_40.volume['cuft_total']
    logs_40 = pd.DataFrame.from_records([l.as_dict() for l in vc_40.logs]).round(1)

    # print(vc.volume)

    # 16 foot maximum log lenght
    mrule['maxlen'] = 16
    mrule['mtopp'] = 5.0
    _mrule = pynvel.init_merchrule(**mrule)

    vc_16 = pynvel.VolumeCalculator(volume_eq=vol_eq, merch_rule=_mrule)
    e = vc_16.calc(dbh_ob=dbh_ob, total_ht=total_ht)

    bdft_16 = vc_16.volume['bdft_gross_prim']
    cuft_16 = vc_16.volume['cuft_total']
    logs_16 = pd.DataFrame.from_records([l.as_dict() for l in vc_16.logs]).round(1)

    # print(vc.volume)

    # print(cuft_40, bdft_40)
    # print(cuft_16, bdft_16)

    # print(logs_40)
    # print(logs_16)

    # The two should be the same for cubic volume
    assert np.isclose(cuft_40, cuft_16, atol=0.001)

    # 40 foot logs should report less Scribner volume
    assert np.less(bdft_40, bdft_16)

if __name__ == "__main__":
  # import sys;sys.argv = ['', 'Test.testName']
  unittest.main()
