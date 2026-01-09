import unittest

import numpy as np
import pandas as pd

import pynvel

class Test(unittest.TestCase):

  def setUp(self):
    """
    Test setup, load datasets
    """

    pynvel.config = pynvel.get_config()

    # TODO: Test multiple equations using nose-parameterized
    self.vol_eq = 'F01FW2W202'
    self.log_len = 16

    # # Default region 6 merchandizing; mrules.f:153
    # self.mrule = pynvel.init_merchrule(
    #     evod=2, opt=23,
    #     maxlen=self.log_len, minlen=2.0, minlent=2.0,
    #     merchl=8.0, mtopp=5, mtops=2, trim=0.5, stump=0.0,
    #     cor='Y', minbfd=8
    #     )

  def test_df(self):
    
    vc = pynvel.VolumeCalculator(
        region=6, forest='12'
        , volume_eq=self.vol_eq
        # , merch_rule=self.mrule
        , calc_products=True)
    
    # # Verify the product classification methods work on arrays
    vol = vc.calc_array(
        np.array([18.0,24.0]),
        np.array([120.0,150.0])
        )
    df = pd.DataFrame(vol, columns=vc.volume_labels)
    print(df)
    
if __name__ == "__main__":
    unittest.main()