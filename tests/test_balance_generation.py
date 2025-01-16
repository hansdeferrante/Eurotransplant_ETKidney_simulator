import sys
import unittest
import os

sys.path.append('./')

import simulator.magic_values.column_names as cn
from simulator.code.BalanceSystem import BalanceSystem
import simulator.code.utils.read_input_files as rdr
import simulator.magic_values.etkidney_simulator_settings as es
import os
from typing import Set, Optional


class TestBalance(unittest.TestCase):
    """This tests whether the blood group rules are correctly applied
    """

    def balance_test_for_gv(self, group_vars: Optional[Set[str]] = None):
        """Test whether ML have acceptable donor/recipient combinations.
        """

        ss = rdr.read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )
        d_bal = rdr.read_historic_donor_balances(
            ss.PATH_BALANCES,
            sim_start_date=ss.SIM_START_DATE,
            max_window_length=ss.WINDOW_FOR_BALANCE
        )

        d_bal.loc[:, 'value_nl'] = (
            -1 * (d_bal.donor_alloc_country == 'Netherlands') +
            1 * (d_bal.recipient_country == 'Netherlands')
        )

        # Add tstart and tstop columns for existing balances
        balances = BalanceSystem.from_balance_df(
            ss=ss,
            df_init_balances=d_bal,
            group_vars=group_vars
        )

        if group_vars is None:
            group_vars = set(('group_var',))
            man_counts = d_bal['value_nl'].sum()
            assert (
                man_counts ==
                balances.return_national_balances(
                    group_values=('1',)

                )['Netherlands']
            )
        else:
            man_counts = d_bal.groupby(
                list(group_vars)
            )['value_nl'].sum().to_dict()
            for group, counts in man_counts.items():
                assert counts == balances.return_national_balances(
                    group_values=(group,)
                )['Netherlands']

    def test_balance(self):
        self.balance_test_for_gv(group_vars=None)
        self.balance_test_for_gv(group_vars=set((cn.D_BLOODGROUP, )))


if __name__ == '__main__':
    unittest.main()
