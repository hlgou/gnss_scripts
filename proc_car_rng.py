import logging
from proc_upd import ProcUpd
from funcs import check_res_sigma, GrtPpplsq


class ProcCarRng(ProcUpd):

    def init_daily(self):
        self._config.carrier_range = False
        self._config.carrier_range_out = True
        return super().init_daily()

    def process_ppp(self):
        logging.info(f"===> Calculate float ambiguities by precise point positioning")
        GrtPpplsq(self._config, 'ppplsq', nmp=self.nthread).run()
        self.basic_check(files=['recover_all', 'ambupd_in'])

        self.editres(bad=80, jump=80, nshort=600, all_sites=True)
        GrtPpplsq(self._config, 'ppplsq', nmp=self.nthread).run()
        self.basic_check(files=['recover_all', 'ambupd_in'])
        check_res_sigma(self._config)
        self.editres(bad=40, jump=40, nshort=600, all_sites=True)

    def ppp_clean(self):
        logging.info(f"===> Detect outliers in carrier-range by PPP")
        self._config.crd_constr = 'FIX'
        self._config.carrier_range = True
        GrtPpplsq(self._config, 'ppplsq', nmp=self.nthread).run()
        self.editres(jump=40, edt_amb=True, all_sites=True)

    def process_daily(self):
        self.process_upd(obs_comb='UC')
        # backup_dir('ambupd', 'ambupd_save')
        # self.ppp_clean()
        # logging.info(f"===> Generate RINEX Carrier-range observation")


if __name__ == '__main__':
    proc = ProcCarRng.from_args()
    proc.process_batch()
