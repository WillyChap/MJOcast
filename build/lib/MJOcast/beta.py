import MJOcast.utils.ProcessForecasts as ProFo 
import MJOcast.utils.ProcessOBS as ProObs
import importlib
importlib.reload(ProFo)
importlib.reload(ProObs)


# Define the path to the YAML file
yaml_file_path = './settings.yaml'
MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
#TODO... it feels like I can remove ""
OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()
MJO_obs.plot_phase_space('2001-01-01',60)

MJO_for = ProFo.MJOforecaster(yaml_file_path,MJO_obs.eof_dict,MJO_obs.MJO_fobs)
DS_CESM_for,OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd = MJO_for.create_forecasts(num_files=1)
MJO_for.plot_phase_space(12,15)
