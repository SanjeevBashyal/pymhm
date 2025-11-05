# -*- coding: utf-8 -*-
"""
Simulation processing module for pymhm
Contains all simulation-related processing methods including NML file generation and mHM execution
"""
import os
import subprocess
from datetime import datetime
from qgis.PyQt.QtWidgets import QMessageBox
from qgis.PyQt.QtCore import QDate
from .utils import DialogUtils


class SimulationProcessor(DialogUtils):
    """
    Handles all simulation processing functionality.
    Expects the dialog instance to have DialogUtils attributes.
    """
    
    def __init__(self, dialog):
        """
        Initialize the simulation processor with reference to the dialog.
        
        Args:
            dialog: The pymhmDialog instance (needs DialogUtils methods)
        """
        self.dialog = dialog
        # Copy DialogUtils methods to self
        self.log_message = dialog.log_message
        self.check_prerequisites = dialog.check_prerequisites
        
    def get_process_case_value(self, combo_box):
        """
        Extract numeric value from combo box text.
        Combo box text format: "1 - Description" or "-1 - Description"
        """
        current_text = combo_box.currentText()
        if current_text:
            # Extract the number before the first dash
            parts = current_text.split(' - ', 1)
            if parts:
                try:
                    return int(parts[0].strip())
                except ValueError:
                    return 1
        return 1
    
    def get_boolean_value(self, checkbox):
        """Get boolean value from checkbox."""
        return checkbox.isChecked()
    
    def get_date_string(self, date_edit):
        """Get date string in YYYY-MM-DD format from QDateEdit."""
        date = date_edit.date()
        return date.toString("yyyy-MM-dd")
    
    def read_ui_inputs(self):
        """
        Read all inputs from the simulation tab UI elements.
        
        Returns:
            dict: Dictionary containing all simulation parameters
        """
        d = self.dialog
        
        # Process Selection
        inputs = {
            'processCase': {
                1: self.get_process_case_value(d.comboBox_interception),
                2: self.get_process_case_value(d.comboBox_snow),
                3: self.get_process_case_value(d.comboBox_soilmoisture),
                4: self.get_process_case_value(d.comboBox_directrunoff),
                5: self.get_process_case_value(d.comboBox_pet),
                6: self.get_process_case_value(d.comboBox_interflow),
                7: self.get_process_case_value(d.comboBox_percolation),
                8: self.get_process_case_value(d.comboBox_routing),
                9: self.get_process_case_value(d.comboBox_baseflow),
                10: self.get_process_case_value(d.comboBox_neutrons),
                11: self.get_process_case_value(d.comboBox_rivertemp),
            },
            # Time Settings
            'timestep': int(d.comboBox_timestep.currentText()),
            'warming_days': d.spinBox_warming_days.value(),
            'sim_start': self.get_date_string(d.dateEdit_sim_start),
            'sim_end': self.get_date_string(d.dateEdit_sim_end),
            # Resolution Settings
            'resolution_hydrology': d.spinBox_resolution_hydrology.value(),
            'resolution_routing': d.spinBox_resolution_routing.value(),
            'nDomains': d.spinBox_nDomains.value(),
            # Soil and LAI Settings
            'nSoilHorizons_mHM': d.spinBox_soil_horizons.value(),
            'tillageDepth': d.spinBox_tillage_depth.value(),
            'lai_timestep': self.get_process_case_value(d.comboBox_lai_timestep),
            'lai_format': d.comboBox_lai_format.currentText(),
            # Parameters
            'canopyInterceptionFactor': d.doubleSpinBox_canopy_interception.value(),
            'snowTreshholdTemperature': d.doubleSpinBox_snow_threshold_temp.value(),
            'degreeDayFactor_forest': d.doubleSpinBox_degday_forest.value(),
            'degreeDayFactor_impervious': d.doubleSpinBox_degday_impervious.value(),
            'degreeDayFactor_pervious': d.doubleSpinBox_degday_previous.value(),
            'pet_min_corr': d.doubleSpinBox_pet_min_corr.value(),
            'pet_max_corr': d.doubleSpinBox_pet_max_corr.value(),
            'pet_aspect_threshold': d.doubleSpinBox_pet_aspect_threshold.value(),
            'hargreaves_coeff': d.doubleSpinBox_hargreaves_coeff.value(),
            'interflow_storage_capacity': d.doubleSpinBox_interflow_storage_capacity.value(),
            'recharge_coeff': d.doubleSpinBox_recharge_coeff.value(),
            'muskingum_traveltime_constant': d.doubleSpinBox_muskingum_traveltime_constant.value(),
            'streamflow_celerity': d.doubleSpinBox_streamflow_celerity.value(),
            'slope_factor': d.doubleSpinBox_slope_factor.value(),
            'impervious_storage_capacity': d.doubleSpinBox_impervious_storage_capacity.value(),
            # Output Settings
            'output_timestep_mode': d.comboBox_output_timestep.currentText(),
            'output_timestep_value': d.spinBox_output_timestep_value.value(),
            'output_time_reference': self.get_process_case_value(d.comboBox_output_time_reference),
            'output_deflate_level': d.spinBox_output_deflate.value(),
            'output_double_precision': d.checkBox_output_double_precision.isChecked(),
            # Output States
            'output_intercept': d.checkBox_output_intercept.isChecked(),
            'output_snowpack': d.checkBox_output_snowpack.isChecked(),
            'output_soilmoist': d.checkBox_output_soilmoist.isChecked(),
            'output_volsoilmoist': d.checkBox_output_volsoilmoist.isChecked(),
            'output_meanvolsoilmoist': d.checkBox_output_meanvolsoilmoist.isChecked(),
            'output_sealstw': d.checkBox_output_sealstw.isChecked(),
            'output_unsatstw': d.checkBox_output_unsatstw.isChecked(),
            'output_satstw': d.checkBox_output_satstw.isChecked(),
            'output_neutrons': d.checkBox_output_neutrons.isChecked(),
            # Output Fluxes
            'output_pet': d.checkBox_output_pet.isChecked(),
            'output_aet': d.checkBox_output_aet.isChecked(),
            'output_total_runoff': d.checkBox_output_total_runoff.isChecked(),
            'output_direct_runoff': d.checkBox_output_direct_runoff.isChecked(),
            'output_fast_interflow': d.checkBox_output_fast_interflow.isChecked(),
            'output_slow_interflow': d.checkBox_output_slow_interflow.isChecked(),
            'output_baseflow': d.checkBox_output_baseflow.isChecked(),
            'output_recharge': d.checkBox_output_recharge.isChecked(),
            'output_infiltration': d.checkBox_output_infiltration.isChecked(),
            'output_aet_soil': d.checkBox_output_aet_soil.isChecked(),
            'output_effective_prec': d.checkBox_output_effective_prec.isChecked(),
            'output_snowmelt': d.checkBox_output_snowmelt.isChecked(),
            'output_routed_streamflow': d.checkBox_output_routed_streamflow.isChecked(),
            'output_routed_temp': d.checkBox_output_routed_temp.isChecked(),
        }
        
        return inputs
    
    def write_mhm_nml(self, inputs, output_dir):
        """
        Write mhm.nml file.
        
        Args:
            inputs: Dictionary of input parameters
            output_dir: Directory where NML files should be written
        """
        nml_path = os.path.join(output_dir, "mhm.nml")
        
        # Parse dates
        start_date = datetime.strptime(inputs['sim_start'], "%Y-%m-%d")
        end_date = datetime.strptime(inputs['sim_end'], "%Y-%m-%d")
        
        with open(nml_path, 'w') as f:
            f.write("&project_description\n")
            f.write('project_details="mHM project"\n')
            f.write('setup_description="model run configured via pymhm"\n')
            f.write('simulation_type="historical simulation"\n')
            f.write('Conventions="XXX"\n')
            f.write('contact="pymhm user"\n')
            f.write('mHM_details="Helmholtz Center for Environmental Research - UFZ, Department Computational Hydrosystems, Stochastic Hydrology Group"\n')
            f.write('history="model run version 1"\n')
            f.write("/\n")
            f.write("\n")
            
            f.write("&mainconfig\n")
            f.write("iFlag_cordinate_sys = 0\n")
            f.write(f"nDomains = {inputs['nDomains']}\n")
            f.write(f"resolution_Hydrology(1) = {inputs['resolution_hydrology']}\n")
            f.write("L0Domain(1) = 1\n")
            f.write("write_restart = .FALSE.\n")
            f.write("read_opt_domain_data(1) = 0\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&mainconfig_mhm_mrm\n")
            f.write('mhm_file_RestartIn(1) = "restart/mHM_restart_001.nc"\n')
            f.write('mrm_file_RestartIn(1) = "restart/mRM_restart_001.nc"\n')
            f.write(f"resolution_Routing(1) = {inputs['resolution_routing']}\n")
            f.write(f"timestep = {inputs['timestep']}\n")
            f.write("read_restart = .FALSE.\n")
            f.write("optimize = .FALSE.\n")
            f.write("optimize_restart = .FALSE.\n")
            f.write("opti_method = 1\n")
            f.write("opti_function = 10\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&mainconfig_mrm\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&directories_general\n")
            f.write('dirConfigOut = "./"\n')
            f.write('dirCommonFiles = "input/morph/"\n')
            f.write('dir_Morpho(1) = "input/morph/"\n')
            f.write('dir_LCover(1) = "input/luse/"\n')
            f.write('mhm_file_RestartOut(1) = "restart/mHM_restart_001.nc"\n')
            f.write('mrm_file_RestartOut(1) = "restart/mRM_restart_001.nc"\n')
            f.write('dir_Out(1) = "output_b1/"\n')
            f.write('file_LatLon(1) = "input/latlon/latlon_1.nc"\n')
            f.write("/\n")
            f.write("\n")
            
            f.write("&directories_mHM\n")
            f.write('inputFormat_meteo_forcings = "nc"\n')
            f.write('dir_Precipitation(1) = "input/meteo/pre/"\n')
            f.write('dir_Temperature(1) = "input/meteo/tavg/"\n')
            f.write('dir_ReferenceET(1) = "input/meteo/pet/"\n')
            f.write("time_step_model_inputs(1) = 0\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&directories_mRM\n")
            f.write('dir_Gauges(1) = "input/gauge/"\n')
            f.write('dir_Total_Runoff(1) = "output_b1/"\n')
            f.write('dir_Bankfull_Runoff(1) = "input/optional_data/"\n')
            f.write("/\n")
            f.write("\n")
            
            f.write("&processSelection\n")
            pc = inputs['processCase']
            f.write(f"processCase(1) = {pc[1]}\n")
            f.write(f"processCase(2) = {pc[2]}\n")
            f.write(f"processCase(3) = {pc[3]}\n")
            f.write(f"processCase(4) = {pc[4]}\n")
            f.write(f"processCase(6) = {pc[6]}\n")
            f.write(f"processCase(7) = {pc[7]}\n")
            f.write(f"processCase(8) = {pc[8]}\n")
            f.write(f"processCase(9) = {pc[9]}\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&LCover\n")
            f.write("nLCoverScene = 1\n")
            f.write(f"LCoverYearStart(1) = {start_date.year}\n")
            f.write(f"LCoverYearEnd(1) = {end_date.year}\n")
            f.write('LCoverfName(1) = "lc_1981.asc"\n')
            f.write("/\n")
            f.write("\n")
            
            f.write("&time_periods\n")
            f.write(f"warming_Days(1) = {inputs['warming_days']}\n")
            f.write(f"eval_Per(1)%yStart = {start_date.year}\n")
            f.write(f"eval_Per(1)%mStart = {start_date.month:02d}\n")
            f.write(f"eval_Per(1)%dStart = {start_date.day:02d}\n")
            f.write(f"eval_Per(1)%yEnd = {end_date.year}\n")
            f.write(f"eval_Per(1)%mEnd = {end_date.month:02d}\n")
            f.write(f"eval_Per(1)%dEnd = {end_date.day:02d}\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&soildata\n")
            f.write("iFlag_soilDB = 0\n")
            f.write(f"tillageDepth = {inputs['tillageDepth']}\n")
            f.write(f"nSoilHorizons_mHM = {inputs['nSoilHorizons_mHM']}\n")
            if inputs['nSoilHorizons_mHM'] > 1:
                f.write("soil_Depth(1) = 200\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&LAI_data_information\n")
            f.write(f"timeStep_LAI_input = {inputs['lai_timestep']}\n")
            f.write(f'inputFormat_gridded_LAI = "{inputs["lai_format"]}"\n')
            f.write("/\n")
            f.write("\n")
            
            f.write("&LCover_MPR\n")
            f.write("fracSealed_cityArea = 0.6\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&evaluation_gauges\n")
            f.write("nGaugesTotal = 1\n")
            f.write("NoGauges_domain(1) = 1\n")
            f.write("Gauge_id(1,1) = 398\n")
            f.write('gauge_filename(1,1) = "00398.txt"\n')
            f.write("/\n")
            f.write("\n")
            
            f.write("&inflow_gauges\n")
            f.write("InflowGauge_Headwater(1,1) = .FALSE.\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&panEvapo\n")
            f.write("evap_coeff = 1.30, 1.20, 0.72, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.50\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&nightDayRatio\n")
            f.write("read_meteo_weights = .FALSE.\n")
            f.write("fnight_prec = 0.46, 0.50, 0.52, 0.51, 0.48, 0.50, 0.49, 0.48, 0.52, 0.56, 0.50, 0.47\n")
            f.write("fnight_pet = 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10\n")
            f.write("fnight_temp = -0.76, -1.30, -1.88, -2.38, -2.72, -2.75, -2.74, -3.04, -2.44, -1.60, -0.94, -0.53\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&Optimization\n")
            f.write("nIterations = 7\n")
            f.write("seed = 1235876\n")
            f.write("sce_ngs = 2\n")
            f.write("/\n")
        
        self.log_message(f"Written: {nml_path}")
    
    def write_mhm_parameter_nml(self, inputs, output_dir):
        """
        Write mhm_parameter.nml file.
        
        Args:
            inputs: Dictionary of input parameters
            output_dir: Directory where NML files should be written
        """
        nml_path = os.path.join(output_dir, "mhm_parameter.nml")
        
        with open(nml_path, 'w') as f:
            f.write("&interception1\n")
            f.write(f"canopyInterceptionFactor = 0.1500, 0.4000, {inputs['canopyInterceptionFactor']:.4f}, 1, 1\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&snow1\n")
            f.write(f"snowTreshholdTemperature = -2.0000, 2.0000, {inputs['snowTreshholdTemperature']:.1f}, 1, 1\n")
            f.write(f"degreeDayFactor_forest = 0.0001, 4.0000, {inputs['degreeDayFactor_forest']:.4f}, 1, 1\n")
            f.write(f"degreeDayFactor_impervious = 0.0000, 1.0000, {inputs['degreeDayFactor_impervious']:.4f}, 1, 1\n")
            f.write(f"degreeDayFactor_pervious = 0.0000, 2.0000, {inputs['degreeDayFactor_pervious']:.4f}, 1, 1\n")
            f.write("increaseDegreeDayFactorByPrecip = 0.1000, 0.9000, 0.5, 1, 1\n")
            f.write("maxDegreeDayFactor_forest = 0.0000, 8.0000, 3.0, 1, 1\n")
            f.write("maxDegreeDayFactor_impervious = 0.0000, 8.0000, 3.5, 1, 1\n")
            f.write("maxDegreeDayFactor_pervious = 0.0000, 8.0000, 4.0, 1, 1\n")
            f.write("/\n")
            f.write("\n")
            
            # Soil moisture parameters (simplified - using case 1)
            f.write("&soilmoisture1\n")
            f.write("orgMatterContent_forest = 0.0000, 20.000, 3.4, 1, 1\n")
            f.write("orgMatterContent_impervious = 0.0000, 1.0000, 0.1, 1, 1\n")
            f.write("orgMatterContent_pervious = 0.0000, 4.0000, 0.6, 1, 1\n")
            f.write("PTF_lower66_5_constant = 0.6462, 0.9506, 0.76, 1, 1\n")
            f.write("PTF_lower66_5_clay = 0.0001, 0.0029, 0.0009, 1, 1\n")
            f.write("PTF_lower66_5_Db = -0.3727, -0.1871, -0.264, 1, 1\n")
            f.write("PTF_higher66_5_constant = 0.5358, 1.1232, 0.89, 1, 1\n")
            f.write("PTF_higher66_5_clay = -0.0055, 0.0049, -0.001, 1, 1\n")
            f.write("PTF_higher66_5_Db = -0.5513, -0.0913, -0.324, 1, 1\n")
            f.write("PTF_Ks_constant = -1.2000, -0.2850, -0.585, 1, 1\n")
            f.write("PTF_Ks_sand = 0.0060, 0.0260, 0.0125, 1, 1\n")
            f.write("PTF_Ks_clay = 0.0030, 0.0130, 0.0063, 1, 1\n")
            f.write("PTF_Ks_curveSlope = 60.960, 60.960, 60.960, 0, 1\n")
            f.write("rootFractionCoefficient_forest = 0.9000, 0.9990, 0.97, 1, 1\n")
            f.write("rootFractionCoefficient_impervious = 0.9000, 0.9500, 0.93, 1, 1\n")
            f.write("rootFractionCoefficient_pervious = 0.0010, 0.0900, 0.02, 1, 1\n")
            f.write("infiltrationShapeFactor = 1.0000, 4.0000, 1.75, 1, 1\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&directRunoff1\n")
            f.write(f"imperviousStorageCapacity = 0.0000, 5.0000, {inputs['impervious_storage_capacity']:.2f}, 1, 1\n")
            f.write("/\n")
            f.write("\n")
            
            # PET parameters based on processCase(5)
            pet_case = inputs['processCase'][5]
            if pet_case == -1:
                f.write("&PETminus1\n")
                f.write("PET_a_forest = 0.3000, 1.3000, 0.3000, 1, 1\n")
                f.write("PET_a_impervious = 0.3000, 1.3000, 0.8000, 1, 1\n")
                f.write("PET_a_pervious = 0.3000, 1.3000, 1.3000, 1, 1\n")
                f.write("PET_b = 0.0000, 1.5000, 1.5000, 1, 1\n")
                f.write("PET_c = -2.000, 0.0000, -0.700, 1, 1\n")
                f.write("/\n")
                f.write("\n")
            elif pet_case == 0:
                f.write("&PET0\n")
                f.write(f"minCorrectionFactorPET = 0.7000, 1.3000, {inputs['pet_min_corr']:.4f}, 1, 1\n")
                f.write(f"maxCorrectionFactorPET = 0.0000, 0.2000, {inputs['pet_max_corr']:.4f}, 1, 1\n")
                f.write(f"aspectTresholdPET = 160.00, 200.00, {inputs['pet_aspect_threshold']:.1f}, 1, 1\n")
                f.write("/\n")
                f.write("\n")
            elif pet_case == 1:
                f.write("&PET1\n")
                f.write(f"minCorrectionFactorPET = 0.7000, 1.3000, {inputs['pet_min_corr']:.4f}, 1, 1\n")
                f.write(f"maxCorrectionFactorPET = 0.0000, 0.2000, {inputs['pet_max_corr']:.4f}, 1, 1\n")
                f.write(f"aspectTresholdPET = 160.00, 200.00, {inputs['pet_aspect_threshold']:.1f}, 1, 1\n")
                f.write(f"HargreavesSamaniCoeff = 0.0016, 0.0030, {inputs['hargreaves_coeff']:.4f}, 1, 1\n")
                f.write("/\n")
                f.write("\n")
            elif pet_case == 2:
                f.write("&PET2\n")
                f.write("PriestleyTaylorCoeff = 0.75, 1.75, 1.1900, 1, 1\n")
                f.write("PriestleyTaylorLAIcorr = -0.50, 0.20, 0.0580, 1, 1\n")
                f.write("/\n")
                f.write("\n")
            elif pet_case == 3:
                f.write("&PET3\n")
                f.write("canopyheigth_forest = 15.00, 40.00, 15.000, 1, 1\n")
                f.write("canopyheigth_impervious = 0.01, 0.50, 0.0200, 1, 1\n")
                f.write("canopyheigth_pervious = 0.10, 5.00, 0.1100, 1, 1\n")
                f.write("displacementheight_coeff = 0.50, 0.85, 0.6400, 1, 1\n")
                f.write("roughnesslength_momentum_coeff = 0.09, 0.16, 0.0950, 1, 1\n")
                f.write("roughnesslength_heat_coeff = 0.07, 0.13, 0.0750, 1, 1\n")
                f.write("stomatal_resistance = 10.00, 200.00, 56.000, 1, 1\n")
                f.write("/\n")
                f.write("\n")
            
            f.write("&interflow1\n")
            f.write(f"interflowStorageCapacityFactor = 75.000, 200.00, {inputs['interflow_storage_capacity']:.2f}, 1, 1\n")
            f.write("interflowRecession_slope = 0.0000, 10.000, 7.0, 1, 1\n")
            f.write("fastInterflowRecession_forest = 1.0000, 3.0000, 1.5, 1, 1\n")
            f.write("slowInterflowRecession_Ks = 1.0000, 30.000, 15.0, 1, 1\n")
            f.write("exponentSlowInterflow = 0.0500, 0.3000, 0.125, 1, 1\n")
            f.write("/\n")
            f.write("\n")
            
            f.write("&percolation1\n")
            f.write(f"rechargeCoefficient = 0.0000, 50.000, {inputs['recharge_coeff']:.2f}, 1, 1\n")
            f.write("rechargeFactor_karstic = -5.0000, 5.0000, -1.0, 1, 1\n")
            f.write("gain_loss_GWreservoir_karstic = 1.0000, 1.0000, 1.0, 0, 1\n")
            f.write("/\n")
            f.write("\n")
            
            # Routing parameters based on processCase(8)
            routing_case = inputs['processCase'][8]
            if routing_case == 1:
                f.write("&routing1\n")
                f.write(f"muskingumTravelTime_constant = 0.3100, 0.3500, {inputs['muskingum_traveltime_constant']:.4f}, 1, 1\n")
                f.write("muskingumTravelTime_riverLength = 0.0700, 0.0800, 0.075, 1, 1\n")
                f.write("muskingumTravelTime_riverSlope = 1.9500, 2.1000, 2.0, 1, 1\n")
                f.write("muskingumTravelTime_impervious = 0.0900, 0.1100, 0.1, 1, 1\n")
                f.write("muskingumAttenuation_riverSlope = 0.0100, 0.5000, 0.3, 1, 1\n")
                f.write("/\n")
                f.write("\n")
            elif routing_case == 2:
                f.write("&routing2\n")
                f.write(f"streamflow_celerity = 0.1, 15., {inputs['streamflow_celerity']:.2f}, 0, 1\n")
                f.write("/\n")
                f.write("\n")
            elif routing_case == 3:
                f.write("&routing3\n")
                f.write(f"slope_factor = 0.1, 100., {inputs['slope_factor']:.2f}, 0, 1\n")
                f.write("/\n")
                f.write("\n")
            
            f.write("&geoparameter\n")
            f.write("GeoParam(1,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("GeoParam(2,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("GeoParam(3,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("GeoParam(4,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("GeoParam(5,:) = 1.000, 1000.00, 100.0, 0, 1\n")
            f.write("GeoParam(6,:) = 1.000, 1000.00, 100.0, 0, 1\n")
            f.write("GeoParam(7,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("GeoParam(8,:) = 1.000, 1000.00, 100.0, 0, 1\n")
            f.write("GeoParam(9,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("GeoParam(10,:) = 1.000, 1000.00, 100.0, 1, 1\n")
            f.write("/\n")
        
        self.log_message(f"Written: {nml_path}")
    
    def write_mhm_outputs_nml(self, inputs, output_dir):
        """
        Write mhm_outputs.nml file.
        
        Args:
            inputs: Dictionary of input parameters
            output_dir: Directory where NML files should be written
        """
        nml_path = os.path.join(output_dir, "mhm_outputs.nml")
        
        # Determine output timestep value
        output_timestep_value = inputs['output_timestep_value']
        
        with open(nml_path, 'w') as f:
            f.write("&NLoutputResults\n")
            f.write(f"output_deflate_level = {inputs['output_deflate_level']}\n")
            f.write(f"output_double_precision = {'.TRUE.' if inputs['output_double_precision'] else '.FALSE.'}\n")
            f.write(f"output_time_reference = {inputs['output_time_reference']}\n")
            f.write(f"timeStep_model_outputs = {output_timestep_value}\n")
            f.write(f"outputFlxState(1) = {'.TRUE.' if inputs['output_intercept'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(2) = {'.TRUE.' if inputs['output_snowpack'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(3) = {'.TRUE.' if inputs['output_soilmoist'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(4) = {'.TRUE.' if inputs['output_volsoilmoist'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(5) = {'.TRUE.' if inputs['output_meanvolsoilmoist'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(6) = {'.TRUE.' if inputs['output_sealstw'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(7) = {'.TRUE.' if inputs['output_unsatstw'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(8) = {'.TRUE.' if inputs['output_satstw'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(9) = {'.TRUE.' if inputs['output_pet'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(10) = {'.TRUE.' if inputs['output_aet'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(11) = {'.TRUE.' if inputs['output_total_runoff'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(12) = {'.TRUE.' if inputs['output_direct_runoff'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(13) = {'.TRUE.' if inputs['output_fast_interflow'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(14) = {'.TRUE.' if inputs['output_slow_interflow'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(15) = {'.TRUE.' if inputs['output_baseflow'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(16) = {'.TRUE.' if inputs['output_recharge'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(17) = {'.TRUE.' if inputs['output_infiltration'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(18) = {'.TRUE.' if inputs['output_neutrons'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(19) = {'.TRUE.' if inputs['output_aet_soil'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(20) = {'.TRUE.' if inputs['output_effective_prec'] else '.FALSE.'}\n")
            f.write(f"outputFlxState(21) = {'.TRUE.' if inputs['output_snowmelt'] else '.FALSE.'}\n")
            f.write("/\n")
        
        self.log_message(f"Written: {nml_path}")
    
    def write_mrm_outputs_nml(self, inputs, output_dir):
        """
        Write mrm_outputs.nml file.
        
        Args:
            inputs: Dictionary of input parameters
            output_dir: Directory where NML files should be written
        """
        nml_path = os.path.join(output_dir, "mrm_outputs.nml")
        
        with open(nml_path, 'w') as f:
            f.write("&NLoutputResults\n")
            f.write(f"output_deflate_level_mrm = {inputs['output_deflate_level']}\n")
            f.write(f"output_double_precision_mrm = {'.true.' if inputs['output_double_precision'] else '.false.'}\n")
            f.write(f"output_time_reference_mrm = {inputs['output_time_reference']}\n")
            f.write("timeStep_model_outputs_mrm = -1\n")
            f.write(f"outputFlxState_mrm(1) = {'.TRUE.' if inputs['output_routed_streamflow'] else '.FALSE.'}\n")
            f.write(f"outputFlxState_mrm(2) = {'.TRUE.' if inputs['output_routed_temp'] else '.FALSE.'}\n")
            f.write("/\n")
        
        self.log_message(f"Written: {nml_path}")
    
    def create_nml_files(self):
        """
        Read UI inputs and create all NML files.
        This function is called by the createNML pushbutton.
        """
        if not self.check_prerequisites():
            return False
        
        self.log_message("\n--- Creating NML files ---")
        
        # Determine output directory
        model_inputs_dir = os.path.join(self.dialog.project_folder, "model_inputs")
        if not os.path.exists(model_inputs_dir):
            os.makedirs(model_inputs_dir)
            self.log_message(f"Created directory: {model_inputs_dir}")
        
        try:
            # Read all UI inputs
            inputs = self.read_ui_inputs()
            
            # Write all NML files
            self.write_mhm_nml(inputs, model_inputs_dir)
            self.write_mhm_parameter_nml(inputs, model_inputs_dir)
            self.write_mhm_outputs_nml(inputs, model_inputs_dir)
            self.write_mrm_outputs_nml(inputs, model_inputs_dir)
            
            self.log_message("NML files created successfully!")
            return True
            
        except Exception as e:
            self.log_message(f"ERROR: Failed to create NML files. Details: {str(e)}")
            QMessageBox.critical(
                self.dialog, "Error", f"Failed to create NML files.\n{str(e)}")
            return False
    
    def run_mhm(self):
        """
        Run mHM command via conda.
        This function is called by the RUN pushbutton.
        First calls create_nml_files, then runs mHM.
        """
        if not self.check_prerequisites():
            return False
        
        self.log_message("\n--- Running mHM ---")
        
        # First create NML files
        self.log_message("Creating NML files...")
        if not self.create_nml_files():
            self.log_message("ERROR: Failed to create NML files. Cannot run mHM.")
            return False
        
        # Determine model inputs directory
        model_inputs_dir = os.path.join(self.dialog.project_folder, "model_inputs")
        
        if not os.path.exists(model_inputs_dir):
            self.log_message(f"ERROR: Model inputs directory not found: {model_inputs_dir}")
            QMessageBox.critical(
                self.dialog, "Error", 
                f"Model inputs directory not found: {model_inputs_dir}")
            return False
        
        try:
            # Change to model_inputs directory and run mHM via conda
            self.log_message(f"Running mHM in directory: {model_inputs_dir}")
            
            # Run mHM command via conda
            # Assuming mHM is installed in a conda environment
            # The command format: conda run -n <env_name> mHM <folder_name>
            # Or if mHM is in PATH: mHM model_inputs
            
            # Try to run mHM directly first (if it's in PATH)
            cmd = ["mHM", "model_inputs"]
            
            self.log_message(f"Executing command: {' '.join(cmd)}")
            
            # Run the command
            process = subprocess.Popen(
                cmd,
                cwd=self.dialog.project_folder,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
            )
            
            # Stream output in real-time
            self.log_message("mHM execution started...")
            output_lines = []
            for line in process.stdout:
                line = line.strip()
                if line:
                    self.log_message(line)
                    output_lines.append(line)
            
            # Wait for process to complete
            process.wait()
            
            if process.returncode == 0:
                self.log_message("mHM execution completed successfully!")
                QMessageBox.information(
                    self.dialog, "Success", 
                    "mHM execution completed successfully!")
                return True
            else:
                self.log_message(f"ERROR: mHM execution failed with return code {process.returncode}")
                QMessageBox.critical(
                    self.dialog, "Error", 
                    f"mHM execution failed.\nReturn code: {process.returncode}\nCheck log for details.")
                return False
                
        except FileNotFoundError:
            # If mHM is not in PATH, try via conda
            self.log_message("mHM not found in PATH. Trying conda...")
            try:
                # Try common conda environment names
                conda_envs = ["mhm", "mHM", "base", "conda"]
                
                for env_name in conda_envs:
                    cmd = ["conda", "run", "-n", env_name, "mHM", "model_inputs"]
                    self.log_message(f"Trying: {' '.join(cmd)}")
                    
                    process = subprocess.Popen(
                        cmd,
                        cwd=self.dialog.project_folder,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        bufsize=1,
                        universal_newlines=True
                    )
                    
                    self.log_message("mHM execution started...")
                    for line in process.stdout:
                        line = line.strip()
                        if line:
                            self.log_message(line)
                    
                    process.wait()
                    
                    if process.returncode == 0:
                        self.log_message("mHM execution completed successfully!")
                        QMessageBox.information(
                            self.dialog, "Success", 
                            "mHM execution completed successfully!")
                        return True
                    else:
                        continue
                
                # If all conda environments failed
                self.log_message("ERROR: Could not find mHM in any conda environment.")
                QMessageBox.critical(
                    self.dialog, "Error", 
                    "Could not find mHM executable.\n"
                    "Please ensure mHM is installed and in your PATH,\n"
                    "or configure the conda environment name.")
                return False
                
            except Exception as e:
                self.log_message(f"ERROR: Failed to run mHM via conda. Details: {str(e)}")
                QMessageBox.critical(
                    self.dialog, "Error", 
                    f"Failed to run mHM.\n{str(e)}\n"
                    "Please ensure mHM is installed and accessible.")
                return False
                
        except Exception as e:
            self.log_message(f"ERROR: Failed to run mHM. Details: {str(e)}")
            QMessageBox.critical(
                self.dialog, "Error", 
                f"Failed to run mHM.\n{str(e)}")
            return False

