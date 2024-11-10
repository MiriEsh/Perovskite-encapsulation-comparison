import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import seaborn
import scipy
import os


if __name__ == '__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__))
    output_path = dir_path + "\\outputs"
    Perovskites = pd.read_csv(dir_path + "\\inputs\\Perovsite database query.csv")

    #adding encapsulation info retrieved from the articles:
    additional_glue_info = pd.read_csv(dir_path + "\\inputs\\additions_to_Encapsulation_info.csv")
    for i, row in additional_glue_info.iterrows():
        Perovskites.loc[Perovskites.Ref_ID == additional_glue_info.Ref_ID[i], ['Encapsulation_stack_sequence', 'Encapsulation_edge_sealing_materials']] = \
        additional_glue_info.new_Encapsulation_stack_sequence[i], additional_glue_info.new_Encapsulation_edge_sealing_materials[i]
    print('added/changed encapsulation info using: additions_to_Encapsulation_info.csv')
    # adding stability info from 5 articles in the DB:
    additional_PCE_info = pd.read_csv(dir_path + "\\inputs\\additions_to_PCE_info.csv")
    for i, row in additional_PCE_info.iterrows():
        Perovskites.loc[Perovskites.Ref_ID == additional_PCE_info.Ref_ID[i], ['Stability_time_total_exposure',
                                                                               'Stability_PCE_end_of_experiment']] =   additional_PCE_info.new_Stability_time_total_exposure[i], additional_PCE_info.new_Stability_PCE_end_of_experiment[i]
    print('added stability info info using: additions_to_PCE_info.csv')
    #taking out the rows with stability atmosphere N2 or Ar:
    num_of_rows_before_taking_out_atmosphereN2 = len(Perovskites)
    Perovskites = Perovskites[(Perovskites.Stability_atmosphere != 'N2') & (Perovskites.Stability_atmosphere != 'Ar') & (Perovskites.Stability_atmosphere != 'Water')]
    num_of_rows = len(Perovskites)
    num_of_rows_with_atmosphere_N2 = num_of_rows_before_taking_out_atmosphereN2 - num_of_rows
    print('took out ', num_of_rows_with_atmosphere_N2 , ' rows with Stability atmosphere N2 or Ar')

    #taking out rows where the articles did not have any encapsulation
    Perovskites = Perovskites[Perovskites['Ref_DOI_number'].str.contains('j.tsf.2020.137786|c8tc05773g|acs.nanolett.8b00025|s41598-018-19612-7|acs.nanolett.8b03685|s40820-017-0140-x|c7ra07579k|adfm.201900466|acs.nanolett.8b01440|adfm.201904684|acs.nanolett.8b01440|c8ra08334g|adfm.201900466|j.orgel.2019.105387',case=False)==False]
    num_after_deleting_nonencapsulted = len(Perovskites)
    print('took out ', num_of_rows- num_after_deleting_nonencapsulted, ' rows with nonencapulated cells')

    #taking out rows with stability mode equal to Short circuit
    Perovskites = Perovskites[Perovskites.Stability_potential_bias_load_condition != 'Short circuit']
    num_after_deleting_short_circuit = len(Perovskites)
    print('took out ',num_after_deleting_nonencapsulted - num_after_deleting_short_circuit, 'rows with stability mode Short circuit')

    #concatenate the Ecapsulation_stack_sequence and the Encapsulation_edge_sealing_materials:
    Perovskites['Encapsulation_materials'] = Perovskites['Encapsulation_stack_sequence'] + ' ' + Perovskites['Encapsulation_edge_sealing_materials']

    #take out all the 'Unknown' and 'none' from the materials:
    Perovskites['Encapsulation_materials'] = Perovskites['Encapsulation_materials'].str.replace('Unknown|none', '')

    #take out rows with clamp as a material
    Perovskites = Perovskites[Perovskites['Encapsulation_materials'].str.contains('Clamp', case=False) == False]
    num_of_rows_with_clamp = num_after_deleting_short_circuit - len(Perovskites)
    print('took out ', num_of_rows_with_clamp, ' rows with Clamp in the encapsulation materials')
    PCE_end_greater_than_100 =Perovskites[Perovskites['Stability_PCE_end_of_experiment'] >= 100]
    PCE_end_greater_than_100.to_csv(output_path+"\\Experiments_with_T80_greater_than_100.csv")
    Perovskites = Perovskites.drop(Perovskites[Perovskites['Stability_PCE_end_of_experiment'] >= 100].index)
    # if Stability_PCE_end_of_experiment=0 take Stability_PCE_end_of_experiment=0.001 in order to prevent converging ln(0) to -infinity
    Perovskites.loc[Perovskites['Stability_PCE_end_of_experiment'] == 0, 'Stability_PCE_end_of_experiment'] = 0.001

    #BURN IN FALSE:
    # TAKE T80
    # T80 doesn't exist take Ts80
    # T80, Ts80 don't exist use T95 via: ln(0.8) * T95/ ln(0.95)
    # T80, Ts80, T95 don't exist ->CALC VIA EXP
    Perovskites.loc[Perovskites['Stability_PCE_burn_in_observed'] == False, 'Calculated_T80'] = np.log(0.8) * Perovskites['Stability_time_total_exposure']/np.log(Perovskites['Stability_PCE_end_of_experiment']/100.0)
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == False) & (Perovskites['Stability_PCE_T95'].isna()==False)), 'Calculated_T80'] = np.log(0.8) * Perovskites['Stability_PCE_T95'] / np.log(0.95)
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == False) & (Perovskites['Stability_PCE_Ts80'].isna()==False)), 'Calculated_T80'] = Perovskites['Stability_PCE_Ts80']
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == False) & (Perovskites['Stability_PCE_T80'].isna()==False)), 'Calculated_T80'] = Perovskites['Stability_PCE_T80']
    #BURN IN TRUE:
    #TAKE TS80
    #TS80 doesn't exist - use T80
    #T80 doesn't exist- use Te80
    #TS80/Te80 don't exist CALC VIA EXP using PCE (1000H)+ END POINT
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == True) & (Perovskites['Stability_PCE_after_1000_h'].isna()==False) &
    (Perovskites['Stability_PCE_after_1000_h'] != Perovskites['Stability_PCE_end_of_experiment'])), 'Calculated_T80'] = \
    -np.log(0.8)*(Perovskites['Stability_time_total_exposure']-1000.0)/np.log(Perovskites['Stability_PCE_after_1000_h']/Perovskites['Stability_PCE_end_of_experiment'])
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == True) & (Perovskites[
        'Stability_PCE_Te80'].isna() == False)), 'Calculated_T80'] = Perovskites['Stability_PCE_Te80']
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == True) & (Perovskites[
        'Stability_PCE_Ts80'].isna() == False)), 'Calculated_T80'] = Perovskites['Stability_PCE_Ts80']

    Perovskites['logT80'] = np.log10(Perovskites['Calculated_T80'])
    Perovskites = Perovskites.drop(Perovskites[Perovskites['Calculated_T80'].isna()].index)

    print('There are a total of', len(Perovskites), 'with stability info')

    experiments_with_T80 = Perovskites[
        ['Ref_ID', 'Ref_DOI_number', 'Perovskite_composition_short_form', 'JV_reverse_scan_PCE', 'JV_forward_scan_PCE',
         'JV_default_PCE',
         'Stability_PCE_initial_value',
         'JV_default_PCE_scan_direction',
         'Stabilised_performance_PCE',
         'Stability_PCE_end_of_experiment', 'Stability_time_total_exposure',
         'Perovskite_thickness', 'Perovskite_deposition_thermal_annealing_temperature',
         'Perovskite_deposition_thermal_annealing_time', 'Calculated_T80']]
    Perovskites.to_csv(output_path+"\\Experiments_with_T80_for_supp.csv")

    # classifying to groups:
    Perovskites['classification'] = ''
    Other = Perovskites[Perovskites['Encapsulation_materials'].str.contains('glass cans|Graphene|Field|particles|Polymer hydrophobic film|Eu(TTA)2(Phen)MAA|Carbon-nt|3M acrylic elastomer (3M VHB 4905) ', case=False) == True]
    Other['classification'] = 'Other'
    Perovskites = Perovskites[~Perovskites.index.isin(Other.index)]
    Al2O3 = Perovskites[Perovskites['Encapsulation_materials'].str.contains('Al2O3', case=False) == True]
    Al2O3['classification'] = 'Al2O3'
    Perovskites = Perovskites[~Perovskites.index.isin(Al2O3.index)]
    Glass = Perovskites[Perovskites['Encapsulation_materials'].str.contains('SLG|Glass', case= False)== True]
    Glass_epoxy = Glass[Glass['Encapsulation_materials'].str.contains('epoxy|Araldite|Ossila',case =False) == True]
    Glass_epoxy['classification'] = 'Glass + Epoxy'
    Glass = Glass[~Glass.index.isin(Glass_epoxy.index)]  # taking only non epoxy glass rows so that UV epoxy is classified as epoxy and not as UV-glue
    Glass_UV_glue = Glass[Glass['Encapsulation_materials'].str.contains('UVCA|UV|ultraviolet|Light| LT',case =False) == True]
    Glass_UV_glue['classification'] = 'Glass + UV glue'
    Glass = Glass[~Glass.index.isin(Glass_UV_glue.index)]  # so that if there is both UV glue and polymer (24 rows), it is cllasified as UV-glue and not as polymer
    Glass_butyl_rubber = Glass[Glass['Encapsulation_materials'].str.contains('Polyisobutene|Polyisobutylene|butyl', case=False) == True]
    Glass_butyl_rubber['classification'] = 'Glass + Butyl rubber'
    Glass = Glass[~Glass.index.isin(Glass_butyl_rubber.index)]
    Glass_polymer = Glass[Glass['Encapsulation_materials'].str.contains('Thermoplastic|poly|surlyn|EVA|Ethylene|Parylene|Kapton|NOA',case =False) == True]
    Glass_polymer['classification'] = 'Glass + polymer'
    Glass_No_sealant = Glass[Glass['Encapsulation_materials'].str.contains('Thermoplastic|poly|surlyn|EVA|Ethylene|Parylene|Kapton|Polyisobutylene|butyl|UV|ultraviolet|Light| LT|NOA|epoxy|Araldite|Ossila', case=False) == False]
    Glass_No_sealant['classification'] = 'Glass no sealant'
    Non_Glass = Perovskites[Perovskites['Encapsulation_materials'].str.contains('SLG|Glass', case=False) == False]
    Foils= Non_Glass[Non_Glass['Encapsulation_materials'].str.contains('tape|plasma polymers|viewbarrier|Barrier foil|sheets|UHPBF|ITO|Kapton|PEN|PET|NOA', case= False)== True]
    Foils['classification'] = 'Foils'
    Non_Glass = Non_Glass[~Non_Glass.index.isin(Foils.index)]
    polymer = Non_Glass[Non_Glass['Encapsulation_materials'].str.contains('PVP |Cyanoacrylate|PCL|Carbon-epoxy|Teflon|PDMS|Paraffin|LDPE|Meltronix|Pattex|Thermoplastic|poly|surlyn|EVA|Ethylene|Parylene|EVOH|PMMA',case=False) == True]
    polymer['classification'] = 'Polymers'
    # finished classifying to groups

    info_for_9_categories = pd.concat(
        [Other, Glass_epoxy,Glass_UV_glue, Glass_polymer, Glass_butyl_rubber, Glass_No_sealant, Al2O3, Foils, polymer])
    all_info_for_analysis = info_for_9_categories[['Ref_ID', 'Ref_DOI_number', 'Encapsulation_stack_sequence', 'Encapsulation_edge_sealing_materials', 'Stability_time_total_exposure',
                        'Stability_PCE_end_of_experiment', 'Stability_PCE_T95','Stability_PCE_Ts95','Stability_PCE_T80','Stability_PCE_Ts80','Stability_PCE_Te80',
                        'Stability_PCE_Tse80','Stability_PCE_after_1000_h','Stability_PCE_burn_in_observed','Stability_light_source_type','Stability_protocol',
                        'Stability_potential_bias_load_condition', 'Calculated_T80','logT80', 'classification']]
    info_for_9_categories.to_csv(output_path+"\\Experiments_with_Encap_and_stability.csv")

    Graph_titles= ["All experiments","Ambient experiments","Experiments at 65 degrees","Dark Open Circuit Experiments","Light Experiments","Light MPPT Experiments"]
    #preparing data for the graphs X=encapsulation_group X[0]-> all experiments, X[1]-> ambient experiments, X[2]-> experiments in 65deg,
    # X[3]-> dark OC experiments, X[4]->light OC experiments, X[5]-> light MPP experiments
    Glass_butyl_rubber_T80=[]
    Glass_butyl_rubber_T80.append(Glass_butyl_rubber['Calculated_T80'].reset_index(drop=True))
    Glass_butyl_rubber_T80.append(Glass_butyl_rubber[Glass_butyl_rubber['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    Glass_butyl_rubber_T80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_protocol'].str.contains('2|3'))&
                               (~Glass_butyl_rubber['Stability_protocol'].str.contains('IEC'))]['Calculated_T80'].reset_index(drop=True))
    Glass_butyl_rubber_T80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_light_source_type']=='Dark')|
                             (Glass_butyl_rubber['Stability_protocol'].str.contains('D',regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_butyl_rubber_T80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_light_source_type']!='Dark')&(~Glass_butyl_rubber['Stability_protocol'].str.contains('D',regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_butyl_rubber_T80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_potential_bias_load_condition'] == 'MPPT') |
                                            (Glass_butyl_rubber['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))

    Glass_butyl_rubber_logT80 = []
    Glass_butyl_rubber_logT80.append(Glass_butyl_rubber['logT80'].reset_index(drop=True))
    Glass_butyl_rubber_logT80.append(
        Glass_butyl_rubber[Glass_butyl_rubber['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(
            drop=True))
    Glass_butyl_rubber_logT80.append(
        Glass_butyl_rubber[(Glass_butyl_rubber['Stability_protocol'].str.contains('2|3')) &
                           (~Glass_butyl_rubber['Stability_protocol'].str.contains('IEC'))]['logT80'].reset_index(
            drop=True))
    Glass_butyl_rubber_logT80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_light_source_type'] == 'Dark') |
                                                         (Glass_butyl_rubber['Stability_protocol'].str.contains('D',regex=False))]['logT80'].reset_index(drop=True))
    Glass_butyl_rubber_logT80.append(Glass_butyl_rubber[
                                          (Glass_butyl_rubber['Stability_light_source_type'] != 'Dark') & (
                                              ~Glass_butyl_rubber['Stability_protocol'].str.contains('D',regex=False))]['logT80'].reset_index(drop=True))
    Glass_butyl_rubber_logT80.append(
        Glass_butyl_rubber[(Glass_butyl_rubber['Stability_potential_bias_load_condition'] == 'MPPT') |
                           (Glass_butyl_rubber['Stability_potential_bias_load_condition'] == 'Passive resistance')]['logT80'].reset_index(drop=True))

    Glass_polymer_T80=[]
    Glass_polymer_T80.append(Glass_polymer['Calculated_T80'].reset_index(drop=True))
    Glass_polymer_T80.append(Glass_polymer[Glass_polymer['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    Glass_polymer_T80.append(Glass_polymer[Glass_polymer['Stability_protocol'].str.contains('2|3')]['Calculated_T80'].reset_index(drop=True))
    Glass_polymer_T80.append(Glass_polymer[(Glass_polymer['Stability_light_source_type'] == 'Dark') |
                                   (Glass_polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_polymer_T80.append(Glass_polymer[(Glass_polymer['Stability_light_source_type'] != 'Dark')&(~Glass_polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_polymer_T80.append(Glass_polymer[(Glass_polymer['Stability_potential_bias_load_condition'] == 'MPPT') |
                                      (Glass_polymer['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))
    Glass_polymer_logT80 = []
    Glass_polymer_logT80.append(Glass_polymer['logT80'].reset_index(drop=True))
    Glass_polymer_logT80.append(
        Glass_polymer[Glass_polymer['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(drop=True))
    Glass_polymer_logT80.append(
        Glass_polymer[Glass_polymer['Stability_protocol'].str.contains('2|3')]['logT80'].reset_index(drop=True))
    Glass_polymer_logT80.append(Glass_polymer[(Glass_polymer['Stability_light_source_type'] == 'Dark') |
                                               (Glass_polymer['Stability_protocol'].str.contains('D', regex=False))][
                                     'logT80'].reset_index(drop=True))
    Glass_polymer_logT80.append(Glass_polymer[(Glass_polymer['Stability_light_source_type'] != 'Dark') & (
        ~Glass_polymer['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Glass_polymer_logT80.append(Glass_polymer[(Glass_polymer['Stability_potential_bias_load_condition'] == 'MPPT') |
                                               (Glass_polymer['Stability_potential_bias_load_condition'] == 'Passive resistance')]['logT80'].reset_index(drop=True))

    Al2O3_T80=[]
    Al2O3_T80.append(Al2O3['Calculated_T80'].reset_index(drop=True))
    Al2O3_T80.append(Al2O3[Al2O3['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    Al2O3_T80.append(Al2O3[Al2O3['Stability_protocol'].str.contains('2|3')]['Calculated_T80'].reset_index(drop=True))
    Al2O3_T80.append(Al2O3[(Al2O3['Stability_light_source_type'] == 'Dark') |
                       (Al2O3['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Al2O3_T80.append(Al2O3[(Al2O3['Stability_light_source_type'] != 'Dark')&(~Al2O3['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Al2O3_T80.append(Al2O3[(Al2O3['Stability_potential_bias_load_condition'] == 'MPPT') |
                      (Al2O3['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))

    Al2O3_logT80 = []
    Al2O3_logT80.append(Al2O3['logT80'].reset_index(drop=True))
    Al2O3_logT80.append(Al2O3[Al2O3['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(drop=True))
    Al2O3_logT80.append(Al2O3[Al2O3['Stability_protocol'].str.contains('2|3')]['logT80'].reset_index(drop=True))
    Al2O3_logT80.append(Al2O3[(Al2O3['Stability_light_source_type'] == 'Dark') |
                               (Al2O3['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Al2O3_logT80.append(Al2O3[(Al2O3['Stability_light_source_type'] != 'Dark') & (
        ~Al2O3['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Al2O3_logT80.append(Al2O3[(Al2O3['Stability_potential_bias_load_condition'] == 'MPPT') |
                               (Al2O3['Stability_potential_bias_load_condition'] == 'Passive resistance')][
                             'logT80'].reset_index(drop=True))

    Foils_T80=[]
    Foils_T80.append(Foils['Calculated_T80'].reset_index(drop=True))
    Foils_T80.append(Foils[Foils['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    Foils_T80.append(Foils[Foils['Stability_protocol'].str.contains('2|3')]['Calculated_T80'].reset_index(drop=True))
    Foils_T80.append(Foils[(Foils['Stability_light_source_type'] == 'Dark') |
                       (Foils['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Foils_T80.append(Foils[(Foils['Stability_light_source_type'] != 'Dark')&(~Foils['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Foils_T80.append(Foils[(Foils['Stability_potential_bias_load_condition'] == 'MPPT') |
                      (Foils['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))

    Foils_logT80 = []
    Foils_logT80.append(Foils['logT80'].reset_index(drop=True))
    Foils_logT80.append(Foils[Foils['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(drop=True))
    Foils_logT80.append(Foils[Foils['Stability_protocol'].str.contains('2|3')]['logT80'].reset_index(drop=True))
    Foils_logT80.append(Foils[(Foils['Stability_light_source_type'] == 'Dark') |
                               (Foils['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Foils_logT80.append(Foils[(Foils['Stability_light_source_type'] != 'Dark') & (
        ~Foils['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Foils_logT80.append(Foils[(Foils['Stability_potential_bias_load_condition'] == 'MPPT') |
                               (Foils['Stability_potential_bias_load_condition'] == 'Passive resistance')]['logT80'].reset_index(drop=True))

    polymer_T80=[]
    polymer_T80.append(polymer['Calculated_T80'].reset_index(drop=True))
    polymer_T80.append(polymer[polymer['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    polymer_T80.append(polymer[polymer['Stability_protocol'].str.contains('2|3')]['Calculated_T80'].reset_index(drop=True))
    polymer_T80.append(polymer[(polymer['Stability_light_source_type'] == 'Dark') |
                           (polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    polymer_T80.append(polymer[(polymer['Stability_light_source_type'] != 'Dark')&(~polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    polymer_T80.append(polymer[(polymer['Stability_potential_bias_load_condition'] == 'MPPT') |
                          (polymer['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))

    polymer_logT80 = []
    polymer_logT80.append(polymer['logT80'].reset_index(drop=True))
    polymer_logT80.append(
        polymer[polymer['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(drop=True))
    polymer_logT80.append(polymer[polymer['Stability_protocol'].str.contains('2|3')]['logT80'].reset_index(drop=True))
    polymer_logT80.append(polymer[(polymer['Stability_light_source_type'] == 'Dark') |
                                   (polymer['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    polymer_logT80.append(polymer[(polymer['Stability_light_source_type'] != 'Dark') & (
        ~polymer['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    polymer_logT80.append(polymer[(polymer['Stability_potential_bias_load_condition'] == 'MPPT') |
                                   (polymer['Stability_potential_bias_load_condition'] == 'Passive resistance')][
                               'logT80'].reset_index(drop=True))

    Glass_UV_glue_T80=[]
    Glass_UV_glue_T80.append(Glass_UV_glue['Calculated_T80'].reset_index(drop=True))
    Glass_UV_glue_T80.append(Glass_UV_glue[Glass_UV_glue['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    Glass_UV_glue_T80.append(Glass_UV_glue[Glass_UV_glue['Stability_protocol'].str.contains('2|3')]['Calculated_T80'].reset_index(drop=True))
    Glass_UV_glue_T80.append(Glass_UV_glue[(Glass_UV_glue['Stability_light_source_type'] == 'Dark') |
                                       (Glass_UV_glue['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_UV_glue_T80.append(Glass_UV_glue[(Glass_UV_glue['Stability_light_source_type'] != 'Dark')&(~Glass_UV_glue['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_UV_glue_T80.append(Glass_UV_glue[(Glass_UV_glue['Stability_potential_bias_load_condition'] == 'MPPT') |
                                      (Glass_UV_glue['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))

    Glass_UV_glue_logT80 = []
    Glass_UV_glue_logT80.append(Glass_UV_glue['logT80'].reset_index(drop=True))
    Glass_UV_glue_logT80.append(
        Glass_UV_glue[Glass_UV_glue['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(drop=True))
    Glass_UV_glue_logT80.append(
        Glass_UV_glue[Glass_UV_glue['Stability_protocol'].str.contains('2|3')]['logT80'].reset_index(drop=True))
    Glass_UV_glue_logT80.append(Glass_UV_glue[(Glass_UV_glue['Stability_light_source_type'] == 'Dark') |
                                               (Glass_UV_glue['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Glass_UV_glue_logT80.append(Glass_UV_glue[(Glass_UV_glue['Stability_light_source_type'] != 'Dark') & (
        ~Glass_UV_glue['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Glass_UV_glue_logT80.append(Glass_UV_glue[(Glass_UV_glue['Stability_potential_bias_load_condition'] == 'MPPT') |
                                               (Glass_UV_glue['Stability_potential_bias_load_condition'] == 'Passive resistance')]['logT80'].reset_index(drop=True))

    Glass_epoxy_T80=[]
    Glass_epoxy_T80.append(Glass_epoxy['Calculated_T80'].reset_index(drop=True))
    Glass_epoxy_T80.append(Glass_epoxy[Glass_epoxy['Stability_protocol'].str.contains('1|Other')]['Calculated_T80'].reset_index(drop=True))
    Glass_epoxy_T80.append(Glass_epoxy[Glass_epoxy['Stability_protocol'].str.contains('2|3')]['Calculated_T80'].reset_index(drop=True))
    Glass_epoxy_T80.append(Glass_epoxy[(Glass_epoxy['Stability_light_source_type'] == 'Dark') |
                                   (Glass_epoxy['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_epoxy_T80.append(Glass_epoxy[(Glass_epoxy['Stability_light_source_type'] != 'Dark')&(~Glass_epoxy['Stability_protocol'].str.contains('D', regex=False))]['Calculated_T80'].reset_index(drop=True))
    Glass_epoxy_T80.append(Glass_epoxy[(Glass_epoxy['Stability_potential_bias_load_condition'] == 'MPPT') |
                                  (Glass_epoxy['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_T80'].reset_index(drop=True))

    Glass_epoxy_logT80 = []
    Glass_epoxy_logT80.append(Glass_epoxy['logT80'].reset_index(drop=True))
    Glass_epoxy_logT80.append(
        Glass_epoxy[Glass_epoxy['Stability_protocol'].str.contains('1|Other')]['logT80'].reset_index(drop=True))
    Glass_epoxy_logT80.append(
        Glass_epoxy[Glass_epoxy['Stability_protocol'].str.contains('2|3')]['logT80'].reset_index(drop=True))
    Glass_epoxy_logT80.append(Glass_epoxy[(Glass_epoxy['Stability_light_source_type'] == 'Dark') |
                                           (Glass_epoxy['Stability_protocol'].str.contains('D', regex=False))][
                                   'logT80'].reset_index(drop=True))
    Glass_epoxy_logT80.append(Glass_epoxy[(Glass_epoxy['Stability_light_source_type'] != 'Dark') & (
        ~Glass_epoxy['Stability_protocol'].str.contains('D', regex=False))]['logT80'].reset_index(drop=True))
    Glass_epoxy_logT80.append(Glass_epoxy[(Glass_epoxy['Stability_potential_bias_load_condition'] == 'MPPT') |
                                           (Glass_epoxy['Stability_potential_bias_load_condition'] == 'Passive resistance')]['logT80'].reset_index(drop=True))

    table = []
    table.append(['', 'Foil','Polymer','Glass UV glue','Glass epoxy','Glass butyl rubber','Glass Polymer','Al2O3', 'Total'])
    my_colors = {'Foil': 'teal', 'Polymer': 'steelblue', 'SLG+UV glue': 'slateblue', 'SLG+epoxy':'darkorchid','SLG+but_rub':'mediumvioletred' ,'SLG+polym': 'sienna', 'Al2O3':'darkkhaki' }
    for i in range(6):
        data_for_box_plot = pd.concat([Foils_T80[i],polymer_T80[i], Glass_UV_glue_T80[i], Glass_epoxy_T80[i],Glass_butyl_rubber_T80[i],Glass_polymer_T80[i],Al2O3_T80[i]], axis =1)  # horizontal_concat
        data_for_box_plot=data_for_box_plot.set_axis(["Foil","Polymer","SLG+UV glue","SLG+epoxy","SLG+but_rub", "SLG+polym", "Al2O3"], axis="columns")
        plt.rcParams.update({'font.size': 24})
        seaborn.catplot(x="variable", y="value", data=pd.melt(data_for_box_plot.dropna(how= "all", axis=1)), kind="box",palette=my_colors)
        seaborn.swarmplot(x='variable', y='value', data=pd.melt(data_for_box_plot.dropna(how= "all", axis=1)), color="gray")
        seaborn.pointplot(data=pd.melt(data_for_box_plot.dropna(how= "all", axis=1)), x="variable", y="value", estimator=np.mean,
                          join=False, ci=None, markers='+', color='black', zorder=3)
        plt.ylabel("T80")
        plt.yscale('log')
        plt.minorticks_on()
        plt.xticks(rotation=75)
        plt.xlabel("")
        plt.title(Graph_titles[i])
        table.append(['num of experiments'+Graph_titles[i],len(Foils_logT80[i]), len(polymer_logT80[i]), len(Glass_UV_glue_logT80[i]), len(Glass_epoxy_logT80[i]),len(Glass_butyl_rubber_logT80[i]), len(Glass_polymer_logT80[i]),
                      len(Al2O3_logT80[i])])
        table.append(['average'+Graph_titles[i], np.average(Foils_logT80[i]), np.average(polymer_logT80[i]),np.average(Glass_UV_glue_logT80[i]),
                      np.average(Glass_epoxy_logT80[i]),np.average(Glass_butyl_rubber_logT80[i]), np.average(Glass_polymer_logT80[i]), np.average(Al2O3_logT80[i])])
        table.append(['stdev'+Graph_titles[i], np.std(Foils_logT80[i]),np.std(polymer_logT80[i]), np.std(Glass_UV_glue_logT80[i]), np.std(Glass_epoxy_logT80[i]),np.std(Glass_butyl_rubber_logT80[i]), np.std(Glass_polymer_logT80[i]), np.std(Al2O3_logT80[i])])

        t_test = []  # currently only printing t-test combinations of groups with more than 5 experiments - if more experiments are added after Sep. 3rd 2024 change this and add all combinations
        t_test.append([None,'Polymer','Glass+UV glue','Glass+epoxy','Glass +butyl rubber','Glass +polymer','Al2O3'])
        t_test.append(['Foil', scipy.stats.ttest_ind(Foils_logT80[i], polymer_logT80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(Foils_logT80[i], Glass_UV_glue_logT80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(Foils_logT80[i], Glass_epoxy_logT80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(Foils_logT80[i], Glass_butyl_rubber_logT80[i], axis=0, equal_var=False).pvalue,
                       scipy.stats.ttest_ind(Foils_logT80[i], Glass_polymer_logT80[i], axis=0, equal_var=False).pvalue,scipy.stats.ttest_ind(Foils_logT80[i], Al2O3_logT80[i], axis=0, equal_var=False).pvalue])
        t_test.append(['Polymer',None, scipy.stats.ttest_ind(polymer_logT80[i], Glass_UV_glue_logT80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(polymer_logT80[i], Glass_epoxy_logT80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(polymer_logT80[i], Glass_butyl_rubber_logT80[i], axis=0, equal_var=False).pvalue,
                       scipy.stats.ttest_ind(polymer_logT80[i], Glass_polymer_logT80[i], axis=0, equal_var=False).pvalue,scipy.stats.ttest_ind(polymer_logT80[i], Al2O3_logT80[i], axis=0, equal_var=False).pvalue])
        t_test.append(['Glass+UV_glue',None, None, scipy.stats.ttest_ind(Glass_UV_glue_logT80[i], Glass_epoxy_logT80[i], axis=0,equal_var=False).pvalue,
                           scipy.stats.ttest_ind(Glass_UV_glue_logT80[i], Glass_butyl_rubber_logT80[i], axis=0,equal_var=False).pvalue,
                           scipy.stats.ttest_ind(Glass_UV_glue_logT80[i], Glass_polymer_logT80[i], axis=0, equal_var=False).pvalue,
                           scipy.stats.ttest_ind(Glass_UV_glue_logT80[i], Al2O3_logT80[i], axis=0, equal_var=False).pvalue])
        t_test.append(['Glass+epoxy',None, None, None,
                                scipy.stats.ttest_ind(Glass_epoxy_logT80[i], Glass_butyl_rubber_logT80[i], axis=0,equal_var=False).pvalue,
                                scipy.stats.ttest_ind(Glass_epoxy_logT80[i], Glass_polymer_logT80[i], axis=0,equal_var=False).pvalue,
                                scipy.stats.ttest_ind(Glass_epoxy_logT80[i], Al2O3_logT80[i], axis=0,equal_var=False).pvalue])
        t_test.append(['Glass+butyl rubber',None, None, None, None,
                              scipy.stats.ttest_ind(Glass_butyl_rubber_logT80[i], Glass_polymer_logT80[i], axis=0,equal_var=False).pvalue,
                              scipy.stats.ttest_ind(Glass_butyl_rubber_logT80[i], Al2O3_logT80[i], axis=0,equal_var=False).pvalue])
        t_test.append(['Glass+polymer',None, None, None, None, None,
                                     scipy.stats.ttest_ind(Glass_polymer_logT80[i], Al2O3_logT80[i], axis=0,equal_var=False).pvalue])
        t_test = pd.DataFrame(t_test)
        t_test.to_csv(output_path +"\\t_test_"+Graph_titles[i] +"_logT80.csv")
    plt.show()
    table= pd.DataFrame(table)
    table.to_csv(output_path +"\\frequency table logT80.csv")
