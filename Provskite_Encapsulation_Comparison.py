import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import seaborn
import scipy

if __name__ == '__main__':
    Perovskites = pd.read_csv('Perovsite database query.csv')

    #adding encapsulation info retrieved from the articles:
    additional_glue_info = pd.read_csv(
        'additions_to_Encapsulation_stack_sequence_and_Encapsulation_sealing_material.csv')
    for i, row in additional_glue_info.iterrows():
        Perovskites.loc[Perovskites.Ref_ID == additional_glue_info.Ref_ID[i], ['Encapsulation_stack_sequence', 'Encapsulation_edge_sealing_materials']] = \
        additional_glue_info.new_Encapsulation_stack_sequence[i], additional_glue_info.new_Encapsulation_edge_sealing_materials[i]
    print('added/changed encapsulation info from: additions_to_Encapsulation_stack_sequence_and_Encapsulation_sealing_material.csv')
    # adding stability info from 5 articles in the DB:
    Perovskites.loc[Perovskites.Ref_ID == 22840, ['Stability_time_total_exposure','Stability_PCE_end_of_experiment']] = 90, 110
    Perovskites.loc[Perovskites.Ref_ID == 43586, ['Stability_time_total_exposure', 'Stability_PCE_end_of_experiment']] = 3260, 90
    Perovskites.loc[Perovskites.Ref_ID == 22379, ['Stability_PCE_end_of_experiment']] = 95
    Perovskites.loc[Perovskites.Ref_ID == 42394, ['Stability_PCE_end_of_experiment']] = 95
    Perovskites.loc[Perovskites.Ref_ID == 38323, ['Stability_PCE_end_of_experiment']] = 86
    print('added stability info from 5 articles: 10.1039/c9ta01859j, 10.1021/acs.jpclett.0c00923, 10.1021/acsami.9b23532,10.1021/acsami.7b07625,10.1016/j.xcrp.2021.100648')

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

    Perovskites = Perovskites.drop(Perovskites[Perovskites['Stability_PCE_end_of_experiment'] >= 100].index)
    # if Stability_PCE_end_of_experiment=0 take Stability_PCE_end_of_experiment=0.1 in order to prevent converging to ln(0)=-infinity
    Perovskites.loc[Perovskites['Stability_PCE_end_of_experiment'] == 0, 'Stability_PCE_end_of_experiment'] = 0.001

    #BURN IN FALSE:
    # TAKE T80
    # T80 DOESN'T EXIST THEN TAKE TS80 (for burn-in= false)
    # T80, Ts80 don't exist USE T95 via: ln(0.8) * T95/ ln(0.95)
    # T80, Ts80, T95 don't exist ->CALC VIA EXP
    Perovskites.loc[Perovskites['Stability_PCE_burn_in_observed'] == False, 'Calculated_Ts80'] = np.log(0.8) * Perovskites['Stability_time_total_exposure']/np.log(Perovskites['Stability_PCE_end_of_experiment']/100.0)
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == False) & (Perovskites['Stability_PCE_T95'].isna()==False)), 'Calculated_Ts80'] = np.log(0.8) * Perovskites['Stability_PCE_T95'] / np.log(0.95)
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == False) & (Perovskites['Stability_PCE_Ts80'].isna()==False)), 'Calculated_Ts80'] = Perovskites['Stability_PCE_Ts80']
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == False) & (Perovskites['Stability_PCE_T80'].isna()==False)), 'Calculated_Ts80'] = Perovskites['Stability_PCE_T80']
    #BURN IN TRUE:
    #TAKE TS80
    #T80 DOESN'T EXIST THEN TAKE Te80 (for burn-in= true row)
    #TS80/Te80 don't exist CALC VIA PCE (1000H)+ END POINT
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == True) & (Perovskites['Stability_PCE_after_1000_h'].isna()==False) &
    (Perovskites['Stability_PCE_after_1000_h'] != Perovskites['Stability_PCE_end_of_experiment'])), 'Calculated_Ts80'] = \
    -np.log(0.8)*(Perovskites['Stability_time_total_exposure']-1000.0)/np.log(Perovskites['Stability_PCE_after_1000_h']/Perovskites['Stability_PCE_end_of_experiment'])
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == True) & (Perovskites[
        'Stability_PCE_Te80'].isna() == False)), 'Calculated_Ts80'] = Perovskites['Stability_PCE_Te80']
    Perovskites.loc[((Perovskites['Stability_PCE_burn_in_observed'] == True) & (Perovskites[
        'Stability_PCE_Ts80'].isna() == False)), 'Calculated_Ts80'] = Perovskites['Stability_PCE_Ts80']

    Perovskites = Perovskites.drop(Perovskites[Perovskites['Calculated_Ts80'].isna()].index)
    Perovskites = Perovskites.drop(Perovskites[Perovskites['Calculated_Ts80'] ==1].index)

    print('There are a total of', len(Perovskites), 'with stability info')

    experiments_with_Ts80 = Perovskites[
        ['Ref_ID', 'Ref_DOI_number', 'Perovskite_composition_short_form', 'JV_reverse_scan_PCE', 'JV_forward_scan_PCE',
         'JV_default_PCE',
         'Stability_PCE_initial_value',
         'JV_default_PCE_scan_direction',
         'Stabilised_performance_PCE',
         'Stability_PCE_end_of_experiment', 'Stability_time_total_exposure',
         'Perovskite_thickness', 'Perovskite_deposition_thermal_annealing_temperature',
         'Perovskite_deposition_thermal_annealing_time', 'Calculated_Ts80']]
    Perovskites.to_csv('422_Experiments_with_Ts80_for_supp.csv')

    # classifying to groups:
    Other = Perovskites[Perovskites['Encapsulation_materials'].str.contains('glass cans|Graphene|Field|particles|polyethylene terephthalate|Polymer hydrophobic film|MAA|Carbon-nt|acrylic elastomer', case=False) == True]
    Perovskites = Perovskites[~Perovskites.index.isin(Other.index)]
    Al2O3 = Perovskites[Perovskites['Encapsulation_materials'].str.contains('Al2O3', case=False) == True]
    Perovskites = Perovskites[~Perovskites.index.isin(Al2O3.index)]
    Glass = Perovskites[Perovskites['Encapsulation_materials'].str.contains('SLG|Glass', case= False)== True]
    Glass_epoxy = Glass[Glass['Encapsulation_materials'].str.contains('epoxy|Araldite|Ossila',case =False) == True]
    Glass = Glass[~Glass.index.isin(Glass_epoxy.index)]  # taking only non epoxy glass rows so that UV epoxy is classified as epoxy and not as UV-glue
    Glass_UV_glue = Glass[Glass['Encapsulation_materials'].str.contains('UVCA|UV|ultraviolet|Light| LT',case =False) == True]
    Glass = Glass[~Glass.index.isin(Glass_UV_glue.index)]  # so that if there is both UV glue and polymer (24 rows), it is cllasified as UV-glue and not as polymer
    Glass_butyl_rubber = Glass[Glass['Encapsulation_materials'].str.contains('Polyisobutene|Polyisobutylene|butyl', case=False) == True]
    Glass = Glass[~Glass.index.isin(Glass_butyl_rubber.index)]
    Glass_polymer = Glass[Glass['Encapsulation_materials'].str.contains('Thermoplastic|poly|surlyn|EVA|Ethylene|Parylene|Kapton|NOA',case =False) == True]
    Glass_No_sealant = Glass[Glass['Encapsulation_materials'].str.contains('Thermoplastic|poly|surlyn|EVA|Ethylene|Parylene|Kapton|Polyisobutylene|butyl|UV|ultraviolet|Light| LT|NOA|epoxy|Araldite|Ossila', case=False) == False]
    Non_Glass = Perovskites[Perovskites['Encapsulation_materials'].str.contains('SLG|Glass', case=False) == False]
    Foils= Non_Glass[Non_Glass['Encapsulation_materials'].str.contains('tape|plasma polymers|viewbarrier|Barrier foil|sheets|UHPBF|ITO|Kapton|PEN|PET|NOA', case= False)== True]
    polymer = Non_Glass[Non_Glass['Encapsulation_materials'].str.contains('PVP |Cyanoacrylate|PCL|Carbon-epoxy|Teflon|PDMS|Paraffin|LDPE|Meltronix|Pattex|Thermoplastic|poly|surlyn|EVA|Ethylene|Parylene|EVOH|PMMA',case=False) == True]
    Non_Glass = Non_Glass[~Non_Glass.index.isin(Foils.index)]
    No_info = Non_Glass[~Non_Glass.index.isin(polymer.index)]
    # finished classifying to groups

    info_for_analysis = pd.concat(
        [Glass_epoxy,Glass_UV_glue, Glass_polymer, Glass_butyl_rubber, Al2O3, Foils, polymer])
    all_info_for_analysis = info_for_analysis[['Ref_ID', 'Ref_DOI_number', 'Encapsulation_stack_sequence', 'Encapsulation_edge_sealing_materials', 'Stability_time_total_exposure',
                        'Stability_PCE_end_of_experiment', 'Stability_PCE_T95','Stability_PCE_Ts95','Stability_PCE_T80','Stability_PCE_Ts80','Stability_PCE_Te80',
                        'Stability_PCE_Tse80','Stability_PCE_after_1000_h','Stability_PCE_burn_in_observed','Stability_light_source_type','Stability_protocol',
                        'Stability_potential_bias_load_condition', 'Calculated_Ts80']]
    info_for_analysis.to_csv('experiments_with_Encap_and_stability.csv')

    Graph_titles= ["All 252 experiments","Ambient experiments","Experiments at 65 degrees","Dark Open Circuit Experiments","Light Experiments","Light MPPT Experiments"]
    #preparing data for the graphs X=encapsulation_group X[0]-> all experiments, X[1]-> ambient experiments, X[2]-> experiments in 65deg,
    # X[3]-> dark OC experiments, X[4]->light OC experiments, X[5]-> light MPP experiments
    Glass_butyl_rubber_Ts80=[]
    Glass_butyl_rubber_Ts80.append(Glass_butyl_rubber['Calculated_Ts80'].reset_index(drop=True))
    Glass_butyl_rubber_Ts80.append(Glass_butyl_rubber[Glass_butyl_rubber['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_butyl_rubber_Ts80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_protocol'].str.contains('2|3'))&
                               (~Glass_butyl_rubber['Stability_protocol'].str.contains('IEC'))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_butyl_rubber_Ts80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_light_source_type']=='Dark')|
                             (Glass_butyl_rubber['Stability_protocol'].str.contains('D',regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_butyl_rubber_Ts80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_light_source_type']!='Dark')&(~Glass_butyl_rubber['Stability_protocol'].str.contains('D',regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_butyl_rubber_Ts80.append(Glass_butyl_rubber[(Glass_butyl_rubber['Stability_potential_bias_load_condition'] == 'MPPT') |
                                            (Glass_butyl_rubber['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_polymer_Ts80=[]
    Glass_polymer_Ts80.append(Glass_polymer['Calculated_Ts80'].reset_index(drop=True))
    Glass_polymer_Ts80.append(Glass_polymer[Glass_polymer['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_polymer_Ts80.append(Glass_polymer[Glass_polymer['Stability_protocol'].str.contains('2|3')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_polymer_Ts80.append(Glass_polymer[(Glass_polymer['Stability_light_source_type'] == 'Dark') |
                                   (Glass_polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_polymer_Ts80.append(Glass_polymer[(Glass_polymer['Stability_light_source_type'] != 'Dark')&(~Glass_polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_polymer_Ts80.append(Glass_polymer[(Glass_polymer['Stability_potential_bias_load_condition'] == 'MPPT') |
                                      (Glass_polymer['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    Al2O3_Ts80=[]
    Al2O3_Ts80.append(Al2O3['Calculated_Ts80'].reset_index(drop=True))
    Al2O3_Ts80.append(Al2O3[Al2O3['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    Al2O3_Ts80.append(Al2O3[Al2O3['Stability_protocol'].str.contains('2|3')]['Calculated_Ts80'].reset_index(drop=True))
    Al2O3_Ts80.append(Al2O3[(Al2O3['Stability_light_source_type'] == 'Dark') |
                       (Al2O3['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Al2O3_Ts80.append(Al2O3[(Al2O3['Stability_light_source_type'] != 'Dark')&(~Al2O3['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Al2O3_Ts80.append(Al2O3[(Al2O3['Stability_potential_bias_load_condition'] == 'MPPT') |
                      (Al2O3['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    Foils_Ts80=[]
    Foils_Ts80.append(Foils['Calculated_Ts80'].reset_index(drop=True))
    Foils_Ts80.append(Foils[Foils['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    Foils_Ts80.append(Foils[Foils['Stability_protocol'].str.contains('2|3')]['Calculated_Ts80'].reset_index(drop=True))
    Foils_Ts80.append(Foils[(Foils['Stability_light_source_type'] == 'Dark') |
                       (Foils['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Foils_Ts80.append(Foils[(Foils['Stability_light_source_type'] != 'Dark')&(~Foils['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Foils_Ts80.append(Foils[(Foils['Stability_potential_bias_load_condition'] == 'MPPT') |
                      (Foils['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    polymer_Ts80=[]
    polymer_Ts80.append(polymer['Calculated_Ts80'].reset_index(drop=True))
    polymer_Ts80.append(polymer[polymer['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    polymer_Ts80.append(polymer[polymer['Stability_protocol'].str.contains('2|3')]['Calculated_Ts80'].reset_index(drop=True))
    polymer_Ts80.append(polymer[(polymer['Stability_light_source_type'] == 'Dark') |
                           (polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    polymer_Ts80.append(polymer[(polymer['Stability_light_source_type'] != 'Dark')&(~polymer['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    polymer_Ts80.append(polymer[(polymer['Stability_potential_bias_load_condition'] == 'MPPT') |
                          (polymer['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_UV_glue_Ts80=[]
    Glass_UV_glue_Ts80.append(Glass_UV_glue['Calculated_Ts80'].reset_index(drop=True))
    Glass_UV_glue_Ts80.append(Glass_UV_glue[Glass_UV_glue['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_UV_glue_Ts80.append(Glass_UV_glue[Glass_UV_glue['Stability_protocol'].str.contains('2|3')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_UV_glue_Ts80.append(Glass_UV_glue[(Glass_UV_glue['Stability_light_source_type'] == 'Dark') |
                                       (Glass_UV_glue['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_UV_glue_Ts80.append(Glass_UV_glue[(Glass_UV_glue['Stability_light_source_type'] != 'Dark')&(~Glass_UV_glue['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_UV_glue_Ts80.append(Glass_UV_glue[(Glass_UV_glue['Stability_potential_bias_load_condition'] == 'MPPT') |
                                      (Glass_UV_glue['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_epoxy_Ts80=[]
    Glass_epoxy_Ts80.append(Glass_epoxy['Calculated_Ts80'].reset_index(drop=True))
    Glass_epoxy_Ts80.append(Glass_epoxy[Glass_epoxy['Stability_protocol'].str.contains('1|Other')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_epoxy_Ts80.append(Glass_epoxy[Glass_epoxy['Stability_protocol'].str.contains('2|3')]['Calculated_Ts80'].reset_index(drop=True))
    Glass_epoxy_Ts80.append(Glass_epoxy[(Glass_epoxy['Stability_light_source_type'] == 'Dark') |
                                   (Glass_epoxy['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_epoxy_Ts80.append(Glass_epoxy[(Glass_epoxy['Stability_light_source_type'] != 'Dark')&(~Glass_epoxy['Stability_protocol'].str.contains('D', regex=False))]['Calculated_Ts80'].reset_index(drop=True))
    Glass_epoxy_Ts80.append(Glass_epoxy[(Glass_epoxy['Stability_potential_bias_load_condition'] == 'MPPT') |
                                  (Glass_epoxy['Stability_potential_bias_load_condition'] == 'Passive resistance')]['Calculated_Ts80'].reset_index(drop=True))
    table = []
    table.append(['', 'Foils','Polymer','Glass UV glue','Glass epoxy','Glass butyl rubber','Glass Polymer','Al2O3'])
    my_colors = {'Foils': 'teal', 'Polymers': 'steelblue', 'Glass+UV glue': 'slateblue', 'Glass+epoxy':'darkorchid','Glass+butyl rubber':'mediumvioletred' ,'Glass+polymer': 'sienna', 'Al2O3':'darkkhaki' }
    for i in range(6):
        data_for_box_plot = pd.concat([Foils_Ts80[i],polymer_Ts80[i], Glass_UV_glue_Ts80[i], Glass_epoxy_Ts80[i],Glass_butyl_rubber_Ts80[i],Glass_polymer_Ts80[i],Al2O3_Ts80[i]], axis =1)  # horizontal_concat
        data_for_box_plot=data_for_box_plot.set_axis(["Foils","Polymers","Glass+UV glue","Glass+epoxy","Glass+butyl rubber", "Glass+polymer", "Al2O3"], axis="columns")
        plt.rcParams.update({'font.size': 16})
        seaborn.catplot(x="variable", y="value", data=pd.melt(data_for_box_plot.dropna(how= "all", axis=1)), kind="box",palette=my_colors)
        seaborn.swarmplot(x='variable', y='value', data=pd.melt(data_for_box_plot.dropna(how= "all", axis=1)), color="gray")
        seaborn.pointplot(data=pd.melt(data_for_box_plot.dropna(how= "all", axis=1)), x="variable", y="value", estimator=np.mean,
                          join=False, ci=None, markers='+', color='black', zorder=3)
        plt.ylabel("Ts80")
        plt.yscale('log')
        plt.minorticks_on()
        plt.xticks(rotation=25)
        plt.xlabel("")
        plt.title(Graph_titles[i])

        table.append(['num of experiments'+Graph_titles[i],len(Foils_Ts80[i]), len(polymer_Ts80[i]), len(Glass_UV_glue_Ts80[i]), len(Glass_epoxy_Ts80[i]),len(Glass_butyl_rubber_Ts80[i]), len(Glass_polymer_Ts80[i]),
                      len(Al2O3_Ts80[i])])
        table.append(['average'+Graph_titles[i], np.average(Foils_Ts80[i]), np.average(polymer_Ts80[i]),np.average(Glass_UV_glue_Ts80[i]),
                      np.average(Glass_epoxy_Ts80[i]),np.average(Glass_butyl_rubber_Ts80[i]), np.average(Glass_polymer_Ts80[i]), np.average(Al2O3_Ts80[i])])
        table.append(['stdev'+Graph_titles[i], np.std(Foils_Ts80[i]),np.std(polymer_Ts80[i]), np.std(Glass_UV_glue_Ts80[i]), np.std(Glass_epoxy_Ts80[i]),np.std(Glass_butyl_rubber_Ts80[i]), np.std(Glass_polymer_Ts80[i]), np.std(Al2O3_Ts80[i])])
        t_test = []
        t_test.append([None,'Polymers','Glass+UV glue','Glass+epoxy','Glass +butyl rubber','Glass +polymer','Al2O3'])
        t_test.append(['Foils', scipy.stats.ttest_ind(Foils_Ts80[i], polymer_Ts80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(Foils_Ts80[i], Glass_UV_glue_Ts80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(Foils_Ts80[i], Glass_epoxy_Ts80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(Foils_Ts80[i], Glass_butyl_rubber_Ts80[i], axis=0, equal_var=False).pvalue,
                       scipy.stats.ttest_ind(Foils_Ts80[i], Glass_polymer_Ts80[i], axis=0, equal_var=False).pvalue,scipy.stats.ttest_ind(Foils_Ts80[i], Al2O3_Ts80[i], axis=0, equal_var=False).pvalue])
        t_test.append(['Polymers',None, scipy.stats.ttest_ind(polymer_Ts80[i], Glass_UV_glue_Ts80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(polymer_Ts80[i], Glass_epoxy_Ts80[i], axis=0, equal_var=False).pvalue, scipy.stats.ttest_ind(polymer_Ts80[i], Glass_butyl_rubber_Ts80[i], axis=0, equal_var=False).pvalue,
                       scipy.stats.ttest_ind(polymer_Ts80[i], Glass_polymer_Ts80[i], axis=0, equal_var=False).pvalue,scipy.stats.ttest_ind(polymer_Ts80[i], Al2O3_Ts80[i], axis=0, equal_var=False).pvalue])
        t_test.append(['Glass+UV_glue',None, None, scipy.stats.ttest_ind(Glass_UV_glue_Ts80[i], Glass_epoxy_Ts80[i], axis=0,equal_var=False).pvalue,
                           scipy.stats.ttest_ind(Glass_UV_glue_Ts80[i], Glass_butyl_rubber_Ts80[i], axis=0,equal_var=False).pvalue,
                           scipy.stats.ttest_ind(Glass_UV_glue_Ts80[i], Glass_polymer_Ts80[i], axis=0, equal_var=False).pvalue,
                           scipy.stats.ttest_ind(Glass_UV_glue_Ts80[i], Al2O3_Ts80[i], axis=0, equal_var=False).pvalue])
        t_test.append(['Glass+epoxy',None, None, None,
                                scipy.stats.ttest_ind(Glass_epoxy_Ts80[i], Glass_butyl_rubber_Ts80[i], axis=0,equal_var=False).pvalue,
                                scipy.stats.ttest_ind(Glass_epoxy_Ts80[i], Glass_polymer_Ts80[i], axis=0,equal_var=False).pvalue,
                                scipy.stats.ttest_ind(Glass_epoxy_Ts80[i], Al2O3_Ts80[i], axis=0,equal_var=False).pvalue])
        t_test.append(['Glass+butyl rubber',None, None, None, None,
                              scipy.stats.ttest_ind(Glass_butyl_rubber_Ts80[i], Glass_polymer_Ts80[i], axis=0,equal_var=False).pvalue,
                              scipy.stats.ttest_ind(Glass_butyl_rubber_Ts80[i], Al2O3_Ts80[i], axis=0,equal_var=False).pvalue])
        t_test.append(['Glass+polymer',None, None, None, None, None,
                                     scipy.stats.ttest_ind(Glass_polymer_Ts80[i], Al2O3_Ts80[i], axis=0,equal_var=False).pvalue])
        t_test = pd.DataFrame(t_test)
        t_test.to_csv('t_test_'+Graph_titles[i] +'_Ts80.csv')
    plt.show()
    table= pd.DataFrame(table)
    table.to_csv('frequency table Ts80.csv')
