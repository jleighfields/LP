
# Author: Justin Fields
# Date: 2020-09-03
# This script turns the jupyter note book 'LP_ortools.ipynb' 
# into a function that can be called from R.

def run_lp(diverse_profile
            ,peak_load
            ,restrict_gas
            ,min_charge_level
            ,init_ch_level
            ,batt_hours
            ,batt_eff
            ,use_outside_energy
            ,outside_energy_cost
            ,gas_mw_cost
            ,gas_mwh_cost
            ,batt_mw_cost
            ,wind_mw_cost
            ,wind_mwh_cost
            ,solar_mw_cost
            ,solar_mwh_cost
            ):


    """Run a linear program that minimizes costs with the constraints:
        1. load must be served
        2. hourly and monthly hydro contraints
        3. battery state of charge limits
    This optimization does not consider an outside market, it only minimizes
    costs to serve native load with native resources.

    Keyword arguments:
    profile_year -- the year to use for load and re generation
    restrict_gas -- the maximum amount of gas generation as a percent of load
    min_charge_level -- minimum charge level of batteries
    init_ch_level -- initial charge level of batteries
    batt_hours -- duration of hours for batteries
    batt_eff -- efficiency of batteries
    use_outside_energy -- use outside energy to meet load 
    outside_energy_cost -- cost of outside energy
    gas_mw_cost -- gas fixed cost $/MW including carbon costs
    gas_mwh_cost -- gas variable cost $/MWh including carbon costs
    batt_mw_cost -- battery fixed cost $/MW including carbon costs
    wind_mw_cost -- wind fixed costs $/MW including carbon costs
    wind_mwh_cost -- wind variable costs $/MWh including carbon costs
    solar_mw_cost -- solar fixed costs $/MW including carbon costs
    solar_mwh_cost -- solar variable costs $/MWh including carbon costs
    """

    import pandas as pd
    import numpy as np
    import time
    import ortools
    from ortools.linear_solver import pywraplp

    # set up some variables

    # boolean for debug printing
    debug_print = False
    
    # boolean for hydro 
    use_hydro = False

    # restrict gas to a portion of total load, 0-1 or None
    # e.g. 0.05 -> 5% limit on gas generation 
    # and  1.0 -> no limit on gas generation
    # divide by 100 since input is in percentages
    restrict_gas = restrict_gas/100


    # read profile data for load and re gen
    if diverse_profile == 'Yes':
        profile_data = 'lp_test_profiles_8760_2018_re_diverse.csv'
    else:
        profile_data = 'lp_test_profiles_8760_2018_re.csv'
        
    # read profiles
    df = pd.read_csv(profile_data, index_col = 'Hour')
    
    # scale load to peak_load input
    df['2030_load'] = df['2030_load'] * peak_load
    
    print('Load and generation data:')
    print(df.head().T)
    print('\n')

    print('Summary of load and generation data:')
    print(df.describe().T)
    print('\n')

    # read hydro parameters
    if use_hydro:
        hydro = pd.read_csv('lp_test_hydro_data.csv')
    
        print('Hydro data:')
        print(hydro.head().T)
        print('\n')


    # start timer
    total_time_0 = time.time()


    # Create the linear solver with the GLOP backend.
    solver = pywraplp.Solver('simple_lp_program', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)


    # build capacity decision variables for the resources
    batt = solver.NumVar(100, solver.infinity(), 'batt')
    solar = solver.NumVar(0, solver.infinity(), 'solar')
    wind = solver.NumVar(0, solver.infinity(), 'wind')
    gas = solver.NumVar(0, solver.infinity(), 'gas')


    # use this to fix the amount of solar and wind and solve for battery MW
    #solver.Add(solar <= 900)
    #solver.Add(solar >= 900)
    #solver.Add(wind <= 850)
    #solver.Add(wind >= 850)

    print('Adding variables for build capacity')
    print('Number of variables =', solver.NumVariables(), '\n')

    print('Adding hourly variables and constraints')
    # generation decision variables
    t0 = time.time()

    # create arrays to hold hourly varaibles
    # we don't need solar and wind since the profiles are fixed
    batt_ch         = [None] * len(df.index)
    batt_disch      = [None] * len(df.index)
    gas_gen         = [None] * len(df.index)
    if use_hydro:
        hydro_gen_crsp  = [None] * len(df.index)
        hydro_gen_lap   = [None] * len(df.index)
    SOC             = [None] * len(df.index)
    if use_outside_energy: outside_energy  = [None] * len(df.index)


    t_h = time.time() # for tracking time in hourly loop
    for h in df.index:
        
        if h % 100 == 0 and debug_print:
                print("h: ", h, "\tET: ", round(time.time()-t_h, 2))
                t_h = time.time()
                
        # add hourly variables
        batt_ch[h] = solver.NumVar(0, solver.infinity(), 'batt_ch[{}]'.format(h))
        batt_disch[h] = solver.NumVar(0, solver.infinity(), 'batt_disch[{}]'.format(h))
        gas_gen[h] = solver.NumVar(0, solver.infinity(), 'gas_gen[{}]'.format(h))
        if use_hydro:
            hydro_gen_crsp[h] = solver.NumVar(0, solver.infinity(), 'hydro_gen_crsp[{}]'.format(h))
            hydro_gen_lap[h] = solver.NumVar(0, solver.infinity(), 'hydro_gen_lap[{}]'.format(h))
        SOC[h] = solver.NumVar(0, solver.infinity(), 'SOC[{}]'.format(h))
        if use_outside_energy: outside_energy[h] = solver.NumVar(0, solver.infinity(), 'outside_energy[{}]'.format(h))
        
        
        # add hourly constraints
        
        # set SOC[h] equal to previous hour SOC 
        # plus the change from charging or discharging
        if h == 0:
            solver.Add(SOC[h] <= init_ch_level*4*batt + batt_ch[h] - batt_disch[h]/batt_eff)
            solver.Add(SOC[h] >= init_ch_level*4*batt + batt_ch[h] - batt_disch[h]/batt_eff)
        else:
            solver.Add(SOC[h] <= SOC[h-1] + batt_ch[h] - batt_disch[h]/batt_eff)
            solver.Add(SOC[h] >= SOC[h-1] + batt_ch[h] - batt_disch[h]/batt_eff)
            

        
        # SOC mwh constraints
        # max mwh constraint
        solver.Add(SOC[h] <= batt_hours*batt)
        
        # min mwh constraint
        solver.Add(SOC[h] >= min_charge_level*batt_hours*batt)
        
        # get battery capacity
        #solver.Add(batt_cap[h] >= batt_hours*batt - SOC[h])
        #solver.Add(batt_cap[h] <= batt_hours*batt - SOC[h])
        
        # hourly demand constraints
        
        if use_outside_energy and use_hydro: 
            solver.Add(hydro_gen_crsp[h] + hydro_gen_lap[h] +
                    (df.solar[h]*solar) +    # based on df hourly profile
                    (df.wind[h]*wind) +      # based on df hourly profile
                    gas_gen[h] + 
                    batt_disch[h] - batt_ch[h] +
                    outside_energy[h]
                    >= df['2030_load'][h])
        elif use_outside_energy: 
            solver.Add(hydro_gen_crsp[h] + hydro_gen_lap[h] +
                    (df.solar[h]*solar) +    # based on df hourly profile
                    (df.wind[h]*wind) +      # based on df hourly profile
                    gas_gen[h] + 
                    batt_disch[h] - batt_ch[h] +
                    outside_energy[h]
                    >= df['2030_load'][h])
        elif use_hydro:
            solver.Add(hydro_gen_crsp[h] + hydro_gen_lap[h] +
            (df.solar[h]*solar) +    # based on df hourly profile
            (df.wind[h]*wind) +      # based on df hourly profile
            gas_gen[h] + 
            batt_disch[h] - batt_ch[h]
            >= df['2030_load'][h])
        else:
            solver.Add((df.solar[h]*solar) +    # based on df hourly profile
            (df.wind[h]*wind) +      # based on df hourly profile
            gas_gen[h] + 
            batt_disch[h] - batt_ch[h]
            >= df['2030_load'][h])
        
        
        # hourly generation constraints
        solver.Add(batt_ch[h] <= batt) 
        solver.Add(batt_disch[h] <= batt)
        solver.Add(gas_gen[h] <= gas)
        
        # hydro hourly constraints by month
        if use_hydro:
            solver.Add(hydro_gen_crsp[h] >= df.MoCRSPMin[h])
            solver.Add(hydro_gen_crsp[h] <= df.MoCRSPMax[h])
        
            solver.Add(hydro_gen_lap[h] >= df.MoLAPMin[h])
            solver.Add(hydro_gen_lap[h] <= df.MoLAPMax[h])
        


    # total gas gen constraint    
    if restrict_gas != None:
        solver.Add(solver.Sum(gas_gen) <= restrict_gas*sum(df['2030_load'])) 
        
        
    # hydro constraints by month
    if use_hydro:
        for i in range(1,13):
            # sum energy from h0 to h1 to get monthly energy
            h0 = int(hydro.HourMin[hydro.Months == i])
            h1 = int(hydro.HourMax[hydro.Months == i] + 1)
            
            # must take all hydro energy
            solver.Add(solver.Sum(hydro_gen_crsp[h0:h1]) <= float(hydro.CRSPMWh[hydro.Months == i]))
            solver.Add(solver.Sum(hydro_gen_crsp[h0:h1]) >= float(hydro.CRSPMWh[hydro.Months == i]))
            
            solver.Add(solver.Sum(hydro_gen_lap[h0:h1]) <= float(hydro.LAPMWh[hydro.Months == i]))
            solver.Add(solver.Sum(hydro_gen_lap[h0:h1]) >= float(hydro.LAPMWh[hydro.Months == i]))

            
    t1 = time.time()

    print('time to build model (seconds): {0:,.2f}\n'.format((t1-t0),1))
    print('Number of variables: {0:,}'.format(solver.NumVariables()))
    print('Number of contraints: {0:,}'.format(solver.NumConstraints()))
    print('\n')


    # look at the coeficients
    solar_coef = (solar_mw_cost + solar_mwh_cost*sum(df.solar))
    wind_coef = (wind_mw_cost + wind_mwh_cost*sum(df.wind))
    names = ['solar_coef', 'wind_coef', 'batt_mw_coef', 'gas_mw_coef', 'gas_mwh_coef']
    coefs = np.round([solar_coef, wind_coef, batt_mw_cost, gas_mw_cost, gas_mwh_cost],0)

    print('Coeficients for costs function:')
    for i in range(len(names)):
        print(names[i] +': \t'+str(coefs[i]))

    print('\n')


    # Objective Function

    objective = solver.Objective()

    # get capacity and energy costs coefficient for wind and solar
    # sum profile and multiply by variable costs to get the energy component of the coefficient
    solar_coef = (solar_mw_cost + solar_mwh_cost*sum(df.solar))
    wind_coef = (wind_mw_cost + wind_mwh_cost*sum(df.wind))

    # set the coefficients in the objective function for the capacity variables
    # noting that the solar and wind coefficients include the variable costs
    objective.SetCoefficient(batt, batt_mw_cost)
    objective.SetCoefficient(solar, solar_coef)
    objective.SetCoefficient(wind, wind_coef)
    objective.SetCoefficient(gas, gas_mw_cost)

    # add energy costs
    for h in df.index:
        objective.SetCoefficient(gas_gen[h], gas_mwh_cost)
        if use_outside_energy: objective.SetCoefficient(outside_energy[h], outside_energy_cost)

        # disincentivize charging and discharging at the same time
        # this removes hours that both charge and discharge
        objective.SetCoefficient(batt_disch[h], 0.02)
        #objective.SetCoefficient(batt_ch[h], -0.01)

        # benefit to keeping the batteries charged
        objective.SetCoefficient(SOC[h], -0.01)


        # minimize the cost to serve the system
        objective.SetMinimization()

    # solve the system
    print('Starting optimization...')
    t0 = time.time()
    status = solver.Solve()
    t1 = time.time()

    print('time to solve (minutes): {0:,.2f}\n'.format((t1-t0)/60,1))

    print('Solution is optimal: ', status == solver.OPTIMAL, '\n')

    obj_val = objective.Value()
    print('Solution:')
    print('Objective value = {0:,.0f}\n'.format(obj_val))


    print('Build variables:')
    batt_mw = batt.solution_value()
    solar_mw = solar.solution_value()
    wind_mw = wind.solution_value()
    gas_mw = gas.solution_value()
    cap_mw = {'batt_mw':batt_mw,
            'solar_mw':solar_mw,
            'wind_mw':wind_mw,
            'gas_mw':gas_mw}


    # print results for build variables
    for r in [batt, solar, wind, gas]:
        print('{0} \t= {1:,.0f}'.format(r, r.solution_value()))

    print('\n')


    # get the solved values to check the solution
    # create a new data frame to hold the final solution values
    print('Gathering hourly data...\n')
    final_df = df[['Date', '2030_load', 'solar', 'wind']].copy()
    final_df['solar'] = final_df['solar'] * solar.solution_value()
    final_df['wind'] = final_df['wind'] * wind.solution_value()
    final_df['gas'] = 0
    final_df['batt_charge'] = 0
    final_df['batt_discharge'] = 0
    final_df['SOC'] = 0
    if use_hydro:
        final_df['crsp'] = 0
        final_df['lap'] = 0
    if use_outside_energy: final_df['outside_energy'] = 0


    # get the battery charge and discharge by hour
    batt_ch_hourly = [None] * len(df.index)
    batt_disch_hourly = [None] * len(df.index)
    for h in df.index:
        batt_ch_hourly[h] = batt_ch[h].solution_value()
        batt_disch_hourly[h] = batt_disch[h].solution_value()
        

    # get cumulative sums for calculating SOC
    batt_ch_hourly = np.cumsum(batt_ch_hourly)
    batt_disch_hourly = np.cumsum(batt_disch_hourly)


    # get the hourly data
    for h in final_df.index:
        if h % 100 == 0 and debug_print:
            print("h: ", h)
            
        final_df.loc[h, 'gas'] = gas_gen[h].solution_value()
        final_df.loc[h, 'batt_charge'] = batt_ch[h].solution_value()
        final_df.loc[h, 'batt_discharge'] = batt_disch[h].solution_value()
        final_df.loc[h, 'SOC'] = SOC[h].solution_value()
        if use_hydro:
            final_df.loc[h, 'crsp'] = hydro_gen_crsp[h].solution_value()
            final_df.loc[h, 'lap'] = hydro_gen_lap[h].solution_value()
        if use_outside_energy: final_df.loc[h, 'outside_energy'] = outside_energy[h].solution_value()
        

    # calc net load for a check on the results
    if use_outside_energy and use_hydro: 
        final_df['net_load'] = round((final_df['crsp'] +
                                    final_df['lap'] +
                                    final_df['solar'] + 
                                    final_df['wind'] + 
                                    final_df['gas'] +
                                    final_df['batt_discharge'] -
                                    final_df['batt_charge'] +
                                    final_df['outside_energy'] - 
                                    final_df['2030_load']), 2)
        
    elif use_hydro:
        final_df['net_load'] = round((final_df['crsp'] +
                                    final_df['lap'] +
                                    final_df['solar'] + 
                                    final_df['wind'] + 
                                    final_df['gas'] +
                                    final_df['batt_discharge'] -
                                    final_df['batt_charge'] -
                                    final_df['2030_load']), 2)
    elif use_outside_energy:
        final_df['net_load'] = round((final_df['solar'] + 
                                    final_df['wind'] + 
                                    final_df['gas'] +
                                    final_df['batt_discharge'] -
                                    final_df['batt_charge'] +
                                    final_df['outside_energy'] -
                                    final_df['2030_load']), 2)
    else:
        final_df['net_load'] = round((final_df['solar'] + 
                                    final_df['wind'] + 
                                    final_df['gas'] +
                                    final_df['batt_discharge'] -
                                    final_df['batt_charge'] -
                                    final_df['2030_load']), 2)
        

    final_df['load_and_charge'] = round((final_df['batt_charge'] +
                                        final_df['2030_load']), 2)


    # set the index to hours in 2030
    final_df.set_index(
        pd.date_range(start='2030-01-01 01:00:00', periods=final_df.shape[0], freq='h'),
        inplace=True
    )


    # take a look
    print(final_df.head(5).T)
    print('\n')


    # summarize the data
    print('Summary of hourly data:\n')
    print(final_df.describe().T)
    print('\n')


    # this should be empty...
    print('Any negative net load? Should be empty...')
    print(final_df[(final_df.net_load < 0)].T)
    print('\n')


    # this should be empty...
    print('Any hours with both charging and discharging? Should be empty...')
    print(final_df[(final_df.batt_discharge > 0) & (final_df.batt_charge > 0)].T)
    print('\n')


    if use_outside_energy:
        outside_energy_percent = 100*final_df.outside_energy.sum()/final_df['2030_load'].sum()
        print('Outside energy as a percentage of load: {0:,.3f}%\n'.format(outside_energy_percent))


    gas_percent = 100*final_df.gas.sum()/final_df['2030_load'].sum()
    print('Gas generation as a percentage of load: {0:,.2f}%\n'.format(gas_percent))


    re_percent = 100*((final_df.solar.sum()+final_df.wind.sum())/final_df['2030_load'].sum())
    print('RE generation as a percentage of load: {0:,.2f}%\n'.format(re_percent))


    excess_percent = 100*(final_df.net_load.sum()/final_df['2030_load'].sum())
    print('Excess generation as a percentage of load: {0:,.2f}%\n'.format(excess_percent))


    excess_re_percent = 100*final_df.batt_discharge.sum()/final_df.batt_charge.sum()
    print('Batt discharge as a percentage of batt charge: {0:,.2f}%\n'.format(excess_re_percent))


    total_time_1 = time.time()


    print('total time to build, solve, and verify (minutes): {0:,.2f}\n'.format((total_time_1-total_time_0)/60))


    # return dictionary for displaying results
    return({'obj_val':obj_val, 'cap_mw':cap_mw, 'final_df':final_df})


# for testing
# note: cannot use __main__ in dashboard use test_function instead
#if __name__ == "__main__":

test_function = False
if test_function:
    print('\n')

    # set up inputs for optimization

    # select year for profiles
    profile_year = '2018'

    use_outside_energy = False
    outside_energy_cost = 10000

    # restrict gas to a portion of total load, 0-1 or None
    # e.g. 0.05 -> 5% limit on gas generation 
    # and  1.0 -> no limit on gas generation
    restrict_gas = 100

    # battery parameters
    min_charge_level = 0.1
    init_ch_level = 0.5
    batt_hours = 4
    batt_eff = 0.85


    # cost of CO2
    # C02 values from CSU study
    gas_co2_ton_per_mwh = (411 + 854)/2000
    # assumed 20 year life
    gas_co2_ton_per_mw = (5981+1000+35566+8210+10165+1425)/(6*18)/20

    wind_co2_ton_per_mwh = 0.2/2000
    # assumed 20 year life
    wind_co2_ton_per_mw = (754+10-241)/20

    solar_co2_ton_per_mwh = 2.1/2000
    # assumed 20 year life
    solar_co2_ton_per_mw = (1202+250-46)/20

    # battery C02 given in lbs
    # assumed 15 year life
    batt_co2_ton_per_mw = (1940400-83481+4903)/2000/15


    # Carbon 2030 $/ton = $9.06 
    #co2_cost = 9.06
    co2_cost = 160

    # calculate costs in $/MWh and $/MW 
    # these costs include the cost of carbon

    cc_gas_mwh = co2_cost * gas_co2_ton_per_mwh
    cc_gas_mw = co2_cost * gas_co2_ton_per_mw

    cc_wind_mwh = co2_cost * wind_co2_ton_per_mwh
    cc_wind_mw = co2_cost * wind_co2_ton_per_mw

    cc_solar_mwh = co2_cost * solar_co2_ton_per_mwh
    cc_solar_mw = co2_cost * solar_co2_ton_per_mw

    cc_batt_mw = co2_cost * batt_co2_ton_per_mw


    # capacity cost in $/kw-mo
    gas_cap_cost = 11.27
    gas_fixed_cost = cc_gas_mw + gas_cap_cost*12*1000 # converted to $/MW-yr

    heat_rate = 8883 # btu/kwh
    vom = 7.16 # $/mwh
    gas_fuel_cost = 4.37 # $/mmbtu
    gas_variable_cost = cc_gas_mwh + gas_fuel_cost*heat_rate/1000 + vom

    batt_cost = 8.25
    batt_fixed_cost = cc_batt_mw + batt_cost*12*1000  # converted to $/MW-yr

    wind_kw_mo = 1.13
    wind_fixed_cost = cc_wind_mw + wind_kw_mo*12*1000  # converted to $/MW-yr
    wind_variable_cost = cc_wind_mwh + 41.01 # $/mwh


    solar_kw_mo = 1.13
    solar_fixed_cost = cc_solar_mw + solar_kw_mo*12*1000  # converted to $/MW-yr
    solar_variable_cost = cc_solar_mwh + 33.51 # $/mwh

    run_lp(profile_year
                    ,restrict_gas
                    ,min_charge_level
                    ,init_ch_level
                    ,batt_hours
                    ,batt_eff
                    ,use_outside_energy
                    ,outside_energy_cost
                    ,gas_fixed_cost
                    ,gas_variable_cost
                    ,batt_fixed_cost
                    ,wind_fixed_cost
                    ,wind_variable_cost
                    ,solar_fixed_cost
                    ,solar_variable_cost
    )

    print('Finished')





