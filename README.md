# Linear Program Overview

This [dashboard](https://jleighfields.shinyapps.io/LP_capacity_model/) tests using OR-Tools to deploy a linear program that selects the optimal combination of the generation resources to serve electric load at the minimal cost.  The model is built for the year 2030 since many utilities have carbon reduction goals for the year 2030.  The linear program logic is developed in this [notebook](https://github.com/jleighfields/LP/blob/master/src/LP_ortools.ipynb) located in `src` folder.  The model is deployed in a dashboard built with `Shiny` and the `Flexdashboard` packages.  This allows the user to run sensitivity tests for various assumptions, including resource and carbon costs.  These cost inputs are collected through the dashboard interface and passed to the OR-Tools code via the `reticulate` package.  These inputs are set by the user on the `Optimization inputs and plots` tab.  After setting the inputs clicking the `Run optimization` will start solving the model run using the input values from the user.  The model run will take about a minute to solve.  Once a solution has been found the value boxes will be populated and the `LP results plot` will be populated.  The plot shows 7 days of hourly data and the plot start date can be set by the user.  All of the hourly results are available on the `Output data` tab.   

## Costs

Costs are allocated in terms of \$/MWh and \$/MW for all resources.  The \$/MWh can be thought of as a variable cost and the \$/MW as a fixed cost.  The life cycle carbon costs are added to these components and are based on a life cycle carbon study performed by Colorado State University.  Generally, carbon emissions associated with buidling, O&M, and decommissioning are included in the \$/MW costs and the operational emissions are included in the \$/MWh costs.

## Resources

Four resources are considered for serving load:

* Wind generation
* Solar generation
* Batteries
* Natural gas generation

Two generation patterns are available for the wind and solar resources, one pattern is from a single location while a second pattern represents a diversified portfolio of wind and solar generation from various locations.  This helps to exemplify the benefits of locational diversity to reduce periods of no generation and excess generation.  The MW capacities for wind and solar units are multiplied by the generation profiles to get the energy production from these resources.

Batteries are modeled as 4-hour lithium-ion batteries.

Gas generation units are modeled as completely flexible units that do not have minimum run times.  This simplifies the analysis and allows the problem to be structured as a linear program.  More complex models can be modeled as mixed-integer programs but they generally require a commercial solver which complicates setting up this model as an interactive application.  As such, Reciprocating Internal Combustion Engines (RICE) are used as the basis for the gas units since they are extremely flexible.

Outside markets are not explicitly considered in this model.  The use of outside energy can be set by the `use_outside_energy` flag and has an associated cost defined by `outside_energy_cost`.  These variables are intended to capture times when emergency may be needed to serve load.  If `use_outside_energy` is set to `False` then load will always be served from the resources listed above.
