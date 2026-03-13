# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 08:42:15 2024

@author: 20192020
"""
import Parameters as P
from Parameters import mesa, random, np, deque, math
from TumorCell import TumorCellAgent
from Lymphocyte import LymphocyteAgent 
from SimResults import SimResults

import numpy as np 

class CancerModel(mesa.Model):
    """Class to initialize ABM and run simulations"""
    
    # save space for objects that will not be changed 
    __slots__ = ("pars",  "grid", "N_tumor_died", "N_tumor_killed",
                 "N_lymp_died", "N_lymp_influx", "N_lymp_exhaust", 
                 "N_tumor", "N_tumor_stem", "N_lymp", "N_tot_tumor", 
                 "N_tot_lymp", "kill_agents", "schedule", "nSteps", "verbose",
                 "therapy_lymph", "therapy_tumor", "therapy_system", "sim_results",
                 "adj_par_tumor", "adj_par_lymph", "adj_par_system")
    
    def __init__(self, N, width, height, pars, verbose = False,
                 init_config = None):
        """
        Parameters
        ----------
        N : dictionary
            Number of cells to randomly place in configuration if no initial 
            configuration is given 
        width, height : float
            Width and height of the grid for the model simulations 
        pars : dictionary 
            Parameters and corresponding values 
        verbose : Boolean, optional
            If True, more information is printed for each action and step done by the 
            model. The default is False.
        init_config : None or DataFrame, optional
            An initial configuration can be given as DataFrame, else a random 
            initial configuration will be used. The default is None.
        """
        super().__init__()
        
        # Get parameter values
        self.pars = pars
        
        # Prepare simulation grid 
        self.grid = mesa.space.SingleGrid(width, height, False) # can only have one cell at one grid cell 
        
        # Variables that are kept for each step 
        self.N_tumor_died = 0 
        self.N_tumor_killed = 0 
        self.N_lymp_died = 0 
        self.N_lymp_influx = 0 
        self.N_lymp_exhaust = 0 
        
        self.N_tumor = 0 # keeps track of number of tumor cells in model
        self.N_tumor_stem = 0 # keeps track of number of tumor cells that are stem cells 
        self.N_lymp = 0 # keeps track of number of lymphocytes in model
        
        # Counter to know how many tumor / lymphocytes have been present in the simulation
        self.N_tot_tumor = 0
        self.N_tot_lymp = 0 
        
        # Initialize additional parameters 
        self.kill_agents = deque()
        
        # Create scheduler and assign it to the model
        # agents are handled randomly in each step of the model 
        self.schedule = mesa.time.RandomActivation(self)
        
        # Other (practical) properties
        self.nSteps = 0 # number of steps this model has performed
        self.verbose = verbose
        self.therapy_lymph = False #whether therapy is active that affects lymphocytes
        self.therapy_tumor = False #whether therapy is active that affect tumors
        self.therapy_system = False #whether therapy is active that affect general system properties 
        
        # list of parameters that need to be adjusted in this step 
        self.adj_par_tumor = []
        self.adj_par_lymph = []
        self.adj_par_system = []
        
        # prepare class to save results
        self.sim_results = SimResults()
        
        # Create random initial configuration if no initial seeding is provided
        if init_config is None:
            # Retrieve properties from given parameters 
            self.N_tumor = N['tumor'] 
            self.N_lymp = N['lymphocytes']  
            self.N_tot_tumor = N['tumor']
            self.N_tot_lymph = N['lymphocytes']
            
            # Create tumor cell agents 
            for i in range(self.N_tumor):
                # probability tumor cell is stem cell 
                is_stem = np.random.binomial(1, P.TUstem)
                
                self.add_tumor_cell(i, is_stem, verbose, rand_loc = True)
                
            # Create lymphocyte agents 
            for i in range(self.N_lymp):
                # probability tumor cell is stem cell 
                self.add_lymphocyte_cell(i, verbose, rand_loc = True)
                
        else: 
            self.import_init_configuration(init_config)
                
        # first snapshot
        self.record_metrics()

    def step(self):
        """Advance the model by one step."""
        # check if therapy has been administered at this step 
        if P.therapy_administration is not None:
            dummy = P.therapy_administration.copy()
            
            # check if one of given therapy needs to be adminstered
            while dummy.count(self.nSteps) > 0:
                current_idx = dummy.index(self.nSteps)
                
                # check in which class changes need to be made 
                if 'tumor' in P.therapy_effect[current_idx]:
                    self.therapy_tumor = True 
                    self.adj_par_tumor.append(P.therapy_par[current_idx])
                elif 'lymphocyte' in P.therapy_effect[current_idx]:
                    self.therapy_lymph = True
                    self.adj_par_lymph.append(P.therapy_par[current_idx])
                elif 'system' in P.therapy_effect[current_idx]:
                    self.therapy_system = True 
                    self.adj_par_system.append(P.therapy_par[current_idx])
                        
                # change value so if there is another therapy at this time step
                # it will also be simulated 
                dummy[current_idx] = dummy[current_idx] + 1                   
                
        # if a therapy is administered that affects a parameters linked to 'system'
        if self.therapy_system:
            for adj_par in self.adj_par_system:
                for (key, value) in adj_par.items():
                    # check if parameter is classified as system parameter
                    if key in P.system_pars: 
                        # change parameter value 
                        old = self.pars[key]
                        
                        if 'Prob' in key:
                            self.pars[key] = min(self.pars[key] * value, 1) 
                        else:
                            self.pars[key] = self.pars[key] * value
                        
                        if self.verbose > 0:
                            print(f'>> A system parameter was affected by therapy! {key} changes from {old} to {self.pars[key]}')
                        
        # remove cells that died 
        self.clear_dying_cells()
        
        # influx of lymphocytes 
        self.lymphocytes_influx()

        # The model's step will go here for now 
        # this will call the step method of each agent and print the agent's unique_id
        self.schedule.step()
        self.nSteps += 1 
        
        # save some properties at this time point 
        self.record_metrics()
        
        # update model parameters so also new T cells / tumor cells etc. are initialized with adjusted values
        # only need to update once for each therapy 
        for i in range(0, len(self.adj_par_tumor)):
            for (key, item) in self.adj_par_tumor[i].items():
                if self.verbose > 0:
                    print('>> Treatment has been administered, tumor cell parameters will be updated')
                
                if 'p' in key:
                    self.pars[key] = min(self.pars[key] * item, 1)
                else:
                    self.pars[key] = self.pars[key] * item
                    
        for i in range(0, len(self.adj_par_lymph)):  
            for (key, item) in self.adj_par_lymph[i].items():
                if self.verbose > 0:
                    print('>> Treatment has been administered, lymphocyte parameters will be updated')
                
                if 'p' in key:
                    self.pars[key] = min(self.pars[key] * item, 1)
                else:
                    self.pars[key] = self.pars[key] * item
        
        # reset some variables 
        self.N_lymp_died = 0 
        self.N_tumor_killed = 0
        self.N_tumor_died = 0 
        self.N_lymp_exhaust = 0
        # note : self.N_lymp_influx gets reset in self.lymphocytes_influx
        
        self.therapy_system = False 
        self.adj_par_system = []
        self.therapy_tumor = False 
        self.adj_par_tumor = []
        self.therapy_lymph = False 
        self.adj_par_lymph = []
        
    def record_metrics(self):
        """Save statistics of current status of the model"""
        self.sim_results.register_lymphocyte_died(self.N_lymp_died)
        self.sim_results.register_tumor_cell_killed(self.N_tumor_killed)
        self.sim_results.register_tumor_cell_died(self.N_tumor_died)
        self.sim_results.register_lymphocyte_influx(self.N_lymp_influx)
        self.sim_results.register_Ncells(self.N_tumor, self.N_lymp)
        self.sim_results.register_lymphocyte_exhausted(self.N_lymp_exhaust)
        self.sim_results.register_tumor_stem(self.N_tumor_stem)
        self.sim_results.register_tumor_total(self.N_tot_tumor)
        self.sim_results.register_lymp_total(self.N_tot_lymp)
        self.sim_results.register_new_round()
        
        
    def clear_dying_cells(self):
        """Remove cells that died"""        
        for _ in range(len(self.kill_agents)):
            # get the fallen agent and remove from list 
            died_agent = self.kill_agents.pop() 
            
            # update counter 
            if died_agent.type == 'tumor':
                self.N_tumor -= 1
                
                if died_agent.is_stem == 1:
                    self.N_tumor_stem -= 1
                
            else:
                self.N_lymp -= 1 
            
            self.schedule.remove(died_agent)
            self.grid.remove_agent(died_agent)
            
            
    def add_tumor_cell(self, i, is_stem, verbose,loc = [], rand_loc=True):
        """
        Add new cell if tumor cell proliferates 

        Parameters
        ----------
        i : int 
            unique_id for tumor cell 
        is_stem : Boolean 
            whether tumor cell will be a stem cell (0 = False, 1 = True)
        loc : list, optional
            x and y coordinates for location of new cell. The default is [].
        rand_loc : Boolean, optional
            whether the new cell is placed randomly or according to given loc. 
            The default is True
        """

        a = TumorCellAgent(i, self, is_stem, verbose) # create agent 
        
        # check if it's a stem cell 
        if is_stem == 1:
            self.N_tumor_stem += 1 
        
        # add agent to the scheduler 
        self.schedule.add(a)
        
        # if agent should be placed at random location 
        if rand_loc : 
            # Add the agent to a random empty grid cell 
            new_loc = random.choice(list(self.grid.empties.copy()))
            self.grid.place_agent(a, new_loc)
        # if agent should be placed at specific location 
        else : 
            self.grid.place_agent(a, (loc[0], loc[1]))
            
            
    def add_lymphocyte_cell(self, i, verbose, loc = [], rand_loc=True, lst_loc=None):
        """
        Add new cell if lymphocyte proliferates 

        Parameters
        ----------
        i : int 
            unique_id for lymphocyte cell 
        loc : list, optional
            x and y coordinates for location of new cell. The default is [].
        rand_loc : Boolean, optional
            whether the new cell is placed randomly or according to given loc. 
            The default is True
        lst_loc : list, optional
            if only a selected few positions can be chosen, a list can be given
        """

        a = LymphocyteAgent(i, self, verbose) # create agent 
        
        # add agent to the scheduler 
        self.schedule.add(a)
        
        # if agent should be placed at random location 
        if rand_loc and not P.IMinfluxEdge: 
            # Add the agent to a random empty grid cell 
            new_loc = random.choice(list(self.grid.empties.copy()))
            self.grid.place_agent(a, new_loc)
        # if agent should be placed at edge 
        elif P.IMinfluxEdge and lst_loc is not None: 
            new_loc = random.choice(lst_loc)
            self.grid.place_agent(a, new_loc)
        # if agent is placed at specific location 
        else : 
            self.grid.place_agent(a, (loc[0], loc[1]))
            new_loc = (loc[0], loc[1])
            
        # register new lymphocyte with corresponding position 
        self.sim_results.register_lymph_movement('LYM'+str(i), new_loc)
        
        return new_loc
            
    def lymphocytes_influx(self):
        """Influx of new lymphocytes"""
        # currently assumption that lymphocytes can appear at any place 
        if random.random() < self.pars['IMinfluxProb']:
            if self.verbose:
                print("We got an influx of lymphocytes!")
            
            # scale number of lymphocytes influx with a factor * number of tumor cells
            n_new_cells = self.pars['IMinfluxRate'] * self.pars['IMrateDynamic'] * self.N_tumor
            
            # check if there is still enough space for all the new cells 
            # else take the minimum 
            n_empty_spaces = len(list(self.grid.empties.copy()))
            possible_pos = None
            
            # adjust possible spaces if lymphocyte can only entire through edge
            if P.IMinfluxEdge:
                possible_pos = list(self.grid.empties.copy())
                possible_pos = [x for x in possible_pos if x[0] == 0 or x[1] == 0 or x[0] == (P.height-1) or x[1] == (P.height-1)]
                n_empty_spaces = len(possible_pos)
            
            if n_new_cells > n_empty_spaces and self.verbose:
                print(f"We got more influx than empty spaces! Influx adjusted from {n_new_cells} to {n_empty_spaces}")
            n_new_cells = min(n_new_cells, n_empty_spaces)
            for _ in range(int(math.ceil(n_new_cells))):
                new_loc = self.add_lymphocyte_cell(self.N_tot_lymp, self.verbose, 
                                         lst_loc=possible_pos)
                
                # remove location of newly added lymphocyte from list of candidates
                if P.IMinfluxEdge:
                    possible_pos.remove(new_loc)
                
                self.N_lymp += 1 
                self.N_tot_lymp += 1 
                
            self.N_lymp_influx += n_new_cells
            
        else:
            self.N_lymp_influx = 0 
            
    
    def import_init_configuration(self, init_config):
        """Use given init_config to create the initial configuration for the model"""
        
        if self.verbose > 0:
            print('>> Loading initial configuration')
        
        for index, row in init_config.iterrows():
            # check cell type, then add cell to simulation 
            if row['cell'] == 'tumor':
                self.add_tumor_cell(self.N_tot_tumor, 
                                    row['stem'], self.verbose, 
                                    loc = [row['x'], row['y']],
                                           rand_loc = False)
                
                self.N_tumor += 1
                self.N_tot_tumor += 1 
                
            elif row['cell'] == 'lymphocyte':
                self.add_lymphocyte_cell(self.N_tot_lymp,
                                         self.verbose, 
                                         loc = [row['x'], row['y']],
                                         rand_loc = False)
                
                self.N_lymp += 1 
                self.N_tot_lymp += 1 
                
    def export_init_configuration(self, filepath):
        """Write initial configuration dataframe to CSV"""
        self.init_config_df.to_csv(filepath)
                   
            
            
            
            
            
                