# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 11:50:32 2024

@author: 20192020
"""
import Parameters as P
from Parameters import mesa, random, deque

class TumorCellAgent(mesa.Agent):
    """Agent to describe behavior of a tumor cell"""
    
    __slots__ = ("is_stem", "model", "unique_id", "verbose", "type", "n_prol", "pars", 
                 "damage", "engaged")
    
    def __init__(self, unique_id, model, is_stem, verbose = False):
        """
        Parameters
        ----------
        unique_id : integer
            A value that is unique to each lymphocyte in the simulation
        model : CancerModel class 
            Inherit properties from the main ABM class so any updates there 
            is also given to the agents 
        is_stem : Boolean
            Value to indicate whether this tumor cell is a stem cell or not 
        verbose : Boolean, optional
            If True, more information is printed for each action and step done by the 
            model. The default is False.

        """
        
        # pass parameters to the parent class 
        super().__init__(unique_id, model)
        
        # given properties of tumor cell
        self.is_stem = is_stem # stem cell? 0 = False, 1 = True 
        self.model = model 
        self.unique_id = 'TC'+str(unique_id)
        self.verbose = verbose 
        self.type = 'tumor'
        self.n_prol = 0 
        
        if self.verbose > 0:
            print(f'Hello! I am {self.unique_id}!')
        
        # add parameters to agent (note these parameters might change due to treatments)
        self.pars = self.model.pars.copy()
        self.set_parameter_values()
        
        # other properties of tumor cell
        self.damage = 0 # how much damage tumor cell sustained 
        self.engaged = False # whether it is being targeted by a T cell 
        

    def set_parameter_values(self):
        """Method to set parameter values, can be called again if parameter value changes"""
        self.TUpprol = self.pars['TUpprol'] # probability of proliferation
        self.TUpmig = self.pars['TUpmig'] # probability of migration
        self.TUpdeath = self.pars['TUpdeath'] # probability of cell death 
        self.TUpmax = self.pars['TUpmax'] # max proliferation capacity 
        self.TUps = self.pars['TUps'] # tumor stem cell probability of symmetric division 
        
    def step(self):
        """Method called when performing a step in the model"""
        
        # remove cells that died 
        self.model.clear_dying_cells()
        
        # agent is on the move : migrates, grows or dies
        self.go_grow_die()
        
        
    def get_possible_steps(self):
        """Get empty (neighboring) grid cells"""
        possible_steps = deque()
        
        for x in self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False):
            
            if self.model.grid.is_cell_empty(x):
                possible_steps.append(x)
                
        return possible_steps 
        
        
        
    def go_grow_die(self):
        """It's the turn for the agent to act: migrate, grow, die or idle"""
        # check if tumor cells get a disadvantage due to therapy
        if self.model.therapy_tumor:
            self.therapy_boost()
        
        # cell is killed by the immune system and dies spontaneously
        # only occurs if cell sustained enough damage OR (cell isn't a stem cell AND 
        # cell undergoes cell death)         
        if self.damage >= P.TUdamageThresh or (self.is_stem == 0 and random.random() < self.TUpdeath):
            self.model.kill_agents.append(self) # add agent to list of cells to remove
            
            # check if cell died on its own or was killed
            if self.engaged:
                self.model.N_tumor_killed += 1 
            else:
                self.model.N_tumor_died += 1 
            
            if self.verbose > 0:
                print(f'{self.unique_id} just died :)')
        
        else:
            # get only empty neighborhood spaces 
            possible_steps = self.get_possible_steps()
            
            # if there is space, check whether cell migrates or proliferates 
            if len(possible_steps) > 0 and not self.engaged:
                # check if cell proliferates 
                if random.random() < self.TUpprol and self.n_prol < self.TUpmax:
                    self.n_prol += 1 
                    
                    # assumption : stem cell divides into stem cell
                    # non stem cell divides into non stem cell 
                    # unless there is an asymmetric division (stem cell divides into non stem cell)
                    if self.is_stem == 0:
                        new_stem = 0 
                    else:
                        new_stem = int(random.random() < self.TUps)
                    
                    # add new tumor cell to model 
                    self.model.add_tumor_cell(self.model.N_tot_tumor, is_stem = new_stem,
                                              verbose = self.verbose, 
                                              loc = random.choice(possible_steps), 
                                              rand_loc = False
                                              )
                    # add cell to counter 
                    self.model.N_tumor += 1 
                    self.model.N_tot_tumor += 1 
                    
                    if self.verbose > 1:
                        print(f'{self.unique_id} just proliferated (stem cell = {new_stem})!')
            
                
                elif random.random() < self.TUpmig:
                    # move agent if it starts to migrate 
                    self.move()
                    
                    if self.verbose > 0:
                        print(f'{self.unique_id} is on the move...')
                
                else: 
                    # do nothing 
                    if self.verbose > 1:
                        print(f'{self.unique_id} was idle')
                    pass 
                
            else:
                # no space 
                if self.verbose > 1 and len(possible_steps) <=0:
                    print(f'{self.unique_id} was idle (no space)')
                elif self.verbose > 0 and self.engaged:
                    print(f'>> {self.unique_id} was idle (being killed)')
                pass  
            
            
    def move(self):
        """Move cell"""
        possible_steps = self.get_possible_steps()
        new_position = random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)
        
        
    def therapy_boost(self):
        """Model effect of therapy on tumor cell"""
        
        for adj_par in self.model.adj_par_tumor:
            for (key, value) in adj_par.items():
                if key in P.tumor_pars: # check if parameter is classified as tumor parameter
                    # change parameter value 
                    old = self.pars[key]
                    
                    if 'p' in key:
                        self.pars[key] = min(self.pars[key] * value, 1) 
                    else:
                        self.pars[key] = self.pars[key] * value
                    
                    if self.verbose > 0:
                        print(f'>> {self.unique_id} is affected by therapy! {key} changes from {old} to {self.pars[key]}')
             
        # (re)initialize parameter values with adjustments 
        self.set_parameter_values()
        
        
                                               
		

        
        
        
        
        