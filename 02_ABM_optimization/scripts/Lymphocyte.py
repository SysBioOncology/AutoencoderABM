# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:34:59 2024

@author: 20192020
"""
import Parameters as P
from Parameters import mesa, random, deque, math, np 

class LymphocyteAgent(mesa.Agent):
    """Agent to describe behavior of a lymphocyte"""
    
    __slots__ = ("model", "unique_id", "verbose", "type", "pars", "IMpkill_orig", 
                 "damage", "engaged", "n_proliferated", "n_killed", "tumor_target", "countdown",
                 "tumor_map", "pos")
  
    def __init__(self, unique_id, model, verbose = False):
        """
        Parameters
        ----------
        unique_id : integer
            A value that is unique to each lymphocyte in the simulation
        model : CancerModel class 
            Inherit properties from the main ABM class so any updates there 
            is also given to the agents 
        verbose : Boolean, optional
            If True, more information is printed for each action and step done by the 
            model. The default is False.

        """
        # pass parameters to the parent class 
        super().__init__(unique_id, model)
        
        # given properties of tumor cell
        self.model = model 
        self.unique_id = 'LYM'+str(unique_id)
        self.verbose = verbose 
        self.type = 'lymphocyte'
        
        if self.verbose > 0:
            print(f'Hello! I am {self.unique_id}!')
        
        # add parameters to agent
        self.pars = self.model.pars.copy()
        self.set_parameter_values()
        
        # keep some original parameter values 
        self.IMpkill_orig = model.pars['IMpkill']
        
        # other properties of tumor cell
        self.damage = 0 # how much damage tumor cell sustained 
        self.engaged = False # whether cell is engaged 
        self.n_proliferated = 0 # keep track of how often cell proliferated
        self.n_killed = 0 # keep track of how many tumor cells it has killed
        self.tumor_target = None # tumor cell that the lymphocyte is currently killing
        self.countdown = 0 # tracker to keep track when immune cell is ready for action again 
        
        # tumor map used if directed migration is active (1 = tumor, 0 = no tumor)
        self.tumor_map = np.zeros((P.width, P.height)) 
        
    def set_parameter_values(self):
        """Method to set parameter values, can be called again if parameter value changes"""
        self.IMkmax = self.pars['IMkmax']
        self.IMpmax = self.pars['IMpmax']
        self.IMpmig = self.pars['IMpmig']
        self.IMpkill = self.pars['IMpkill']
        self.IMpprol = self.pars['IMpprol'] 
        self.IMpdeath = self.pars['IMpdeath']
        self.IMrwalk = self.pars['IMrwalk']
        self.IMinfluxEdge = self.pars['IMinfluxEdge']
        self.IMdirectedWidth = self.pars['IMdirectedWidth']
        self.IMdecayWeight = self.pars['IMdecayWeight']
        self.IMdecay = self.pars['IMdecay']
        
    def step(self):
        """Method called when performing a step in the model"""
        # remove cells that died 
        self.model.clear_dying_cells()
        
        # agent is on the move : migrates, attacks, grows or dies
        self.action()
        
    def get_possible_steps(self):
        """Get empty (neighboring) grid cells"""
        possible_steps = deque()
        
        for x in self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False):
            
            if self.model.grid.is_cell_empty(x):
                possible_steps.append(x)
                
        return possible_steps 
        
    def action(self):
        """It's the turn for the agent to act: migrate, attack, grow, die or idle"""
        skipAttack = False
        
        # check if immune cell receives boost from immunotherapy 
        if self.model.therapy_lymph:
            self.therapy_boost()
        
        # check if engaged with a tumor cell 
        if self.engaged:
            # countdown until lymphocyte has dealt with the tumor cell accordingly 
            self.countdown -= 1
            
            # update damage done to the tumor cell 
            self.tumor_target.damage += P.TUdamageThresh / (
                math.ceil(P.engagementDuration / P.oneStepDuration))
            
            if self.verbose > 0:
                print(f'>> {self.unique_id} is currently killing {self.tumor_target.unique_id}!')
            
            # check if cell should be done with its engagement 
            if self.countdown == 0:
                
                if self.verbose > 0:
                    print(f'>> {self.unique_id} successfully killed {self.tumor_target.unique_id}!')
                
                self.engaged = False 
                self.tumor_target = None 
        
        else: 
            # cell dies 
            if random.random() < self.IMpdeath:
                self.model.kill_agents.append(self) # add agent to list of cells to remove
                self.model.N_lymp_died += 1 
                
                if self.verbose > 0:
                    print(f'{self.unique_id} just died :(')
                skipAttack = True 
                
            else: 
                # get only empty neighborhood spaces 
                possible_steps = self.get_possible_steps()
                
                # if there is space, check whether cell migrates or proliferates 
                if len(possible_steps) > 0:
                    # check if cell proliferates and proliferation capacity has not been exceeded
                    if random.random() < self.IMpprol and self.n_proliferated <= self.IMpmax:
                        
                        # add new lymphocyte to model 
                        self.model.add_lymphocyte_cell(self.model.N_tot_lymp,
                                                  verbose = self.verbose, 
                                                  loc = random.choice(possible_steps), 
                                                  rand_loc = False
                                                  )
                        # add cell to counter 
                        self.model.N_lymp += 1 
                        self.model.N_tot_lymp += 1 
                        self.n_proliferated += 1 
                        
                        if self.verbose > 1:
                            print(f'{self.unique_id} just proliferated!')
                            
                        skipAttack = True 
                        
                    elif random.random() < self.IMpmig:
                        # move agent if it starts to migrate 
                        self.move()
                        
                        if self.verbose > 0:
                            print(f'{self.unique_id} is on the move...')
                    
                    else: 
                        # do nothing / idle 
                        if self.verbose > 1:
                            print(f'{self.unique_id} was idle')
                        pass 
                    
                else:
                    # no space so idle 
                    if self.verbose > 1:
                        print(f'{self.unique_id} was idle (no space)')
                    pass  
                
        
            # lymphocyte can kill cell if it didn't proliferate or die 
            if not skipAttack:
                # if cell is ready to attack (aka has not exceeded kill capacity)
                if self.n_killed < self.IMkmax and random.random() < self.IMpkill:
                    
                    # check if there is a tumor cell that can be killed as neighbor
                    target = self.find_target()
                    
                    if target != 'nothing':
                        if self.verbose > 0: 
                            print(f'>> {self.unique_id} started targeting {target.unique_id}!')
                        
                        self.n_killed += 1  
                        
                        # check if immune cell is now exhausted 
                        if self.n_killed >= self.IMkmax:
                            self.model.N_lymp_exhaust += 1 
                            
                        self.engaged = True 
                        self.tumor_target = target
                        self.countdown = math.ceil(P.engagementDuration / P.oneStepDuration)
                        
                        # update damage done to the tumor cell 
                        self.tumor_target.damage += P.TUdamageThresh / (
                            math.ceil(P.engagementDuration / P.oneStepDuration))
                        self.tumor_target.engaged = True 
                    
                    else:
                        if self.verbose > 0:
                            print(f'>> {self.unique_id} wants to kill but there are no targets')
                
                
            
        
    def move(self):
        """Move cell"""
        possible_steps = self.get_possible_steps()
        
        # decide if movement will be completely random or towards area with most tumor cells
        if random.random() < self.IMrwalk: # random movement
            new_position = random.choice(possible_steps)
        else: # directed movement 
            # update tumor map 
            self.update_tumor_map()
        
            # use tumor count in each direction to get new pos 
            new_position, tumor_count = self.get_best_direction(possible_steps)
            
            if self.verbose > 0:
                print(f'>> {self.unique_id} is moving with purpose! Going to {new_position} with {tumor_count} tumors in sight')
                
        self.model.grid.move_agent(self, new_position)
        
        # register new position 
        self.model.sim_results.register_lymph_movement(self.unique_id, self.pos)
        
    def update_tumor_map(self):
        """Update Boolean mask on where tumor cells are on the grid"""
        # reset map 
        self.tumor_map.fill(0)
        
        # update map with tumor movements 
        for cell_content, (x, y) in self.model.grid.coord_iter():
            if cell_content != None and cell_content.type == 'tumor':
                self.tumor_map[x, y] = 1
        
    def get_best_direction(self, possible_steps):
        """Get tumor count in a straight line in each direction with possible movement"""
        best_pos = possible_steps[0] 
        best_count = 0 
        
        # shuffle possible_steps so avoid bias to first position having highest 
        # possibility of being chosen if all equal 
        random.shuffle(possible_steps)
        
        # go over all directions of possible steps 
        for p in possible_steps:
            
            # check if neighbor is on diagonal line
            if p[0] != self.pos[0] and p[1] != self.pos[1]:
                # check in which diagonal direction neighbor is 
                if p[0] < self.pos[0] and p[1] < self.pos[1]: #top, left  
                    line_pos_x = [xi for xi in reversed(range(0, p[0]))] 
                    line_pos_y = [yi for yi in reversed(range(0, p[1]))]
                    direction = ('diagonal', 'top', 'left')
                    
                elif p[0] > self.pos[0] and p[1] < self.pos[1]: #bottom, left
                    line_pos_x = [xi for xi in range(p[0], P.height)] 
                    line_pos_y = [yi for yi in reversed(range(0, p[1]))]
                    direction = ('diagonal', 'bottom', 'left')
                    
                elif p[0] < self.pos[0] and p[1] > self.pos[1]: #top, right
                    line_pos_x = [xi for xi in reversed(range(0, p[0]))] 
                    line_pos_y = [yi for yi in range(p[1], P.width)]
                    direction = ('diagonal', 'top', 'right')
                    
                else: #bottom, right
                    line_pos_x = [xi for xi in range(p[0], P.height)] 
                    line_pos_y = [yi for yi in range(p[1], P.width)]
                    direction = ('diagonal', 'bottom', 'right')
                    
                    
            # check if neighbor is on straight line 
            else:
                # check if neighbor is top/bottom/left/right
                if p[0] > self.pos[0]:  # bottom
                    line_pos_x = [xi for xi in range(p[0], P.height)] 
                    line_pos_y = [self.pos[1] for _ in range(p[0], P.height)]
                    direction = ('straight', 'bottom')
                    
                elif p[0] < self.pos[0]: # top
                    line_pos_x = [xi for xi in reversed(range(0, p[0]))] 
                    line_pos_y = [self.pos[1] for _ in reversed(range(0, p[0]))] 
                    direction = ('straight', 'top')
                    
                else: 
                    if p[1] > self.pos[1]: # right
                        line_pos_x = [self.pos[0] for _ in range(p[1], P.width)] 
                        line_pos_y = [yi for yi in range(p[1], P.width)]
                        direction = ('straight', 'right')

                    else: # left
                        line_pos_x = [self.pos[0] for _ in reversed(range(0, p[1]))] 
                        line_pos_y = [yi for yi in reversed(range(0, p[1]))]
                        direction = ('straight', 'left')
                        
            # get tumor count 
            # ensure x and y lists are same length 
            line_pos_x = line_pos_x[0:len(line_pos_y)]
            line_pos_y = line_pos_y[0:len(line_pos_x)]
            
            tumor_counts = self.tumor_map[line_pos_x, line_pos_y]
                
            # multiply with weights depending on distance 
            # ensure also that P.IMdecay has same length as tumor_counts
            if P.IMdirectedWidth > 0: # if more neighbors are taken into account to avoid blind spots
               tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
               current_count = sum(tumor_counts)
               
               # get neighboring cells of central line 
               for step_size in range(1, P.IMdirectedWidth+1):
                   current_count += self.get_wider_direction(direction, line_pos_x, line_pos_y, 
                                                        tumor_counts, step_size)
                         
            else:
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                current_count = sum(tumor_counts)
                
            # set cutoff value so current_count < threshold is just 0 
            # avoids something E-30 being seen as largest 
            if current_count < 0.01: 
                current_count = 0 
            
            if current_count > best_count: # found new best 
                best_count = current_count
                best_pos = p
            elif current_count == best_count: # if equal, choose random one
                if random.choice([1, 2]) == 2:
                    best_count = current_count
                    best_pos = p 
                    
        return best_pos, best_count
        
        
    def get_wider_direction(self, direction, line_pos_x, line_pos_y, tumor_counts, step_size):
        """Widen view of tumors to take into account for directed migration"""
        additional_count = 0 
        
        # skip a few coordinates if resulting position end up behind current cell position 
        skip = max(step_size-2, 0)
        
        # define direction
        # take into account that side steps might fall out of grid 
        if direction[0] == 'diagonal':
            # get both left and right neighbor of central line 
            if direction[1] == 'bottom':
                side_pos_x = [x-step_size for x in line_pos_x if x-step_size >= 0][skip:]
                tumor_counts = self.tumor_map[side_pos_x, line_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
            else: # top 
                side_pos_x = [x+step_size for x in line_pos_x if x+step_size < P.height][skip:]
                tumor_counts = self.tumor_map[side_pos_x, line_pos_y[:len(side_pos_x)]][skip:]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
            
            if direction[2] =='left':
                side_pos_y = [y+step_size for y in line_pos_y if y+step_size < P.width][skip:]
                tumor_counts = self.tumor_map[line_pos_x[:len(side_pos_y)], side_pos_y]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
            else: #right
                side_pos_y = [y-step_size for y in line_pos_y if y-step_size >= 0][skip:]
                tumor_counts = self.tumor_map[line_pos_x[:len(side_pos_y)], side_pos_y]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
        else: # straight
            if direction[1] == 'bottom':
                side_pos_x = [x-step_size for x in line_pos_x if x-step_size >= 0][skip:]
                side_pos_y = [y-step_size for y in line_pos_y if y-step_size >= 0][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                side_pos_y = [y+step_size for y in line_pos_y if y+step_size < P.width][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                
            elif direction[1] == 'left':
                side_pos_x = [x-step_size for x in line_pos_x if x-step_size >= 0][skip:]
                side_pos_y = [y+step_size for y in line_pos_y if y+step_size < P.width][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                side_pos_x = [x+step_size for x in line_pos_x if x+step_size < P.height][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                
            elif direction[1] == 'right':
                side_pos_x = [x-step_size for x in line_pos_x if x-step_size >= 0][skip:]
                side_pos_y = [y-step_size for y in line_pos_y if y-step_size >= 0][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                side_pos_x = [x+step_size for x in line_pos_x if x+step_size < P.height][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                
            else: # top 
                side_pos_x = [x+step_size for x in line_pos_x if x+step_size < P.height][skip:]
                side_pos_y = [y-step_size for y in line_pos_y if y-step_size >= 0][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
                side_pos_y = [y+step_size for y in line_pos_y if y+step_size < P.width][skip:]
                tumor_counts = self.tumor_map[side_pos_x[:len(side_pos_y)], side_pos_y[:len(side_pos_x)]]
                tumor_counts = np.multiply(tumor_counts, P.IMdecay[:len(tumor_counts)])
                additional_count += sum(tumor_counts)
                
        return additional_count
        
    def find_target(self):
        """Find tumor cell to kill"""
        target = 'nothing' 
        
        # get neighbors of this lymphocyte 
        neighbors = self.model.grid.get_neighbors(self.pos, moore=True, include_center=False)
        
        if len(neighbors) > 0:
            # randomly shuffle list 
            random.shuffle(neighbors)
            
            # get the first random neighbor that is a tumor cell 
            # and the tumor cell is not being killed 
            for one_neighbor in neighbors:
                if one_neighbor.type == 'tumor' and not one_neighbor.engaged:
                    target = one_neighbor 
                    break 
            
        return target
    
    def therapy_boost(self):
        """Improve lymphocyte killing of tumor cells after immunotherapy"""
        
        for adj_par in self.model.adj_par_lymph:
            for (key, value) in adj_par.items():
                if key in P.lymph_pars: # check if parameter is classified as lymphocyte parameter
                    # change parameter value 
                    old = self.pars[key]
                    
                    if 'p' in key:
                        self.pars[key] = min(self.pars[key] * value, 1) 
                    else:
                        self.pars[key] = self.pars[key] * value
                    
                    if self.verbose > 0:
                        print(f'>> {self.unique_id} received a boost in performance! {key} improves from {old} to {self.pars[key]}')
             
        # (re)initialize parameter values with adjustments 
        self.set_parameter_values()

        
        
        
        
        
    
