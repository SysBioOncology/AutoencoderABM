# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:15:29 2024

@author: 20192020
"""
import Parameters as P
from Parameters import mesa, random, np, deque, math
import matplotlib.pyplot as plt 
import seaborn as sns 
sns.set_style("white")

class SimResults:
    """Class to save some simulation statistics"""
    
    __slots__ = ("Ntumor_alive", "Nlymp_alive", "Ntumor_died", "Nlymp_died", 
                 "Ntumor_killed", "Nlymp_influx", "Nlymp_exhausted", 
                 "Ntumor_stem", "Ntumor_total", "Nlymp_total", "time", "current_time", 
                 "lymph_pos_x", "lymph_pos_y")
    
    def __init__(self):
        # create objects to save simulation stats
        self.Ntumor_alive = deque([0])
        self.Nlymp_alive = deque([0])
        self.Ntumor_died = deque([0]) # tumor cells dying on their own 
        self.Nlymp_died = deque([0])
        self.Ntumor_killed = deque([0])
        self.Nlymp_influx = deque([0])
        self.Nlymp_exhausted = deque([0])
        self.Ntumor_stem = deque([0])
        self.Ntumor_total = 0 
        self.Nlymp_total = 0 
        self.time = deque([0])
        self.current_time = 0
        self.lymph_pos_x = {}
        self.lymph_pos_y = {}
        
    def register_new_round(self):
        """Keep track of how long we are into the simulation"""
        self.time.append(self.current_time + P.oneStepDuration)
        self.current_time += P.oneStepDuration
        
    def register_Ncells(self, N_tumor, N_lymp):
        """Update how many cells are currently alive in the model"""
        self.Ntumor_alive.append(N_tumor)
        self.Nlymp_alive.append(N_lymp)
        
    def register_tumor_cell_died(self, died):
        """Add how many tumor cells died (not due to lymphocytes) this round"""
        self.Ntumor_died.append(died)
        
    def register_lymphocyte_died(self, died):
        """Add how many lymphocytes died this round"""
        self.Nlymp_died.append(died)
        
    def register_tumor_cell_killed(self, killed):
        """Add how many tumor cells were killed this round"""
        self.Ntumor_killed.append(killed)
        
    def register_lymphocyte_influx(self, influx):
        """Add how many lymphocytes were added this round"""
        self.Nlymp_influx.append(influx)
        
    def register_lymphocyte_exhausted(self, exhaust):
        """Add how many lymphocytes were exhausted this round"""
        self.Nlymp_exhausted.append(exhaust)
        
    def register_tumor_stem(self, n_stem):
        """Add how many of the tumor cells that are alive are stem cells"""
        self.Ntumor_stem.append(n_stem)
        
    def register_tumor_total(self, total):
        """Update number of tumor cells that have been in simulation"""
        self.Ntumor_total = total 
        
    def register_lymp_total(self, total):
        """Update number of lymphocytes that have been in simulation"""
        self.Nlymp_total = total 
        
    def register_lymph_movement(self, lymph_id, pos):
        """Update tracker for lymphocyte movement"""
        
        if lymph_id not in self.lymph_pos_x.keys():
            self.lymph_pos_x[lymph_id] = deque([pos[0]])
            self.lymph_pos_y[lymph_id] = deque([pos[1]])
        else:
            self.lymph_pos_x[lymph_id].append(pos[0])
            self.lymph_pos_y[lymph_id].append(pos[1])
            
        
        
    def get_total_tumor(self):
        """Return total number of tumor cells that are / have been in simulation"""
        return(self.Ntumor_total)
    
    def get_total_lymp(self):
        """Return total number of lymphocytes that are / have been in simulation"""
        return(self.Nlymp_total)
        
    def get_final_tumor(self):
        """Return final number of tumor cells in simulation"""
        return(self.Ntumor_alive[-1])
    
    def get_final_lymp(self):
        """Return final number of lymphocytes in simulation"""
        return(self.Nlymp_alive[-1])
        
    def calculate_lymp_tumor_ratio(self):
        """Determine ratio lymphocytes : tumor cells"""
        if self.Ntumor_alive[-1] > 0:
            output = self.Nlymp_alive[-1] / self.Ntumor_alive[-1]
        else:
            output = 0
        
        return(output)
    
    def calculate_killed_tumor_ratio(self):
        """Determine ratio total number of killed tumor : total number of tumor"""
        return sum(self.Ntumor_killed) / self.Ntumor_total
    
    def calculate_exhaust_lymp_ratio(self):
        """Determine ratio exhausted lymphocytes : lymphocytes"""
        if self.Nlymp_alive[-1] > 0:
            output = self.Nlymp_exhausted[-1] / self.Nlymp_alive[-1]
        else:
            output = 0
            
        return(output)
    
    def calculate_stem_tumor_ratio(self):
        """Determine ratio tumor stem cells : tumor cells"""
        if self.Ntumor_alive[-1] > 0:
            output = self.Ntumor_stem[-1] / self.Ntumor_alive[-1]
        else:
            output = 0
        
        return(output)
        

    def visualize_cell_counts(self, plot_title = None, maxTime = 1200):
        """Visualize saved properties from simulation"""
        
        # six plots to visualize some simulation results 
        fig, axes = plt.subplots(2, 4)
        fig.set_size_inches(10, 6)
    
    
        if P.therapy_administration is not None:
            plot_therapy = [x * P.oneStepDuration for x in P.therapy_administration]
            plot_therapy = [x for x in plot_therapy if x <= max(self.time)]
            
            axes[0,0].vlines(plot_therapy,
                             min(self.Ntumor_alive), max(self.Ntumor_alive),
                             color='black', linestyle='--', alpha=0.7)
            axes[0,1].vlines(plot_therapy,
                             min(self.Ntumor_died), max(self.Ntumor_died),
                             color='black', linestyle='--', alpha=0.7)
            axes[0,2].vlines(plot_therapy,
                             min(self.Ntumor_killed), max(self.Ntumor_killed),
                             color='black', linestyle='--', alpha=0.7)
            axes[0,3].vlines(plot_therapy,
                             min(self.Ntumor_stem), max(self.Ntumor_stem),
                             color='black', linestyle='--', alpha=0.7)
            axes[1,0].vlines(plot_therapy,
                             min(self.Nlymp_alive), max(self.Nlymp_alive),
                             color='black', linestyle='--', alpha=0.7)
            axes[1,1].vlines(plot_therapy,
                             min(self.Nlymp_died), max(self.Nlymp_died),
                             color='black', linestyle='--', alpha=0.7)
            axes[1,2].vlines(plot_therapy,
                             min(self.Nlymp_influx), max(self.Nlymp_influx),
                             color='black', linestyle='--', alpha=0.7)
            axes[1,3].vlines(plot_therapy,
                             min(self.Nlymp_exhausted), max(self.Nlymp_exhausted),
                             color='black', linestyle='--', alpha=0.7)
            
        # Ntumor - alive 
        axes[0,0].set_title('Alive')
        axes[0,0].set_ylabel('Tumor cells', weight = 'bold', fontsize=18)
        #axes[0,0].set_xlim(0, maxTime)
        sns.lineplot(x=self.time, y=self.Ntumor_alive, ax=axes[0,0],
                     color = 'black')
        
        # Ntumor - died 
        axes[0,1].set_title('Died')
        #axes[0,1].set_xlim(0, maxTime)
        sns.lineplot(x=self.time, y=self.Ntumor_died, ax=axes[0,1],
                     color = 'black')
        
        # Ntumor - killed 
        axes[0,2].set_title('Killed')
        #axes[0,2].set_xlim(0, maxTime)
        sns.lineplot(x=self.time, y=self.Ntumor_killed, ax=axes[0,2],
                     color = 'black')
        
        # Ntumor - stem cells 
        axes[0,3].set_title('Stem cells')
        sns.lineplot(x=self.time, y=self.Ntumor_stem, ax=axes[0,3],
                     color = 'black')
        
        # Nlymphocytes - alive 
        axes[1,0].set_title('Alive')
        axes[1,0].set_ylabel('Lymphocytes', weight = 'bold', fontsize=18)
        axes[1,0].set_xlabel('Time (h)')
        #axes[1,0].set_xlim(0, maxTime)
        sns.lineplot(x=self.time, y=self.Nlymp_alive, ax=axes[1,0],
                     color = 'black')
        
        # Nlymphocytes - died 
        axes[1,1].set_title('Died')
        axes[1,1].set_xlabel('Time (h)')
        #axes[1,1].set_xlim(0, maxTime)
        sns.lineplot(x=self.time, y=self.Nlymp_died, ax=axes[1,1],
                     color = 'black')
        
        # Nlymphocytes - influx 
        axes[1,2].set_title('Influx')
        axes[1,2].set_xlabel('Time (h)')
        #axes[1,2].set_xlim(0, maxTime)
        sns.lineplot(x=self.time, y=self.Nlymp_influx, ax=axes[1,2],
                     color = 'black')
        
        # Nlymphocytes - exhausted 
        axes[1,3].set_title('Exhausted')
        axes[1,3].set_xlabel('Time (h)')
        sns.lineplot(x=self.time, y=self.Nlymp_exhausted, ax=axes[1,3],
                     color = 'black')
        
        plt.tight_layout()
        
        if plot_title != None:
            plt.savefig(plot_title,dpi=600)
        else:
            plt.savefig(f'../output/SimulationCounts.png',dpi=600)
            
    def visualize_lymph_movement(self, filepath = None):
        """Visualize tracking lines of each lymphocyte"""
        # prepare plot 
        plt.figure()
        
        # plot every cell 
        for cell in self.lymph_pos_x.keys():
            x = self.lymph_pos_x[cell]
            y = self.lymph_pos_y[cell]
        
            plt.plot(y, x, 'k-', linewidth=0.4, alpha=0.4)
            
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Lymphocyte tracking data')
        plt.xlim(0, P.width)
        plt.ylim(P.height, 0)
        
        if filepath != None:
            plt.savefig(filepath,dpi=600)
        else:
            plt.savefig(f'../output/SimulationCounts.png',dpi=600)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
