# -*- coding: utf-8 -*-
"""
Create on Fri Oct 22 13:58:05 2021

@author: BenJanevic
"""

class Source:
    def __init__(self, index, x, y, Q, fixed=False):
        self.index = index
        self.x = x
        self.y = y
        self.Q = Q
        self.fixed = fixed
        

    def __repr__(self):
        return f"({self.x}, {self.y}): {round(self.Q, 2)}"
    
    
    def is_fixed(self):
        return self.fixed