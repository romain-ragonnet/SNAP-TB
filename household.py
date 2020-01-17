from numpy import random

class household:
    def __init__(self, id):
        self.id = id
        self.size = 0
        self.individual_ids = [] # list containing the ids of all individuals currently living in the household
        self.school_id = None
        self.repopulate_date = 0.  # date at which a new couple moves in
        self.minimum_time_to_next_baby = random.uniform(low=0., high=365.25*4)
        self.last_baby_time = 0  # most recent date at which a last baby was born in the household
