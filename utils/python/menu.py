import sys

def menu():
    message = '''
    ----------------------------------------------------------------------------------------
                            MENU: CHOOSE ONE OF THE FOLLOWING OPTIONS
    ----------------------------------------------------------------------------------------
    
    1   Run simulation
    2   Extract and plot data from simulation
    3   Compute topological charge and visualize spin vectors
    '''
    
    try:
        usrChoice = int(input(message))
        return usrChoice
    
    except ValueError:
        sys.exit('Exiting... The input was not an integer. Run and try again.')
     