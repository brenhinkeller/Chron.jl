# Chron.jl

A two-part framework for (1) estimating eruption/deposition age          
distributions from complex mineral age spectra and (2) subsequently building  
a stratigraphic age model based on those distributions. Each step relies on   
a Markov-Chain Monte Carlo model.                                             

The first model uses an informative prior distribution to estimate the      
times of first (i.e., saturation) and last  mineral crystallization (i.e.,    
eruption/deposition).                                                         

The second model uses the estimated (posterior) eruption/deposition age
distributions along with the constraint of stratigraphic superposition to     
produce an age-depth model                                                      
