The original dataset is ds_filtered conducted by DisGeNet group.

Constructing_Graphs.py generates two networks based on Formula 1 in the manuscript including disease network and gene networks. It reads the file “ds_filtered” as an input. The output for the disease network is disease_network_1.
To preserve the multi-scale backbone of the weighted human disease network (WHDN) while removing less relevant and meaningful edges we use a multi-scale filtering method proposed by Serrano et al. [2009].

The result is the file “r_disease_network_1” and the giant component of the deducted network is given in “Giantr_disease_network_1”. 

Giantr_disease_network_1 is the file we apply the measures including DC, CC, BC, and W-DIL methods on it. In this file there are three columns that shows the given weight between two diseases.  

DIL_Method_W.py calculate the node importance scores calculated by DIL-W method. It reads file “Giantr_disease_network_1” as an input and gives scores to the nodes in the disease network. If you want to compare the DIL-W result with DC, CC, and BC the proper lines in the file must be uncommented. 
