# DIAS-Summer-Internship
I added my c++ files for example 3. I tried to put some descriptive comments in there. If you want to run any of the code then I think the only thing
you will need to change is the file paths that are specified because they are all written to refer to my local files. Other than that the code should run fine. 
The update function might look a bit weird. I originally wrote it where I computed the hamiltonian twice (one at t and one at t + dt) for each update. 
Then I realised I could make the code more efficient instead of computing the hamiltonian at t and t+dt I could just save the value of H(t+dt) into 
the variable holding H(t) so that on each run I only have to compute H(t+dt) and then H(t) will be assigned to whatever H(t+dt) was in the last iteration. 
If this doesn't make sense I can explain it more. Maybe there's an even better way to do it with less function calls again, but this method made the code twice
as fast as it was when I was calling the function twice per iteration which makes sense because now it's only half as many function calls. Let me know if anything
is weird or unclear and I'll be happy to help as much as I can if you're gonna start writing in c++. I haven't really written much for the D0 brane stuff, I've only
defined all my matrix functions like adding, multiplying, commutator, etc. When I have the simulation going I'll upload it into this file. 
