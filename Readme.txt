PACKAGE CONTENTS
Soft-Gluons is a project to compute the eikonal contribution of soft gluon radiation in (up to now) 2 -> 3 parton scattering. To validate the result, the ratio R_s of the eikonal contribution divided by the actual matrix element (taken from Sherpa) is calculated.

To this date, Soft-Gluons contains the following files:
> interface.py
> eikonal.py
> eikonal.c

USAGE
The 'interface.py' is an altered version of the Sherpa Python interface contained in the Sherpa package (sherpa.hepforge.org) which is used to compute the colour summed, squared matrix element of a given process.
There are currently two codes to compute the eikonal contribution of soft gluon radiation: 'eikonal.py' and 'eikonal.c'. While 'eikonal.c' isn't interfaced to Sherpa yet, the python one is. Both codes use given hard and soft anomalous dimension matrices to compute the eikonal contribution.

Usage of eikonal.py:
In 'eikonal.py', two functions are defined. The function 'writeRunCard' saves a file 'Run.dat' that is needed by Sherpa and further defines the process and configurations. Its arguments are the process identifiers (see below) and the centre of mass frame energy of one parton beam, where a symmetric collides is considered. 
The function 'calcTrace' calculates the the trace of the altered hard matrix of the underlying 2 -> 2 parton process when a soft gluon is emitted. Its arguments are the momenta p_1, p_2, p_3, p_4 of the hard 2 -> 2 process, the soft momentum p_s of the soft gluon and the process identifier (see below). 
The function 'calcMatEl' is defined in 'interface.py' and returns the actual 2 -> 3 squared matrix element taken from Sherpa. Its arguments are again the particle momenta p_1, p_2, p_3, p_4 and p_s.
In 'eikonal.py', an array of lambda_s and R_s is constructed, which can be printed or used for plotting.

Process identifiers:
Process identifiers are used to distinguish various hard 2 -> 2 processes. The process variable is a tuple [process id, subprocess id], where process id and subprocess id are integers with the following identification:
process id      subprocess id       parton process
1               1                   q qbar -> q qbar (same flavour overall)
1               2                   q qbar -> q' qbar' (same flavours in initial and final states respectively)
1               3                   q qbar' -> q qbar' (mixed flavours)
2               1                   q q -> q q (same flavour)
2               2                   q q' -> q q' (mixed flavour)
3               not needed          q g -> q g
4               1                   q qbar -> g g
4               2                   g g -> q qbar
5               not needed          g g -> g g
As it can be seen from the table, the process ids group similar processes and the subprocess ids are used to distinguish minor differences among them. The subprocess ids which are 'not needed' may just be set to 1 or any other integer value.

INSTALLATION
To work properly, Sherpa needs to be set up with '--enable-pyext' in the ./configure prompt and the files 'eikonal.py' and 'interface.py' need to be placed in the ME2-Python folder created by Sherpa: /path-to-sherpa/share/SHERPA-MC/Examples/API/ME2-Python/, where 'path-to-sherpa denotes the path to which Sherpa was installed.

AUTHORS
Jose Llanes Jurado, Christian Tobias Preuss
